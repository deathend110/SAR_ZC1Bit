function out = zc1bit_screen_impl(cfg)
%ZC1BIT_SCREEN_IMPL 阶段一（置零规则快筛）单核内存优先工作流实现。
% 严格仅在 35 个校准样本（奇数 LocalSampleIdx）上执行，绝不触碰测试集。
% 设计要点：
%   1. 单核顺序执行（不使用 MATLAB 并行池），避免进程间大矩阵拷贝与内存带宽争抢；
%   2. 先扫描扁平检查点 checkpoints/screen/<ThresholdKey>.mat，全命中时 0 次读取回波文件；
%   3. 存在未完成项时对回波文件仅顺序读 1 轮（35 样本约 403 MB 常驻内存），并外层按 (Range_q, Azimuth_q)
%      在 cfg.runtime.cache_gb（推荐 4-6 GB）预算内复用 signal_up，最后按 ThresholdKey 原子落盘 -v7。

if nargin == 0
    cfg = zc1bit_config();
end
directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>

record = zc1bit_prepare_run(cfg);
p = zc1bit_paths(cfg);
S60 = load(cfg.parameter_file_60);
[threshold_groups, variant_specs] = zc1bit_build_screen_tasks(cfg, S60);

manifest = record.manifest;
cal_manifest = manifest(manifest.Split == "calibration", :);
if isempty(cal_manifest)
    error("zc1bit:EmptyCalibration", "样本清单中未找到校准样本。");
end

[cache_gb, preload_signals] = resolve_runtime_options(cfg);
cache_budget_bytes = cache_gb * (1024 ^ 3);

use_fast = true;
if isfield(cfg.screen, "use_fast_path")
    use_fast = logical(cfg.screen.use_fast_path);
end

n_tg = height(threshold_groups);
n_cal = height(cal_manifest);
n_var = height(variant_specs);

% 1. 先扫描全部扁平 ThresholdKey 检查点，仅汇总未完成任务
detail_by_tg = cell(n_tg, 1);
completed_tg = false(n_tg, 1);
cp_dir = fullfile(p.checkpoints, "screen");

for tg_idx = 1:n_tg
    tg = threshold_groups(tg_idx, :);
    cp_path = fullfile(cp_dir, tg.ThresholdKey + ".mat");
    if ~isfile(cp_path)
        continue;
    end
    state = zc1bit_load_state(cp_path, struct());
    if isfield(state, "completed") && isequal(state.completed, true) && isfield(state, "detail_rows")
        expected_params = build_tg_expected_params(tg, variant_specs, cal_manifest, S60, cfg);
        zc1bit_check_parameters(state.parameters, expected_params, ...
            sprintf("阶段一检查点(%s)", tg.ThresholdKey));
        detail_by_tg{tg_idx} = state.detail_rows;
        completed_tg(tg_idx) = true;
    end
end

% 2. 按需加载或构建校准集 GT 缓存；仅当存在未完成项或 GT 缓存缺失时顺序读盘 1 轮
gt_cache_path = fullfile(cp_dir, "gt_cal_cache.mat");
gt_cache = try_load_cal_gt_cache(cfg, S60, cal_manifest, gt_cache_path);
cal_signals = {};

if any(~completed_tg) || isempty(gt_cache)
    cal_signals = preload_all_calibration_signals(cal_manifest, S60, cfg.raw_decimation_step);
    if isempty(gt_cache)
        gt_cache = build_and_save_cal_gt_cache(cfg, S60, cal_manifest, cal_signals, gt_cache_path);
    end
end

if ~isfolder(p.screen)
    mkdir(p.screen);
end
writetable(gt_cache.roi_manifest, fullfile(p.screen, "screen_cal_enl_roi_manifest.csv"));

fprintf("\n=== 阶段一规则快筛：%d 阈值配置 × %d 变体（%d 校准样本 | 已命中 %d/%d 检查点 | cache=%.1f GB）===\n", ...
    n_tg, n_var, n_cal, sum(completed_tg), n_tg, cache_gb);

% 3. 仅对未完成的 ThresholdKey 按 (Range_q, Azimuth_q) 外层分组并复用上采样回波
pending_indices = find(~completed_tg);
if ~isempty(pending_indices)
    fast_ctx = [];
    if use_fast
        fast_ctx = build_fast_focus_context( ...
            S60, threshold_groups.Range_q(pending_indices), threshold_groups.Azimuth_q(pending_indices));
    end

    base_info = whos("cal_signals");
    base_bytes = base_info.bytes;
    if isa(cal_signals{1}, "single")
        elem_bytes = 8;
    else
        elem_bytes = 16;
    end

    pending_allocs = unique([ ...
        threshold_groups.Range_q(pending_indices), ...
        threshold_groups.Azimuth_q(pending_indices)], "rows", "stable");

    for a_idx = 1:size(pending_allocs, 1)
        rq = pending_allocs(a_idx, 1);
        aq = pending_allocs(a_idx, 2);
        group_tg_idx = pending_indices( ...
            threshold_groups.Range_q(pending_indices) == rq & ...
            threshold_groups.Azimuth_q(pending_indices) == aq);
        grid = resolve_sampling_grid(rq, aq, [S60.nrn, S60.nan], S60);

        if rq == 1 && aq == 1
            % Q=1 零上采样开销，直接复用内存中的基带回波
            for tg_idx = group_tg_idx(:).'
                detail_rows = run_tg_with_up_signals( ...
                    cfg, S60, grid, threshold_groups(tg_idx, :), variant_specs, ...
                    cal_manifest, cal_signals, gt_cache, fast_ctx, use_fast);
                save_tg_checkpoint(cp_dir, threshold_groups(tg_idx, :), variant_specs, ...
                    cal_manifest, S60, cfg, detail_rows);
                detail_by_tg{tg_idx} = detail_rows;
            end
            continue;
        end

        est_up_bytes = n_cal * prod(grid.upsampled_size) * elem_bytes;
        if preload_signals && (base_bytes + est_up_bytes <= cache_budget_bytes)
            fprintf("  [R%dA%d] 内存缓存模式：预构建 %d 个上采样回波（约 %.2f GB）...\n", ...
                rq, aq, n_cal, est_up_bytes / (1024 ^ 3));
            up_signals = cell(n_cal, 1);
            for k = 1:n_cal
                up_signals{k} = two_dim_upsample_fft(cal_signals{k}, rq, aq);
            end
            for tg_idx = group_tg_idx(:).'
                detail_rows = run_tg_with_up_signals( ...
                    cfg, S60, grid, threshold_groups(tg_idx, :), variant_specs, ...
                    cal_manifest, up_signals, gt_cache, fast_ctx, use_fast);
                save_tg_checkpoint(cp_dir, threshold_groups(tg_idx, :), variant_specs, ...
                    cal_manifest, S60, cfg, detail_rows);
                detail_by_tg{tg_idx} = detail_rows;
            end
            clear up_signals
        else
            fprintf("  [R%dA%d] 流式低内存模式：预计整组上采样 %.2f GB 超出剩余预算，逐样本复用（%d 组）...\n", ...
                rq, aq, est_up_bytes / (1024 ^ 3), numel(group_tg_idx));
            n_g = numel(group_tg_idx);
            tg_parts = cell(n_g, n_cal);
            for k = 1:n_cal
                sig_up = two_dim_upsample_fft(cal_signals{k}, rq, aq);
                sample_row = cal_manifest(k, :);
                for m = 1:n_g
                    tg = threshold_groups(group_tg_idx(m), :);
                    tg_parts{m, k} = evaluate_one_sample_for_tg( ...
                        cfg, S60, grid, tg, variant_specs, sample_row, ...
                        sig_up, gt_cache, [], fast_ctx, use_fast);
                end
                clear sig_up
            end
            for m = 1:n_g
                tg_idx = group_tg_idx(m);
                tg = threshold_groups(tg_idx, :);
                detail_rows = vertcat(tg_parts{m, :});
                save_tg_checkpoint(cp_dir, tg, variant_specs, cal_manifest, S60, cfg, detail_rows);
                detail_by_tg{tg_idx} = detail_rows;
            end
            clear tg_parts
        end
    end
end

% 4. 组装全量明细并一次性写入最终报表与汇总 MAT
detail_table = vertcat(detail_by_tg{:});
validate_detail_table(detail_table);

summary_table = build_screen_summary(detail_table, threshold_groups, variant_specs);
ranking_table = build_screen_rule_ranking(summary_table, variant_specs);

writetable(detail_table, fullfile(p.screen, "screen_detail.csv"));
writetable(summary_table, fullfile(p.screen, "screen_summary.csv"));
writetable(ranking_table, fullfile(p.screen, "screen_rule_ranking.csv"));
zc1bit_atomic_save(fullfile(p.screen, "screen_results.mat"), struct( ...
    "threshold_groups", threshold_groups, ...
    "variant_specs", variant_specs, ...
    "cal_manifest", cal_manifest, ...
    "roi_manifest", gt_cache.roi_manifest, ...
    "detail_table", detail_table, ...
    "summary_table", summary_table, ...
    "ranking_table", ranking_table));
write_screen_metadata(p.screen, cfg, cal_manifest, threshold_groups, variant_specs, ranking_table, cache_gb);

fprintf("阶段一规则快筛完成，输出目录：%s\n", p.screen);
out = struct( ...
    "detail_table", detail_table, ...
    "summary_table", summary_table, ...
    "ranking_table", ranking_table, ...
    "output_dir", p.screen);
end

function [cache_gb, preload_signals] = resolve_runtime_options(cfg)
cache_gb = 5;
preload_signals = true;
if isfield(cfg, "runtime") && isstruct(cfg.runtime)
    if isfield(cfg.runtime, "cache_gb") && ~isempty(cfg.runtime.cache_gb)
        cache_gb = double(cfg.runtime.cache_gb);
    end
    if isfield(cfg.runtime, "preload_signals") && ~isempty(cfg.runtime.preload_signals)
        preload_signals = logical(cfg.runtime.preload_signals);
    end
end
if ~isscalar(cache_gb) || ~isfinite(cache_gb) || cache_gb <= 0
    error("zc1bit:InvalidCacheBudget", "cfg.runtime.cache_gb 必须为正有限标量。");
end
end

function expected_params = build_tg_expected_params(tg, variant_specs, cal_manifest, S60, cfg)
expected_params = struct( ...
    "threshold", zc1bit_parameters(tg), ...
    "variants", zc1bit_parameters(variant_specs), ...
    "manifest", zc1bit_parameters(cal_manifest), ...
    "S60", zc1bit_parameters(S60), ...
    "enl_window_size", cfg.enl_window_size, ...
    "enl_stride", cfg.enl_stride, ...
    "entropy_num_bins", cfg.entropy_num_bins);
end

function save_tg_checkpoint(cp_dir, tg, variant_specs, cal_manifest, S60, cfg, detail_rows)
cp_path = fullfile(cp_dir, tg.ThresholdKey + ".mat");
expected_params = build_tg_expected_params(tg, variant_specs, cal_manifest, S60, cfg);
state = struct( ...
    "parameters", expected_params, ...
    "detail_rows", detail_rows, ...
    "completed", true);
zc1bit_atomic_save(cp_path, struct("state", state));
fprintf("  [完成] %s（已原子落盘 %d 条校准明细）\n", tg.ThresholdKey, height(detail_rows));
end

function gt_cache = try_load_cal_gt_cache(cfg, S60, cal_manifest, cache_path)
gt_cache = [];
if ~isfile(cache_path)
    return;
end
state = zc1bit_load_state(cache_path, struct());
if ~isfield(state, "completed") || ~isequal(state.completed, true)
    return;
end
expected_params = build_gt_expected_params(cfg, S60, cal_manifest);
zc1bit_check_parameters(state.parameters, expected_params, "校准集 GT 缓存");
gt_cache = state.gt_cache;
end

function gt_cache = build_and_save_cal_gt_cache(cfg, S60, cal_manifest, cal_signals, cache_path)
n_cal = height(cal_manifest);
img_gts = cell(n_cal, 1);
gt_entropy = zeros(n_cal, 1);
for idx = 1:n_cal
    img_gt = generate_gt_image(cal_signals{idx}, S60);
    img_gts{idx} = img_gt;
    gt_entropy(idx) = compute_entropy(img_gt, cfg.entropy_num_bins);
end

[is_enl, enl_top, enl_left, ~, roi_manifest] = build_enl_regions(img_gts, cal_manifest, cfg);
gt_enl = nan(n_cal, 1);
for idx = find(is_enl).'
    gt_enl(idx) = compute_enl(img_gts{idx}, enl_top(idx), enl_left(idx), cfg.enl_window_size);
end

gt_cache = struct( ...
    "SampleID", cal_manifest.SampleID, ...
    "img_gts", {img_gts}, ...
    "gt_entropy", gt_entropy, ...
    "is_enl", is_enl, ...
    "enl_top", enl_top, ...
    "enl_left", enl_left, ...
    "gt_enl", gt_enl, ...
    "roi_manifest", roi_manifest);
expected_params = build_gt_expected_params(cfg, S60, cal_manifest);
zc1bit_atomic_save(cache_path, struct("state", struct( ...
    "parameters", expected_params, ...
    "gt_cache", gt_cache, ...
    "completed", true)));
end

function expected_params = build_gt_expected_params(cfg, S60, cal_manifest)
expected_params = struct( ...
    "manifest", zc1bit_parameters(cal_manifest), ...
    "S60", zc1bit_parameters(S60), ...
    "enl_dataset_names", cfg.enl_dataset_names, ...
    "enl_window_size", cfg.enl_window_size, ...
    "enl_stride", cfg.enl_stride, ...
    "entropy_num_bins", cfg.entropy_num_bins, ...
    "raw_decimation_step", cfg.raw_decimation_step);
end

function cal_signals = preload_all_calibration_signals(cal_manifest, S60, decim_step)
% 对每个唯一回波 MAT 文件仅顺序读取一次，提取全部校准信号常驻内存。
n_cal = height(cal_manifest);
cal_signals = cell(n_cal, 1);
[unique_files, ~, file_group] = unique(cal_manifest.FilePath, "stable");
row_idx = 1 : decim_step : (1 + decim_step * (S60.nrn - 1));

fprintf("顺序加载 %d 个回波文件中的 %d 个校准样本到内存...\n", numel(unique_files), n_cal);
for f_idx = 1:numel(unique_files)
    indices = find(file_group == f_idx);
    var_name = char(cal_manifest.EchoVariable(indices(1)));
    loaded = load(unique_files(f_idx), var_name);
    raw = loaded.(var_name);
    for idx = indices(:).'
        c_start = cal_manifest.CStart(idx);
        sig60 = raw(row_idx, c_start : c_start + S60.nan - 1);
        if ~isequal(size(sig60), [S60.nrn, S60.nan])
            error("zc1bit:SignalSize", "回波块尺寸 [%d %d] 与 S60 [%d %d] 不一致。", ...
                size(sig60, 1), size(sig60, 2), S60.nrn, S60.nan);
        end
        cal_signals{idx} = sig60;
    end
    clear raw loaded
end
end

function detail_rows = run_tg_with_up_signals( ...
        cfg, S60, grid, tg, variant_specs, cal_manifest, up_signals, gt_cache, fast_ctx, use_fast)
n_cal = height(cal_manifest);
sft_carrier = [];
if ismember(tg.ThresholdType, ["RSFT", "ASFT", "2D-SFT"])
    cparams = struct("fr_over_Br", tg.FrOverBr, "fa_over_Ba", tg.FaOverBa, "phi0", tg.Phi0);
    sft_carrier = build_sft_carrier(grid, cparams, S60);
end

parts = cell(n_cal, 1);
for k = 1:n_cal
    parts{k} = evaluate_one_sample_for_tg( ...
        cfg, S60, grid, tg, variant_specs, cal_manifest(k, :), ...
        up_signals{k}, gt_cache, sft_carrier, fast_ctx, use_fast);
end
detail_rows = vertcat(parts{:});
end

function rows = evaluate_one_sample_for_tg( ...
        cfg, S60, grid, tg, variant_specs, sample_row, ...
        signal_up, gt_cache, sft_carrier, fast_ctx, use_fast)

n_var = height(variant_specs);
sid = sample_row.SampleID;
gt_idx = find(gt_cache.SampleID == sid, 1);
img_gt = gt_cache.img_gts{gt_idx};
is_enl = gt_cache.is_enl(gt_idx);
top_r = gt_cache.enl_top(gt_idx);
left_c = gt_cache.enl_left(gt_idx);
gt_enl_val = gt_cache.gt_enl(gt_idx);

ttype = tg.ThresholdType;
rq = tg.Range_q;
aq = tg.Azimuth_q;
axis_mode = tg.Axis;
ent_bins = cfg.entropy_num_bins;
enl_win = cfg.enl_window_size;

t_spec = struct("ThresholdType", ttype);
t_opts = struct("nyquist_margin", cfg.nyquist_margin, "carrier", sft_carrier);
if ttype == "RT"
    t_spec.As = tg.As;
    rt_seed = stable_rt_seed(cfg.seed, tg.ThresholdKey, sid);
    t_opts.rand_stream = RandStream("mt19937ar", "Seed", rt_seed);
elseif ismember(ttype, ["RSFT", "ASFT", "2D-SFT"])
    t_spec.STR_dB = tg.STR_dB;
    t_spec.fr_over_Br = tg.FrOverBr;
    t_spec.fa_over_Ba = tg.FaOverBa;
    t_spec.phi0 = tg.Phi0;
end

U = build_threshold(signal_up, grid, t_spec, S60, t_opts);
channel_1bit = quantize_1bit_with_u(signal_up, U);
clear U

psnr_vec = zeros(n_var, 1);
ssim_vec = zeros(n_var, 1);
ent_vec = zeros(n_var, 1);
enl_vec = nan(n_var, 1);
abs_enl_diff_vec = nan(n_var, 1);
zr_vec = zeros(n_var, 1);
mbr_vec = zeros(n_var, 1);
mbi_vec = zeros(n_var, 1);
abm_vec = zeros(n_var, 1);

v_rules = variant_specs.RuleKey;
v_modes = variant_specs.DecisionMode;
v_dcs = variant_specs.RemoveDC;
v_keys = variant_specs.VariantKey;

for v_idx = 1:n_var
    [psnr_vec(v_idx), ssim_vec(v_idx), ent_vec(v_idx), enl_vec(v_idx), ...
        zr_vec(v_idx), mbr_vec(v_idx), mbi_vec(v_idx), abm_vec(v_idx)] = ...
        evaluate_one_variant(channel_1bit, img_gt, v_rules(v_idx), axis_mode, ...
        v_modes(v_idx), v_dcs(v_idx), rq, aq, S60, fast_ctx, use_fast, ...
        ent_bins, is_enl, top_r, left_c, enl_win);
end
if is_enl
    abs_enl_diff_vec = abs(enl_vec - gt_enl_val);
end

rows = table( ...
    repmat(sid, n_var, 1), ...
    repmat(sample_row.LocalSampleIdx, n_var, 1), ...
    repmat(sample_row.Dataset, n_var, 1), ...
    repmat(tg.ThresholdKey, n_var, 1), ...
    repmat(tg.Branch, n_var, 1), ...
    repmat(axis_mode, n_var, 1), ...
    repmat(tg.GroupName, n_var, 1), ...
    repmat(tg.Q, n_var, 1), ...
    repmat(rq, n_var, 1), ...
    repmat(aq, n_var, 1), ...
    repmat(ttype, n_var, 1), ...
    repmat(tg.ParamTag, n_var, 1), ...
    repmat(tg.As, n_var, 1), ...
    repmat(tg.STR_dB, n_var, 1), ...
    repmat(tg.FrOverBr, n_var, 1), ...
    repmat(tg.FaOverBa, n_var, 1), ...
    v_rules, v_modes, v_dcs, v_keys, tg.ThresholdKey + "__" + v_keys, ...
    psnr_vec, ssim_vec, ent_vec, enl_vec, repmat(gt_enl_val, n_var, 1), abs_enl_diff_vec, ...
    zr_vec, mbr_vec, mbi_vec, abm_vec, ...
    'VariableNames', { ...
    'SampleID', 'LocalSampleIdx', 'Dataset', ...
    'ThresholdKey', 'Branch', 'Axis', 'GroupName', 'Q', 'Range_q', 'Azimuth_q', ...
    'ThresholdType', 'ParamTag', 'As', 'STR_dB', 'FrOverBr', 'FaOverBa', ...
    'RuleKey', 'DecisionMode', 'RemoveDC', 'VariantKey', 'TaskKey', ...
    'PSNR', 'SSIM', 'Entropy', 'ENL', 'GT_ENL', 'AbsENLDiff', ...
    'ZeroRatio', 'MeanBiasReal', 'MeanBiasImag', 'AbsMeanBias'});
end

function [psnr_v, ssim_v, ent_v, enl_v, zr_v, mbr_v, mbi_v, abm_v] = ...
        evaluate_one_variant(channel_1bit, img_gt, rule_key, axis_mode, ...
        decision_mode, remove_dc, rq, aq, S60, fast_ctx, use_fast, ...
        ent_bins, is_enl, top_r, left_c, enl_win)

qopts = struct("Axis", axis_mode, "DecisionMode", decision_mode, "RemoveDC", remove_dc);
[channel_p15, qdiag] = apply_zc_pseudo15bit(channel_1bit, rule_key, qopts);
if use_fast
    img_out = fast_focus_complex_roi(channel_p15, fast_ctx, rq, aq);
else
    img_out = focus_upsampled_channel(channel_p15, S60, rq, aq);
end
m = compute_image_metrics(img_out, img_gt);
psnr_v = m.PSNR;
ssim_v = m.SSIM;
ent_v = compute_entropy(img_out, ent_bins);
if is_enl
    enl_v = compute_enl(img_out, top_r, left_c, enl_win);
else
    enl_v = NaN;
end
zr_v = qdiag.ZeroRatio;
mbr_v = qdiag.MeanBiasReal;
mbi_v = qdiag.MeanBiasImag;
abm_v = qdiag.AbsMeanBias;
end

function summary_table = build_screen_summary(detail_table, threshold_groups, variant_specs)
n_tg = height(threshold_groups);
n_var = height(variant_specs);
n_rows = n_tg * n_var;

TaskKey = strings(n_rows, 1);
ThresholdKey = strings(n_rows, 1);
Branch = strings(n_rows, 1);
Axis = strings(n_rows, 1);
GroupName = strings(n_rows, 1);
Q = zeros(n_rows, 1);
Range_q = zeros(n_rows, 1);
Azimuth_q = zeros(n_rows, 1);
ThresholdType = strings(n_rows, 1);
ParamTag = strings(n_rows, 1);
RuleKey = strings(n_rows, 1);
DecisionMode = strings(n_rows, 1);
RemoveDC = false(n_rows, 1);
VariantKey = strings(n_rows, 1);
SampleCount = zeros(n_rows, 1);
ENLSampleCount = zeros(n_rows, 1);
PSNR_Mean = zeros(n_rows, 1);
PSNR_Std = zeros(n_rows, 1);
SSIM_Mean = zeros(n_rows, 1);
SSIM_Std = zeros(n_rows, 1);
Entropy_Mean = zeros(n_rows, 1);
ENL_Mean = nan(n_rows, 1);
GT_ENL_Mean = nan(n_rows, 1);
AbsENLDiff_Mean = nan(n_rows, 1);
ZeroRatio_Mean = zeros(n_rows, 1);
AbsMeanBias_Mean = zeros(n_rows, 1);
DeltaPSNR_vs_B0 = nan(n_rows, 1);
P_PSNR_vs_B0 = nan(n_rows, 1);
DeltaSSIM_vs_B0 = nan(n_rows, 1);
P_SSIM_vs_B0 = nan(n_rows, 1);
DeltaAbsENLDiff_vs_B0 = nan(n_rows, 1);
P_AbsENLDiff_vs_B0 = nan(n_rows, 1);

ptr = 0;
for tg_idx = 1:n_tg
    tg = threshold_groups(tg_idx, :);
    tg_rows = detail_table(detail_table.ThresholdKey == tg.ThresholdKey, :);
    b0_rows = tg_rows(tg_rows.VariantKey == "B0_par_raw", :);
    b0_rows = sortrows(b0_rows, "SampleID");

    for v_idx = 1:n_var
        vs = variant_specs(v_idx, :);
        ptr = ptr + 1;
        sub = tg_rows(tg_rows.VariantKey == vs.VariantKey, :);
        sub = sortrows(sub, "SampleID");
        if isempty(sub)
            error("zc1bit:MissingScreenRows", "汇总阶段缺少任务结果：%s / %s", tg.ThresholdKey, vs.VariantKey);
        end

        TaskKey(ptr) = tg.ThresholdKey + "__" + vs.VariantKey;
        ThresholdKey(ptr) = tg.ThresholdKey;
        Branch(ptr) = tg.Branch;
        Axis(ptr) = tg.Axis;
        GroupName(ptr) = tg.GroupName;
        Q(ptr) = tg.Q;
        Range_q(ptr) = tg.Range_q;
        Azimuth_q(ptr) = tg.Azimuth_q;
        ThresholdType(ptr) = tg.ThresholdType;
        ParamTag(ptr) = tg.ParamTag;
        RuleKey(ptr) = vs.RuleKey;
        DecisionMode(ptr) = vs.DecisionMode;
        RemoveDC(ptr) = vs.RemoveDC;
        VariantKey(ptr) = vs.VariantKey;
        SampleCount(ptr) = height(sub);
        enl_mask = ~isnan(sub.ENL);
        ENLSampleCount(ptr) = sum(enl_mask);

        PSNR_Mean(ptr) = mean(sub.PSNR);
        PSNR_Std(ptr) = std(sub.PSNR);
        SSIM_Mean(ptr) = mean(sub.SSIM);
        SSIM_Std(ptr) = std(sub.SSIM);
        Entropy_Mean(ptr) = mean(sub.Entropy);
        if any(enl_mask)
            ENL_Mean(ptr) = mean(sub.ENL(enl_mask));
            GT_ENL_Mean(ptr) = mean(sub.GT_ENL(enl_mask));
            AbsENLDiff_Mean(ptr) = mean(sub.AbsENLDiff(enl_mask));
        end
        ZeroRatio_Mean(ptr) = mean(sub.ZeroRatio);
        AbsMeanBias_Mean(ptr) = mean(sub.AbsMeanBias);

        if ~isempty(b0_rows) && isequal(sub.SampleID, b0_rows.SampleID)
            wp = compute_paired_wilcoxon(sub.PSNR, b0_rows.PSNR);
            ws = compute_paired_wilcoxon(sub.SSIM, b0_rows.SSIM);
            DeltaPSNR_vs_B0(ptr) = wp.MeanDiff;
            P_PSNR_vs_B0(ptr) = wp.PValue;
            DeltaSSIM_vs_B0(ptr) = ws.MeanDiff;
            P_SSIM_vs_B0(ptr) = ws.PValue;
            b0_enl_mask = ~isnan(b0_rows.AbsENLDiff);
            if any(enl_mask & b0_enl_mask)
                we = compute_paired_wilcoxon( ...
                    sub.AbsENLDiff(enl_mask & b0_enl_mask), ...
                    b0_rows.AbsENLDiff(enl_mask & b0_enl_mask));
                DeltaAbsENLDiff_vs_B0(ptr) = we.MeanDiff;
                P_AbsENLDiff_vs_B0(ptr) = we.PValue;
            end
        end
    end
end

summary_table = table( ...
    TaskKey, ThresholdKey, Branch, Axis, GroupName, Q, Range_q, Azimuth_q, ...
    ThresholdType, ParamTag, RuleKey, DecisionMode, RemoveDC, VariantKey, ...
    SampleCount, ENLSampleCount, ...
    PSNR_Mean, PSNR_Std, SSIM_Mean, SSIM_Std, Entropy_Mean, ...
    ENL_Mean, GT_ENL_Mean, AbsENLDiff_Mean, ZeroRatio_Mean, AbsMeanBias_Mean, ...
    DeltaPSNR_vs_B0, P_PSNR_vs_B0, DeltaSSIM_vs_B0, P_SSIM_vs_B0, ...
    DeltaAbsENLDiff_vs_B0, P_AbsENLDiff_vs_B0);
end

function ranking_table = build_screen_rule_ranking(summary_table, variant_specs)
n_var = height(variant_specs);
RuleKey = variant_specs.RuleKey;
DecisionMode = variant_specs.DecisionMode;
RemoveDC = variant_specs.RemoveDC;
VariantKey = variant_specs.VariantKey;

ConfigCount = zeros(n_var, 1);
MeanPSNR = zeros(n_var, 1);
MeanSSIM = zeros(n_var, 1);
MeanAbsENLDiff = nan(n_var, 1);
MeanZeroRatio = zeros(n_var, 1);
MeanAbsBias = zeros(n_var, 1);
MeanDeltaPSNR = nan(n_var, 1);
MeanDeltaSSIM = nan(n_var, 1);
MeanDeltaAbsENLDiff = nan(n_var, 1);
SSIMWinCount = zeros(n_var, 1);
SSIMSigWinCount = zeros(n_var, 1);

for idx = 1:n_var
    sub = summary_table(summary_table.VariantKey == VariantKey(idx), :);
    ConfigCount(idx) = height(sub);
    MeanPSNR(idx) = mean(sub.PSNR_Mean);
    MeanSSIM(idx) = mean(sub.SSIM_Mean);
    if any(~isnan(sub.AbsENLDiff_Mean))
        MeanAbsENLDiff(idx) = mean(sub.AbsENLDiff_Mean, "omitnan");
    end
    MeanZeroRatio(idx) = mean(sub.ZeroRatio_Mean);
    MeanAbsBias(idx) = mean(sub.AbsMeanBias_Mean);
    MeanDeltaPSNR(idx) = mean(sub.DeltaPSNR_vs_B0, "omitnan");
    MeanDeltaSSIM(idx) = mean(sub.DeltaSSIM_vs_B0, "omitnan");
    if any(~isnan(sub.DeltaAbsENLDiff_vs_B0))
        MeanDeltaAbsENLDiff(idx) = mean(sub.DeltaAbsENLDiff_vs_B0, "omitnan");
    end
    SSIMWinCount(idx) = sum(sub.DeltaSSIM_vs_B0 > 0);
    SSIMSigWinCount(idx) = sum(sub.DeltaSSIM_vs_B0 > 0 & sub.P_SSIM_vs_B0 < 0.05);
end

ranking_table = table( ...
    RuleKey, DecisionMode, RemoveDC, VariantKey, ConfigCount, ...
    MeanSSIM, MeanPSNR, MeanAbsENLDiff, MeanZeroRatio, MeanAbsBias, ...
    MeanDeltaSSIM, MeanDeltaPSNR, MeanDeltaAbsENLDiff, SSIMWinCount, SSIMSigWinCount);
ranking_table = sortrows(ranking_table, ...
    {'MeanSSIM', 'MeanPSNR', 'MeanAbsENLDiff', 'VariantKey'}, ...
    {'descend', 'descend', 'ascend', 'ascend'});
ranking_table.Rank = (1:height(ranking_table)).';
ranking_table = movevars(ranking_table, "Rank", "Before", "RuleKey");
end

function validate_detail_table(detail_table)
if isempty(detail_table) || height(detail_table) == 0
    error("zc1bit:EmptyDetailTable", "阶段一明细表为空。");
end
for col = ["PSNR", "SSIM", "Entropy", "ZeroRatio", "MeanBiasReal", "MeanBiasImag", "AbsMeanBias"]
    vals = detail_table.(col);
    if any(isnan(vals)) || any(isinf(vals) & vals < 0)
        error("zc1bit:InvalidDetailValue", "明细列 %s 包含非法值。", col);
    end
end
if any(~isfinite(detail_table.SSIM)) || any(~isfinite(detail_table.Entropy))
    error("zc1bit:NonFiniteDetail", "明细表 SSIM 或 Entropy 存在非有限值。");
end
end

function seed = stable_rt_seed(base_seed, threshold_key, sample_id)
chars = double(char(threshold_key));
h = 0;
for k = 1:numel(chars)
    h = mod(h * 131 + chars(k), 1000003);
end
seed = mod(double(base_seed) + h * 1000 + double(sample_id), 2^31 - 1);
end

function write_screen_metadata(screen_dir, cfg, cal_manifest, threshold_groups, variant_specs, ranking_table, cache_gb)
fid = fopen(fullfile(screen_dir, "screen_metadata.txt"), "w", "n", "UTF-8");
if fid < 0
    error("zc1bit:MetadataWrite", "无法写入 screen_metadata.txt。");
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "Experiment=SAR-ZC1Bit Stage 1 Screen\n");
fprintf(fid, "ExecutionMode=single_core\n");
fprintf(fid, "CacheBudgetGB=%.6g\n", cache_gb);
fprintf(fid, "Seed=%d\n", cfg.seed);
fprintf(fid, "CalibrationSampleCount=%d\n", height(cal_manifest));
fprintf(fid, "ThresholdConfigCount=%d\n", height(threshold_groups));
fprintf(fid, "VariantCount=%d\n", height(variant_specs));
fprintf(fid, "Rules=%s\n", join(string(cfg.screen.rules), ","));
fprintf(fid, "DecisionModes=%s\n", join(string(cfg.screen.decision_modes), ","));
fprintf(fid, "Top1Variant=%s (MeanSSIM=%.6f, MeanPSNR=%.4f)\n", ...
    ranking_table.VariantKey(1), ranking_table.MeanSSIM(1), ranking_table.MeanPSNR(1));
end
