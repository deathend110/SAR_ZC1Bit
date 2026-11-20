function out = zc1bit_screen_impl(cfg)
%ZC1BIT_SCREEN_IMPL 阶段一（置零规则快筛）核心工作流实现。
% 严格仅在 35 个校准样本（奇数 LocalSampleIdx）上执行，绝不触碰测试集。
% 按“阈值配置 × 场景”原子保存 checkpoint，中断恢复仅核对实际计算参数。

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

% 1. 准备或恢复校准集 GT 与固定 64x64 ENL ROI（不读取测试集）
gt_cache = load_or_build_cal_gt_cache(cfg, S60, cal_manifest, p);
if ~isfolder(p.screen)
    mkdir(p.screen);
end
writetable(gt_cache.roi_manifest, fullfile(p.screen, "screen_cal_enl_roi_manifest.csv"));

% 2. 构建频域直裁快路径上下文（若启用）
use_fast = true;
if isfield(cfg.screen, "use_fast_path")
    use_fast = logical(cfg.screen.use_fast_path);
end
fast_ctx = [];
if use_fast
    fast_ctx = build_fast_focus_context(S60, threshold_groups.Range_q, threshold_groups.Azimuth_q);
end

% 3. 解析并行 worker（失败时自动降级串行，数值结果完全一致）
workers = resolve_parallel_workers(cfg.screen.workers);

% 4. 按“阈值配置 × 场景”执行并原子落盘
scenes = unique(cal_manifest.Dataset, "stable");
n_tg = height(threshold_groups);
n_scenes = numel(scenes);
n_var = height(variant_specs);

all_detail_parts = cell(n_tg * n_scenes, 1);
part_ptr = 0;

fprintf("\n=== 阶段一规则快筛：%d 阈值配置 × %d 场景 × %d 变体（仅 %d 校准样本）===\n", ...
    n_tg, n_scenes, n_var, height(cal_manifest));

scene_echo_cache = struct("scene", "", "signals", {{}}, "sample_ids", []);

for tg_idx = 1:n_tg
    tg = threshold_groups(tg_idx, :);
    grid = resolve_sampling_grid(tg.Range_q, tg.Azimuth_q, [S60.nrn, S60.nan], S60);
    sft_carrier = [];
    if ismember(tg.ThresholdType, ["RSFT", "ASFT", "2D-SFT"])
        cparams = struct("fr_over_Br", tg.FrOverBr, "fa_over_Ba", tg.FaOverBa, "phi0", tg.Phi0);
        sft_carrier = build_sft_carrier(grid, cparams, S60);
    end

    for s_idx = 1:n_scenes
        scene_name = scenes(s_idx);
        scene_key = make_scene_key(scene_name);
        scene_mask = (cal_manifest.Dataset == scene_name);
        scene_manifest = cal_manifest(scene_mask, :);

        cp_path = fullfile(p.checkpoints, "screen", tg.ThresholdKey, scene_key + ".mat");
        expected_params = struct( ...
            "threshold", zc1bit_parameters(tg), ...
            "variants", zc1bit_parameters(variant_specs), ...
            "manifest", zc1bit_parameters(scene_manifest), ...
            "S60", zc1bit_parameters(S60), ...
            "enl_window_size", cfg.enl_window_size, ...
            "enl_stride", cfg.enl_stride, ...
            "entropy_num_bins", cfg.entropy_num_bins);

        part_ptr = part_ptr + 1;
        if isfile(cp_path)
            state = zc1bit_load_state(cp_path, struct());
            if isfield(state, "completed") && isequal(state.completed, true) && isfield(state, "detail_rows")
                zc1bit_check_parameters(state.parameters, expected_params, ...
                    sprintf("阶段一检查点(%s/%s)", tg.ThresholdKey, scene_key));
                all_detail_parts{part_ptr} = state.detail_rows;
                continue;
            end
        end

        % 按需加载当前场景的 5 个校准回波（避免断点全跳过时读取磁盘）
        if scene_echo_cache.scene ~= scene_name
            [scene_echo_cache.sample_ids, scene_echo_cache.signals] = ...
                load_scene_calibration_signals(scene_manifest, S60, cfg.raw_decimation_step);
            scene_echo_cache.scene = scene_name;
        end

        fprintf("  [%02d/%02d] %s | 场景 %s (%d 样本)\n", ...
            tg_idx, n_tg, tg.ThresholdKey, scene_name, height(scene_manifest));

        detail_rows = evaluate_scene_threshold_group( ...
            cfg, S60, grid, tg, variant_specs, scene_manifest, ...
            scene_echo_cache.signals, gt_cache, sft_carrier, fast_ctx, use_fast, workers);

        state = struct( ...
            "parameters", expected_params, ...
            "detail_rows", detail_rows, ...
            "completed", true);
        zc1bit_atomic_save(cp_path, struct("state", state));
        all_detail_parts{part_ptr} = detail_rows;
    end
end

% 5. 组装明细表并生成汇总、规则排名与配对统计
detail_table = vertcat(all_detail_parts{:});
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
write_screen_metadata(p.screen, cfg, cal_manifest, threshold_groups, variant_specs, ranking_table);

fprintf("阶段一规则快筛完成，输出目录：%s\n", p.screen);
out = struct( ...
    "detail_table", detail_table, ...
    "summary_table", summary_table, ...
    "ranking_table", ranking_table, ...
    "output_dir", p.screen);
end

function gt_cache = load_or_build_cal_gt_cache(cfg, S60, cal_manifest, p)
cache_path = fullfile(p.checkpoints, "screen", "gt_cal_cache.mat");
expected_params = struct( ...
    "manifest", zc1bit_parameters(cal_manifest), ...
    "S60", zc1bit_parameters(S60), ...
    "enl_dataset_names", cfg.enl_dataset_names, ...
    "enl_window_size", cfg.enl_window_size, ...
    "enl_stride", cfg.enl_stride, ...
    "entropy_num_bins", cfg.entropy_num_bins, ...
    "raw_decimation_step", cfg.raw_decimation_step);

if isfile(cache_path)
    state = zc1bit_load_state(cache_path, struct());
    if isfield(state, "completed") && isequal(state.completed, true)
        zc1bit_check_parameters(state.parameters, expected_params, "校准集 GT 缓存");
        gt_cache = state.gt_cache;
        return;
    end
end

n_cal = height(cal_manifest);
img_gts = cell(n_cal, 1);
gt_entropy = zeros(n_cal, 1);
scenes = unique(cal_manifest.Dataset, "stable");

for s_idx = 1:numel(scenes)
    scene_mask = (cal_manifest.Dataset == scenes(s_idx));
    scene_rows = cal_manifest(scene_mask, :);
    [~, signals] = load_scene_calibration_signals(scene_rows, S60, cfg.raw_decimation_step);
    indices = find(scene_mask);
    for k = 1:numel(indices)
        idx = indices(k);
        img_gt = generate_gt_image(signals{k}, S60);
        img_gts{idx} = img_gt;
        gt_entropy(idx) = compute_entropy(img_gt, cfg.entropy_num_bins);
    end
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
zc1bit_atomic_save(cache_path, struct("state", struct( ...
    "parameters", expected_params, ...
    "gt_cache", gt_cache, ...
    "completed", true)));
end

function [sample_ids, signals] = load_scene_calibration_signals(scene_manifest, S60, decim_step)
n = height(scene_manifest);
sample_ids = scene_manifest.SampleID;
signals = cell(n, 1);
loaded = load(scene_manifest.FilePath(1), char(scene_manifest.EchoVariable(1)));
raw = loaded.(char(scene_manifest.EchoVariable(1)));
row_idx = 1 : decim_step : (1 + decim_step * (S60.nrn - 1));
for k = 1:n
    c_start = scene_manifest.CStart(k);
    sig60 = raw(row_idx, c_start : c_start + S60.nan - 1);
    if ~isequal(size(sig60), [S60.nrn, S60.nan])
        error("zc1bit:SignalSize", "回波块尺寸 [%d %d] 与 S60 [%d %d] 不一致。", ...
            size(sig60, 1), size(sig60, 2), S60.nrn, S60.nan);
    end
    signals{k} = sig60;
end
end

function rows = evaluate_scene_threshold_group( ...
        cfg, S60, grid, tg, variant_specs, scene_manifest, ...
        signals, gt_cache, sft_carrier, fast_ctx, use_fast, workers)

n_samples = height(scene_manifest);
n_var = height(variant_specs);
n_total = n_samples * n_var;

psnr_vec = zeros(n_total, 1);
ssim_vec = zeros(n_total, 1);
ent_vec = zeros(n_total, 1);
enl_vec = nan(n_total, 1);
gt_enl_vec = nan(n_total, 1);
abs_enl_diff_vec = nan(n_total, 1);
zr_vec = zeros(n_total, 1);
mbr_vec = zeros(n_total, 1);
mbi_vec = zeros(n_total, 1);
abm_vec = zeros(n_total, 1);

sample_id_vec = zeros(n_total, 1);
local_idx_vec = zeros(n_total, 1);
dataset_vec = strings(n_total, 1);
rule_vec = strings(n_total, 1);
mode_vec = strings(n_total, 1);
dc_vec = false(n_total, 1);
vkey_vec = strings(n_total, 1);
tkey_vec = strings(n_total, 1);

v_rules = variant_specs.RuleKey;
v_modes = variant_specs.DecisionMode;
v_dcs = variant_specs.RemoveDC;
v_keys = variant_specs.VariantKey;
axis_mode = tg.Axis;
rq = tg.Range_q;
aq = tg.Azimuth_q;
ttype = tg.ThresholdType;
ent_bins = cfg.entropy_num_bins;
enl_win = cfg.enl_window_size;

for k = 1:n_samples
    sid = scene_manifest.SampleID(k);
    gt_idx = find(gt_cache.SampleID == sid, 1);
    img_gt = gt_cache.img_gts{gt_idx};
    is_enl = gt_cache.is_enl(gt_idx);
    top_r = gt_cache.enl_top(gt_idx);
    left_c = gt_cache.enl_left(gt_idx);
    gt_enl_val = gt_cache.gt_enl(gt_idx);

    signal_up = two_dim_upsample_fft(signals{k}, rq, aq);

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
    clear signal_up U

    p_arr = zeros(n_var, 1);
    s_arr = zeros(n_var, 1);
    e_arr = zeros(n_var, 1);
    enl_arr = nan(n_var, 1);
    zr_arr = zeros(n_var, 1);
    mbr_arr = zeros(n_var, 1);
    mbi_arr = zeros(n_var, 1);
    abm_arr = zeros(n_var, 1);

    if workers > 0
        parfor (v_idx = 1:n_var, workers)
            [p_arr(v_idx), s_arr(v_idx), e_arr(v_idx), enl_arr(v_idx), ...
                zr_arr(v_idx), mbr_arr(v_idx), mbi_arr(v_idx), abm_arr(v_idx)] = ...
                evaluate_one_variant(channel_1bit, img_gt, v_rules(v_idx), axis_mode, ...
                v_modes(v_idx), v_dcs(v_idx), rq, aq, S60, fast_ctx, use_fast, ...
                ent_bins, is_enl, top_r, left_c, enl_win);
        end
    else
        for v_idx = 1:n_var
            [p_arr(v_idx), s_arr(v_idx), e_arr(v_idx), enl_arr(v_idx), ...
                zr_arr(v_idx), mbr_arr(v_idx), mbi_arr(v_idx), abm_arr(v_idx)] = ...
                evaluate_one_variant(channel_1bit, img_gt, v_rules(v_idx), axis_mode, ...
                v_modes(v_idx), v_dcs(v_idx), rq, aq, S60, fast_ctx, use_fast, ...
                ent_bins, is_enl, top_r, left_c, enl_win);
        end
    end

    dest = (k - 1) * n_var + (1:n_var);
    psnr_vec(dest) = p_arr;
    ssim_vec(dest) = s_arr;
    ent_vec(dest) = e_arr;
    enl_vec(dest) = enl_arr;
    gt_enl_vec(dest) = repmat(gt_enl_val, n_var, 1);
    if is_enl
        abs_enl_diff_vec(dest) = abs(enl_arr - gt_enl_val);
    end
    zr_vec(dest) = zr_arr;
    mbr_vec(dest) = mbr_arr;
    mbi_vec(dest) = mbi_arr;
    abm_vec(dest) = abm_arr;

    sample_id_vec(dest) = repmat(sid, n_var, 1);
    local_idx_vec(dest) = repmat(scene_manifest.LocalSampleIdx(k), n_var, 1);
    dataset_vec(dest) = repmat(scene_manifest.Dataset(k), n_var, 1);
    rule_vec(dest) = v_rules;
    mode_vec(dest) = v_modes;
    dc_vec(dest) = v_dcs;
    vkey_vec(dest) = v_keys;
    tkey_vec(dest) = tg.ThresholdKey + "__" + v_keys;
end

rows = table( ...
    sample_id_vec, local_idx_vec, dataset_vec, ...
    repmat(tg.ThresholdKey, n_total, 1), ...
    repmat(tg.Branch, n_total, 1), ...
    repmat(tg.Axis, n_total, 1), ...
    repmat(tg.GroupName, n_total, 1), ...
    repmat(tg.Q, n_total, 1), ...
    repmat(tg.Range_q, n_total, 1), ...
    repmat(tg.Azimuth_q, n_total, 1), ...
    repmat(tg.ThresholdType, n_total, 1), ...
    repmat(tg.ParamTag, n_total, 1), ...
    repmat(tg.As, n_total, 1), ...
    repmat(tg.STR_dB, n_total, 1), ...
    repmat(tg.FrOverBr, n_total, 1), ...
    repmat(tg.FaOverBa, n_total, 1), ...
    rule_vec, mode_vec, dc_vec, vkey_vec, tkey_vec, ...
    psnr_vec, ssim_vec, ent_vec, enl_vec, gt_enl_vec, abs_enl_diff_vec, ...
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

function workers = resolve_parallel_workers(requested_workers)
workers = 0;
if requested_workers <= 0
    return;
end
pool = gcp("nocreate");
if isempty(pool)
    try
        max_cores = max(1, feature("numcores"));
        pool = parpool("threads", min(requested_workers, max_cores));
        workers = pool.NumWorkers;
    catch err
        warning("zc1bit:ParallelFallback", ...
            "并行池启动失败，已自动降级为串行执行 (workers=0)：%s", err.message);
        workers = 0;
    end
else
    workers = min(requested_workers, pool.NumWorkers);
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

function key = make_scene_key(scene_name)
key = lower(regexprep(string(scene_name), "[^a-zA-Z0-9]+", "_"));
key = regexprep(key, "^_+|_+$", "");
end

function write_screen_metadata(screen_dir, cfg, cal_manifest, threshold_groups, variant_specs, ranking_table)
fid = fopen(fullfile(screen_dir, "screen_metadata.txt"), "w", "n", "UTF-8");
if fid < 0
    error("zc1bit:MetadataWrite", "无法写入 screen_metadata.txt。");
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "Experiment=SAR-ZC1Bit Stage 1 Screen\n");
fprintf(fid, "Seed=%d\n", cfg.seed);
fprintf(fid, "CalibrationSampleCount=%d\n", height(cal_manifest));
fprintf(fid, "ThresholdConfigCount=%d\n", height(threshold_groups));
fprintf(fid, "VariantCount=%d\n", height(variant_specs));
fprintf(fid, "Rules=%s\n", join(string(cfg.screen.rules), ","));
fprintf(fid, "DecisionModes=%s\n", join(string(cfg.screen.decision_modes), ","));
fprintf(fid, "Top1Variant=%s (MeanSSIM=%.6f, MeanPSNR=%.4f)\n", ...
    ranking_table.VariantKey(1), ranking_table.MeanSSIM(1), ranking_table.MeanPSNR(1));
end
