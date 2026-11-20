function out = zc1bit_search_impl(cfg)
%ZC1BIT_SEARCH_IMPL 阶段二（35 校准样本独立阈值参数搜索）单核内存优先工作流实现。
% 仅覆盖单向距离 R1A1..R9A1 × {ZT,RT,RSFT} 与单向方位 R1A1..R1A9 × {ZT,RT,ASFT}，
% 对基线 B0 与阶段一自动选出的 Top 5 非 B0 变体各自独立寻优；严禁使用任何测试样本。

if nargin == 0
    cfg = zc1bit_config();
end
directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>

record = zc1bit_prepare_run(cfg);
p = zc1bit_paths(cfg);
S60 = load(cfg.parameter_file_60);

selected_variants = zc1bit_resolve_stage2_variants(cfg, p);
units = zc1bit_build_stage2_units(cfg, S60);

manifest = record.manifest;
cal_manifest = manifest(manifest.Split == "calibration", :);
if isempty(cal_manifest)
    error("zc1bit:EmptyCalibration", "样本清单中未找到校准样本。");
end

[cache_gb, preload_signals] = resolve_runtime_options(cfg);
cache_budget_bytes = cache_gb * (1024 ^ 3);
use_fast = true;
if isfield(cfg.stage2, "use_fast_path")
    use_fast = logical(cfg.stage2.use_fast_path);
end

n_units = height(units);
n_cal = height(cal_manifest);
n_var = height(selected_variants);

cp_dir = fullfile(p.checkpoints, "search");
summary_by_unit = cell(n_units, 1);
best_by_unit = cell(n_units, 1);
completed_units = false(n_units, 1);

for u_idx = 1:n_units
    unit = units(u_idx, :);
    cp_path = fullfile(cp_dir, unit.UnitKey + ".mat");
    if ~isfile(cp_path)
        continue;
    end
    state = zc1bit_load_state(cp_path, struct());
    if isfield(state, "completed") && isequal(state.completed, true) && ...
            isfield(state, "summary_rows") && isfield(state, "best_rows")
        expected_params = build_unit_expected_params(unit, selected_variants, cal_manifest, S60, cfg);
        zc1bit_check_parameters(state.parameters, expected_params, ...
            sprintf("阶段二搜索检查点(%s)", unit.UnitKey));
        summary_by_unit{u_idx} = state.summary_rows;
        best_by_unit{u_idx} = state.best_rows;
        completed_units(u_idx) = true;
    end
end

gt_cache_path = fullfile(cp_dir, "gt_cal_cache.mat");
screen_gt_path = fullfile(p.checkpoints, "screen", "gt_cal_cache.mat");
gt_cache = zc1bit_stage2_gt_cache(cfg, S60, cal_manifest, {}, gt_cache_path, screen_gt_path);
cal_signals = {};

if any(~completed_units) || isempty(gt_cache)
    cal_signals = zc1bit_stage2_load_signals(cal_manifest, S60, cfg.raw_decimation_step, "校准样本");
    if isempty(gt_cache)
        gt_cache = zc1bit_stage2_gt_cache(cfg, S60, cal_manifest, cal_signals, gt_cache_path, screen_gt_path);
    end
end

if ~isfolder(p.search)
    mkdir(p.search);
end

fprintf("\n=== 阶段二参数搜索：%d 单元 × %d 变体（%d 校准样本 | 已命中 %d/%d 检查点 | cache=%.1f GB）===\n", ...
    n_units, n_var, n_cal, sum(completed_units), n_units, cache_gb);

pending_indices = find(~completed_units);
if ~isempty(pending_indices)
    fast_ctx = [];
    if use_fast
        fast_ctx = build_fast_focus_context( ...
            S60, units.Range_q(pending_indices), units.Azimuth_q(pending_indices));
    end

    base_info = whos("cal_signals");
    base_bytes = base_info.bytes;
    if isa(cal_signals{1}, "single")
        elem_bytes = 8;
    else
        elem_bytes = 16;
    end

    pending_allocs = unique([units.Range_q(pending_indices), units.Azimuth_q(pending_indices)], "rows", "stable");
    for a_idx = 1:size(pending_allocs, 1)
        rq = pending_allocs(a_idx, 1);
        aq = pending_allocs(a_idx, 2);
        group_u_idx = pending_indices(units.Range_q(pending_indices) == rq & units.Azimuth_q(pending_indices) == aq);
        sgrid = resolve_sampling_grid(rq, aq, [S60.nrn, S60.nan], S60);

        if rq == 1 && aq == 1
            for u_idx = group_u_idx(:).'
                [srows, brows] = run_unit_with_cached_up( ...
                    cfg, S60, sgrid, units(u_idx, :), selected_variants, ...
                    cal_manifest, cal_signals, gt_cache, fast_ctx, use_fast);
                save_unit_checkpoint(cp_dir, units(u_idx, :), selected_variants, cal_manifest, S60, cfg, srows, brows);
                summary_by_unit{u_idx} = srows;
                best_by_unit{u_idx} = brows;
            end
            continue;
        end

        est_up_bytes = n_cal * prod(sgrid.upsampled_size) * elem_bytes;
        if preload_signals && (base_bytes + est_up_bytes <= cache_budget_bytes)
            fprintf("  [R%dA%d] 内存缓存模式：预构建 %d 个上采样回波（约 %.2f GB）...\n", ...
                rq, aq, n_cal, est_up_bytes / (1024 ^ 3));
            up_signals = cell(n_cal, 1);
            for k = 1:n_cal
                up_signals{k} = two_dim_upsample_fft(cal_signals{k}, rq, aq);
            end
            for u_idx = group_u_idx(:).'
                [srows, brows] = run_unit_with_cached_up( ...
                    cfg, S60, sgrid, units(u_idx, :), selected_variants, ...
                    cal_manifest, up_signals, gt_cache, fast_ctx, use_fast);
                save_unit_checkpoint(cp_dir, units(u_idx, :), selected_variants, cal_manifest, S60, cfg, srows, brows);
                summary_by_unit{u_idx} = srows;
                best_by_unit{u_idx} = brows;
            end
            clear up_signals
        else
            fprintf("  [R%dA%d] 流式低内存模式：预计整组上采样 %.2f GB 超出预算，逐样本复用（%d 单元）...\n", ...
                rq, aq, est_up_bytes / (1024 ^ 3), numel(group_u_idx));
            [s_parts, b_parts] = run_alloc_streaming( ...
                cfg, S60, sgrid, units(group_u_idx, :), selected_variants, ...
                cal_manifest, cal_signals, gt_cache, fast_ctx, use_fast);
            for m = 1:numel(group_u_idx)
                u_idx = group_u_idx(m);
                save_unit_checkpoint(cp_dir, units(u_idx, :), selected_variants, cal_manifest, S60, cfg, s_parts{m}, b_parts{m});
                summary_by_unit{u_idx} = s_parts{m};
                best_by_unit{u_idx} = b_parts{m};
            end
        end
    end
end

summary_table = vertcat(summary_by_unit{:});
best_params_table = vertcat(best_by_unit{:});
validate_search_outputs(summary_table, best_params_table, units, selected_variants);

writetable(selected_variants, fullfile(p.search, "search_selected_variants.csv"));
writetable(summary_table, fullfile(p.search, "search_summary.csv"));
writetable(best_params_table, fullfile(p.search, "search_best_params.csv"));

signature = struct( ...
    "stage2", zc1bit_parameters(cfg.stage2), ...
    "nyquist_margin", cfg.nyquist_margin, ...
    "seed", cfg.seed, ...
    "S60", zc1bit_parameters(S60), ...
    "cal_manifest", zc1bit_parameters(cal_manifest));

zc1bit_atomic_save(fullfile(p.search, "search_results.mat"), struct( ...
    "signature", signature, ...
    "units", units, ...
    "selected_variants", selected_variants, ...
    "cal_manifest", cal_manifest, ...
    "summary_table", summary_table, ...
    "best_params_table", best_params_table));
write_search_metadata(p.search, cfg, cal_manifest, units, selected_variants, best_params_table, cache_gb);

fprintf("阶段二参数搜索完成，输出目录：%s\n", p.search);
out = struct( ...
    "selected_variants", selected_variants, ...
    "units", units, ...
    "cal_manifest", cal_manifest, ...
    "summary_table", summary_table, ...
    "best_params_table", best_params_table, ...
    "output_dir", p.search);
end

function [summary_rows, best_rows] = run_unit_with_cached_up( ...
    cfg, S60, sgrid, unit, selected_variants, cal_manifest, up_signals, gt_cache, fast_ctx, use_fast)

[coarse_cands, cand_meta] = zc1bit_stage2_candidates(unit, cfg, S60, "coarse");
coarse_details = evaluate_cands_over_cached_up( ...
    cfg, S60, sgrid, unit, coarse_cands, selected_variants, cal_manifest, up_signals, gt_cache, fast_ctx, use_fast);
coarse_summary = aggregate_unit_summary(coarse_details, unit, coarse_cands, selected_variants);

summary_rows = coarse_summary;
if ismember(unit.ThresholdType, ["RSFT", "ASFT"]) && logical(cfg.stage2.sft_refine)
    coarse_best = select_best_per_variant(coarse_summary, selected_variants, cand_meta);
    refine_cands = zc1bit_stage2_candidates(unit, cfg, S60, "refine", coarse_best, coarse_cands.ParamTag);
    if ~isempty(refine_cands)
        refine_details = evaluate_cands_over_cached_up( ...
            cfg, S60, sgrid, unit, refine_cands, selected_variants, cal_manifest, up_signals, gt_cache, fast_ctx, use_fast);
        refine_summary = aggregate_unit_summary(refine_details, unit, refine_cands, selected_variants);
        summary_rows = [coarse_summary; refine_summary];
    end
end

best_rows = select_best_per_variant(summary_rows, selected_variants, cand_meta);
end

function [s_parts, b_parts] = run_alloc_streaming( ...
    cfg, S60, sgrid, group_units, selected_variants, cal_manifest, cal_signals, gt_cache, fast_ctx, use_fast)

n_g = height(group_units);
n_cal = height(cal_manifest);
rq = double(group_units.Range_q(1));
aq = double(group_units.Azimuth_q(1));

coarse_cands = cell(n_g, 1);
cand_metas = cell(n_g, 1);
coarse_carriers = cell(n_g, 1);
for m = 1:n_g
    [coarse_cands{m}, cand_metas{m}] = zc1bit_stage2_candidates(group_units(m, :), cfg, S60, "coarse");
    coarse_carriers{m} = build_unit_carrier_map(sgrid, S60, group_units(m, :), coarse_cands{m});
end

coarse_details_by_u = cell(n_g, n_cal);
for k = 1:n_cal
    sig_up = two_dim_upsample_fft(cal_signals{k}, rq, aq);
    sample_row = cal_manifest(k, :);
    for m = 1:n_g
        coarse_details_by_u{m, k} = evaluate_cands_one_sample( ...
            cfg, S60, sgrid, group_units(m, :), coarse_cands{m}, coarse_carriers{m}, ...
            selected_variants, sample_row, sig_up, gt_cache, fast_ctx, use_fast);
    end
    clear sig_up
end

summary_by_u = cell(n_g, 1);
refine_cands = cell(n_g, 1);
refine_carriers = cell(n_g, 1);
need_refine = false(n_g, 1);
for m = 1:n_g
    unit = group_units(m, :);
    csum = aggregate_unit_summary(vertcat(coarse_details_by_u{m, :}), unit, coarse_cands{m}, selected_variants);
    summary_by_u{m} = csum;
    if ismember(unit.ThresholdType, ["RSFT", "ASFT"]) && logical(cfg.stage2.sft_refine)
        cbest = select_best_per_variant(csum, selected_variants, cand_metas{m});
        rc = zc1bit_stage2_candidates(unit, cfg, S60, "refine", cbest, coarse_cands{m}.ParamTag);
        if ~isempty(rc)
            refine_cands{m} = rc;
            refine_carriers{m} = build_unit_carrier_map(sgrid, S60, unit, rc);
            need_refine(m) = true;
        end
    end
end

if any(need_refine)
    refine_details_by_u = cell(n_g, n_cal);
    for k = 1:n_cal
        sig_up = two_dim_upsample_fft(cal_signals{k}, rq, aq);
        sample_row = cal_manifest(k, :);
        for m = find(need_refine).'
            refine_details_by_u{m, k} = evaluate_cands_one_sample( ...
                cfg, S60, sgrid, group_units(m, :), refine_cands{m}, refine_carriers{m}, ...
                selected_variants, sample_row, sig_up, gt_cache, fast_ctx, use_fast);
        end
        clear sig_up
    end
    for m = find(need_refine).'
        rsum = aggregate_unit_summary( ...
            vertcat(refine_details_by_u{m, :}), group_units(m, :), refine_cands{m}, selected_variants);
        summary_by_u{m} = [summary_by_u{m}; rsum];
    end
end

s_parts = summary_by_u;
b_parts = cell(n_g, 1);
for m = 1:n_g
    b_parts{m} = select_best_per_variant(s_parts{m}, selected_variants, cand_metas{m});
end
end

function details = evaluate_cands_over_cached_up( ...
    cfg, S60, sgrid, unit, cands, selected_variants, cal_manifest, up_signals, gt_cache, fast_ctx, use_fast)

carriers = build_unit_carrier_map(sgrid, S60, unit, cands);
n_cal = height(cal_manifest);
parts = cell(n_cal, 1);
for k = 1:n_cal
    parts{k} = evaluate_cands_one_sample( ...
        cfg, S60, sgrid, unit, cands, carriers, selected_variants, ...
        cal_manifest(k, :), up_signals{k}, gt_cache, fast_ctx, use_fast);
end
details = vertcat(parts{:});
end

function carriers = build_unit_carrier_map(sgrid, S60, unit, cands)
n_c = height(cands);
carriers = cell(n_c, 1);
if ~ismember(unit.ThresholdType, ["RSFT", "ASFT"]) || n_c == 0
    return;
end
freq_pairs = [cands.FrOverBr, cands.FaOverBa, cands.Phi0];
[uniq_pairs, ~, idx_map] = unique(freq_pairs, "rows", "stable");
uniq_carriers = cell(size(uniq_pairs, 1), 1);
for k = 1:size(uniq_pairs, 1)
    cparams = struct("fr_over_Br", uniq_pairs(k, 1), "fa_over_Ba", uniq_pairs(k, 2), "phi0", uniq_pairs(k, 3));
    uniq_carriers{k} = build_sft_carrier(sgrid, cparams, S60);
end
for k = 1:n_c
    carriers{k} = uniq_carriers{idx_map(k)};
end
end

function rows = evaluate_cands_one_sample( ...
    cfg, S60, sgrid, unit, cands, carriers, selected_variants, sample_row, sig_up, gt_cache, fast_ctx, use_fast)

n_c = height(cands);
parts = cell(n_c, 1);
ttype = string(unit.ThresholdType);
axis_mode = string(unit.Axis);

for c_idx = 1:n_c
    crow = cands(c_idx, :);
    t_spec = struct("ThresholdType", ttype);
    if ttype == "RT"
        t_spec.As = crow.As;
    elseif ismember(ttype, ["RSFT", "ASFT"])
        t_spec.STR_dB = crow.STR_dB;
        t_spec.fr_over_Br = crow.FrOverBr;
        t_spec.fa_over_Ba = crow.FaOverBa;
        t_spec.phi0 = crow.Phi0;
    end
    sub = zc1bit_stage2_eval_threshold( ...
        cfg, S60, sgrid, t_spec, sample_row, sig_up, selected_variants, axis_mode, ...
        gt_cache, carriers{c_idx}, fast_ctx, use_fast, crow.ParamKey);
    n_sub = height(sub);
    sub.ParamTag = repmat(crow.ParamTag, n_sub, 1);
    sub.ParamKey = repmat(crow.ParamKey, n_sub, 1);
    parts{c_idx} = sub;
end
rows = vertcat(parts{:});
end

function sum_tbl = aggregate_unit_summary(details, unit, cands, selected_variants)
n_c = height(cands);
n_v = height(selected_variants);
n_rows = n_c * n_v;

UnitKey = repmat(string(unit.UnitKey), n_rows, 1);
Branch = repmat(string(unit.Branch), n_rows, 1);
Axis = repmat(string(unit.Axis), n_rows, 1);
GroupName = repmat(string(unit.GroupName), n_rows, 1);
Q = repmat(double(unit.Q), n_rows, 1);
Range_q = repmat(double(unit.Range_q), n_rows, 1);
Azimuth_q = repmat(double(unit.Azimuth_q), n_rows, 1);
ThresholdType = repmat(string(unit.ThresholdType), n_rows, 1);

ParamTag = strings(n_rows, 1);
ParamKey = strings(n_rows, 1);
As = nan(n_rows, 1);
STR_dB = nan(n_rows, 1);
FrOverBr = nan(n_rows, 1);
FaOverBa = nan(n_rows, 1);
Phi0 = zeros(n_rows, 1);
PhaseTag = strings(n_rows, 1);

VariantOrder = zeros(n_rows, 1);
VariantRole = strings(n_rows, 1);
VariantKey = strings(n_rows, 1);
RuleKey = strings(n_rows, 1);
DecisionMode = strings(n_rows, 1);
RemoveDC = false(n_rows, 1);

SampleCount = zeros(n_rows, 1);
ENLSampleCount = zeros(n_rows, 1);
Calibration_SSIM_Mean = zeros(n_rows, 1);
Calibration_SSIM_Std = zeros(n_rows, 1);
Calibration_PSNR_Mean = zeros(n_rows, 1);
Calibration_PSNR_Std = zeros(n_rows, 1);
Calibration_Entropy_Mean = zeros(n_rows, 1);
Calibration_ENL_Mean = nan(n_rows, 1);
Calibration_GT_ENL_Mean = nan(n_rows, 1);
Calibration_AbsENLDiff_Mean = nan(n_rows, 1);
Calibration_ZeroRatio_Mean = zeros(n_rows, 1);
Calibration_AbsMeanBias_Mean = zeros(n_rows, 1);

ptr = 0;
for c_idx = 1:n_c
    crow = cands(c_idx, :);
    cmask = details.ParamTag == crow.ParamTag;
    cdet = details(cmask, :);
    for v_idx = 1:n_v
        vrow = selected_variants(v_idx, :);
        sub = cdet(cdet.VariantKey == vrow.VariantKey, :);
        if isempty(sub)
            error("zc1bit:MissingSearchDetail", "搜索汇总缺少单元 %s / 参数 %s / 变体 %s。", ...
                unit.UnitKey, crow.ParamTag, vrow.VariantKey);
        end
        ptr = ptr + 1;
        ParamTag(ptr) = crow.ParamTag;
        ParamKey(ptr) = crow.ParamKey;
        As(ptr) = crow.As;
        STR_dB(ptr) = crow.STR_dB;
        FrOverBr(ptr) = crow.FrOverBr;
        FaOverBa(ptr) = crow.FaOverBa;
        Phi0(ptr) = crow.Phi0;
        PhaseTag(ptr) = crow.PhaseTag;

        VariantOrder(ptr) = vrow.VariantOrder;
        VariantRole(ptr) = vrow.VariantRole;
        VariantKey(ptr) = vrow.VariantKey;
        RuleKey(ptr) = vrow.RuleKey;
        DecisionMode(ptr) = vrow.DecisionMode;
        RemoveDC(ptr) = logical(vrow.RemoveDC);

        SampleCount(ptr) = height(sub);
        enl_mask = ~isnan(sub.ENL);
        ENLSampleCount(ptr) = sum(enl_mask);
        Calibration_SSIM_Mean(ptr) = mean(sub.SSIM);
        Calibration_SSIM_Std(ptr) = std(sub.SSIM);
        Calibration_PSNR_Mean(ptr) = mean(sub.PSNR);
        Calibration_PSNR_Std(ptr) = std(sub.PSNR);
        Calibration_Entropy_Mean(ptr) = mean(sub.Entropy);
        if any(enl_mask)
            Calibration_ENL_Mean(ptr) = mean(sub.ENL(enl_mask));
            Calibration_GT_ENL_Mean(ptr) = mean(sub.GT_ENL(enl_mask));
            Calibration_AbsENLDiff_Mean(ptr) = mean(sub.AbsENLDiff(enl_mask));
        end
        Calibration_ZeroRatio_Mean(ptr) = mean(sub.ZeroRatio);
        Calibration_AbsMeanBias_Mean(ptr) = mean(sub.AbsMeanBias);
    end
end

sum_tbl = table( ...
    UnitKey, Branch, Axis, GroupName, Q, Range_q, Azimuth_q, ThresholdType, ...
    ParamTag, ParamKey, As, STR_dB, FrOverBr, FaOverBa, Phi0, PhaseTag, ...
    VariantOrder, VariantRole, VariantKey, RuleKey, DecisionMode, RemoveDC, ...
    SampleCount, ENLSampleCount, ...
    Calibration_SSIM_Mean, Calibration_SSIM_Std, Calibration_PSNR_Mean, Calibration_PSNR_Std, ...
    Calibration_Entropy_Mean, Calibration_ENL_Mean, Calibration_GT_ENL_Mean, ...
    Calibration_AbsENLDiff_Mean, Calibration_ZeroRatio_Mean, Calibration_AbsMeanBias_Mean);
end

function best_tbl = select_best_per_variant(sum_tbl, selected_variants, cand_meta)
n_v = height(selected_variants);
parts = cell(n_v, 1);
for v_idx = 1:n_v
    vkey = selected_variants.VariantKey(v_idx);
    sub = sum_tbl(sum_tbl.VariantKey == vkey, :);
    if isempty(sub)
        error("zc1bit:MissingVariantSummary", "未找到变体 %s 的搜索结果。", vkey);
    end
    sub = sortrows(sub, ...
        {'Calibration_SSIM_Mean', 'Calibration_PSNR_Mean', 'Calibration_AbsENLDiff_Mean', 'Calibration_ZeroRatio_Mean', 'ParamTag'}, ...
        {'descend', 'descend', 'ascend', 'ascend', 'ascend'});
    winner = sub(1, :);
    [t_rt, t_str, t_freq] = check_boundary_touch(winner, cand_meta);
    winner.TouchesRTBoundary = t_rt;
    winner.TouchesSTRBoundary = t_str;
    winner.TouchesFreqBoundary = t_freq;
    winner.TouchesAnyBoundary = t_rt || t_str || t_freq;
    parts{v_idx} = winner;
end
best_tbl = vertcat(parts{:});
end

function [t_rt, t_str, t_freq] = check_boundary_touch(row, cand_meta)
t_rt = false;
t_str = false;
t_freq = false;
ttype = string(row.ThresholdType);
tol = 1e-9;
if ttype == "RT" && ~isempty(cand_meta.as_grid)
    t_rt = abs(row.As - min(cand_meta.as_grid)) < tol || abs(row.As - max(cand_meta.as_grid)) < tol;
elseif ismember(ttype, ["RSFT", "ASFT"])
    if ~isempty(cand_meta.str_grid)
        t_str = row.STR_dB <= min(cand_meta.str_grid) + tol || row.STR_dB >= max(cand_meta.str_grid) - tol;
    end
    if ttype == "RSFT"
        f_val = row.FrOverBr;
    else
        f_val = row.FaOverBa;
    end
    if ~isempty(cand_meta.freq_grid)
        t_freq = f_val <= min(cand_meta.freq_grid) + tol || f_val >= max(cand_meta.freq_grid) - tol;
    end
end
end

function save_unit_checkpoint(cp_dir, unit, selected_variants, cal_manifest, S60, cfg, summary_rows, best_rows)
cp_path = fullfile(cp_dir, unit.UnitKey + ".mat");
expected_params = build_unit_expected_params(unit, selected_variants, cal_manifest, S60, cfg);
state = struct( ...
    "parameters", expected_params, ...
    "summary_rows", summary_rows, ...
    "best_rows", best_rows, ...
    "completed", true);
zc1bit_atomic_save(cp_path, struct("state", state));
fprintf("  [完成] %s（已原子落盘 %d 条校准搜索汇总 / %d 条最优参数）\n", ...
    unit.UnitKey, height(summary_rows), height(best_rows));
end

function expected_params = build_unit_expected_params(unit, selected_variants, cal_manifest, S60, cfg)
expected_params = struct( ...
    "unit", zc1bit_parameters(unit), ...
    "selected_variants", zc1bit_parameters(selected_variants), ...
    "stage2", zc1bit_parameters(cfg.stage2), ...
    "nyquist_margin", cfg.nyquist_margin, ...
    "seed", cfg.seed, ...
    "manifest", zc1bit_parameters(cal_manifest), ...
    "S60", zc1bit_parameters(S60), ...
    "enl_window_size", cfg.enl_window_size, ...
    "enl_stride", cfg.enl_stride, ...
    "entropy_num_bins", cfg.entropy_num_bins);
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

function validate_search_outputs(summary_table, best_params_table, units, selected_variants)
if isempty(summary_table) || height(summary_table) == 0
    error("zc1bit:EmptySearchSummary", "阶段二搜索汇总表为空。");
end
expected_best = height(units) * height(selected_variants);
if height(best_params_table) ~= expected_best
    error("zc1bit:IncompleteSearchBest", ...
        "阶段二最优参数数量 (%d) 与预期 (%d) 不一致。", height(best_params_table), expected_best);
end
for col = ["Calibration_SSIM_Mean", "Calibration_PSNR_Mean", "Calibration_Entropy_Mean", "Calibration_ZeroRatio_Mean"]
    vals = summary_table.(col);
    if any(~isfinite(vals))
        error("zc1bit:NonFiniteSearchSummary", "搜索汇总列 %s 存在非有限值。", col);
    end
end
end

function write_search_metadata(search_dir, cfg, cal_manifest, units, selected_variants, best_params_table, cache_gb)
fid = fopen(fullfile(search_dir, "search_metadata.txt"), "w", "n", "UTF-8");
if fid < 0
    error("zc1bit:MetadataWrite", "无法写入 search_metadata.txt。");
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "Experiment=SAR-ZC1Bit Stage 2 Search\n");
fprintf(fid, "ExecutionMode=single_core\n");
fprintf(fid, "CacheBudgetGB=%.6g\n", cache_gb);
fprintf(fid, "CalibrationSampleCount=%d\n", height(cal_manifest));
fprintf(fid, "UnitCount=%d\n", height(units));
fprintf(fid, "SelectedVariants=%s\n", join(selected_variants.VariantKey, ","));
fprintf(fid, "BestParamRowCount=%d\n", height(best_params_table));
fprintf(fid, "SFTRefine=%d\n", logical(cfg.stage2.sft_refine));
end
