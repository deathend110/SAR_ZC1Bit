function out = zc1bit_evaluate_impl(cfg)
%ZC1BIT_EVALUATE_IMPL 阶段二（35 测试样本最终评估、带外泄漏诊断与双侧配对 Wilcoxon 检验）。
% 严格仅使用确认锁定参数与 35 个测试样本（偶数 LocalSampleIdx），绝不交叉污染。

if nargin == 0
    cfg = zc1bit_config();
end
directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>

record = zc1bit_prepare_run(cfg);
p = zc1bit_paths(cfg);
S60 = load(cfg.parameter_file_60);

confirm_mat = fullfile(p.locked, "confirm_results.mat");
if ~isfile(confirm_mat)
    error("zc1bit:MissingLockedConfirm", ...
        "未找到锁参确认结果：%s。请先运行 run_zc1bit_confirm(cfg)。", confirm_mat);
end
loaded = load(confirm_mat, "signature", "units", "selected_variants", "locked_confirm_plan");
if ~isfield(loaded, "locked_confirm_plan")
    error("zc1bit:InvalidConfirmResults", "confirm_results.mat 缺少 locked_confirm_plan。");
end

cal_manifest = record.manifest(record.manifest.Split == "calibration", :);
test_manifest = record.manifest(record.manifest.Split == "test", :);
if isempty(test_manifest)
    error("zc1bit:EmptyTestSplit", "样本清单中未找到测试样本。");
end

expected_sig = struct( ...
    "stage2", zc1bit_parameters(cfg.stage2), ...
    "nyquist_margin", cfg.nyquist_margin, ...
    "seed", cfg.seed, ...
    "S60", zc1bit_parameters(S60), ...
    "cal_manifest", zc1bit_parameters(cal_manifest));
if isfield(loaded, "signature")
    zc1bit_check_parameters(loaded.signature, expected_sig, "锁参确认与当前配置");
end

units = loaded.units;
selected_variants = loaded.selected_variants;
locked_plan = loaded.locked_confirm_plan;

[cache_gb, preload_signals] = resolve_runtime_options(cfg);
cache_budget_bytes = cache_gb * (1024 ^ 3);
use_fast = true;
if isfield(cfg.stage2, "use_fast_path")
    use_fast = logical(cfg.stage2.use_fast_path);
end

n_units = height(units);
n_test = height(test_manifest);
cp_dir = fullfile(p.checkpoints, "evaluate");

detail_by_unit = cell(n_units, 1);
leakage_by_unit = cell(n_units, 1);
completed_units = false(n_units, 1);

for u_idx = 1:n_units
    unit = units(u_idx, :);
    cp_path = fullfile(cp_dir, unit.UnitKey + ".mat");
    if ~isfile(cp_path)
        continue;
    end
    state = zc1bit_load_state(cp_path, struct());
    if isfield(state, "completed") && isequal(state.completed, true) && ...
            isfield(state, "detail_rows") && isfield(state, "leakage_rows")
        unit_plan = locked_plan(locked_plan.UnitKey == unit.UnitKey, :);
        expected_params = build_eval_expected_params(unit, unit_plan, test_manifest, S60, cfg);
        zc1bit_check_parameters(state.parameters, expected_params, ...
            sprintf("阶段二评估检查点(%s)", unit.UnitKey));
        detail_by_unit{u_idx} = state.detail_rows;
        leakage_by_unit{u_idx} = state.leakage_rows;
        completed_units(u_idx) = true;
    end
end

gt_cache_path = fullfile(cp_dir, "gt_test_cache.mat");
gt_cache = zc1bit_stage2_gt_cache(cfg, S60, test_manifest, {}, gt_cache_path);
test_signals = {};

if any(~completed_units) || isempty(gt_cache)
    test_signals = zc1bit_stage2_load_signals(test_manifest, S60, cfg.raw_decimation_step, "测试样本");
    if isempty(gt_cache)
        gt_cache = zc1bit_stage2_gt_cache(cfg, S60, test_manifest, test_signals, gt_cache_path);
    end
end

if ~isfolder(p.evaluate)
    mkdir(p.evaluate);
end
if ~isfolder(p.stats)
    mkdir(p.stats);
end
writetable(gt_cache.roi_manifest, fullfile(p.evaluate, "evaluate_test_enl_roi_manifest.csv"));

fprintf("\n=== 阶段二测试集评估：%d 单元 × %d 测试样本（计划 %d 条 | 已命中 %d/%d 检查点 | cache=%.1f GB）===\n", ...
    n_units, n_test, height(locked_plan), sum(completed_units), n_units, cache_gb);

pending_indices = find(~completed_units);
if ~isempty(pending_indices)
    fast_ctx = [];
    if use_fast
        fast_ctx = build_fast_focus_context( ...
            S60, units.Range_q(pending_indices), units.Azimuth_q(pending_indices));
    end

    [diag_sample_indices, diag_rc_gts] = build_leakage_gt_refs(test_manifest, test_signals, S60);
    diag_groups = select_diagnostic_groups(units);

    base_info = whos("test_signals");
    base_bytes = base_info.bytes;
    if isa(test_signals{1}, "single")
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
                unit = units(u_idx, :);
                uplan = locked_plan(locked_plan.UnitKey == unit.UnitKey, :);
                [drows, lrows] = evaluate_unit_cached_up( ...
                    cfg, S60, sgrid, unit, uplan, test_manifest, test_signals, ...
                    gt_cache, fast_ctx, use_fast, diag_groups, diag_sample_indices, diag_rc_gts);
                save_eval_checkpoint(cp_dir, unit, uplan, test_manifest, S60, cfg, drows, lrows);
                detail_by_unit{u_idx} = drows;
                leakage_by_unit{u_idx} = lrows;
            end
            continue;
        end

        est_up_bytes = n_test * prod(sgrid.upsampled_size) * elem_bytes;
        if preload_signals && (base_bytes + est_up_bytes <= cache_budget_bytes)
            fprintf("  [R%dA%d] 内存缓存模式：预构建 %d 个测试上采样回波（约 %.2f GB）...\n", ...
                rq, aq, n_test, est_up_bytes / (1024 ^ 3));
            up_signals = cell(n_test, 1);
            for k = 1:n_test
                up_signals{k} = two_dim_upsample_fft(test_signals{k}, rq, aq);
            end
            for u_idx = group_u_idx(:).'
                unit = units(u_idx, :);
                uplan = locked_plan(locked_plan.UnitKey == unit.UnitKey, :);
                [drows, lrows] = evaluate_unit_cached_up( ...
                    cfg, S60, sgrid, unit, uplan, test_manifest, up_signals, ...
                    gt_cache, fast_ctx, use_fast, diag_groups, diag_sample_indices, diag_rc_gts);
                save_eval_checkpoint(cp_dir, unit, uplan, test_manifest, S60, cfg, drows, lrows);
                detail_by_unit{u_idx} = drows;
                leakage_by_unit{u_idx} = lrows;
            end
            clear up_signals
        else
            fprintf("  [R%dA%d] 流式低内存模式：预计整组上采样 %.2f GB 超出预算，逐样本复用（%d 单元）...\n", ...
                rq, aq, est_up_bytes / (1024 ^ 3), numel(group_u_idx));
            n_g = numel(group_u_idx);
            d_parts = cell(n_g, n_test);
            l_parts = cell(n_g, numel(diag_sample_indices));
            for k = 1:n_test
                sig_up = two_dim_upsample_fft(test_signals{k}, rq, aq);
                srow = test_manifest(k, :);
                diag_pos = find(diag_sample_indices == k, 1);
                for m = 1:n_g
                    unit = units(group_u_idx(m), :);
                    uplan = locked_plan(locked_plan.UnitKey == unit.UnitKey, :);
                    [d_parts{m, k}, l_row] = evaluate_unit_one_sample( ...
                        cfg, S60, sgrid, unit, uplan, [], srow, sig_up, gt_cache, fast_ctx, use_fast, ...
                        diag_groups, ~isempty(diag_pos), diag_rc_gts, diag_pos);
                    if ~isempty(diag_pos)
                        l_parts{m, diag_pos} = l_row;
                    end
                end
                clear sig_up
            end
            for m = 1:n_g
                u_idx = group_u_idx(m);
                unit = units(u_idx, :);
                uplan = locked_plan(locked_plan.UnitKey == unit.UnitKey, :);
                drows = vertcat(d_parts{m, :});
                lrows = vertcat_nonempty(l_parts(m, :), empty_leakage_table());
                save_eval_checkpoint(cp_dir, unit, uplan, test_manifest, S60, cfg, drows, lrows);
                detail_by_unit{u_idx} = drows;
                leakage_by_unit{u_idx} = lrows;
            end
        end
    end
end

detail_table = vertcat(detail_by_unit{:});
leakage_table = vertcat_nonempty(leakage_by_unit, empty_leakage_table());
validate_eval_detail(detail_table);

summary_table = build_evaluate_summary(detail_table, locked_plan);
wilcoxon_table = build_paired_wilcoxon_table(detail_table, units, selected_variants, locked_plan, cfg);

writetable(detail_table, fullfile(p.evaluate, "evaluate_detail.csv"));
writetable(summary_table, fullfile(p.evaluate, "evaluate_summary.csv"));
writetable(leakage_table, fullfile(p.evaluate, "evaluate_leakage_summary.csv"));
writetable(wilcoxon_table, fullfile(p.stats, "paired_wilcoxon_vs_b0.csv"));

zc1bit_atomic_save(fullfile(p.evaluate, "evaluate_results.mat"), struct( ...
    "units", units, ...
    "selected_variants", selected_variants, ...
    "locked_confirm_plan", locked_plan, ...
    "test_manifest", test_manifest, ...
    "roi_manifest", gt_cache.roi_manifest, ...
    "detail_table", detail_table, ...
    "summary_table", summary_table, ...
    "leakage_table", leakage_table, ...
    "wilcoxon_table", wilcoxon_table));
write_evaluate_metadata(p.evaluate, test_manifest, units, selected_variants, summary_table, wilcoxon_table, cache_gb);

fprintf("阶段二测试集评估完成，输出目录：%s\n", p.evaluate);
out = struct( ...
    "detail_table", detail_table, ...
    "summary_table", summary_table, ...
    "leakage_table", leakage_table, ...
    "wilcoxon_table", wilcoxon_table, ...
    "test_manifest", test_manifest, ...
    "output_dir", p.evaluate, ...
    "stats_dir", p.stats);
end

function [diag_indices, diag_rc_gts] = build_leakage_gt_refs(test_manifest, test_signals, S60)
[~, diag_indices] = unique(test_manifest.Dataset, "stable");
diag_indices = sort(diag_indices(:).');
n_d = numel(diag_indices);
diag_rc_gts = cell(n_d, 1);
for d_idx = 1:n_d
    sig = test_signals{diag_indices(d_idx)};
    diag_rc_gts{d_idx} = Range_Compress(sig, S60.fc, S60.tnrn, S60.gama, S60.R0, S60.C, S60.Fs, S60.Tp);
end
end

function diag_groups = select_diagnostic_groups(units)
preferred = ["R1A1", "R4A1", "R9A1", "R1A4", "R1A9"];
avail = unique(units.GroupName, "stable");
diag_groups = avail(ismember(avail, preferred));
if isempty(diag_groups)
    diag_groups = avail([1, end]);
    diag_groups = unique(diag_groups, "stable");
end
end

function [detail_rows, leakage_rows] = evaluate_unit_cached_up( ...
    cfg, S60, sgrid, unit, uplan, test_manifest, up_signals, ...
    gt_cache, fast_ctx, use_fast, diag_groups, diag_sample_indices, diag_rc_gts)

carriers = build_plan_carrier_map(sgrid, S60, unit, uplan);
n_test = height(test_manifest);
d_parts = cell(n_test, 1);
l_parts = cell(numel(diag_sample_indices), 1);

for k = 1:n_test
    diag_pos = find(diag_sample_indices == k, 1);
    [d_parts{k}, lrow] = evaluate_unit_one_sample( ...
        cfg, S60, sgrid, unit, uplan, carriers, test_manifest(k, :), up_signals{k}, ...
        gt_cache, fast_ctx, use_fast, diag_groups, ~isempty(diag_pos), diag_rc_gts, diag_pos);
    if ~isempty(diag_pos)
        l_parts{diag_pos} = lrow;
    end
end

detail_rows = vertcat(d_parts{:});
leakage_rows = vertcat_nonempty(l_parts, empty_leakage_table());
end

function carriers = build_plan_carrier_map(sgrid, S60, unit, uplan)
[uniq_tags, first_idx] = unique(uplan.ParamTag, "stable");
n_u = numel(uniq_tags);
carriers = containers.Map("KeyType", "char", "ValueType", "any");
if ~ismember(unit.ThresholdType, ["RSFT", "ASFT"])
    return;
end
for k = 1:n_u
    row = uplan(first_idx(k), :);
    cparams = struct("fr_over_Br", row.FrOverBr, "fa_over_Ba", row.FaOverBa, "phi0", row.Phi0);
    carriers(char(uniq_tags(k))) = build_sft_carrier(sgrid, cparams, S60);
end
end

function [rows, lrows] = evaluate_unit_one_sample( ...
    cfg, S60, sgrid, unit, uplan, carrier_map, sample_row, sig_up, ...
    gt_cache, fast_ctx, use_fast, diag_groups, is_diag_sample, diag_rc_gts, diag_pos)

[uniq_tags, ~, pgroup] = unique(uplan.ParamTag, "stable");
n_p = numel(uniq_tags);
parts = cell(n_p, 1);
ttype = string(unit.ThresholdType);
axis_mode = string(unit.Axis);

need_leakage = is_diag_sample && ismember(string(unit.GroupName), diag_groups);
l_parts = cell(n_p, 1);
leak_thr = 0.35;
if isfield(cfg.stage2, "leakage_threshold_ratio")
    leak_thr = double(cfg.stage2.leakage_threshold_ratio);
end

for p_idx = 1:n_p
    sub_plan = uplan(pgroup == p_idx, :);
    p0 = sub_plan(1, :);
    t_spec = struct("ThresholdType", ttype);
    carrier = [];
    if ttype == "RT"
        t_spec.As = p0.As;
    elseif ismember(ttype, ["RSFT", "ASFT"])
        t_spec.STR_dB = p0.STR_dB;
        t_spec.fr_over_Br = p0.FrOverBr;
        t_spec.fa_over_Ba = p0.FaOverBa;
        t_spec.phi0 = p0.Phi0;
        if ~isempty(carrier_map) && isKey(carrier_map, char(p0.ParamTag))
            carrier = carrier_map(char(p0.ParamTag));
        end
    end

    res_tbl = zc1bit_stage2_eval_threshold( ...
        cfg, S60, sgrid, t_spec, sample_row, sig_up, sub_plan, axis_mode, ...
        gt_cache, carrier, fast_ctx, use_fast, p0.ParamKey);

    n_sub = height(res_tbl);
    res_tbl.LockedConfigKey = sub_plan.LockedConfigKey;
    res_tbl.Track = sub_plan.Track;
    res_tbl.VariantRole = sub_plan.VariantRole;
    res_tbl.UnitKey = repmat(string(unit.UnitKey), n_sub, 1);
    res_tbl.Branch = repmat(string(unit.Branch), n_sub, 1);
    res_tbl.Axis = repmat(axis_mode, n_sub, 1);
    res_tbl.GroupName = repmat(string(unit.GroupName), n_sub, 1);
    res_tbl.Q = repmat(double(unit.Q), n_sub, 1);
    res_tbl.Range_q = repmat(double(unit.Range_q), n_sub, 1);
    res_tbl.Azimuth_q = repmat(double(unit.Azimuth_q), n_sub, 1);
    res_tbl.ThresholdType = repmat(ttype, n_sub, 1);
    res_tbl.ParamTag = sub_plan.ParamTag;
    res_tbl.ParamKey = sub_plan.ParamKey;
    res_tbl.As = sub_plan.As;
    res_tbl.STR_dB = sub_plan.STR_dB;
    res_tbl.FrOverBr = sub_plan.FrOverBr;
    res_tbl.FaOverBa = sub_plan.FaOverBa;
    parts{p_idx} = res_tbl;

    if need_leakage
        l_parts{p_idx} = compute_group_leakage( ...
            cfg, S60, sgrid, t_spec, carrier, p0.ParamKey, sig_up, sub_plan, unit, ...
            sample_row, diag_rc_gts{diag_pos}, leak_thr);
    end
end

rows = vertcat(parts{:});
lrows = vertcat_nonempty(l_parts, empty_leakage_table());
end

function ltbl = compute_group_leakage( ...
    cfg, S60, sgrid, t_spec, carrier, seed_tag, sig_up, sub_plan, unit, sample_row, rc_gt, leak_thr)

ttype = string(unit.ThresholdType);
t_opts = struct("nyquist_margin", cfg.nyquist_margin, "carrier", carrier);
if ttype == "RT"
    rt_seed = zc1bit_stage2_rt_seed(cfg.seed, seed_tag, sample_row.SampleID);
    t_opts.rand_stream = RandStream("mt19937ar", "Seed", rt_seed);
end
U = build_threshold(sig_up, sgrid, t_spec, S60, t_opts);
c1 = quantize_1bit_with_u(sig_up, U);
clear U

n_sub = height(sub_plan);
off_vec = zeros(n_sub, 1);
rng_vec = zeros(n_sub, 1);
az_vec = zeros(n_sub, 1);
rq = double(unit.Range_q);
aq = double(unit.Azimuth_q);
axis_mode = string(unit.Axis);

for k = 1:n_sub
    qopts = struct("Axis", axis_mode, "DecisionMode", sub_plan.DecisionMode(k), "RemoveDC", logical(sub_plan.RemoveDC(k)));
    cp15 = apply_zc_pseudo15bit(c1, sub_plan.RuleKey(k), qopts);
    [~, rc_crop] = focus_upsampled_channel(cp15, S60, rq, aq);
    [off_vec(k), rng_vec(k), az_vec(k)] = compute_leakage_metrics(rc_crop, rc_gt, leak_thr);
end

ltbl = table( ...
    repmat(sample_row.SampleID, n_sub, 1), ...
    repmat(sample_row.Dataset, n_sub, 1), ...
    sub_plan.LockedConfigKey, ...
    repmat(string(unit.UnitKey), n_sub, 1), ...
    repmat(string(unit.GroupName), n_sub, 1), ...
    repmat(ttype, n_sub, 1), ...
    sub_plan.Track, ...
    sub_plan.VariantKey, ...
    sub_plan.ParamTag, ...
    off_vec, rng_vec, az_vec, ...
    'VariableNames', { ...
    'SampleID', 'Dataset', 'LockedConfigKey', 'UnitKey', 'GroupName', ...
    'ThresholdType', 'Track', 'VariantKey', 'ParamTag', ...
    'OffSupportRatio', 'RangeOffRatio', 'AzimuthOffRatio'});
end

function sum_tbl = build_evaluate_summary(detail_table, locked_plan)
n_p = height(locked_plan);
LockedConfigKey = locked_plan.LockedConfigKey;
UnitKey = locked_plan.UnitKey;
Branch = locked_plan.Branch;
Axis = locked_plan.Axis;
GroupName = locked_plan.GroupName;
Q = locked_plan.Q;
Range_q = locked_plan.Range_q;
Azimuth_q = locked_plan.Azimuth_q;
ThresholdType = locked_plan.ThresholdType;
Track = locked_plan.Track;
VariantOrder = locked_plan.VariantOrder;
VariantRole = locked_plan.VariantRole;
VariantKey = locked_plan.VariantKey;
RuleKey = locked_plan.RuleKey;
DecisionMode = locked_plan.DecisionMode;
RemoveDC = logical(locked_plan.RemoveDC);
ParamTag = locked_plan.ParamTag;
As = locked_plan.As;
STR_dB = locked_plan.STR_dB;
FrOverBr = locked_plan.FrOverBr;
FaOverBa = locked_plan.FaOverBa;

TestSampleCount = zeros(n_p, 1);
ENLSampleCount = zeros(n_p, 1);
PSNR_Mean = zeros(n_p, 1);
PSNR_Std = zeros(n_p, 1);
SSIM_Mean = zeros(n_p, 1);
SSIM_Std = zeros(n_p, 1);
Entropy_Mean = zeros(n_p, 1);
ENL_Mean = nan(n_p, 1);
GT_ENL_Mean = nan(n_p, 1);
AbsENLDiff_Mean = nan(n_p, 1);
ZeroRatio_Mean = zeros(n_p, 1);
AbsMeanBias_Mean = zeros(n_p, 1);

for k = 1:n_p
    sub = detail_table(detail_table.LockedConfigKey == LockedConfigKey(k), :);
    if isempty(sub)
        error("zc1bit:MissingEvaluateRows", "评估汇总缺少配置：%s", LockedConfigKey(k));
    end
    TestSampleCount(k) = height(sub);
    enl_mask = ~isnan(sub.ENL);
    ENLSampleCount(k) = sum(enl_mask);
    PSNR_Mean(k) = mean(sub.PSNR);
    PSNR_Std(k) = std(sub.PSNR);
    SSIM_Mean(k) = mean(sub.SSIM);
    SSIM_Std(k) = std(sub.SSIM);
    Entropy_Mean(k) = mean(sub.Entropy);
    if any(enl_mask)
        ENL_Mean(k) = mean(sub.ENL(enl_mask));
        GT_ENL_Mean(k) = mean(sub.GT_ENL(enl_mask));
        AbsENLDiff_Mean(k) = mean(sub.AbsENLDiff(enl_mask));
    end
    ZeroRatio_Mean(k) = mean(sub.ZeroRatio);
    AbsMeanBias_Mean(k) = mean(sub.AbsMeanBias);
end

sum_tbl = table( ...
    LockedConfigKey, UnitKey, Branch, Axis, GroupName, Q, Range_q, Azimuth_q, ...
    ThresholdType, Track, VariantOrder, VariantRole, VariantKey, RuleKey, DecisionMode, RemoveDC, ...
    ParamTag, As, STR_dB, FrOverBr, FaOverBa, ...
    TestSampleCount, ENLSampleCount, ...
    PSNR_Mean, PSNR_Std, SSIM_Mean, SSIM_Std, Entropy_Mean, ...
    ENL_Mean, GT_ENL_Mean, AbsENLDiff_Mean, ZeroRatio_Mean, AbsMeanBias_Mean);
end

function w_tbl = build_paired_wilcoxon_table(detail_table, units, selected_variants, locked_plan, cfg)
base_key = string(cfg.stage2.baseline_variant);
cand_variants = selected_variants(selected_variants.VariantKey ~= base_key, :);

rows = cell(0, 1);
for u_idx = 1:height(units)
    unit = units(u_idx, :);
    ukey = unit.UnitKey;
    b0_rows = sortrows(detail_table( ...
        detail_table.UnitKey == ukey & detail_table.Track == "independent_recalibrated" & ...
        detail_table.VariantKey == base_key, :), "SampleID");

    for v_idx = 1:height(cand_variants)
        vrow = cand_variants(v_idx, :);
        vkey = vrow.VariantKey;
        ind_rows = sortrows(detail_table( ...
            detail_table.UnitKey == ukey & detail_table.Track == "independent_recalibrated" & ...
            detail_table.VariantKey == vkey, :), "SampleID");

        if ~isempty(ind_rows) && isequal(ind_rows.SampleID, b0_rows.SampleID)
            rows{end + 1, 1} = make_wilcoxon_row("vs_B0", unit, "independent_recalibrated", vrow, ind_rows, b0_rows); %#ok<AGROW>
        end

        has_tr = any(locked_plan.UnitKey == ukey & locked_plan.Track == "b0_transferred" & locked_plan.VariantKey == vkey);
        if has_tr
            tr_rows = sortrows(detail_table( ...
                detail_table.UnitKey == ukey & detail_table.Track == "b0_transferred" & ...
                detail_table.VariantKey == vkey, :), "SampleID");
            if ~isempty(tr_rows) && isequal(tr_rows.SampleID, b0_rows.SampleID)
                rows{end + 1, 1} = make_wilcoxon_row("vs_B0", unit, "b0_transferred", vrow, tr_rows, b0_rows); %#ok<AGROW>
            end
            if ~isempty(ind_rows) && ~isempty(tr_rows) && isequal(ind_rows.SampleID, tr_rows.SampleID)
                rows{end + 1, 1} = make_wilcoxon_row("recal_vs_transfer", unit, "ablation", vrow, ind_rows, tr_rows); %#ok<AGROW>
            end
        end
    end
end

if isempty(rows)
    error("zc1bit:EmptyWilcoxon", "未生成任何配对 Wilcoxon 检验结果。");
end
w_tbl = vertcat(rows{:});
end

function row = make_wilcoxon_row(comp_type, unit, track_label, vrow, a_rows, b_rows)
wp = compute_paired_wilcoxon(a_rows.PSNR, b_rows.PSNR);
ws = compute_paired_wilcoxon(a_rows.SSIM, b_rows.SSIM);

enl_mask = ~isnan(a_rows.AbsENLDiff) & ~isnan(b_rows.AbsENLDiff);
n_enl = sum(enl_mask);
if n_enl > 0
    we = compute_paired_wilcoxon(a_rows.AbsENLDiff(enl_mask), b_rows.AbsENLDiff(enl_mask));
    d_enl_mean = we.MeanDiff;
    d_enl_std = we.StdDiff;
    d_enl_med = we.MedianDiff;
    p_enl = we.PValue;
else
    d_enl_mean = NaN;
    d_enl_std = NaN;
    d_enl_med = NaN;
    p_enl = NaN;
end

row = table( ...
    string(comp_type), string(unit.UnitKey), string(unit.Branch), string(unit.Axis), ...
    string(unit.GroupName), double(unit.Q), string(unit.ThresholdType), string(track_label), ...
    string(vrow.VariantRole), string(vrow.VariantKey), string(vrow.RuleKey), ...
    wp.N, n_enl, ...
    wp.MeanDiff, wp.StdDiff, wp.MedianDiff, wp.PValue, wp.SignedRank, ...
    ws.MeanDiff, ws.StdDiff, ws.MedianDiff, ws.PValue, ws.SignedRank, ...
    d_enl_mean, d_enl_std, d_enl_med, p_enl, ...
    'VariableNames', { ...
    'ComparisonType', 'UnitKey', 'Branch', 'Axis', 'GroupName', 'Q', 'ThresholdType', 'Track', ...
    'VariantRole', 'VariantKey', 'RuleKey', 'SampleCount', 'ENLSampleCount', ...
    'DeltaPSNR_Mean', 'DeltaPSNR_Std', 'DeltaPSNR_Median', 'P_PSNR', 'SignedRank_PSNR', ...
    'DeltaSSIM_Mean', 'DeltaSSIM_Std', 'DeltaSSIM_Median', 'P_SSIM', 'SignedRank_SSIM', ...
    'DeltaAbsENLDiff_Mean', 'DeltaAbsENLDiff_Std', 'DeltaAbsENLDiff_Median', 'P_AbsENLDiff'});
end

function tbl = empty_leakage_table()
tbl = table(zeros(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
    strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    'VariableNames', { ...
    'SampleID', 'Dataset', 'LockedConfigKey', 'UnitKey', 'GroupName', ...
    'ThresholdType', 'Track', 'VariantKey', 'ParamTag', ...
    'OffSupportRatio', 'RangeOffRatio', 'AzimuthOffRatio'});
end

function out_tbl = vertcat_nonempty(cells_in, fallback_empty)
non_empty = false(numel(cells_in), 1);
for k = 1:numel(cells_in)
    non_empty(k) = ~isempty(cells_in{k}) && height(cells_in{k}) > 0;
end
if any(non_empty)
    out_tbl = vertcat(cells_in{non_empty});
else
    out_tbl = fallback_empty;
end
end

function save_eval_checkpoint(cp_dir, unit, uplan, test_manifest, S60, cfg, detail_rows, leakage_rows)
cp_path = fullfile(cp_dir, unit.UnitKey + ".mat");
expected_params = build_eval_expected_params(unit, uplan, test_manifest, S60, cfg);
state = struct( ...
    "parameters", expected_params, ...
    "detail_rows", detail_rows, ...
    "leakage_rows", leakage_rows, ...
    "completed", true);
zc1bit_atomic_save(cp_path, struct("state", state));
fprintf("  [完成] %s（已原子落盘 %d 条测试评估明细）\n", unit.UnitKey, height(detail_rows));
end

function expected_params = build_eval_expected_params(unit, unit_plan, test_manifest, S60, cfg)
expected_params = struct( ...
    "unit", zc1bit_parameters(unit), ...
    "locked_plan", zc1bit_parameters(unit_plan), ...
    "stage2", zc1bit_parameters(cfg.stage2), ...
    "nyquist_margin", cfg.nyquist_margin, ...
    "seed", cfg.seed, ...
    "manifest", zc1bit_parameters(test_manifest), ...
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

function validate_eval_detail(detail_table)
if isempty(detail_table) || height(detail_table) == 0
    error("zc1bit:EmptyEvaluateDetail", "阶段二测试评估明细表为空。");
end
for col = ["PSNR", "SSIM", "Entropy", "ZeroRatio", "MeanBiasReal", "MeanBiasImag", "AbsMeanBias"]
    vals = detail_table.(col);
    if any(isnan(vals)) || any(isinf(vals) & vals < 0)
        error("zc1bit:InvalidEvaluateValue", "测试评估明细列 %s 包含非法值。", col);
    end
end
if any(~isfinite(detail_table.SSIM)) || any(~isfinite(detail_table.Entropy))
    error("zc1bit:NonFiniteEvaluateDetail", "测试评估明细表 SSIM 或 Entropy 存在非有限值。");
end
end

function write_evaluate_metadata(eval_dir, test_manifest, units, selected_variants, summary_table, wilcoxon_table, cache_gb)
fid = fopen(fullfile(eval_dir, "evaluate_metadata.txt"), "w", "n", "UTF-8");
if fid < 0
    error("zc1bit:MetadataWrite", "无法写入 evaluate_metadata.txt。");
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "Experiment=SAR-ZC1Bit Stage 2 Evaluate\n");
fprintf(fid, "ExecutionMode=single_core\n");
fprintf(fid, "CacheBudgetGB=%.6g\n", cache_gb);
fprintf(fid, "TestSampleCount=%d\n", height(test_manifest));
fprintf(fid, "UnitCount=%d\n", height(units));
fprintf(fid, "SelectedVariants=%s\n", join(selected_variants.VariantKey, ","));
fprintf(fid, "SummaryRowCount=%d\n", height(summary_table));
fprintf(fid, "WilcoxonRowCount=%d\n", height(wilcoxon_table));
end
