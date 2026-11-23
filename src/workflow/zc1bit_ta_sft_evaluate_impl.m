function out = zc1bit_ta_sft_evaluate_impl(cfg)
%ZC1BIT_TA_SFT_EVALUATE_IMPL TA-SFT 35 测试样本评估、带外泄漏、Wilcoxon 检验与初相失配消融。
% 严格仅使用 35 个测试样本（偶数 LocalSampleIdx），绝不修改或覆盖阶段一/二已有产物。

if nargin == 0
    cfg = zc1bit_config();
end
directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>

record = zc1bit_prepare_run(cfg);
p = zc1bit_paths(cfg);
S60 = load(cfg.parameter_file_60);

search_dir = fullfile(p.ta_sft, "search");
search_mat = fullfile(search_dir, "ta_sft_search_results.mat");
if ~isfile(search_mat)
    error("zc1bit:MissingTaSftSearch", ...
        "未找到 TA-SFT 搜索结果：%s。请先运行 run_zc1bit_ta_sft(cfg, ""search"")。", search_mat);
end
loaded = load(search_mat, "signature", "units", "b0_sft_params", "method_order", "locked_params");
if ~isfield(loaded, "locked_params") || ~isfield(loaded, "units") || ~isfield(loaded, "b0_sft_params")
    error("zc1bit:InvalidTaSftSearch", "ta_sft_search_results.mat 缺少必要字段。");
end

cal_manifest = record.manifest(record.manifest.Split == "calibration", :);
test_manifest = record.manifest(record.manifest.Split == "test", :);
if isempty(test_manifest)
    error("zc1bit:EmptyTestSplit", "样本清单中未找到测试样本。");
end

units = loaded.units;
b0_sft_params = loaded.b0_sft_params;
method_order = string(loaded.method_order(:));
locked_params = loaded.locked_params;

if isfield(loaded, "signature")
    expected_search_sig = struct( ...
        "ta_sft", build_ta_sft_search_sig(cfg), ...
        "nyquist_margin", cfg.nyquist_margin, ...
        "seed", cfg.seed, ...
        "S60", zc1bit_parameters(S60), ...
        "cal_manifest", zc1bit_parameters(cal_manifest), ...
        "b0_sft_params", zc1bit_parameters(b0_sft_params));
    check_ta_sft_parameters(loaded.signature, expected_search_sig, "TA-SFT 搜索结果与当前配置");
end

if height(locked_params) ~= height(units) * numel(method_order)
    error("zc1bit:IncompleteTaSftLocked", "TA-SFT 锁参行数不完整。");
end

[cache_gb, preload_signals] = resolve_runtime_options(cfg);
cache_budget_bytes = cache_gb * (1024 ^ 3);
use_fast = true;
if isfield(cfg.ta_sft, "use_fast_path")
    use_fast = logical(cfg.ta_sft.use_fast_path);
end

eval_cp_dir = fullfile(p.checkpoints, "ta_sft", "evaluate");
rob_cp_dir = fullfile(p.checkpoints, "ta_sft", "robustness");
eval_out_dir = fullfile(p.ta_sft, "evaluate");
stats_out_dir = fullfile(p.ta_sft, "stats");
rob_out_dir = fullfile(p.ta_sft, "robustness");
ensure_dir(eval_cp_dir);
ensure_dir(rob_cp_dir);
ensure_dir(eval_out_dir);
ensure_dir(stats_out_dir);
ensure_dir(rob_out_dir);

n_units = height(units);
n_test = height(test_manifest);
detail_by_unit = cell(n_units, 1);
leakage_by_unit = cell(n_units, 1);
completed_eval = false(n_units, 1);

for u_idx = 1:n_units
    unit = units(u_idx, :);
    cp_path = fullfile(eval_cp_dir, unit.UnitKey + ".mat");
    if ~isfile(cp_path)
        continue;
    end
    state = zc1bit_load_state(cp_path, struct());
    if isfield(state, "completed") && isequal(state.completed, true) && ...
            isfield(state, "detail_rows") && isfield(state, "leakage_rows")
        u_locked = locked_params(locked_params.UnitKey == unit.UnitKey, :);
        b0_row = b0_sft_params(b0_sft_params.UnitKey == unit.UnitKey, :);
        expected_params = build_eval_expected_params(unit, b0_row, u_locked, test_manifest, S60, cfg);
        check_ta_sft_parameters(state.parameters, expected_params, ...
            sprintf("TA-SFT 评估检查点(%s)", unit.UnitKey));
        detail_by_unit{u_idx} = state.detail_rows;
        leakage_by_unit{u_idx} = state.leakage_rows;
        completed_eval(u_idx) = true;
    end
end

rob_units = units(ismember(units.UnitKey, string(cfg.ta_sft.phase_mismatch_units)), :);
n_rob = height(rob_units);
rob_by_unit = cell(n_rob, 1);
completed_rob = false(n_rob, 1);
for r_idx = 1:n_rob
    unit = rob_units(r_idx, :);
    cp_path = fullfile(rob_cp_dir, unit.UnitKey + ".mat");
    if ~isfile(cp_path)
        continue;
    end
    state = zc1bit_load_state(cp_path, struct());
    if isfield(state, "completed") && isequal(state.completed, true) && ...
            isfield(state, "robustness_rows")
        u_locked = locked_params(locked_params.UnitKey == unit.UnitKey, :);
        b0_row = b0_sft_params(b0_sft_params.UnitKey == unit.UnitKey, :);
        expected_params = build_rob_expected_params(unit, b0_row, u_locked, test_manifest, S60, cfg);
        check_ta_sft_parameters(state.parameters, expected_params, ...
            sprintf("TA-SFT 相偏鲁棒性检查点(%s)", unit.UnitKey));
        rob_by_unit{r_idx} = state.robustness_rows;
        completed_rob(r_idx) = true;
    end
end

gt_cache_path = fullfile(eval_cp_dir, "gt_test_cache.mat");
stage2_gt_path = fullfile(p.checkpoints, "evaluate", "gt_test_cache.mat");
gt_cache = zc1bit_stage2_gt_cache(cfg, S60, test_manifest, {}, gt_cache_path, stage2_gt_path);
test_signals = {};

if any(~completed_eval) || any(~completed_rob) || isempty(gt_cache)
    test_signals = zc1bit_stage2_load_signals(test_manifest, S60, cfg.raw_decimation_step, "测试样本");
    if isempty(gt_cache)
        gt_cache = zc1bit_stage2_gt_cache(cfg, S60, test_manifest, test_signals, gt_cache_path, stage2_gt_path);
    end
end

writetable(gt_cache.roi_manifest, fullfile(eval_out_dir, "ta_sft_test_enl_roi_manifest.csv"));

fprintf("\n=== TA-SFT 测试集评估：%d 单元 × %d 方法（%d 测试样本 | 主评估命中 %d/%d | 相偏命中 %d/%d | cache=%.1f GB）===\n", ...
    n_units, numel(method_order), n_test, sum(completed_eval), n_units, sum(completed_rob), n_rob, cache_gb);

pending_eval = find(~completed_eval);
if ~isempty(pending_eval)
    fast_ctx = [];
    if use_fast
        fast_ctx = build_fast_focus_context( ...
            S60, units.Range_q(pending_eval), units.Azimuth_q(pending_eval));
    end
    [diag_sample_indices, diag_rc_gts] = build_leakage_gt_refs(test_manifest, test_signals, S60);
    diag_groups = select_diag_groups(units, cfg);

    base_info = whos("test_signals");
    base_bytes = base_info.bytes;
    if isa(test_signals{1}, "single")
        elem_bytes = 8;
    else
        elem_bytes = 16;
    end

    pending_allocs = unique([units.Range_q(pending_eval), units.Azimuth_q(pending_eval)], "rows", "stable");
    for a_idx = 1:size(pending_allocs, 1)
        rq = pending_allocs(a_idx, 1);
        aq = pending_allocs(a_idx, 2);
        group_u_idx = pending_eval(units.Range_q(pending_eval) == rq & units.Azimuth_q(pending_eval) == aq);
        sgrid = resolve_sampling_grid(rq, aq, [S60.nrn, S60.nan], S60);

        if rq == 1 && aq == 1
            for u_idx = group_u_idx(:).'
                unit = units(u_idx, :);
                b0_row = b0_sft_params(b0_sft_params.UnitKey == unit.UnitKey, :);
                u_locked = order_unit_locked_rows(locked_params(locked_params.UnitKey == unit.UnitKey, :), method_order);
                [drows, lrows] = evaluate_unit_cached_up( ...
                    cfg, S60, sgrid, unit, b0_row, u_locked, test_manifest, test_signals, ...
                    gt_cache, fast_ctx, use_fast, diag_groups, diag_sample_indices, diag_rc_gts);
                save_eval_checkpoint(eval_cp_dir, unit, b0_row, u_locked, test_manifest, S60, cfg, drows, lrows);
                detail_by_unit{u_idx} = drows;
                leakage_by_unit{u_idx} = lrows;
            end
            continue;
        end

        est_up_bytes = n_test * prod(sgrid.upsampled_size) * elem_bytes;
        if preload_signals && (base_bytes + est_up_bytes <= cache_budget_bytes)
            fprintf("  [评估 R%dA%d] 内存缓存模式：预构建 %d 个测试上采样回波（约 %.2f GB）...\n", ...
                rq, aq, n_test, est_up_bytes / (1024 ^ 3));
            up_signals = cell(n_test, 1);
            for k = 1:n_test
                up_signals{k} = two_dim_upsample_fft(test_signals{k}, rq, aq);
            end
            for u_idx = group_u_idx(:).'
                unit = units(u_idx, :);
                b0_row = b0_sft_params(b0_sft_params.UnitKey == unit.UnitKey, :);
                u_locked = order_unit_locked_rows(locked_params(locked_params.UnitKey == unit.UnitKey, :), method_order);
                [drows, lrows] = evaluate_unit_cached_up( ...
                    cfg, S60, sgrid, unit, b0_row, u_locked, test_manifest, up_signals, ...
                    gt_cache, fast_ctx, use_fast, diag_groups, diag_sample_indices, diag_rc_gts);
                save_eval_checkpoint(eval_cp_dir, unit, b0_row, u_locked, test_manifest, S60, cfg, drows, lrows);
                detail_by_unit{u_idx} = drows;
                leakage_by_unit{u_idx} = lrows;
            end
            clear up_signals
        else
            fprintf("  [评估 R%dA%d] 流式低内存模式：逐样本复用上采样（%d 单元）...\n", rq, aq, numel(group_u_idx));
            n_g = numel(group_u_idx);
            d_parts = cell(n_g, n_test);
            l_parts = cell(n_g, numel(diag_sample_indices));
            b0_list = cell(n_g, 1);
            locked_list = cell(n_g, 1);
            carrier_list = cell(n_g, 1);
            for m = 1:n_g
                unit = units(group_u_idx(m), :);
                b0_r = b0_sft_params(b0_sft_params.UnitKey == unit.UnitKey, :);
                b0_list{m} = b0_r;
                locked_list{m} = order_unit_locked_rows(locked_params(locked_params.UnitKey == unit.UnitKey, :), method_order);
                carrier_list{m} = build_sft_carrier(sgrid, struct( ...
                    "fr_over_Br", b0_r.FrOverBr, "fa_over_Ba", b0_r.FaOverBa, "phi0", b0_r.Phi0), S60);
            end
            for k = 1:n_test
                sig_up = two_dim_upsample_fft(test_signals{k}, rq, aq);
                srow = test_manifest(k, :);
                diag_pos = find(diag_sample_indices == k, 1);
                for m = 1:n_g
                    unit = units(group_u_idx(m), :);
                    [d_parts{m, k}, lrow] = evaluate_unit_one_sample( ...
                        cfg, S60, sgrid, unit, b0_list{m}, carrier_list{m}, locked_list{m}, ...
                        srow, sig_up, gt_cache, fast_ctx, use_fast, ...
                        diag_groups, ~isempty(diag_pos), diag_rc_gts, diag_pos);
                    if ~isempty(diag_pos)
                        l_parts{m, diag_pos} = lrow;
                    end
                end
                clear sig_up
            end
            for m = 1:n_g
                u_idx = group_u_idx(m);
                unit = units(u_idx, :);
                drows = vertcat(d_parts{m, :});
                lrows = vertcat_nonempty(l_parts(m, :), empty_leakage_table());
                save_eval_checkpoint(eval_cp_dir, unit, b0_list{m}, locked_list{m}, test_manifest, S60, cfg, drows, lrows);
                detail_by_unit{u_idx} = drows;
                leakage_by_unit{u_idx} = lrows;
            end
        end
    end
end

pending_rob = find(~completed_rob);
if ~isempty(pending_rob)
    fast_ctx_rob = [];
    if use_fast
        fast_ctx_rob = build_fast_focus_context( ...
            S60, rob_units.Range_q(pending_rob), rob_units.Azimuth_q(pending_rob));
    end
    base_info = whos("test_signals");
    base_bytes = base_info.bytes;
    if isa(test_signals{1}, "single")
        elem_bytes = 8;
    else
        elem_bytes = 16;
    end

    for r_idx = pending_rob(:).'
        unit = rob_units(r_idx, :);
        rq = double(unit.Range_q);
        aq = double(unit.Azimuth_q);
        sgrid = resolve_sampling_grid(rq, aq, [S60.nrn, S60.nan], S60);
        b0_row = b0_sft_params(b0_sft_params.UnitKey == unit.UnitKey, :);
        u_locked = order_unit_locked_rows(locked_params(locked_params.UnitKey == unit.UnitKey, :), method_order);

        est_up_bytes = n_test * prod(sgrid.upsampled_size) * elem_bytes;
        use_up_cache = (rq == 1 && aq == 1) || (preload_signals && (base_bytes + est_up_bytes <= cache_budget_bytes));
        if rq == 1 && aq == 1
            rob_rows = run_robustness_unit( ...
                cfg, S60, sgrid, unit, b0_row, u_locked, method_order, test_manifest, test_signals, ...
                gt_cache, fast_ctx_rob, use_fast, false);
        elseif use_up_cache
            fprintf("  [相偏消融 %s] 内存缓存模式...\n", unit.UnitKey);
            up_signals = cell(n_test, 1);
            for k = 1:n_test
                up_signals{k} = two_dim_upsample_fft(test_signals{k}, rq, aq);
            end
            rob_rows = run_robustness_unit( ...
                cfg, S60, sgrid, unit, b0_row, u_locked, method_order, test_manifest, up_signals, ...
                gt_cache, fast_ctx_rob, use_fast, false);
            clear up_signals
        else
            fprintf("  [相偏消融 %s] 流式低内存模式...\n", unit.UnitKey);
            rob_rows = run_robustness_unit( ...
                cfg, S60, sgrid, unit, b0_row, u_locked, method_order, test_manifest, test_signals, ...
                gt_cache, fast_ctx_rob, use_fast, true);
        end
        save_rob_checkpoint(rob_cp_dir, unit, b0_row, u_locked, test_manifest, S60, cfg, rob_rows);
        rob_by_unit{r_idx} = rob_rows;
    end
end

detail_table = vertcat(detail_by_unit{:});
leakage_table = vertcat_nonempty(leakage_by_unit, empty_leakage_table());
robustness_table = vertcat_nonempty(rob_by_unit, empty_robustness_table());
validate_eval_detail(detail_table);

summary_table = build_evaluate_summary(detail_table, units, b0_sft_params, locked_params, method_order);
wilcoxon_table = build_paired_wilcoxon_table(detail_table, units, method_order);

writetable(detail_table, fullfile(eval_out_dir, "ta_sft_evaluate_detail.csv"));
writetable(summary_table, fullfile(eval_out_dir, "ta_sft_evaluate_summary.csv"));
writetable(leakage_table, fullfile(eval_out_dir, "ta_sft_leakage_summary.csv"));
writetable(wilcoxon_table, fullfile(stats_out_dir, "ta_sft_paired_wilcoxon.csv"));
writetable(robustness_table, fullfile(rob_out_dir, "ta_sft_phase_robustness_summary.csv"));

zc1bit_atomic_save(fullfile(eval_out_dir, "ta_sft_evaluate_results.mat"), struct( ...
    "units", units, ...
    "b0_sft_params", b0_sft_params, ...
    "method_order", method_order, ...
    "locked_params", locked_params, ...
    "test_manifest", test_manifest, ...
    "roi_manifest", gt_cache.roi_manifest, ...
    "detail_table", detail_table, ...
    "summary_table", summary_table, ...
    "leakage_table", leakage_table, ...
    "wilcoxon_table", wilcoxon_table, ...
    "robustness_table", robustness_table));
write_evaluate_metadata(eval_out_dir, test_manifest, units, method_order, summary_table, wilcoxon_table, robustness_table, cache_gb);

fprintf("TA-SFT 测试集评估与相偏消融完成，输出目录：%s\n", p.ta_sft);
out = struct( ...
    "detail_table", detail_table, ...
    "summary_table", summary_table, ...
    "leakage_table", leakage_table, ...
    "wilcoxon_table", wilcoxon_table, ...
    "robustness_table", robustness_table, ...
    "test_manifest", test_manifest, ...
    "output_dir", eval_out_dir, ...
    "stats_dir", stats_out_dir, ...
    "robustness_dir", rob_out_dir);
end

function u_locked = order_unit_locked_rows(raw_rows, method_order)
n_m = numel(method_order);
parts = cell(n_m, 1);
for m = 1:n_m
    hit = raw_rows(raw_rows.MethodKey == method_order(m), :);
    if height(hit) ~= 1
        error("zc1bit:MissingLockedMethod", "单元缺失或重复方法锁参：%s。", method_order(m));
    end
    parts{m} = hit;
end
u_locked = vertcat(parts{:});
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

function diag_groups = select_diag_groups(units, cfg)
preferred = string(cfg.ta_sft.leakage_groups);
avail = unique(units.GroupName, "stable");
diag_groups = avail(ismember(avail, preferred));
if isempty(diag_groups)
    diag_groups = unique(avail([1, end]), "stable");
end
end

function [detail_rows, leakage_rows] = evaluate_unit_cached_up( ...
    cfg, S60, sgrid, unit, b0_row, u_locked, test_manifest, up_signals, ...
    gt_cache, fast_ctx, use_fast, diag_groups, diag_sample_indices, diag_rc_gts)

cparams = struct("fr_over_Br", b0_row.FrOverBr, "fa_over_Ba", b0_row.FaOverBa, "phi0", b0_row.Phi0);
carrier = build_sft_carrier(sgrid, cparams, S60);

n_test = height(test_manifest);
d_parts = cell(n_test, 1);
l_parts = cell(numel(diag_sample_indices), 1);

for k = 1:n_test
    diag_pos = find(diag_sample_indices == k, 1);
    [d_parts{k}, lrow] = evaluate_unit_one_sample( ...
        cfg, S60, sgrid, unit, b0_row, carrier, u_locked, test_manifest(k, :), up_signals{k}, ...
        gt_cache, fast_ctx, use_fast, diag_groups, ~isempty(diag_pos), diag_rc_gts, diag_pos);
    if ~isempty(diag_pos)
        l_parts{diag_pos} = lrow;
    end
end

detail_rows = vertcat(d_parts{:});
leakage_rows = vertcat_nonempty(l_parts, empty_leakage_table());
end

function [rows, lrow] = evaluate_unit_one_sample( ...
    cfg, S60, sgrid, unit, b0_row, carrier, u_locked, sample_row, sig_up, ...
    gt_cache, fast_ctx, use_fast, diag_groups, is_diag_sample, diag_rc_gts, diag_pos)

sid = sample_row.SampleID;
gt_idx = find(gt_cache.SampleID == sid, 1);
if isempty(gt_idx)
    error("zc1bit:MissingSampleGT", "GT 缓存中缺少 SampleID=%d。", sid);
end
img_gt = gt_cache.img_gts{gt_idx};
is_enl = gt_cache.is_enl(gt_idx);
top_r = gt_cache.enl_top(gt_idx);
left_c = gt_cache.enl_left(gt_idx);
gt_enl_val = gt_cache.gt_enl(gt_idx);

rq = double(sgrid.q_range);
aq = double(sgrid.q_azimuth);
axis_mode = string(unit.Axis);
ttype = string(unit.ThresholdType);

t_spec = struct( ...
    "ThresholdType", ttype, ...
    "STR_dB", double(b0_row.STR_dB), ...
    "fr_over_Br", double(b0_row.FrOverBr), ...
    "fa_over_Ba", double(b0_row.FaOverBa), ...
    "phi0", double(b0_row.Phi0));
t_opts = struct("nyquist_margin", cfg.nyquist_margin, "carrier", carrier);
U = build_threshold(sig_up, sgrid, t_spec, S60, t_opts);
channel_1bit = quantize_1bit_with_u(sig_up, U);
clear U

n_m = height(u_locked);
psnr_vec = zeros(n_m, 1);
ssim_vec = zeros(n_m, 1);
ent_vec = zeros(n_m, 1);
enl_vec = nan(n_m, 1);
abs_enl_diff_vec = nan(n_m, 1);
zr_vec = zeros(n_m, 1);
mr_vec = zeros(n_m, 1);
eg_vec = nan(n_m, 1);
abm_vec = zeros(n_m, 1);

need_leakage = is_diag_sample && ismember(string(unit.GroupName), diag_groups);
leak_thr = 0.35;
if isfield(cfg.ta_sft, "leakage_threshold_ratio")
    leak_thr = double(cfg.ta_sft.leakage_threshold_ratio);
end
off_vec = nan(n_m, 1);
rng_vec = nan(n_m, 1);
az_vec = nan(n_m, 1);

ent_bins = cfg.entropy_num_bins;
enl_win = cfg.enl_window_size;

for m_idx = 1:n_m
    mkey = u_locked.MethodKey(m_idx);
    hp_name = u_locked.HyperparamName(m_idx);
    hp_val = double(u_locked.HyperparamValue(m_idx));
    ropts = struct("Axis", axis_mode);
    if hp_name ~= "none" && hp_name ~= "blind"
        ropts.(char(hp_name)) = hp_val;
    end
    [channel_rec, qdiag] = apply_ta_sft_reconstruct(channel_1bit, carrier, mkey, ropts);
    if use_fast
        img_out = fast_focus_complex_roi(channel_rec, fast_ctx, rq, aq);
    else
        img_out = focus_upsampled_channel(channel_rec, S60, rq, aq);
    end
    m = compute_image_metrics(img_out, img_gt);
    psnr_vec(m_idx) = m.PSNR;
    ssim_vec(m_idx) = m.SSIM;
    ent_vec(m_idx) = compute_entropy(img_out, ent_bins);
    if is_enl
        enl_vec(m_idx) = compute_enl(img_out, top_r, left_c, enl_win);
        abs_enl_diff_vec(m_idx) = abs(enl_vec(m_idx) - gt_enl_val);
    end
    zr_vec(m_idx) = qdiag.ZeroRatio;
    mr_vec(m_idx) = qdiag.ModifiedRatio;
    eg_vec(m_idx) = qdiag.EstGamma;
    abm_vec(m_idx) = qdiag.AbsMeanBias;

    if need_leakage
        [~, rc_crop] = focus_upsampled_channel(channel_rec, S60, rq, aq);
        [off_vec(m_idx), rng_vec(m_idx), az_vec(m_idx)] = compute_leakage_metrics(rc_crop, diag_rc_gts{diag_pos}, leak_thr);
    end
end

rows = table( ...
    repmat(sid, n_m, 1), ...
    repmat(sample_row.LocalSampleIdx, n_m, 1), ...
    repmat(sample_row.Dataset, n_m, 1), ...
    repmat(string(unit.UnitKey), n_m, 1), ...
    repmat(string(unit.Branch), n_m, 1), ...
    repmat(axis_mode, n_m, 1), ...
    repmat(string(unit.GroupName), n_m, 1), ...
    repmat(double(unit.Q), n_m, 1), ...
    repmat(rq, n_m, 1), ...
    repmat(aq, n_m, 1), ...
    repmat(ttype, n_m, 1), ...
    u_locked.MethodOrder, ...
    u_locked.MethodKey, ...
    u_locked.HyperparamName, ...
    u_locked.HyperparamValue, ...
    u_locked.HyperparamTag, ...
    psnr_vec, ssim_vec, ent_vec, enl_vec, repmat(gt_enl_val, n_m, 1), abs_enl_diff_vec, ...
    zr_vec, mr_vec, eg_vec, abm_vec, ...
    'VariableNames', { ...
    'SampleID', 'LocalSampleIdx', 'Dataset', 'UnitKey', 'Branch', 'Axis', 'GroupName', ...
    'Q', 'Range_q', 'Azimuth_q', 'ThresholdType', ...
    'MethodOrder', 'MethodKey', 'HyperparamName', 'HyperparamValue', 'HyperparamTag', ...
    'PSNR', 'SSIM', 'Entropy', 'ENL', 'GT_ENL', 'AbsENLDiff', ...
    'ZeroRatio', 'ModifiedRatio', 'EstGamma', 'AbsMeanBias'});

if need_leakage
    lrow = table( ...
        repmat(sid, n_m, 1), ...
        repmat(sample_row.Dataset, n_m, 1), ...
        repmat(string(unit.UnitKey), n_m, 1), ...
        repmat(string(unit.GroupName), n_m, 1), ...
        repmat(ttype, n_m, 1), ...
        u_locked.MethodOrder, ...
        u_locked.MethodKey, ...
        u_locked.HyperparamTag, ...
        off_vec, rng_vec, az_vec, ...
        'VariableNames', { ...
        'SampleID', 'Dataset', 'UnitKey', 'GroupName', 'ThresholdType', ...
        'MethodOrder', 'MethodKey', 'HyperparamTag', ...
        'OffSupportRatio', 'RangeOffRatio', 'AzimuthOffRatio'});
else
    lrow = empty_leakage_table();
end
end

function rob_rows = run_robustness_unit( ...
    cfg, S60, sgrid, unit, b0_row, u_locked, method_order, test_manifest, signals_in, ...
    gt_cache, fast_ctx, use_fast, streaming)

rq = double(sgrid.q_range);
aq = double(sgrid.q_azimuth);
axis_mode = string(unit.Axis);
ttype = string(unit.ThresholdType);
n_test = height(test_manifest);
n_m = numel(method_order);

cparams = struct("fr_over_Br", b0_row.FrOverBr, "fa_over_Ba", b0_row.FaOverBa, "phi0", b0_row.Phi0);
carrier_true = build_sft_carrier(sgrid, cparams, S60);
t_spec = struct( ...
    "ThresholdType", ttype, ...
    "STR_dB", double(b0_row.STR_dB), ...
    "fr_over_Br", double(b0_row.FrOverBr), ...
    "fa_over_Ba", double(b0_row.FaOverBa), ...
    "phi0", double(b0_row.Phi0));

dphi_list = double(cfg.ta_sft.phase_mismatch_deg(:).');
n_ang = numel(dphi_list);

psnr_mat = zeros(n_test, n_m, n_ang);
ssim_mat = zeros(n_test, n_m, n_ang);
zr_mat = zeros(n_test, n_m, n_ang);
mr_mat = zeros(n_test, n_m, n_ang);
eg_mat = nan(n_test, n_m, n_ang);

rx_carriers = cell(n_ang, 1);
for a_idx = 1:n_ang
    rx_carriers{a_idx} = carrier_true * exp(1i * deg2rad(dphi_list(a_idx)));
end

for k = 1:n_test
    sid = test_manifest.SampleID(k);
    gt_idx = find(gt_cache.SampleID == sid, 1);
    img_gt = gt_cache.img_gts{gt_idx};

    if streaming
        sig_up = two_dim_upsample_fft(signals_in{k}, rq, aq);
    else
        sig_up = signals_in{k};
    end
    t_opts = struct("nyquist_margin", cfg.nyquist_margin, "carrier", carrier_true);
    U_true = build_threshold(sig_up, sgrid, t_spec, S60, t_opts);
    channel_1bit = quantize_1bit_with_u(sig_up, U_true);
    clear U_true
    if streaming
        clear sig_up
    end

    % B0_par_raw 与 R1_par_raw 与接收载波相偏完全无关，每个样本仅算 1 次后复用至全部角度。
    for m_idx = 1:n_m
        mkey = u_locked.MethodKey(m_idx);
        if ismember(mkey, ["B0_par_raw", "R1_par_raw"])
            ropts = struct("Axis", axis_mode);
            [ch_rec, qdiag] = apply_ta_sft_reconstruct(channel_1bit, [], mkey, ropts);
            if use_fast
                img_out = fast_focus_complex_roi(ch_rec, fast_ctx, rq, aq);
            else
                img_out = focus_upsampled_channel(ch_rec, S60, rq, aq);
            end
            mval = compute_image_metrics(img_out, img_gt);
            psnr_mat(k, m_idx, :) = mval.PSNR;
            ssim_mat(k, m_idx, :) = mval.SSIM;
            zr_mat(k, m_idx, :) = qdiag.ZeroRatio;
            mr_mat(k, m_idx, :) = qdiag.ModifiedRatio;
        end
    end

    for a_idx = 1:n_ang
        c_rx = rx_carriers{a_idx};
        for m_idx = 1:n_m
            mkey = u_locked.MethodKey(m_idx);
            if ismember(mkey, ["B0_par_raw", "R1_par_raw"])
                continue;
            end
            hp_name = u_locked.HyperparamName(m_idx);
            hp_val = double(u_locked.HyperparamValue(m_idx));
            ropts = struct("Axis", axis_mode);
            if hp_name ~= "none" && hp_name ~= "blind"
                ropts.(char(hp_name)) = hp_val;
            end
            [ch_rec, qdiag] = apply_ta_sft_reconstruct(channel_1bit, c_rx, mkey, ropts);
            if use_fast
                img_out = fast_focus_complex_roi(ch_rec, fast_ctx, rq, aq);
            else
                img_out = focus_upsampled_channel(ch_rec, S60, rq, aq);
            end
            mval = compute_image_metrics(img_out, img_gt);
            psnr_mat(k, m_idx, a_idx) = mval.PSNR;
            ssim_mat(k, m_idx, a_idx) = mval.SSIM;
            zr_mat(k, m_idx, a_idx) = qdiag.ZeroRatio;
            mr_mat(k, m_idx, a_idx) = qdiag.ModifiedRatio;
            eg_mat(k, m_idx, a_idx) = qdiag.EstGamma;
        end
    end
end

b0_col = find(method_order == "B0_par_raw", 1);
r1_col = find(method_order == "R1_par_raw", 1);

total_rows = n_ang * n_m;
UnitKey = repmat(string(unit.UnitKey), total_rows, 1);
Branch = repmat(string(unit.Branch), total_rows, 1);
Axis = repmat(axis_mode, total_rows, 1);
GroupName = repmat(string(unit.GroupName), total_rows, 1);
Q = repmat(double(unit.Q), total_rows, 1);
ThresholdType = repmat(ttype, total_rows, 1);
PhaseMismatchDeg = zeros(total_rows, 1);
MethodOrder = zeros(total_rows, 1);
MethodKey = strings(total_rows, 1);
PSNR_Mean = zeros(total_rows, 1);
PSNR_Std = zeros(total_rows, 1);
SSIM_Mean = zeros(total_rows, 1);
SSIM_Std = zeros(total_rows, 1);
DeltaPSNR_vs_B0 = zeros(total_rows, 1);
P_PSNR_vs_B0 = ones(total_rows, 1);
DeltaSSIM_vs_B0 = zeros(total_rows, 1);
P_SSIM_vs_B0 = ones(total_rows, 1);
DeltaPSNR_vs_R1 = zeros(total_rows, 1);
P_PSNR_vs_R1 = ones(total_rows, 1);
DeltaSSIM_vs_R1 = zeros(total_rows, 1);
P_SSIM_vs_R1 = ones(total_rows, 1);
ZeroRatio_Mean = zeros(total_rows, 1);
ModifiedRatio_Mean = zeros(total_rows, 1);
EstGamma_Mean = nan(total_rows, 1);

ptr = 0;
for a_idx = 1:n_ang
    dphi = dphi_list(a_idx);
    b0_p = psnr_mat(:, b0_col, a_idx);
    b0_s = ssim_mat(:, b0_col, a_idx);
    r1_p = psnr_mat(:, r1_col, a_idx);
    r1_s = ssim_mat(:, r1_col, a_idx);
    for m_idx = 1:n_m
        ptr = ptr + 1;
        p_vec = psnr_mat(:, m_idx, a_idx);
        s_vec = ssim_mat(:, m_idx, a_idx);
        PhaseMismatchDeg(ptr) = dphi;
        MethodOrder(ptr) = u_locked.MethodOrder(m_idx);
        MethodKey(ptr) = u_locked.MethodKey(m_idx);
        PSNR_Mean(ptr) = mean(p_vec);
        PSNR_Std(ptr) = std(p_vec);
        SSIM_Mean(ptr) = mean(s_vec);
        SSIM_Std(ptr) = std(s_vec);
        ZeroRatio_Mean(ptr) = mean(zr_mat(:, m_idx, a_idx));
        ModifiedRatio_Mean(ptr) = mean(mr_mat(:, m_idx, a_idx));
        eg_vec = eg_mat(:, m_idx, a_idx);
        if any(~isnan(eg_vec))
            EstGamma_Mean(ptr) = mean(eg_vec(~isnan(eg_vec)));
        end
        if m_idx ~= b0_col
            wp_b0 = compute_paired_wilcoxon(p_vec(:), b0_p(:));
            ws_b0 = compute_paired_wilcoxon(s_vec(:), b0_s(:));
            DeltaPSNR_vs_B0(ptr) = wp_b0.MeanDiff;
            P_PSNR_vs_B0(ptr) = wp_b0.PValue;
            DeltaSSIM_vs_B0(ptr) = ws_b0.MeanDiff;
            P_SSIM_vs_B0(ptr) = ws_b0.PValue;
        end
        if m_idx ~= r1_col
            wp_r1 = compute_paired_wilcoxon(p_vec(:), r1_p(:));
            ws_r1 = compute_paired_wilcoxon(s_vec(:), r1_s(:));
            DeltaPSNR_vs_R1(ptr) = wp_r1.MeanDiff;
            P_PSNR_vs_R1(ptr) = wp_r1.PValue;
            DeltaSSIM_vs_R1(ptr) = ws_r1.MeanDiff;
            P_SSIM_vs_R1(ptr) = ws_r1.PValue;
        end
    end
end

rob_rows = table( ...
    UnitKey, Branch, Axis, GroupName, Q, ThresholdType, PhaseMismatchDeg, ...
    MethodOrder, MethodKey, PSNR_Mean, PSNR_Std, SSIM_Mean, SSIM_Std, ...
    DeltaPSNR_vs_B0, P_PSNR_vs_B0, DeltaSSIM_vs_B0, P_SSIM_vs_B0, ...
    DeltaPSNR_vs_R1, P_PSNR_vs_R1, DeltaSSIM_vs_R1, P_SSIM_vs_R1, ...
    ZeroRatio_Mean, ModifiedRatio_Mean, EstGamma_Mean);
end

function sum_tbl = build_evaluate_summary(detail_table, units, b0_sft_params, locked_params, method_order)
n_u = height(units);
n_m = numel(method_order);
n_rows = n_u * n_m;

UnitKey = strings(n_rows, 1);
Branch = strings(n_rows, 1);
Axis = strings(n_rows, 1);
GroupName = strings(n_rows, 1);
Q = zeros(n_rows, 1);
Range_q = zeros(n_rows, 1);
Azimuth_q = zeros(n_rows, 1);
ThresholdType = strings(n_rows, 1);
B0_ParamTag = strings(n_rows, 1);
B0_STR_dB = zeros(n_rows, 1);
B0_FrOverBr = zeros(n_rows, 1);
B0_FaOverBa = zeros(n_rows, 1);
MethodOrder = zeros(n_rows, 1);
MethodKey = strings(n_rows, 1);
HyperparamName = strings(n_rows, 1);
HyperparamValue = nan(n_rows, 1);
HyperparamTag = strings(n_rows, 1);
TestSampleCount = zeros(n_rows, 1);
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
ModifiedRatio_Mean = zeros(n_rows, 1);
EstGamma_Mean = nan(n_rows, 1);
AbsMeanBias_Mean = zeros(n_rows, 1);

ptr = 0;
for u_idx = 1:n_u
    unit = units(u_idx, :);
    uk = unit.UnitKey;
    b0_row = b0_sft_params(b0_sft_params.UnitKey == uk, :);
    u_locked = order_unit_locked_rows(locked_params(locked_params.UnitKey == uk, :), method_order);
    for m_idx = 1:n_m
        ptr = ptr + 1;
        mkey = method_order(m_idx);
        sub = detail_table(detail_table.UnitKey == uk & detail_table.MethodKey == mkey, :);
        if isempty(sub)
            error("zc1bit:MissingEvaluateRows", "测试汇总缺少单元/方法：%s / %s", uk, mkey);
        end
        lr = u_locked(m_idx, :);
        UnitKey(ptr) = uk;
        Branch(ptr) = string(unit.Branch);
        Axis(ptr) = string(unit.Axis);
        GroupName(ptr) = string(unit.GroupName);
        Q(ptr) = double(unit.Q);
        Range_q(ptr) = double(unit.Range_q);
        Azimuth_q(ptr) = double(unit.Azimuth_q);
        ThresholdType(ptr) = string(unit.ThresholdType);
        B0_ParamTag(ptr) = string(b0_row.ParamTag);
        B0_STR_dB(ptr) = double(b0_row.STR_dB);
        B0_FrOverBr(ptr) = double(b0_row.FrOverBr);
        B0_FaOverBa(ptr) = double(b0_row.FaOverBa);
        MethodOrder(ptr) = double(lr.MethodOrder);
        MethodKey(ptr) = mkey;
        HyperparamName(ptr) = string(lr.HyperparamName);
        HyperparamValue(ptr) = double(lr.HyperparamValue);
        HyperparamTag(ptr) = string(lr.HyperparamTag);
        TestSampleCount(ptr) = height(sub);
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
        ModifiedRatio_Mean(ptr) = mean(sub.ModifiedRatio);
        eg_mask = ~isnan(sub.EstGamma);
        if any(eg_mask)
            EstGamma_Mean(ptr) = mean(sub.EstGamma(eg_mask));
        end
        AbsMeanBias_Mean(ptr) = mean(sub.AbsMeanBias);
    end
end

sum_tbl = table( ...
    UnitKey, Branch, Axis, GroupName, Q, Range_q, Azimuth_q, ThresholdType, ...
    B0_ParamTag, B0_STR_dB, B0_FrOverBr, B0_FaOverBa, ...
    MethodOrder, MethodKey, HyperparamName, HyperparamValue, HyperparamTag, ...
    TestSampleCount, ENLSampleCount, ...
    PSNR_Mean, PSNR_Std, SSIM_Mean, SSIM_Std, Entropy_Mean, ...
    ENL_Mean, GT_ENL_Mean, AbsENLDiff_Mean, ...
    ZeroRatio_Mean, ModifiedRatio_Mean, EstGamma_Mean, AbsMeanBias_Mean);
end

function w_tbl = build_paired_wilcoxon_table(detail_table, units, method_order)
rows = cell(0, 1);
ta_methods = method_order(~ismember(method_order, ["B0_par_raw", "R1_par_raw"]));
non_b0_methods = method_order(method_order ~= "B0_par_raw");

for u_idx = 1:height(units)
    unit = units(u_idx, :);
    uk = unit.UnitKey;
    b0_rows = sortrows(detail_table(detail_table.UnitKey == uk & detail_table.MethodKey == "B0_par_raw", :), "SampleID");
    r1_rows = sortrows(detail_table(detail_table.UnitKey == uk & detail_table.MethodKey == "R1_par_raw", :), "SampleID");

    for m_idx = 1:numel(non_b0_methods)
        mkey = non_b0_methods(m_idx);
        m_order = find(method_order == mkey, 1);
        a_rows = sortrows(detail_table(detail_table.UnitKey == uk & detail_table.MethodKey == mkey, :), "SampleID");
        if ~isempty(a_rows) && isequal(a_rows.SampleID, b0_rows.SampleID)
            rows{end + 1, 1} = make_wilcoxon_row("vs_B0", unit, m_order, mkey, a_rows, b0_rows); %#ok<AGROW>
        end
    end

    for m_idx = 1:numel(ta_methods)
        mkey = ta_methods(m_idx);
        m_order = find(method_order == mkey, 1);
        a_rows = sortrows(detail_table(detail_table.UnitKey == uk & detail_table.MethodKey == mkey, :), "SampleID");
        if ~isempty(a_rows) && isequal(a_rows.SampleID, r1_rows.SampleID)
            rows{end + 1, 1} = make_wilcoxon_row("vs_R1", unit, m_order, mkey, a_rows, r1_rows); %#ok<AGROW>
        end
    end
end

if isempty(rows)
    error("zc1bit:EmptyTaSftWilcoxon", "未生成任何 TA-SFT Wilcoxon 配对检验行。");
end
w_tbl = vertcat(rows{:});
end

function row = make_wilcoxon_row(comp_type, unit, m_order, mkey, a_rows, b_rows)
wp = compute_paired_wilcoxon(a_rows.PSNR(:), b_rows.PSNR(:));
ws = compute_paired_wilcoxon(a_rows.SSIM(:), b_rows.SSIM(:));

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
    string(unit.GroupName), double(unit.Q), string(unit.ThresholdType), ...
    double(m_order), string(mkey), wp.N, n_enl, ...
    wp.MeanDiff, wp.StdDiff, wp.MedianDiff, wp.PValue, wp.SignedRank, ...
    ws.MeanDiff, ws.StdDiff, ws.MedianDiff, ws.PValue, ws.SignedRank, ...
    d_enl_mean, d_enl_std, d_enl_med, p_enl, ...
    'VariableNames', { ...
    'ComparisonType', 'UnitKey', 'Branch', 'Axis', 'GroupName', 'Q', 'ThresholdType', ...
    'MethodOrder', 'MethodKey', 'SampleCount', 'ENLSampleCount', ...
    'DeltaPSNR_Mean', 'DeltaPSNR_Std', 'DeltaPSNR_Median', 'P_PSNR', 'SignedRank_PSNR', ...
    'DeltaSSIM_Mean', 'DeltaSSIM_Std', 'DeltaSSIM_Median', 'P_SSIM', 'SignedRank_SSIM', ...
    'DeltaAbsENLDiff_Mean', 'DeltaAbsENLDiff_Std', 'DeltaAbsENLDiff_Median', 'P_AbsENLDiff'});
end

function tbl = empty_leakage_table()
tbl = table(zeros(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
    zeros(0, 1), strings(0, 1), strings(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    'VariableNames', { ...
    'SampleID', 'Dataset', 'UnitKey', 'GroupName', 'ThresholdType', ...
    'MethodOrder', 'MethodKey', 'HyperparamTag', ...
    'OffSupportRatio', 'RangeOffRatio', 'AzimuthOffRatio'});
end

function tbl = empty_robustness_table()
tbl = table(strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), zeros(0, 1), strings(0, 1), ...
    zeros(0, 1), zeros(0, 1), strings(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    'VariableNames', { ...
    'UnitKey', 'Branch', 'Axis', 'GroupName', 'Q', 'ThresholdType', 'PhaseMismatchDeg', ...
    'MethodOrder', 'MethodKey', 'PSNR_Mean', 'PSNR_Std', 'SSIM_Mean', 'SSIM_Std', ...
    'DeltaPSNR_vs_B0', 'P_PSNR_vs_B0', 'DeltaSSIM_vs_B0', 'P_SSIM_vs_B0', ...
    'DeltaPSNR_vs_R1', 'P_PSNR_vs_R1', 'DeltaSSIM_vs_R1', 'P_SSIM_vs_R1', ...
    'ZeroRatio_Mean', 'ModifiedRatio_Mean', 'EstGamma_Mean'});
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

function save_eval_checkpoint(cp_dir, unit, b0_row, u_locked, test_manifest, S60, cfg, detail_rows, leakage_rows)
cp_path = fullfile(cp_dir, unit.UnitKey + ".mat");
expected_params = build_eval_expected_params(unit, b0_row, u_locked, test_manifest, S60, cfg);
state = struct( ...
    "parameters", expected_params, ...
    "detail_rows", detail_rows, ...
    "leakage_rows", leakage_rows, ...
    "completed", true);
zc1bit_atomic_save(cp_path, struct("state", state));
fprintf("  [评估完成] %s（已原子落盘 %d 条测试评估明细）\n", unit.UnitKey, height(detail_rows));
end

function save_rob_checkpoint(cp_dir, unit, b0_row, u_locked, test_manifest, S60, cfg, rob_rows)
cp_path = fullfile(cp_dir, unit.UnitKey + ".mat");
expected_params = build_rob_expected_params(unit, b0_row, u_locked, test_manifest, S60, cfg);
state = struct( ...
    "parameters", expected_params, ...
    "robustness_rows", rob_rows, ...
    "completed", true);
zc1bit_atomic_save(cp_path, struct("state", state));
fprintf("  [相偏完成] %s（已原子落盘 %d 条相偏鲁棒性记录）\n", unit.UnitKey, height(rob_rows));
end

function expected_params = build_eval_expected_params(unit, b0_row, u_locked, test_manifest, S60, cfg)
expected_params = struct( ...
    "unit", zc1bit_parameters(unit), ...
    "b0_sft_row", zc1bit_parameters(b0_row), ...
    "locked_unit_rows", u_locked(:, ["MethodKey", "HyperparamName", "HyperparamValue", "HyperparamTag"]), ...
    "ta_sft_eval", struct( ...
        "leakage_groups", string(cfg.ta_sft.leakage_groups), ...
        "leakage_threshold_ratio", double(cfg.ta_sft.leakage_threshold_ratio)), ...
    "nyquist_margin", cfg.nyquist_margin, ...
    "seed", cfg.seed, ...
    "manifest", zc1bit_parameters(test_manifest), ...
    "S60", zc1bit_parameters(S60), ...
    "enl_window_size", cfg.enl_window_size, ...
    "enl_stride", cfg.enl_stride, ...
    "entropy_num_bins", cfg.entropy_num_bins);
end

function expected_params = build_rob_expected_params(unit, b0_row, u_locked, test_manifest, S60, cfg)
expected_params = struct( ...
    "unit", zc1bit_parameters(unit), ...
    "b0_sft_row", zc1bit_parameters(b0_row), ...
    "locked_unit_rows", u_locked(:, ["MethodKey", "HyperparamName", "HyperparamValue", "HyperparamTag"]), ...
    "ta_sft_rob", struct("phase_mismatch_deg", double(cfg.ta_sft.phase_mismatch_deg)), ...
    "nyquist_margin", cfg.nyquist_margin, ...
    "seed", cfg.seed, ...
    "manifest", zc1bit_parameters(test_manifest), ...
    "S60", zc1bit_parameters(S60));
end

function sig = build_ta_sft_search_sig(cfg)
ts = cfg.ta_sft;
sig = struct( ...
    "branches", string(ts.branches), ...
    "range_allocations", double(ts.range_allocations), ...
    "azimuth_allocations", double(ts.azimuth_allocations), ...
    "range_threshold_types", string(ts.range_threshold_types), ...
    "azimuth_threshold_types", string(ts.azimuth_threshold_types), ...
    "tau_grid", double(ts.tau_grid), ...
    "alpha_grid", double(ts.alpha_grid), ...
    "gamma_grid", double(ts.gamma_grid));
end

function check_ta_sft_parameters(saved, current, context)
zc1bit_check_parameters(saved, current, context);
for sub_name = ["ta_sft", "ta_sft_eval", "ta_sft_rob"]
    if isfield(saved, sub_name) && isfield(current, sub_name)
        compare_struct_fields_strict(saved.(sub_name), current.(sub_name), context + "." + sub_name);
    end
end
if isfield(saved, "b0_sft_row") && isfield(current, "b0_sft_row")
    zc1bit_check_parameters(saved.b0_sft_row, current.b0_sft_row, context + ".b0_sft_row");
end
if isfield(saved, "locked_unit_rows") && isfield(current, "locked_unit_rows")
    if ~isequaln(saved.locked_unit_rows, current.locked_unit_rows)
        error("zc1bit:ParameterMismatch", "参数 %s.locked_unit_rows 与检查点记录不一致。", context);
    end
end
end

function compare_struct_fields_strict(saved_s, curr_s, context)
fnames = string(fieldnames(curr_s)).';
for fn = fnames
    if ~isfield(saved_s, fn)
        continue;
    end
    if ~isequaln(saved_s.(fn), curr_s.(fn))
        error("zc1bit:ParameterMismatch", "参数 %s.%s 与检查点记录不一致。", context, fn);
    end
end
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
    error("zc1bit:EmptyTaSftEvaluateDetail", "TA-SFT 测试评估明细表为空。");
end
for col = ["PSNR", "SSIM", "Entropy", "ZeroRatio", "ModifiedRatio", "AbsMeanBias"]
    vals = detail_table.(col);
    if any(~isfinite(vals))
        error("zc1bit:NonFiniteTaSftEvaluateDetail", "TA-SFT 测试评估明细列 %s 存在非有限值。", col);
    end
end
end

function write_evaluate_metadata(eval_dir, test_manifest, units, method_order, summary_table, wilcoxon_table, robustness_table, cache_gb)
fid = fopen(fullfile(eval_dir, "ta_sft_evaluate_metadata.txt"), "w", "n", "UTF-8");
if fid < 0
    error("zc1bit:MetadataWrite", "无法写入 ta_sft_evaluate_metadata.txt。");
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "Experiment=SAR-ZC1Bit TA-SFT Evaluate\n");
fprintf(fid, "ExecutionMode=single_core\n");
fprintf(fid, "CacheBudgetGB=%.6g\n", cache_gb);
fprintf(fid, "TestSampleCount=%d\n", height(test_manifest));
fprintf(fid, "UnitCount=%d\n", height(units));
fprintf(fid, "Methods=%s\n", join(method_order, ","));
fprintf(fid, "SummaryRowCount=%d\n", height(summary_table));
fprintf(fid, "WilcoxonRowCount=%d\n", height(wilcoxon_table));
fprintf(fid, "RobustnessRowCount=%d\n", height(robustness_table));
end

function ensure_dir(d)
if ~isfolder(d)
    mkdir(d);
end
end
