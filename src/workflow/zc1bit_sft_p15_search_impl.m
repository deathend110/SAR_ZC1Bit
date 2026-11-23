function out = zc1bit_sft_p15_search_impl(cfg)
%ZC1BIT_SFT_P15_SEARCH_IMPL 严格伪 1.5-bit（{+1,0,-1}）SFT 专题校准寻优与锁参。
% 仅使用 35 个校准样本，前端 SFT 参数固定为阶段二 B0_par_raw 锁定参数，
% 对 P15_1_carrier_agree / P15_2_rlp_she / P15_3_arcsin_lin / P15_4_dual_domain 进行超参寻优，
% 并同步评估 B0_par_raw 与 R1_par_raw 作为校准参照基线。

if nargin == 0
    cfg = zc1bit_config();
end
directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>

record = zc1bit_prepare_run(cfg);
p = zc1bit_paths(cfg);
S60 = load(cfg.parameter_file_60);

units = build_sft_p15_units(cfg, S60);
b0_sft_params = resolve_b0_sft_params(cfg, p, units);
cand_table = build_sft_p15_candidates(cfg);
method_order = [ ...
    "B0_par_raw", "R1_par_raw", ...
    "P15_1_carrier_agree", "P15_2_rlp_she", ...
    "P15_3_arcsin_lin", "P15_4_dual_domain"].';

manifest = record.manifest;
cal_manifest = manifest(manifest.Split == "calibration", :);
if isempty(cal_manifest)
    error("zc1bit:EmptyCalibration", "样本清单中未找到校准样本。");
end

[cache_gb, preload_signals] = resolve_runtime_options(cfg);
cache_budget_bytes = cache_gb * (1024 ^ 3);
use_fast = true;
if isfield(cfg.sft_p15, "use_fast_path")
    use_fast = logical(cfg.sft_p15.use_fast_path);
end

n_units = height(units);
n_cal = height(cal_manifest);
cp_dir = fullfile(p.checkpoints, "sft_p15", "search");
out_dir = fullfile(p.sft_p15, "search");
ensure_dir(cp_dir);
ensure_dir(out_dir);

summary_by_unit = cell(n_units, 1);
locked_by_unit = cell(n_units, 1);
completed_units = false(n_units, 1);

for u_idx = 1:n_units
    unit = units(u_idx, :);
    b0_row = b0_sft_params(b0_sft_params.UnitKey == unit.UnitKey, :);
    cp_path = fullfile(cp_dir, unit.UnitKey + ".mat");
    if ~isfile(cp_path)
        continue;
    end
    state = zc1bit_load_state(cp_path, struct());
    if isfield(state, "completed") && isequal(state.completed, true) && ...
            isfield(state, "summary_rows") && isfield(state, "locked_rows")
        expected_params = build_search_expected_params(unit, b0_row, cal_manifest, S60, cfg);
        check_sft_p15_parameters(state.parameters, expected_params, ...
            sprintf("SFT-P15 搜索检查点(%s)", unit.UnitKey));
        summary_by_unit{u_idx} = state.summary_rows;
        locked_by_unit{u_idx} = state.locked_rows;
        completed_units(u_idx) = true;
    end
end

gt_cache_path = fullfile(cp_dir, "gt_cal_cache.mat");
stage2_gt_path = fullfile(p.checkpoints, "search", "gt_cal_cache.mat");
gt_cache = zc1bit_stage2_gt_cache(cfg, S60, cal_manifest, {}, gt_cache_path, stage2_gt_path);
cal_signals = {};

if any(~completed_units) || isempty(gt_cache)
    cal_signals = zc1bit_stage2_load_signals(cal_manifest, S60, cfg.raw_decimation_step, "校准样本");
    if isempty(gt_cache)
        gt_cache = zc1bit_stage2_gt_cache(cfg, S60, cal_manifest, cal_signals, gt_cache_path, stage2_gt_path);
    end
end

fprintf("\n=== 严格伪 1.5-bit SFT 校准寻优：%d 单元 × %d 候选配置（%d 校准样本 | 已命中 %d/%d | cache=%.1f GB）===\n", ...
    n_units, height(cand_table), n_cal, sum(completed_units), n_units, cache_gb);

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
                [srows, lrows] = run_unit_cached_up( ...
                    cfg, S60, sgrid, units(u_idx, :), b0_sft_params, cand_table, method_order, ...
                    cal_manifest, cal_signals, gt_cache, fast_ctx, use_fast);
                save_search_checkpoint(cp_dir, units(u_idx, :), b0_sft_params, cal_manifest, S60, cfg, srows, lrows);
                summary_by_unit{u_idx} = srows;
                locked_by_unit{u_idx} = lrows;
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
                [srows, lrows] = run_unit_cached_up( ...
                    cfg, S60, sgrid, units(u_idx, :), b0_sft_params, cand_table, method_order, ...
                    cal_manifest, up_signals, gt_cache, fast_ctx, use_fast);
                save_search_checkpoint(cp_dir, units(u_idx, :), b0_sft_params, cal_manifest, S60, cfg, srows, lrows);
                summary_by_unit{u_idx} = srows;
                locked_by_unit{u_idx} = lrows;
            end
            clear up_signals
        else
            fprintf("  [R%dA%d] 流式低内存模式：预计整组上采样 %.2f GB 超出预算，逐样本复用（%d 单元）...\n", ...
                rq, aq, est_up_bytes / (1024 ^ 3), numel(group_u_idx));
            [s_parts, l_parts] = run_alloc_streaming( ...
                cfg, S60, sgrid, units(group_u_idx, :), b0_sft_params, cand_table, method_order, ...
                cal_manifest, cal_signals, gt_cache, fast_ctx, use_fast);
            for m = 1:numel(group_u_idx)
                u_idx = group_u_idx(m);
                save_search_checkpoint(cp_dir, units(u_idx, :), b0_sft_params, cal_manifest, S60, cfg, s_parts{m}, l_parts{m});
                summary_by_unit{u_idx} = s_parts{m};
                locked_by_unit{u_idx} = l_parts{m};
            end
        end
    end
end

summary_table = vertcat(summary_by_unit{:});
locked_params = vertcat(locked_by_unit{:});
validate_search_outputs(summary_table, locked_params, units, method_order);

writetable(summary_table, fullfile(out_dir, "sft_p15_search_summary.csv"));
writetable(locked_params, fullfile(out_dir, "sft_p15_locked_params.csv"));

signature = struct( ...
    "sft_p15", build_sft_p15_cfg_signature(cfg), ...
    "nyquist_margin", cfg.nyquist_margin, ...
    "seed", cfg.seed, ...
    "S60", zc1bit_parameters(S60), ...
    "cal_manifest", zc1bit_parameters(cal_manifest), ...
    "b0_sft_params", zc1bit_parameters(b0_sft_params));

zc1bit_atomic_save(fullfile(out_dir, "sft_p15_search_results.mat"), struct( ...
    "signature", signature, ...
    "units", units, ...
    "b0_sft_params", b0_sft_params, ...
    "cand_table", cand_table, ...
    "method_order", method_order, ...
    "cal_manifest", cal_manifest, ...
    "summary_table", summary_table, ...
    "locked_params", locked_params));
write_search_metadata(out_dir, cfg, cal_manifest, units, method_order, locked_params, cache_gb);

fprintf("严格伪 1.5-bit SFT 校准搜索完成，输出目录：%s\n", out_dir);
out = struct( ...
    "units", units, ...
    "b0_sft_params", b0_sft_params, ...
    "method_order", method_order, ...
    "cal_manifest", cal_manifest, ...
    "summary_table", summary_table, ...
    "locked_params", locked_params, ...
    "output_dir", out_dir);
end

function units = build_sft_p15_units(cfg, S60)
sp = cfg.sft_p15;
branches = string(sp.branches);

est_count = 0;
for b_idx = 1:numel(branches)
    br = branches(b_idx);
    if br == "range_branch"
        est_count = est_count + size(sp.range_allocations, 1) * numel(string(sp.range_threshold_types));
    elseif br == "azimuth_branch"
        est_count = est_count + size(sp.azimuth_allocations, 1) * numel(string(sp.azimuth_threshold_types));
    else
        error("zc1bit:UnknownBranch", "未知 SFT-P15 分支：%s", br);
    end
end
if est_count == 0
    error("zc1bit:EmptySftP15Units", "SFT-P15 未配置任何实验单元。");
end

UnitIdx = (1:est_count).';
Branch = strings(est_count, 1);
Axis = strings(est_count, 1);
GroupName = strings(est_count, 1);
Q = zeros(est_count, 1);
Range_q = zeros(est_count, 1);
Azimuth_q = zeros(est_count, 1);
ThresholdType = strings(est_count, 1);
UnitKey = strings(est_count, 1);

ptr = 0;
for b_idx = 1:numel(branches)
    br = branches(b_idx);
    if br == "range_branch"
        axis_mode = "range";
        allocs = double(sp.range_allocations);
        ttypes = string(sp.range_threshold_types);
        allowed_tt = "RSFT";
        btag = "rb";
        if any(allocs(:, 2) ~= 1)
            error("zc1bit:SftP15AxisMismatch", "range_branch 要求 Azimuth_q == 1。");
        end
    else
        axis_mode = "azimuth";
        allocs = double(sp.azimuth_allocations);
        ttypes = string(sp.azimuth_threshold_types);
        allowed_tt = "ASFT";
        btag = "ab";
        if any(allocs(:, 1) ~= 1)
            error("zc1bit:SftP15AxisMismatch", "azimuth_branch 要求 Range_q == 1。");
        end
    end
    if any(allocs(:, 1) > 1 & allocs(:, 2) > 1)
        error("zc1bit:SftP15BaruNotAllowed", "SFT-P15 仅允许单向上采样，禁止 BARU。");
    end
    if any(~ismember(ttypes, allowed_tt))
        error("zc1bit:SftP15ThresholdNotAllowed", "分支 %s 仅支持 [%s]。", br, allowed_tt);
    end

    for a_idx = 1:size(allocs, 1)
        rq = allocs(a_idx, 1);
        aq = allocs(a_idx, 2);
        resolve_sampling_grid(rq, aq, [S60.nrn, S60.nan], S60);
        gname = sprintf("R%dA%d", rq, aq);
        for t_idx = 1:numel(ttypes)
            ttype = ttypes(t_idx);
            type_key = lower(replace(ttype, "-", ""));
            ptr = ptr + 1;
            Branch(ptr) = br;
            Axis(ptr) = axis_mode;
            GroupName(ptr) = gname;
            Q(ptr) = rq * aq;
            Range_q(ptr) = rq;
            Azimuth_q(ptr) = aq;
            ThresholdType(ptr) = ttype;
            UnitKey(ptr) = sprintf("%s_%s_%s_%s", btag, gname, axis_mode, type_key);
        end
    end
end

units = table(UnitIdx, Branch, Axis, GroupName, Q, Range_q, Azimuth_q, ThresholdType, UnitKey);
if numel(unique(units.UnitKey)) ~= height(units)
    error("zc1bit:SftP15UnitKeyCollision", "SFT-P15 存在重复 UnitKey。");
end
end

function b0_tbl = resolve_b0_sft_params(cfg, p, units)
if isfield(cfg.sft_p15, "b0_locked_table") && istable(cfg.sft_p15.b0_locked_table) && ...
        ~isempty(cfg.sft_p15.b0_locked_table)
    raw_tbl = cfg.sft_p15.b0_locked_table;
else
    src_file = "";
    if isfield(cfg.sft_p15, "baseline_locked_file")
        src_file = string(cfg.sft_p15.baseline_locked_file);
    end
    if strlength(src_file) == 0
        src_file = fullfile(p.locked, "locked_params.csv");
    end
    if ~isfile(src_file)
        error("zc1bit:MissingLockedParams", ...
            "未找到阶段二 B0 锁参表：%s。请先运行阶段二 confirm 或配置 cfg.sft_p15.baseline_locked_file。", src_file);
    end
    raw_tbl = readtable(src_file, "TextType", "string");
end

raw_tbl.UnitKey = string(raw_tbl.UnitKey);
raw_tbl.VariantKey = string(raw_tbl.VariantKey);
raw_tbl.ParamTag = string(raw_tbl.ParamTag);
raw_tbl.ParamKey = string(raw_tbl.ParamKey);

n_u = height(units);
parts = cell(n_u, 1);
for k = 1:n_u
    uk = units.UnitKey(k);
    hit = raw_tbl(raw_tbl.UnitKey == uk & raw_tbl.VariantKey == "B0_par_raw", :);
    if height(hit) ~= 1
        error("zc1bit:MissingB0SftUnit", "B0_par_raw 锁参表中缺失或重复单元：%s（找到 %d 行）。", uk, height(hit));
    end
    if ~ismember("Phi0", string(hit.Properties.VariableNames)) || ~isfinite(double(hit.Phi0))
        phi0_val = 0;
    else
        phi0_val = double(hit.Phi0);
    end
    parts{k} = table( ...
        uk, string(units.ThresholdType(k)), hit.ParamTag(1), hit.ParamKey(1), ...
        double(hit.STR_dB(1)), double(hit.FrOverBr(1)), double(hit.FaOverBa(1)), phi0_val, ...
        'VariableNames', { ...
        'UnitKey', 'ThresholdType', 'ParamTag', 'ParamKey', ...
        'STR_dB', 'FrOverBr', 'FaOverBa', 'Phi0'});
end
b0_tbl = vertcat(parts{:});
end

function cand_tbl = build_sft_p15_candidates(cfg)
sp = cfg.sft_p15;
tau1 = sort(unique(round(double(sp.p1_tau_grid(:)), 6)));
lmin_w1 = double(sp.p2_w1_lmin_grid(:));
lmin_w2 = double(sp.p2_w2_lmin_grid(:));
tau3 = sort(unique(round(double(sp.p3_tau_grid(:)), 6)));
tau4 = sort(unique(round(double(sp.p4_tau_grid(:)), 6)));
lmin4 = double(sp.p4_lmin_grid(:));

rows = cell(0, 1);
rows{end + 1, 1} = make_cand_row(1, "B0_par_raw", "b0", NaN, NaN, NaN, NaN);
rows{end + 1, 1} = make_cand_row(2, "R1_par_raw", "r1", NaN, NaN, NaN, NaN);

for k = 1:numel(tau1)
    v = tau1(k);
    tag = "tau" + format_num_tag(v);
    rows{end + 1, 1} = make_cand_row(3, "P15_1_carrier_agree", tag, v, NaN, NaN, NaN); %#ok<AGROW>
end

for k = 1:numel(lmin_w1)
    lv = lmin_w1(k);
    tag = "w1_lmin" + format_num_tag(lv);
    rows{end + 1, 1} = make_cand_row(4, "P15_2_rlp_she", tag, NaN, 1, lv, NaN); %#ok<AGROW>
end

for k = 1:numel(lmin_w2)
    lv = lmin_w2(k);
    tag = "w2_lmin" + format_num_tag(lv);
    rows{end + 1, 1} = make_cand_row(4, "P15_2_rlp_she", tag, NaN, 2, lv, NaN); %#ok<AGROW>
end

for k = 1:numel(tau3)
    v = tau3(k);
    tag = "taulin" + format_num_tag(v);
    rows{end + 1, 1} = make_cand_row(5, "P15_3_arcsin_lin", tag, NaN, NaN, NaN, v); %#ok<AGROW>
end

for i = 1:numel(tau4)
    tv = tau4(i);
    for j = 1:numel(lmin4)
        lv = lmin4(j);
        tag = sprintf("tau%s_lmin%s", format_num_tag(tv), format_num_tag(lv));
        rows{end + 1, 1} = make_cand_row(6, "P15_4_dual_domain", tag, tv, 1, lv, NaN); %#ok<AGROW>
    end
end

cand_tbl = vertcat(rows{:});
cand_tbl.CandIdx = (1:height(cand_tbl)).';
cand_tbl = cand_tbl(:, [ ...
    "CandIdx", "MethodOrder", "MethodKey", "HyperparamTag", ...
    "Tau", "HalfWidth", "MinRunLength", "TauLin"]);
end

function row = make_cand_row(m_order, m_key, tag, tau, hw, lmin, taulin)
row = table( ...
    double(m_order), string(m_key), string(tag), ...
    double(tau), double(hw), double(lmin), double(taulin), ...
    'VariableNames', { ...
    'MethodOrder', 'MethodKey', 'HyperparamTag', ...
    'Tau', 'HalfWidth', 'MinRunLength', 'TauLin'});
end

function tag = format_num_tag(val)
if isinf(val)
    tag = "Inf";
    return;
end
s = sprintf("%.6g", val);
tag = replace(string(s), ".", "p");
tag = replace(tag, "-", "m");
end

function [summary_rows, locked_rows] = run_unit_cached_up( ...
    cfg, S60, sgrid, unit, b0_sft_params, cand_table, method_order, ...
    cal_manifest, up_signals, gt_cache, fast_ctx, use_fast)

b0_row = b0_sft_params(b0_sft_params.UnitKey == unit.UnitKey, :);
carrier = build_sft_carrier(sgrid, struct( ...
    "fr_over_Br", b0_row.FrOverBr, ...
    "fa_over_Ba", b0_row.FaOverBa, ...
    "phi0", b0_row.Phi0), S60);

n_cal = height(cal_manifest);
parts = cell(n_cal, 1);
for k = 1:n_cal
    parts{k} = evaluate_unit_one_sample( ...
        cfg, S60, sgrid, unit, b0_row, carrier, cand_table, ...
        cal_manifest(k, :), up_signals{k}, gt_cache, fast_ctx, use_fast);
end
details = vertcat(parts{:});
summary_rows = aggregate_unit_search(details, unit, b0_row, cand_table);
locked_rows = select_best_per_method(summary_rows, method_order);
end

function [s_parts, l_parts] = run_alloc_streaming( ...
    cfg, S60, sgrid, group_units, b0_sft_params, cand_table, method_order, ...
    cal_manifest, cal_signals, gt_cache, fast_ctx, use_fast)

n_g = height(group_units);
n_cal = height(cal_manifest);
rq = double(group_units.Range_q(1));
aq = double(group_units.Azimuth_q(1));

b0_rows = cell(n_g, 1);
carriers = cell(n_g, 1);
for m = 1:n_g
    uk = group_units.UnitKey(m);
    b0_r = b0_sft_params(b0_sft_params.UnitKey == uk, :);
    b0_rows{m} = b0_r;
    carriers{m} = build_sft_carrier(sgrid, struct( ...
        "fr_over_Br", b0_r.FrOverBr, ...
        "fa_over_Ba", b0_r.FaOverBa, ...
        "phi0", b0_r.Phi0), S60);
end

det_by_u = cell(n_g, n_cal);
for k = 1:n_cal
    sig_up = two_dim_upsample_fft(cal_signals{k}, rq, aq);
    srow = cal_manifest(k, :);
    for m = 1:n_g
        det_by_u{m, k} = evaluate_unit_one_sample( ...
            cfg, S60, sgrid, group_units(m, :), b0_rows{m}, carriers{m}, cand_table, ...
            srow, sig_up, gt_cache, fast_ctx, use_fast);
    end
    clear sig_up
end

s_parts = cell(n_g, 1);
l_parts = cell(n_g, 1);
for m = 1:n_g
    details_u = vertcat(det_by_u{m, :});
    s_parts{m} = aggregate_unit_search(details_u, group_units(m, :), b0_rows{m}, cand_table);
    l_parts{m} = select_best_per_method(s_parts{m}, method_order);
end
end

function rows = evaluate_unit_one_sample( ...
    cfg, S60, sgrid, unit, b0_row, carrier, cand_table, ...
    sample_row, sig_up, gt_cache, fast_ctx, use_fast)

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

t_spec = struct( ...
    "ThresholdType", string(unit.ThresholdType), ...
    "STR_dB", double(b0_row.STR_dB), ...
    "fr_over_Br", double(b0_row.FrOverBr), ...
    "fa_over_Ba", double(b0_row.FaOverBa), ...
    "phi0", double(b0_row.Phi0));
t_opts = struct("nyquist_margin", cfg.nyquist_margin, "carrier", carrier);
U = build_threshold(sig_up, sgrid, t_spec, S60, t_opts);
c1 = quantize_1bit_with_u(sig_up, U);
clear U

n_c = height(cand_table);
psnr_vec = zeros(n_c, 1);
ssim_vec = zeros(n_c, 1);
ent_vec = zeros(n_c, 1);
enl_vec = nan(n_c, 1);
abs_enl_vec = nan(n_c, 1);
zr_vec = zeros(n_c, 1);
prot_vec = zeros(n_c, 1);
abm_vec = zeros(n_c, 1);

ent_bins = cfg.entropy_num_bins;
enl_win = cfg.enl_window_size;

for c_idx = 1:n_c
    mkey = cand_table.MethodKey(c_idx);
    ropts = build_reconstruct_options(axis_mode, cand_table(c_idx, :));
    [ch_out, qdiag] = apply_sft_p15_reconstruct(c1, carrier, mkey, ropts);
    if use_fast
        img_out = fast_focus_complex_roi(ch_out, fast_ctx, rq, aq);
    else
        img_out = focus_upsampled_channel(ch_out, S60, rq, aq);
    end
    m = compute_image_metrics(img_out, img_gt);
    psnr_vec(c_idx) = m.PSNR;
    ssim_vec(c_idx) = m.SSIM;
    ent_vec(c_idx) = compute_entropy(img_out, ent_bins);
    if is_enl
        enl_vec(c_idx) = compute_enl(img_out, top_r, left_c, enl_win);
        abs_enl_vec(c_idx) = abs(enl_vec(c_idx) - gt_enl_val);
    end
    zr_vec(c_idx) = qdiag.ZeroRatio;
    prot_vec(c_idx) = qdiag.ProtectedEdgeRatio;
    abm_vec(c_idx) = qdiag.AbsMeanBias;
end

rows = table( ...
    repmat(sid, n_c, 1), ...
    cand_table.CandIdx, ...
    psnr_vec, ssim_vec, ent_vec, enl_vec, repmat(gt_enl_val, n_c, 1), abs_enl_vec, ...
    zr_vec, prot_vec, abm_vec, ...
    'VariableNames', { ...
    'SampleID', 'CandIdx', ...
    'PSNR', 'SSIM', 'Entropy', 'ENL', 'GT_ENL', 'AbsENLDiff', ...
    'ZeroRatio', 'ProtectedEdgeRatio', 'AbsMeanBias'});
end

function ropts = build_reconstruct_options(axis_mode, crow)
ropts = struct("Axis", axis_mode);
if ~isnan(crow.Tau)
    ropts.Tau = double(crow.Tau);
end
if ~isnan(crow.HalfWidth)
    ropts.HalfWidth = double(crow.HalfWidth);
end
if ~isnan(crow.MinRunLength)
    ropts.MinRunLength = double(crow.MinRunLength);
end
if ~isnan(crow.TauLin)
    ropts.TauLin = double(crow.TauLin);
end
end

function sum_tbl = aggregate_unit_search(details, unit, b0_row, cand_table)
n_c = height(cand_table);

UnitKey = repmat(string(unit.UnitKey), n_c, 1);
Branch = repmat(string(unit.Branch), n_c, 1);
Axis = repmat(string(unit.Axis), n_c, 1);
GroupName = repmat(string(unit.GroupName), n_c, 1);
Q = repmat(double(unit.Q), n_c, 1);
Range_q = repmat(double(unit.Range_q), n_c, 1);
Azimuth_q = repmat(double(unit.Azimuth_q), n_c, 1);
ThresholdType = repmat(string(unit.ThresholdType), n_c, 1);

B0_ParamTag = repmat(string(b0_row.ParamTag), n_c, 1);
B0_STR_dB = repmat(double(b0_row.STR_dB), n_c, 1);
B0_FrOverBr = repmat(double(b0_row.FrOverBr), n_c, 1);
B0_FaOverBa = repmat(double(b0_row.FaOverBa), n_c, 1);
B0_Phi0 = repmat(double(b0_row.Phi0), n_c, 1);

MethodOrder = cand_table.MethodOrder;
MethodKey = cand_table.MethodKey;
HyperparamTag = cand_table.HyperparamTag;
Tau = cand_table.Tau;
HalfWidth = cand_table.HalfWidth;
MinRunLength = cand_table.MinRunLength;
TauLin = cand_table.TauLin;

SampleCount = zeros(n_c, 1);
ENLSampleCount = zeros(n_c, 1);
Calibration_SSIM_Mean = zeros(n_c, 1);
Calibration_SSIM_Std = zeros(n_c, 1);
Calibration_PSNR_Mean = zeros(n_c, 1);
Calibration_PSNR_Std = zeros(n_c, 1);
Calibration_Entropy_Mean = zeros(n_c, 1);
Calibration_ENL_Mean = nan(n_c, 1);
Calibration_GT_ENL_Mean = nan(n_c, 1);
Calibration_AbsENLDiff_Mean = nan(n_c, 1);
Calibration_ZeroRatio_Mean = zeros(n_c, 1);
Calibration_ProtectedEdgeRatio_Mean = zeros(n_c, 1);
Calibration_AbsMeanBias_Mean = zeros(n_c, 1);

for c_idx = 1:n_c
    sub = details(details.CandIdx == cand_table.CandIdx(c_idx), :);
    SampleCount(c_idx) = height(sub);
    enl_mask = ~isnan(sub.ENL);
    ENLSampleCount(c_idx) = sum(enl_mask);
    Calibration_SSIM_Mean(c_idx) = mean(sub.SSIM);
    Calibration_SSIM_Std(c_idx) = std(sub.SSIM);
    Calibration_PSNR_Mean(c_idx) = mean(sub.PSNR);
    Calibration_PSNR_Std(c_idx) = std(sub.PSNR);
    Calibration_Entropy_Mean(c_idx) = mean(sub.Entropy);
    if any(enl_mask)
        Calibration_ENL_Mean(c_idx) = mean(sub.ENL(enl_mask));
        Calibration_GT_ENL_Mean(c_idx) = mean(sub.GT_ENL(enl_mask));
        Calibration_AbsENLDiff_Mean(c_idx) = mean(sub.AbsENLDiff(enl_mask));
    end
    Calibration_ZeroRatio_Mean(c_idx) = mean(sub.ZeroRatio);
    Calibration_ProtectedEdgeRatio_Mean(c_idx) = mean(sub.ProtectedEdgeRatio);
    Calibration_AbsMeanBias_Mean(c_idx) = mean(sub.AbsMeanBias);
end

sum_tbl = table( ...
    UnitKey, Branch, Axis, GroupName, Q, Range_q, Azimuth_q, ThresholdType, ...
    B0_ParamTag, B0_STR_dB, B0_FrOverBr, B0_FaOverBa, B0_Phi0, ...
    MethodOrder, MethodKey, HyperparamTag, Tau, HalfWidth, MinRunLength, TauLin, ...
    SampleCount, ENLSampleCount, ...
    Calibration_SSIM_Mean, Calibration_SSIM_Std, Calibration_PSNR_Mean, Calibration_PSNR_Std, ...
    Calibration_Entropy_Mean, Calibration_ENL_Mean, Calibration_GT_ENL_Mean, ...
    Calibration_AbsENLDiff_Mean, Calibration_ZeroRatio_Mean, Calibration_ProtectedEdgeRatio_Mean, ...
    Calibration_AbsMeanBias_Mean);
end

function locked_tbl = select_best_per_method(sum_tbl, method_order)
n_m = numel(method_order);
parts = cell(n_m, 1);
for m_idx = 1:n_m
    mkey = method_order(m_idx);
    sub = sum_tbl(sum_tbl.MethodKey == mkey, :);
    if isempty(sub)
        error("zc1bit:MissingSftP15MethodSummary", "未找到方法 %s 的寻优汇总行。", mkey);
    end
    sub = sortrows(sub, ...
        {'Calibration_SSIM_Mean', 'Calibration_PSNR_Mean', 'Calibration_AbsENLDiff_Mean', 'Calibration_ZeroRatio_Mean', 'HyperparamTag'}, ...
        {'descend', 'descend', 'ascend', 'ascend', 'ascend'});
    parts{m_idx} = sub(1, :);
end
locked_tbl = vertcat(parts{:});
end

function save_search_checkpoint(cp_dir, unit, b0_sft_params, cal_manifest, S60, cfg, summary_rows, locked_rows)
b0_row = b0_sft_params(b0_sft_params.UnitKey == unit.UnitKey, :);
cp_path = fullfile(cp_dir, unit.UnitKey + ".mat");
expected_params = build_search_expected_params(unit, b0_row, cal_manifest, S60, cfg);
state = struct( ...
    "parameters", expected_params, ...
    "summary_rows", summary_rows, ...
    "locked_rows", locked_rows, ...
    "completed", true);
zc1bit_atomic_save(cp_path, struct("state", state));
fprintf("  [完成] %s（已原子落盘 %d 条校准寻优汇总 / %d 条方法锁参）\n", ...
    unit.UnitKey, height(summary_rows), height(locked_rows));
end

function expected_params = build_search_expected_params(unit, b0_row, cal_manifest, S60, cfg)
expected_params = struct( ...
    "unit", zc1bit_parameters(unit), ...
    "b0_sft_row", zc1bit_parameters(b0_row), ...
    "sft_p15", build_sft_p15_cfg_signature(cfg), ...
    "nyquist_margin", cfg.nyquist_margin, ...
    "seed", cfg.seed, ...
    "manifest", zc1bit_parameters(cal_manifest), ...
    "S60", zc1bit_parameters(S60), ...
    "enl_window_size", cfg.enl_window_size, ...
    "enl_stride", cfg.enl_stride, ...
    "entropy_num_bins", cfg.entropy_num_bins);
end

function sig = build_sft_p15_cfg_signature(cfg)
sp = cfg.sft_p15;
sig = struct( ...
    "branches", string(sp.branches), ...
    "range_allocations", double(sp.range_allocations), ...
    "azimuth_allocations", double(sp.azimuth_allocations), ...
    "range_threshold_types", string(sp.range_threshold_types), ...
    "azimuth_threshold_types", string(sp.azimuth_threshold_types), ...
    "p1_tau_grid", double(sp.p1_tau_grid), ...
    "p2_w1_lmin_grid", double(sp.p2_w1_lmin_grid), ...
    "p2_w2_lmin_grid", double(sp.p2_w2_lmin_grid), ...
    "p3_tau_grid", double(sp.p3_tau_grid), ...
    "p4_tau_grid", double(sp.p4_tau_grid), ...
    "p4_lmin_grid", double(sp.p4_lmin_grid));
end

function check_sft_p15_parameters(saved, current, context)
zc1bit_check_parameters(saved, current, context);
if isfield(saved, "sft_p15") && isfield(current, "sft_p15")
    compare_struct_fields_strict(saved.sft_p15, current.sft_p15, context + ".sft_p15");
end
if isfield(saved, "b0_sft_row") && isfield(current, "b0_sft_row")
    zc1bit_check_parameters(saved.b0_sft_row, current.b0_sft_row, context + ".b0_sft_row");
end
end

function compare_struct_fields_strict(saved_s, curr_s, context)
fnames = string(fieldnames(curr_s)).';
for fn = fnames
    if ~isfield(saved_s, fn)
        continue;
    end
    a = saved_s.(fn);
    b = curr_s.(fn);
    if ~isequaln(a, b)
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

function validate_search_outputs(summary_table, locked_params, units, method_order)
if isempty(summary_table) || height(summary_table) == 0
    error("zc1bit:EmptySftP15SearchSummary", "SFT-P15 搜索汇总表为空。");
end
expected_locked = height(units) * numel(method_order);
if height(locked_params) ~= expected_locked
    error("zc1bit:IncompleteSftP15Locked", ...
        "SFT-P15 锁参行数为 %d，与预期 %d 不一致。", height(locked_params), expected_locked);
end
for col = ["Calibration_SSIM_Mean", "Calibration_PSNR_Mean", "Calibration_Entropy_Mean", ...
        "Calibration_ZeroRatio_Mean", "Calibration_ProtectedEdgeRatio_Mean", "Calibration_AbsMeanBias_Mean"]
    if any(~isfinite(summary_table.(col)))
        error("zc1bit:NonFiniteSftP15Search", "搜索汇总列 %s 存在非有限值。", col);
    end
end
end

function write_search_metadata(out_dir, cfg, cal_manifest, units, method_order, locked_params, cache_gb)
fid = fopen(fullfile(out_dir, "sft_p15_search_metadata.txt"), "w", "n", "UTF-8");
if fid < 0
    error("zc1bit:MetadataWrite", "无法写入 sft_p15_search_metadata.txt。");
end
cleanup = onCleanup(@() fclose(fid));
sp = cfg.sft_p15;
fprintf(fid, "Experiment=SAR-ZC1Bit Strict Pseudo-1.5-Bit SFT Search\n");
fprintf(fid, "ExecutionMode=single_core\n");
fprintf(fid, "CacheBudgetGB=%.6g\n", cache_gb);
fprintf(fid, "CalibrationSampleCount=%d\n", height(cal_manifest));
fprintf(fid, "UnitCount=%d\n", height(units));
fprintf(fid, "Methods=%s\n", join(method_order, ","));
fprintf(fid, "LockedRowCount=%d\n", height(locked_params));
fprintf(fid, "P1TauGrid=%s\n", mat2str(sp.p1_tau_grid));
fprintf(fid, "P2W1LminGrid=%s\n", mat2str(sp.p2_w1_lmin_grid));
fprintf(fid, "P2W2LminGrid=%s\n", mat2str(sp.p2_w2_lmin_grid));
fprintf(fid, "P3TauGrid=%s\n", mat2str(sp.p3_tau_grid));
fprintf(fid, "P4TauGrid=%s\n", mat2str(sp.p4_tau_grid));
fprintf(fid, "P4LminGrid=%s\n", mat2str(sp.p4_lmin_grid));
end

function ensure_dir(d)
if ~isfolder(d)
    mkdir(d);
end
end
