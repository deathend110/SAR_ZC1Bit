function rows = zc1bit_stage2_eval_threshold( ...
    cfg, S60, sgrid, threshold_spec, sample_row, signal_up, variant_tbl, axis_mode, ...
    gt_cache, carrier, fast_ctx, use_fast, seed_tag)
%ZC1BIT_STAGE2_EVAL_THRESHOLD 在单样本与单个阈值参数下构造一次 1-bit 通道并批量评估多变体。

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

ttype = string(threshold_spec.ThresholdType);
t_opts = struct("nyquist_margin", cfg.nyquist_margin, "carrier", carrier);
if ttype == "RT"
    rt_seed = zc1bit_stage2_rt_seed(cfg.seed, seed_tag, sid);
    t_opts.rand_stream = RandStream("mt19937ar", "Seed", rt_seed);
end

U = build_threshold(signal_up, sgrid, threshold_spec, S60, t_opts);
channel_1bit = quantize_1bit_with_u(signal_up, U);
clear U

n_var = height(variant_tbl);
psnr_vec = zeros(n_var, 1);
ssim_vec = zeros(n_var, 1);
ent_vec = zeros(n_var, 1);
enl_vec = nan(n_var, 1);
abs_enl_diff_vec = nan(n_var, 1);
zr_vec = zeros(n_var, 1);
mbr_vec = zeros(n_var, 1);
mbi_vec = zeros(n_var, 1);
abm_vec = zeros(n_var, 1);

ent_bins = cfg.entropy_num_bins;
enl_win = cfg.enl_window_size;

for v_idx = 1:n_var
    qopts = struct( ...
        "Axis", axis_mode, ...
        "DecisionMode", variant_tbl.DecisionMode(v_idx), ...
        "RemoveDC", logical(variant_tbl.RemoveDC(v_idx)));
    [channel_p15, qdiag] = apply_zc_pseudo15bit(channel_1bit, variant_tbl.RuleKey(v_idx), qopts);
    if use_fast
        img_out = fast_focus_complex_roi(channel_p15, fast_ctx, rq, aq);
    else
        img_out = focus_upsampled_channel(channel_p15, S60, rq, aq);
    end
    m = compute_image_metrics(img_out, img_gt);
    psnr_vec(v_idx) = m.PSNR;
    ssim_vec(v_idx) = m.SSIM;
    ent_vec(v_idx) = compute_entropy(img_out, ent_bins);
    if is_enl
        enl_vec(v_idx) = compute_enl(img_out, top_r, left_c, enl_win);
        abs_enl_diff_vec(v_idx) = abs(enl_vec(v_idx) - gt_enl_val);
    end
    zr_vec(v_idx) = qdiag.ZeroRatio;
    mbr_vec(v_idx) = qdiag.MeanBiasReal;
    mbi_vec(v_idx) = qdiag.MeanBiasImag;
    abm_vec(v_idx) = qdiag.AbsMeanBias;
end

rows = table( ...
    repmat(sid, n_var, 1), ...
    repmat(sample_row.LocalSampleIdx, n_var, 1), ...
    repmat(sample_row.Dataset, n_var, 1), ...
    variant_tbl.VariantKey, ...
    variant_tbl.RuleKey, ...
    variant_tbl.DecisionMode, ...
    logical(variant_tbl.RemoveDC), ...
    psnr_vec, ssim_vec, ent_vec, enl_vec, repmat(gt_enl_val, n_var, 1), abs_enl_diff_vec, ...
    zr_vec, mbr_vec, mbi_vec, abm_vec, ...
    'VariableNames', { ...
    'SampleID', 'LocalSampleIdx', 'Dataset', ...
    'VariantKey', 'RuleKey', 'DecisionMode', 'RemoveDC', ...
    'PSNR', 'SSIM', 'Entropy', 'ENL', 'GT_ENL', 'AbsENLDiff', ...
    'ZeroRatio', 'MeanBiasReal', 'MeanBiasImag', 'AbsMeanBias'});
end
