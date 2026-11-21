function [img_out, img_gt, metrics, diag] = zc1bit_render_locked_image( ...
    cfg, S60, sample_row, signal60, plan_row, fast_ctx, img_gt, expected_row)
%ZC1BIT_RENDER_LOCKED_IMAGE 按阶段二锁定参数重渲染单样本、单单元、单变体图像。
% plan_row 来自 locked_confirm_plan.csv，必须是一行 table。
% 与 evaluate 共用同一 RT 种子、同一阈值/量化/伪 1.5-bit/快路径成像函数，
% 因此可用于验证绘图与 evaluate_detail 的 PSNR/SSIM 是否一致。

if nargin < 6
    fast_ctx = [];
end
if nargin < 7
    img_gt = [];
end
if nargin < 8
    expected_row = table();
end
if height(plan_row) ~= 1
    error("zc1bit:PlotPlanRow", "plan_row 必须恰好包含一行，实际 %d 行。", height(plan_row));
end

rq = double(plan_row.Range_q);
aq = double(plan_row.Azimuth_q);
grid = resolve_sampling_grid(rq, aq, [S60.nrn, S60.nan], S60);
if isempty(fast_ctx)
    fast_ctx = build_fast_focus_context(S60, rq, aq);
end
if isempty(img_gt)
    img_gt = generate_gt_image(signal60, S60);
end

signal_up = two_dim_upsample_fft(signal60, rq, aq);
ttype = string(plan_row.ThresholdType);
t_spec = struct("ThresholdType", ttype);
t_opts = struct("nyquist_margin", cfg.nyquist_margin);

if ttype == "RT"
    t_spec.As = double(plan_row.As);
    rt_seed = zc1bit_stage2_rt_seed(cfg.seed, string(plan_row.ParamKey), ...
        double(sample_row.SampleID));
    t_opts.rand_stream = RandStream("mt19937ar", "Seed", rt_seed);
elseif ismember(ttype, ["RSFT", "ASFT", "2D-SFT"])
    t_spec.STR_dB = double(plan_row.STR_dB);
    t_spec.fr_over_Br = double(plan_row.FrOverBr);
    t_spec.fa_over_Ba = double(plan_row.FaOverBa);
    t_spec.phi0 = double(plan_row.Phi0);
    cparams = struct("fr_over_Br", t_spec.fr_over_Br, ...
        "fa_over_Ba", t_spec.fa_over_Ba, "phi0", t_spec.phi0);
    t_opts.carrier = build_sft_carrier(grid, cparams, S60);
end

U = build_threshold(signal_up, grid, t_spec, S60, t_opts);
channel_1bit = quantize_1bit_with_u(signal_up, U);
clear U

qopts = struct( ...
    "Axis", string(plan_row.Axis), ...
    "DecisionMode", string(plan_row.DecisionMode), ...
    "RemoveDC", logical(plan_row.RemoveDC));
[channel_p15, qdiag] = apply_zc_pseudo15bit( ...
    channel_1bit, string(plan_row.RuleKey), qopts);
img_out = fast_focus_complex_roi(channel_p15, fast_ctx, rq, aq);
metrics = compute_image_metrics(img_out, img_gt);
metrics.MatchExpected = NaN;
metrics.ExpectedPSNR = NaN;
metrics.ExpectedSSIM = NaN;

if ~isempty(expected_row) && ismember("PSNR", string(expected_row.Properties.VariableNames)) && ...
        ismember("SSIM", string(expected_row.Properties.VariableNames))
    exp_psnr = double(expected_row.PSNR);
    exp_ssim = double(expected_row.SSIM);
    if isinf(metrics.PSNR) && isinf(exp_psnr)
        d_psnr = 0;
    else
        d_psnr = abs(metrics.PSNR - exp_psnr);
    end
    d_ssim = abs(metrics.SSIM - exp_ssim);
    metrics.ExpectedPSNR = exp_psnr;
    metrics.ExpectedSSIM = exp_ssim;
    metrics.MatchExpected = (d_psnr <= 1e-6) && (d_ssim <= 1e-8);
    if ~metrics.MatchExpected
        warning("zc1bit:PlotMetricMismatch", ...
            "重渲染指标与 evaluate_detail 不一致：%s / %s / SampleID=%d；dPSNR=%.3g，dSSIM=%.3g", ...
            char(string(plan_row.UnitKey)), char(string(plan_row.VariantKey)), ...
            double(sample_row.SampleID), d_psnr, d_ssim);
    end
end

diag = struct( ...
    "UnitKey", string(plan_row.UnitKey), ...
    "VariantKey", string(plan_row.VariantKey), ...
    "ThresholdType", ttype, ...
    "Q", double(plan_row.Q), ...
    "Range_q", rq, ...
    "Azimuth_q", aq, ...
    "ZeroRatio", qdiag.ZeroRatio, ...
    "MeanBiasReal", qdiag.MeanBiasReal, ...
    "MeanBiasImag", qdiag.MeanBiasImag, ...
    "AbsMeanBias", qdiag.AbsMeanBias);
end
