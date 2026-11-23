function [channel_out, diag] = apply_ta_sft_reconstruct(channel_1bit, carrier, method_key, options)
%APPLY_TA_SFT_RECONSTRUCT 阈值感知 SFT 1-bit 接收端重构（TA-SFT）。
% 严格仅接受标准 1-bit（{+1, -1}）量化复矩阵与接收端已知单位模 SFT 载波，
% 物理隔绝未量化信号与 sigma_hat，I/Q 各自独立处理且矩阵尺寸不变。
%
% 支持 method_key：
%   "B0_par_raw"     - 原始 1-bit（±1 不变）
%   "R1_par_raw"     - 盲双沿双边并行置零（伪 1.5-bit 对照）
%   "TA1_gate_zc"    - 方案 1：阈值门控选择性置零（严格三电平 {+1, 0, -1}，超参 Tau）
%   "TA2_anchor_zc"  - 方案 2（连续）：过零锚点 -alpha*u_mid 替换（超参 Alpha）
%   "TA2_anchor_q5"  - 方案 2（5 电平）：先节点平均再对称量化至 {-1,-0.5,0,+0.5,+1}（超参 Alpha）
%   "TA3_mmse_fixed" - 方案 3A：标称阈信比逐点截断高斯 MMSE（超参 Gamma）
%   "TA3_mmse_est"   - 方案 3B：纯 1-bit 相关盲估计 gamma_hat + 逐点 MMSE（0 超参）

if nargin < 4 || isempty(options)
    options = struct();
end
axis_mode = "range";
if isfield(options, "Axis") && ~isempty(options.Axis)
    axis_mode = string(options.Axis);
end
if ~isscalar(axis_mode) || ~ismember(axis_mode, ["range", "azimuth"])
    error("zc1bit:AxisMode", "TA-SFT 仅支持 range 或 azimuth 单向轴，收到：%s", axis_mode);
end

method_key = string(method_key);
allowed_methods = [ ...
    "B0_par_raw", "R1_par_raw", ...
    "TA1_gate_zc", "TA2_anchor_zc", "TA2_anchor_q5", ...
    "TA3_mmse_fixed", "TA3_mmse_est"];
if ~isscalar(method_key) || ~ismember(method_key, allowed_methods)
    error("zc1bit:UnknownTaSftMethod", "未知 TA-SFT 重构方法：%s", method_key);
end

if isempty(channel_1bit) || ~isnumeric(channel_1bit)
    error("zc1bit:Invalid1BitInput", "输入必须为非空数值复矩阵。");
end
re_in = real(channel_1bit);
im_in = imag(channel_1bit);
if ~all(re_in == 1 | re_in == -1, "all") || ~all(im_in == 1 | im_in == -1, "all")
    error("zc1bit:Invalid1BitInput", "TA-SFT 输入实部与虚部必须严格属于 {+1, -1}。");
end

requires_carrier = ~ismember(method_key, ["B0_par_raw", "R1_par_raw"]);
if requires_carrier
    if isempty(carrier) || ~isnumeric(carrier) || ~isequal(size(carrier), size(channel_1bit)) || ...
            any(~isfinite(real(carrier)), "all") || any(~isfinite(imag(carrier)), "all")
        error("zc1bit:InvalidCarrier", "方法 %s 需要与输入同尺寸且全有限的复载波。", method_key);
    end
elseif ~isempty(carrier)
    if ~isnumeric(carrier) || ~isequal(size(carrier), size(channel_1bit)) || ...
            any(~isfinite(real(carrier)), "all") || any(~isfinite(imag(carrier)), "all")
        error("zc1bit:InvalidCarrier", "传入的复载波尺寸不匹配或包含非有限值。");
    end
end

param_val = NaN;
est_gamma = NaN;

switch method_key
    case "B0_par_raw"
        channel_out = channel_1bit;

    case "R1_par_raw"
        qopts = struct("Axis", axis_mode, "DecisionMode", "parallel", "RemoveDC", false);
        channel_out = apply_zc_pseudo15bit(channel_1bit, "R1", qopts);

    case "TA1_gate_zc"
        tau = resolve_nonneg_scalar(options, "Tau", 0.5);
        param_val = tau;
        u_re = double(real(carrier));
        u_im = double(imag(carrier));
        re_out = apply_ta1_gate(re_in, u_re, tau, axis_mode);
        im_out = apply_ta1_gate(im_in, u_im, tau, axis_mode);
        channel_out = complex(re_out, im_out);

    case {"TA2_anchor_zc", "TA2_anchor_q5"}
        alpha = resolve_nonneg_scalar(options, "Alpha", 1.0);
        param_val = alpha;
        quantize_q5 = (method_key == "TA2_anchor_q5");
        u_re = double(real(carrier));
        u_im = double(imag(carrier));
        re_out = apply_ta2_anchor(re_in, u_re, alpha, axis_mode, quantize_q5);
        im_out = apply_ta2_anchor(im_in, u_im, alpha, axis_mode, quantize_q5);
        channel_out = complex(re_out, im_out);

    case "TA3_mmse_fixed"
        gamma = resolve_nonneg_scalar(options, "Gamma", 1.0);
        param_val = gamma;
        u_re = double(real(carrier));
        u_im = double(imag(carrier));
        re_out = apply_ta3_mmse(re_in, u_re, gamma);
        im_out = apply_ta3_mmse(im_in, u_im, gamma);
        channel_out = complex(re_out, im_out);

    case "TA3_mmse_est"
        u_re = double(real(carrier));
        u_im = double(imag(carrier));
        est_gamma = estimate_gamma_from_1bit(re_in, im_in, u_re, u_im);
        param_val = est_gamma;
        re_out = apply_ta3_mmse(re_in, u_re, est_gamma);
        im_out = apply_ta3_mmse(im_in, u_im, est_gamma);
        channel_out = complex(re_out, im_out);
end

re_final = real(channel_out);
im_final = imag(channel_out);
zero_ratio_real = mean(re_final(:) == 0);
zero_ratio_imag = mean(im_final(:) == 0);
mod_ratio_real = mean(re_final(:) ~= re_in(:));
mod_ratio_imag = mean(im_final(:) ~= im_in(:));
mean_bias_real = mean(double(re_final(:)));
mean_bias_imag = mean(double(im_final(:)));

diag = struct( ...
    "MethodKey", method_key, ...
    "Axis", axis_mode, ...
    "ParamValue", param_val, ...
    "EstGamma", est_gamma, ...
    "ZeroRatioReal", zero_ratio_real, ...
    "ZeroRatioImag", zero_ratio_imag, ...
    "ZeroRatio", 0.5 * (zero_ratio_real + zero_ratio_imag), ...
    "ModifiedRatio", 0.5 * (mod_ratio_real + mod_ratio_imag), ...
    "MeanBiasReal", mean_bias_real, ...
    "MeanBiasImag", mean_bias_imag, ...
    "AbsMeanBias", hypot(mean_bias_real, mean_bias_imag));
end

function val = resolve_nonneg_scalar(options, field_name, default_val)
val = default_val;
if isfield(options, field_name) && ~isempty(options.(field_name))
    val = double(options.(field_name));
end
if ~isscalar(val) || ~isnumeric(val) || ~isfinite(val) || val < 0
    error("zc1bit:InvalidHyperparam", "参数 %s 必须为非负有限标量。", field_name);
end
end

function y = apply_ta1_gate(x, u, tau, axis_mode)
y = x;
if tau == 0
    return;
end
[nr, nc] = size(x);
if axis_mode == "range"
    if nr < 2 || nc == 0
        return;
    end
    zc_edge = (x(1:end - 1, :) ~= x(2:end, :));
    u_mid = 0.5 * (u(1:end - 1, :) + u(2:end, :));
    gate = zc_edge & (abs(u_mid) <= tau);
    if ~any(gate, "all")
        return;
    end
    mask = false(nr, nc);
    mask(1:end - 1, :) = mask(1:end - 1, :) | gate;
    mask(2:end, :) = mask(2:end, :) | gate;
    y(mask) = 0;
else
    if nc < 2 || nr == 0
        return;
    end
    zc_edge = (x(:, 1:end - 1) ~= x(:, 2:end));
    u_mid = 0.5 * (u(:, 1:end - 1) + u(:, 2:end));
    gate = zc_edge & (abs(u_mid) <= tau);
    if ~any(gate, "all")
        return;
    end
    mask = false(nr, nc);
    mask(:, 1:end - 1) = mask(:, 1:end - 1) | gate;
    mask(:, 2:end) = mask(:, 2:end) | gate;
    y(mask) = 0;
end
end

function y = apply_ta2_anchor(x, u, alpha, axis_mode, quantize_q5)
[nr, nc] = size(x);
if (axis_mode == "range" && (nr < 2 || nc == 0)) || ...
        (axis_mode == "azimuth" && (nc < 2 || nr == 0))
    y = x;
    return;
end

val_sum = zeros(nr, nc);
hit_cnt = zeros(nr, nc);

if axis_mode == "range"
    zc_edge = (x(1:end - 1, :) ~= x(2:end, :));
    if any(zc_edge, "all")
        u_mid = 0.5 * (u(1:end - 1, :) + u(2:end, :));
        v_edge = zeros(nr - 1, nc);
        v_edge(zc_edge) = -alpha * u_mid(zc_edge);
        val_sum(1:end - 1, :) = val_sum(1:end - 1, :) + v_edge;
        val_sum(2:end, :) = val_sum(2:end, :) + v_edge;
        cnt_edge = double(zc_edge);
        hit_cnt(1:end - 1, :) = hit_cnt(1:end - 1, :) + cnt_edge;
        hit_cnt(2:end, :) = hit_cnt(2:end, :) + cnt_edge;
    end
else
    zc_edge = (x(:, 1:end - 1) ~= x(:, 2:end));
    if any(zc_edge, "all")
        u_mid = 0.5 * (u(:, 1:end - 1) + u(:, 2:end));
        v_edge = zeros(nr, nc - 1);
        v_edge(zc_edge) = -alpha * u_mid(zc_edge);
        val_sum(:, 1:end - 1) = val_sum(:, 1:end - 1) + v_edge;
        val_sum(:, 2:end) = val_sum(:, 2:end) + v_edge;
        cnt_edge = double(zc_edge);
        hit_cnt(:, 1:end - 1) = hit_cnt(:, 1:end - 1) + cnt_edge;
        hit_cnt(:, 2:end) = hit_cnt(:, 2:end) + cnt_edge;
    end
end

mask = hit_cnt > 0;
if ~any(mask, "all")
    y = x;
    return;
end

y_dbl = double(x);
node_val = val_sum(mask) ./ hit_cnt(mask);
if quantize_q5
    % 必须先节点平均再统一对称 5 电平量化，确保元素严格属于 {-1, -0.5, 0, +0.5, +1}。
    node_val = round(2 * max(-1, min(1, node_val))) / 2;
end
y_dbl(mask) = node_val;
y = cast(y_dbl, "like", x);
end

function y = apply_ta3_mmse(x, u, gamma)
if gamma == 0
    y = x;
    return;
end
x_dbl = double(x);
z = x_dbl .* (gamma * u);
imr = inv_mills_ratio(z);
y = cast(x_dbl .* imr, "like", x);
end

function out = inv_mills_ratio(z)
%INV_MILLS_RATIO 计算标准正态逆 Mills 比 phi(z) / Phi(z)，全实轴数值稳定。
out = zeros(size(z));
pos = z >= 0;
if any(pos, "all")
    zp = z(pos);
    out(pos) = sqrt(2 / pi) * exp(-0.5 * (zp .^ 2)) ./ erfc(-zp / sqrt(2));
end
neg = ~pos;
if any(neg, "all")
    zn = z(neg);
    out(neg) = sqrt(2 / pi) ./ erfcx(-zn / sqrt(2));
end
if any(~isfinite(out), "all")
    error("zc1bit:NonFiniteMillsRatio", "逆 Mills 比计算出现非有限值。");
end
end

function gamma_hat = estimate_gamma_from_1bit(re_in, im_in, u_re, u_im)
%ESTIMATE_GAMMA_FROM_1BIT 仅利用 1-bit 数据与单位模载波内积反演全局阈信比 gamma_hat。
persistent gamma_lut g_lut
if isempty(gamma_lut)
    gamma_lut = [linspace(0, 4, 801), linspace(4.02, 20, 800)].';
    theta = (0:1023) * (2 * pi / 1024);
    cos_theta = cos(theta);
    g_lut = zeros(size(gamma_lut));
    inv_sqrt2 = 1 / sqrt(2);
    for k = 1:numel(gamma_lut)
        g_lut(k) = mean(erf((gamma_lut(k) * inv_sqrt2) * cos_theta) .* cos_theta);
    end
end

rho_re = mean(double(re_in(:)) .* u_re(:));
rho_im = mean(double(im_in(:)) .* u_im(:));
rho = 0.5 * (rho_re + rho_im);

if ~isfinite(rho) || rho <= 0
    gamma_hat = 0;
elseif rho >= g_lut(end)
    gamma_hat = gamma_lut(end);
else
    gamma_hat = interp1(g_lut, gamma_lut, rho, "pchip");
end
end
