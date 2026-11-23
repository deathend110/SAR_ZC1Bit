function [channel_out, diag] = apply_sft_p15_reconstruct(channel_1bit, carrier, method_key, options)
%APPLY_SFT_P15_RECONSTRUCT 面向 SFT 阈值的严格伪 1.5-bit（{+1, 0, -1}）接收端重构。
% 硬性契约：
%   1. 输入 channel_1bit 实部与虚部必须严格属于 {+1, -1}；
%   2. 输出 channel_out 实部与虚部必须严格属于 {+1, 0, -1}，尺寸与精度保持不变；
%   3. 绝不接收未量化回波或能量尺度 sigma_hat。
%
% 支持方法（method_key）：
%   "B0_par_raw"          - 原始 1-bit 基线（全盲）
%   "R1_par_raw"          - 并行双沿双边过零置零（全盲，原版对照）
%   "P15_1_carrier_agree" - 方案 1：载波置信度门控置零 (b .* u > Tau -> 0)
%   "P15_2_rlp_she"       - 方案 2：游程保护 + 可调半宽对称置零（全盲 SHE-ZC）
%   "P15_3_arcsin_lin"    - 方案 3：arcsin 反翘曲单侧置零 (|u_mid| > TauLin 削多数侧)
%   "P15_4_dual_domain"   - 方案 4：波峰置信度 + 游程保护边沿双域联合置零

if nargin < 4 || isempty(options)
    options = struct();
end

method_key = string(method_key);
allowed_methods = [ ...
    "B0_par_raw", "R1_par_raw", ...
    "P15_1_carrier_agree", "P15_2_rlp_she", ...
    "P15_3_arcsin_lin", "P15_4_dual_domain"];
if ~isscalar(method_key) || ~ismember(method_key, allowed_methods)
    error("zc1bit:UnknownMethod", "未知伪 1.5-bit SFT 重构方法：%s", method_key);
end

axis_mode = "range";
if isfield(options, "Axis") && ~isempty(options.Axis)
    axis_mode = string(options.Axis);
end
if ~isscalar(axis_mode) || ~ismember(axis_mode, ["range", "azimuth"])
    error("zc1bit:AxisMode", "未知工作主轴：%s（仅支持 range / azimuth）。", axis_mode);
end

re_in = real(channel_1bit);
im_in = imag(channel_1bit);
if isempty(channel_1bit) || ~isnumeric(channel_1bit) || ...
        ~all(re_in == 1 | re_in == -1, "all") || ~all(im_in == 1 | im_in == -1, "all")
    error("zc1bit:Invalid1BitInput", "输入 channel_1bit 实部与虚部必须严格属于 {+1, -1}。");
end

needs_carrier = ismember(method_key, ["P15_1_carrier_agree", "P15_3_arcsin_lin", "P15_4_dual_domain"]);
if needs_carrier
    if nargin < 2 || isempty(carrier) || ~isnumeric(carrier) || ...
            ~isequal(size(carrier), size(channel_1bit)) || ~all(isfinite(carrier), "all")
        error("zc1bit:InvalidCarrier", "方法 %s 需要与输入同尺寸且有限的单位模 SFT 载波。", method_key);
    end
    u_re = real(carrier);
    u_im = imag(carrier);
else
    u_re = [];
    u_im = [];
end

tau_val = NaN;
hw_val = NaN;
lmin_val = NaN;
taulin_val = NaN;
prot_re = 0;
prot_im = 0;

switch method_key
    case "B0_par_raw"
        re_out = re_in;
        im_out = im_in;

    case "R1_par_raw"
        qopts = struct("Axis", axis_mode, "DecisionMode", "parallel", "RemoveDC", false);
        [channel_out, qdiag] = apply_zc_pseudo15bit(channel_1bit, "R1", qopts);
        re_out = real(channel_out);
        im_out = imag(channel_out);
        assert_strict_three_level(re_out, im_out);
        diag = build_diag_struct(method_key, axis_mode, tau_val, hw_val, lmin_val, taulin_val, ...
            re_out, im_out, 0, qdiag);
        return;

    case "P15_1_carrier_agree"
        tau_val = parse_nonneg_finite_scalar(options, "Tau", "P15_1_carrier_agree");
        re_out = apply_p1_channel(re_in, u_re, tau_val);
        im_out = apply_p1_channel(im_in, u_im, tau_val);

    case "P15_2_rlp_she"
        if isfield(options, "HalfWidth") && ~isempty(options.HalfWidth)
            hw_val = double(options.HalfWidth);
        else
            hw_val = 1;
        end
        if ~isscalar(hw_val) || ~isfinite(hw_val) || hw_val < 1 || mod(hw_val, 1) ~= 0
            error("zc1bit:InvalidHyperparam", "HalfWidth 必须为正整数标量。");
        end
        lmin_val = parse_min_run_length(options, hw_val);
        [re_out, prot_re] = apply_p2_channel(re_in, axis_mode, hw_val, lmin_val);
        [im_out, prot_im] = apply_p2_channel(im_in, axis_mode, hw_val, lmin_val);

    case "P15_3_arcsin_lin"
        taulin_val = parse_nonneg_finite_scalar(options, "TauLin", "P15_3_arcsin_lin");
        re_out = apply_p3_channel(re_in, u_re, axis_mode, taulin_val);
        im_out = apply_p3_channel(im_in, u_im, axis_mode, taulin_val);

    case "P15_4_dual_domain"
        tau_val = parse_nonneg_finite_scalar(options, "Tau", "P15_4_dual_domain");
        hw_val = 1;
        lmin_val = parse_min_run_length(options, hw_val);
        [re_out, prot_re] = apply_p4_channel(re_in, u_re, axis_mode, tau_val, lmin_val);
        [im_out, prot_im] = apply_p4_channel(im_in, u_im, axis_mode, tau_val, lmin_val);
end

re_out = cast(re_out, "like", re_in);
im_out = cast(im_out, "like", im_in);
assert_strict_three_level(re_out, im_out);
channel_out = complex(re_out, im_out);
prot_ratio = 0.5 * (prot_re + prot_im);
diag = build_diag_struct(method_key, axis_mode, tau_val, hw_val, lmin_val, taulin_val, ...
    re_out, im_out, prot_ratio, []);
end

function assert_strict_three_level(re_out, im_out)
valid_re = all(re_out == -1 | re_out == 0 | re_out == 1, "all");
valid_im = all(im_out == -1 | im_out == 0 | im_out == 1, "all");
if ~valid_re || ~valid_im
    error("zc1bit:InvalidPseudo15Output", "严格伪 1.5-bit 输出实部与虚部必须全部属于 {-1, 0, +1}。");
end
end

function val = parse_nonneg_finite_scalar(options, field_name, method_key)
if ~isfield(options, field_name) || isempty(options.(field_name))
    error("zc1bit:MissingHyperparam", "方法 %s 必须提供参数 options.%s。", method_key, field_name);
end
val = double(options.(field_name));
if ~isscalar(val) || ~isfinite(val) || val < 0
    error("zc1bit:InvalidHyperparam", "方法 %s 的参数 %s 必须为非负有限标量。", method_key, field_name);
end
end

function l_min = parse_min_run_length(options, half_width)
if ~isfield(options, "MinRunLength") || isempty(options.MinRunLength)
    error("zc1bit:MissingHyperparam", "游程保护方案必须提供参数 options.MinRunLength。");
end
l_min = double(options.MinRunLength);
if ~isscalar(l_min) || isnan(l_min) || l_min < 1
    error("zc1bit:InvalidHyperparam", "MinRunLength 必须为 >= 1 或 Inf 的标量。");
end
if ~isinf(l_min) && mod(l_min, 1) ~= 0
    error("zc1bit:InvalidHyperparam", "有限 MinRunLength 必须为整数。");
end
if half_width > 1 && l_min < half_width
    error("zc1bit:InvalidHyperparam", "当 HalfWidth > 1 时要求 MinRunLength >= HalfWidth。");
end
end

function y = apply_p1_channel(x, u, tau)
y = x;
if tau >= 1
    return;
end
mask = (double(x) .* double(u)) > tau;
y(mask) = 0;
end

function [y, prot_ratio] = apply_p2_channel(x, axis_mode, half_width, l_min)
if axis_mode == "azimuth"
    [mask, prot_ratio] = compute_rlp_she_mask_dim1(x.', half_width, l_min);
    mask = mask.';
else
    [mask, prot_ratio] = compute_rlp_she_mask_dim1(x, half_width, l_min);
end
y = x;
y(mask) = 0;
end

function [mask, prot_ratio] = compute_rlp_she_mask_dim1(x, half_width, l_min)
[nr, nc] = size(x);
mask = false(nr, nc);
if nr < 2 || nc == 0
    prot_ratio = 0;
    return;
end

zc_edge = (x(1:end - 1, :) ~= x(2:end, :));
n_zc = sum(zc_edge(:));
if n_zc == 0
    prot_ratio = 0;
    return;
end
if isinf(l_min)
    prot_ratio = 1;
    return;
end

% 按列内切断计算每个位置的同号游程长度，防止跨列误合并
is_start = false(nr, nc);
is_start(1, :) = true;
is_start(2:end, :) = (x(2:end, :) ~= x(1:end - 1, :));
run_id = cumsum(is_start(:));
run_starts = find(is_start(:));
run_ends = [run_starts(2:end) - 1; nr * nc];
run_lens = run_ends - run_starts + 1;
elem_run_len = reshape(run_lens(run_id), nr, nc);

valid_edge = zc_edge & (elem_run_len(1:end - 1, :) >= l_min) & (elem_run_len(2:end, :) >= l_min);
prot_ratio = sum(zc_edge(:) & ~valid_edge(:)) / n_zc;

for d = 1:half_width
    if d > nr - 1
        continue;
    end
    e_d = valid_edge & (elem_run_len(1:end - 1, :) >= d) & (elem_run_len(2:end, :) >= d);
    mask(1 : nr - d, :) = mask(1 : nr - d, :) | e_d(d : nr - 1, :);
    mask(1 + d : nr, :) = mask(1 + d : nr, :) | e_d(1 : nr - d, :);
end
end

function y = apply_p3_channel(x, u, axis_mode, tau_lin)
if axis_mode == "azimuth"
    mask = compute_p3_mask_dim1(x.', u.', tau_lin).';
else
    mask = compute_p3_mask_dim1(x, u, tau_lin);
end
y = x;
y(mask) = 0;
end

function mask = compute_p3_mask_dim1(x, u, tau_lin)
[nr, nc] = size(x);
mask = false(nr, nc);
if nr < 2 || nc == 0 || tau_lin >= 1
    return;
end

zc_edge = (x(1:end - 1, :) ~= x(2:end, :));
u_mid = 0.5 * (double(u(1:end - 1, :)) + double(u(2:end, :)));
hi_edge = zc_edge & (abs(u_mid) > tau_lin);
if ~any(hi_edge(:))
    return;
end

% u_mid < 0 表示真实信号 S > 0（削 +1 侧）；u_mid > 0 表示真实信号 S < 0（削 -1 侧）
target_pos = hi_edge & (u_mid < 0);
target_neg = hi_edge & (u_mid > 0);
mask(1:end - 1, :) = mask(1:end - 1, :) | ...
    (target_pos & (x(1:end - 1, :) == 1)) | (target_neg & (x(1:end - 1, :) == -1));
mask(2:end, :) = mask(2:end, :) | ...
    (target_pos & (x(2:end, :) == 1)) | (target_neg & (x(2:end, :) == -1));
end

function [y, prot_ratio] = apply_p4_channel(x, u, axis_mode, tau, l_min)
if tau >= 1
    mask_agree = false(size(x));
else
    mask_agree = (double(x) .* double(u)) > tau;
end
if axis_mode == "azimuth"
    [mask_rlp, prot_ratio] = compute_rlp_she_mask_dim1(x.', 1, l_min);
    mask_rlp = mask_rlp.';
else
    [mask_rlp, prot_ratio] = compute_rlp_she_mask_dim1(x, 1, l_min);
end
y = x;
y(mask_agree | mask_rlp) = 0;
end

function diag = build_diag_struct(method_key, axis_mode, tau_val, hw_val, lmin_val, taulin_val, ...
    re_out, im_out, prot_ratio, r1_diag)

if ~isempty(r1_diag)
    zr_re = r1_diag.ZeroRatioReal;
    zr_im = r1_diag.ZeroRatioImag;
    zr = r1_diag.ZeroRatio;
    mb_re = r1_diag.MeanBiasReal;
    mb_im = r1_diag.MeanBiasImag;
    abm = r1_diag.AbsMeanBias;
else
    zr_re = mean(re_out(:) == 0);
    zr_im = mean(im_out(:) == 0);
    zr = 0.5 * (zr_re + zr_im);
    mb_re = mean(double(re_out(:)));
    mb_im = mean(double(im_out(:)));
    abm = hypot(mb_re, mb_im);
end

diag = struct( ...
    "MethodKey", method_key, ...
    "Axis", axis_mode, ...
    "Tau", tau_val, ...
    "HalfWidth", hw_val, ...
    "MinRunLength", lmin_val, ...
    "TauLin", taulin_val, ...
    "ZeroRatioReal", zr_re, ...
    "ZeroRatioImag", zr_im, ...
    "ZeroRatio", zr, ...
    "ProtectedEdgeRatio", prot_ratio, ...
    "MeanBiasReal", mb_re, ...
    "MeanBiasImag", mb_im, ...
    "AbsMeanBias", abm);
end
