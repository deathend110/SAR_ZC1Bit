function [channel_out, diag] = apply_zc_pseudo15bit(channel_1bit, rule_key, options)
%APPLY_ZC_PSEUDO15BIT 对标准 1-bit 复回波施加伪 1.5-bit 过零置零重构（{+1, 0, -1}）。
% 仅作用于 1-bit 量化之后、距离压缩之前，实部与虚部各自独立处理，矩阵尺寸保持不变。
%
% 支持规则：B0, R1..R9, L1
% options 字段：
%   Axis         - "range" (dim 1, 默认) | "azimuth" (dim 2) | "both_union" | "both_seq_ra" | "both_seq_ar"
%   DecisionMode - "parallel" (默认) | "serial"
%   RemoveDC     - false (默认) | true

if nargin < 3 || isempty(options)
    options = struct();
end
axis_mode = "range";
if isfield(options, "Axis") && ~isempty(options.Axis)
    axis_mode = string(options.Axis);
end
decision_mode = "parallel";
if isfield(options, "DecisionMode") && ~isempty(options.DecisionMode)
    decision_mode = string(options.DecisionMode);
end
remove_dc = false;
if isfield(options, "RemoveDC") && ~isempty(options.RemoveDC)
    remove_dc = logical(options.RemoveDC);
end

rule_key = string(rule_key);
spec = resolve_rule_spec(rule_key);
if ~ismember(decision_mode, ["parallel", "serial"])
    error("zc1bit:DecisionMode", "未知判定模式：%s（仅支持 parallel / serial）。", decision_mode);
end
if ~ismember(axis_mode, ["range", "azimuth", "both_union", "both_seq_ra", "both_seq_ar"])
    error("zc1bit:AxisMode", "未知置零方向：%s。", axis_mode);
end

re_in = real(channel_1bit);
im_in = imag(channel_1bit);
if ~all(re_in == 1 | re_in == -1, "all") || ~all(im_in == 1 | im_in == -1, "all")
    error("zc1bit:Invalid1BitInput", "伪 1.5-bit 输入实部与虚部必须严格属于 {+1, -1}。");
end

re_out = apply_channel_axis(re_in, spec, axis_mode, decision_mode);
im_out = apply_channel_axis(im_in, spec, axis_mode, decision_mode);

% 置零率与天然直流偏置统一在去 DC 之前统计，反映真实三电平输出特征。
zero_ratio_real = mean(re_out(:) == 0);
zero_ratio_imag = mean(im_out(:) == 0);
mean_bias_real = mean(double(re_out(:)));
mean_bias_imag = mean(double(im_out(:)));

if remove_dc
    re_out = re_out - cast(mean_bias_real, "like", re_out);
    im_out = im_out - cast(mean_bias_imag, "like", im_out);
end

channel_out = complex(re_out, im_out);

diag = struct( ...
    "RuleKey", rule_key, ...
    "Axis", axis_mode, ...
    "DecisionMode", decision_mode, ...
    "RemoveDC", remove_dc, ...
    "ZeroRatioReal", zero_ratio_real, ...
    "ZeroRatioImag", zero_ratio_imag, ...
    "ZeroRatio", 0.5 * (zero_ratio_real + zero_ratio_imag), ...
    "MeanBiasReal", mean_bias_real, ...
    "MeanBiasImag", mean_bias_imag, ...
    "AbsMeanBias", hypot(mean_bias_real, mean_bias_imag), ...
    "PostDCMeanBiasReal", mean(double(re_out(:))), ...
    "PostDCMeanBiasImag", mean(double(im_out(:))));
end

function spec = resolve_rule_spec(rule_key)
% 统一四布尔映射表：[zero_left_p2n, zero_right_p2n, zero_left_n2p, zero_right_n2p]
switch rule_key
    case "B0", flags = [false, false, false, false];
    case "R1", flags = [true,  true,  true,  true];
    case "R2", flags = [true,  false, false, false];
    case "R3", flags = [false, true,  false, false];
    case "R4", flags = [false, false, false, true];
    case "R5", flags = [false, false, true,  false];
    case "R6", flags = [true,  false, false, true];
    case "R7", flags = [false, true,  true,  false];
    case "L1", flags = [true,  false, true,  false];
    case "R8", flags = [true,  true,  false, false];
    case "R9", flags = [false, false, true,  true];
    otherwise
        error("zc1bit:UnknownRule", "未知过零置零规则：%s", rule_key);
end
spec = struct( ...
    "RuleKey", rule_key, ...
    "ZeroLeftP2N", flags(1), ...
    "ZeroRightP2N", flags(2), ...
    "ZeroLeftN2P", flags(3), ...
    "ZeroRightN2P", flags(4), ...
    "IsNoOp", ~any(flags));
end

function y = apply_channel_axis(x, spec, axis_mode, decision_mode)
if spec.IsNoOp
    y = x;
    return;
end
switch axis_mode
    case "range"
        y = apply_1d_dim1(x, spec, decision_mode);
    case "azimuth"
        y = apply_1d_dim1(x.', spec, decision_mode).';
    case "both_union"
        yr = apply_1d_dim1(x, spec, decision_mode);
        ya = apply_1d_dim1(x.', spec, decision_mode).';
        y = x;
        y(yr == 0 | ya == 0) = 0;
    case "both_seq_ra"
        y = apply_1d_dim1(x, spec, decision_mode);
        y = apply_1d_dim1(y.', spec, decision_mode).';
    case "both_seq_ar"
        y = apply_1d_dim1(x.', spec, decision_mode).';
        y = apply_1d_dim1(y, spec, decision_mode);
end
end

function y = apply_1d_dim1(x, spec, decision_mode)
% 沿第 1 维执行过零置零；支持列向量与二维矩阵（向量化并行列处理）。
[nr, nc] = size(x);
y = x;
if nr < 2 || nc == 0
    return;
end

if decision_mode == "parallel"
    p2n = (x(1:end - 1, :) == 1) & (x(2:end, :) == -1);
    n2p = (x(1:end - 1, :) == -1) & (x(2:end, :) == 1);
    mask = false(nr, nc);
    if spec.ZeroLeftP2N,  mask(1:end - 1, :) = mask(1:end - 1, :) | p2n; end
    if spec.ZeroRightP2N, mask(2:end, :)     = mask(2:end, :)     | p2n; end
    if spec.ZeroLeftN2P,  mask(1:end - 1, :) = mask(1:end - 1, :) | n2p; end
    if spec.ZeroRightN2P, mask(2:end, :)     = mask(2:end, :)     | n2p; end
    y(mask) = 0;
else
    for k = 1:nr - 1
        top = y(k, :);
        bot = y(k + 1, :);
        p2n = (top == 1) & (bot == -1);
        n2p = (top == -1) & (bot == 1);
        if spec.ZeroLeftP2N && any(p2n),  y(k, p2n) = 0; end
        if spec.ZeroRightP2N && any(p2n), y(k + 1, p2n) = 0; end
        if spec.ZeroLeftN2P && any(n2p),  y(k, n2p) = 0; end
        if spec.ZeroRightN2P && any(n2p), y(k + 1, n2p) = 0; end
    end
end
end
