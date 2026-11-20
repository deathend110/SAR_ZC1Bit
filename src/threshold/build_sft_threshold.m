function [threshold, sigma_hat, amplitude] = build_sft_threshold(signal_up, grid, params, S60, carrier)
%BUILD_SFT_THRESHOLD 构造单频阈值矩阵 U（覆盖 RSFT / ASFT / 2D-SFT）。

sigma_hat = estimate_sigma_hat(signal_up);
amplitude = sigma_hat / (10 ^ (params.STR_dB / 20));

if nargin < 5 || isempty(carrier)
    carrier = build_sft_carrier(grid, params, S60);
elseif ~isequal(size(carrier), size(signal_up))
    error("zc1bit:CarrierSize", "预计算 SFT 载波尺寸必须与上采样回波一致。");
end

threshold = amplitude * carrier;
end
