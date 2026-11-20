function [threshold, sigma_hat] = build_splitrt_threshold(signal_up, As, rand_stream)
%BUILD_SPLITRT_THRESHOLD 构造可分离二维随机阈值（SplitRT）。
% 先分别生成距离列相位 phi_r 与方位行相位 phi_a，再相位相加取复指数。

if ~isscalar(As) || ~isfinite(As) || As < 0
    error("zc1bit:InvalidAs", "SplitRT 幅度系数 As 必须为非负有限标量。");
end

[nr_up, na_up] = size(signal_up);
if nargin < 3 || isempty(rand_stream)
    phi_r = 2 * pi * rand(nr_up, 1);
    phi_a = 2 * pi * rand(1, na_up);
else
    phi_r = 2 * pi * rand(rand_stream, nr_up, 1);
    phi_a = 2 * pi * rand(rand_stream, 1, na_up);
end

sigma_hat = estimate_sigma_hat(signal_up);
threshold = (As * sigma_hat) * exp(1i * (phi_r + phi_a));
end
