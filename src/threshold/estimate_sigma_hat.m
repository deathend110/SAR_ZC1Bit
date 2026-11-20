function sigma_hat = estimate_sigma_hat(signal_up)
%ESTIMATE_SIGMA_HAT 基于上采样后复回波估计瑞利幅度尺度 sigma_hat。
% 锁定公式：sigma_hat = sqrt(2/pi) * mean(|S^up|)。

if isempty(signal_up) || ~isnumeric(signal_up)
    error("zc1bit:SignalInput", "上采样回波必须为非空数值矩阵。");
end
sigma_hat = sqrt(2 / pi) * mean(abs(signal_up(:)));
end
