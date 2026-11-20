function [Noise, stats] = gaussian(signal, SNR_dB, verbose)
%GAUSSIAN 生成指定输入 SNR_dB 的零均值复高斯加性白噪声。
% 仅返回噪声矩阵与有限样本统计，不直接返回 signal + Noise。

if nargin < 3
    verbose = true;
end

% 以 1-bit 量化前原始复回波的平均功率定义输入信噪比
signal_power = mean(abs(signal(:)) .^ 2);
target_noise_power = signal_power / (10 ^ (SNR_dB / 10));
target_noise_std = sqrt(target_noise_power);

[rows, cols] = size(signal);
scale_factor = target_noise_std / sqrt(2);

Noise = scale_factor * (randn(rows, cols) + 1j * randn(rows, cols));
Noise = Noise - mean(Noise, "all");

actual_noise_power = mean(abs(Noise(:)) .^ 2);
actual_SNR_dB = 10 * log10(signal_power / actual_noise_power);
stats = struct( ...
    "TargetSNRdB", SNR_dB, ...
    "ActualSNRdB", actual_SNR_dB, ...
    "SignalPower", signal_power, ...
    "TargetNoisePower", target_noise_power, ...
    "ActualNoisePower", actual_noise_power, ...
    "TargetNoiseStd", target_noise_std, ...
    "ActualNoiseStd", std(Noise(:)), ...
    "NoiseMean", mean(Noise(:)));

if verbose
    fprintf("目标 SNR: %.3f dB | 实际 SNR: %.3f dB | 均值模: %.3e\n", ...
        SNR_dB, actual_SNR_dB, abs(stats.NoiseMean));
end
end
