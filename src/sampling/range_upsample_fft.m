function signal_up = range_upsample_fft(signal, range_q)
%RANGE_UPSAMPLE_FFT 沿第 1 维（距离向）执行频域中心补零带限插值。

if ~isscalar(range_q) || ~isfinite(range_q) || range_q < 1
    error("zc1bit:SamplingFactor", "距离上采样倍率必须是不小于 1 的有限数。");
end

[nr, na] = size(signal);
nr_up = round(range_q * nr);
if nr_up == nr
    signal_up = signal;
    return;
end

spectrum = fftshift(fft(signal, [], 1), 1);
pad_total = nr_up - nr;
pad_top = floor(pad_total / 2);
pad_bottom = pad_total - pad_top;
spectrum_up = [ ...
    zeros(pad_top, na, "like", spectrum); ...
    spectrum; ...
    zeros(pad_bottom, na, "like", spectrum)];
signal_up = ifft(ifftshift(spectrum_up, 1), [], 1) * (nr_up / nr);
end
