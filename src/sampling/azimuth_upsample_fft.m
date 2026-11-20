function signal_up = azimuth_upsample_fft(signal, azimuth_q)
%AZIMUTH_UPSAMPLE_FFT 沿第 2 维（方位向）执行频域中心补零带限插值。

if ~isscalar(azimuth_q) || ~isfinite(azimuth_q) || azimuth_q < 1
    error("zc1bit:SamplingFactor", "方位上采样倍率必须是不小于 1 的有限数。");
end

[nr, na] = size(signal);
na_up = round(azimuth_q * na);
if na_up == na
    signal_up = signal;
    return;
end

spectrum = fftshift(fft(signal, [], 2), 2);
pad_total = na_up - na;
pad_left = floor(pad_total / 2);
pad_right = pad_total - pad_left;
spectrum_up = [ ...
    zeros(nr, pad_left, "like", spectrum), ...
    spectrum, ...
    zeros(nr, pad_right, "like", spectrum)];
signal_up = ifft(ifftshift(spectrum_up, 2), [], 2) * (na_up / na);
end
