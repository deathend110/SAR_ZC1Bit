function signal_up = two_dim_upsample_fft(signal, range_q, azimuth_q)
%TWO_DIM_UPSAMPLE_FFT 先距离向（dim 1）后方位向（dim 2）执行二维 FFT 补零插值。
% 参数顺序严格与 R{range_q}A{azimuth_q} 一致：先 range_q，后 azimuth_q。

signal_up = signal;
if range_q > 1
    signal_up = range_upsample_fft(signal_up, range_q);
end
if azimuth_q > 1
    signal_up = azimuth_upsample_fft(signal_up, azimuth_q);
end
end
