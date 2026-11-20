function signal_down = two_dim_downsample_fft(signal, range_q, azimuth_q, S60)
%TWO_DIM_DOWNSAMPLE_FFT 距离压缩后先裁方位向（dim 2），再裁距离向（dim 1）回到基网格。
% 参数顺序严格与 R{range_q}A{azimuth_q} 一致：先 range_q，后 azimuth_q。

signal_down = signal;
if azimuth_q > 1
    signal_down = crop_azimuth_spectrum(signal_down, S60.nan);
end
if range_q > 1
    signal_down = crop_range_spectrum(signal_down, S60.nrn);
end
end
