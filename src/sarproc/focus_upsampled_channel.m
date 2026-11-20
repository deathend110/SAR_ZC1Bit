function [img_out, RC_crop, RC_raw] = focus_upsampled_channel(channel_in, S60, range_q, azimuth_q)
%FOCUS_UPSAMPLED_CHANNEL 标准慢路径：上采样网格距离压缩 → 频谱裁剪 → RCMC → 方位成像 → ROI 归一化。
% 参数顺序严格为 (channel_in, S60, range_q, azimuth_q)。

grid = resolve_sampling_grid(range_q, azimuth_q, [S60.nrn, S60.nan], S60);
if ~isequal(size(channel_in), grid.upsampled_size)
    error("zc1bit:ChannelSize", "输入回波尺寸 [%d %d] 与网格 [%d %d] 不一致。", ...
        size(channel_in, 1), size(channel_in, 2), grid.upsampled_size(1), grid.upsampled_size(2));
end

[tnrn_up, Fs_up] = range_axis(size(channel_in, 1), grid.q_range_eff, S60);
RC_raw = Range_Compress(channel_in, S60.fc, tnrn_up, S60.gama, ...
    S60.R0, S60.C, Fs_up, S60.Tp);
RC_crop = two_dim_downsample_fft(RC_raw, range_q, azimuth_q, S60);
RCMC_out = RCMC(RC_crop, S60.lambda, S60.fnrn, S60.fnan, ...
    S60.R0, S60.C, S60.v);
IMG = SAR_Imaging(RCMC_out, S60.lambda, S60.Fs, ...
    S60.R0, S60.C, S60.v, S60.tnan, S60.Ta, S60.prf);
img_out = normalize_image(extract_roi(IMG, S60));
end
