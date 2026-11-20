function [img_out, complex_roi] = fast_focus_complex_roi(channel, context, range_q, azimuth_q)
%FAST_FOCUS_COMPLEX_ROI 使用预计算频域核直接完成裁剪、RCMC 与方位聚焦，输出归一化 ROI。

r_idx = find(abs(context.range_factors - range_q) < 1e-12, 1);
a_idx = find(abs(context.azimuth_factors - azimuth_q) < 1e-12, 1);
if isempty(r_idx) || isempty(a_idx)
    error("zc1bit:FastContextFactor", "快路径上下文未预计算倍率 R=%g, A=%g。", range_q, azimuth_q);
end

nr_up = context.range_nr_ups(r_idx);
na_up = context.azimuth_na_ups(a_idx);
if ~isequal(size(channel), [nr_up, na_up])
    error("zc1bit:FastInputSize", "输入尺寸 [%d %d] 与快路径网格 [%d %d] 不一致。", ...
        size(channel, 1), size(channel, 2), nr_up, na_up);
end

range_spectrum = fftshift(fft(fftshift(channel, 1), [], 1), 1) .* context.range_kernels{r_idx};
if nr_up > context.nrn
    range_spectrum = range_spectrum(context.range_crop_indices{r_idx}, :);
end
range_doppler = fftshift(fft(fftshift(range_spectrum, 2), [], 2), 2);
if na_up > context.nan
    range_doppler = range_doppler(:, context.azimuth_crop_indices{a_idx});
end

range_azimuth = fftshift(ifft(fftshift(range_doppler .* context.rcmc_kernel, 1), [], 1), 1);
image = fftshift(ifft(fftshift(range_azimuth .* context.azimuth_kernel, 2), [], 2), 2);
complex_roi = image(context.roi_rows, context.roi_cols);
img_out = normalize_image(abs(complex_roi));
end
