function output = crop_azimuth_spectrum(input, target_cols)
%CROP_AZIMUTH_SPECTRUM 沿第 2 维（方位向）裁剪中心多普勒频带回目标列数。

current_cols = size(input, 2);
if target_cols > current_cols
    error("zc1bit:CropSize", "目标方位点数 (%d) 不能大于当前尺寸 (%d)。", target_cols, current_cols);
elseif target_cols == current_cols
    output = input;
    return;
end

spectrum = fftshift(fft(input, [], 2), 2);
center = floor(current_cols / 2) + 1;
half_width = floor(target_cols / 2);
if mod(target_cols, 2) == 0
    indices = center - half_width : center + half_width - 1;
else
    indices = center - half_width : center + half_width;
end
output = ifft(ifftshift(spectrum(:, indices), 2), [], 2);
end
