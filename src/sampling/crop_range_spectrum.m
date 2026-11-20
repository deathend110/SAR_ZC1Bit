function output = crop_range_spectrum(input, target_rows)
%CROP_RANGE_SPECTRUM 沿第 1 维（距离向）裁剪中心频带回目标行数。

current_rows = size(input, 1);
if target_rows > current_rows
    error("zc1bit:CropSize", "目标距离点数 (%d) 不能大于当前尺寸 (%d)。", target_rows, current_rows);
elseif target_rows == current_rows
    output = input;
    return;
end

spectrum = fftshift(fft(input, [], 1), 1);
center = floor(current_rows / 2) + 1;
half_width = floor(target_rows / 2);
if mod(target_rows, 2) == 0
    indices = center - half_width : center + half_width - 1;
else
    indices = center - half_width : center + half_width;
end
output = ifft(ifftshift(spectrum(indices, :), 1), [], 1);
end
