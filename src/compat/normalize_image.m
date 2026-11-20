function y = normalize_image(x)
%NORMALIZE_IMAGE 将实/复图像幅值统一归一化到 [0, 1]。

mag = abs(x);
peak = max(mag(:));
if peak == 0
    y = mag;
else
    y = mag / peak;
end
end
