function value = compute_enl(img, top, left, window_size)
%COMPUTE_ENL 在固定均匀区 ROI 上计算等效视数 ENL = mu^2 / var（基于强度 patch.^2）。

patch = double(img(top : top + window_size - 1, left : left + window_size - 1));
intensity = patch .^ 2;
mu = mean(intensity(:));
variance = var(intensity(:), 0);
if variance <= eps(max(mu ^ 2, 1))
    value = Inf;
else
    value = mu ^ 2 / variance;
end
end
