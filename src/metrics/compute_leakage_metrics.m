function [off_ratio, range_ratio, azimuth_ratio] = compute_leakage_metrics(X, reference_matrix, threshold_ratio)
%COMPUTE_LEAKAGE_METRICS 计算距离压缩频谱带外能量泄漏比。

if nargin < 3 || isempty(threshold_ratio)
    threshold_ratio = 0.35;
end

ref_spec = abs(fftshift(fft2(reference_matrix)));
support_mask = ref_spec >= threshold_ratio * max(ref_spec(:));
spec = abs(fftshift(fft2(X))) .^ 2;
total_energy = sum(spec(:)) + eps;
off_ratio = sum(spec(~support_mask), "all") / total_energy;

range_profile = sum(spec, 2);
azimuth_profile = sum(spec, 1).';
range_mask = any(support_mask, 2);
azimuth_mask = any(support_mask, 1).';
range_ratio = sum(range_profile(~range_mask)) / (sum(range_profile) + eps);
azimuth_ratio = sum(azimuth_profile(~azimuth_mask)) / (sum(azimuth_profile) + eps);
end
