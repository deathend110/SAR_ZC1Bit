function value = compute_entropy(img, num_bins)
%COMPUTE_ENTROPY 计算归一化幅度图的 256-bin 信息熵（bit）。

if nargin < 2 || isempty(num_bins)
    num_bins = 256;
end

values = min(max(double(img(:)), 0), 1);
counts = histcounts(values, linspace(0, 1, num_bins + 1));
probabilities = counts / max(sum(counts), 1);
probabilities = probabilities(probabilities > 0);
value = -sum(probabilities .* log2(probabilities));
end
