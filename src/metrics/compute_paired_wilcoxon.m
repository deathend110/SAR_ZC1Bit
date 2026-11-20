function stats_out = compute_paired_wilcoxon(values_a, values_b)
%COMPUTE_PAIRED_WILCOXON 执行双侧配对 Wilcoxon 符号秩检验（values_a - values_b）。
% 输入先统一转为双精度列向量，避免行列混用导致索引越界。

a = double(values_a(:));
b = double(values_b(:));
if numel(a) ~= numel(b) || isempty(a)
    error("zc1bit:PairedLength", "配对样本向量长度必须非空且相等。");
end
if any(~isfinite(a)) || any(~isfinite(b))
    error("zc1bit:NonFinitePairedValues", "Wilcoxon 配对检验输入存在非有限值。");
end

diff_vals = a - b;
mean_diff = mean(diff_vals);
std_diff = std(diff_vals);
median_diff = median(diff_vals);

has_signrank = exist("signrank", "file") == 2 || exist("signrank", "file") == 5;
if ~has_signrank
    p_value = NaN;
    signed_rank = NaN;
elseif all(diff_vals == 0)
    p_value = 1;
    signed_rank = 0;
else
    [p_value, ~, sr_stats] = signrank(a, b);
    if isfield(sr_stats, "signedrank")
        signed_rank = double(sr_stats.signedrank);
    else
        signed_rank = NaN;
    end
end

stats_out = struct( ...
    "N", numel(diff_vals), ...
    "MeanDiff", mean_diff, ...
    "StdDiff", std_diff, ...
    "MedianDiff", median_diff, ...
    "PValue", double(p_value), ...
    "SignedRank", signed_rank);
end
