function sample_id = zc1bit_b0_median_sample_id(detail_table, baseline_key)
%ZC1BIT_B0_MEDIAN_SAMPLE_ID 选择 B0 跨单元平均 SSIM 位于中位数的测试样本。
% detail_table 应已过滤为 Track=="independent_recalibrated"。

baseline_key = string(baseline_key);
sub = detail_table(detail_table.VariantKey == baseline_key, :);
if isempty(sub)
    error("zc1bit:PlotNoBaselineRows", "明细表没有基线 %s 的逐样本行。", baseline_key);
end
ids = unique(double(sub.SampleID));
means = nan(numel(ids), 1);
for k = 1:numel(ids)
    mask = double(sub.SampleID) == ids(k);
    means(k) = mean(double(sub.SSIM(mask)), "omitnan");
end
valid = isfinite(means);
if ~any(valid)
    error("zc1bit:PlotNoFiniteBaseline", "基线逐样本 SSIM 没有有限值。");
end
med = median(means(valid), "omitnan");
[~, local_idx] = min(abs(means - med));
sample_id = ids(local_idx);
end
