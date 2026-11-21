function T = zc1bit_build_calibration_delta(S, baseline_key)
%ZC1BIT_BUILD_CALIBRATION_DELTA 由 search_best_params 构造校准集 Δ 指标表。
% 输出列：Branch, ThresholdType, Q, VariantKey, DeltaSSIM_Mean, DeltaPSNR_Mean。

baseline_key = string(baseline_key);
rows = height(S);
mask = S.VariantKey ~= baseline_key;
idx = find(mask);
Branch = strings(numel(idx), 1);
ThresholdType = strings(numel(idx), 1);
Q = zeros(numel(idx), 1);
VariantKey = strings(numel(idx), 1);
DeltaSSIM_Mean = zeros(numel(idx), 1);
DeltaPSNR_Mean = zeros(numel(idx), 1);

for k = 1:numel(idx)
    i = idx(k);
    b0_mask = S.UnitKey == S.UnitKey(i) & S.VariantKey == baseline_key;
    if ~any(b0_mask)
        error("zc1bit:PlotMissingCalibrationB0", "校准搜索缺少 B0 行：%s", S.UnitKey(i));
    end
    j = find(b0_mask, 1);
    Branch(k) = S.Branch(i);
    ThresholdType(k) = S.ThresholdType(i);
    Q(k) = S.Q(i);
    VariantKey(k) = S.VariantKey(i);
    DeltaSSIM_Mean(k) = S.Calibration_SSIM_Mean(i) - S.Calibration_SSIM_Mean(j);
    DeltaPSNR_Mean(k) = S.Calibration_PSNR_Mean(i) - S.Calibration_PSNR_Mean(j);
end
T = table(Branch, ThresholdType, Q, VariantKey, DeltaSSIM_Mean, DeltaPSNR_Mean);
end
