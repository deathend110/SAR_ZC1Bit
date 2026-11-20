function metrics = compute_image_metrics(img_out, img_gt)
%COMPUTE_IMAGE_METRICS 计算单幅归一化幅度图的 PSNR 与 SSIM，并拦截非法非有限值。

if ~isequal(size(img_out), size(img_gt))
    error("zc1bit:MetricSize", "待评价图像与 GT 尺寸不一致。");
end

psnr_val = double(psnr(img_out, img_gt));
ssim_val = double(ssim(img_out, img_gt));

if isnan(psnr_val) || (isinf(psnr_val) && psnr_val < 0)
    error("zc1bit:InvalidMetric", "PSNR 出现非法值：%g。", psnr_val);
end
if ~isfinite(ssim_val)
    error("zc1bit:InvalidMetric", "SSIM 出现非有限值：%g。", ssim_val);
end

metrics = struct("PSNR", psnr_val, "SSIM", ssim_val);
end
