function [is_enl, tops, lefts, scores, roi_manifest] = build_enl_regions(img_gts, manifest, cfg)
%BUILD_ENL_REGIONS 基于 clean GT 为 ENL 场景样本选定并锁定固定 64x64 均匀区 ROI。

n_samples = height(manifest);
if numel(img_gts) ~= n_samples
    error("zc1bit:ENLSampleCount", "GT 图像元胞长度必须与样本清单行数一致。");
end

is_enl = ismember(string(manifest.Dataset), string(cfg.enl_dataset_names));
tops = zeros(n_samples, 1);
lefts = zeros(n_samples, 1);
scores = nan(n_samples, 1);

for idx = find(is_enl).'
    [tops(idx), lefts(idx), scores(idx)] = select_uniform_roi( ...
        img_gts{idx}, cfg.enl_window_size, cfg.enl_stride);
end

roi_manifest = manifest(is_enl, :);
n_roi = height(roi_manifest);
roi_manifest.ROITop = tops(is_enl);
roi_manifest.ROILeft = lefts(is_enl);
roi_manifest.ROIHeight = repmat(cfg.enl_window_size, n_roi, 1);
roi_manifest.ROIWidth = repmat(cfg.enl_window_size, n_roi, 1);
roi_manifest.UniformityScore = scores(is_enl);
end
