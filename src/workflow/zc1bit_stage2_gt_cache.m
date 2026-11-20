function gt_cache = zc1bit_stage2_gt_cache(cfg, S60, split_manifest, signals, cache_path, reuse_path)
%ZC1BIT_STAGE2_GT_CACHE 加载、复用或构建指定样本划分的 GT 与固定 64x64 ENL ROI 缓存。

if nargin < 6
    reuse_path = "";
end

expected_params = build_gt_expected_params(cfg, S60, split_manifest);
gt_cache = try_load_one(cache_path, expected_params);
if ~isempty(gt_cache)
    return;
end

if strlength(string(reuse_path)) > 0
    gt_cache = try_load_one(reuse_path, expected_params);
    if ~isempty(gt_cache)
        zc1bit_atomic_save(cache_path, struct("state", struct( ...
            "parameters", expected_params, ...
            "gt_cache", gt_cache, ...
            "completed", true)));
        return;
    end
end

if isempty(signals)
    gt_cache = [];
    return;
end

n_samples = height(split_manifest);
img_gts = cell(n_samples, 1);
gt_entropy = zeros(n_samples, 1);
for idx = 1:n_samples
    img_gt = generate_gt_image(signals{idx}, S60);
    img_gts{idx} = img_gt;
    gt_entropy(idx) = compute_entropy(img_gt, cfg.entropy_num_bins);
end

[is_enl, enl_top, enl_left, ~, roi_manifest] = build_enl_regions(img_gts, split_manifest, cfg);
gt_enl = nan(n_samples, 1);
for idx = find(is_enl).'
    gt_enl(idx) = compute_enl(img_gts{idx}, enl_top(idx), enl_left(idx), cfg.enl_window_size);
end

gt_cache = struct( ...
    "SampleID", split_manifest.SampleID, ...
    "img_gts", {img_gts}, ...
    "gt_entropy", gt_entropy, ...
    "is_enl", is_enl, ...
    "enl_top", enl_top, ...
    "enl_left", enl_left, ...
    "gt_enl", gt_enl, ...
    "roi_manifest", roi_manifest);

zc1bit_atomic_save(cache_path, struct("state", struct( ...
    "parameters", expected_params, ...
    "gt_cache", gt_cache, ...
    "completed", true)));
end

function gt_cache = try_load_one(path_val, expected_params)
gt_cache = [];
if strlength(string(path_val)) == 0 || ~isfile(path_val)
    return;
end
state = zc1bit_load_state(path_val, struct());
if ~isfield(state, "completed") || ~isequal(state.completed, true) || ~isfield(state, "gt_cache")
    return;
end
zc1bit_check_parameters(state.parameters, expected_params, "GT 缓存");
gt_cache = state.gt_cache;
end

function expected_params = build_gt_expected_params(cfg, S60, split_manifest)
expected_params = struct( ...
    "manifest", zc1bit_parameters(split_manifest), ...
    "S60", zc1bit_parameters(S60), ...
    "enl_dataset_names", cfg.enl_dataset_names, ...
    "enl_window_size", cfg.enl_window_size, ...
    "enl_stride", cfg.enl_stride, ...
    "entropy_num_bins", cfg.entropy_num_bins, ...
    "raw_decimation_step", cfg.raw_decimation_step);
end
