function [sample_row, signal60, img_gt, manifest] = zc1bit_select_render_sample( ...
    cfg, D, record, S60, sample_id)
%ZC1BIT_SELECT_RENDER_SAMPLE 选择成像脚本使用的统一测试样本并加载回波。
% sample_id 为空时，自动取 B0 跨单元平均 SSIM 的中位数样本，避免挑选最好样本。

if nargin < 5
    sample_id = [];
end
if ~isfield(cfg, "data_root") || strlength(string(cfg.data_root)) == 0 || ...
        ~isfolder(cfg.data_root)
    error("zc1bit:PlotMissingDataRoot", ...
        "原始回波目录不存在：%s。请检查 inputs/local_paths.m。", ...
        string(cfg.data_root));
end

manifest = record.manifest;
file_paths = strings(height(manifest), 1);
for i = 1:height(manifest)
    file_paths(i) = fullfile(cfg.data_root, manifest.Dataset(i), manifest.File(i));
end
manifest.FilePath = file_paths;

test_manifest = manifest(manifest.Split == "test", :);
if isempty(test_manifest)
    error("zc1bit:PlotEmptyTestManifest", "运行记录中没有测试样本。");
end
if isempty(sample_id)
    sample_id = zc1bit_b0_median_sample_id(D.eval_detail_ind, D.baseline_key);
end
sample_row = test_manifest(test_manifest.SampleID == double(sample_id), :);
if height(sample_row) ~= 1
    error("zc1bit:PlotSampleID", ...
        "测试样本 SampleID=%s 不存在或重复（找到 %d 行）。", ...
        string(sample_id), height(sample_row));
end
if ~isfile(sample_row.FilePath)
    error("zc1bit:PlotMissingEcho", "回波文件不存在：%s", sample_row.FilePath);
end

signals = zc1bit_stage2_load_signals(sample_row, S60, cfg.raw_decimation_step, "绘图样本");
signal60 = signals{1};
img_gt = generate_gt_image(signal60, S60);
end
