function record = zc1bit_prepare_run(cfg)
%ZC1BIT_PREPARE_RUN 校验配置、成像参数与回波文件头，构建或复用 70 样本 35/35 分层清单。
% 仅读取 rstart*.mat 的变量元数据（whos），不在 preflight 阶段加载大回波矩阵。

validate_config_schema(cfg);
zc1bit_validate_roots(cfg);

S60 = load(cfg.parameter_file_60);
required_s60 = ["nrn", "nan", "fc", "B", "Fs", "prf", "R0", "C", "v", ...
    "Tp", "Ta", "gama", "lambda", "fnrn", "fnan", "tnrn", "tnan", "R_total", "A_num"];
if ~all(isfield(S60, required_s60))
    error("zc1bit:ImagingFields", "FS60 参数文件缺少必要成像字段。");
end
resolve_azimuth_bandwidth(S60);

p = zc1bit_paths(cfg);
record_path = fullfile(p.protocol, "run_record.mat");
resuming = isfile(record_path);

if resuming
    loaded = load(record_path, "record");
    if ~isfield(loaded, "record")
        error("zc1bit:ProtocolRecord", "现有 run_record.mat 缺少 record 字段。");
    end
    record = loaded.record;
    zc1bit_check_parameters(record.signature.config, cfg, "工作流配置");
    zc1bit_check_parameters(record.signature.S60, S60, "成像参数 S60");
    manifest = record.manifest;
    manifest.FilePath = fullfile(cfg.data_root, manifest.Dataset, manifest.File);
    verify_echo_files_exist(manifest);
    record.manifest = manifest;
else
    manifest = build_sample_manifest(cfg, S60);
    record = struct();
    record.signature = struct( ...
        "config", zc1bit_parameters(cfg), ...
        "S60", zc1bit_parameters(S60));
    record.manifest = manifest;
    zc1bit_atomic_save(record_path, struct("record", record));
    writetable(manifest(:, ["SampleID", "DatasetIdx", "LocalSampleIdx", ...
        "Dataset", "File", "CStart", "Split"]), ...
        fullfile(p.protocol, "manifest_70.csv"));
end

[threshold_groups, variant_specs] = zc1bit_build_screen_tasks(cfg, S60);
zc1bit_atomic_save(fullfile(p.protocol, "preflight.mat"), struct( ...
    "record", record, ...
    "threshold_groups", threshold_groups, ...
    "variant_specs", variant_specs));

n_cal = sum(record.manifest.Split == "calibration");
n_test = sum(record.manifest.Split == "test");
fprintf("预检完成：%d 校准样本、%d 测试样本（测试集隔离），阶段一 %d 个阈值配置 × %d 个规则变体。\n", ...
    n_cal, n_test, height(threshold_groups), height(variant_specs));
end

function validate_config_schema(cfg)
if exist("psnr", "file") ~= 2 || exist("ssim", "file") ~= 2
    error("zc1bit:Toolbox", "需要 Image Processing Toolbox 的 psnr 与 ssim 函数。");
end
if ~isstring(cfg.dataset_names) || ~isrow(cfg.dataset_names) || isempty(cfg.dataset_names)
    error("zc1bit:DatasetNames", "dataset_names 必须为非空 string 行向量。");
end
if numel(unique(cfg.dataset_names)) ~= numel(cfg.dataset_names)
    error("zc1bit:DatasetNames", "dataset_names 不能包含重复场景名。");
end
if ~isscalar(cfg.num_samples_per_dataset) || cfg.num_samples_per_dataset < 2 || ...
        mod(cfg.num_samples_per_dataset, 1) ~= 0
    error("zc1bit:SamplesPerDataset", "num_samples_per_dataset 必须为不小于 2 的整数。");
end
if ~isscalar(cfg.nyquist_margin) || cfg.nyquist_margin <= 0 || cfg.nyquist_margin >= 1
    error("zc1bit:NyquistMargin", "nyquist_margin 必须位于 (0, 1) 之间。");
end
if ~isscalar(cfg.raw_decimation_step) || cfg.raw_decimation_step < 1 || ...
        mod(cfg.raw_decimation_step, 1) ~= 0
    error("zc1bit:DecimationStep", "raw_decimation_step 必须为正整数。");
end
end

function manifest = build_sample_manifest(cfg, S60)
num_ds = numel(cfg.dataset_names);
n_per_ds = cfg.num_samples_per_dataset;
total = num_ds * n_per_ds;

SampleID = (1:total).';
DatasetIdx = zeros(total, 1);
LocalSampleIdx = zeros(total, 1);
Dataset = strings(total, 1);
File = strings(total, 1);
FilePath = strings(total, 1);
EchoVariable = strings(total, 1);
CStart = zeros(total, 1);
RawHeight = zeros(total, 1);
RawWidth = zeros(total, 1);
Split = strings(total, 1);

min_rows = 1 + cfg.raw_decimation_step * (S60.nrn - 1);
ptr = 0;
for ds_idx = 1:num_ds
    ds_name = cfg.dataset_names(ds_idx);
    ds_dir = fullfile(cfg.data_root, ds_name);
    if ~isfolder(ds_dir)
        error("zc1bit:MissingDatasetDir", "场景目录不存在：%s", ds_dir);
    end
    mat_files = dir(fullfile(ds_dir, "rstart*.mat"));
    mat_names = sort(string({mat_files.name}));
    if isempty(mat_names)
        error("zc1bit:NoEchoMat", "数据集 %s 中未找到 rstart*.mat。", ds_name);
    end
    pick_idx = mod(cfg.seed, numel(mat_names)) + 1;
    picked_name = mat_names(pick_idx);
    mat_path = fullfile(ds_dir, picked_name);

    info = whos("-file", mat_path);
    if isempty(info)
        error("zc1bit:EmptyMat", "回波 MAT 文件为空：%s", mat_path);
    end
    % 仅选取第一个二维复数 single/double 矩阵作为回波变量。
    valid_var = false(numel(info), 1);
    for v_idx = 1:numel(info)
        valid_var(v_idx) = ismember(string(info(v_idx).class), ["double", "single"]) && ...
            info(v_idx).complex && numel(info(v_idx).size) == 2 && ...
            info(v_idx).size(1) >= min_rows && info(v_idx).size(2) >= S60.nan;
    end
    if ~any(valid_var)
        error("zc1bit:EchoSchema", "回波文件缺少尺寸与复数类型匹配的矩阵：%s", mat_path);
    end
    first_idx = find(valid_var, 1);
    var_info = info(first_idx);

    starts = stratified_starts(var_info.size(2), S60.nan, n_per_ds);
    for local_idx = 1:n_per_ds
        ptr = ptr + 1;
        DatasetIdx(ptr) = ds_idx;
        LocalSampleIdx(ptr) = local_idx;
        Dataset(ptr) = ds_name;
        File(ptr) = picked_name;
        FilePath(ptr) = mat_path;
        EchoVariable(ptr) = string(var_info.name);
        CStart(ptr) = starts(local_idx);
        RawHeight(ptr) = var_info.size(1);
        RawWidth(ptr) = var_info.size(2);
        if mod(local_idx, 2) == 1
            Split(ptr) = "calibration";
        else
            Split(ptr) = "test";
        end
    end
end

manifest = table(SampleID, DatasetIdx, LocalSampleIdx, Dataset, File, ...
    FilePath, EchoVariable, CStart, RawHeight, RawWidth, Split);
end

function starts = stratified_starts(raw_width, window_width, num_samples)
max_start = raw_width - window_width + 1;
if max_start < 1
    error("zc1bit:RawWidthTooShort", "原始回波列数 %d 不足以裁出宽度 %d。", raw_width, window_width);
end
starts = zeros(num_samples, 1);
for idx = 1:num_samples
    starts(idx) = min(max(round((idx - 0.5) / num_samples * max_start), 1), max_start);
end
end

function verify_echo_files_exist(manifest)
unique_paths = unique(manifest.FilePath);
for k = 1:numel(unique_paths)
    if ~isfile(unique_paths(k))
        error("zc1bit:MissingEchoFile", "回波文件不存在：%s", unique_paths(k));
    end
end
end
