function cfg = zc1bit_config()
%ZC1BIT_CONFIG SAR-ZC1Bit 集中配置；默认仅执行 preflight 预检，不启动耗时实验。

repo_root = string(fileparts(mfilename("fullpath")));
cfg = struct();
cfg.repo_root = repo_root;
cfg.version = "zc1bit_workflow_v1";
cfg.stage = "preflight";

% 本机路径解耦：优先读取未跟踪的 inputs/local_paths.m
local_paths_file = fullfile(repo_root, "inputs", "local_paths.m");
if isfile(local_paths_file)
    lp = read_local_paths(local_paths_file);
    cfg.data_root = string(lp.data_root);
    cfg.parameter_file_60 = string(lp.parameter_file_60);
    if isfield(lp, "results_root") && strlength(string(lp.results_root)) > 0
        cfg.results_root = string(lp.results_root);
    else
        cfg.results_root = fullfile(repo_root, "results_zc1bit");
    end
    if isfield(lp, "tmp_root") && strlength(string(lp.tmp_root)) > 0
        cfg.tmp_root = string(lp.tmp_root);
    else
        cfg.tmp_root = fullfile(repo_root, "tmp");
    end
else
    cfg.data_root = "";
    cfg.parameter_file_60 = fullfile(repo_root, "inputs", "parameters", "FS60_params.mat");
    cfg.results_root = fullfile(repo_root, "results_zc1bit");
    cfg.tmp_root = fullfile(repo_root, "tmp");
end

% 统一样本与指标基准口径（对齐 EXPERIMENT_PROTOCOL.md）
cfg.seed = 2026;
cfg.dataset_names = [ ...
    "SAR_Dataset_Bangkok_1", "SAR_Dataset_city1_histeq", ...
    "SAR_Dataset_city2_histeq", "SAR_Dataset_SAR_figure", ...
    "SAR_Dataset_filed", "SAR_Dataset_port", "SAR_Dataset_suburb"];
cfg.num_samples_per_dataset = 10;
cfg.raw_decimation_step = 3;
cfg.nyquist_margin = 0.98;

cfg.enl_dataset_names = ["SAR_Dataset_filed", "SAR_Dataset_port"];
cfg.enl_window_size = 64;
cfg.enl_stride = 16;
cfg.entropy_num_bins = 256;

% 阶段一：过零置零规则快筛配置（默认单向距离/方位 × {ZT, RT, RSFT, ASFT}；2D-SFT 与 BARU 为后续可选）
cfg.screen = struct();
cfg.screen.branches = ["range_branch", "azimuth_branch"];
cfg.screen.range_allocations = [1, 1; 4, 1; 9, 1];
cfg.screen.azimuth_allocations = [1, 1; 1, 4; 1, 9];
cfg.screen.range_threshold_types = ["ZT", "RT", "RSFT"];
cfg.screen.azimuth_threshold_types = ["ZT", "RT", "ASFT"];
cfg.screen.rt_as_values = [0.3, 0.6, 0.9];
cfg.screen.sft_str_offsets = [0, 2];
cfg.screen.phi0 = 0;
cfg.screen.rules = ["B0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "L1", "R8", "R9"];
cfg.screen.decision_modes = ["parallel", "serial"];
cfg.screen.remove_dc_options = [false, true];
cfg.screen.include_b0_dc = true;
cfg.screen.sft_baselines = build_default_sft_baselines();
cfg.screen.use_fast_path = true;

% 运行期内存与 I/O 优化设置（单核执行，不启用并行池；不参与断点计算参数校验）
cfg.runtime = struct("cache_gb", 5, "preload_signals", true);

% 阶段二：单向距离/方位 × 四类阈值独立寻优与评估（自动读取阶段一 Top 5；不含 BARU 与 2D-SFT）
cfg.stage2 = struct();
cfg.stage2.branches = ["range_branch", "azimuth_branch"];
cfg.stage2.range_allocations = [(1:9).', ones(9, 1)];
cfg.stage2.azimuth_allocations = [ones(9, 1), (1:9).'];
cfg.stage2.range_threshold_types = ["ZT", "RT", "RSFT"];
cfg.stage2.azimuth_threshold_types = ["ZT", "RT", "ASFT"];
cfg.stage2.baseline_variant = "B0_par_raw";
cfg.stage2.top_k = 5;
cfg.stage2.distinct_rules_only = false;
cfg.stage2.dedup_equivalent = true;
cfg.stage2.selected_variants = string.empty(0, 1);
cfg.stage2.ranking_file = "";
cfg.stage2.rt_as_grid = round([0.05, 0.1:0.1:1.5], 6);
cfg.stage2.sft_str_grid = [-6, -4, -2, 0, 2, 4];
cfg.stage2.sft_freq_grid = round([0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.2], 6);
cfg.stage2.sft_refine = true;
cfg.stage2.sft_refine_str_step = 1;
cfg.stage2.sft_refine_freq_step = 0.1;
cfg.stage2.phi0 = 0;
cfg.stage2.include_b0_transfer = true;
cfg.stage2.leakage_threshold_ratio = 0.35;
cfg.stage2.use_fast_path = true;
end

function lp = read_local_paths(local_paths_file)
inputs_dir = fileparts(local_paths_file);
prev_dir = pwd;
guard = onCleanup(@() cd(prev_dir));
cd(inputs_dir);
lp = local_paths();
end

function baselines = build_default_sft_baselines()
% 基于历史参考表构造阶段一严格满足 Nyquist 边界的单轴代表性基准参数（RSFT / ASFT）。
entries = { ...
    "range",   1, 1, "RSFT",  4, 0.5, 0.0; ...
    "range",   4, 1, "RSFT", -4, 2.0, 0.0; ...
    "range",   9, 1, "RSFT", -6, 1.6, 0.0; ...
    "azimuth", 1, 1, "ASFT",  4, 0.0, 0.5; ...
    "azimuth", 1, 4, "ASFT", -3, 0.0, 1.8; ...
    "azimuth", 1, 9, "ASFT", -5, 0.0, 1.5};

n = size(entries, 1);
baselines = repmat(struct( ...
    "Axis", "", "Range_q", 0, "Azimuth_q", 0, ...
    "ThresholdType", "", "STR_dB", 0, "FrOverBr", 0, "FaOverBa", 0), n, 1);
for k = 1:n
    baselines(k).Axis = string(entries{k, 1});
    baselines(k).Range_q = entries{k, 2};
    baselines(k).Azimuth_q = entries{k, 3};
    baselines(k).ThresholdType = string(entries{k, 4});
    baselines(k).STR_dB = entries{k, 5};
    baselines(k).FrOverBr = entries{k, 6};
    baselines(k).FaOverBa = entries{k, 7};
end
end
