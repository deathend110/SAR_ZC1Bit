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
else
    cfg.data_root = "";
    cfg.parameter_file_60 = fullfile(repo_root, "inputs", "parameters", "FS60_params.mat");
    cfg.results_root = fullfile(repo_root, "results_zc1bit");
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

% 阶段一：过零置零规则快筛配置
cfg.screen = struct();
cfg.screen.branches = ["range_branch", "azimuth_branch"];
cfg.screen.range_allocations = [1, 1; 4, 1; 9, 1];
cfg.screen.azimuth_allocations = [1, 1; 1, 4; 1, 9];
cfg.screen.range_threshold_types = ["ZT", "RT", "RSFT", "2D-SFT"];
cfg.screen.azimuth_threshold_types = ["ZT", "RT", "ASFT", "2D-SFT"];
cfg.screen.rt_as_values = [0.3, 0.6, 0.9];
cfg.screen.sft_str_offsets = [0, 2];
cfg.screen.phi0 = 0;
cfg.screen.rules = ["B0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "L1", "R8", "R9"];
cfg.screen.decision_modes = ["parallel", "serial"];
cfg.screen.remove_dc_options = [false, true];
cfg.screen.include_b0_dc = true;
cfg.screen.sft_baselines = build_default_sft_baselines();
cfg.screen.use_fast_path = true;
cfg.screen.workers = 0;
end

function lp = read_local_paths(local_paths_file)
inputs_dir = fileparts(local_paths_file);
prev_dir = pwd;
guard = onCleanup(@() cd(prev_dir));
cd(inputs_dir);
lp = local_paths();
end

function baselines = build_default_sft_baselines()
% 基于历史参考表（V4 RSFT / 2D-SFT）构造阶段一严格满足 Nyquist 边界的代表性基准参数。
entries = { ...
    "range",   1, 1, "RSFT",   4,  0.5, 0.0; ...
    "range",   1, 1, "2D-SFT", 4,  0.5, 0.5; ...
    "range",   4, 1, "RSFT",  -4,  2.0, 0.0; ...
    "range",   4, 1, "2D-SFT",-3,  2.0, 0.5; ...
    "range",   9, 1, "RSFT",  -6,  1.6, 0.0; ...
    "range",   9, 1, "2D-SFT",-5,  1.6, 0.5; ...
    "azimuth", 1, 1, "ASFT",   4,  0.0, 0.5; ...
    "azimuth", 1, 1, "2D-SFT", 4,  0.5, 0.5; ...
    "azimuth", 1, 4, "ASFT",  -3,  0.0, 1.8; ...
    "azimuth", 1, 4, "2D-SFT",-3,  0.5, 1.8; ...
    "azimuth", 1, 9, "ASFT",  -5,  0.0, 1.5; ...
    "azimuth", 1, 9, "2D-SFT",-5,  0.5, 1.5};

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
