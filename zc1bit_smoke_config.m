function cfg = zc1bit_smoke_config(fixture_root, parameter_file)
%ZC1BIT_SMOKE_CONFIG 轻量冒烟闭环配置（单场景、少量规则与阈值档位，默认 preflight）。
% 用于测试与快速闭环校验，不依赖外部 70 样本大数据集。

cfg = zc1bit_config();
cfg.data_root = string(fixture_root);
cfg.results_root = fullfile(string(fixture_root), "results_zc1bit");
cfg.parameter_file_60 = string(parameter_file);
cfg.dataset_names = "SmokeScene";
cfg.enl_dataset_names = "SmokeScene";
cfg.num_samples_per_dataset = 2;
cfg.raw_decimation_step = 1;
cfg.enl_window_size = 4;
cfg.enl_stride = 2;

cfg.screen.branches = ["range_branch", "azimuth_branch"];
cfg.screen.range_allocations = [1, 1; 2, 1];
cfg.screen.azimuth_allocations = [1, 1; 1, 2];
cfg.screen.range_threshold_types = ["ZT", "RT", "RSFT"];
cfg.screen.azimuth_threshold_types = ["ZT", "ASFT", "2D-SFT"];
cfg.screen.rt_as_values = 0.6;
cfg.screen.sft_str_offsets = 0;
cfg.screen.rules = ["B0", "R1", "R4", "L1"];
cfg.screen.decision_modes = ["parallel", "serial"];
cfg.screen.remove_dc_options = [false, true];
cfg.screen.include_b0_dc = false;
cfg.screen.sft_baselines = build_smoke_sft_baselines();
cfg.stage = "preflight";
end

function baselines = build_smoke_sft_baselines()
entries = { ...
    "range",   1, 1, "RSFT",   0, 0.2, 0.0; ...
    "range",   2, 1, "RSFT",   0, 0.2, 0.0; ...
    "azimuth", 1, 1, "ASFT",   0, 0.0, 0.2; ...
    "azimuth", 1, 1, "2D-SFT", 0, 0.2, 0.2; ...
    "azimuth", 1, 2, "ASFT",   0, 0.0, 0.2; ...
    "azimuth", 1, 2, "2D-SFT", 0, 0.2, 0.2};
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
