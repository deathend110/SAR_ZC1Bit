%% SAR-ZC1Bit 一键交互式启动脚本
% 在 MATLAB 编辑器中打开本文件，修改下方参数区后点击“运行”。
% 默认仅执行 preflight（预检配置与 70 样本清单）；确认无误后将 cfg.stage 改为 "screen" 启动阶段一快筛。
clc;

script_directory = string(fileparts(mfilename("fullpath")));
addpath(script_directory);
cfg = zc1bit_config();

%% ============================= 参数区 ==================================
% preflight：只核验路径、FS60 参数、70 样本 35/35 划分与阶段一任务列表（不执行成像）。
% screen   ：执行阶段一过零置零规则快筛（仅使用 35 个校准样本，支持断点续算）。
cfg.stage = "preflight";

% 阶段一并行线程数（0 表示串行；正整数尝试启用 thread parpool，失败时自动退化串行）。
cfg.screen.workers = 0;
cfg.screen.use_fast_path = true;

% 如需缩小/调整阶段一快筛范围，可取消以下注释按需定制：
% cfg.screen.rules = ["B0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "L1", "R8", "R9"];
% cfg.screen.rt_as_values = [0.3, 0.6, 0.9];
% cfg.screen.sft_str_offsets = [0, 2];

%% ============================ 启动流程 =================================
fprintf("SAR-ZC1Bit 执行阶段：%s\n", cfg.stage);
fprintf("回波根目录：%s\n", cfg.data_root);
fprintf("成像参数文件：%s\n", cfg.parameter_file_60);
fprintf("结果输出目录：%s\n", cfg.results_root);

workflow_timer = tic;
workflow_result = run_zc1bit_workflow(cfg);
fprintf("所选阶段执行完成，用时 %.2f 分钟。\n", toc(workflow_timer) / 60);
