%% SAR-ZC1Bit 一键交互式启动脚本
% 在 MATLAB 编辑器中打开本文件，确认下方参数区与待执行阶段序列后点击“运行”。
% 说明：
%   1. 默认一阶段已完成，依次执行 ["search", "confirm", "evaluate"]；
%      若尚未运行或需重跑一阶段，可改为 ["screen", "search", "confirm", "evaluate"]（已有断点自动复用）。
%   2. 阶段二默认采用跨 5 大过零机理族去重口径（SSIMSigWinCount -> SSIMWinCount -> Q>1 过采样 MeanSSIM），
%      自动选出覆盖全部 10 条规则机理的 Top 6 变体（R1_par/ser + R9 + R3 + L1 + R7）与基线 B0_par_raw 对比。
clc;

script_directory = string(fileparts(mfilename("fullpath")));
addpath(script_directory);
cfg = zc1bit_config();
setup_zc1bit_paths(cfg);

%% ============================= 参数区 ==================================
% 按需编辑待执行阶段序列（可选："preflight", "screen", "search", "confirm", "evaluate"）
stages_to_run = ["search", "confirm", "evaluate"];

% 单核运行（不使用并行池），通过 4–6 GB 内存缓存减少机械盘读盘与同倍率上采样 FFT 次数。
cfg.runtime.cache_gb = 6;
cfg.runtime.preload_signals = true;
cfg.screen.use_fast_path = true;
cfg.stage2.use_fast_path = true;

% 阶段二自动跨机理族 Top 6 口径参数（按需取消注释微调）：
cfg.stage2.top_k = 6;
cfg.stage2.dedup_equivalent = true;
cfg.stage2.family_diversified = true;
% cfg.stage2.distinct_rules_only = false;
% cfg.stage2.sft_refine = true;
% cfg.stage2.selected_variants = ["R1_par_raw", "R1_ser_raw", "R9_par_raw", "R3_par_dc", "L1_par_raw", "R7_par_dc"];

%% ============================ 启动流程 =================================
fprintf("回波根目录：%s\n", cfg.data_root);
fprintf("成像参数文件：%s\n", cfg.parameter_file_60);
fprintf("结果输出目录：%s\n", cfg.results_root);

if any(ismember(stages_to_run, ["search", "confirm", "evaluate", "all"]))
    preview_tbl = zc1bit_resolve_stage2_variants(cfg);
    fprintf("\n=== 阶段二已自动锁定参评变体（B0 + Top %d）===\n", cfg.stage2.top_k);
    disp(preview_tbl(:, ["VariantOrder", "VariantRole", "VariantKey", ...
        "Stage1SSIMSigWinCount", "Stage1SSIMWinCount", "Stage1UpsampledSSIM", "Stage1MeanSSIM"]));
end

stage_results = struct();
for s_idx = 1:numel(stages_to_run)
    stage_name = stages_to_run(s_idx);
    cfg.stage = stage_name;
    fprintf("\nSAR-ZC1Bit 执行阶段：%s\n", cfg.stage);
    stage_timer = tic;
    stage_results.(stage_name) = run_zc1bit_workflow(cfg);
    fprintf("阶段 %s 执行完成，用时 %.2f 分钟。\n", stage_name, toc(stage_timer) / 60);
end
