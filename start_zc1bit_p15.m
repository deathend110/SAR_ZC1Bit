% 1. 代码静态检查
checkcode("zc1bit_config.m", "run_zc1bit_sft_p15.m", ...
    "src/workflow/zc1bit_paths.m", "src/quantization/apply_sft_p15_reconstruct.m", ...
    "src/workflow/zc1bit_sft_p15_search_impl.m", "src/workflow/zc1bit_sft_p15_evaluate_impl.m", ...
    "scripts/plot_sft_p15_results.m", "tests/TestSftP15Reconstruction.m")

% 2. 运行全部单元测试（含原阶段一/二、ta_sft 与 TestSftP15Reconstruction）
assertSuccess(runtests("tests"))

% 3. 运行严格伪 1.5-bit SFT 全流程（18 单元校准寻优 -> 测试评估、载波抑制 -> 4 单元初相失配消融）
cfg = zc1bit_config();
cfg.runtime.cache_gb = 5;              % 建议 4–6 GB；修改路径或缓存预算可断点续算
out = run_zc1bit_sft_p15(cfg, "all");  % 也可分步传入 "search" 或 "evaluate"

% 4. 生成全部对比曲线、超参演化、载波抑制柱状图、相偏衰减与典型样本 SAR 图像
run("scripts/plot_sft_p15_results.m");