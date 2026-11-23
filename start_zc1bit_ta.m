% 1. 代码静态检查
checkcode("zc1bit_config.m", "run_zc1bit_ta_sft.m", ...
    "src/workflow/zc1bit_paths.m", "src/quantization/apply_ta_sft_reconstruct.m", ...
    "src/workflow/zc1bit_ta_sft_search_impl.m", "src/workflow/zc1bit_ta_sft_evaluate_impl.m", ...
    "scripts/plot_ta_sft_results.m", "tests/TestTaSftReconstruction.m")

% 2. 运行全部单元测试（含原阶段一/二测试与 TestTaSftReconstruction）
assertSuccess(runtests("tests"))

% 3. 运行 TA-SFT 全流程（18 单元校准寻优 -> 测试评估、带外泄漏 -> 4 单元初相失配消融）
cfg = zc1bit_config();
cfg.runtime.cache_gb = 6;              % 建议 4–6 GB；修改路径或缓存预算可断点续算
out = run_zc1bit_ta_sft(cfg, "all");   % 也可分步传入 "search" 或 "evaluate"

% 4. 生成全部对比曲线、超参演化、相偏衰减与典型样本 SAR 图像
run("scripts/plot_ta_sft_results.m");