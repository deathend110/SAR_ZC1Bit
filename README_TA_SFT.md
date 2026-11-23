# SAR-ZC1Bit — TA-SFT 专题使用说明

**阈值感知 SFT 1-bit 接收端重构（Threshold-Aware SFT Reconstruction）**是面向 `SAR-ZC1Bit` 正弦频率阈值（`RSFT` / `ASFT`）分支的独立实验扩展：在保持前端采集与下行传输为 **100% 严格 1-bit（$\pm 1$）** 且绝不接触未量化幅度尺度 $\hat\sigma$ 的前提下，利用接收端已知单位模 SFT 载波，对比 5 个阈值感知重构方法与 `B0_par_raw`、`R1_par_raw` 的成像表现与初相失配鲁棒性。

详细数学推导与实验协议见 [docs/ta_sft_experiment_design.md](docs/ta_sft_experiment_design.md)，前期主实验结果见 [docs/experiment_report.md](docs/experiment_report.md)。

---

## 1. 方法一览（7 组同场对比）

| MethodKey | 类型 | 送入 RC 的数值形式 | 寻优超参 |
| :--- | :--- | :--- | :--- |
| `B0_par_raw` | 原始 1-bit 基线 | $\{+1, -1\}$ | 无 |
| `R1_par_raw` | 盲过零伪 1.5-bit 对照 | $\{+1, 0, -1\}$ | 无 |
| `TA1_gate_zc` | 方案 1：载波门控选择性置零 | **严格三电平** $\{+1, 0, -1\}$ | 门限 $\tau \in [0.05, 1.0]$ |
| `TA2_anchor_zc` | 方案 2（连续）：过零锚点 $-\alpha\bar{u}$ 替换 | 非过零 $\pm 1$ + 过零浮点 | 权重 $\alpha \in [0, 3.0]$ |
| `TA2_anchor_q5` | 方案 2（量化）：先节点平均再对称 5 电平量化 | **严格 5 电平** $\{-1, -0.5, 0, +0.5, +1\}$ | 权重 $\alpha \in [0, 3.0]$ |
| `TA3_mmse_fixed` | 方案 3A：固定标称阈信比逐点 MMSE | 全点连续浮点 | 标称比 $\gamma \in [0.2, 2.6]$ |
| `TA3_mmse_est` | 方案 3B：纯 1-bit 相关盲估计 $\hat\gamma$ + MMSE | 全点连续浮点 | **0 超参**（自适应盲估计） |

---

## 2. 前置条件

1. 已完成本地路径配置 `inputs/local_paths.m`（参考 [inputs/local_paths.example.m](inputs/local_paths.example.m)）。
2. 已完成阶段二 `confirm` 步骤，存在 `results_zc1bit/locked/locked_params.csv`（用于读取 18 个 SFT 单元上的 `B0_par_raw` 锁定前端参数；也可通过 `cfg.ta_sft.baseline_locked_file` 指定其他文件）。

---

## 3. 相关文件与输出目录

### 源码与入口
- 配置与路径：[zc1bit_config.m](zc1bit_config.m)（`cfg.ta_sft`）、[src/workflow/zc1bit_paths.m](src/workflow/zc1bit_paths.m)（`p.ta_sft`）
- 核心重构算子：[src/quantization/apply_ta_sft_reconstruct.m](src/quantization/apply_ta_sft_reconstruct.m)
- 校准寻优与评估引擎：[src/workflow/zc1bit_ta_sft_search_impl.m](src/workflow/zc1bit_ta_sft_search_impl.m)、[src/workflow/zc1bit_ta_sft_evaluate_impl.m](src/workflow/zc1bit_ta_sft_evaluate_impl.m)
- 一键运行与绘图：[run_zc1bit_ta_sft.m](run_zc1bit_ta_sft.m)、[scripts/plot_ta_sft_results.m](scripts/plot_ta_sft_results.m)
- 单元测试：[tests/TestTaSftReconstruction.m](tests/TestTaSftReconstruction.m)

### 隔离产出目录（均不覆盖阶段一/二原有结果）
- 断点缓存：`results_zc1bit/checkpoints/ta_sft/{search,evaluate,robustness}/`（按 `UnitKey` 扁平 `-v7` 原子保存，自动复用阶段二 GT 缓存）
- 校准寻优：`results_zc1bit/ta_sft/search/{ta_sft_search_summary.csv, ta_sft_locked_params.csv, ta_sft_search_results.mat}`
- 测试评估：`results_zc1bit/ta_sft/evaluate/{ta_sft_evaluate_detail.csv, ta_sft_evaluate_summary.csv, ta_sft_leakage_summary.csv, ta_sft_evaluate_results.mat}`
- 显著性检验：`results_zc1bit/ta_sft/stats/ta_sft_paired_wilcoxon.csv`（含 `vs_B0` 与 `vs_R1`）
- 相偏鲁棒性：`results_zc1bit/ta_sft/robustness/ta_sft_phase_robustness_summary.csv`
- 图表与渲染：`results_zc1bit/ta_sft/figures/{curves,images,tables}/`

---

## 4. 快速运行命令（MATLAB）

```matlab
% 1. 代码静态检查
checkcode("zc1bit_config.m", "run_zc1bit_ta_sft.m", ...
    "src/workflow/zc1bit_paths.m", "src/quantization/apply_ta_sft_reconstruct.m", ...
    "src/workflow/zc1bit_ta_sft_search_impl.m", "src/workflow/zc1bit_ta_sft_evaluate_impl.m", ...
    "scripts/plot_ta_sft_results.m", "tests/TestTaSftReconstruction.m")

% 2. 运行全部单元测试（含原阶段一/二测试与 TestTaSftReconstruction）
assertSuccess(runtests("tests"))

% 3. 运行 TA-SFT 全流程（18 单元校准寻优 -> 测试评估、带外泄漏 -> 4 单元初相失配消融）
cfg = zc1bit_config();
cfg.runtime.cache_gb = 5;              % 建议 4–6 GB；修改路径或缓存预算可断点续算
out = run_zc1bit_ta_sft(cfg, "all");   % 也可分步传入 "search" 或 "evaluate"

% 4. 生成全部对比曲线、超参演化、相偏衰减与典型样本 SAR 图像
run("scripts/plot_ta_sft_results.m");
```
