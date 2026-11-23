# SAR-ZC1Bit — 严格伪 1.5-bit SFT 重构专题使用说明（SFT-P15）

**严格伪 1.5-bit SFT 重构（Strict Pseudo-1.5-Bit SFT Reconstruction，`SFT-P15`）**是面向 `SAR-ZC1Bit` 正弦频率阈值（`RSFT` / `ASFT`）分支的独立专题实验扩展：在保持前端采集与下行传输为 **100% 严格 1-bit（$\pm 1$）**、且送入距离压缩（RC）前的输出强制严格属于 **三电平伪 1.5-bit $\{+1, 0, -1\}$**（绝不接触未量化回波 `signal_up` 与幅度尺度 $\hat\sigma$）的前提下，针对 SFT 脉宽调制载波特性对比 4 套改进方案与 `B0_par_raw`、`R1_par_raw`，探索如何在严守 $\{+1, 0, -1\}$ 约束下系统性超越 `B0` + SFT。

详细数学推导、防泄漏红线与实验协议见 [docs/sft_p15_experiment_design.md](docs/sft_p15_experiment_design.md)，前期主实验结果见 [docs/experiment_report.md](docs/experiment_report.md)。

---

## 1. 方法一览（6 组同场对比）

| MethodKey | 类型与机理 | 送入 RC 的数值集合 | 需 $\phi_0$/载波 | 寻优超参网格（共 68 候选/单元） |
| :--- | :--- | :---: | :---: | :--- |
| `B0_par_raw` | 原始 1-bit 基线 | $\{+1, -1\}$ | 否 | 无（1 点） |
| `R1_par_raw` | 原版盲过零伪 1.5-bit 对照 | **严格三电平** $\{+1, 0, -1\}$ | 否 | 无（1 点） |
| `P15_1_carrier_agree` | **方案 1**：载波置信度门控置零（$b\cdot u > \tau \to 0$） | **严格三电平** $\{+1, 0, -1\}$ | 是 | $\tau \in \{0.1\text{–}0.8, 0.85, 0.9, 0.95, 1.0\}$（12 点） |
| `P15_2_rlp_she` | **方案 2**：游程保护 + 可调半宽对称置零（全盲 SHE-ZC） | **严格三电平** $\{+1, 0, -1\}$ | **否（全盲）** | $(w=1, L_{\min}\in\{1\text{–}6,\infty\}) + (w=2, L_{\min}\in\{5,6,7\})$（10 点） |
| `P15_3_arcsin_lin` | **方案 3**：$\arcsin$ 反翘曲单侧置零 | **严格三电平** $\{+1, 0, -1\}$ | 是 | $\tau_{\mathrm{lin}} \in \{0.2, 0.3, \dots, 1.0\}$（9 点） |
| `P15_4_dual_domain` | **方案 4**：波峰置信度 + 游程保护边沿双域联合置零 | **严格三电平** $\{+1, 0, -1\}$ | 是 | 二维 $\tau\in\{0.5,0.7,0.8,0.85,0.9,0.95,1.0\} \times L_{\min}\in\{2,3,4,5,\infty\}$（35 点） |

---

## 2. 前置条件

1. 已完成本地路径配置 `inputs/local_paths.m`（参考 [inputs/local_paths.example.m](inputs/local_paths.example.m)）。
2. 已完成阶段二 `confirm` 步骤，存在 `results_zc1bit/locked/locked_params.csv`（用于读取 18 个 SFT 单元上的 `B0_par_raw` 锁定前端 SFT 参数；也可通过 `cfg.sft_p15.baseline_locked_file` 指定其他锁参文件）。

---

## 3. 相关文件与输出目录

### 源码与入口
- 配置与路径：[zc1bit_config.m](zc1bit_config.m)（`cfg.sft_p15`）、[src/workflow/zc1bit_paths.m](src/workflow/zc1bit_paths.m)（`p.sft_p15`）
- 核心重构算子：[src/quantization/apply_sft_p15_reconstruct.m](src/quantization/apply_sft_p15_reconstruct.m)
- 校准寻优与评估引擎：[src/workflow/zc1bit_sft_p15_search_impl.m](src/workflow/zc1bit_sft_p15_search_impl.m)、[src/workflow/zc1bit_sft_p15_evaluate_impl.m](src/workflow/zc1bit_sft_p15_evaluate_impl.m)
- 一键运行与绘图：[run_zc1bit_sft_p15.m](run_zc1bit_sft_p15.m)、[scripts/plot_sft_p15_results.m](scripts/plot_sft_p15_results.m)
- 单元测试：[tests/TestSftP15Reconstruction.m](tests/TestSftP15Reconstruction.m)

### 隔离产出目录（均不覆盖阶段一/二与 `ta_sft` 原有结果）
- 断点缓存：`results_zc1bit/checkpoints/sft_p15/{search,evaluate,robustness}/`（按 `UnitKey` 扁平 `-v7` 原子保存，自动复用阶段二 GT 缓存）
- 校准寻优：`results_zc1bit/sft_p15/search/{sft_p15_search_summary.csv, sft_p15_locked_params.csv, sft_p15_search_results.mat, sft_p15_search_metadata.txt}`
- 测试评估与诊断：`results_zc1bit/sft_p15/evaluate/{sft_p15_evaluate_detail.csv, sft_p15_evaluate_summary.csv, sft_p15_leakage_summary.csv, sft_p15_test_enl_roi_manifest.csv, sft_p15_evaluate_results.mat, sft_p15_evaluate_metadata.txt}`（含精确频域载波基频抑制比 `CarrierToneRatio` / `CarrierSuppression_dB`）
- 显著性检验：`results_zc1bit/sft_p15/stats/sft_p15_paired_wilcoxon.csv`（含 `vs_B0` 与 `vs_R1`）
- 相偏鲁棒性：`results_zc1bit/sft_p15/robustness/sft_p15_phase_robustness_summary.csv`
- 图表与渲染：`results_zc1bit/sft_p15/figures/{curves,images,tables}/`（含 `tables/sft_p15_render_metrics.csv`）

---

## 4. 快速运行命令（MATLAB）

```matlab
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
```
