# SAR-ZC1Bit 工作流使用说明

## 1. 路径配置与环境准备

1. 本仓库通过 `inputs/local_paths.m` 解耦本机路径。初次使用请复制 [local_paths.example.m](../inputs/local_paths.example.m) 为 `inputs/local_paths.m`（已被 `.gitignore` 忽略），配置：
   - `paths.data_root`：外部 180 MHz 原始回波目录（例如 `G:\MATLAB-G\SAR Full PSF`）。
   - `paths.parameter_file_60`：`FS60_params.mat` 文件路径（默认指向 `inputs/parameters/FS60_params.mat`）。
   - `paths.results_root`：结果、断点与报告输出目录（默认 `results_zc1bit/`）。
   - `paths.tmp_root`：工作区临时目录（默认 `tmp/`，按 `tmp/tests/`、`tmp/smoke/`、`tmp/scratch/` 分类保存，Git 忽略）。
2. 所有阶段入口均通过 `setup_zc1bit_paths(cfg)` 仅注入本仓库目录及 `src/` 子目录，并持有 `onCleanup` 目录守卫，支持从任意 MATLAB 当前目录调用且在异常退出时自动恢复工作目录。

## 2. 阶段执行说明

### 2.1 预检（`preflight`，默认）

```matlab
cfg = zc1bit_config();
cfg.stage = "preflight";
out = run_zc1bit_workflow(cfg);
```

- 校验本机目录、工具箱（`psnr`、`ssim`）、`FS60_params.mat` 字段及 Nyquist 安全边界。
- 仅通过 `whos` 读取各场景选定 `rstart*.mat` 的元数据头，生成或核对 70 样本（35 奇数校准 / 35 偶数测试）清单并保存 `results_zc1bit/protocol/manifest_70.csv`。
- **不加载大回波矩阵，不执行任何成像计算。**

### 2.2 阶段一：过零置零规则快筛（`screen`）

```matlab
cfg = zc1bit_config();
cfg.stage = "screen";
cfg.runtime.cache_gb = 5;         % 4–6 GB 内存缓存预算（复用单轮读盘与同倍率 signal_up）
cfg.runtime.preload_signals = true;
out = run_zc1bit_workflow(cfg);   % 或直接调用 run_zc1bit_screen(cfg)
```

- **数据隔离**：阶段一**仅在 35 个校准样本**上执行，35 个测试样本保持完全隔离。
- **实验覆盖（36 阈值配置 × 42 规则变体）**：
  - 距离单向分支（dim 1）：`R1A1`、`R4A1`、`R9A1` × `{ZT, RT, RSFT}`；
  - 方位单向分支（dim 2）：`R1A1`、`R1A4`、`R1A9` × `{ZT, RT, ASFT}`；
  - 候选规则：`B0` + `R1`–`R9` + `L1`；
  - 正交开关：并行判定（`parallel`）vs. 串行原地判定（`serial`）、原始输出（`raw`）vs. 去 DC（`dc_removed`）。
  - *注：双向 BARU 与 `2D-SFT` 保留底层源码接口作后续可选扩展，默认不运行。*
- **单核内存优先与断点恢复**：全程保持单核无并行池运行；存在未完成项时对 7 个 `.mat` 仅顺序读 1 轮并按 `cfg.runtime.cache_gb`（推荐 4–6 GB）复用上采样结果；按 `ThresholdKey` 扁平原子保存于 `results_zc1bit/checkpoints/screen/<ThresholdKey>.mat`（默认 `-v7` 格式），全命中时 0 次读盘。恢复时**仅核对实际计算参数**，修改路径、`runtime.cache_gb`、`preload_signals` 或重启 MATLAB 均可直接复用已完成进度。
- **阶段一产出**（`results_zc1bit/screen/`）：
  - `screen_detail.csv`：逐样本、逐任务明细（含 PSNR、SSIM、Entropy、ENL、`|ENL - GT_ENL|`、ZeroRatio、I/Q 均值偏置）。
  - `screen_summary.csv`：按任务汇总均值/标准差及相对同阈值 `B0` 的配对差值与双侧 Wilcoxon $p$-value。
  - `screen_rule_ranking.csv`：按全局平均 SSIM 降序（PSNR 降序破局）对各规则与开关组合排序，用于筛选进入阶段二的 Top 2–3 优胜规则。
  - `screen_cal_enl_roi_manifest.csv`：校准集中 `filed` + `port` 样本的 GT 自动 64×64 均匀区 ROI 坐标。
  - `screen_results.mat` 与 `screen_metadata.txt`：完整结构化存档与运行元数据。

### 2.3 阶段二：自动 Top 5 独立寻优、锁参与测试集评估（`search → confirm → evaluate` / `all`）

阶段一生成 `results_zc1bit/screen/screen_rule_ranking.csv` 后，阶段二默认**自动排除 `B0` 并读取排名前 5 的最优伪 1.5-bit 变体**，与基线 `B0_par_raw` 共 6 个变体在两条单向主线（距离 `R1A1`–`R9A1` × `{ZT, RT, RSFT}` 与方位 `R1A1`–`R1A9` × `{ZT, RT, ASFT}`，共 54 单元，不含 BARU 与 `2D-SFT`）上执行：

1. **参数寻优（`search`）**：仅使用 35 个校准样本各自独立寻优（SSIM 优先、PSNR 破局），产出 `results_zc1bit/search/{search_selected_variants.csv, search_summary.csv, search_best_params.csv, search_results.mat}`。
2. **锁参确认（`confirm`）**：校验 54 单元搜索网格完整后锁定独立最优参数，并生成 `b0_transferred` 原参数直迁消融计划，产出 `results_zc1bit/locked/{locked_params.csv, locked_confirm_plan.csv, confirm_results.mat}`。
3. **测试集评估与统计（`evaluate`）**：严格隔离使用 35 个测试样本评估 PSNR/SSIM/ENL/Entropy/频谱泄漏并执行双侧配对 Wilcoxon 检验，产出 `results_zc1bit/evaluate/` 与 `results_zc1bit/stats/paired_wilcoxon_vs_b0.csv`。

可通过 `cfg.stage = "search"` / `"confirm"` / `"evaluate"` 分步调用（或 `cfg.stage = "all"` 顺序串行复用断点执行），也可按需调整 `cfg.stage2.top_k`、`cfg.stage2.distinct_rules_only`、`cfg.stage2.sft_refine` 与 `cfg.stage2.selected_variants`。

## 3. 轻量单元测试

`tests/` 下的所有测试不依赖外部 `SAR Full PSF` 大数据集，使用小型合成信号并在工作区 `<repo_root>/tmp/tests/` 下创建隔离临时目录（测试退出自动清理，不写入 C 盘），验证核心数学性质、快速成像等价性及断点恢复逻辑：

```matlab
results = runtests("tests");
assertSuccess(results);
```
