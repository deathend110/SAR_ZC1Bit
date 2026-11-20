# SAR-ZC1Bit 工作流使用说明

## 1. 路径配置与环境准备

1. 本仓库通过 `inputs/local_paths.m` 解耦本机路径。初次使用请复制 [local_paths.example.m](../inputs/local_paths.example.m) 为 `inputs/local_paths.m`（已被 `.gitignore` 忽略），配置：
   - `paths.data_root`：外部 180 MHz 原始回波目录（例如 `G:\MATLAB-G\SAR Full PSF`）。
   - `paths.parameter_file_60`：`FS60_params.mat` 文件路径（默认指向 `inputs/parameters/FS60_params.mat`）。
   - `paths.results_root`：结果、断点与报告输出目录（默认 `results_zc1bit/`）。
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
cfg.screen.workers = 0; % 可设为正整数以启用 thread parpool，失败时自动退化串行
out = run_zc1bit_workflow(cfg); % 或直接调用 run_zc1bit_screen(cfg)
```

- **数据隔离**：阶段一**仅在 35 个校准样本**上执行，35 个测试样本保持完全隔离。
- **实验覆盖**：
  - 距离单向分支（dim 1）：`R1A1`、`R4A1`、`R9A1` × `{ZT, RT, RSFT, 2D-SFT}`；
  - 方位单向分支（dim 2）：`R1A1`、`R1A4`、`R1A9` × `{ZT, RT, ASFT, 2D-SFT}`；
  - 候选规则：`B0` + `R1`–`R9` + `L1`；
  - 正交开关：并行判定（`parallel`）vs. 串行原地判定（`serial`）、原始输出（`raw`）vs. 去 DC（`dc_removed`）。
- **断点恢复**：按“阈值配置 × 场景”原子保存于 `results_zc1bit/checkpoints/screen/<ThresholdKey>/<SceneKey>.mat`。恢复时**仅核对实际计算参数**，修改路径、worker 数或重启 MATLAB 均可直接复用已完成进度。
- **阶段一产出**（`results_zc1bit/screen/`）：
  - `screen_detail.csv`：逐样本、逐任务明细（含 PSNR、SSIM、Entropy、ENL、`|ENL - GT_ENL|`、ZeroRatio、I/Q 均值偏置）。
  - `screen_summary.csv`：按任务汇总均值/标准差及相对同阈值 `B0` 的配对差值与双侧 Wilcoxon $p$-value。
  - `screen_rule_ranking.csv`：按全局平均 SSIM 降序（PSNR 降序破局）对各规则与开关组合排序，用于筛选进入阶段二的 Top 2–3 优胜规则。
  - `screen_cal_enl_roi_manifest.csv`：校准集中 `filed` + `port` 样本的 GT 自动 64×64 均匀区 ROI 坐标。
  - `screen_results.mat` 与 `screen_metadata.txt`：完整结构化存档与运行元数据。

### 2.3 阶段二入口（预留）

`run_zc1bit_search`、`run_zc1bit_confirm`、`run_zc1bit_evaluate` 已作为根目录公共入口预留；在阶段一选定优胜规则并按协议授权实现阶段二前，调用这些阶段会显式抛出 `zc1bit:StageNotImplemented`。

## 3. 轻量单元测试

`tests/` 下的所有测试不依赖外部 `SAR Full PSF` 大数据集，使用小型合成信号与虚拟临时目录验证核心数学性质、快速成像等价性及断点恢复逻辑：

```matlab
results = runtests("tests");
assertSuccess(results);
```
