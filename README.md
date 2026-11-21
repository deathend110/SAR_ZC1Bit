# SAR-ZC1Bit

基于纯 1-bit（$\pm 1$）SAR 回波的伪 1.5-bit 过零置零重构（$\{+1, 0, -1\}$）成像研究仓库。实验协议以 [experiment_protocol.md](docs/experiment_protocol.md) 与 [sar_zc1bit_experiment_design.md](docs/sar_zc1bit_experiment_design.md) 为唯一真相源，工程结构对齐 `SAR Dataset Transformer`。

## 目录结构

- 根目录：交互入口（`start_zc1bit.m`）、集中配置（`zc1bit_config.m` / `zc1bit_smoke_config.m`）、路径守卫（`setup_zc1bit_paths.m`）及阶段入口（`run_zc1bit_*.m`）。
- `src/`：独立维护的 MATLAB 源码目录，禁止跨仓库 `addpath` 强耦合：
  - `src/compat/`：`gaussian.m`、`normalize_image.m`
  - `src/sampling/`：距离与方位 FFT 上采样、频谱裁剪及网格解析（底层保留双向 BARU 接口作后续可选扩展）
  - `src/threshold/`：ZT、SplitRT、RSFT、ASFT 四类主阈值（保留 2D-SFT 接口作可选扩展）与 $\hat\sigma$ 估计、Nyquist 校验
  - `src/quantization/`：非对称 $\operatorname{sgn}_+(0)=+1$ 1-bit 量化与伪 1.5-bit 规则族（`B0`、`R1`–`R9`、`L1`）
  - `src/sarproc/`：`Range_Compress`、`RCMC`、`SAR_Imaging`、标准路径与频域直裁快路径
  - `src/metrics/`：PSNR、SSIM、固定 64×64 GT ROI ENL、256-bin Entropy、频谱泄漏比与配对 Wilcoxon 检验
  - `src/workflow/`：路径管理、原子 checkpoint 落盘、仅实际参数核对、70 样本 35/35 划分、阶段一规则快筛引擎及阶段二（`search` / `confirm` / `evaluate`）引擎
- `inputs/`：
  - `local_paths.example.m`：本机数据与参数路径配置模板；本地使用前复制为未跟踪的 `inputs/local_paths.m`。
  - `parameters/`：存放本机 `FS60_params.mat`（Git 忽略）。
  - `reference/`：存放只读历史参考参数表与 70/20 样本/ROI 清单 CSV。
- `tests/`：不依赖外部大数据集的轻量单元测试。
- `docs/`：协议文档、技术手册与 [工作流使用说明](docs/WORKFLOW_USAGE.md)。
- `results_zc1bit/`：实验产物与断点目录（按 `protocol/`、`screen/`、`checkpoints/` 等子目录存放，Git 忽略）。
- `tmp/`：工作区临时目录（`tmp/tests/`、`tmp/smoke/`、`tmp/scratch/`，Git 忽略），单元测试与临时中间产物均保存于工作区内，严禁写入 C 盘。

## 快速开始

1. 复制 `inputs/local_paths.example.m` 为 `inputs/local_paths.m`，确认 `data_root`（如 `G:\MATLAB-G\SAR Full PSF`）与 `parameter_file_60` 指向正确路径。
2. 在 MATLAB 中打开 `start_zc1bit.m`，保持 `cfg.stage = "preflight"` 点击运行，核对预检输出的校准/测试样本数与阶段一任务规模。
3. 确认就绪后将 `cfg.stage` 改为 `"screen"` 运行阶段一规则快筛（默认单核运行、`4–6 GB` 内存缓存单轮读盘，按 `ThresholdKey` 扁平断点续算）；快筛完成后可依次运行 `"search"` → `"confirm"` → `"evaluate"`（或 `"all"`）执行阶段二单向四阈值自动 Top 5 寻优、锁参与测试集评估。
4. 阶段二完成后，可直接打开 `scripts/` 下的独立脚本生成 SSIM/PSNR vs Q 曲线、ENL/机理诊断、最优参数汇总与最优指标 SAR 图像；详见 [scripts/README.md](scripts/README.md)。
5. 运行轻量单元测试：

```matlab
results = runtests("tests");
assertSuccess(results);
```
