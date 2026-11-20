# inputs 目录说明

- `local_paths.example.m`：Git 跟踪的路径配置模板。初次克隆仓库后，请复制为 `local_paths.m` 并修改为本机路径。
- `local_paths.m`：本地路径配置文件（已被 `.gitignore` 忽略，不提交到 Git）。
- `parameters/`：存放本机 `FS60_params.mat`（已被 `*.mat` 规则忽略）。
- `reference/`：存放从 `SAR_Azimuth_Range_Upsample` 导入的只读参考 CSV（受 `!/inputs/reference/*.csv` 规则跟踪），用于阶段一快筛基准、原参数直迁消融与 70/20 样本/ROI 口径核对，不替代阶段二的独立校准搜索：
  - `v5_sample_manifest_70.csv`：70 样本清单参考表
  - `v5_enl_roi_manifest_20.csv`：`filed` + `port` 共 20 样本的 64×64 ENL 固定 ROI 参考表
  - `v4_rsft2d_group_summary.csv`：2D-SFT 全倍率（$Q=1\text{–}10$）历史锁定参考表
  - `v4_rsft_calibration_summary.csv`：1D RSFT 历史校准参考表
  - `v5_rt_as_summary.csv`：SplitRT 幅度系数 $A_s$ 敏感性历史参考表
