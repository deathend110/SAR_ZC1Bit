# SAR-ZC1Bit 阶段二成像与绘图脚本

本目录脚本用于阶段二 `search → confirm → evaluate` 完成后的结果可视化。每个脚本都可以在 MATLAB 编辑器中单独打开并点击“运行”，不依赖工作区变量，也不需要先运行其他脚本。

## 前置条件

请先完成阶段二，并确认以下文件存在：

```text
results_zc1bit/search/search_selected_variants.csv
results_zc1bit/search/search_best_params.csv
results_zc1bit/locked/confirm_results.mat
results_zc1bit/locked/locked_confirm_plan.csv
results_zc1bit/evaluate/evaluate_summary.csv
results_zc1bit/evaluate/evaluate_detail.csv
results_zc1bit/evaluate/evaluate_leakage_summary.csv   # 可选
results_zc1bit/stats/paired_wilcoxon_vs_b0.csv
results_zc1bit/protocol/run_record.mat
```

缺失时脚本会给出中文错误提示。所有输出写到 `results_zc1bit/figures/`，不会写入 C 盘临时目录。

## 脚本一览

| 脚本 | 作用 | 主要输出 |
| --- | --- | --- |
| `plot_ssim_vs_q.m` | SSIM vs Q；测试集 + 校准集；两线图 + 总览图 | `figures/curves/{test,calibration}/ssim_vs_q/` |
| `plot_psnr_vs_q.m` | PSNR vs Q；同样两套来源与布局 | `figures/curves/{test,calibration}/psnr_vs_q/` |
| `plot_delta_metric_vs_q.m` | ΔSSIM / ΔPSNR vs Q，测试集标 p<0.05 | `figures/curves/{test,calibration}/delta_metric_vs_q/` |
| `plot_enl_vs_q.m` | ENL 与 \|ENL−ENL_GT\| vs Q；测试集 + 校准集各 10 个 filed/port 样本（以结果列实际值为准） | `figures/curves/{test,calibration}/enl_vs_q/` |
| `plot_mechanism_vs_q.m` | ZeroRatio / DC bias vs Q（测试集 + 校准集）；诊断组泄漏柱状图（测试集） | `figures/curves/{test,calibration}/mechanism_vs_q/` |
| `plot_optimal_params_vs_q.m` | RT As、RSFT/ASFT STR 与频率 vs Q，标触界点 | `figures/curves/ablation/optimal_params_vs_q/` |
| `plot_final_gain_bars.m` | Q>1 或 Q=9 的平均 ΔSSIM/ΔPSNR 柱状图 | `figures/curves/{test,calibration}/final_gain_bars/` |
| `plot_transfer_ablation.m` | 独立重校准 − B0 参数直迁（测试集） | `figures/curves/ablation/transfer_vs_recalibration/` |
| `render_best_unit_images.m` | 每个方法各自测试集最优单元的 SAR 图像 | `figures/images/best_unit/` |
| `render_common_unit_images.m` | 固定共同单元、同一样本的 7 方法对比 | `figures/images/common_unit/` |
| `plot_ta_sft_results.m` | TA-SFT 专题曲线、相偏鲁棒性与 SAR 拼图 | `ta_sft/figures/` |

共享辅助函数（不需要直接运行）：

- `zc1bit_script_begin.m`：路径注入、读取 `run_record`/`S60`、创建输出目录。
- `zc1bit_load_stage2_results.m`：读取并校验阶段二全部表格。
- `zc1bit_plot_metric_vs_q.m`：指标 vs Q 的两线/总览绘图。
- `zc1bit_plot_delta_vs_q.m`：Δ 指标 vs Q 绘图。
- `zc1bit_render_locked_image.m`：按锁参精确重渲染单样本图像。
- `zc1bit_select_render_sample.m`：选择代表样本并加载回波/GT。
- `zc1bit_save_image_montage.m`：GT + 7 方法图像/误差图拼图。

## 常用可调参数

在每个脚本顶部：

- `sample_id = []`：成像脚本默认取测试集 B0 跨单元平均 SSIM 的中位数样本；填入具体 `SampleID` 可覆盖。
- `unit_keys = string.empty(0, 1)`（仅 `render_common_unit_images.m`）：为空时每个 branch 自动选 B0 测试 SSIM 最高的单元；也可指定 `["rb_R9A1_range_rsft"; "ab_R1A9_azimuth_asft"]`。
- `q_mode = "upsampled"`（仅 `plot_final_gain_bars.m`）：`"upsampled"` 用 Q>1，`"q9"` 仅用 Q=9。

## 数据口径

- 测试集曲线：`evaluate_summary.csv` 中 `Track == "independent_recalibrated"` 的行；参数来自校准集锁参，**不在测试集上重新选优**。
- 校准集曲线：`search_best_params.csv` 的 `Calibration_SSIM_Mean / Calibration_PSNR_Mean`。
- 基线：`VariantKey == "B0_par_raw"`，所有图统一黑色虚线。
- 显著性：测试集使用 `paired_wilcoxon_vs_b0.csv` 的 `P_SSIM / P_PSNR`；校准集没有逐样本搜索明细，因此不标 p 值。
- `Q=1` 是无上采样对照点，通常会出现明显掉点；图不做裁剪。
- “各自最优单元”图像展示每个方法的最佳测试单元，不是同一单元对比；受控对比请看 `render_common_unit_images.m`。

## 图像重渲染一致性

`zc1bit_render_locked_image.m` 复用了 evaluate 的同一套函数，并通过 `src/workflow/zc1bit_stage2_rt_seed.m` 复现 RT 随机阈值。渲染完成后会把单样本 PSNR/SSIM 与 `evaluate_detail.csv` 对比：

- SSIM 绝对误差应小于 `1e-8`；
- PSNR 绝对误差应小于 `1e-6`。

如果出现 `zc1bit:PlotMetricMismatch` 警告，优先检查 `ParamKey`、`SampleID` 或阶段二结果是否已过期。
