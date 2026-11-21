%% ΔSSIM / ΔPSNR vs Q（测试集 + 校准集）
% 测试集：paired_wilcoxon_vs_b0.csv 中的 vs_B0 独立重校准行；p<0.05 标实心。
% 校准集：由 search_best_params.csv 的逐单元最优值减去同单元 B0 最优值。
% 输出目录：results_zc1bit/figures/curves/{test,calibration}/delta_metric_vs_q/

script_dir = string(fileparts(mfilename("fullpath")));
addpath(char(script_dir));
addpath(char(fileparts(script_dir)));

[cfg, p, S60, record] = zc1bit_script_begin(); %#ok<ASGLU>
D = zc1bit_load_stage2_results(cfg, p, record);
style = zc1bit_plot_style();
variants = D.pseudo_keys;

%% 测试集 Δ
test_table = D.wilcoxon_ind;
test_dir = fullfile(p.figures, "curves", "test", "delta_metric_vs_q");
zc1bit_plot_delta_vs_q(test_table, "DeltaSSIM_Mean", "P_SSIM", ...
    "Test ΔSSIM", test_dir, variants, style, ...
    struct("metric_label", "SSIM", "title_prefix", "Test ΔSSIM"));
zc1bit_plot_delta_vs_q(test_table, "DeltaPSNR_Mean", "P_PSNR", ...
    "Test ΔPSNR", test_dir, variants, style, ...
    struct("metric_label", "PSNR (dB)", "title_prefix", "Test ΔPSNR"));
fprintf("已输出测试集 Δ 指标曲线：%s\n", test_dir);

%% 校准集 Δ
cal_table = zc1bit_build_calibration_delta(D.search_best, D.baseline_key);
cal_dir = fullfile(p.figures, "curves", "calibration", "delta_metric_vs_q");
zc1bit_plot_delta_vs_q(cal_table, "DeltaSSIM_Mean", "", ...
    "Calibration ΔSSIM", cal_dir, variants, style, ...
    struct("metric_label", "SSIM", "title_prefix", "Calibration ΔSSIM"));
zc1bit_plot_delta_vs_q(cal_table, "DeltaPSNR_Mean", "", ...
    "Calibration ΔPSNR", cal_dir, variants, style, ...
    struct("metric_label", "PSNR (dB)", "title_prefix", "Calibration ΔPSNR"));
fprintf("已输出校准集 Δ 指标曲线：%s\n", cal_dir);
