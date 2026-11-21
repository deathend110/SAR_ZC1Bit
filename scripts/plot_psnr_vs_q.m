%% PSNR vs Q 曲线（测试集 + 校准集，两线图 + 总览图）
% 直接运行本脚本即可；需要先完成阶段二 search -> confirm -> evaluate。
% 输出目录：results_zc1bit/figures/curves/{test,calibration}/psnr_vs_q/

script_dir = string(fileparts(mfilename("fullpath")));
addpath(char(script_dir));
addpath(char(fileparts(script_dir)));

[cfg, p, S60, record] = zc1bit_script_begin(); %#ok<ASGLU>
D = zc1bit_load_stage2_results(cfg, p, record);
style = zc1bit_plot_style();
variants = D.pseudo_keys;

test_table = zc1bit_attach_pvalues(D.eval_ind, D.wilcoxon, ...
    "P_PSNR", "vs_B0", "independent_recalibrated");
test_dir = fullfile(p.figures, "curves", "test", "psnr_vs_q");
zc1bit_plot_metric_vs_q(test_table, "PSNR_Mean", "PSNR_Std", "P_PSNR", ...
    "Test PSNR", test_dir, variants, D.baseline_key, style, ...
    struct("metric_label", "PSNR (dB)", "show_two_line", true, "show_overview", true));
fprintf("已输出测试集 PSNR-vs-Q 曲线：%s\n", test_dir);

cal_table = D.search_best;
cal_table.PSNR_Mean = cal_table.Calibration_PSNR_Mean;
if ismember("Calibration_PSNR_Std", string(cal_table.Properties.VariableNames))
    cal_table.PSNR_Std = cal_table.Calibration_PSNR_Std;
else
    cal_table.PSNR_Std = nan(height(cal_table), 1);
end
cal_dir = fullfile(p.figures, "curves", "calibration", "psnr_vs_q");
zc1bit_plot_metric_vs_q(cal_table, "PSNR_Mean", "PSNR_Std", "", ...
    "Calibration PSNR", cal_dir, variants, D.baseline_key, style, ...
    struct("metric_label", "PSNR (dB)", "show_two_line", true, "show_overview", true));
fprintf("已输出校准集 PSNR-vs-Q 曲线：%s\n", cal_dir);
