%% SSIM vs Q 曲线（测试集 + 校准集，两线图 + 总览图）
% 直接运行本脚本即可；需要先完成阶段二 search -> confirm -> evaluate。
% 输出目录：results_zc1bit/figures/curves/{test,calibration}/ssim_vs_q/

script_dir = string(fileparts(mfilename("fullpath")));
addpath(char(script_dir));
addpath(char(fileparts(script_dir)));

[cfg, p, S60, record] = zc1bit_script_begin(); %#ok<ASGLU>
D = zc1bit_load_stage2_results(cfg, p, record);
style = zc1bit_plot_style();
variants = D.pseudo_keys;

% 测试集：锁参后的 35 样本 SSIM；p 值来自配对 Wilcoxon。
test_table = zc1bit_attach_pvalues(D.eval_ind, D.wilcoxon, ...
    "P_SSIM", "vs_B0", "independent_recalibrated");
test_dir = fullfile(p.figures, "curves", "test", "ssim_vs_q");
zc1bit_plot_metric_vs_q(test_table, "SSIM_Mean", "SSIM_Std", "P_SSIM", ...
    "Test SSIM", test_dir, variants, D.baseline_key, style, ...
    struct("metric_label", "SSIM", "show_two_line", true, "show_overview", true));
fprintf("已输出测试集 SSIM-vs-Q 曲线：%s\n", test_dir);

% 校准集：search_best_params 中该单元最优参数对应的校准 SSIM。
cal_table = D.search_best;
cal_table.SSIM_Mean = cal_table.Calibration_SSIM_Mean;
if ismember("Calibration_SSIM_Std", string(cal_table.Properties.VariableNames))
    cal_table.SSIM_Std = cal_table.Calibration_SSIM_Std;
else
    cal_table.SSIM_Std = nan(height(cal_table), 1);
end
cal_dir = fullfile(p.figures, "curves", "calibration", "ssim_vs_q");
zc1bit_plot_metric_vs_q(cal_table, "SSIM_Mean", "SSIM_Std", "", ...
    "Calibration SSIM", cal_dir, variants, D.baseline_key, style, ...
    struct("metric_label", "SSIM", "show_two_line", true, "show_overview", true));
fprintf("已输出校准集 SSIM-vs-Q 曲线：%s\n", cal_dir);
