%% ENL 与 |ENL - ENL_GT| vs Q（测试集 + 校准集）
% 每个划分中 filed + port 的样本数不同（当前 35/35 划分下各 10 个；
% 以 ENLSampleCount / Calibration_ENL_Mean 实际值为准）。
% 输出目录：results_zc1bit/figures/curves/{test,calibration}/enl_vs_q/

script_dir = string(fileparts(mfilename("fullpath")));
addpath(char(script_dir));
addpath(char(fileparts(script_dir)));

[cfg, p, S60, record] = zc1bit_script_begin(); %#ok<ASGLU>
D = zc1bit_load_stage2_results(cfg, p, record);
style = zc1bit_plot_style();
variants = D.pseudo_keys;

test_dir = fullfile(p.figures, "curves", "test", "enl_vs_q");
gt_enl_ref = mean(double(D.eval_ind.GT_ENL_Mean), "omitnan");
if isfinite(gt_enl_ref)
    enl_ref = struct("metric_label", "ENL", "show_two_line", false, ...
        "show_overview", true, "reference_value", gt_enl_ref, ...
        "reference_label", "GT ENL");
else
    enl_ref = struct("metric_label", "ENL", "show_two_line", false, ...
        "show_overview", true);
end
zc1bit_plot_metric_vs_q(D.eval_ind, "ENL_Mean", "", "", ...
    "Test ENL", test_dir, variants, D.baseline_key, style, enl_ref);

zc1bit_plot_metric_vs_q(D.eval_ind, "AbsENLDiff_Mean", "", "", ...
    "Test |ENL-ENL_GT|", test_dir, variants, D.baseline_key, style, ...
    struct("metric_label", "|ENL - ENL_GT|", "show_two_line", false, ...
    "show_overview", true, "reference_value", 0, "reference_label", "ideal"));
fprintf("已输出测试集 ENL 曲线：%s\n", test_dir);

%% 校准集 ENL（与测试集合并即为 filed + port 的 20 个样本）
cal_table = D.search_best;
cal_table.ENL_Mean = cal_table.Calibration_ENL_Mean;
cal_table.GT_ENL_Mean = cal_table.Calibration_GT_ENL_Mean;
cal_table.AbsENLDiff_Mean = cal_table.Calibration_AbsENLDiff_Mean;
cal_dir = fullfile(p.figures, "curves", "calibration", "enl_vs_q");
gt_enl_cal = mean(double(cal_table.GT_ENL_Mean), "omitnan");
if isfinite(gt_enl_cal)
    enl_cal_ref = struct("metric_label", "ENL", "show_two_line", false, ...
        "show_overview", true, "reference_value", gt_enl_cal, ...
        "reference_label", "GT ENL");
else
    enl_cal_ref = struct("metric_label", "ENL", "show_two_line", false, ...
        "show_overview", true);
end
zc1bit_plot_metric_vs_q(cal_table, "ENL_Mean", "", "", ...
    "Calibration ENL", cal_dir, variants, D.baseline_key, style, enl_cal_ref);
zc1bit_plot_metric_vs_q(cal_table, "AbsENLDiff_Mean", "", "", ...
    "Calibration |ENL-ENL_GT|", cal_dir, variants, D.baseline_key, style, ...
    struct("metric_label", "|ENL - ENL_GT|", "show_two_line", false, ...
    "show_overview", true, "reference_value", 0, "reference_label", "ideal"));
fprintf("已输出校准集 ENL 曲线：%s\n", cal_dir);
