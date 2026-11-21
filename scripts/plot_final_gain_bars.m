%% 汇总增益柱状图：Q>1（或 Q=9）平均 ΔSSIM / ΔPSNR
% 测试集与校准集分别输出；每个 branch × threshold 一张双面板图。
% 输出目录：results_zc1bit/figures/curves/{test,calibration}/final_gain_bars/

%% 可调参数
q_mode = "upsampled";   % "upsampled" = Q>1；"q9" = 仅 Q=9

script_dir = string(fileparts(mfilename("fullpath")));
addpath(char(script_dir));
addpath(char(fileparts(script_dir)));

[cfg, p, S60, record] = zc1bit_script_begin(); %#ok<ASGLU>
D = zc1bit_load_stage2_results(cfg, p, record);
style = zc1bit_plot_style();
variants = D.pseudo_keys;

sources = { ...
    struct("name", "test", "table", D.wilcoxon_ind, ...
        "p_ssim", "P_SSIM", "p_psnr", "P_PSNR"), ...
    struct("name", "calibration", ...
        "table", zc1bit_build_calibration_delta(D.search_best, D.baseline_key), ...
        "p_ssim", "", "p_psnr", "")};

for s_idx = 1:numel(sources)
    src = sources{s_idx};
    T = src.table;
    out_dir = fullfile(p.figures, "curves", src.name, "final_gain_bars");
    zc1bit_ensure_dir(out_dir);
    combos = unique(T(:, ["Branch", "ThresholdType"]), "rows", "stable");
    for c_idx = 1:height(combos)
        br = combos.Branch(c_idx);
        tt = combos.ThresholdType(c_idx);
        [vals_ssim, err_ssim] = aggregate_gain(T, br, tt, variants, ...
            "DeltaSSIM_Mean", src.p_ssim, q_mode);
        [vals_psnr, err_psnr] = aggregate_gain(T, br, tt, variants, ...
            "DeltaPSNR_Mean", src.p_psnr, q_mode);

        fig = figure("Visible", "off", "Color", "w", ...
            "Position", [100, 100, 1200, 520]);
        tiledlayout(fig, 1, 2, "TileSpacing", "compact", "Padding", "compact");
        nexttile;
        plot_bar_panel(vals_ssim, err_ssim, variants, "ΔSSIM vs B0", ...
            sprintf("%s | %s", br, tt), style);
        nexttile;
        plot_bar_panel(vals_psnr, err_psnr, variants, "ΔPSNR vs B0 (dB)", ...
            sprintf("%s | %s", br, tt), style);

        file_name = sprintf("%s_%s_%s_gain_bars.png", ...
            src.name, lower(br), lower(string(tt)));
        exportgraphics(fig, fullfile(out_dir, file_name), ...
            "Resolution", style.export_resolution);
        close(fig);
    end
    fprintf("已输出 %s 汇总增益柱状图：%s\n", src.name, out_dir);
end

function [vals, errs] = aggregate_gain(T, branch_key, threshold_key, variants, metric_col, p_col, q_mode)
vals = nan(numel(variants), 1);
errs = zeros(numel(variants), 1);
for v_idx = 1:numel(variants)
    vkey = variants(v_idx);
    sub = T(T.Branch == branch_key & T.ThresholdType == threshold_key & ...
        T.VariantKey == vkey, :);
    if string(q_mode) == "q9"
        sub = sub(sub.Q == 9, :);
    else
        sub = sub(sub.Q > 1, :);
    end
    if isempty(sub)
        continue;
    end
    vals(v_idx) = mean(double(sub.(metric_col)), "omitnan");
    errs(v_idx) = std(double(sub.(metric_col)), "omitnan");
    if ~isfinite(errs(v_idx))
        errs(v_idx) = 0;
    end
end
end

function plot_bar_panel(vals, errs, variants, y_label, title_text, style)
n = numel(variants);
hold on;
grid on;
box on;
for k = 1:n
    if isfinite(vals(k))
        b = bar(k, vals(k), 0.7);
        b.FaceColor = style.pseudo_colors(mod(k - 1, size(style.pseudo_colors, 1)) + 1, :);
        b.EdgeColor = "none";
    end
end
errorbar(1:n, vals, errs, "Color", "k", "LineStyle", "none", "LineWidth", 1.0);
plot([0.4, n + 0.6], [0, 0], "--", "Color", style.baseline_color, "LineWidth", 1.2);

max_abs = max(abs([vals(isfinite(vals)); 0]));
if max_abs == 0
    max_abs = 1;
end
pad = 0.12 * max_abs;
ylim([-max_abs - pad, max_abs + pad]);
xlim([0.4, n + 0.6]);
xticks(1:n);
xticklabels(variants);
xtickangle(25);
set(gca, "TickLabelInterpreter", "none", ...
    "FontName", style.font_name, "FontSize", style.font_size - 1);
ylabel(y_label, "FontSize", style.label_size, "Interpreter", "none");
title(title_text, "Interpreter", "none", "FontSize", style.title_size);
end
