%% 机理诊断：ZeroRatio / |DC bias| vs Q，以及诊断组泄漏指标
% 输出目录：results_zc1bit/figures/curves/test/mechanism_vs_q/
% 泄漏指标只覆盖 evaluate_impl 选定的诊断组（R1A1/R4A1/R9A1/R1A4/R1A9），
% 图标题会明确标注“diagnostic groups”，不代表全部 54 个单元。

script_dir = string(fileparts(mfilename("fullpath")));
addpath(char(script_dir));
addpath(char(fileparts(script_dir)));

[cfg, p, S60, record] = zc1bit_script_begin(); %#ok<ASGLU>
D = zc1bit_load_stage2_results(cfg, p, record);
style = zc1bit_plot_style();
variants = D.pseudo_keys;
out_dir = fullfile(p.figures, "curves", "test", "mechanism_vs_q");

zc1bit_plot_metric_vs_q(D.eval_ind, "ZeroRatio_Mean", "", "", ...
    "Test ZeroRatio", out_dir, variants, D.baseline_key, style, ...
    struct("metric_label", "Zero ratio", "show_two_line", false, ...
    "show_overview", true, "reference_value", 0, "reference_label", "no zeroing"));
zc1bit_plot_metric_vs_q(D.eval_ind, "AbsMeanBias_Mean", "", "", ...
    "Test |mean bias|", out_dir, variants, D.baseline_key, style, ...
    struct("metric_label", "|mean bias|", "show_two_line", false, ...
    "show_overview", true, "reference_value", 0, "reference_label", "zero bias"));
fprintf("已输出测试集 ZeroRatio / DC bias 曲线：%s\n", out_dir);

%% 校准集 ZeroRatio / DC bias
cal_table = D.search_best;
cal_table.ZeroRatio_Mean = cal_table.Calibration_ZeroRatio_Mean;
cal_table.AbsMeanBias_Mean = cal_table.Calibration_AbsMeanBias_Mean;
cal_dir = fullfile(p.figures, "curves", "calibration", "mechanism_vs_q");
zc1bit_plot_metric_vs_q(cal_table, "ZeroRatio_Mean", "", "", ...
    "Calibration ZeroRatio", cal_dir, variants, D.baseline_key, style, ...
    struct("metric_label", "Zero ratio", "show_two_line", false, ...
    "show_overview", true, "reference_value", 0, "reference_label", "no zeroing"));
zc1bit_plot_metric_vs_q(cal_table, "AbsMeanBias_Mean", "", "", ...
    "Calibration |mean bias|", cal_dir, variants, D.baseline_key, style, ...
    struct("metric_label", "|mean bias|", "show_two_line", false, ...
    "show_overview", true, "reference_value", 0, "reference_label", "zero bias"));
fprintf("已输出校准集 ZeroRatio / DC bias 曲线：%s\n", cal_dir);

%% 诊断组泄漏指标（按单元/变体对诊断样本取均值）
L = D.leakage;
if isempty(L) || height(L) == 0
    fprintf("evaluate_leakage_summary.csv 为空，跳过泄漏指标图。\n");
    return;
end
vars = string(L.Properties.VariableNames);
if ismember("Track", vars)
    L = L(L.Track == "independent_recalibrated", :);
end
if isempty(L) || height(L) == 0
    fprintf("泄漏表没有 independent_recalibrated 行，跳过泄漏指标图。\n");
    return;
end

[found, loc] = ismember(L.UnitKey, D.units.UnitKey);
if any(~found)
    L = L(found, :);
    loc = loc(found);
end
L.Q = D.units.Q(loc);
L.Axis = D.units.Axis(loc);
L.GroupName = string(L.GroupName);
if isempty(L) || height(L) == 0
    fprintf("泄漏表与 units 无匹配行，跳过泄漏指标图。\n");
    return;
end

unit_keys = unique(L.UnitKey, "stable");
variant_keys = [D.baseline_key; D.pseudo_keys];
metric_names = ["OffSupportRatio", "RangeOffRatio", "AzimuthOffRatio"];
metric_labels = ["off-support ratio", "range off-support", "azimuth off-support"];
available = metric_names(ismember(metric_names, string(L.Properties.VariableNames)));
labels_available = metric_labels(ismember(metric_names, string(L.Properties.VariableNames)));
if isempty(available)
    fprintf("泄漏表缺少 OffSupportRatio/RangeOffRatio/AzimuthOffRatio，跳过。\n");
    return;
end

n_u = numel(unit_keys);
n_m = numel(available);
n_v = numel(variant_keys);
values = nan(n_u, n_v, n_m);
for m_idx = 1:n_m
    col = available(m_idx);
    col_vals = double(L.(col));
    for u_idx = 1:n_u
        ukey = unit_keys(u_idx);
        for v_idx = 1:n_v
            vkey = variant_keys(v_idx);
            mask = L.UnitKey == ukey & L.VariantKey == vkey;
            if any(mask)
                values(u_idx, v_idx, m_idx) = mean(col_vals(mask), "omitnan");
            end
        end
    end
end

labels = strings(n_u, 1);
for u_idx = 1:n_u
    mask = L.UnitKey == unit_keys(u_idx);
    first = find(mask, 1);
    labels(u_idx) = sprintf("%s_%s", ...
        L.GroupName(first), L.Axis(first));
end

fig = figure("Visible", "off", "Color", "w", ...
    "Position", [100, 100, 1200, 760]);
tiledlayout(fig, 1, n_m, "TileSpacing", "compact", "Padding", "compact");
for m_idx = 1:n_m
    nexttile;
    b = bar(values(:, :, m_idx), "grouped");
    for v_idx = 1:n_v
        b(v_idx).FaceColor = style.pseudo_colors( ...
            mod(v_idx - 1, size(style.pseudo_colors, 1)) + 1, :);
        b(v_idx).EdgeColor = "none";
    end
    grid on;
    box on;
    ylim([0, 1]);
    title(sprintf("Leakage (%s, diagnostic groups)", labels_available(m_idx)), ...
        "Interpreter", "none", "FontSize", style.title_size);
    xlabel("Diagnostic unit", "FontSize", style.label_size);
    ylabel(metric_labels(m_idx), "FontSize", style.label_size, "Interpreter", "none");
    xticks(1:n_u);
    xticklabels(labels);
    xtickangle(30);
    set(gca, "FontName", style.font_name, "FontSize", style.font_size - 1, ...
        "TickLabelInterpreter", "none");
end
legend(variant_keys, "Location", "eastoutside", "Interpreter", "none");
exportgraphics(fig, fullfile(out_dir, "leakage_diagnostic_bars.png"), ...
    "Resolution", style.export_resolution);
close(fig);
fprintf("已输出诊断组泄漏指标：%s\n", fullfile(out_dir, "leakage_diagnostic_bars.png"));
