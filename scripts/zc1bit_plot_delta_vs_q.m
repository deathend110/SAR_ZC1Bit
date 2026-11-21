function zc1bit_plot_delta_vs_q(T, metric_col, p_col, source_name, out_dir, variants, style, options)
%ZC1BIT_PLOT_DELTA_VS_Q 绘制相对 B0 的 Δ指标随 Q 变化的总览曲线。
% T 必须包含 Branch, ThresholdType, Q, VariantKey 与 metric_col（Δ 指标）。
% 零线代表 B0 基线；测试集可对 p<0.05 的点加实心标记。

if nargin < 7 || isempty(style)
    style = zc1bit_plot_style();
end
if nargin < 8 || isempty(options)
    options = struct();
end
if ~isfield(options, "metric_label") || isempty(options.metric_label)
    options.metric_label = metric_col;
end
if ~isfield(options, "title_prefix") || isempty(options.title_prefix)
    options.title_prefix = string(source_name);
end

vars = string(T.Properties.VariableNames);
required = ["Branch", "ThresholdType", "Q", "VariantKey", metric_col];
missing = required(~ismember(required, vars));
if ~isempty(missing)
    error("zc1bit:PlotMissingColumn", "Δ曲线表缺少必要列：%s", join(missing, ", "));
end

variants = string(variants(:));
zc1bit_ensure_dir(out_dir);
combos = unique(T(:, ["Branch", "ThresholdType"]), "rows", "stable");

for c_idx = 1:height(combos)
    br = combos.Branch(c_idx);
    tt = combos.ThresholdType(c_idx);
    sub = T(T.Branch == br & T.ThresholdType == tt, :);
    branch_tag = regexprep(string(br), "_branch$", "");
    source_tag = lower(regexprep(string(source_name), "[^A-Za-z0-9]+", "_"));
    source_tag = regexprep(source_tag, "^_+|_+$", "");

    vals = double(sub.(metric_col));
    vals = vals(isfinite(vals));
    if isempty(vals)
        error("zc1bit:PlotEmptyMetric", "Δ指标没有有限值：%s / %s", br, tt);
    end
    max_abs = max(abs([vals; 0]));
    if max_abs == 0
        max_abs = 1;
    end
    pad = 0.10 * max_abs;

    fig = figure("Visible", "off", "Color", "w", "Position", style.figure_position);
    hold on;
    grid on;
    box on;

    h_zero = plot([1, 9], [0, 0], style.baseline_style, ...
        "Color", style.baseline_color, "LineWidth", style.line_width + 0.2);

    h = gobjects(numel(variants), 1);
    names = strings(numel(variants), 1);
    for v_idx = 1:numel(variants)
        vkey = variants(v_idx);
        vs = sub(sub.VariantKey == vkey, :);
        vs = sortrows(vs, "Q");
        if height(vs) ~= 9
            error("zc1bit:PlotMissingCurve", ...
                "Δ曲线不完整：%s / %s / %s，Q 点数为 %d（预期 9）。", ...
                br, tt, vkey, height(vs));
        end
        color = style.pseudo_colors(mod(v_idx - 1, size(style.pseudo_colors, 1)) + 1, :);
        var_y = vs.(metric_col);
        h(v_idx) = plot(vs.Q, var_y, "-o", ...
            "Color", color, "LineWidth", style.line_width, ...
            "MarkerSize", style.marker_size, "MarkerFaceColor", "w");
        names(v_idx) = vkey;

        if ~isempty(p_col) && ismember(p_col, vars)
            p_vals = vs.(p_col);
            sig = isfinite(p_vals) & p_vals < 0.05;
            if any(sig)
                plot(vs.Q(sig), var_y(sig), "o", ...
                    "Color", color, "MarkerSize", style.sig_marker_size, ...
                    "MarkerFaceColor", color, "MarkerEdgeColor", color, ...
                    "HandleVisibility", "off");
            end
        end
    end

    title(sprintf("%s | %s | %s | Δ%s", ...
        string(options.title_prefix), branch_tag, string(tt), options.metric_label), ...
        "Interpreter", "none", "FontSize", style.title_size);
    xlabel("Upsampling factor Q", "FontSize", style.label_size);
    ylabel("Δ" + string(options.metric_label) + " vs B0", "FontSize", style.label_size, "Interpreter", "none");
    xlim([1, 9]);
    xticks(1:9);
    ylim([-max_abs - pad, max_abs + pad]);
    set(gca, "FontName", style.font_name, "FontSize", style.font_size);
    legend([h_zero; h], ["B0 (zero line)"; names], ...
        "Location", "eastoutside", "Interpreter", "none");

    file_name = sprintf("%s_%s_%s_delta_%s.png", source_tag, ...
        lower(branch_tag), lower(string(tt)), lower(string(metric_col)));
    exportgraphics(fig, fullfile(out_dir, file_name), ...
        "Resolution", style.export_resolution);
    close(fig);
end
end
