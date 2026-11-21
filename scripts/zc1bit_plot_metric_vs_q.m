function zc1bit_plot_metric_vs_q(T, metric_col, std_col, p_col, source_name, out_dir, variants, baseline_key, style, options)
%ZC1BIT_PLOT_METRIC_VS_Q 绘制指标随上采样倍率 Q 变化的曲线。
% 支持“每个变体与 B0 的两线图”和“每个 branch × threshold 的多线总览图”。
% T 必须包含 Branch, ThresholdType, Q, VariantKey 与 metric_col。

if nargin < 9 || isempty(style)
    style = zc1bit_plot_style();
end
if nargin < 10 || isempty(options)
    options = struct();
end
options = set_default(options, "show_two_line", true);
options = set_default(options, "show_overview", true);
options = set_default(options, "show_errorbar", false);
options = set_default(options, "metric_label", metric_col);
options = set_default(options, "reference_value", NaN);
options = set_default(options, "reference_label", "");
options = set_default(options, "sample_count", NaN);

vars = string(T.Properties.VariableNames);
required = ["Branch", "ThresholdType", "Q", "VariantKey", metric_col];
missing = required(~ismember(required, vars));
if ~isempty(missing)
    error("zc1bit:PlotMissingColumn", "曲线表缺少必要列：%s", join(missing, ", "));
end

variants = string(variants(:));
two_dir = fullfile(out_dir, "two_line");
over_dir = fullfile(out_dir, "overview");
if options.show_two_line, zc1bit_ensure_dir(two_dir); end
if options.show_overview, zc1bit_ensure_dir(over_dir); end

combos = unique(T(:, ["Branch", "ThresholdType"]), "rows", "stable");
for c_idx = 1:height(combos)
    br = combos.Branch(c_idx);
    tt = combos.ThresholdType(c_idx);
    sub = T(T.Branch == br & T.ThresholdType == tt, :);

    base = sub(sub.VariantKey == baseline_key, :);
    base = sortrows(base, "Q");
    if height(base) ~= 9
        error("zc1bit:PlotMissingCurve", ...
            "基线曲线不完整：%s / %s，Q 点数为 %d（预期 9）。", br, tt, height(base));
    end

    y_limits = resolve_y_limits(sub.(metric_col), options.reference_value, options.metric_label);
    branch_tag = regexprep(string(br), "_branch$", "");
    branch_label = branch_tag;
    source_tag = lower(regexprep(string(source_name), "[^A-Za-z0-9]+", "_"));
    source_tag = regexprep(source_tag, "^_+|_+$", "");

    if options.show_two_line
        for v_idx = 1:numel(variants)
            vkey = variants(v_idx);
            vs = sub(sub.VariantKey == vkey, :);
            vs = sortrows(vs, "Q");
            if height(vs) ~= 9
                error("zc1bit:PlotMissingCurve", ...
                    "变体曲线不完整：%s / %s / %s，Q 点数为 %d（预期 9）。", ...
                    br, tt, vkey, height(vs));
            end

            color = style.pseudo_colors(mod(v_idx - 1, size(style.pseudo_colors, 1)) + 1, :);
            fig = figure("Visible", "off", "Color", "w", "Position", style.figure_position);
            hold on;
            grid on;
            box on;

            base_y = base.(metric_col);
            var_y = vs.(metric_col);
            h_base = plot(base.Q, base_y, style.baseline_style, ...
                "Color", style.baseline_color, "LineWidth", style.line_width + 0.2);
            h_var = plot(vs.Q, var_y, "-o", ...
                "Color", color, "LineWidth", style.line_width, ...
                "MarkerSize", style.marker_size, "MarkerFaceColor", "w");
            if options.show_errorbar && ~isempty(std_col) && ismember(std_col, vars)
                var_e = vs.(std_col);
                errorbar(vs.Q, var_y, var_e, ...
                    "Color", color, "LineWidth", 1.0, ...
                    "LineStyle", "none", "Marker", "none", ...
                    "CapSize", 5, "HandleVisibility", "off");
            end

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
            draw_reference_line(1, 9, options);

            title(sprintf("%s | %s | %s | %s", ...
                string(source_name), branch_label, string(tt), vkey), ...
                "Interpreter", "none", "FontSize", style.title_size);
            xlabel("Upsampling factor Q", "FontSize", style.label_size);
            ylabel(options.metric_label, "FontSize", style.label_size, "Interpreter", "none");
            xlim([1, 9]);
            xticks(1:9);
            ylim(y_limits);
            set(gca, "FontName", style.font_name, "FontSize", style.font_size);
            legend([h_base, h_var], ...
                {"B0 (1-bit baseline)", sprintf("%s (pseudo 1.5-bit)", vkey)}, ...
                "Location", "best", "Interpreter", "none");

            file_name = sprintf("%s_%s_%s_%s.png", source_tag, ...
                lower(branch_tag), lower(string(tt)), vkey);
            exportgraphics(fig, fullfile(two_dir, file_name), ...
                "Resolution", style.export_resolution);
            close(fig);
        end
    end

    if options.show_overview
        fig = figure("Visible", "off", "Color", "w", "Position", style.figure_position);
        hold on;
        grid on;
        box on;
        h = gobjects(numel(variants) + 1, 1);
        names = strings(numel(variants) + 1, 1);
        base_y = base.(metric_col);
        h(1) = plot(base.Q, base_y, style.baseline_style, ...
            "Color", style.baseline_color, "LineWidth", style.line_width + 0.2);
        names(1) = "B0 (1-bit baseline)";
        for v_idx = 1:numel(variants)
            vkey = variants(v_idx);
            vs = sub(sub.VariantKey == vkey, :);
            vs = sortrows(vs, "Q");
            if height(vs) ~= 9
                error("zc1bit:PlotMissingCurve", ...
                    "总览曲线不完整：%s / %s / %s，Q 点数为 %d（预期 9）。", ...
                    br, tt, vkey, height(vs));
            end
            color = style.pseudo_colors(mod(v_idx - 1, size(style.pseudo_colors, 1)) + 1, :);
            var_y = vs.(metric_col);
            h(v_idx + 1) = plot(vs.Q, var_y, "-o", ...
                "Color", color, "LineWidth", style.line_width, ...
                "MarkerSize", style.marker_size, "MarkerFaceColor", "w");
            names(v_idx + 1) = vkey;
        end
        draw_reference_line(1, 9, options);

        title(sprintf("%s | %s | %s | overview", ...
            string(source_name), branch_label, string(tt)), ...
            "Interpreter", "none", "FontSize", style.title_size);
        xlabel("Upsampling factor Q", "FontSize", style.label_size);
        ylabel(options.metric_label, "FontSize", style.label_size, "Interpreter", "none");
        xlim([1, 9]);
        xticks(1:9);
        ylim(y_limits);
        set(gca, "FontName", style.font_name, "FontSize", style.font_size);
        legend(h, names, "Location", "eastoutside", "Interpreter", "none");

        file_name = sprintf("%s_%s_%s_overview.png", source_tag, ...
            lower(branch_tag), lower(string(tt)));
        exportgraphics(fig, fullfile(over_dir, file_name), ...
            "Resolution", style.export_resolution);
        close(fig);
    end
end
end

function y_limits = resolve_y_limits(values, reference_value, metric_label)
values = double(values(:));
values = values(isfinite(values));
if isempty(values)
    error("zc1bit:PlotEmptyMetric", "指标列没有有限值：%s", metric_label);
end
if isfinite(reference_value)
    values = [values; double(reference_value)];
end
y_min = min(values);
y_max = max(values);
if y_max <= y_min
    pad = max(abs(y_min), 1) * 0.01;
else
    pad = 0.06 * (y_max - y_min);
end
y_limits = [y_min - pad, y_max + pad];
if string(metric_label) == "SSIM"
    y_limits(1) = max(y_limits(1), 0);
    y_limits(2) = min(y_limits(2), 1);
end
end

function draw_reference_line(x_min, x_max, options)
if ~isfinite(options.reference_value)
    return;
end
plot([x_min, x_max], [options.reference_value, options.reference_value], ":", ...
    "Color", [0.35, 0.35, 0.35], "LineWidth", 1.2, "HandleVisibility", "off");
if strlength(string(options.reference_label)) > 0
    text(x_min + 0.1, options.reference_value, options.reference_label, ...
        "Color", [0.35, 0.35, 0.35], "FontSize", 9, ...
        "VerticalAlignment", "bottom", "Interpreter", "none");
end
end

function options = set_default(options, name, value)
if ~isfield(options, name) || isempty(options.(name))
    options.(name) = value;
end
end
