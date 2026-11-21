%% 最优阈值参数 vs Q（校准集）
% RT: As(Q)；RSFT: STR(Q), fr/Br(Q)；ASFT: STR(Q), fa/Ba(Q)。
% 每个伪规则一条实线，B0 黑色虚线；触界点为空心方形。
% 输出目录：results_zc1bit/figures/curves/ablation/optimal_params_vs_q/

script_dir = string(fileparts(mfilename("fullpath")));
addpath(char(script_dir));
addpath(char(fileparts(script_dir)));

[cfg, p, S60, record] = zc1bit_script_begin(); %#ok<ASGLU>
D = zc1bit_load_stage2_results(cfg, p, record);
style = zc1bit_plot_style();
variants = D.pseudo_keys;
out_dir = fullfile(p.figures, "curves", "ablation", "optimal_params_vs_q");
zc1bit_ensure_dir(out_dir);
S = D.search_best;

%% RT: As(Q)
T = S(S.ThresholdType == "RT", :);
branches = unique(T.Branch, "stable");
if ~isempty(branches)
    fig = figure("Visible", "off", "Color", "w", ...
        "Position", [100, 100, 1200, 520]);
    tiledlayout(fig, 1, numel(branches), "TileSpacing", "compact", "Padding", "compact");
    for b_idx = 1:numel(branches)
        nexttile;
        plot_parameter_panel(T, branches(b_idx), "As", "A_s", ...
            variants, D.baseline_key, style);
    end
    exportgraphics(fig, fullfile(out_dir, "optimal_rt_as.png"), ...
        "Resolution", style.export_resolution);
    close(fig);
    fprintf("已输出 RT 最优 As：%s\n", fullfile(out_dir, "optimal_rt_as.png"));
end

%% RSFT: STR(Q), fr/Br(Q)
T = S(S.ThresholdType == "RSFT", :);
if ~isempty(T)
    fig = figure("Visible", "off", "Color", "w", ...
        "Position", [100, 100, 1200, 520]);
    tiledlayout(fig, 1, 2, "TileSpacing", "compact", "Padding", "compact");
    nexttile;
    plot_parameter_panel(T, "range_branch", "STR_dB", "STR (dB)", ...
        variants, D.baseline_key, style);
    nexttile;
    plot_parameter_panel(T, "range_branch", "FrOverBr", "f_r / B_r", ...
        variants, D.baseline_key, style);
    exportgraphics(fig, fullfile(out_dir, "optimal_rsft_params.png"), ...
        "Resolution", style.export_resolution);
    close(fig);
    fprintf("已输出 RSFT 最优参数：%s\n", fullfile(out_dir, "optimal_rsft_params.png"));
end

%% ASFT: STR(Q), fa/Ba(Q)
T = S(S.ThresholdType == "ASFT", :);
if ~isempty(T)
    fig = figure("Visible", "off", "Color", "w", ...
        "Position", [100, 100, 1200, 520]);
    tiledlayout(fig, 1, 2, "TileSpacing", "compact", "Padding", "compact");
    nexttile;
    plot_parameter_panel(T, "azimuth_branch", "STR_dB", "STR (dB)", ...
        variants, D.baseline_key, style);
    nexttile;
    plot_parameter_panel(T, "azimuth_branch", "FaOverBa", "f_a / B_a", ...
        variants, D.baseline_key, style);
    exportgraphics(fig, fullfile(out_dir, "optimal_asft_params.png"), ...
        "Resolution", style.export_resolution);
    close(fig);
    fprintf("已输出 ASFT 最优参数：%s\n", fullfile(out_dir, "optimal_asft_params.png"));
end

function plot_parameter_panel(T, branch_key, param_col, y_label, variants, baseline_key, style)
sub = T(T.Branch == string(branch_key), :);
if isempty(sub)
    title(sprintf("missing %s", branch_key), "Interpreter", "none");
    return;
end
vars = string(T.Properties.VariableNames);
has_boundary = ismember("TouchesAnyBoundary", vars);
hold on;
grid on;
box on;
h = gobjects(numel(variants) + 1, 1);
names = strings(numel(variants) + 1, 1);

base = sub(sub.VariantKey == baseline_key, :);
base = sortrows(base, "Q");
base_y = base.(param_col);
h(1) = plot(base.Q, base_y, style.baseline_style, ...
    "Color", style.baseline_color, "LineWidth", style.line_width + 0.2);
names(1) = "B0";
for v_idx = 1:numel(variants)
    vkey = variants(v_idx);
    vs = sub(sub.VariantKey == vkey, :);
    vs = sortrows(vs, "Q");
    color = style.pseudo_colors(mod(v_idx - 1, size(style.pseudo_colors, 1)) + 1, :);
    var_y = vs.(param_col);
    h(v_idx + 1) = plot(vs.Q, var_y, "-o", ...
        "Color", color, "LineWidth", style.line_width, ...
        "MarkerSize", style.marker_size, "MarkerFaceColor", "w");
    names(v_idx + 1) = vkey;
    if has_boundary
        b = logical(vs.TouchesAnyBoundary);
        if any(b)
            plot(vs.Q(b), var_y(b), "s", ...
                "Color", color, "MarkerSize", style.sig_marker_size + 2, ...
                "LineWidth", 1.2, "HandleVisibility", "off");
        end
    end
end

title(sprintf("%s | %s", string(branch_key), param_col), ...
    "Interpreter", "none", "FontSize", style.title_size);
xlabel("Upsampling factor Q", "FontSize", style.label_size);
ylabel(y_label, "FontSize", style.label_size, "Interpreter", "none");
xlim([1, 9]);
xticks(1:9);
set(gca, "FontName", style.font_name, "FontSize", style.font_size);
legend(h, names, "Location", "eastoutside", "Interpreter", "none");
end
