%% 严格伪 1.5-bit（{+1,0,-1}）SFT 专题实验绘图与 SAR 成像对比脚本
% 直接运行即可；需要先完成 run_zc1bit_sft_p15(cfg, "all")。
% 全部输出保存至 results_zc1bit/sft_p15/figures/，绝不修改阶段一/二与 ta_sft 已有图表。

script_dir = string(fileparts(mfilename("fullpath")));
addpath(char(script_dir));
addpath(char(fileparts(script_dir)));

[cfg, p, S60, record] = zc1bit_script_begin();
eval_mat = fullfile(p.sft_p15, "evaluate", "sft_p15_evaluate_results.mat");
if ~isfile(eval_mat)
    error("zc1bit:MissingSftP15EvalResults", ...
        "未找到 SFT-P15 评估结果：%s。请先运行 run_zc1bit_sft_p15(cfg, ""all"")。", eval_mat);
end
E = load(eval_mat, "units", "b0_sft_params", "method_order", "locked_params", ...
    "detail_table", "summary_table", "leakage_table", "wilcoxon_table", "robustness_table");

style = zc1bit_plot_style();
fig_root = fullfile(p.sft_p15, "figures");
curves_dir = fullfile(fig_root, "curves");
images_dir = fullfile(fig_root, "images");
tables_dir = fullfile(fig_root, "tables");
zc1bit_ensure_dir(fig_root);
zc1bit_ensure_dir(curves_dir);
zc1bit_ensure_dir(images_dir);
zc1bit_ensure_dir(tables_dir);

method_order = string(E.method_order(:));
method_styles = build_sft_p15_styles(method_order);

% 1. SSIM & PSNR vs Q（RSFT / ASFT 总览，标出 vs B0 统计显著正增益）
plot_metric_overview(E.summary_table, E.wilcoxon_table, "RSFT", "SSIM", ...
    method_order, method_styles, style, fullfile(curves_dir, "sft_p15_ssim_vs_q_rsft.png"));
plot_metric_overview(E.summary_table, E.wilcoxon_table, "ASFT", "SSIM", ...
    method_order, method_styles, style, fullfile(curves_dir, "sft_p15_ssim_vs_q_asft.png"));
plot_metric_overview(E.summary_table, E.wilcoxon_table, "RSFT", "PSNR", ...
    method_order, method_styles, style, fullfile(curves_dir, "sft_p15_psnr_vs_q_rsft.png"));
plot_metric_overview(E.summary_table, E.wilcoxon_table, "ASFT", "PSNR", ...
    method_order, method_styles, style, fullfile(curves_dir, "sft_p15_psnr_vs_q_asft.png"));

% 2. Delta SSIM & Delta PSNR vs B0 曲线
plot_delta_vs_b0(E.wilcoxon_table, "RSFT", method_order, method_styles, style, ...
    fullfile(curves_dir, "sft_p15_delta_vs_b0_rsft.png"));
plot_delta_vs_b0(E.wilcoxon_table, "ASFT", method_order, method_styles, style, ...
    fullfile(curves_dir, "sft_p15_delta_vs_b0_asft.png"));

% 3. 锁定超参随 Q 演化
plot_optimal_hyperparams(E.summary_table, style, ...
    fullfile(curves_dir, "sft_p15_optimal_hyperparams.png"));

% 4. 带外频谱泄漏 + SFT 载波基频抑制比柱状图
if ~isempty(E.leakage_table) && height(E.leakage_table) > 0
    plot_leakage_and_carrier_bars(E.leakage_table, method_order, style, ...
        fullfile(curves_dir, "sft_p15_leakage_and_carrier_bars.png"));
end

% 5. 初相失配鲁棒性曲线（突出全盲方案 2 的水平抗相偏优势）
if isfield(E, "robustness_table") && ~isempty(E.robustness_table) && height(E.robustness_table) > 0
    plot_phase_robustness(E.robustness_table, method_order, method_styles, style, "SSIM", ...
        fullfile(curves_dir, "sft_p15_phase_robustness_ssim.png"));
    plot_phase_robustness(E.robustness_table, method_order, method_styles, style, "PSNR", ...
        fullfile(curves_dir, "sft_p15_phase_robustness_psnr.png"));
end

% 6. 典型样本（测试集 B0_par_raw SSIM 中位数样本）SAR 成像与误差拼图
render_metrics = render_sft_p15_montages(cfg, S60, record, E, method_order, style, images_dir);
if ~isempty(render_metrics) && height(render_metrics) > 0
    writetable(render_metrics, fullfile(tables_dir, "sft_p15_render_metrics.csv"));
end

fprintf("严格伪 1.5-bit SFT 专题图表已输出至：%s\n", fig_root);

function st = build_sft_p15_styles(method_order)
n_m = numel(method_order);
st = containers.Map("KeyType", "char", "ValueType", "any");
palette = [ ...
    0.00, 0.00, 0.00; ... % B0_par_raw: 黑虚线
    0.45, 0.45, 0.45; ... % R1_par_raw: 灰点划线
    0.85, 0.33, 0.10; ... % P15_1_carrier_agree: 橙红
    0.00, 0.45, 0.74; ... % P15_2_rlp_she: 深蓝（全盲）
    0.47, 0.67, 0.19; ... % P15_3_arcsin_lin: 绿
    0.49, 0.18, 0.56];    % P15_4_dual_domain: 紫
line_styles = ["--", "-.", "-", "-", "-", "-"];
markers = ["none", "s", "o", "^", "d", "p"];

for k = 1:n_m
    mkey = method_order(k);
    idx = min(k, size(palette, 1));
    st(char(mkey)) = struct( ...
        "Color", palette(idx, :), ...
        "LineStyle", line_styles(idx), ...
        "Marker", markers(idx), ...
        "LineWidth", 1.8 + 0.3 * double(mkey == "B0_par_raw"));
end
end

function plot_metric_overview(sum_tbl, wx_tbl, ttype, metric_name, method_order, mstyles, style, out_path)
sub = sum_tbl(sum_tbl.ThresholdType == ttype, :);
if isempty(sub)
    return;
end
mean_col = metric_name + "_Mean";
p_col = "P_" + metric_name;
delta_col = "Delta" + metric_name + "_Mean";
wx_sub = wx_tbl(wx_tbl.ComparisonType == "vs_B0" & wx_tbl.ThresholdType == ttype, :);

fig = figure("Visible", "off", "Color", "w", "Position", style.figure_position);
hold on;
grid on;
box on;

handles = gobjects(0, 1);
labels = strings(0, 1);
for k = 1:numel(method_order)
    mkey = method_order(k);
    ms = sortrows(sub(sub.MethodKey == mkey, :), "Q");
    if isempty(ms)
        continue;
    end
    st = mstyles(char(mkey));
    h = plot(ms.Q, ms.(mean_col), ...
        "LineStyle", st.LineStyle, "Color", st.Color, ...
        "LineWidth", st.LineWidth, "Marker", st.Marker, ...
        "MarkerSize", style.marker_size, "MarkerFaceColor", "w");
    handles(end + 1, 1) = h; %#ok<AGROW>
    labels(end + 1, 1) = mkey; %#ok<AGROW>

    if mkey ~= "B0_par_raw" && ~isempty(wx_sub)
        wm = sortrows(wx_sub(wx_sub.MethodKey == mkey, :), "Q");
        if ~isempty(wm)
            sig_mask = (wm.(delta_col) > 0) & (wm.(p_col) < 0.05);
            if any(sig_mask)
                sig_q = wm.Q(sig_mask);
                [tf, loc] = ismember(sig_q, ms.Q);
                loc = loc(tf);
                if ~isempty(loc)
                    plot(ms.Q(loc), ms.(mean_col)(loc), "o", ...
                        "Color", st.Color, "MarkerFaceColor", st.Color, ...
                        "MarkerSize", style.sig_marker_size, "HandleVisibility", "off");
                end
            end
        end
    end
end

xlabel("Upsampling factor Q", "FontSize", style.label_size);
if metric_name == "PSNR"
    ylabel("PSNR (dB)", "FontSize", style.label_size);
else
    ylabel(metric_name, "FontSize", style.label_size);
end
title(sprintf("Strict P1.5-bit SFT Test %s | %s (filled = sig > B0, p<0.05)", metric_name, ttype), ...
    "Interpreter", "none", "FontSize", style.title_size);
xticks(unique(sort(sub.Q)));
if ~isempty(handles)
    legend(handles, cellstr(labels), "Interpreter", "none", "Location", "eastoutside");
end
exportgraphics(fig, char(out_path), "Resolution", style.export_resolution);
close(fig);
end

function plot_delta_vs_b0(wx_tbl, ttype, method_order, mstyles, style, out_path)
sub = wx_tbl(wx_tbl.ComparisonType == "vs_B0" & wx_tbl.ThresholdType == ttype, :);
if isempty(sub)
    return;
end
fig = figure("Visible", "off", "Color", "w", "Position", [100, 100, 980, 740]);
tl = tiledlayout(fig, 2, 1, "TileSpacing", "compact", "Padding", "compact");

non_b0 = method_order(method_order ~= "B0_par_raw");
for panel = 1:2
    ax = nexttile(tl);
    hold(ax, "on");
    grid(ax, "on");
    box(ax, "on");
    yline(ax, 0, "--k", "LineWidth", 1.2, "HandleVisibility", "off");
    if panel == 1
        dcol = "DeltaSSIM_Mean";
        pcol = "P_SSIM";
        ylab = "\DeltaSSIM vs B0";
    else
        dcol = "DeltaPSNR_Mean";
        pcol = "P_PSNR";
        ylab = "\DeltaPSNR vs B0 (dB)";
    end

    handles = gobjects(0, 1);
    labels = strings(0, 1);
    for k = 1:numel(non_b0)
        mkey = non_b0(k);
        wm = sortrows(sub(sub.MethodKey == mkey, :), "Q");
        if isempty(wm)
            continue;
        end
        st = mstyles(char(mkey));
        h = plot(ax, wm.Q, wm.(dcol), ...
            "LineStyle", st.LineStyle, "Color", st.Color, ...
            "LineWidth", st.LineWidth, "Marker", st.Marker, ...
            "MarkerSize", style.marker_size, "MarkerFaceColor", "w");
        handles(end + 1, 1) = h; %#ok<AGROW>
        labels(end + 1, 1) = mkey; %#ok<AGROW>
        sig = (wm.(dcol) > 0) & (wm.(pcol) < 0.05);
        if any(sig)
            plot(ax, wm.Q(sig), wm.(dcol)(sig), "o", ...
                "Color", st.Color, "MarkerFaceColor", st.Color, ...
                "MarkerSize", style.sig_marker_size, "HandleVisibility", "off");
        end
    end
    ylabel(ax, ylab, "FontSize", style.label_size);
    xticks(ax, unique(sort(sub.Q)));
    if panel == 1
        title(ax, sprintf("Strict P1.5-bit SFT Gain vs B0 | %s", ttype), ...
            "Interpreter", "none", "FontSize", style.title_size);
        if ~isempty(handles)
            legend(ax, handles, cellstr(labels), "Interpreter", "none", "Location", "eastoutside");
        end
    else
        xlabel(ax, "Upsampling factor Q", "FontSize", style.label_size);
    end
end

exportgraphics(fig, char(out_path), "Resolution", style.export_resolution);
close(fig);
end

function plot_optimal_hyperparams(sum_tbl, style, out_path)
fig = figure("Visible", "off", "Color", "w", "Position", [100, 100, 1120, 760]);
tl = tiledlayout(fig, 2, 2, "TileSpacing", "compact", "Padding", "compact");
q_ticks = unique(sort(sum_tbl.Q));

% 子图 1：P15_1 Tau
ax1 = nexttile(tl);
hold(ax1, "on"); grid(ax1, "on"); box(ax1, "on");
for ttype = ["RSFT", "ASFT"]
    sub = sortrows(sum_tbl(sum_tbl.MethodKey == "P15_1_carrier_agree" & sum_tbl.ThresholdType == ttype, :), "Q");
    if ~isempty(sub)
        plot(ax1, sub.Q, sub.Tau, "-o", "LineWidth", 1.8, "DisplayName", char(ttype));
    end
end
title(ax1, "P15_1 Carrier-Agree \tau", "FontSize", style.title_size);
xlabel(ax1, "Q", "FontSize", style.label_size);
ylabel(ax1, "\tau (1 = B0)", "FontSize", style.label_size);
xticks(ax1, q_ticks);
ylim(ax1, [0, 1.08]);
legend(ax1, "Location", "best", "Interpreter", "none");

% 子图 2：P15_2 MinRunLength (Inf 映射为 8) 与 HalfWidth
ax2 = nexttile(tl);
hold(ax2, "on"); grid(ax2, "on"); box(ax2, "on");
inf_cap = 8;
for ttype = ["RSFT", "ASFT"]
    sub = sortrows(sum_tbl(sum_tbl.MethodKey == "P15_2_rlp_she" & sum_tbl.ThresholdType == ttype, :), "Q");
    if isempty(sub)
        continue;
    end
    y_vals = double(sub.MinRunLength);
    y_vals(isinf(y_vals)) = inf_cap;
    plot(ax2, sub.Q, y_vals, "-o", "LineWidth", 1.8, "DisplayName", char(ttype));
    w2_mask = double(sub.HalfWidth) > 1;
    if any(w2_mask)
        plot(ax2, sub.Q(w2_mask), y_vals(w2_mask), "s", "MarkerSize", 8, ...
            "MarkerEdgeColor", "k", "DisplayName", char(ttype + " (w=2)"));
    end
end
title(ax2, "P15_2 Blind RLP-SHE (L_{min}, s = w=2)", "FontSize", style.title_size);
xlabel(ax2, "Q", "FontSize", style.label_size);
ylabel(ax2, "L_{min}", "FontSize", style.label_size);
xticks(ax2, q_ticks);
yticks(ax2, [1:7, inf_cap]);
yticklabels(ax2, ["1", "2", "3", "4", "5", "6", "7", "Inf (B0)"]);
ylim(ax2, [0.5, inf_cap + 0.5]);
legend(ax2, "Location", "best", "Interpreter", "none");

% 子图 3：P15_3 TauLin
ax3 = nexttile(tl);
hold(ax3, "on"); grid(ax3, "on"); box(ax3, "on");
for ttype = ["RSFT", "ASFT"]
    sub = sortrows(sum_tbl(sum_tbl.MethodKey == "P15_3_arcsin_lin" & sum_tbl.ThresholdType == ttype, :), "Q");
    if ~isempty(sub)
        plot(ax3, sub.Q, sub.TauLin, "-o", "LineWidth", 1.8, "DisplayName", char(ttype));
    end
end
title(ax3, "P15_3 Arcsin Linearization \tau_{lin}", "FontSize", style.title_size);
xlabel(ax3, "Q", "FontSize", style.label_size);
ylabel(ax3, "\tau_{lin} (1 = B0)", "FontSize", style.label_size);
xticks(ax3, q_ticks);
ylim(ax3, [0.1, 1.08]);
legend(ax3, "Location", "best", "Interpreter", "none");

% 子图 4：P15_4 联合 Tau（点旁标 Lmin）
ax4 = nexttile(tl);
hold(ax4, "on"); grid(ax4, "on"); box(ax4, "on");
for t_idx = 1:2
    ttype = ["RSFT", "ASFT"];
    ttype = ttype(t_idx);
    sub = sortrows(sum_tbl(sum_tbl.MethodKey == "P15_4_dual_domain" & sum_tbl.ThresholdType == ttype, :), "Q");
    if isempty(sub)
        continue;
    end
    plot(ax4, sub.Q, sub.Tau, "-o", "LineWidth", 1.8, "DisplayName", char(ttype));
    for k = 1:height(sub)
        if isinf(sub.MinRunLength(k))
            lbl = "Inf";
        else
            lbl = sprintf("L%g", sub.MinRunLength(k));
        end
        y_off = 0.03 * (3 - 2 * t_idx);
        text(ax4, sub.Q(k), sub.Tau(k) + y_off, lbl, "FontSize", 8, ...
            "HorizontalAlignment", "center", "Interpreter", "none");
    end
end
title(ax4, "P15_4 Dual-Domain (\tau & L_{min})", "FontSize", style.title_size);
xlabel(ax4, "Q", "FontSize", style.label_size);
ylabel(ax4, "\tau", "FontSize", style.label_size);
xticks(ax4, q_ticks);
ylim(ax4, [0.4, 1.12]);
legend(ax4, "Location", "best", "Interpreter", "none");

exportgraphics(fig, char(out_path), "Resolution", style.export_resolution);
close(fig);
end

function plot_leakage_and_carrier_bars(leak_tbl, method_order, style, out_path)
units_u = unique(leak_tbl(:, ["UnitKey", "GroupName", "ThresholdType"]), "rows", "stable");
if isempty(units_u)
    return;
end
n_u = height(units_u);
n_m = numel(method_order);

metrics = ["OffSupportRatio", "RangeOffRatio", "AzimuthOffRatio", "CarrierSuppression_dB"];
titles = [ ...
    "Total Off-Support Ratio", ...
    "Range Off-Support Ratio", ...
    "Azimuth Off-Support Ratio", ...
    "Carrier Tone Suppression vs B0 (dB)"];

fig = figure("Visible", "off", "Color", "w", "Position", [100, 100, 1280, 760]);
tl = tiledlayout(fig, 2, 2, "TileSpacing", "compact", "Padding", "compact");
x_labels = units_u.GroupName + " " + units_u.ThresholdType;

for p_idx = 1:4
    ax = nexttile(tl);
    mat = nan(n_u, n_m);
    mcol = metrics(p_idx);
    for u = 1:n_u
        for m = 1:n_m
            mask = leak_tbl.UnitKey == units_u.UnitKey(u) & leak_tbl.MethodKey == method_order(m);
            if any(mask)
                mat(u, m) = mean(leak_tbl.(mcol)(mask), "omitnan");
            end
        end
    end
    bar(ax, mat, "grouped");
    grid(ax, "on");
    box(ax, "on");
    xticks(ax, 1:n_u);
    xticklabels(ax, cellstr(x_labels));
    xtickangle(ax, 25);
    if p_idx <= 3
        ylim(ax, [0, 1]);
    else
        yline(ax, 0, "--k", "LineWidth", 1.1, "HandleVisibility", "off");
    end
    title(ax, titles(p_idx), "FontSize", style.title_size);
    if p_idx == 4
        legend(ax, cellstr(method_order), "Interpreter", "none", "Location", "best");
    end
end

exportgraphics(fig, char(out_path), "Resolution", style.export_resolution);
close(fig);
end

function plot_phase_robustness(rob_tbl, method_order, mstyles, style, metric_name, out_path)
units_u = unique(rob_tbl(:, ["UnitKey", "GroupName", "ThresholdType"]), "rows", "stable");
n_u = height(units_u);
if n_u == 0
    return;
end
n_cols = min(2, n_u);
n_rows = ceil(n_u / n_cols);
fig = figure("Visible", "off", "Color", "w", "Position", [100, 100, 560 * n_cols, 400 * n_rows]);
tl = tiledlayout(fig, n_rows, n_cols, "TileSpacing", "compact", "Padding", "compact");

mean_col = metric_name + "_Mean";
non_ref = method_order(~ismember(method_order, ["B0_par_raw", "R1_par_raw"]));

for u = 1:n_u
    ax = nexttile(tl);
    hold(ax, "on");
    grid(ax, "on");
    box(ax, "on");
    sub = rob_tbl(rob_tbl.UnitKey == units_u.UnitKey(u), :);
    angs = unique(sort(sub.PhaseMismatchDeg));

    for ref_key = ["B0_par_raw", "R1_par_raw"]
        rsub = sub(sub.MethodKey == ref_key, :);
        if ~isempty(rsub)
            st = mstyles(char(ref_key));
            yline(ax, mean(rsub.(mean_col)), st.LineStyle, ...
                "Color", st.Color, "LineWidth", st.LineWidth, "DisplayName", char(ref_key));
        end
    end

    for k = 1:numel(non_ref)
        mkey = non_ref(k);
        msub = sortrows(sub(sub.MethodKey == mkey, :), "PhaseMismatchDeg");
        if isempty(msub)
            continue;
        end
        st = mstyles(char(mkey));
        plot(ax, msub.PhaseMismatchDeg, msub.(mean_col), ...
            "LineStyle", st.LineStyle, "Color", st.Color, ...
            "LineWidth", st.LineWidth, "Marker", st.Marker, ...
            "MarkerSize", style.marker_size, "MarkerFaceColor", "w", ...
            "DisplayName", char(mkey));
    end

    title(ax, sprintf("%s (%s)", units_u.GroupName(u), units_u.ThresholdType(u)), ...
        "Interpreter", "none", "FontSize", style.title_size);
    xlabel(ax, "Receiver Phase Mismatch \Delta\phi_0 (deg)", "FontSize", style.label_size);
    ylabel(ax, metric_name, "FontSize", style.label_size);
    xticks(ax, angs);
    if u == n_u
        legend(ax, "Interpreter", "none", "Location", "best");
    end
end

exportgraphics(fig, char(out_path), "Resolution", style.export_resolution);
close(fig);
end

function render_tbl = render_sft_p15_montages(cfg, S60, record, E, method_order, style, images_dir)
render_tbl = table();
if ~isfield(cfg, "data_root") || strlength(string(cfg.data_root)) == 0 || ~isfolder(cfg.data_root)
    warning("zc1bit:MissingDataRootForRender", "未找到原始回波目录，跳过 SFT-P15 SAR 成像拼图渲染。");
    return;
end

det = E.detail_table;
b0_sub = det(det.MethodKey == "B0_par_raw", :);
ids = unique(double(b0_sub.SampleID));
means = nan(numel(ids), 1);
for k = 1:numel(ids)
    means(k) = mean(double(b0_sub.SSIM(double(b0_sub.SampleID) == ids(k))), "omitnan");
end
[~, loc] = min(abs(means - median(means, "omitnan")));
sample_id = ids(loc);

manifest = record.manifest;
for k = 1:height(manifest)
    manifest.FilePath(k) = fullfile(cfg.data_root, manifest.Dataset(k), manifest.File(k));
end
sample_row = manifest(manifest.Split == "test" & double(manifest.SampleID) == sample_id, :);
if height(sample_row) ~= 1 || ~isfile(sample_row.FilePath)
    warning("zc1bit:RenderSampleUnavailable", "SampleID=%d 回波不可用，跳过 SAR 渲染。", sample_id);
    return;
end

signals = zc1bit_stage2_load_signals(sample_row, S60, cfg.raw_decimation_step, "SFT-P15 渲染样本");
sig60 = signals{1};
img_gt = generate_gt_image(sig60, S60);

target_units = ["rb_R9A1_range_rsft", "ab_R1A9_azimuth_asft"];
avail_units = target_units(ismember(target_units, E.units.UnitKey));
if isempty(avail_units)
    avail_units = E.units.UnitKey(1);
end

rows = cell(0, 1);
use_fast = true;
if isfield(cfg.sft_p15, "use_fast_path")
    use_fast = logical(cfg.sft_p15.use_fast_path);
end

for u_idx = 1:numel(avail_units)
    uk = avail_units(u_idx);
    unit = E.units(E.units.UnitKey == uk, :);
    b0_row = E.b0_sft_params(E.b0_sft_params.UnitKey == uk, :);
    rq = double(unit.Range_q);
    aq = double(unit.Azimuth_q);
    sgrid = resolve_sampling_grid(rq, aq, [S60.nrn, S60.nan], S60);
    sig_up = two_dim_upsample_fft(sig60, rq, aq);

    cparams = struct("fr_over_Br", b0_row.FrOverBr, "fa_over_Ba", b0_row.FaOverBa, "phi0", b0_row.Phi0);
    carrier = build_sft_carrier(sgrid, cparams, S60);
    t_spec = struct( ...
        "ThresholdType", string(unit.ThresholdType), ...
        "STR_dB", double(b0_row.STR_dB), ...
        "fr_over_Br", double(b0_row.FrOverBr), ...
        "fa_over_Ba", double(b0_row.FaOverBa), ...
        "phi0", double(b0_row.Phi0));
    U = build_threshold(sig_up, sgrid, t_spec, S60, struct("nyquist_margin", cfg.nyquist_margin, "carrier", carrier));
    channel_1bit = quantize_1bit_with_u(sig_up, U);
    clear U sig_up

    fast_ctx = [];
    if use_fast
        fast_ctx = build_fast_focus_context(S60, rq, aq);
    end

    n_m = numel(method_order);
    img_list = cell(n_m + 1, 1);
    title_list = strings(n_m + 1, 1);
    img_list{1} = img_gt;
    title_list(1) = sprintf("GT (Sample %d)", sample_id);

    for m_idx = 1:n_m
        mkey = method_order(m_idx);
        lr = E.locked_params(E.locked_params.UnitKey == uk & E.locked_params.MethodKey == mkey, :);
        ropts = struct("Axis", string(unit.Axis));
        if ~isnan(lr.Tau), ropts.Tau = double(lr.Tau); end
        if ~isnan(lr.HalfWidth), ropts.HalfWidth = double(lr.HalfWidth); end
        if ~isnan(lr.MinRunLength), ropts.MinRunLength = double(lr.MinRunLength); end
        if ~isnan(lr.TauLin), ropts.TauLin = double(lr.TauLin); end
        ch_rec = apply_sft_p15_reconstruct(channel_1bit, carrier, mkey, ropts);
        if use_fast
            img_out = fast_focus_complex_roi(ch_rec, fast_ctx, rq, aq);
        else
            img_out = focus_upsampled_channel(ch_rec, S60, rq, aq);
        end
        mval = compute_image_metrics(img_out, img_gt);
        img_list{m_idx + 1} = img_out;
        title_list(m_idx + 1) = sprintf("%s (%s)\nPSNR %.2f | SSIM %.4f", ...
            mkey, string(lr.HyperparamTag), mval.PSNR, mval.SSIM);

        rows{end + 1, 1} = table( ...
            double(sample_id), uk, string(unit.ThresholdType), double(unit.Q), ...
            mkey, string(lr.HyperparamTag), mval.PSNR, mval.SSIM, ...
            'VariableNames', {'SampleID', 'UnitKey', 'ThresholdType', 'Q', 'MethodKey', 'HyperparamTag', 'PSNR', 'SSIM'}); %#ok<AGROW>
    end

    zc1bit_save_image_montage(img_list, title_list, ...
        fullfile(images_dir, "sft_p15_" + uk + "_montage.png"), style, "image");
    zc1bit_save_image_montage(img_list, title_list, ...
        fullfile(images_dir, "sft_p15_" + uk + "_error_montage.png"), style, "error");
end

if ~isempty(rows)
    render_tbl = vertcat(rows{:});
end
end
