%% TA-SFT 专题实验绘图与 SAR 成像对比脚本
% 直接运行即可；需要先完成 run_zc1bit_ta_sft(cfg, "all")。
% 全部输出保存至 results_zc1bit/ta_sft/figures/，绝不修改阶段一/二已有图表。

script_dir = string(fileparts(mfilename("fullpath")));
addpath(char(script_dir));
addpath(char(fileparts(script_dir)));

[cfg, p, S60, record] = zc1bit_script_begin();
eval_mat = fullfile(p.ta_sft, "evaluate", "ta_sft_evaluate_results.mat");
if ~isfile(eval_mat)
    error("zc1bit:MissingTaSftEvalResults", ...
        "未找到 TA-SFT 评估结果：%s。请先运行 run_zc1bit_ta_sft(cfg, ""all"")。", eval_mat);
end
E = load(eval_mat, "units", "b0_sft_params", "method_order", "locked_params", ...
    "detail_table", "summary_table", "leakage_table", "wilcoxon_table", "robustness_table");

style = zc1bit_plot_style();
fig_root = fullfile(p.ta_sft, "figures");
curves_dir = fullfile(fig_root, "curves");
images_dir = fullfile(fig_root, "images");
tables_dir = fullfile(fig_root, "tables");
zc1bit_ensure_dir(fig_root);
zc1bit_ensure_dir(curves_dir);
zc1bit_ensure_dir(images_dir);
zc1bit_ensure_dir(tables_dir);

method_order = string(E.method_order(:));
method_styles = build_ta_sft_styles(method_order);

% 1. SSIM & PSNR vs Q（RSFT / ASFT 总览，标出 vs B0 统计显著正增益）
plot_metric_overview(E.summary_table, E.wilcoxon_table, "RSFT", "SSIM", ...
    method_order, method_styles, style, fullfile(curves_dir, "ta_sft_ssim_vs_q_rsft.png"));
plot_metric_overview(E.summary_table, E.wilcoxon_table, "ASFT", "SSIM", ...
    method_order, method_styles, style, fullfile(curves_dir, "ta_sft_ssim_vs_q_asft.png"));
plot_metric_overview(E.summary_table, E.wilcoxon_table, "RSFT", "PSNR", ...
    method_order, method_styles, style, fullfile(curves_dir, "ta_sft_psnr_vs_q_rsft.png"));
plot_metric_overview(E.summary_table, E.wilcoxon_table, "ASFT", "PSNR", ...
    method_order, method_styles, style, fullfile(curves_dir, "ta_sft_psnr_vs_q_asft.png"));

% 2. Delta SSIM & Delta PSNR vs B0 曲线
plot_delta_vs_b0(E.wilcoxon_table, "RSFT", method_order, method_styles, style, ...
    fullfile(curves_dir, "ta_sft_delta_vs_b0_rsft.png"));
plot_delta_vs_b0(E.wilcoxon_table, "ASFT", method_order, method_styles, style, ...
    fullfile(curves_dir, "ta_sft_delta_vs_b0_asft.png"));

% 3. 锁定超参（tau, alpha, gamma 与盲估计 gamma_hat）随 Q 变化
plot_optimal_hyperparams(E.summary_table, style, ...
    fullfile(curves_dir, "ta_sft_optimal_hyperparams.png"));

% 4. 带外频谱泄漏诊断柱状图
if ~isempty(E.leakage_table) && height(E.leakage_table) > 0
    plot_leakage_bars(E.leakage_table, method_order, style, ...
        fullfile(curves_dir, "ta_sft_leakage_bars.png"));
end

% 5. 初相失配鲁棒性衰减曲线
if isfield(E, "robustness_table") && ~isempty(E.robustness_table) && height(E.robustness_table) > 0
    plot_phase_robustness(E.robustness_table, method_order, method_styles, style, "SSIM", ...
        fullfile(curves_dir, "ta_sft_phase_robustness_ssim.png"));
    plot_phase_robustness(E.robustness_table, method_order, method_styles, style, "PSNR", ...
        fullfile(curves_dir, "ta_sft_phase_robustness_psnr.png"));
end

% 6. 典型样本（测试集 B0_par_raw SSIM 中位数样本）SAR 图像与误差对比
render_metrics = render_ta_sft_montages(cfg, S60, record, E, method_order, style, images_dir);
if ~isempty(render_metrics) && height(render_metrics) > 0
    writetable(render_metrics, fullfile(tables_dir, "ta_sft_render_metrics.csv"));
end

fprintf("TA-SFT 专题图表已输出至：%s\n", fig_root);

function st = build_ta_sft_styles(method_order)
n_m = numel(method_order);
st = containers.Map("KeyType", "char", "ValueType", "any");
palette = [ ...
    0.00, 0.00, 0.00; ... % B0: 黑虚线
    0.45, 0.45, 0.45; ... % R1: 灰点划线
    0.85, 0.33, 0.10; ... % TA1_gate_zc: 橙红
    0.00, 0.45, 0.74; ... % TA2_anchor_zc: 深蓝
    0.30, 0.75, 0.93; ... % TA2_anchor_q5: 青蓝
    0.47, 0.67, 0.19; ... % TA3_mmse_fixed: 绿
    0.49, 0.18, 0.56];    % TA3_mmse_est: 紫
line_styles = ["--", "-.", "-", "-", "--", "-", "-"];
markers = ["none", "s", "o", "^", "v", "d", "p"];

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
title(sprintf("TA-SFT Test %s | %s (filled circle = sig > B0, p<0.05)", metric_name, ttype), ...
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
        title(ax, sprintf("TA-SFT Gain vs B0 | %s", ttype), "Interpreter", "none", "FontSize", style.title_size);
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
fig = figure("Visible", "off", "Color", "w", "Position", [100, 100, 1060, 740]);
tl = tiledlayout(fig, 2, 2, "TileSpacing", "compact", "Padding", "compact");

configs = { ...
    "TA1_gate_zc", "HyperparamValue", "TA1 Gate Threshold \tau"; ...
    "TA2_anchor_zc", "HyperparamValue", "TA2 Continuous Scale \alpha"; ...
    "TA2_anchor_q5", "HyperparamValue", "TA2 5-Level Scale \alpha"; ...
    "TA3_mmse_fixed", "HyperparamValue", "TA3 \gamma (Locked vs Blind \hat{\gamma})"};

for idx = 1:4
    ax = nexttile(tl);
    hold(ax, "on");
    grid(ax, "on");
    box(ax, "on");
    mkey = string(configs{idx, 1});
    col_name = string(configs{idx, 2});
    ttl = string(configs{idx, 3});

    for ttype = ["RSFT", "ASFT"]
        sub = sortrows(sum_tbl(sum_tbl.MethodKey == mkey & sum_tbl.ThresholdType == ttype, :), "Q");
        if isempty(sub)
            continue;
        end
        plot(ax, sub.Q, sub.(col_name), "-o", "LineWidth", 1.8, "DisplayName", char(ttype));
        if mkey == "TA3_mmse_fixed"
            est_sub = sortrows(sum_tbl(sum_tbl.MethodKey == "TA3_mmse_est" & sum_tbl.ThresholdType == ttype, :), "Q");
            if ~isempty(est_sub)
                plot(ax, est_sub.Q, est_sub.EstGamma_Mean, "--s", "LineWidth", 1.6, ...
                    "DisplayName", char(ttype + " (blind est)"));
            end
        end
    end
    title(ax, ttl, "FontSize", style.title_size);
    xlabel(ax, "Q", "FontSize", style.label_size);
    xticks(ax, unique(sort(sum_tbl.Q)));
    legend(ax, "Location", "best", "Interpreter", "none");
end

exportgraphics(fig, char(out_path), "Resolution", style.export_resolution);
close(fig);
end

function plot_leakage_bars(leak_tbl, method_order, style, out_path)
units_u = unique(leak_tbl(:, ["UnitKey", "GroupName", "ThresholdType"]), "rows", "stable");
if isempty(units_u)
    return;
end
n_u = height(units_u);
n_m = numel(method_order);

metrics = ["OffSupportRatio", "RangeOffRatio", "AzimuthOffRatio"];
titles = ["Total Off-Support Ratio", "Range Off-Support Ratio", "Azimuth Off-Support Ratio"];

fig = figure("Visible", "off", "Color", "w", "Position", [100, 100, 1380, 480]);
tl = tiledlayout(fig, 1, 3, "TileSpacing", "compact", "Padding", "compact");
x_labels = units_u.GroupName + " " + units_u.ThresholdType;

for p_idx = 1:3
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
    xtickangle(ax, 30);
    ylim(ax, [0, 1]);
    title(ax, titles(p_idx), "FontSize", style.title_size);
    if p_idx == 3
        legend(ax, cellstr(method_order), "Interpreter", "none", "Location", "eastoutside");
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
ta_methods = method_order(~ismember(method_order, ["B0_par_raw", "R1_par_raw"]));

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

    for k = 1:numel(ta_methods)
        mkey = ta_methods(k);
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

function render_tbl = render_ta_sft_montages(cfg, S60, record, E, method_order, style, images_dir)
render_tbl = table();
if ~isfield(cfg, "data_root") || strlength(string(cfg.data_root)) == 0 || ~isfolder(cfg.data_root)
    warning("zc1bit:MissingDataRootForRender", "未找到原始回波目录，跳过 TA-SFT SAR 成像拼图渲染。");
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

signals = zc1bit_stage2_load_signals(sample_row, S60, cfg.raw_decimation_step, "TA-SFT 渲染样本");
sig60 = signals{1};
img_gt = generate_gt_image(sig60, S60);

target_units = ["rb_R9A1_range_rsft", "ab_R1A9_azimuth_asft"];
avail_units = target_units(ismember(target_units, E.units.UnitKey));
if isempty(avail_units)
    avail_units = E.units.UnitKey(1);
end

rows = cell(0, 1);
use_fast = true;
if isfield(cfg.ta_sft, "use_fast_path")
    use_fast = logical(cfg.ta_sft.use_fast_path);
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
        hp_name = string(lr.HyperparamName);
        if hp_name ~= "none" && hp_name ~= "blind"
            ropts.(char(hp_name)) = double(lr.HyperparamValue);
        end
        ch_rec = apply_ta_sft_reconstruct(channel_1bit, carrier, mkey, ropts);
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
        fullfile(images_dir, "ta_sft_" + uk + "_montage.png"), style, "image");
    zc1bit_save_image_montage(img_list, title_list, ...
        fullfile(images_dir, "ta_sft_" + uk + "_error_montage.png"), style, "error");
end

if ~isempty(rows)
    render_tbl = vertcat(rows{:});
end
end
