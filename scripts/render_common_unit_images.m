%% 共同单元 SAR 图像：B0 + 6 个过零变体
% 默认每个 branch 自动选择 B0 测试 SSIM 最高的单元（common_unit_mode="b0_best"），
% 也可在顶部用 unit_keys 指定。同一单元、同一样本下 7 个方法各用该单元的锁参，
% 是严格可比的视觉对比。
% 输出目录：results_zc1bit/figures/images/common_unit/

%% 可调参数
unit_keys = string.empty(0, 1);   % 例：["rb_R9A1_range_rsft"; "ab_R1A9_azimuth_asft"]
sample_id = [];                    % [] = 自动取测试集 B0 SSIM 中位数样本

script_dir = string(fileparts(mfilename("fullpath")));
addpath(char(script_dir));
addpath(char(fileparts(script_dir)));

[cfg, p, S60, record] = zc1bit_script_begin(); %#ok<ASGLU>
D = zc1bit_load_stage2_results(cfg, p, record);
style = zc1bit_plot_style();
[sample_row, signal60, img_gt] = zc1bit_select_render_sample( ...
    cfg, D, record, S60, sample_id);
out_dir = fullfile(p.figures, "images", "common_unit");
zc1bit_ensure_dir(out_dir);

selection_mode = "user";
unit_keys = string(unit_keys(:));
if isempty(unit_keys)
    selection_mode = "b0_best";
    branches = unique(D.eval_ind.Branch, "stable");
    for b_idx = 1:numel(branches)
        br = branches(b_idx);
        b0 = D.eval_ind(D.eval_ind.VariantKey == D.baseline_key & ...
            D.eval_ind.Branch == br, :);
        b0 = sortrows(b0, {'SSIM_Mean', 'PSNR_Mean'}, {'descend', 'descend'});
        if isempty(b0)
            error("zc1bit:PlotMissingBaseline", "branch %s 缺少 B0 评估行。", br);
        end
        unit_keys(end + 1, 1) = b0.UnitKey(1); %#ok<AGROW>
    end
end
unit_keys = unique(unit_keys, "stable");
if isempty(unit_keys)
    error("zc1bit:PlotEmptyCommonUnits", "没有可用的共同单元。");
end

variant_keys = [D.baseline_key; D.pseudo_keys];
n_variants = numel(variant_keys);
metric_parts = cell(numel(unit_keys), 1);
unit_rows = cell(numel(unit_keys), 1);

for u_idx = 1:numel(unit_keys)
    ukey = unit_keys(u_idx);
    unit_row = D.units(D.units.UnitKey == ukey, :);
    if height(unit_row) ~= 1
        error("zc1bit:PlotUnitKey", "单元 %s 在 units 中应有 1 行，实际 %d 行。", ...
            ukey, height(unit_row));
    end
    unit_rows{u_idx} = unit_row;
    fast_ctx = build_fast_focus_context(S60, ...
        double(unit_row.Range_q), double(unit_row.Azimuth_q));

    imgs = cell(n_variants, 1);
    titles = strings(n_variants, 1);
    VariantKey = strings(n_variants, 1);
    VariantRole = strings(n_variants, 1);
    UnitKey = repmat(ukey, n_variants, 1);
    Branch = repmat(unit_row.Branch, n_variants, 1);
    ThresholdType = repmat(unit_row.ThresholdType, n_variants, 1);
    Q = repmat(double(unit_row.Q), n_variants, 1);
    PSNR = nan(n_variants, 1);
    SSIM = nan(n_variants, 1);
    ExpectedPSNR = nan(n_variants, 1);
    ExpectedSSIM = nan(n_variants, 1);
    MatchExpected = false(n_variants, 1);
    ImageFile = strings(n_variants, 1);

    for v_idx = 1:n_variants
        vkey = variant_keys(v_idx);
        plan = D.locked_plan(D.locked_plan.Track == "independent_recalibrated" & ...
            D.locked_plan.UnitKey == ukey & D.locked_plan.VariantKey == vkey, :);
        if height(plan) ~= 1
            error("zc1bit:PlotPlanRow", ...
                "锁定计划中 %s / %s 应有 1 行，实际 %d 行。", ukey, vkey, height(plan));
        end
        expected = D.eval_detail_ind( ...
            D.eval_detail_ind.SampleID == sample_row.SampleID & ...
            D.eval_detail_ind.UnitKey == ukey & ...
            D.eval_detail_ind.VariantKey == vkey, :);
        if height(expected) ~= 1
            error("zc1bit:PlotExpectedRow", ...
                "evaluate_detail 中 SampleID=%d / %s / %s 应有 1 行，实际 %d 行。", ...
                sample_row.SampleID, ukey, vkey, height(expected));
        end

        [img, ~, m] = zc1bit_render_locked_image( ...
            cfg, S60, sample_row, signal60, plan, fast_ctx, img_gt, expected);
        imgs{v_idx} = img;
        role = D.variants.VariantRole(D.variants.VariantKey == vkey);
        VariantKey(v_idx) = vkey;
        VariantRole(v_idx) = role(1);
        PSNR(v_idx) = m.PSNR;
        SSIM(v_idx) = m.SSIM;
        ExpectedPSNR(v_idx) = m.ExpectedPSNR;
        ExpectedSSIM(v_idx) = m.ExpectedSSIM;
        MatchExpected(v_idx) = logical(m.MatchExpected);
        ImageFile(v_idx) = fullfile(out_dir, ...
            "common_unit_" + ukey + "_" + vkey + ".png");
        imwrite(uint8(round(255 * min(max(img, 0), 1))), char(ImageFile(v_idx)));
        titles(v_idx) = sprintf("%s\n%s | %s | Q=%d\nPSNR %.2f | SSIM %.4f", ...
            vkey, ukey, string(unit_row.ThresholdType), double(unit_row.Q), ...
            m.PSNR, m.SSIM);
    end

    all_imgs = [{img_gt}; imgs];
    all_titles = ["GT (normalized amplitude)"; titles];
    zc1bit_save_image_montage(all_imgs, all_titles, ...
        fullfile(out_dir, "common_unit_" + ukey + "_montage.png"), style, "image");
    zc1bit_save_image_montage(all_imgs, all_titles, ...
        fullfile(out_dir, "common_unit_" + ukey + "_error_montage.png"), style, "error");

    metric_parts{u_idx} = table(VariantKey, VariantRole, UnitKey, Branch, ...
        ThresholdType, Q, PSNR, SSIM, ExpectedPSNR, ExpectedSSIM, ...
        MatchExpected, ImageFile);
end

metrics = vertcat(metric_parts{:});
writetable(metrics, fullfile(p.figures, "tables", "render_metrics_common_unit.csv"));

unit_table = vertcat(unit_rows{:});
unit_table.SelectionMode = repmat(string(selection_mode), height(unit_table), 1);
writetable(unit_table(:, ["UnitKey", "Branch", "ThresholdType", "Q", "SelectionMode"]), ...
    fullfile(p.figures, "tables", "common_units.csv"));
fprintf("已输出共同单元图像：%s\n", out_dir);
fprintf("样本 SampleID=%d；共同单元见 %s\n", sample_row.SampleID, ...
    fullfile(p.figures, "tables", "common_units.csv"));
