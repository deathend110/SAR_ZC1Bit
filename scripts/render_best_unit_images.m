%% 各自最优单元 SAR 图像：B0 + 6 个过零变体
% 每个方法在 evaluate_summary 中按测试集 SSIM 最大、PSNR 破局选择自己的最优单元，
% 固定同一样本渲染。注意：不同方法的最优单元可能不同，这不是受控同单元对比；
% 受控对比请运行 render_common_unit_images.m。
% 输出目录：results_zc1bit/figures/images/best_unit/

%% 可调参数
sample_id = [];   % [] = 自动取测试集 B0 SSIM 中位数样本；也可填入 SampleID

script_dir = string(fileparts(mfilename("fullpath")));
addpath(char(script_dir));
addpath(char(fileparts(script_dir)));

[cfg, p, S60, record] = zc1bit_script_begin(); %#ok<ASGLU>
D = zc1bit_load_stage2_results(cfg, p, record);
style = zc1bit_plot_style();
[sample_row, signal60, img_gt] = zc1bit_select_render_sample( ...
    cfg, D, record, S60, sample_id);
out_dir = fullfile(p.figures, "images", "best_unit");
zc1bit_ensure_dir(out_dir);

variant_keys = [D.baseline_key; D.pseudo_keys];
n_variants = numel(variant_keys);
imgs = cell(n_variants, 1);
titles = strings(n_variants, 1);
VariantKey = strings(n_variants, 1);
VariantRole = strings(n_variants, 1);
UnitKey = strings(n_variants, 1);
Branch = strings(n_variants, 1);
ThresholdType = strings(n_variants, 1);
Q = zeros(n_variants, 1);
PSNR = nan(n_variants, 1);
SSIM = nan(n_variants, 1);
ExpectedPSNR = nan(n_variants, 1);
ExpectedSSIM = nan(n_variants, 1);
MatchExpected = false(n_variants, 1);
ImageFile = strings(n_variants, 1);

for v_idx = 1:n_variants
    vkey = variant_keys(v_idx);
    sub = D.eval_ind(D.eval_ind.VariantKey == vkey, :);
    best = zc1bit_select_best_unit(sub, vkey);
    ukey = best.UnitKey;
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

    fast_ctx = build_fast_focus_context(S60, ...
        double(plan.Range_q), double(plan.Azimuth_q));
    [img, ~, m] = zc1bit_render_locked_image( ...
        cfg, S60, sample_row, signal60, plan, fast_ctx, img_gt, expected);
    imgs{v_idx} = img;

    role = D.variants.VariantRole(D.variants.VariantKey == vkey);
    VariantKey(v_idx) = vkey;
    VariantRole(v_idx) = role(1);
    UnitKey(v_idx) = ukey;
    Branch(v_idx) = best.Branch;
    ThresholdType(v_idx) = best.ThresholdType;
    Q(v_idx) = best.Q;
    PSNR(v_idx) = m.PSNR;
    SSIM(v_idx) = m.SSIM;
    ExpectedPSNR(v_idx) = m.ExpectedPSNR;
    ExpectedSSIM(v_idx) = m.ExpectedSSIM;
    MatchExpected(v_idx) = logical(m.MatchExpected);
    ImageFile(v_idx) = fullfile(out_dir, "best_unit_" + vkey + ".png");
    imwrite(uint8(round(255 * min(max(img, 0), 1))), char(ImageFile(v_idx)));

    titles(v_idx) = sprintf("%s\n%s | %s | Q=%d\nPSNR %.2f | SSIM %.4f", ...
        vkey, ukey, string(best.ThresholdType), double(best.Q), m.PSNR, m.SSIM);
end

all_imgs = [{img_gt}; imgs];
all_titles = ["GT (normalized amplitude)"; titles];
zc1bit_save_image_montage(all_imgs, all_titles, ...
    fullfile(out_dir, "best_unit_montage.png"), style, "image");
zc1bit_save_image_montage(all_imgs, all_titles, ...
    fullfile(out_dir, "best_unit_error_montage.png"), style, "error");

metrics = table(VariantKey, VariantRole, UnitKey, Branch, ThresholdType, Q, ...
    PSNR, SSIM, ExpectedPSNR, ExpectedSSIM, MatchExpected, ImageFile);
writetable(metrics, fullfile(p.figures, "tables", "render_metrics_best_unit.csv"));
fprintf("已输出各自最优单元图像：%s\n", out_dir);
fprintf("样本 SampleID=%d；指标复核见 %s\n", sample_row.SampleID, ...
    fullfile(p.figures, "tables", "render_metrics_best_unit.csv"));
