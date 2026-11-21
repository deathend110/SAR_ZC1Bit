%% 重校准消融：独立重校准 − B0 参数直迁（测试集）
% 仅覆盖 RT/RSFT/ASFT（ZT 无自由参数，不生成 transfer 行）。
% 输出目录：results_zc1bit/figures/curves/ablation/transfer_vs_recalibration/

script_dir = string(fileparts(mfilename("fullpath")));
addpath(char(script_dir));
addpath(char(fileparts(script_dir)));

[cfg, p, S60, record] = zc1bit_script_begin(); %#ok<ASGLU>
D = zc1bit_load_stage2_results(cfg, p, record);
style = zc1bit_plot_style();
variants = D.pseudo_keys;

if ~any(D.wilcoxon.ComparisonType == "recal_vs_transfer")
    fprintf("paired_wilcoxon_vs_b0.csv 中没有 recal_vs_transfer 行，跳过消融图。\n");
    return;
end
T = D.wilcoxon(D.wilcoxon.ComparisonType == "recal_vs_transfer", :);
out_dir = fullfile(p.figures, "curves", "ablation", "transfer_vs_recalibration");
zc1bit_plot_delta_vs_q(T, "DeltaSSIM_Mean", "P_SSIM", ...
    "Recal - Transfer ΔSSIM", out_dir, variants, style, ...
    struct("metric_label", "SSIM", "title_prefix", "Recal - Transfer ΔSSIM"));
zc1bit_plot_delta_vs_q(T, "DeltaPSNR_Mean", "P_PSNR", ...
    "Recal - Transfer ΔPSNR", out_dir, variants, style, ...
    struct("metric_label", "PSNR (dB)", "title_prefix", "Recal - Transfer ΔPSNR"));
fprintf("已输出重校准消融曲线：%s\n", out_dir);
