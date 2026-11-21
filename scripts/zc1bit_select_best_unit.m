function row = zc1bit_select_best_unit(T, variant_key)
%ZC1BIT_SELECT_BEST_UNIT 在评估汇总表中按 SSIM 最大、PSNR 破局选择单元。

variant_key = string(variant_key);
sub = T(T.VariantKey == variant_key, :);
if isempty(sub)
    error("zc1bit:PlotNoVariantRows", "评估汇总中没有变体：%s", variant_key);
end
sub = sortrows(sub, {'SSIM_Mean', 'PSNR_Mean'}, {'descend', 'descend'});
row = sub(1, :);
end
