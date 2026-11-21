function T = zc1bit_attach_pvalues(T, W, p_col, comparison_type, track)
%ZC1BIT_ATTACH_PVALUES 将配对 Wilcoxon 的 p 值按 UnitKey+VariantKey 附着到曲线表。
% W 为 paired_wilcoxon_vs_b0.csv；p_col 例如 P_SSIM / P_PSNR。

T = T;
p_col = string(p_col);
p_vals = nan(height(T), 1);
if isempty(W) || height(W) == 0
    T.(p_col) = p_vals;
    return;
end
vars = string(W.Properties.VariableNames);
if ~ismember("UnitKey", vars) || ~ismember("VariantKey", vars) || ~ismember(p_col, vars)
    T.(p_col) = p_vals;
    return;
end
mask = true(height(W), 1);
if ~isempty(comparison_type) && ismember("ComparisonType", vars)
    mask = mask & string(W.ComparisonType) == string(comparison_type);
end
if ~isempty(track) && ismember("Track", vars)
    mask = mask & string(W.Track) == string(track);
end
W = W(mask, :);
w_vals = W.(p_col);
for i = 1:height(T)
    hit = W.UnitKey == T.UnitKey(i) & W.VariantKey == T.VariantKey(i);
    if any(hit)
        p_vals(i) = w_vals(find(hit, 1));
    end
end
T.(p_col) = p_vals;
end
