function [fr_limit, fa_limit] = check_sft_nyquist(grid, params, S60, nyquist_margin)
%CHECK_SFT_NYQUIST 校验 SFT 归一化频率是否满足采样网格 Nyquist 安全边界。

if nargin < 4 || isempty(nyquist_margin)
    nyquist_margin = 0.98;
end
if ~isscalar(nyquist_margin) || ~isfinite(nyquist_margin) || nyquist_margin <= 0 || nyquist_margin >= 1
    error("zc1bit:NyquistMargin", "nyquist_margin 必须位于 (0, 1) 区间内。");
end

ba = resolve_azimuth_bandwidth(S60);
fr_limit = nyquist_margin * grid.Fs_up / (2 * S60.B);
fa_limit = nyquist_margin * grid.PRF_up / (2 * ba);

if params.fr_over_Br < 0 || params.fa_over_Ba < 0
    error("zc1bit:NegativeFrequency", "SFT 归一化频率 fr/Br 与 fa/Ba 不能为负。");
end
if params.fr_over_Br >= fr_limit || params.fa_over_Ba >= fa_limit
    error("zc1bit:SFTNyquistViolation", ...
        "SFT 频率超出 Nyquist 安全边界：要求 fr/Br < %.12g 且 fa/Ba < %.12g；当前为 %.12g、%.12g。", ...
        fr_limit, fa_limit, params.fr_over_Br, params.fa_over_Ba);
end
end
