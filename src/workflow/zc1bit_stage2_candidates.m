function [cand_table, meta] = zc1bit_stage2_candidates(unit_row, cfg, S60, mode, centers_tbl, existing_tags)
%ZC1BIT_STAGE2_CANDIDATES 生成单向四阈值（ZT/RT/RSFT/ASFT）的搜索网格与局部细化点集。
% 所有 SFT 频率点均经过 check_sft_nyquist 严格校验与动态安全裁剪。

if nargin < 4 || strlength(string(mode)) == 0
    mode = "coarse";
else
    mode = string(mode);
end
if nargin < 5
    centers_tbl = table();
end
if nargin < 6 || isempty(existing_tags)
    existing_tags = string.empty(0, 1);
else
    existing_tags = string(existing_tags(:));
end

s2 = cfg.stage2;
ttype = string(unit_row.ThresholdType);
rq = double(unit_row.Range_q);
aq = double(unit_row.Azimuth_q);
sgrid = resolve_sampling_grid(rq, aq, [S60.nrn, S60.nan], S60);
ba = resolve_azimuth_bandwidth(S60);
fr_limit = cfg.nyquist_margin * sgrid.Fs_up / (2 * S60.B);
fa_limit = cfg.nyquist_margin * sgrid.PRF_up / (2 * ba);

as_vals = sort(unique(round(double(s2.rt_as_grid(:)), 6)));
str_vals = sort(unique(round(double(s2.sft_str_grid(:)), 6)));

if ttype == "RSFT"
    freq_limit = fr_limit;
elseif ttype == "ASFT"
    freq_limit = fa_limit;
else
    freq_limit = NaN;
end

if isnan(freq_limit)
    legal_freqs = double.empty(0, 1);
else
    raw_f = sort(unique(round(double(s2.sft_freq_grid(:)), 6)));
    legal_freqs = raw_f(raw_f > 0 & raw_f < freq_limit * (1 - 1e-9));
    if isempty(legal_freqs)
        legal_freqs = unique(round([0.25, 0.5, 0.75] * freq_limit, 4));
        legal_freqs = legal_freqs(legal_freqs > 0 & legal_freqs < freq_limit * (1 - 1e-9));
    end
end

meta = struct( ...
    "as_grid", as_vals, ...
    "str_grid", str_vals, ...
    "freq_grid", legal_freqs, ...
    "fr_limit", fr_limit, ...
    "fa_limit", fa_limit, ...
    "freq_limit", freq_limit);

switch ttype
    case "ZT"
        if mode == "refine"
            cand_table = empty_cand_table();
            return;
        end
        cand_table = build_rows_from_vectors(unit_row, "zt", NaN, NaN, NaN, NaN, s2.phi0, "coarse");

    case "RT"
        if mode == "refine"
            cand_table = empty_cand_table();
            return;
        end
        n = numel(as_vals);
        ptags = "as" + arrayfun(@num_to_key, as_vals);
        cand_table = build_rows_from_vectors( ...
            unit_row, ptags, as_vals, nan(n, 1), nan(n, 1), nan(n, 1), repmat(s2.phi0, n, 1), "coarse");

    case {"RSFT", "ASFT"}
        if mode == "coarse"
            [str_mesh, freq_mesh] = ndgrid(str_vals, legal_freqs);
            cand_table = build_sft_rows(unit_row, sgrid, S60, cfg.nyquist_margin, ...
                str_mesh(:), freq_mesh(:), s2.phi0, "coarse");
        elseif mode == "refine"
            if ~isfield(s2, "sft_refine") || ~logical(s2.sft_refine) || isempty(centers_tbl)
                cand_table = empty_cand_table();
                return;
            end
            cand_table = build_sft_refine_rows(unit_row, sgrid, S60, cfg, freq_limit, centers_tbl, existing_tags);
        else
            error("zc1bit:UnknownCandidateMode", "未知候选生成模式：%s", mode);
        end

    otherwise
        error("zc1bit:Stage2ThresholdNotAllowed", "阶段二不支持阈值类型：%s", ttype);
end

if ~isempty(cand_table) && ~isempty(existing_tags)
    cand_table = cand_table(~ismember(cand_table.ParamTag, existing_tags), :);
    cand_table.ParamIdx = (1:height(cand_table)).';
end
end

function tbl = build_sft_refine_rows(unit_row, sgrid, S60, cfg, freq_limit, centers_tbl, existing_tags)
s2 = cfg.stage2;
ttype = string(unit_row.ThresholdType);
d_str = double(s2.sft_refine_str_step);
d_freq = double(s2.sft_refine_freq_step);

pairs = zeros(0, 2);
for k = 1:height(centers_tbl)
    c_str = double(centers_tbl.STR_dB(k));
    if ttype == "RSFT"
        c_freq = double(centers_tbl.FrOverBr(k));
    else
        c_freq = double(centers_tbl.FaOverBa(k));
    end
    s_list = unique(round(c_str + [-d_str, 0, d_str], 6));
    f_list = unique(round(c_freq + [-d_freq, 0, d_freq], 6));
    f_list = f_list(f_list > 0 & f_list < freq_limit * (1 - 1e-9));
    if isempty(f_list)
        continue;
    end
    [sm, fm] = ndgrid(s_list, f_list);
    pairs = [pairs; [sm(:), fm(:)]]; %#ok<AGROW>
end

if isempty(pairs)
    tbl = empty_cand_table();
    return;
end
pairs = unique(pairs, "rows", "stable");
tbl = build_sft_rows(unit_row, sgrid, S60, cfg.nyquist_margin, pairs(:, 1), pairs(:, 2), s2.phi0, "refine");
if ~isempty(existing_tags) && ~isempty(tbl)
    tbl = tbl(~ismember(tbl.ParamTag, existing_tags), :);
    tbl.ParamIdx = (1:height(tbl)).';
end
end

function tbl = build_sft_rows(unit_row, sgrid, S60, nyquist_margin, str_vec, freq_vec, phi0, phase_tag)
n = numel(str_vec);
if n == 0
    tbl = empty_cand_table();
    return;
end
if string(unit_row.ThresholdType) == "RSFT"
    fr_vec = freq_vec(:);
    fa_vec = zeros(n, 1);
else
    fr_vec = zeros(n, 1);
    fa_vec = freq_vec(:);
end

ptags = strings(n, 1);
for k = 1:n
    check_sft_nyquist(sgrid, struct("fr_over_Br", fr_vec(k), "fa_over_Ba", fa_vec(k)), S60, nyquist_margin);
    ptags(k) = sprintf("s%s_fr%s_fa%s", num_to_key(str_vec(k)), num_to_key(fr_vec(k)), num_to_key(fa_vec(k)));
end
tbl = build_rows_from_vectors(unit_row, ptags, nan(n, 1), str_vec(:), fr_vec, fa_vec, repmat(phi0, n, 1), phase_tag);
end

function tbl = build_rows_from_vectors(unit_row, ptags, as_vec, str_vec, fr_vec, fa_vec, phi0_vec, phase_tag)
n = numel(ptags);
ParamIdx = (1:n).';
ParamTag = string(ptags(:));
ParamKey = repmat(string(unit_row.UnitKey), n, 1) + "_" + ParamTag;
As = double(as_vec(:));
STR_dB = double(str_vec(:));
FrOverBr = double(fr_vec(:));
FaOverBa = double(fa_vec(:));
Phi0 = double(phi0_vec(:));
PhaseTag = repmat(string(phase_tag), n, 1);
tbl = table(ParamIdx, ParamTag, ParamKey, As, STR_dB, FrOverBr, FaOverBa, Phi0, PhaseTag);
end

function tbl = empty_cand_table()
tbl = table(zeros(0, 1), strings(0, 1), strings(0, 1), zeros(0, 1), zeros(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), strings(0, 1), ...
    'VariableNames', {'ParamIdx', 'ParamTag', 'ParamKey', 'As', 'STR_dB', 'FrOverBr', 'FaOverBa', 'Phi0', 'PhaseTag'});
end

function text = num_to_key(value)
text = string(sprintf("%.6g", value));
text = replace(text, "-", "m");
text = replace(text, ".", "p");
end
