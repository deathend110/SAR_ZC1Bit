function selected_variants = zc1bit_resolve_stage2_variants(cfg, p)
%ZC1BIT_RESOLVE_STAGE2_VARIANTS 解析阶段二参评变体（固定基线 B0 + 优化口径自动读取 Top K 非 B0 变体）。
% 口径优先按 SSIMSigWinCount -> SSIMWinCount -> 过采样(Q>1) MeanSSIM/MeanPSNR -> 全局均值排序，
% 并默认合并同一规则下数值等价的 raw/dc 与 parallel/serial 冗余变体。

if nargin < 2 || isempty(p)
    p = zc1bit_paths(cfg);
end
s2 = cfg.stage2;
base_key = string(s2.baseline_variant);

if isfield(s2, "ranking_file") && strlength(string(s2.ranking_file)) > 0
    ranking_path = string(s2.ranking_file);
else
    ranking_path = fullfile(p.screen, "screen_rule_ranking.csv");
end
summary_path = fullfile(fileparts(ranking_path), "screen_summary.csv");

has_ranking = isfile(ranking_path);
if has_ranking
    ranking_tbl = load_and_score_ranking(ranking_path, summary_path);
else
    ranking_tbl = table();
end

manual_keys = string.empty(0, 1);
if isfield(s2, "selected_variants") && ~isempty(s2.selected_variants)
    manual_keys = string(s2.selected_variants(:));
    manual_keys = manual_keys(strlength(manual_keys) > 0 & manual_keys ~= base_key);
    manual_keys = unique(manual_keys, "stable");
end

if isempty(manual_keys)
    if ~has_ranking
        error("zc1bit:MissingScreenRanking", ...
            "未找到阶段一规则排名表：%s。请先运行 run_zc1bit_screen(cfg) 或设置 cfg.stage2.selected_variants。", ...
            ranking_path);
    end
    non_b0 = ranking_tbl(ranking_tbl.RuleKey ~= "B0", :);
    if isempty(non_b0)
        error("zc1bit:EmptyStage2Candidates", "阶段一排名表 %s 中不包含任何非 B0 候选规则。", ranking_path);
    end
    if ~isfield(s2, "dedup_equivalent") || logical(s2.dedup_equivalent)
        non_b0 = deduplicate_equivalent_variants(non_b0);
    end
    if isfield(s2, "family_diversified") && logical(s2.family_diversified)
        non_b0 = pick_family_representatives(non_b0);
    end
    if isfield(s2, "distinct_rules_only") && logical(s2.distinct_rules_only)
        non_b0 = pick_first_per_rule(non_b0);
    end
    top_k = double(s2.top_k);
    if ~isscalar(top_k) || ~isfinite(top_k) || top_k < 1 || mod(top_k, 1) ~= 0
        error("zc1bit:InvalidTopK", "cfg.stage2.top_k 必须为正整数。");
    end
    n_pick = min(top_k, height(non_b0));
    cand_tbl = non_b0(1:n_pick, :);
else
    cand_tbl = build_candidates_from_keys(manual_keys, ranking_tbl);
end

n_cand = height(cand_tbl);
if n_cand == 0
    error("zc1bit:EmptyStage2Candidates", "阶段二未选出任何伪 1.5-bit 候选变体。");
end

n_total = 1 + n_cand;
VariantOrder = (1:n_total).';
VariantRole = strings(n_total, 1);
VariantKey = strings(n_total, 1);
RuleKey = strings(n_total, 1);
DecisionMode = strings(n_total, 1);
RemoveDC = false(n_total, 1);
Stage1Rank = nan(n_total, 1);
Stage1MeanSSIM = nan(n_total, 1);
Stage1MeanPSNR = nan(n_total, 1);
Stage1UpsampledSSIM = nan(n_total, 1);
Stage1SSIMWinCount = nan(n_total, 1);
Stage1SSIMSigWinCount = nan(n_total, 1);

base_info = lookup_or_parse_variant(base_key, ranking_tbl);
VariantRole(1) = "baseline";
VariantKey(1) = base_info.VariantKey;
RuleKey(1) = base_info.RuleKey;
DecisionMode(1) = base_info.DecisionMode;
RemoveDC(1) = base_info.RemoveDC;
Stage1Rank(1) = base_info.Stage1Rank;
Stage1MeanSSIM(1) = base_info.Stage1MeanSSIM;
Stage1MeanPSNR(1) = base_info.Stage1MeanPSNR;
Stage1UpsampledSSIM(1) = base_info.Stage1UpsampledSSIM;
Stage1SSIMWinCount(1) = base_info.Stage1SSIMWinCount;
Stage1SSIMSigWinCount(1) = base_info.Stage1SSIMSigWinCount;

for idx = 1:n_cand
    row_idx = idx + 1;
    VariantRole(row_idx) = "top" + string(idx);
    VariantKey(row_idx) = cand_tbl.VariantKey(idx);
    RuleKey(row_idx) = cand_tbl.RuleKey(idx);
    DecisionMode(row_idx) = cand_tbl.DecisionMode(idx);
    RemoveDC(row_idx) = cand_tbl.RemoveDC(idx);
    Stage1Rank(row_idx) = cand_tbl.Rank(idx);
    Stage1MeanSSIM(row_idx) = cand_tbl.MeanSSIM(idx);
    Stage1MeanPSNR(row_idx) = cand_tbl.MeanPSNR(idx);
    Stage1UpsampledSSIM(row_idx) = cand_tbl.UpsampledMeanSSIM(idx);
    Stage1SSIMWinCount(row_idx) = cand_tbl.SSIMWinCount(idx);
    Stage1SSIMSigWinCount(row_idx) = cand_tbl.SSIMSigWinCount(idx);
end

selected_variants = table( ...
    VariantOrder, VariantRole, VariantKey, RuleKey, DecisionMode, RemoveDC, ...
    Stage1Rank, Stage1SSIMSigWinCount, Stage1SSIMWinCount, ...
    Stage1UpsampledSSIM, Stage1MeanSSIM, Stage1MeanPSNR);
end

function tbl = load_and_score_ranking(ranking_path, summary_path)
tbl = readtable(ranking_path, "TextType", "string");
required = ["Rank", "RuleKey", "DecisionMode", "RemoveDC", "VariantKey", "MeanSSIM", "MeanPSNR"];
vars = string(tbl.Properties.VariableNames);
if ~all(ismember(required, vars))
    error("zc1bit:InvalidScreenRanking", "阶段一排名表 %s 缺少必要列。", ranking_path);
end
tbl.RuleKey = string(tbl.RuleKey);
tbl.DecisionMode = string(tbl.DecisionMode);
tbl.VariantKey = string(tbl.VariantKey);
tbl.RemoveDC = to_logical_col(tbl.RemoveDC);
tbl.Rank = double(tbl.Rank);
tbl.MeanSSIM = double(tbl.MeanSSIM);
tbl.MeanPSNR = double(tbl.MeanPSNR);

if ~ismember("SSIMWinCount", vars)
    tbl.SSIMWinCount = zeros(height(tbl), 1);
else
    tbl.SSIMWinCount = double(tbl.SSIMWinCount);
end
if ~ismember("SSIMSigWinCount", vars)
    tbl.SSIMSigWinCount = zeros(height(tbl), 1);
else
    tbl.SSIMSigWinCount = double(tbl.SSIMSigWinCount);
end
if ~ismember("MeanAbsENLDiff", vars)
    tbl.MeanAbsENLDiff = nan(height(tbl), 1);
else
    tbl.MeanAbsENLDiff = double(tbl.MeanAbsENLDiff);
end

if ismember("UpsampledMeanSSIM", vars) && ismember("UpsampledMeanPSNR", vars)
    tbl.UpsampledMeanSSIM = double(tbl.UpsampledMeanSSIM);
    tbl.UpsampledMeanPSNR = double(tbl.UpsampledMeanPSNR);
elseif isfile(summary_path)
    [up_ssim, up_psnr] = compute_upsampled_metrics_from_summary(tbl.VariantKey, tbl.MeanSSIM, tbl.MeanPSNR, summary_path);
    tbl.UpsampledMeanSSIM = up_ssim;
    tbl.UpsampledMeanPSNR = up_psnr;
else
    tbl.UpsampledMeanSSIM = tbl.MeanSSIM;
    tbl.UpsampledMeanPSNR = tbl.MeanPSNR;
end

tbl = sortrows(tbl, ...
    {'SSIMSigWinCount', 'SSIMWinCount', 'UpsampledMeanSSIM', 'UpsampledMeanPSNR', ...
    'MeanSSIM', 'MeanAbsENLDiff', 'RemoveDC', 'DecisionMode', 'VariantKey'}, ...
    {'descend', 'descend', 'descend', 'descend', 'descend', 'ascend', 'ascend', 'ascend', 'ascend'});
tbl.Rank = (1:height(tbl)).';
end

function [up_ssim, up_psnr] = compute_upsampled_metrics_from_summary(vkeys, fallback_ssim, fallback_psnr, summary_path)
sum_tbl = readtable(summary_path, "TextType", "string");
n = numel(vkeys);
up_ssim = fallback_ssim;
up_psnr = fallback_psnr;
if ~all(ismember(["VariantKey", "Q", "SSIM_Mean", "PSNR_Mean"], string(sum_tbl.Properties.VariableNames)))
    return;
end
sum_tbl.VariantKey = string(sum_tbl.VariantKey);
for k = 1:n
    sub = sum_tbl(sum_tbl.VariantKey == vkeys(k) & double(sum_tbl.Q) > 1, :);
    if ~isempty(sub)
        up_ssim(k) = mean(double(sub.SSIM_Mean));
        up_psnr(k) = mean(double(sub.PSNR_Mean));
    end
end
end

function out_tbl = deduplicate_equivalent_variants(tbl)
% 合并数值等价冗余变体：
% 1) 同 RuleKey 下，仅当判定模式或去DC产生实质差异时才多占一个席位；
% 2) 同 (RuleKey, DecisionMode) 下，若 RemoveDC=true 未带来实质提升，保留严谨三值 raw。
n = height(tbl);
keep = false(n, 1);
kept_idx = zeros(0, 1);

for k = 1:n
    rk = tbl.RuleKey(k);
    dm = tbl.DecisionMode(k);
    dc = tbl.RemoveDC(k);
    is_dup = false;

    for j = kept_idx(:).'
        if tbl.RuleKey(j) ~= rk
            continue;
        end
        d_up = abs(tbl.UpsampledMeanSSIM(k) - tbl.UpsampledMeanSSIM(j));
        d_all = abs(tbl.MeanSSIM(k) - tbl.MeanSSIM(j));
        same_win = tbl.SSIMSigWinCount(k) == tbl.SSIMSigWinCount(j) && ...
            tbl.SSIMWinCount(k) == tbl.SSIMWinCount(j);

        if tbl.DecisionMode(j) == dm && dc && ~tbl.RemoveDC(j) && ...
                tbl.SSIMSigWinCount(k) <= tbl.SSIMSigWinCount(j) && d_up < 5e-4 && d_all < 5e-4
            is_dup = true;
            break;
        end
        if same_win && d_up < 1e-4 && d_all < 1e-4
            is_dup = true;
            break;
        end
    end

    if ~is_dup
        keep(k) = true;
        kept_idx(end + 1, 1) = k; %#ok<AGROW>
    end
end
out_tbl = tbl(keep, :);
end

function out_tbl = pick_family_representatives(tbl)
% 按物理过零机理族去重覆盖：
%   {R1} 保留非等价的 parallel / serial（各自最多 1 席）；
%   {R2,R3,R4,R5}、{R6,R7}、{R8,R9}、{L1} 每族仅保留阶段一排名最优的 1 名代表。
n = height(tbl);
keep = false(n, 1);
seen_slots = string.empty(0, 1);
for k = 1:n
    fam = rule_mechanism_family(tbl.RuleKey(k));
    if fam == "F_BOTH_BOTH"
        slot_key = fam + "_" + tbl.DecisionMode(k);
    else
        slot_key = fam;
    end
    if ~ismember(slot_key, seen_slots)
        seen_slots(end + 1, 1) = slot_key; %#ok<AGROW>
        keep(k) = true;
    end
end
out_tbl = tbl(keep, :);
end

function fam = rule_mechanism_family(rule_key)
switch string(rule_key)
    case "R1"
        fam = "F_BOTH_BOTH";
    case {"R2", "R3", "R4", "R5"}
        fam = "F_ONE_ONE";
    case {"R6", "R7"}
        fam = "F_BOTH_ONE";
    case {"R8", "R9"}
        fam = "F_ONE_BOTH";
    case "L1"
        fam = "F_L1";
    otherwise
        fam = "F_OTHER_" + string(rule_key);
end
end

function out_tbl = pick_first_per_rule(tbl)
n = height(tbl);
keep = false(n, 1);
seen = string.empty(0, 1);
for k = 1:n
    rk = tbl.RuleKey(k);
    if ~ismember(rk, seen)
        seen(end + 1, 1) = rk; %#ok<AGROW>
        keep(k) = true;
    end
end
out_tbl = tbl(keep, :);
end

function cand_tbl = build_candidates_from_keys(keys, ranking_tbl)
n = numel(keys);
Rank = nan(n, 1);
RuleKey = strings(n, 1);
DecisionMode = strings(n, 1);
RemoveDC = false(n, 1);
VariantKey = strings(n, 1);
MeanSSIM = nan(n, 1);
MeanPSNR = nan(n, 1);
UpsampledMeanSSIM = nan(n, 1);
SSIMWinCount = nan(n, 1);
SSIMSigWinCount = nan(n, 1);
for k = 1:n
    info = lookup_or_parse_variant(keys(k), ranking_tbl);
    Rank(k) = info.Stage1Rank;
    RuleKey(k) = info.RuleKey;
    DecisionMode(k) = info.DecisionMode;
    RemoveDC(k) = info.RemoveDC;
    VariantKey(k) = info.VariantKey;
    MeanSSIM(k) = info.Stage1MeanSSIM;
    MeanPSNR(k) = info.Stage1MeanPSNR;
    UpsampledMeanSSIM(k) = info.Stage1UpsampledSSIM;
    SSIMWinCount(k) = info.Stage1SSIMWinCount;
    SSIMSigWinCount(k) = info.Stage1SSIMSigWinCount;
end
cand_tbl = table(Rank, RuleKey, DecisionMode, RemoveDC, VariantKey, ...
    MeanSSIM, MeanPSNR, UpsampledMeanSSIM, SSIMWinCount, SSIMSigWinCount);
end

function info = lookup_or_parse_variant(vkey, ranking_tbl)
vkey = string(vkey);
if ~isempty(ranking_tbl) && any(ranking_tbl.VariantKey == vkey)
    idx = find(ranking_tbl.VariantKey == vkey, 1);
    info = struct( ...
        "VariantKey", vkey, ...
        "RuleKey", ranking_tbl.RuleKey(idx), ...
        "DecisionMode", ranking_tbl.DecisionMode(idx), ...
        "RemoveDC", logical(ranking_tbl.RemoveDC(idx)), ...
        "Stage1Rank", double(ranking_tbl.Rank(idx)), ...
        "Stage1MeanSSIM", double(ranking_tbl.MeanSSIM(idx)), ...
        "Stage1MeanPSNR", double(ranking_tbl.MeanPSNR(idx)), ...
        "Stage1UpsampledSSIM", double(ranking_tbl.UpsampledMeanSSIM(idx)), ...
        "Stage1SSIMWinCount", double(ranking_tbl.SSIMWinCount(idx)), ...
        "Stage1SSIMSigWinCount", double(ranking_tbl.SSIMSigWinCount(idx)));
    return;
end
[rk, dm, dc] = parse_variant_key(vkey);
info = struct( ...
    "VariantKey", vkey, ...
    "RuleKey", rk, ...
    "DecisionMode", dm, ...
    "RemoveDC", dc, ...
    "Stage1Rank", NaN, ...
    "Stage1MeanSSIM", NaN, ...
    "Stage1MeanPSNR", NaN, ...
    "Stage1UpsampledSSIM", NaN, ...
    "Stage1SSIMWinCount", NaN, ...
    "Stage1SSIMSigWinCount", NaN);
end

function [rule_key, decision_mode, remove_dc] = parse_variant_key(vkey)
parts = split(string(vkey), "_");
if numel(parts) ~= 3
    error("zc1bit:InvalidVariantKey", "无法解析 VariantKey：%s（期望格式 <Rule>_<par|ser>_<raw|dc>）。", vkey);
end
rule_key = parts(1);
if parts(2) == "par"
    decision_mode = "parallel";
elseif parts(2) == "ser"
    decision_mode = "serial";
else
    error("zc1bit:InvalidVariantKey", "VariantKey 判定模式非法：%s", vkey);
end
if parts(3) == "raw"
    remove_dc = false;
elseif parts(3) == "dc"
    remove_dc = true;
else
    error("zc1bit:InvalidVariantKey", "VariantKey 去 DC 标识非法：%s", vkey);
end
end

function vals = to_logical_col(raw)
if islogical(raw)
    vals = logical(raw(:));
elseif isnumeric(raw)
    vals = logical(raw(:) ~= 0);
else
    s = lower(strtrim(string(raw(:))));
    vals = (s == "true" | s == "1");
end
end
