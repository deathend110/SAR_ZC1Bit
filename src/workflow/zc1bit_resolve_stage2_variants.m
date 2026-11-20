function selected_variants = zc1bit_resolve_stage2_variants(cfg, p)
%ZC1BIT_RESOLVE_STAGE2_VARIANTS 解析阶段二参评变体（固定基线 B0 + 自动读取阶段一 Top K 非 B0 变体）。

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

has_ranking = isfile(ranking_path);
if has_ranking
    ranking_tbl = load_ranking_table(ranking_path);
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

base_info = lookup_or_parse_variant(base_key, ranking_tbl);
VariantRole(1) = "baseline";
VariantKey(1) = base_info.VariantKey;
RuleKey(1) = base_info.RuleKey;
DecisionMode(1) = base_info.DecisionMode;
RemoveDC(1) = base_info.RemoveDC;
Stage1Rank(1) = base_info.Stage1Rank;
Stage1MeanSSIM(1) = base_info.Stage1MeanSSIM;
Stage1MeanPSNR(1) = base_info.Stage1MeanPSNR;

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
end

selected_variants = table( ...
    VariantOrder, VariantRole, VariantKey, RuleKey, DecisionMode, RemoveDC, ...
    Stage1Rank, Stage1MeanSSIM, Stage1MeanPSNR);
end

function tbl = load_ranking_table(ranking_path)
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
tbl = sortrows(tbl, {'Rank', 'MeanSSIM', 'MeanPSNR'}, {'ascend', 'descend', 'descend'});
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
for k = 1:n
    info = lookup_or_parse_variant(keys(k), ranking_tbl);
    Rank(k) = info.Stage1Rank;
    RuleKey(k) = info.RuleKey;
    DecisionMode(k) = info.DecisionMode;
    RemoveDC(k) = info.RemoveDC;
    VariantKey(k) = info.VariantKey;
    MeanSSIM(k) = info.Stage1MeanSSIM;
    MeanPSNR(k) = info.Stage1MeanPSNR;
end
cand_tbl = table(Rank, RuleKey, DecisionMode, RemoveDC, VariantKey, MeanSSIM, MeanPSNR);
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
        "Stage1MeanPSNR", double(ranking_tbl.MeanPSNR(idx)));
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
    "Stage1MeanPSNR", NaN);
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
