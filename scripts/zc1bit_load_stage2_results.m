function D = zc1bit_load_stage2_results(cfg, p, record)
%ZC1BIT_LOAD_STAGE2_RESULTS 读取阶段二绘图所需的全部表格并做完整性校验。
% 只读 search/locked/evaluate/stats 下的 CSV 与 confirm_results.mat，不加载回波。

if nargin < 1 || isempty(cfg)
    cfg = zc1bit_config();
end
if nargin < 2 || isempty(p)
    p = zc1bit_paths(cfg);
end
if nargin < 3 || isempty(record)
    record_path = fullfile(p.protocol, "run_record.mat");
    loaded_record = load(record_path, "record");
    record = loaded_record.record;
end

paths = struct();
paths.selected = fullfile(p.search, "search_selected_variants.csv");
paths.search_best = fullfile(p.search, "search_best_params.csv");
paths.confirm_mat = fullfile(p.locked, "confirm_results.mat");
paths.locked_plan = fullfile(p.locked, "locked_confirm_plan.csv");
paths.eval_summary = fullfile(p.evaluate, "evaluate_summary.csv");
paths.eval_detail = fullfile(p.evaluate, "evaluate_detail.csv");
paths.wilcoxon = fullfile(p.stats, "paired_wilcoxon_vs_b0.csv");
paths.leakage = fullfile(p.evaluate, "evaluate_leakage_summary.csv");

required_files = [ ...
    paths.search_best, paths.confirm_mat, paths.locked_plan, ...
    paths.eval_summary, paths.eval_detail, paths.wilcoxon];
for k = 1:numel(required_files)
    if ~isfile(required_files(k))
        error("zc1bit:PlotMissingStage2", ...
            "缺少阶段二产物：%s。请先完成 run_zc1bit_search → run_zc1bit_confirm → run_zc1bit_evaluate。", ...
            required_files(k));
    end
end

confirm_loaded = load(paths.confirm_mat, "units", "selected_variants");
if ~isfield(confirm_loaded, "units") || ~isfield(confirm_loaded, "selected_variants")
    error("zc1bit:PlotInvalidConfirm", "confirm_results.mat 缺少 units 或 selected_variants：%s", paths.confirm_mat);
end

D = struct();
D.cfg = cfg;
D.p = p;
D.record = record;
D.paths = paths;
D.units = normalize_table(confirm_loaded.units);
D.units = sortrows(D.units, "UnitIdx");

if isfile(paths.selected)
    D.variants = read_csv_table(paths.selected, ...
        ["VariantOrder", "VariantRole", "VariantKey", "RuleKey", "DecisionMode", "RemoveDC"], ...
        "search_selected_variants.csv");
else
    D.variants = normalize_table(confirm_loaded.selected_variants);
end
D.variants = sortrows(D.variants, "VariantOrder");

D.search_best = read_csv_table(paths.search_best, ...
    ["UnitKey", "Branch", "Axis", "GroupName", "Q", "ThresholdType", ...
    "VariantKey", "RuleKey", "DecisionMode", "RemoveDC", ...
    "ParamTag", "As", "STR_dB", "FrOverBr", "FaOverBa", ...
    "Calibration_SSIM_Mean", "Calibration_PSNR_Mean"], ...
    "search_best_params.csv");

D.locked_plan = read_csv_table(paths.locked_plan, ...
    ["LockedConfigKey", "UnitKey", "Branch", "Axis", "GroupName", "Q", ...
    "ThresholdType", "Track", "VariantKey", "RuleKey", "DecisionMode", "RemoveDC", ...
    "ParamTag", "ParamKey", "As", "STR_dB", "FrOverBr", "FaOverBa", "Phi0"], ...
    "locked_confirm_plan.csv");

D.eval_summary = read_csv_table(paths.eval_summary, ...
    ["LockedConfigKey", "UnitKey", "Branch", "Axis", "GroupName", "Q", ...
    "ThresholdType", "Track", "VariantKey", "RuleKey", "DecisionMode", "RemoveDC", ...
    "ParamTag", "As", "STR_dB", "FrOverBr", "FaOverBa", ...
    "PSNR_Mean", "SSIM_Mean", ...
    "ENL_Mean", "GT_ENL_Mean", "AbsENLDiff_Mean", ...
    "ZeroRatio_Mean", "AbsMeanBias_Mean"], ...
    "evaluate_summary.csv");

D.eval_detail = read_csv_table(paths.eval_detail, ...
    ["SampleID", "LocalSampleIdx", "Dataset", "UnitKey", "Branch", "Q", ...
    "ThresholdType", "Track", "VariantKey", "RuleKey", "DecisionMode", "RemoveDC", ...
    "PSNR", "SSIM", "ENL", "GT_ENL", "AbsENLDiff", "ZeroRatio", "AbsMeanBias"], ...
    "evaluate_detail.csv");

D.wilcoxon = read_csv_table(paths.wilcoxon, ...
    ["ComparisonType", "UnitKey", "Branch", "Axis", "GroupName", "Q", ...
    "ThresholdType", "Track", "VariantKey", "RuleKey", ...
    "DeltaPSNR_Mean", "P_PSNR", "DeltaSSIM_Mean", "P_SSIM", ...
    "DeltaAbsENLDiff_Mean", "P_AbsENLDiff"], ...
    "paired_wilcoxon_vs_b0.csv");

if isfile(paths.leakage)
    D.leakage = read_csv_table(paths.leakage, ...
        ["SampleID", "Dataset", "UnitKey", "GroupName", "ThresholdType", ...
        "Track", "VariantKey", "OffSupportRatio"], ...
        "evaluate_leakage_summary.csv");
else
    D.leakage = table();
end

D.baseline_key = "B0_par_raw";
if any(D.variants.VariantRole == "baseline")
    idx = find(D.variants.VariantRole == "baseline", 1);
    D.baseline_key = D.variants.VariantKey(idx);
end
if ~any(D.variants.VariantKey == D.baseline_key)
    error("zc1bit:PlotMissingBaseline", "参评变体中缺少 baseline 变体：%s", D.baseline_key);
end
D.pseudo_keys = D.variants.VariantKey(D.variants.VariantRole ~= "baseline");
D.pseudo_keys = D.pseudo_keys(:);
D.variant_keys = [D.baseline_key; D.pseudo_keys];

D.eval_ind = D.eval_summary(D.eval_summary.Track == "independent_recalibrated", :);
D.eval_detail_ind = D.eval_detail(D.eval_detail.Track == "independent_recalibrated", :);
D.wilcoxon_ind = D.wilcoxon(D.wilcoxon.Track == "independent_recalibrated" & ...
    D.wilcoxon.ComparisonType == "vs_B0", :);

validate_stage2_combinations(D);
end

function T = read_csv_table(path_value, required_cols, label)
if ~isfile(path_value)
    error("zc1bit:PlotMissingFile", "缺少绘图输入文件：%s（%s）。", path_value, label);
end
T = readtable(path_value, "TextType", "string");
T = normalize_table(T);
vars = string(T.Properties.VariableNames);
missing = required_cols(~ismember(required_cols, vars));
if ~isempty(missing)
    error("zc1bit:PlotMissingColumn", "%s 缺少必要列：%s", label, join(missing, ", "));
end
end

function T = normalize_table(T)
if ~istable(T)
    error("zc1bit:PlotInvalidTable", "绘图输入必须是 table。");
end
vars = string(T.Properties.VariableNames);
string_cols = ["UnitKey", "Branch", "Axis", "GroupName", "ThresholdType", ...
    "VariantRole", "VariantKey", "RuleKey", "DecisionMode", "Track", ...
    "ComparisonType", "ParamTag", "ParamKey", "LockedConfigKey", ...
    "Dataset", "EchoVariable", "File", "Split"];
for name = string_cols
    if ismember(name, vars) && ~isstring(T.(name))
        T.(name) = string(T.(name));
    end
end
if ismember("RemoveDC", vars) && ~islogical(T.RemoveDC)
    T.RemoveDC = logical(T.RemoveDC);
end
bool_cols = ["TouchesRTBoundary", "TouchesSTRBoundary", "TouchesFreqBoundary", ...
    "TouchesAnyBoundary", "SameAsB0Optimal", "IsNoOp"];
for name = bool_cols
    if ismember(name, vars) && ~islogical(T.(name))
        T.(name) = logical(T.(name));
    end
end
if ismember("VariantOrder", vars) && ~isnumeric(T.VariantOrder)
    T.VariantOrder = double(T.VariantOrder);
end
end

function validate_stage2_combinations(D)
n_expected = height(D.units) * height(D.variants);
if height(D.search_best) ~= n_expected
    error("zc1bit:PlotIncompleteSearch", ...
        "search_best_params.csv 行数为 %d，预期 %d（%d 单元 × %d 变体），阶段二搜索可能未完成。", ...
        height(D.search_best), n_expected, height(D.units), height(D.variants));
end
if height(D.eval_ind) ~= n_expected
    error("zc1bit:PlotIncompleteEvaluate", ...
        "evaluate_summary.csv 独立评估行数为 %d，预期 %d，阶段二评估可能未完成。", ...
        height(D.eval_ind), n_expected);
end
if isempty(D.eval_detail_ind)
    error("zc1bit:PlotEmptyEvaluateDetail", "evaluate_detail.csv 不含 independent_recalibrated 行。");
end
for u_idx = 1:height(D.units)
    ukey = D.units.UnitKey(u_idx);
    for v_idx = 1:height(D.variants)
        vkey = D.variants.VariantKey(v_idx);
        if ~any(D.search_best.UnitKey == ukey & D.search_best.VariantKey == vkey)
            error("zc1bit:PlotIncompleteSearch", "搜索结果缺少组合：%s / %s", ukey, vkey);
        end
        if ~any(D.eval_ind.UnitKey == ukey & D.eval_ind.VariantKey == vkey)
            error("zc1bit:PlotIncompleteEvaluate", "测试评估缺少组合：%s / %s", ukey, vkey);
        end
    end
end
end
