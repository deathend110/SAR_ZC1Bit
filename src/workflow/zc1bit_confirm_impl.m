function out = zc1bit_confirm_impl(cfg)
%ZC1BIT_CONFIRM_IMPL 阶段二（校准锁参确认与原参数直迁消融表构造）。
% 严格要求阶段二 search cube 已 100% 完成，绝不从未完成搜索结果中提前锁参。

if nargin == 0
    cfg = zc1bit_config();
end
directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>

record = zc1bit_prepare_run(cfg);
p = zc1bit_paths(cfg);
S60 = load(cfg.parameter_file_60);

search_mat = fullfile(p.search, "search_results.mat");
if ~isfile(search_mat)
    error("zc1bit:MissingSearchResults", ...
        "未找到阶段二搜索结果：%s。请先完成 run_zc1bit_search(cfg)。", search_mat);
end
loaded = load(search_mat, "signature", "units", "selected_variants", "best_params_table");
if ~isfield(loaded, "best_params_table") || ~isfield(loaded, "selected_variants")
    error("zc1bit:InvalidSearchResults", "search_results.mat 缺少必要字段。");
end

cal_manifest = record.manifest(record.manifest.Split == "calibration", :);
expected_sig = struct( ...
    "stage2", zc1bit_parameters(cfg.stage2), ...
    "nyquist_margin", cfg.nyquist_margin, ...
    "seed", cfg.seed, ...
    "S60", zc1bit_parameters(S60), ...
    "cal_manifest", zc1bit_parameters(cal_manifest));
if isfield(loaded, "signature")
    zc1bit_check_parameters(loaded.signature, expected_sig, "阶段二搜索结果与当前配置");
end

units = zc1bit_build_stage2_units(cfg, S60);
selected_variants = loaded.selected_variants;
best_tbl = loaded.best_params_table;

verify_full_cube_completeness(units, selected_variants, best_tbl);

locked_params = best_tbl;
locked_params.Track = repmat("independent_recalibrated", height(locked_params), 1);

plan_parts = cell(2, 1);
plan_parts{1} = add_confirm_keys(locked_params, best_tbl, selected_variants);

include_transfer = true;
if isfield(cfg.stage2, "include_b0_transfer")
    include_transfer = logical(cfg.stage2.include_b0_transfer);
end

base_key = string(cfg.stage2.baseline_variant);
if include_transfer
    plan_parts{2} = build_b0_transfer_rows(units, selected_variants, best_tbl, base_key);
    locked_confirm_plan = vertcat(plan_parts{:});
else
    locked_confirm_plan = plan_parts{1};
end

if numel(unique(locked_confirm_plan.LockedConfigKey)) ~= height(locked_confirm_plan)
    error("zc1bit:LockedKeyCollision", "锁参评估计划存在重复的 LockedConfigKey。");
end

if ~isfolder(p.locked)
    mkdir(p.locked);
end
writetable(locked_params, fullfile(p.locked, "locked_params.csv"));
writetable(locked_confirm_plan, fullfile(p.locked, "locked_confirm_plan.csv"));

signature = struct( ...
    "stage2", zc1bit_parameters(cfg.stage2), ...
    "nyquist_margin", cfg.nyquist_margin, ...
    "seed", cfg.seed, ...
    "S60", zc1bit_parameters(S60), ...
    "cal_manifest", zc1bit_parameters(cal_manifest), ...
    "selected_variants", zc1bit_parameters(selected_variants));

zc1bit_atomic_save(fullfile(p.locked, "confirm_results.mat"), struct( ...
    "signature", signature, ...
    "units", units, ...
    "selected_variants", selected_variants, ...
    "locked_params", locked_params, ...
    "locked_confirm_plan", locked_confirm_plan));
write_confirm_metadata(p.locked, units, selected_variants, locked_params, locked_confirm_plan);

fprintf("阶段二参数锁定确认完成（独立寻优 %d 条，全评估计划 %d 条），输出目录：%s\n", ...
    height(locked_params), height(locked_confirm_plan), p.locked);

out = struct( ...
    "selected_variants", selected_variants, ...
    "units", units, ...
    "locked_params", locked_params, ...
    "locked_confirm_plan", locked_confirm_plan, ...
    "output_dir", p.locked);
end

function verify_full_cube_completeness(units, selected_variants, best_tbl)
n_u = height(units);
n_v = height(selected_variants);
if height(best_tbl) ~= n_u * n_v
    error("zc1bit:IncompleteSearchCube", ...
        "搜索最优参数行数为 %d，未覆盖完整的 %d×%d=%d 单元组合，禁止提前锁参。", ...
        height(best_tbl), n_u, n_v, n_u * n_v);
end
for u_idx = 1:n_u
    ukey = units.UnitKey(u_idx);
    for v_idx = 1:n_v
        vkey = selected_variants.VariantKey(v_idx);
        mask = best_tbl.UnitKey == ukey & best_tbl.VariantKey == vkey;
        if sum(mask) ~= 1
            error("zc1bit:IncompleteSearchCube", "缺失或重复组合：%s / %s。", ukey, vkey);
        end
        if ~isfinite(best_tbl.Calibration_SSIM_Mean(mask)) || ~isfinite(best_tbl.Calibration_PSNR_Mean(mask))
            error("zc1bit:IncompleteSearchCube", "组合 %s / %s 包含非有限校准指标。", ukey, vkey);
        end
    end
end
end

function tbl = add_confirm_keys(tbl, b0_ref_tbl, selected_variants)
n = height(tbl);
tbl.LockedConfigKey = tbl.UnitKey + "__" + tbl.Track + "__" + tbl.VariantKey;
BaselineParamTag = strings(n, 1);
SameAsB0Optimal = false(n, 1);

base_row = selected_variants(selected_variants.VariantRole == "baseline", :);
base_key = base_row.VariantKey(1);

for k = 1:n
    bmask = b0_ref_tbl.UnitKey == tbl.UnitKey(k) & b0_ref_tbl.VariantKey == base_key;
    b_tag = b0_ref_tbl.ParamTag(find(bmask, 1));
    BaselineParamTag(k) = b_tag;
    SameAsB0Optimal(k) = (tbl.ParamTag(k) == b_tag);
end
tbl.BaselineParamTag = BaselineParamTag;
tbl.SameAsB0Optimal = SameAsB0Optimal;
end

function tr_tbl = build_b0_transfer_rows(units, selected_variants, best_tbl, base_key)
cand_rows = selected_variants(selected_variants.VariantKey ~= base_key, :);
param_units = units(units.ThresholdType ~= "ZT", :);
n_pu = height(param_units);
n_c = height(cand_rows);
if n_pu == 0 || n_c == 0
    tr_tbl = best_tbl([], :);
    tr_tbl.Track = strings(0, 1);
    tr_tbl.LockedConfigKey = strings(0, 1);
    tr_tbl.BaselineParamTag = strings(0, 1);
    tr_tbl.SameAsB0Optimal = false(0, 1);
    return;
end

parts = cell(n_pu * n_c, 1);
ptr = 0;
for u_idx = 1:n_pu
    ukey = param_units.UnitKey(u_idx);
    b0_row = best_tbl(best_tbl.UnitKey == ukey & best_tbl.VariantKey == base_key, :);
    b0_row = b0_row(1, :);
    for v_idx = 1:n_c
        vrow = cand_rows(v_idx, :);
        row = b0_row;
        row.VariantOrder = vrow.VariantOrder;
        row.VariantRole = vrow.VariantRole;
        row.VariantKey = vrow.VariantKey;
        row.RuleKey = vrow.RuleKey;
        row.DecisionMode = vrow.DecisionMode;
        row.RemoveDC = logical(vrow.RemoveDC);
        row.Track = "b0_transferred";
        row.LockedConfigKey = ukey + "__b0_transferred__" + vrow.VariantKey;
        row.BaselineParamTag = b0_row.ParamTag;
        row.SameAsB0Optimal = true;
        ptr = ptr + 1;
        parts{ptr} = row;
    end
end
tr_tbl = vertcat(parts{:});
end

function write_confirm_metadata(locked_dir, units, selected_variants, locked_params, locked_confirm_plan)
fid = fopen(fullfile(locked_dir, "confirm_metadata.txt"), "w", "n", "UTF-8");
if fid < 0
    error("zc1bit:MetadataWrite", "无法写入 confirm_metadata.txt。");
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "Experiment=SAR-ZC1Bit Stage 2 Confirm\n");
fprintf(fid, "UnitCount=%d\n", height(units));
fprintf(fid, "SelectedVariants=%s\n", join(selected_variants.VariantKey, ","));
fprintf(fid, "IndependentLockedCount=%d\n", height(locked_params));
fprintf(fid, "FullEvalPlanCount=%d\n", height(locked_confirm_plan));
end
