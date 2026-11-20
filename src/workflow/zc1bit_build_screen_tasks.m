function [threshold_groups, variant_specs] = zc1bit_build_screen_tasks(cfg, S60)
%ZC1BIT_BUILD_SCREEN_TASKS 构建阶段一规则快筛的阈值配置表与规则/开关变体表。
% 覆盖：
%   - range_branch（dim 1）：R1A1, R4A1, R9A1 × {ZT, RT, RSFT, 2D-SFT}
%   - azimuth_branch（dim 2）：R1A1, R1A4, R1A9 × {ZT, RT, ASFT, 2D-SFT}
%   - 规则变体：B0 + R1..R9 + L1 × {parallel, serial} × {raw, dc_removed}

sc = cfg.screen;
rows_tg = cell(0, 14);

for b_idx = 1:numel(sc.branches)
    branch = string(sc.branches(b_idx));
    if branch == "range_branch"
        axis_mode = "range";
        allocs = sc.range_allocations;
        ttypes = string(sc.range_threshold_types);
        branch_tag = "rb";
    elseif branch == "azimuth_branch"
        axis_mode = "azimuth";
        allocs = sc.azimuth_allocations;
        ttypes = string(sc.azimuth_threshold_types);
        branch_tag = "ab";
    else
        error("zc1bit:UnknownBranch", "未知阶段一分支：%s", branch);
    end

    for a_idx = 1:size(allocs, 1)
        rq = allocs(a_idx, 1);
        aq = allocs(a_idx, 2);
        Q = rq * aq;
        group_name = sprintf("R%dA%d", rq, aq);
        grid = resolve_sampling_grid(rq, aq, [S60.nrn, S60.nan], S60);

        for t_idx = 1:numel(ttypes)
            ttype = ttypes(t_idx);
            switch ttype
                case "ZT"
                    tkey = sprintf("%s_%s_%s_zt", branch_tag, group_name, axis_mode);
                    rows_tg(end + 1, :) = {branch, axis_mode, group_name, Q, rq, aq, ...
                        ttype, "zt", NaN, NaN, NaN, NaN, sc.phi0, string(tkey)}; %#ok<AGROW>

                case "RT"
                    for as_val = sc.rt_as_values(:).'
                        as_tag = num_to_key(as_val);
                        param_tag = "as" + as_tag;
                        tkey = sprintf("%s_%s_%s_rt_%s", branch_tag, group_name, axis_mode, param_tag);
                        rows_tg(end + 1, :) = {branch, axis_mode, group_name, Q, rq, aq, ...
                            ttype, param_tag, as_val, NaN, NaN, NaN, sc.phi0, string(tkey)}; %#ok<AGROW>
                    end

                case {"RSFT", "ASFT", "2D-SFT"}
                    base = lookup_sft_baseline(sc.sft_baselines, axis_mode, rq, aq, ttype);
                    for d_str = sc.sft_str_offsets(:).'
                        str_val = base.STR_dB + d_str;
                        fr_val = base.FrOverBr;
                        fa_val = base.FaOverBa;
                        spec = struct("fr_over_Br", fr_val, "fa_over_Ba", fa_val);
                        check_sft_nyquist(grid, spec, S60, cfg.nyquist_margin);

                        param_tag = sprintf("s%s_fr%s_fa%s", ...
                            num_to_key(str_val), num_to_key(fr_val), num_to_key(fa_val));
                        type_key = lower(replace(ttype, "-", ""));
                        tkey = sprintf("%s_%s_%s_%s_%s", ...
                            branch_tag, group_name, axis_mode, type_key, param_tag);
                        rows_tg(end + 1, :) = {branch, axis_mode, group_name, Q, rq, aq, ...
                            ttype, string(param_tag), NaN, str_val, fr_val, fa_val, sc.phi0, string(tkey)}; %#ok<AGROW>
                    end

                otherwise
                    error("zc1bit:ThresholdType", "不支持的阶段一阈值类型：%s", ttype);
            end
        end
    end
end

if isempty(rows_tg)
    error("zc1bit:EmptyThresholdGroups", "阶段一未生成任何阈值配置。");
end

threshold_groups = cell2table(rows_tg, "VariableNames", { ...
    'Branch', 'Axis', 'GroupName', 'Q', 'Range_q', 'Azimuth_q', ...
    'ThresholdType', 'ParamTag', 'As', 'STR_dB', 'FrOverBr', 'FaOverBa', 'Phi0', 'ThresholdKey'});
threshold_groups.Branch = string(threshold_groups.Branch);
threshold_groups.Axis = string(threshold_groups.Axis);
threshold_groups.GroupName = string(threshold_groups.GroupName);
threshold_groups.Q = double(threshold_groups.Q);
threshold_groups.Range_q = double(threshold_groups.Range_q);
threshold_groups.Azimuth_q = double(threshold_groups.Azimuth_q);
threshold_groups.ThresholdType = string(threshold_groups.ThresholdType);
threshold_groups.ParamTag = string(threshold_groups.ParamTag);
threshold_groups.As = double(threshold_groups.As);
threshold_groups.STR_dB = double(threshold_groups.STR_dB);
threshold_groups.FrOverBr = double(threshold_groups.FrOverBr);
threshold_groups.FaOverBa = double(threshold_groups.FaOverBa);
threshold_groups.Phi0 = double(threshold_groups.Phi0);
threshold_groups.ThresholdKey = string(threshold_groups.ThresholdKey);

if numel(unique(threshold_groups.ThresholdKey)) ~= height(threshold_groups)
    error("zc1bit:ThresholdKeyCollision", "阶段一存在重复的 ThresholdKey。");
end

% 规则 × 正交开关变体表
rows_var = cell(0, 4);
rules = string(sc.rules);
modes = string(sc.decision_modes);
dc_opts = logical(sc.remove_dc_options);

for r_idx = 1:numel(rules)
    rkey = rules(r_idx);
    if rkey == "B0"
        % B0 无过零置零动作，并行与串行完全等价，仅保留 parallel。
        b0_dc_opts = false;
        if isfield(sc, "include_b0_dc") && sc.include_b0_dc
            b0_dc_opts = dc_opts;
        end
        for dc = b0_dc_opts(:).'
            vkey = make_variant_key(rkey, "parallel", dc);
            rows_var(end + 1, :) = {rkey, "parallel", dc, vkey}; %#ok<AGROW>
        end
        continue;
    end

    for m_idx = 1:numel(modes)
        mode = modes(m_idx);
        for dc = dc_opts(:).'
            vkey = make_variant_key(rkey, mode, dc);
            rows_var(end + 1, :) = {rkey, mode, dc, vkey}; %#ok<AGROW>
        end
    end
end

if isempty(rows_var)
    error("zc1bit:EmptyVariants", "阶段一未生成任何规则变体。");
end

variant_specs = cell2table(rows_var, "VariableNames", {'RuleKey', 'DecisionMode', 'RemoveDC', 'VariantKey'});
variant_specs.RuleKey = string(variant_specs.RuleKey);
variant_specs.DecisionMode = string(variant_specs.DecisionMode);
variant_specs.RemoveDC = logical(variant_specs.RemoveDC);
variant_specs.VariantKey = string(variant_specs.VariantKey);

if numel(unique(variant_specs.VariantKey)) ~= height(variant_specs)
    error("zc1bit:VariantKeyCollision", "阶段一存在重复的 VariantKey。");
end
end

function base = lookup_sft_baseline(baselines, axis_mode, rq, aq, ttype)
mask = string({baselines.Axis}).' == axis_mode & ...
    [baselines.Range_q].' == rq & ...
    [baselines.Azimuth_q].' == aq & ...
    string({baselines.ThresholdType}).' == ttype;
if sum(mask) ~= 1
    error("zc1bit:MissingSFTBaseline", ...
        "阶段一缺少唯一 SFT 基准配置：Axis=%s, R=%d, A=%d, Type=%s。", axis_mode, rq, aq, ttype);
end
base = baselines(mask);
end

function key = make_variant_key(rule_key, mode, remove_dc)
mode_tag = "par";
if mode == "serial"
    mode_tag = "ser";
end
dc_tag = "raw";
if remove_dc
    dc_tag = "dc";
end
key = string(sprintf("%s_%s_%s", rule_key, mode_tag, dc_tag));
end

function text = num_to_key(value)
text = string(sprintf("%.6g", value));
text = replace(text, "-", "m");
text = replace(text, ".", "p");
end
