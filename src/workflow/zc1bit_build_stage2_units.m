function units = zc1bit_build_stage2_units(cfg, S60)
%ZC1BIT_BUILD_STAGE2_UNITS 构建阶段二单向主线单元表（range R1A1..R9A1 与 azimuth R1A1..R1A9 × 4 类阈值）。
% 严格拒绝双向 BARU（Range_q>1 && Azimuth_q>1）与 2D-SFT。

s2 = cfg.stage2;
branches = string(s2.branches);

est_count = 0;
for b_idx = 1:numel(branches)
    br = branches(b_idx);
    if br == "range_branch"
        est_count = est_count + size(s2.range_allocations, 1) * numel(s2.range_threshold_types);
    elseif br == "azimuth_branch"
        est_count = est_count + size(s2.azimuth_allocations, 1) * numel(s2.azimuth_threshold_types);
    else
        error("zc1bit:UnknownBranch", "未知阶段二分支：%s", br);
    end
end
if est_count == 0
    error("zc1bit:EmptyStage2Units", "阶段二未配置任何实验单元。");
end

UnitIdx = (1:est_count).';
Branch = strings(est_count, 1);
Axis = strings(est_count, 1);
GroupName = strings(est_count, 1);
Q = zeros(est_count, 1);
Range_q = zeros(est_count, 1);
Azimuth_q = zeros(est_count, 1);
ThresholdType = strings(est_count, 1);
UnitKey = strings(est_count, 1);

ptr = 0;
for b_idx = 1:numel(branches)
    br = branches(b_idx);
    if br == "range_branch"
        axis_mode = "range";
        allocs = double(s2.range_allocations);
        ttypes = string(s2.range_threshold_types);
        allowed_tt = ["ZT", "RT", "RSFT"];
        btag = "rb";
    else
        axis_mode = "azimuth";
        allocs = double(s2.azimuth_allocations);
        ttypes = string(s2.azimuth_threshold_types);
        allowed_tt = ["ZT", "RT", "ASFT"];
        btag = "ab";
    end

    if any(allocs(:, 1) > 1 & allocs(:, 2) > 1)
        error("zc1bit:Stage2BarNotAllowed", "阶段二主线禁止双向 BARU 倍率对，仅允许单向距离或单向方位上采样。");
    end
    if br == "range_branch" && any(allocs(:, 2) ~= 1)
        error("zc1bit:Stage2AxisMismatch", "range_branch 要求 Azimuth_q == 1。");
    end
    if br == "azimuth_branch" && any(allocs(:, 1) ~= 1)
        error("zc1bit:Stage2AxisMismatch", "azimuth_branch 要求 Range_q == 1。");
    end
    if any(~ismember(ttypes, allowed_tt))
        error("zc1bit:Stage2ThresholdNotAllowed", ...
            "阶段二 %s 分支仅支持 [%s]，禁止 2D-SFT。", br, join(allowed_tt, ", "));
    end

    for a_idx = 1:size(allocs, 1)
        rq = allocs(a_idx, 1);
        aq = allocs(a_idx, 2);
        resolve_sampling_grid(rq, aq, [S60.nrn, S60.nan], S60);
        gname = sprintf("R%dA%d", rq, aq);

        for t_idx = 1:numel(ttypes)
            ttype = ttypes(t_idx);
            type_key = lower(replace(ttype, "-", ""));
            ptr = ptr + 1;
            Branch(ptr) = br;
            Axis(ptr) = axis_mode;
            GroupName(ptr) = gname;
            Q(ptr) = rq * aq;
            Range_q(ptr) = rq;
            Azimuth_q(ptr) = aq;
            ThresholdType(ptr) = ttype;
            UnitKey(ptr) = sprintf("%s_%s_%s_%s", btag, gname, axis_mode, type_key);
        end
    end
end

units = table(UnitIdx, Branch, Axis, GroupName, Q, Range_q, Azimuth_q, ThresholdType, UnitKey);
if numel(unique(units.UnitKey)) ~= height(units)
    error("zc1bit:Stage2UnitKeyCollision", "阶段二存在重复的 UnitKey。");
end
end
