function output = zc1bit_workflow_impl(cfg)
%ZC1BIT_WORKFLOW_IMPL 统一工作流调度入口：默认仅预检，阶段一支持 preflight / screen。

if nargin == 0
    cfg = zc1bit_config();
end
allowed = ["preflight", "screen", "search", "confirm", "evaluate", "all"];
stage = string(cfg.stage);
if ~isscalar(stage) || ~ismember(stage, allowed)
    error("zc1bit:UnknownStage", "未知执行阶段：%s", stage);
end

directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>
record = zc1bit_prepare_run(cfg);
p = zc1bit_paths(cfg);
output = struct("paths", p, "record", record);

switch stage
    case "preflight"
        return;
    case "screen"
        output.screen = zc1bit_screen_impl(cfg);
    case {"search", "confirm", "evaluate", "all"}
        error("zc1bit:StageNotImplemented", ...
            "阶段 %s 属于阶段二全量校准与评估，待阶段一快筛选出 Top 2-3 规则后按协议推进实现。", stage);
end
end
