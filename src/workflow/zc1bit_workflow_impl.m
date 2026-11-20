function output = zc1bit_workflow_impl(cfg)
%ZC1BIT_WORKFLOW_IMPL 统一工作流调度入口：支持 preflight / screen / search / confirm / evaluate / all。

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
    case "search"
        output.search = zc1bit_search_impl(cfg);
    case "confirm"
        output.confirm = zc1bit_confirm_impl(cfg);
    case "evaluate"
        output.evaluate = zc1bit_evaluate_impl(cfg);
    case "all"
        output.screen = zc1bit_screen_impl(cfg);
        output.search = zc1bit_search_impl(cfg);
        output.confirm = zc1bit_confirm_impl(cfg);
        output.evaluate = zc1bit_evaluate_impl(cfg);
end
end
