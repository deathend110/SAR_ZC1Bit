function varargout = run_zc1bit_sft_p15(cfg, action)
%RUN_ZC1BIT_SFT_P15 严格伪 1.5-bit（{+1,0,-1}）SFT 专题实验一键入口。
% 用法：
%   out = run_zc1bit_sft_p15();                  % 默认使用 zc1bit_config() 执行 search + evaluate
%   out = run_zc1bit_sft_p15(cfg, "search");     % 仅执行 35 校准样本超参寻优
%   out = run_zc1bit_sft_p15(cfg, "evaluate");   % 仅执行 35 测试样本评估、泄漏/载波抑制与相偏消融
%   out = run_zc1bit_sft_p15(cfg, "all");        % 顺序执行 search -> evaluate

if nargin < 1 || isempty(cfg)
    cfg = zc1bit_config();
end
if nargin < 2 || isempty(action)
    action = "all";
end
action = string(action);
allowed = ["search", "evaluate", "all"];
if ~isscalar(action) || ~ismember(action, allowed)
    error("zc1bit:UnknownSftP15Action", "未知 SFT-P15 执行动作：%s（支持 search / evaluate / all）。", action);
end

directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>
out = struct();
switch action
    case "search"
        out.search = zc1bit_sft_p15_search_impl(cfg);
    case "evaluate"
        out.evaluate = zc1bit_sft_p15_evaluate_impl(cfg);
    case "all"
        out.search = zc1bit_sft_p15_search_impl(cfg);
        out.evaluate = zc1bit_sft_p15_evaluate_impl(cfg);
end
if nargout > 0
    varargout{1} = out;
end
end
