function varargout = run_zc1bit_ta_sft(cfg, action)
%RUN_ZC1BIT_TA_SFT 阈值感知 SFT 1-bit 接收端重构（TA-SFT）专题一键入口。
% 用法：
%   out = run_zc1bit_ta_sft();                  % 默认使用 zc1bit_config() 执行 search + evaluate
%   out = run_zc1bit_ta_sft(cfg, "search");     % 仅执行 35 校准样本超参寻优
%   out = run_zc1bit_ta_sft(cfg, "evaluate");   % 仅执行 35 测试样本评估、泄漏与相偏鲁棒性
%   out = run_zc1bit_ta_sft(cfg, "all");        % 顺序执行 search -> evaluate

if nargin < 1 || isempty(cfg)
    cfg = zc1bit_config();
end
if nargin < 2 || isempty(action)
    action = "all";
end
action = string(action);
allowed = ["search", "evaluate", "all"];
if ~isscalar(action) || ~ismember(action, allowed)
    error("zc1bit:UnknownTaSftAction", "未知 TA-SFT 执行动作：%s（支持 search / evaluate / all）。", action);
end

directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>
out = struct();
switch action
    case "search"
        out.search = zc1bit_ta_sft_search_impl(cfg);
    case "evaluate"
        out.evaluate = zc1bit_ta_sft_evaluate_impl(cfg);
    case "all"
        out.search = zc1bit_ta_sft_search_impl(cfg);
        out.evaluate = zc1bit_ta_sft_evaluate_impl(cfg);
end
if nargout > 0
    varargout{1} = out;
end
end
