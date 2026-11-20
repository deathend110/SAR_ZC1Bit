function varargout = run_zc1bit_search(varargin) %#ok<STOUT>
%RUN_ZC1BIT_SEARCH 阶段二（校准集独立阈值参数搜索）公共入口。
% 当前仓库已完成阶段一快筛闭环；阶段二待阶段一选出 Top 2-3 规则后按协议推进实现。

if nargin >= 1
    cfg = varargin{1};
    directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>
end
error("zc1bit:StageNotImplemented", ...
    "阶段二参数搜索 (search) 尚未启用，请先完成阶段一快筛 (run_zc1bit_screen) 并锁定 Top 2-3 规则。");
end
