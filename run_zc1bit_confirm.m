function varargout = run_zc1bit_confirm(varargin) %#ok<STOUT>
%RUN_ZC1BIT_CONFIRM 阶段二（参数锁定确认）公共入口。

if nargin >= 1
    cfg = varargin{1};
    directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>
end
error("zc1bit:StageNotImplemented", ...
    "阶段二参数确认 (confirm) 尚未启用，请先完成阶段一快筛与阶段二参数搜索。");
end
