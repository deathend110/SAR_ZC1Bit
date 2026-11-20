function varargout = run_zc1bit_search(varargin)
%RUN_ZC1BIT_SEARCH 阶段二（校准集独立阈值参数搜索）公共入口，仅使用 35 个校准样本。

if nargin == 0
    cfg = zc1bit_config();
else
    cfg = varargin{1};
end
directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>
[varargout{1:nargout}] = zc1bit_search_impl(cfg);
end
