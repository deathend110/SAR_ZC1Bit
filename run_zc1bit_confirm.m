function varargout = run_zc1bit_confirm(varargin)
%RUN_ZC1BIT_CONFIRM 阶段二（参数锁定确认与原参数直迁消融表）公共入口。

if nargin == 0
    cfg = zc1bit_config();
else
    cfg = varargin{1};
end
directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>
[varargout{1:nargout}] = zc1bit_confirm_impl(cfg);
end
