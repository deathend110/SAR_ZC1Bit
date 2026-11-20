function varargout = run_zc1bit_workflow(varargin)
%RUN_ZC1BIT_WORKFLOW SAR-ZC1Bit 统一工作流公共入口，具体实现位于 src/workflow。

if nargin == 0
    cfg = zc1bit_config();
else
    cfg = varargin{1};
end
directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>
[varargout{1:nargout}] = zc1bit_workflow_impl(cfg);
end
