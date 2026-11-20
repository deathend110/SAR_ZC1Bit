function varargout = run_zc1bit_evaluate(varargin)
%RUN_ZC1BIT_EVALUATE 阶段二（35 测试样本最终评估与配对统计）公共入口。

if nargin == 0
    cfg = zc1bit_config();
else
    cfg = varargin{1};
end
directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>
[varargout{1:nargout}] = zc1bit_evaluate_impl(cfg);
end
