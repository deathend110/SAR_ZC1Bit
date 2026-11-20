function varargout = run_zc1bit_evaluate(varargin) %#ok<STOUT>
%RUN_ZC1BIT_EVALUATE 阶段二（35 测试样本最终评估与配对统计）公共入口。

if nargin >= 1
    cfg = varargin{1};
    directory_guard = setup_zc1bit_paths(cfg); %#ok<NASGU>
end
error("zc1bit:StageNotImplemented", ...
    "阶段二测试集评估 (evaluate) 尚未启用，严禁未锁定参数前直接评估测试集。");
end
