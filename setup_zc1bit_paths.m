function directory_guard = setup_zc1bit_paths(cfg)
%SETUP_ZC1BIT_PATHS 仅注入本仓库根目录与 src/ 各模块路径，禁止跨仓库 addpath。
% 若接收返回值 directory_guard，则在作用期间将工作目录切换至本仓根路径，并在退出或报错后恢复原目录。

if nargin == 0 || ~isfield(cfg, "repo_root")
    repo_root = string(fileparts(mfilename("fullpath")));
else
    repo_root = string(cfg.repo_root);
end

if nargout > 0
    previous_directory = pwd;
    cd(repo_root);
    directory_guard = onCleanup(@() cd(previous_directory));
end

addpath(repo_root);
for name = ["workflow", "quantization", "threshold", "sampling", "sarproc", "metrics", "compat"]
    addpath(fullfile(repo_root, "src", name));
end
end
