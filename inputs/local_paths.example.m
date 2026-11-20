function paths = local_paths()
%LOCAL_PATHS 本机 SAR 回波数据与成像参数路径配置模板。
% 使用前请复制本文件为 inputs/local_paths.m（已被 .gitignore 忽略），
% 并修改为本机真实路径；禁止在核心源码中硬编码个人绝对路径。

repo_root = string(fileparts(fileparts(mfilename("fullpath"))));

paths = struct();
paths.data_root = "G:\MATLAB-G\SAR Full PSF";
paths.parameter_file_60 = fullfile(repo_root, "inputs", "parameters", "FS60_params.mat");
paths.results_root = fullfile(repo_root, "results_zc1bit");
end
