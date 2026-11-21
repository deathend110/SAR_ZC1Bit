function [cfg, p, S60, record] = zc1bit_script_begin()
%ZC1BIT_SCRIPT_BEGIN 绘图/成像脚本公共初始化。
% 只注入本仓库路径、读取 run_record 与 FS60 参数，不强制加载外部回波；
% render_* 脚本在真正读取回波前再校验数据根目录并刷新 FilePath。

script_dir = string(fileparts(mfilename("fullpath")));
repo_root = string(fileparts(script_dir));
addpath(char(repo_root));
addpath(char(script_dir));

cfg = zc1bit_config();
setup_zc1bit_paths(cfg);
p = zc1bit_paths(cfg);

record_path = fullfile(p.protocol, "run_record.mat");
if ~isfile(record_path)
    error("zc1bit:PlotMissingRunRecord", ...
        "未找到运行记录：%s。请先完成 preflight 与阶段二。", record_path);
end
loaded = load(record_path, "record");
if ~isfield(loaded, "record") || ~isfield(loaded.record, "manifest")
    error("zc1bit:PlotInvalidRunRecord", "运行记录缺少 record.manifest：%s", record_path);
end
record = loaded.record;

if ~isfile(cfg.parameter_file_60)
    error("zc1bit:PlotMissingS60", "未找到成像参数文件：%s", cfg.parameter_file_60);
end
S60 = load(cfg.parameter_file_60);

zc1bit_ensure_dir(p.figures);
zc1bit_ensure_dir(fullfile(p.figures, "curves"));
zc1bit_ensure_dir(fullfile(p.figures, "images"));
zc1bit_ensure_dir(fullfile(p.figures, "tables"));
end
