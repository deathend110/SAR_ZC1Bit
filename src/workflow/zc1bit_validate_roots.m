function zc1bit_validate_roots(cfg)
%ZC1BIT_VALIDATE_ROOTS 校验本机回波根目录与 FS60 参数文件是否有效。

if ~isfield(cfg, "data_root") || strlength(string(cfg.data_root)) == 0 || ~isfolder(cfg.data_root)
    error("zc1bit:MissingDataRoot", ...
        "原始 SAR 回波目录不存在：%s。请复制 inputs/local_paths.example.m 为 inputs/local_paths.m 并配置正确路径。", ...
        string(cfg.data_root));
end

if ~isfield(cfg, "parameter_file_60") || strlength(string(cfg.parameter_file_60)) == 0 || ~isfile(cfg.parameter_file_60)
    error("zc1bit:MissingFS60", ...
        "FS60 成像参数文件不存在：%s。请将 FS60_params.mat 放入 inputs/parameters/ 或在 inputs/local_paths.m 中指定。", ...
        string(cfg.parameter_file_60));
end
end
