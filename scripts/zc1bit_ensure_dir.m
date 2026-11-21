function zc1bit_ensure_dir(path_value)
%ZC1BIT_ENSURE_DIR 确保目录存在；用于绘图脚本输出目录。

path_str = string(path_value);
if strlength(path_str) == 0
    return;
end
if ~isfolder(path_str)
    [ok, msg] = mkdir(path_str);
    if ~ok
        error("zc1bit:PlotMkdirFailed", "无法创建输出目录 %s：%s", path_str, msg);
    end
end
end
