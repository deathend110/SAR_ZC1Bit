function zc1bit_atomic_save(path_value, payload, mat_version)
%ZC1BIT_ATOMIC_SAVE 先创建父目录并写入同目录临时文件，成功后再原子替换。
% 默认对小于 1.5 GB 的载荷使用 -v7（减少机械盘 HDF5 元数据寻道开销），超限自动切换 -v7.3。

path_str = string(path_value);
folder = string(fileparts(path_str));
if strlength(folder) > 0 && ~isfolder(folder)
    [ok, msg] = mkdir(folder);
    if ~ok
        error("zc1bit:MkdirFailed", "无法创建目录 %s：%s", folder, msg);
    end
end

if nargin < 3 || strlength(string(mat_version)) == 0
    info = whos("payload");
    if ~isempty(info) && info.bytes >= 1.5e9
        mat_version = "-v7.3";
    else
        mat_version = "-v7";
    end
else
    mat_version = string(mat_version);
end

scratch = path_str + ".tmp.mat";
save(scratch, "-struct", "payload", char(mat_version));
[ok, msg] = movefile(scratch, path_str, "f");
if ~ok
    error("zc1bit:AtomicSave", "原子保存检查点失败 (%s)：%s", path_str, msg);
end
end
