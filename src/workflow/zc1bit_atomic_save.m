function zc1bit_atomic_save(path_value, payload)
%ZC1BIT_ATOMIC_SAVE 先创建父目录并写入同目录临时文件，成功后再原子替换。

path_str = string(path_value);
folder = string(fileparts(path_str));
if strlength(folder) > 0 && ~isfolder(folder)
    [ok, msg] = mkdir(folder);
    if ~ok
        error("zc1bit:MkdirFailed", "无法创建目录 %s：%s", folder, msg);
    end
end

scratch = path_str + ".tmp.mat";
save(scratch, "-struct", "payload", "-v7.3");
[ok, msg] = movefile(scratch, path_str, "f");
if ~ok
    error("zc1bit:AtomicSave", "原子保存检查点失败 (%s)：%s", path_str, msg);
end
end
