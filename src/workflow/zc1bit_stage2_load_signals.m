function signals = zc1bit_stage2_load_signals(split_manifest, S60, decim_step, split_label)
%ZC1BIT_STAGE2_LOAD_SIGNALS 对指定样本集唯一回波 MAT 文件仅顺序读取 1 轮并常驻内存。

if nargin < 4 || strlength(string(split_label)) == 0
    split_label = "样本";
end
n_samples = height(split_manifest);
signals = cell(n_samples, 1);
if n_samples == 0
    return;
end

[unique_files, ~, file_group] = unique(split_manifest.FilePath, "stable");
row_idx = 1 : decim_step : (1 + decim_step * (S60.nrn - 1));

fprintf("顺序加载 %d 个回波文件中的 %d 个%s到内存...\n", numel(unique_files), n_samples, string(split_label));
for f_idx = 1:numel(unique_files)
    indices = find(file_group == f_idx);
    var_name = char(split_manifest.EchoVariable(indices(1)));
    loaded = load(unique_files(f_idx), var_name);
    raw = loaded.(var_name);
    for idx = indices(:).'
        c_start = split_manifest.CStart(idx);
        sig60 = raw(row_idx, c_start : c_start + S60.nan - 1);
        if ~isequal(size(sig60), [S60.nrn, S60.nan])
            error("zc1bit:SignalSize", "回波块尺寸 [%d %d] 与 S60 [%d %d] 不一致。", ...
                size(sig60, 1), size(sig60, 2), S60.nrn, S60.nan);
        end
        signals{idx} = sig60;
    end
    clear raw loaded
end
end
