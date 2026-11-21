function zc1bit_save_image_montage(img_list, title_list, out_path, style, mode)
%ZC1BIT_SAVE_IMAGE_MONTAGE 保存 GT + 多方法 SAR 图像/误差图拼图。
% img_list{1} 必须是 GT；其余为待比较方法。mode: "image" | "error"。

if nargin < 5 || isempty(mode)
    mode = "image";
end
mode = string(mode);
img_list = img_list(:);
title_list = string(title_list(:));
if numel(img_list) ~= numel(title_list)
    error("zc1bit:PlotMontageSize", "图像数量 (%d) 与标题数量 (%d) 不一致。", ...
        numel(img_list), numel(title_list));
end
zc1bit_ensure_dir(fileparts(string(out_path)));

n = numel(img_list);
n_cols = 4;
n_rows = ceil(n / n_cols);
fig = figure("Visible", "off", "Color", "w", ...
    "Position", [100, 100, 1500, 380 * n_rows]);
tl = tiledlayout(fig, n_rows, n_cols, "TileSpacing", "compact", "Padding", "compact");
gt_img = img_list{1};

if mode == "error"
    err_max = 0;
    err_cells = cell(n, 1);
    for k = 1:n
        err_cells{k} = abs(double(img_list{k}) - double(gt_img));
        err_max = max(err_max, max(err_cells{k}(:)));
    end
    if ~isfinite(err_max) || err_max <= 0
        err_max = 1e-12;
    end
end

last_ax = [];
for k = 1:n
    ax = nexttile(tl);
    last_ax = ax;
    if mode == "error"
        imagesc(ax, err_cells{k});
        colormap(ax, hot(256));
        caxis(ax, [0, err_max]);
    else
        imagesc(ax, double(img_list{k}));
        colormap(ax, gray(256));
        caxis(ax, [0, 1]);
    end
    axis(ax, "image");
    axis(ax, "off");
    title(ax, title_list(k), "Interpreter", "none", "FontSize", 9);
end

if ~isempty(last_ax)
    cb = colorbar(last_ax);
    if mode == "error"
        cb.Label.String = "|error|";
    else
        cb.Label.String = "normalized amplitude";
    end
    if isprop(cb, 'Layout')
        cb.Layout.Tile = "east";
    end
end
exportgraphics(fig, char(out_path), "Resolution", style.export_resolution);
close(fig);
end
