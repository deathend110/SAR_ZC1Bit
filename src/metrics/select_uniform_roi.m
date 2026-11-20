function [top, left, score] = select_uniform_roi(img_gt, window_size, stride)
%SELECT_UNIFORM_ROI 在 clean GT 上自动选定最均匀的 64x64 均匀区 ROI，供全部分支共享复用。

[height, width] = size(img_gt);
if height < window_size || width < window_size
    error("zc1bit:ROIWindow", "图像尺寸 [%d %d] 小于 ENL 窗口 %d。", height, width, window_size);
end

tops = 1:stride:(height - window_size + 1);
lefts = 1:stride:(width - window_size + 1);
count = numel(tops) * numel(lefts);
means = zeros(count, 1);
scores = inf(count, 1);
top_values = zeros(count, 1);
left_values = zeros(count, 1);

ptr = 0;
intensity = double(img_gt) .^ 2;
for row = tops
    for col = lefts
        ptr = ptr + 1;
        patch = intensity(row : row + window_size - 1, col : col + window_size - 1);
        mu = mean(patch(:));
        means(ptr) = mu;
        scores(ptr) = var(patch(:), 0) / (mu ^ 2 + eps);
        top_values(ptr) = row;
        left_values(ptr) = col;
    end
end

lower = interp_percentile(means, 20);
upper = interp_percentile(means, 80);
valid = means >= lower & means <= upper & means > eps;
if ~any(valid)
    valid = means > eps;
end
if ~any(valid)
    valid = true(size(means));
end
valid_idx = find(valid);
[score, local_idx] = min(scores(valid_idx));
best_idx = valid_idx(local_idx);
top = top_values(best_idx);
left = left_values(best_idx);
end

function value = interp_percentile(values, p)
values = sort(values(~isnan(values)));
if isempty(values)
    value = NaN;
    return;
end
pos = 1 + (numel(values) - 1) * p / 100;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    value = values(lo);
else
    value = values(lo) + (values(hi) - values(lo)) * (pos - lo);
end
end
