function context = build_fast_focus_context(S60, range_factors, azimuth_factors)
%BUILD_FAST_FOCUS_CONTEXT 预计算与回波和阈值无关的频域聚焦核，加速搜索与快筛。
% 在中心频谱上直接合并距离压缩、双向裁剪、RCMC 与方位压缩，与标准路径严格等价。

range_factors = [1; range_factors(:)];
azimuth_factors = [1; azimuth_factors(:)];
for k = 1:numel(range_factors)
    if abs(range_factors(k) - round(range_factors(k))) < 1e-12
        range_factors(k) = round(range_factors(k));
    end
end
for k = 1:numel(azimuth_factors)
    if abs(azimuth_factors(k) - round(azimuth_factors(k))) < 1e-12
        azimuth_factors(k) = round(azimuth_factors(k));
    end
end
range_factors = unique(range_factors);
azimuth_factors = unique(azimuth_factors);

context = struct();
context.nrn = S60.nrn;
context.nan = S60.nan;
context.range_factors = range_factors(:).';
context.azimuth_factors = azimuth_factors(:).';
context.range_nr_ups = zeros(1, numel(range_factors));
context.range_kernels = cell(1, numel(range_factors));
context.range_crop_indices = cell(1, numel(range_factors));
context.azimuth_na_ups = zeros(1, numel(azimuth_factors));
context.azimuth_crop_indices = cell(1, numel(azimuth_factors));

for idx = 1:numel(range_factors)
    rq = range_factors(idx);
    grid = resolve_sampling_grid(rq, 1, [S60.nrn, S60.nan], S60);
    nrn_up = grid.upsampled_size(1);
    [tnrn_up, Fs_up] = range_axis(nrn_up, grid.q_range_eff, S60);
    pulse_sample_count = fix(S60.Tp * Fs_up);
    center_idx = floor(nrn_up / 2) + 1;
    start_idx = center_idx - floor(pulse_sample_count / 2);
    end_idx = start_idx + pulse_sample_count - 1;
    assert(pulse_sample_count >= 1 && start_idx >= 1 && end_idx <= nrn_up, ...
        "距离压缩脉冲窗口超出回波范围。");
    range_filter = zeros(nrn_up, 1);
    pulse_time = tnrn_up(start_idx:end_idx);
    range_filter(start_idx:end_idx) = exp(1i * pi * S60.gama * ...
        (pulse_time(:) - 2 * S60.R0 / S60.C) .^ 2);
    context.range_nr_ups(idx) = nrn_up;
    context.range_kernels{idx} = conj(centered_fft(range_filter, 1));
    context.range_crop_indices{idx} = centered_indices(nrn_up, S60.nrn);
end

for idx = 1:numel(azimuth_factors)
    aq = azimuth_factors(idx);
    grid = resolve_sampling_grid(1, aq, [S60.nrn, S60.nan], S60);
    nan_up = grid.upsampled_size(2);
    context.azimuth_na_ups(idx) = nan_up;
    context.azimuth_crop_indices{idx} = centered_indices(nan_up, S60.nan);
end

fnan_row = S60.fnan(:).';
migration = S60.lambda ^ 2 * S60.R0 * fnan_row .^ 2 / (8 * S60.v ^ 2);
context.rcmc_kernel = exp(1i * 2 * pi * (S60.fnrn(:) * (2 * migration / S60.C)));

delta_range = S60.C / (2 * S60.Fs);
range_scale = (-S60.nrn / 2 : S60.nrn / 2 - 1).' * delta_range;
azimuth_rate = -2 * S60.v ^ 2 ./ (S60.lambda * (range_scale + S60.R0));
azimuth_filter = zeros(S60.nrn, S60.nan);
aperture_samples = fix(S60.Ta * S60.prf);
start_idx = -aperture_samples / 2 + S60.nan / 2 + 1;
end_idx = aperture_samples / 2 + S60.nan / 2;
assert(start_idx == round(start_idx) && end_idx == round(end_idx), ...
    "方位压缩窗口索引不是整数。");
aperture_idx = round(start_idx):round(end_idx);
azimuth_filter(:, aperture_idx) = exp(-1i * pi * azimuth_rate * (S60.tnan(aperture_idx) .^ 2));
context.azimuth_kernel = centered_fft(azimuth_filter, 2);

context.roi_rows = S60.nrn / 2 - S60.R_total / 2 + 1 : S60.nrn / 2 + S60.R_total / 2;
context.roi_cols = S60.nan / 2 - S60.A_num / 2 : S60.nan / 2 + S60.A_num / 2 - 1;
end

function Xf = centered_fft(X, dim)
Xf = fftshift(fft(fftshift(X, dim), [], dim), dim);
end

function indices = centered_indices(current_size, target_size)
if target_size > current_size
    error("zc1bit:CropSize", "目标频谱尺寸大于当前尺寸。");
end
center = floor(current_size / 2) + 1;
half_size = floor(target_size / 2);
if mod(target_size, 2) == 0
    indices = center - half_size : center + half_size - 1;
else
    indices = center - half_size : center + half_size;
end
end
