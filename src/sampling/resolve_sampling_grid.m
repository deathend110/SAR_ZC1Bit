function grid = resolve_sampling_grid(range_q, azimuth_q, input_size, S60)
%RESOLVE_SAMPLING_GRID 将距离与方位名义倍率映射为实际整数上采样网格。
% 维度约定：dim 1 为距离向（快时间），dim 2 为方位向（慢时间）。
% 对近整数浮点误差（例如 sqrt(3)^2）执行 1e-12 阈值吸附。

range_q = snap_near_integer(range_q);
azimuth_q = snap_near_integer(azimuth_q);

if ~isscalar(range_q) || ~isfinite(range_q) || range_q < 1
    error("zc1bit:SamplingFactor", "距离上采样倍率必须是不小于 1 的有限实数。");
end
if ~isscalar(azimuth_q) || ~isfinite(azimuth_q) || azimuth_q < 1
    error("zc1bit:SamplingFactor", "方位上采样倍率必须是不小于 1 的有限实数。");
end
if ~isnumeric(input_size) || numel(input_size) ~= 2 || any(input_size < 1) || any(mod(input_size, 1) ~= 0)
    error("zc1bit:InputSize", "input_size 必须为正整数二元向量 [Nr, Na]。");
end

nr_up = round(range_q * input_size(1));
na_up = round(azimuth_q * input_size(2));
q_range_eff = nr_up / input_size(1);
q_azimuth_eff = na_up / input_size(2);
q_total_eff = snap_near_integer(q_range_eff * q_azimuth_eff);

grid = struct( ...
    "q_range", range_q, ...
    "q_azimuth", azimuth_q, ...
    "q_total", snap_near_integer(range_q * azimuth_q), ...
    "q_range_eff", q_range_eff, ...
    "q_azimuth_eff", q_azimuth_eff, ...
    "q_total_eff", q_total_eff, ...
    "input_size", double(input_size(:).'), ...
    "upsampled_size", [nr_up, na_up], ...
    "Fs_up", q_range_eff * S60.Fs, ...
    "PRF_up", q_azimuth_eff * S60.prf);
end

function value = snap_near_integer(value)
if isnumeric(value) && isscalar(value) && isfinite(value) && abs(value - round(value)) < 1e-12
    value = round(value);
end
end
