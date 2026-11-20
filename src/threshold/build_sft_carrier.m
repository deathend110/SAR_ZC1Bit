function carrier = build_sft_carrier(grid, params, S60)
%BUILD_SFT_CARRIER 构造与信号幅度无关的单位模二维 SFT 复载波。
% 统一覆盖 RSFT（fa=0）、ASFT（fr=0）与 2D-SFT，严格先将距离和方位相位相加再取 exp。

phi0 = 0;
if isfield(params, "phi0") && ~isempty(params.phi0)
    phi0 = params.phi0;
end

nr_up = grid.upsampled_size(1);
na_up = grid.upsampled_size(2);
fast_time = ((0:nr_up - 1).' - floor(nr_up / 2)) / grid.Fs_up;
slow_time = ((0:na_up - 1) - floor(na_up / 2)) / grid.PRF_up;

ba = resolve_azimuth_bandwidth(S60);
phase_range = 2 * pi * (params.fr_over_Br * S60.B) * fast_time;
phase_azimuth = 2 * pi * (params.fa_over_Ba * ba) * slow_time;
carrier = exp(1i * (phase_range + phase_azimuth + phi0));
end
