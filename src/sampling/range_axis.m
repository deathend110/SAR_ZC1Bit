function [tnrn_up, Fs_up] = range_axis(nrn_up, range_q_eff, S60)
%RANGE_AXIS 构造距离上采样后的快时间轴（列向量）与等效采样率。

Fs_up = range_q_eff * S60.Fs;
if abs(range_q_eff - 1) < 1e-12 && nrn_up == S60.nrn
    tnrn_up = S60.tnrn;
    return;
end
start_time = 2 * S60.R0 / S60.C - floor(nrn_up / 2) / Fs_up;
tnrn_up = start_time + (0:nrn_up - 1).' / Fs_up;
end
