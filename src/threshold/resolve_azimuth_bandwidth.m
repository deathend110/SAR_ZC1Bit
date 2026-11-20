function bandwidth = resolve_azimuth_bandwidth(S60)
%RESOLVE_AZIMUTH_BANDWIDTH 按 Ba、Bd、2*v/Da 优先级解析方位多普勒带宽。

if isfield(S60, "Ba")
    bandwidth = S60.Ba;
elseif isfield(S60, "Bd")
    bandwidth = S60.Bd;
elseif isfield(S60, "Da")
    bandwidth = 2 * S60.v / S60.Da;
else
    error("zc1bit:AzimuthBandwidth", "FS60 参数必须包含 Ba、Bd 或 Da 字段。");
end

if ~isscalar(bandwidth) || ~isfinite(bandwidth) || bandwidth <= 0
    error("zc1bit:AzimuthBandwidth", "解析得到的方位带宽 Ba 必须为正有限标量。");
end
end
