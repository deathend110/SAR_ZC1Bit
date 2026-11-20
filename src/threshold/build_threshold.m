function [threshold, meta] = build_threshold(signal_up, grid, threshold_spec, S60, options)
%BUILD_THRESHOLD 统一构造 ZT / RT / RSFT / ASFT / 2D-SFT 五类复阈值矩阵。

if nargin < 5 || isempty(options)
    options = struct();
end
nyquist_margin = 0.98;
if isfield(options, "nyquist_margin") && ~isempty(options.nyquist_margin)
    nyquist_margin = options.nyquist_margin;
end
enforce_nyquist = true;
if isfield(options, "enforce_nyquist") && ~isempty(options.enforce_nyquist)
    enforce_nyquist = logical(options.enforce_nyquist);
end
rand_stream = [];
if isfield(options, "rand_stream") && ~isempty(options.rand_stream)
    rand_stream = options.rand_stream;
end
carrier = [];
if isfield(options, "carrier") && ~isempty(options.carrier)
    carrier = options.carrier;
end

ttype = string(threshold_spec.ThresholdType);
sigma_hat = estimate_sigma_hat(signal_up);
meta = struct("ThresholdType", ttype, "SigmaHat", sigma_hat, "Amplitude", 0);

switch ttype
    case "ZT"
        threshold = complex(zeros(size(signal_up), "like", real(signal_up)));
    case "RT"
        threshold = build_splitrt_threshold(signal_up, threshold_spec.As, rand_stream);
        meta.Amplitude = threshold_spec.As * sigma_hat;
    case {"RSFT", "ASFT", "2D-SFT"}
        if ttype == "RSFT" && abs(threshold_spec.fa_over_Ba) > 0
            error("zc1bit:RSFTAzimuthFreq", "RSFT 的 fa_over_Ba 必须严格为 0。");
        end
        if ttype == "ASFT" && abs(threshold_spec.fr_over_Br) > 0
            error("zc1bit:ASFTRangeFreq", "ASFT 的 fr_over_Br 必须严格为 0。");
        end
        if enforce_nyquist
            check_sft_nyquist(grid, threshold_spec, S60, nyquist_margin);
        end
        [threshold, ~, amplitude] = build_sft_threshold(signal_up, grid, threshold_spec, S60, carrier);
        meta.Amplitude = amplitude;
    otherwise
        error("zc1bit:ThresholdType", "未知阈值类型：%s", ttype);
end
end
