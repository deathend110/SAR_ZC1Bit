function [img_gt, complex_roi, RC_gt] = generate_gt_image(signal60, S60)
%GENERATE_GT_IMAGE 对未量化 60MHz 复回波执行标准全精度成像与统一归一化。

RC_gt = Range_Compress(signal60, S60.fc, S60.tnrn, S60.gama, ...
    S60.R0, S60.C, S60.Fs, S60.Tp);
RCMC_gt = RCMC(RC_gt, S60.lambda, S60.fnrn, S60.fnan, ...
    S60.R0, S60.C, S60.v);
IMG_gt = SAR_Imaging(RCMC_gt, S60.lambda, S60.Fs, ...
    S60.R0, S60.C, S60.v, S60.tnan, S60.Ta, S60.prf);

rows = S60.nrn / 2 - S60.R_total / 2 + 1 : S60.nrn / 2 + S60.R_total / 2;
cols = S60.nan / 2 - S60.A_num / 2 : S60.nan / 2 + S60.A_num / 2 - 1;
complex_roi = IMG_gt(rows, cols);
img_gt = normalize_image(abs(complex_roi));
end
