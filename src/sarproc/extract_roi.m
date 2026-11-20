function roi = extract_roi(IMG, S60)
%EXTRACT_ROI 从复数聚焦图像提取与历史 V4/V5 完全一致的线性幅度有效 ROI。

rows = S60.nrn / 2 - S60.R_total / 2 + 1 : S60.nrn / 2 + S60.R_total / 2;
cols = S60.nan / 2 - S60.A_num / 2 : S60.nan / 2 + S60.A_num / 2 - 1;
roi = abs(IMG(rows, cols));
end
