function S1 = quantize_1bit_with_u(S, U)
%QUANTIZE_1BIT_WITH_U 复数 1-bit 加阈值量化（非对称符号 sgn_+(0)=+1）。
% 严格采用加阈值约定 S + U，仅将严格小于 0 的位置置 -1，禁用内置 sign。

if nargin < 2 || isempty(U)
    U = complex(zeros(size(S), "like", real(S)));
end
if ~isequal(size(S), size(U))
    error("zc1bit:ThresholdSizeMismatch", "复回波与阈值矩阵尺寸必须完全一致。");
end

re = ones(size(S), "like", real(S));
im = ones(size(S), "like", real(S));
re(real(S) + real(U) < 0) = -1;
im(imag(S) + imag(U) < 0) = -1;
S1 = complex(re, im);
end
