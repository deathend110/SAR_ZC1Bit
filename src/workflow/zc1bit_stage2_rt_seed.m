function seed = zc1bit_stage2_rt_seed(base_seed, seed_tag, sample_id)
%ZC1BIT_STAGE2_RT_SEED 统一阶段二 RT 随机阈值种子推导。
% 搜索、评估与绘图重渲染必须使用完全相同的哈希规则，才能复现同一阈值场。
% 规则：对 ParamKey 做 131 进制滚动哈希，再叠加基础种子与 SampleID。

chars = double(char(string(seed_tag)));
h = 0;
for k = 1:numel(chars)
    h = mod(h * 131 + chars(k), 1000003);
end
seed = mod(double(base_seed) + h * 1000 + double(sample_id), 2^31 - 1);
end
