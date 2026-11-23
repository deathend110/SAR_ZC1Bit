function p = zc1bit_paths(cfg)
%ZC1BIT_PATHS 统一生成结果、断点、报告与工作区临时目录结构，禁止在各脚本中手写路径拼接。

p = struct();
p.root = string(cfg.results_root);
p.protocol = fullfile(p.root, "protocol");
p.screen = fullfile(p.root, "screen");
p.search = fullfile(p.root, "search");
p.locked = fullfile(p.root, "locked");
p.evaluate = fullfile(p.root, "evaluate");
p.stats = fullfile(p.root, "stats");
p.figures = fullfile(p.root, "figures");
p.checkpoints = fullfile(p.root, "checkpoints");
p.ta_sft = fullfile(p.root, "ta_sft");
p.sft_p15 = fullfile(p.root, "sft_p15");
if isfield(cfg, "tmp_root") && strlength(string(cfg.tmp_root)) > 0
    p.tmp = string(cfg.tmp_root);
elseif isfield(cfg, "repo_root") && strlength(string(cfg.repo_root)) > 0
    p.tmp = fullfile(string(cfg.repo_root), "tmp");
else
    p.tmp = fullfile(p.root, "tmp");
end
end
