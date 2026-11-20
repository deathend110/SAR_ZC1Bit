function p = zc1bit_paths(cfg)
%ZC1BIT_PATHS 统一生成结果、断点与报告目录结构，禁止在各脚本中手写路径拼接。

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
end
