# AGENTS.md — SAR-ZC1Bit 仓库防回归与协作规则

本文件与 [CLAUDE.md](CLAUDE.md) 保持一致，记录 MATLAB + SAR 工程防回归硬规则：

1. **唯一真相源与只读红线**：实验范围以 [experiment_protocol.md](docs/experiment_protocol.md) 为准，外部参考仓库（`G:\MATLAB-G\SAR Dataset Transformer`、`G:\MATLAB-G\SAR_Azimuth_Range_Upsample`、`G:\MATLAB-G\SAR Full PSF`）严格只读；复用函数必须独立维护在本仓 `src/` 内，禁止跨仓 `addpath`。
2. **只修改代码，不启动 MATLAB**：代理仅负责修改代码、文档与非 MATLAB 静态检查（如编码、无 BOM、`git diff --check`），严禁启动 MATLAB 进程（包括 `matlab -batch`）；`checkcode`、`tests/` 单元测试及各阶段实验均由用户在本地 MATLAB 中自行运行。
3. **路径解耦与工作区产物约束**：个人数据路径写入本地未跟踪文件 `inputs/local_paths.m`（模板为 `inputs/local_paths.example.m`），严禁硬编码进核心代码；全部实验与测试产物均须留在仓内工作区（正式产物放 `results_zc1bit/`，测试/冒烟/临时文件放 `tmp/tests/`、`tmp/smoke/`、`tmp/scratch/` 并由 Git 忽略），严禁使用系统 C 盘临时目录（禁用无参 `tempname`/`tempdir` 或 `$env:TEMP`）。
4. **SAR 信号与量化约定**：
   - 第 1 维 = 距离向（快时间），第 2 维 = 方位向（慢时间）；双轴函数参数顺序统一为 `(range_q, azimuth_q)`。
   - 1-bit 量化采用加阈值 $S^\uparrow + U$ 与非对称 $\operatorname{sgn}_+(0)=+1$（先 `ones` 初始化，再将 `< 0` 置 `-1`，禁用内置 `sign`）。
   - 伪 1.5-bit 仅作用于 1-bit 量化后、距离压缩前，实部与虚部独立执行且不改变矩阵尺寸。
   - 二维阈值（SplitRT、2D-SFT）必须先相位相加，再统一取复指数。
5. **MATLAB 防回归**：
   - 表格累加器用 `table()`，禁用 `[]`；列索引使用 `T(:, ["A","B"])` 或 `T(:, {'A','B'})`，禁用 `{"A","B"}`；常量列用 `repmat` 展开。
   - `signrank` 输入前强制 `(:)` 列向量并接收第一返回值 `p`。
   - string cell 比较前先转 `string(...)`；绘图/赋值前检查 `isempty` / `any(mask)`。
   - 落盘前显式建目录，写出 CSV/MAT 前校验 `height`、变量名与非有限值。
6. **Checkpoint 与防泄漏**：
   - 保持单核执行（不启用 `parpool`），阶段一按 `ThresholdKey`（35 样本整体）扁平 `-v7` 原子落盘；断点恢复仅核对实际计算参数（忽略路径、时间戳、版本号、`runtime.cache_gb` 与 `preload_signals`）。
   - 阶段一快筛与阶段二寻优严禁使用 35 个测试样本；测试集仅在阶段二锁参后用于最终评价。
7. **编码与 Git**：所有 `.m`、`.md`、`.csv` 使用 UTF-8（无 BOM）；改动后运行 UTF-8 无 BOM 检查与 `git diff --check`（不启动 MATLAB）。
