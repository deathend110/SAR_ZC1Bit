# AGENTS.md — SAR-ZC1Bit 仓库防回归与协作规则

本文件与 [CLAUDE.md](CLAUDE.md) 保持一致，记录 MATLAB + SAR 工程防回归硬规则：

1. **唯一真相源与只读红线**：实验范围以 [experiment_protocol.md](docs/experiment_protocol.md) 为准，外部参考仓库（`G:\MATLAB-G\SAR Dataset Transformer`、`G:\MATLAB-G\SAR_Azimuth_Range_Upsample`、`G:\MATLAB-G\SAR Full PSF`）严格只读；复用函数必须独立维护在本仓 `src/` 内，禁止跨仓 `addpath`。
2. **实现授权与运行授权分离**：编写或修改代码后，默认只执行 `checkcode` 静态检查与不依赖外部大数据集的 `tests/` 轻量单元测试；未经用户明确要求，严禁启动完整阶段一快筛或阶段二参数搜索。
3. **路径解耦**：个人数据路径写入本地未跟踪文件 `inputs/local_paths.m`（模板为 `inputs/local_paths.example.m`），严禁硬编码进核心代码。
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
   - 断点恢复仅核对实际计算参数（忽略路径、时间戳、版本号、worker 数与缓存大小）。
   - 阶段一快筛与阶段二寻优严禁使用 35 个测试样本；测试集仅在阶段二锁参后用于最终评价。
7. **编码与 Git**：所有 `.m`、`.md`、`.csv` 使用 UTF-8（无 BOM）；改动后运行 `checkcode` 与 `git diff --check`。
