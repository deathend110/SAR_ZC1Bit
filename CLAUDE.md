# CLAUDE.md — SAR-ZC1Bit 仓库协作与工程约束

## 1. 仓库目标与唯一真相源

- 本仓库为 **`SAR-ZC1Bit`**，研究在**保持采集与传输端严格 1-bit（±1）不变**的前提下，于成像前施加伪 1.5-bit 过零置零重构（`{+1, 0, -1}`），在 ZT / RT / RSFT / ASFT 四类阈值与 1–9 倍单向上采样（距离单向与方位单向）下提升 SAR 成像的 PSNR、SSIM 与 ENL；双向 BARU 与 2D-SFT 保留底层接口作为后续可选扩展。
- 实验主题、规则定义、数据划分与指标口径以 `docs/EXPERIMENT_PROTOCOL.md` 为唯一真相源，不得擅自偏离、删减核心规则或扩展无关功能（例如真 1.5-bit 硬件仿真、深度学习去噪等）。
- 仓库架构、工作流分层与代码风格对齐 `G:\MATLAB-G\SAR-Dataset-Transformer`。

## 2. 本机参考目录与只读红线

本机存在以下参考资源，**全部严格只读，严禁修改、删除或写入任何文件**：
- 架构参考：`G:\MATLAB-G\SAR-Dataset-Transformer`
- 算法与历史参数参考：`G:\MATLAB-G\SAR_Azimuth_Range_Upsample`
- 阈值与噪声知识库：`G:\Zhe_KnowledgeBase`

**依赖隔离原则**：
- 需要复用的基础函数（`Range_Compress`、`RCMC`、`SAR_Imaging`、`gaussian`、`normalize_image` 等）必须拷贝到本仓库 `src/` 内独立维护，禁止跨仓库 `addpath` 耦合。
- 本机原始回波数据与参数路径统一写入本地未跟踪文件 `inputs/local_paths.m`（仓库仅提交模板 `inputs/local_paths.example.m`），严禁在核心代码中硬编码个人绝对路径。
- **工作区产物约束（禁止写 C 盘）**：所有实验产物（默认 `results_zc1bit/`）与测试/冒烟/临时产物（`tmp/tests/`、`tmp/smoke/`、`tmp/scratch/`，已由 `.gitignore` 忽略）均必须保存在本仓库工作区内；严禁写入系统 C 盘临时目录（禁止使用无参 `tempname`、`tempdir` 或 `$env:TEMP`）。

## 3. 授权边界（实现授权与运行授权分离）

1. **先分析/确认，再改代码**：用户未明确要求“写代码/修改代码/实现功能/修复问题”前，默认不主动修改代码；若任务主要是讨论、分析或设计，优先输出结论。
2. **只负责修改代码，不启动 MATLAB**：代理仅负责编写与修改代码/文档及非 MATLAB 静态校验（UTF-8 无 BOM、`git diff --check`、文本规则扫描），**严禁命令行启动 MATLAB（禁用 `matlab -batch` 等任何 MATLAB 启动调用）**。`checkcode`、`tests/` 单元测试、冒烟测试与完整实验均由用户在本地 MATLAB 中自行运行。
3. **Windows + MATLAB 进程与运行约束**：
   - 严禁启动新的 MATLAB 进程，也严禁干扰或终止用户已打开的 MATLAB GUI 会话。
   - 代码保持单核单会话设计（不使用 MATLAB 并行池 `parpool`），优先利用内存缓存（推荐 `4–6 GB` 预算）避免机械硬盘重复读取与高频小文件写入。

## 4. 编码、注释与 Git 规范

- **编码**：所有 `.m`、`.md`、`.csv` 文件统一使用 **UTF-8（无 BOM）** 读写。
- **中文注释**：说明关键步骤、物理约定、维度方向与非直观设计意图；不写无意义的逐行注释。
- **Git 习惯**：
  - 默认直接在主工作区开发，不创建额外 `git worktree`；改动前后检查 `git status`。
  - 不覆盖用户未提交改动；严禁误删 `results_zc1bit/`、checkpoint 或已生成图表。
  - 每次完成独立功能或文档改动后，检查 `git diff --check`，并用清晰中文提交一次 `git commit`。

## 5. MATLAB 与 SAR 防回归硬规则

### 5.1 SAR 信号与量化约定
1. **维度锁定**：第 1 维（行）= 距离向（快时间），第 2 维（列）= 方位向（慢时间）。
2. **量化约定**：统一加阈值 $S^\uparrow + U$；非对称符号 $\operatorname{sgn}_+(0)=+1$（初始化为 `+1`，仅将 `< 0` 置 `-1`，禁用内置 `sign`）。
3. **伪 1.5-bit 位置**：仅作用于标准 1-bit 量化输出之后、距离压缩之前，对实部和虚部各自独立执行，不改变矩阵尺寸。
4. **二维阈值构造**：必须先将距离与方位相位相加，再统一取复指数；禁止把两个复阈值直接相加。
5. **同协议对比**：所有分支（1-bit 与伪 1.5-bit、各上采样分配）必须共用同一 GT、同一 ROI 与统一归一化协议。
6. **随机性**：RT 阈值统一固定 `rng(seed, "twister")`，按组名固定 RNG 编号。

### 5.2 MATLAB 常见坑防回归
1. **Table 构造**：累加器初始化为 `table()`（禁用 `[]` 纵向拼接）；列选择用 `T(:, ["A","B"])` 或 `T(:, {'A','B'})`，禁用 `T(:, {"A","B"})`；常量列用 `repmat` 展开到相同行数。
2. **数组方向与统计检验**：输入 `signrank`（Wilcoxon）或配对差值前，先统一 `(:)` 列向量并核对配对顺序与返回值接收方式。
3. **字符串与空结果**：cell 中 string 先统一转 string 数组再用 `==`；筛选绘图前先检查 `isempty` / `any(mask)`。
4. **小数倍率与文件键**：若涉及小数倍率，使用 `round(q * N)` 与近整数吸附（`abs(p2 - round(p2)) < 1e-12`）；checkpoint 与输出文件名统一使用稳定 ASCII `FileKey`。
5. **落盘前建目录**：保存 CSV/MAT/PNG/metadata 前先确保目标目录存在；写出前检查 `height`、变量名与非有限值。

### 5.3 Checkpoint 与防数据泄漏
1. 阶段一按 `ThresholdKey`（35 样本整体）在扁平目录原子保存 checkpoint，`< 1.5 GB` 载荷默认使用 `-v7` 格式（避免逐场景/逐样本高频小文件写入造成机械盘寻道）；后续阶段二搜索同样保持粗粒度原子落盘，减少小文件写盘频率。
2. 断点恢复**仅核对实际计算参数**（倍率、阈值类型与参数、置零规则 Key、判定/去 DC 开关、样本划分、归一化口径），禁止校验完整路径、时间戳、`runtime.cache_gb`、`preload_signals`、`use_fast_path` 或版本号导致无意义重算。
3. 严禁从未完成的搜索 cube 中提前选优；严禁测试集参与任何参数搜索或规则筛选。
