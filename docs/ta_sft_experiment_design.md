# SAR-ZC1Bit 阈值感知 SFT 1-bit 接收端重构（TA-SFT）实验设计书

**版本**：v1.0  
**所属项目**：`SAR-ZC1Bit`  
**文档定位**：阐述基于已知正弦频率阈值（RSFT / ASFT）先验的 1-bit 接收端重构专题实验动机、数学定义、防信息泄漏协议、鲁棒性消融与实现验收标准。

---

## 1. 研究背景与问题提出

### 1.1 前期盲过零实验暴露的核心瓶颈
在阶段一与阶段二主实验（见 [experiment_report.md](experiment_report.md)）中，基于纯符号跳变的盲双沿置零规则 **`R1_par_raw`** 在随机阈值（RT）高倍区间（距离 $Q=4\text{–}9$、方位 $Q=3\text{–}9$）实现了稳定的 PSNR/SSIM 显著提升（$Q=9$ 时 $\Delta\text{SSIM} = +0.0134\text{–}+0.0146$，$\Delta\text{PSNR} = +0.50\text{–}+0.57\text{ dB}$）。然而，在成像质量天花板远高于 RT 的正弦频率阈值（RSFT / ASFT）下，盲过零策略遭遇明显瓶颈：
- **方位 ASFT**：$Q=1\text{–}4$ 明显弱于基线 `B0_par_raw`；虽然在 $Q=7\text{–}9$ 统计显著反超且于 `ab_R1A9_azimuth_asft` 取得全场最高测试均值（SSIM `0.9465` / PSNR `33.01 dB`），但增益幅度仅 `+0.0006` / `+0.13 dB`。
- **距离 RSFT**：$Q=1\text{–}8$ 全线弱于 `B0_par_raw`（$Q=4$ 落差达 `-0.0340` SSIM / `-1.37 dB`），至 $Q=9$ 才仅在 PSNR 上微弱领先（`+0.06 dB`）。

### 1.2 物理根因：$S + U = 0$ 而非 $S = 0$
在 SFT 体制下，叠加于上采样复回波 $S^\uparrow$ 的阈值具有高幅、高频的解析形式（最优 $\mathrm{STR} = -4 \sim -5\text{ dB}$ 对应包络幅值 $A_{\mathrm{SFT}} \approx 1.58\text{–}1.78\,\hat\sigma$）：

$$U[n] = A_{\mathrm{SFT}}\,e^{j(\phi_r[n_r] + \phi_a[n_a] + \phi_0)}$$

经过加阈值 1-bit 量化 $b = \operatorname{sgn}_+(S^\uparrow + U)$ 后，相邻采样点发生过零（$b[k] \neq b[k+1]$）并不代表真实 SAR 回波接近零，而是满足：

$$S_I[n] + U_I[n] \approx 0 \quad \Longleftrightarrow \quad S_I[n] \approx -U_I[n] \qquad \text{（虚部 Q 路同理）}$$

* **在 RT 下**：阈值零均值独立随机，强信号难以被翻转，过零统计集中发生在小幅区域 $|S| < |U|$，故置 `0` 有效。
* **在 SFT 下**：高幅正弦阈值主导了 52%–71% 的邻点翻转；盲目将所有过零点置 `0`，实质上是把大量由正弦载波切过的高幅度真实回波强行抹零，破坏了 SFT 的带外谐波搬移结构。

因此，如果合理利用接收端**已经已知**的 SFT 阈值波形先验，将“盲过零置零”升级为**阈值感知（Threshold-Aware）接收端重构**，便有望克服盲置零的根本缺陷。

---

## 2. 工程假设、系统接口与防信息泄漏红线

### 2.1 硬件与工程合理性
1. **严格 1-bit 前端采集与传输**：星上/采集端依然执行最简单的 $S^\uparrow + U \to \operatorname{sgn}_+$ 比较器，输出与下行数据流保持 **100% 纯 1-bit（$\pm 1$）**，零额外 ADC 硬件复杂度、零额外存储与通信比特开销。
2. **已知 SFT 载波的合法先验**：SFT 参考阈值是由雷达系统主振/DDS 与 ADC 采样时钟相干锁相生成的确定性单频波。频率 $(f_r/B_r, f_a/B_a)$、物理电平与标定初相 $\phi_0$ 属于雷达已知工作模式参数，随常规帧头遥测参数即可在地面成像端精确重构单位模解析载波 $u[n] = \exp(j\theta[n])$。

### 2.2 防信息泄漏接口红线
为确保实验结论绝对可信、不存在任何“偷看未量化回波”的漏洞，核心函数 [apply_ta_sft_reconstruct.m](../src/quantization/apply_ta_sft_reconstruct.m) 在架构层锁死以下红线：
* **禁止传入 `signal_up` 与 $\hat\sigma$**：函数签名固定为 `[channel_out, diag] = apply_ta_sft_reconstruct(channel_1bit, carrier, method_key, options)`。
* **严格输入检验**：`channel_1bit` 实部与虚部必须严格属于 $\{+1, -1\}$；`carrier` 仅为单位模复载波（$|u| = 1$，实虚部取值严格在 $[-1, 1]$，绝不包含任何未量化幅度尺度信息）。
* **样本严格隔离**：全部超参仅在 35 奇数校准样本上寻优锁定；35 偶数测试样本严禁参与任何参数选择。

---

## 3. 重构方案数学定义（7 方法对照）

记单个通道（I 路或 Q 路，各自独立处理）的 1-bit 序列为 $b[n] \in \{+1, -1\}$，对应的接收端归一化载波分量为 $u[n] \in [-1, 1]$（I 路为 $\mathrm{Re}(u)$，Q 路为 $\mathrm{Im}(u)$）。沿当前主轴（RSFT 沿 dim 1，ASFT 沿 dim 2）定义相邻过零边指示 $e[k] = \mathbb{I}(b[k] \neq b[k+1])$ 与边中点归一化载波 $\bar{u}[k] = \frac{1}{2}(u[k] + u[k+1])$。

| 序号 | 方法 Key | 类别定位 | RC 前输出数据集合 | 核心数学规则 | 超参 | 退化边界 |
| :---: | :--- | :--- | :---: | :--- | :---: | :--- |
| 1 | `B0_par_raw` | 原始 1-bit 基线 | $\{+1, -1\}$ | 保持 $b[n]$ 不变 | 无 | — |
| 2 | `R1_par_raw` | 盲过零伪 1.5-bit 对照 | $\{+1, 0, -1\}$ | 全部 $e[k]=1$ 边两端并行置 `0` | 无 | — |
| 3 | `TA1_gate_zc` | **方案 1**：阈值门控选择性置零 | **严格三电平** $\{+1, 0, -1\}$ | 仅对 $e[k]=1$ 且 $|\bar{u}[k]| \le \tau$ 的边两端并行置 `0`，强阈值处过零保留 $\pm 1$ | $\tau \ge 0$ | $\tau=0 \Rightarrow B_0$<br>$\tau=1 \Rightarrow R_1$ |
| 4 | `TA2_anchor_zc` | **方案 2A**：过零锚点连续替换 | 非过零 $\pm 1$ + 过零浮点 | 过零边赋估值 $v[k] = -\alpha\bar{u}[k]$；被单边命中的节点取 $v$，被双边同时命中的节点取均值 $\frac{1}{2}(v[n-1]+v[n])$ | $\alpha \ge 0$ | $\alpha=0 \Rightarrow R_1$ |
| 5 | `TA2_anchor_q5` | **方案 2B**：过零锚点 5 电平量化 | **严格 5 电平** $\{-1, -0.5, 0, +0.5, +1\}$ | **先**按方案 2A 计算过零节点均值 $y_{\mathrm{node}}$，**再**对节点均值执行对称 5 电平量化：$\operatorname{round}(2\operatorname{clip}(y_{\mathrm{node}}, -1, 1))/2$（非过零点保持 $\pm 1$） | $\alpha \ge 0$ | $\alpha=0 \Rightarrow R_1$ |
| 6 | `TA3_mmse_fixed` | **方案 3A**：固定标称 $\gamma$ 逐点 MMSE | 全点连续浮点 | $\hat{y}[n] = b[n]\,\frac{\phi(z_n)}{\Phi(z_n)}$，其中 $z_n = b[n]\cdot \gamma u[n]$ | $\gamma \ge 0$ | $\gamma=0 \Rightarrow B_0$ |
| 7 | `TA3_mmse_est` | **方案 3B**：1-bit 自适应盲估计 MMSE | 全点连续浮点 | 仅从 1-bit 内积相关 $\hat\rho = \frac{1}{2}(\overline{b_I u_I} + \overline{b_Q u_Q})$ 查单调函数 $g(\gamma)$ 反演 $\hat\gamma = g^{-1}(\hat\rho)$，再执行同方案 3A 逐点 MMSE | **0 超参** | 无未量化先验 |

### 关键数值实现保障
1. **先平均再量化（`TA2_anchor_q5`）**：若先对每条边独立量化至 5 电平再求双边均值，会产生 $\pm 0.25$、$\pm 0.75$ 等中间电平破坏 5 电平集合约束；本实现严格先求节点均值再统一做对称 5 电平量化，确保全矩阵所有元素严格属于 $\{-1, -0.5, 0, +0.5, +1\}$。
2. **全实轴无溢出逆 Mills 比（`TA3_mmse_*`）**：
   $$\frac{\phi(z)}{\Phi(z)} = \begin{cases} \sqrt{\dfrac{2}{\pi}}\,\dfrac{e^{-z^2/2}}{\operatorname{erfc}(-z/\sqrt{2})}, & z \ge 0 \\ \sqrt{\dfrac{2}{\pi}}\,\dfrac{1}{\operatorname{erfcx}(-z/\sqrt{2})}, & z < 0 \end{cases}$$
   利用 scaled complementary error function $\operatorname{erfcx}(w) = e^{w^2}\operatorname{erfc}(w)$ 彻底消除 $z \ll -1$ 时的指数下溢与除零风险。
3. **纯 1-bit $\hat\gamma$ 无偏估计原理（`TA3_mmse_est`）**：对于单路零均值高斯回波 $S_I \sim \mathcal{N}(0, \sigma^2)$ 与相干载波 $u_I = \cos\theta$、$U_I = \gamma\sigma\cos\theta$，条件期望 $\mathbb{E}[b_I \mid \theta] = \operatorname{erf}(\gamma\cos\theta/\sqrt{2})$，故交叉相关满足严格单调确定性映射：
   $$g(\gamma) = \mathbb{E}[b_I u_I] = \frac{1}{2\pi}\int_0^{2\pi} \operatorname{erf}\!\left(\frac{\gamma\cos\theta}{\sqrt{2}}\right)\cos\theta\,d\theta \in \left[0, \frac{2}{\pi}\right)$$
   由实虚两路 1-bit 与单位模载波内积均值 $\hat\rho$ 经 PCHIP 分段插值反演即可无偏恢复有效阈信比 $\hat\gamma$，**无需任何先验 $\hat\sigma$ 也不依赖离线超参网格**。

---

## 4. 实验协议与轨道设置

### 4.1 轨道 A：固定阶段二 `B0_par_raw` 最优 SFT 参数（主实验，当前默认）
* **单元范围（共 18 单元）**：
  * 距离单向 RSFT：`rb_R1A1_range_rsft` ~ `rb_R9A1_range_rsft`（9 单元）；
  * 方位单向 ASFT：`ab_R1A1_azimuth_asft` ~ `ab_R1A9_azimuth_asft`（9 单元）。
* **前端 SFT 参数锁定**：直接读取阶段二 [locked_params.csv](../results_zc1bit/locked/locked_params.csv) 中 `VariantKey == "B0_par_raw"` 在对应单元的锁定 $(S_{\mathrm{dB}}, f_r/B_r, f_a/B_a, \phi_0)$，保证前端物理调制完全不变，在 `B0` 最优主场参数下正面检验接收端重构的净收益。
* **35 校准样本超参寻优网格（每单元 44 候选配置）**：
  * `TA1_gate_zc`：$\tau \in \{0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0\}$（11 点）；
  * `TA2_anchor_zc` 与 `TA2_anchor_q5`：$\alpha \in \{0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0\}$（各 10 点）；
  * `TA3_mmse_fixed`：$\gamma \in \{0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.2, 2.6\}$（10 点）；
  * `B0_par_raw`、`R1_par_raw`、`TA3_mmse_est`：无网格超参（各 1 点）。
  * 寻优准则：校准集 35 样本平均 SSIM 最大（平均 PSNR 破局、`AbsENLDiff` 次破局）。

### 4.2 工程消融：接收端初相失配鲁棒性测试
为量化评估非理想硬件锁相（存在固定全局相偏 $\Delta\phi_0$）时的退化边界，在 4 个代表单元（`rb_R4A1_range_rsft`、`rb_R9A1_range_rsft`、`ab_R1A4_azimuth_asft`、`ab_R1A9_azimuth_asft`）上，保持前端生成 `channel_1bit` 使用真实载波 $u_{\mathrm{true}}$，而接收端注入相偏载波 $u_{\mathrm{rx}} = u_{\mathrm{true}}\,e^{j\Delta\phi_0}$：
$$\Delta\phi_0 \in \{0^\circ, 5^\circ, 10^\circ, 15^\circ, 30^\circ, 45^\circ, 90^\circ\}$$
在 35 测试样本上对比 5 个 `TA*` 方法与相偏免疫基线（`B0_par_raw`、`R1_par_raw`）的性能衰减曲线与配对显著性。

### 4.3 后续可选扩展（轨道 B）
若轨道 A 显示某方案（如 `TA2` 或 `TA3`）与更深/不同频率的 SFT 调制存在协同空间，后续可针对代表倍率（$Q\in\{4,9\}$）扩展 SFT 前端参数 $(S, f)$ 与重构超参的联合二维搜索，不纳入本次主实验默认耗时。

---

## 5. 评价指标、统计检验与隔离产出架构

* **全量指标**：复用 Stage 2 同款 GT 缓存与固定 $64\times 64$ 均匀区 ROI（`filed` + `port` 10 测试样本），统计 `PSNR`、`SSIM`、`Entropy`、`ENL`、`AbsENLDiff`、`ZeroRatio`、`ModifiedRatio`（重构后数值发生改变的样本比例）、`EstGamma`、`AbsMeanBias` 及代表性倍率（`R1A1`/`R4A1`/`R9A1`/`R1A4`/`R1A9`）频谱带外泄漏比（`OffSupportRatio`, `RangeOffRatio`, `AzimuthOffRatio`）。
* **双对照显著性检验**：对全部 35 测试样本配对差值执行双侧 Wilcoxon 符号秩检验（`signrank`），同时输出 **`vs_B0`**（6 非基线方法 vs `B0_par_raw`）与 **`vs_R1`**（5 `TA*` 方法 vs `R1_par_raw`）。
* **产物目录严格隔离**：所有断点与产出存于 `results_zc1bit/checkpoints/ta_sft/` 与 `results_zc1bit/ta_sft/`（含 `search/`、`evaluate/`、`stats/`、`robustness/`、`figures/`），绝不改写 Stage 1 / Stage 2 既有目录。
* **验收标准**：
  1. 静态语法与防回归检查（`checkcode`、UTF-8 无 BOM、`git diff --check`）；
  2. 单元测试 [TestTaSftReconstruction.m](../tests/TestTaSftReconstruction.m) 与原有全量 `tests/` 全部通过；
  3. 主实验输出维度完整性校验（锁参 18×7 = 126 行、测试汇总 126 行、Wilcoxon 198 行、相偏消融 4×7×7 = 196 行、全核心指标无 NaN/Inf）。
