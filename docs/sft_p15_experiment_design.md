# SAR-ZC1Bit 严格伪 1.5-bit SFT 重构（SFT-P15）实验设计书

- **版本**：v1.0
- **所属项目**：`SAR-ZC1Bit`
- **文档定位**：阐述在保持“纯 1-bit（$\pm 1$）采集/传输 + 距离压缩前严格三电平伪 1.5-bit（$\{+1, 0, -1\}$）”硬性契约的前提下，针对正弦频率阈值（RSFT / ASFT）物理特性设计的四套 SFT-P15 重构方案之研究动机、数学推导、防信息泄漏机制、载波抑制与相偏鲁棒性实验设计。

---

## 1. 研究动机与 SFT 1-bit PWM 机理分析

### 1.1 前期结果与硬约束下的新问题
在 [experiment_report.md](experiment_report.md) 总结的前期阶段一/二主实验中，基于纯邻点跳变的盲过零置零规则 **`R1_par_raw`**（输出严格属于 $\{+1, 0, -1\}$）在随机阈值（RT）中/高倍区间取得了显著收益，但在正弦频率阈值（RSFT / ASFT）上受到显著制约：
- **距离 `RSFT`**：在 $Q=1\text{–}8$ 全线弱于 `B0_par_raw`（$Q=1$ 落差达 `-0.4535` SSIM / `-6.21 dB`，$Q=4$ 落差 `-0.0340` SSIM / `-1.37 dB`），直到 $Q=9$ 才在 PSNR 上略微领先 `+0.06 dB`（SSIM 仍持平 `-0.0002`）。
- **方位 `ASFT`**：在 $Q=1\text{–}4$ 明显弱于 `B0_par_raw`；在 $Q=9$（`ab_R1A9_azimuth_asft`）虽然取得全场最高均值（SSIM `0.9465` / PSNR `33.01 dB`），但相比 `B0_par_raw` 增益仅 `+0.0006` / `+0.13 dB`。

同时，并行课题 [ta_sft_experiment_design.md](ta_sft_experiment_design.md) 探索了接收端 5 电平量化与连续浮点 MMSE 重构。本专题回到 `SAR-ZC1Bit` 最核心的表示约束——**送入距离压缩前的信号必须严格为三电平伪 1.5-bit $\{+1, 0, -1\}$（仅对原始 1-bit 序列做选择性置零，不引入多电平幅度拟合与连续浮点）**——专门回答一个核心科学问题：
> **在严格 $\{+1, 0, -1\}$ 三电平约束下，能否针对 SFT 确定性正弦阈值（或 SFT 周期性过零）的机理，设计出在低/中/高全倍率上系统性超越 `B0_par_raw` + SFT 的伪 1.5-bit 算子？**

### 1.2 SFT 1-bit PWM 基带传递函数与三大误差源
在 SFT 体制下，叠加于上采样复回波 $S^\uparrow$ 的阈值具有高幅、高频的解析形式（在阶段二 `B0_par_raw` 锁定参数中，于 $Q=5\text{–}9$ 收敛至 $\mathrm{STR} \in \{-4, -5\}\text{ dB}$、$f_r/B_r \in [1.5, 1.8]$、$f_a/B_a \in [1.2, 1.6]$）：

$$U[n] = A_{\mathrm{SFT}}\,e^{j(\phi_r[n_r] + \phi_a[n_a] + \phi_0)}$$

单通道（I 路或 Q 路）执行加阈值 1-bit 量化 $b[n] = \operatorname{sgn}_+(S_c[n] + A_{\mathrm{SFT}}\cos\theta[n])$。对于准静态信号分量 $s \in [-A_{\mathrm{SFT}}, A_{\mathrm{SFT}}]$，该 1-bit 脉宽调制（PWM）输出一个载波周期内的平均基带传递函数与一次谐波（载波基频）幅值为：

$$c_0(s) = \frac{2}{\pi}\arcsin\!\left(\frac{s}{A_{\mathrm{SFT}}}\right) = \frac{2}{\pi}\left[\frac{s}{A_{\mathrm{SFT}}} + \frac{1}{6}\left(\frac{s}{A_{\mathrm{SFT}}}\right)^3 + \frac{3}{40}\left(\frac{s}{A_{\mathrm{SFT}}}\right)^5 + \cdots\right], \qquad a_1^{\mathrm{B0}}(s) = \frac{4}{\pi}\sqrt{1 - \left(\frac{s}{A_{\mathrm{SFT}}}\right)^2}$$

由此可直接定位 **`B0_par_raw` + SFT** 与盲置零 **`R1_par_raw`** 的三大性能损耗根源：
1. **大信号区 $\arcsin$ 上翘非线性畸变**：当 $|s|/A_{\mathrm{SFT}}$ 增大时，$\frac{2}{\pi}\arcsin(s/A_{\mathrm{SFT}})$ 呈三次及更高奇次项上翘，造成动态范围与高散射区相对响应失真；
2. **载波基频（$m=1$）与高次谐波（$m \ge 2$）泄漏与混叠**：SFT 将大量能量集中在载波基频与谐波上；在有限过采样倍率 $Q$ 与给定 $f/B$ 下（阶段二 `B0` 锁频 $f/B$ 在 $Q \ge 2$ 达 $1.0\text{–}2.0$），高次谐波或未被彻底压制的载波边带极易混叠或残留进成像支持域；
3. **盲 `R1` 在低/中 $Q$ 下的短游程过度抹零**：高频正弦载波使每个载波半周期仅跨越少数采样点；`R1_par_raw` 对所有邻点跳变边两端无差别置 `0`，在 `RSFT`/`ASFT` 下导致 $Q=1$ 置零率高达 **76.5%–93.5%**、$Q=2\text{–}4$ 仍高达 **62.0%–71.1%**，大量有效窄脉冲（1–2 点游程）被整段抹零。

---

## 2. 严格 1-bit / 伪 1.5-bit 契约与防信息泄漏红线

### 2.1 两级表示契约
1. **前端采集与下行传输（100% 纯 1-bit）**：保持 `S_up + U -> sgn_+`，传输比特流严格为 $\{+1, -1\}$，零额外 ADC 与通信代价。
2. **距离压缩前重构（100% 严格三电平伪 1.5-bit）**：I/Q 实部与虚部各自独立处理，矩阵尺寸与数值精度（`single`/`double`）保持不变，所有输出元素严格属于 $\{+1, 0, -1\}$（由运行时断言 `zc1bit:InvalidPseudo15Output` 强制保证）。

### 2.2 工程先验分级与接口隔离红线
本专题同时纳入**全盲算子**与**单位模载波辅助算子**两类工作体制：
- **全盲算子**（`B0_par_raw`、`R1_par_raw`、`P15_2_rlp_she`）：不需要任何载波频率、幅度或初相 $\phi_0$ 先验，直接仅从 `{+1, -1}` 序列的局部游程结构计算三电平输出，对接收端初相失配天然完全免疫。
- **单位模载波辅助算子**（`P15_1_carrier_agree`、`P15_3_arcsin_lin`、`P15_4_dual_domain`）：利用雷达已知工作频率/初相遥测生成的单位模解析载波 $u[n] = \exp(j\theta[n])$（实虚部仅提供符号/相位指向，严格在 $[-1, 1]$ 内），绝不涉及任何与回波强度有关的量。

核心算子 [apply_sft_p15_reconstruct.m](../src/quantization/apply_sft_p15_reconstruct.m) 在代码结构上锁定四条防泄漏红线：
- **签名禁区**：`[channel_out, diag] = apply_sft_p15_reconstruct(channel_1bit, carrier, method_key, options)` **严禁**传入未量化上采样回波 `signal_up` 与幅度尺度 $\hat\sigma$。
- **入口 1-bit 硬校验**：检验 `channel_1bit` 实虚部所有元素严格属于 $\{+1, -1\}$（失败抛出 `zc1bit:Invalid1BitInput`）。
- **出口三电平硬断言**：检验 `channel_out` 实虚部所有元素严格属于 $\{-1, 0, +1\}$（失败抛出 `zc1bit:InvalidPseudo15Output`）。
- **样本严格隔离**：全部超参寻优与锁定仅在 35 奇数校准样本上进行；35 偶数测试样本严禁参与任何超参或 SFT 参数决策。

---

## 3. 重构方案数学定义（6 方法同场对比）

记单通道（I 路或 Q 路，各自独立处理）的 1-bit 输入序列为 $b[n] \in \{+1, -1\}$，对应的归一化载波实部或虚部分量为 $u[n] \in [-1, 1]$。沿当前工作主轴（`RSFT` 沿 dim 1 距离向，`ASFT` 沿 dim 2 方位向），定义：
- 相邻过零边指示：$e[k] = \mathbb{I}(b[k] \neq b[k+1])$；
- 边中点归一化载波：$\bar{u}[k] = \frac{1}{2}(u[k] + u[k+1])$；
- 位置 $n$ 所属的连续同号游程长度：$L[n] \ge 1$（沿主轴独立统计，列间/行间严格切断，防止二维边界误合并）。

| 序号 | 方法 Key | 类别与定位 | 需 $\phi_0$/载波 | 核心数学定义 | 寻优超参 | 精确退化边界 |
| :---: | :--- | :--- | :---: | :--- | :---: | :--- |
| 1 | `B0_par_raw` | 原始 1-bit 基线 | 否 | 保持 $b[n] \in \{+1, -1\}$ 不变 | 无 | — |
| 2 | `R1_par_raw` | 盲双沿置零对照 | 否 | 对所有 $e[k]=1$ 的过零边两端并行置 `0` | 无 | — |
| 3 | `P15_1_carrier_agree` | **方案 1**：载波置信度门控置零 | 是 | 若 $b[n]\cdot u[n] > \tau$ 则 $y[n]=0$，否则 $y[n]=b[n]$ | $\tau \in [0.1, 1.0]$ | $\tau=1 \Rightarrow B_0$ |
| 4 | `P15_2_rlp_she` | **方案 2**：游程保护 + 可调半宽对称置零（全盲 SHE-ZC） | **否（全盲）** | 仅对满足 $L[k] \ge L_{\min}$ 且 $L[k+1] \ge L_{\min}$ 的过零边，向两侧各对称置零至多 $w$ 点 | 半宽 $w \in \{1, 2\}$<br>最小游程 $L_{\min}$ | $L_{\min}=\infty \Rightarrow B_0$<br>$(w=1, L_{\min}=1) \Rightarrow R_1$ |
| 5 | `P15_3_arcsin_lin` | **方案 3**：$\arcsin$ 反翘曲单侧置零 | 是 | 对 $|\bar{u}[k]| > \tau_{\mathrm{lin}}$ 的过零边，仅将与 $-\operatorname{sgn}(\bar{u}[k])$ 同号的多数侧置 `0` | $\tau_{\mathrm{lin}} \in [0.2, 1.0]$ | $\tau_{\mathrm{lin}}=1 \Rightarrow B_0$ |
| 6 | `P15_4_dual_domain` | **方案 4**：波峰置信度 + 游程保护边沿双域联合置零 | 是 | 取 `P15_1`（参数 $\tau$）与 `P15_2`（固定 $w=1$，参数 $L_{\min}$）置零掩膜的逻辑并集 | 二维 $(\tau, L_{\min})$ | $(1,\infty)\Rightarrow B_0$<br>$(\tau,\infty)\Rightarrow \text{P15\_1}$<br>$(1,L_{\min})\Rightarrow \text{P15\_2}(w=1)$ |

### 3.1 各方案机理剖析

1. **方案 1：`P15_1_carrier_agree`（载波置信度门控置零 / 逐点 MMSE 最优三电平投影）**

   由统计判决理论，当已知 $b[n]$ 与 $u[n]$ 时，逐点条件期望修正的振幅与偏置取决于 $b[n]\cdot u[n]$：若 $b[n]$ 与载波波峰强同号（$b[n]\cdot u[n] > \tau$），则大概率仅仅是由阈值极性本身切换产生的，而包含的信号信息极易受量化误差与高能载波泄漏主导；将其置 `0` 是连续逐点 MMSE 权重在 $\{+1, 0, -1\}$ 上的最优三电平硬投影。它在小信号区 $|s| < A_{\mathrm{SFT}}\tau$ 上完整保留了 100% 无偏基带线性系数 $c_0(s)$，同时在 $m=1$ 载波峰谷处精确挖孔，大幅压低载波基频能量与高幅端 $\arcsin$ 上翘响应。

2. **方案 2：`P15_2_rlp_she`（游程保护 + 可调半宽对称置零 / 全盲 SHE-ZC）**

   针对低/中 $Q$ 下短游程窄脉冲被盲 `R1` 抹杀、以及高 $Q$ 下需抑制奇次高次谐波的问题，本全盲算子沿工作主轴仅允许左右游程同时满足 $L[k] \ge L_{\min}$ 与 $L[k+1] \ge L_{\min}$ 的“稳健过零边”触发置零，且在边沿左右各向内对称扩展置零 $w \in \{1, 2\}$ 点（对距离边界偏移 $d \in \{1,\dots,w\}$ 亦要求相邻游程 $\ge d$，确保每个游程内只覆盖到有效过零邻域）。在理想双电平周期脉冲调制下，在每个过零边界两侧形成宽度对称的零电平平台 $2\delta$，使其傅里叶谐波系数满足三电平**选择性谐波消除（Selective Harmonic Elimination, SHE）**规律：
   $$a_m^{\text{3-level}} = a_m^{\mathrm{B0}}\cos(m\delta)$$
   既保持信号极性对称与近零直流偏置，又能按 $w$ 与 $L_{\min}$ 在不同倍率下压制高频开关谐波（同时以 $L_{\min}$ 保护 1–2 点窄脉冲不被破坏）。

3. **方案 3：`P15_3_arcsin_lin`（$\arcsin$ 反翘曲单侧置零）**

   SFT 在过零点满足 $S + U \approx 0$，即高载波幅值 $|\bar{u}[k]| > \tau_{\mathrm{lin}}$ 处发生过零时，真实信号极大且具有明确先验符号 $\operatorname{sgn}(S) = -\operatorname{sgn}(\bar{u}[k])$：
   - 当 $\bar{u}[k] < -\tau_{\mathrm{lin}}$（真实信号显著为正 $S > 0$）时，仅将过零边中的 **`+1` 侧** 置 `0`；
   - 当 $\bar{u}[k] > +\tau_{\mathrm{lin}}$（真实信号显著为负 $S < 0$）时，仅将过零边中的 **`-1` 侧** 置 `0`。

   与盲 `R6`/`R7` 会产生高达 `0.24–0.39` 的严重单向直流偏置不同，`P15_3` 因依赖正弦奇对称的 $\bar{u}[k]$ 符号，在正负信号侧严格反对称削峰，天然保持零均值偏置，并专门在 $|s|/A_{\mathrm{SFT}} > \tau_{\mathrm{lin}}$ 的大信号端下压传输曲线以抵消 $\arcsin$ 的三次上翘。

4. **方案 4：`P15_4_dual_domain`（波峰置信度 + 游程保护边沿双域联合置零）**

   联合“相位波峰域”的 `P15_1`（载波基频抑制 + 高端线性化，参数 $\tau$）与“时间游程域”的 `P15_2`（稳健过零边对称窄零槽谐波消除，固定 $w=1$、参数 $L_{\min}$），对二者的置零掩膜直接取逻辑并集，并在 2D 正交网格上进行直接联合寻优。其网格直接内嵌了 `B0_par_raw`、单域 `P15_1` 与单域 `P15_2(w=1)` 三重边界，可检验双域互补协同的净增益。

### 3.2 关键数值实现与反面设计警示
- **二维边界与向量化游程统计**：[apply_sft_p15_reconstruct.m](../src/quantization/apply_sft_p15_reconstruct.m) 在计算游程长度时于每列首行强制重置 `is_start(1, :) = true`（方位向先转置后处理），彻底避免列尾与下一列首行的同号点误合并为一个跨列超长游程。全局禁用 MATLAB 内置 `sign`。
- **反面设计警示（为什么不采用“先低通滤波再硬量化回 $\{+1, 0, -1\}$”）**：对 1-bit + SFT 序列直接低通滤波再按固定阈值重新硬量化为三电平，会抹除高频正弦抖动线性化效应，在通带内重新引入近似 ZT 的三值硬截断奇次谐波；因此本专题全部 4 套方法均直接在 1-bit/三电平离散符号域内作机理驱动的门控或游程变换。

---

## 4. 实验协议、载波抑制诊断与初相失配消融

### 4.1 轨道 A：固定阶段二 `B0_par_raw` 锁定 SFT 参数（主实验）
* **单元范围（共 18 单元）**：
  * 距离单向 `RSFT`：`rb_R1A1_range_rsft` ~ `rb_R9A1_range_rsft`（9 单元）；
  * 方位单向 `ASFT`：`ab_R1A1_azimuth_asft` ~ `ab_R1A9_azimuth_asft`（9 单元）。
* **前端 SFT 参数读取与冻结**：直接读取 [locked_params.csv](../results_zc1bit/locked/locked_params.csv) 中 `VariantKey == "B0_par_raw"` 在上述 18 单元上的锁定 $(S_{\mathrm{dB}}, f_r/B_r, f_a/B_a, \phi_0)$，保证采集端 1-bit 流与基线完全一致。
* **35 校准样本超参寻优网格（每单元共 68 候选）**：
  * `B0_par_raw`、`R1_par_raw`：各 **1** 点（无超参）；
  * `P15_1_carrier_agree`：$\tau \in \{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0\}$（**12** 点）；
  * `P15_2_rlp_she`：$(w=1, L_{\min} \in \{1, 2, 3, 4, 5, 6, \infty\})$（7 点）+ $(w=2, L_{\min} \in \{5, 6, 7\})$（3 点），共 **10** 点；
  * `P15_3_arcsin_lin`：$\tau_{\mathrm{lin}} \in \{0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0\}$（**9** 点）；
  * `P15_4_dual_domain`：直接 2D 网格 $\tau \in \{0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0\}$（7 点）× $L_{\min} \in \{2, 3, 4, 5, \infty\}$（5 点，固定 $w=1$），共 **35** 点。
* **锁参准则**：按 35 校准样本平均 `SSIM` 降序首选，依次以平均 `PSNR` 降序、平均 `AbsENLDiff` 升序、平均 `ZeroRatio` 升序与 `HyperparamTag` 字典序破局。

### 4.2 频域机理诊断：带外泄漏与精确载波基频抑制比
在 5 个代表倍率组 `R1A1`/`R4A1`/`R9A1`/`R1A4`/`R1A9` 对应的 **6 个单元**（`R1A1` 含距离和方位各 1 单元）× 7 个典型场景测试样本（每场景首张偶数样本：`SampleID` 为 `2, 12, 22, 32, 42, 52, 62`）上计算：
1. **距离压缩域带外能量泄漏**：`OffSupportRatio`、`RangeOffRatio`、`AzimuthOffRatio`（支持域阈值 `leakage_threshold_ratio = 0.35`）；
2. **距离压缩前上采样序列载波基频抑制比**：直接对进入 RC 前的严格三电平序列 $y = y_I + j\,y_Q$ 计算精确单频载波 $u[n] = e^{j\theta[n]}$ 投影能量占比及相对 `B0_par_raw` 的抑制增益（dB）：
   $$C_I = \mathrm{mean}\!\left(y_I\,e^{-j\theta}\right), \quad C_Q = \mathrm{mean}\!\left(y_Q\,e^{-j\theta}\right), \quad \mathrm{CarrierToneRatio} = \frac{2\left(|C_I|^2 + |C_Q|^2\right)}{\mathrm{mean}\!\left(y_I^2 + y_Q^2\right) + \varepsilon}$$
   $$\mathrm{CarrierSuppression\_dB} = 10\log_{10}\frac{\mathrm{CarrierToneRatio}_{\mathrm{B0}} + \varepsilon}{\mathrm{CarrierToneRatio} + \varepsilon}$$

### 4.3 工程消融：接收端初相失配鲁棒性测试
在 4 个代表单元（`rb_R4A1_range_rsft`、`rb_R9A1_range_rsft`、`ab_R1A4_azimuth_asft`、`ab_R1A9_azimuth_asft`）上，前端保持真实载波生成 1-bit 流，接收端注入全局初相失配 $\Delta\phi_0$：

$$\Delta\phi_0 \in \{0^\circ, 5^\circ, 10^\circ, 15^\circ, 30^\circ, 45^\circ, 90^\circ\}$$

- **全盲方法免疫机制**：`B0_par_raw`、`R1_par_raw` 与 **`P15_2_rlp_she`** 对接收端载波与 $\Delta\phi_0$ 天然完全免疫（实现中每样本仅计算 1 次后直接复用到全部 7 个相偏角度），形成“相偏无关鲁棒线”，与相干方法 `P15_1`/`P15_3`/`P15_4` 在 $\Delta\phi_0$ 增大时的衰减边界形成直接工程对比。

### 4.4 后续可选扩展（轨道 B）
若轨道 A 显示某严格伪 1.5-bit 方法（如 `P15_2` 或 `P15_4`）在改变 SFT 调制深度/频偏时具备更高上限，后续可在代表倍率扩展 $(S_{\mathrm{dB}}, f)$ 与重构超参的联合寻优，不计入本次默认单核执行预算。

---

## 5. 评价指标、统计检验、隔离架构与验收标准

* **全量测试指标（35 测试样本）**：自动复用 Stage 2 同款 GT 缓存与固定 $64\times 64$ 均匀区 ROI（`filed` + `port` 共 10 个测试样本），统计 `PSNR`、`SSIM`、`Entropy`、`ENL`、`AbsENLDiff`、`ZeroRatio`、`ProtectedEdgeRatio`（被游程条件拦截保护的过零边比例）与 `AbsMeanBias`。
* **双对照配对 Wilcoxon 符号秩检验（`signrank`，强制列向量）**：对全部 18 单元输出两种显著性检验对照表：
  * **`vs_B0`**：5 个非基线方法 vs `B0_par_raw`（$18 \times 5 = 90$ 行）；
  * **`vs_R1`**：4 个 `P15_*` 新方案 vs `R1_par_raw`（$18 \times 4 = 72$ 行），共计 **162** 行。
* **目录与断点严格隔离**：全部缓存、锁参、统计结果与图表存于 `results_zc1bit/checkpoints/sft_p15/` 与 `results_zc1bit/sft_p15/`（含 `search/`、`evaluate/`、`stats/`、`robustness/`、`figures/`），采用单核确定性运行与按 `UnitKey` 扁平 `-v7` 原子落盘，绝不修改 Stage 1 / Stage 2 或 `ta_sft` 任何产物。
* **本地 MATLAB 验收标准**：
  1. **静态检查与规范**：全部修改/新增代码满足 `checkcode` 零严重警告、UTF-8 无 BOM 与 `git diff --check` 清洁；
  2. **单元测试全通过**：[TestSftP15Reconstruction.m](../tests/TestSftP15Reconstruction.m) 6 项独立测试（含三电平硬契约、精确退化边界、真值表、轴对称性、反翘曲奇对称、双域并集与独立小规模端到端寻优/评估/断点恢复）及全仓 `tests/` 全部通过；
  3. **主实验产物维度完整性**：
     * `sft_p15_search_summary.csv`：$18 \times 68 = 1224$ 行；`sft_p15_locked_params.csv`：$18 \times 6 = 108$ 行；
     * `sft_p15_evaluate_detail.csv`：$18 \times 6 \times 35 = 3780$ 行；`sft_p15_evaluate_summary.csv`：$18 \times 6 = 108$ 行；
     * `sft_p15_leakage_summary.csv`：$6 \times 7 \times 6 = 252$ 行；
     * `sft_p15_paired_wilcoxon.csv`：162 行；
     * `sft_p15_phase_robustness_summary.csv`：$4 \times 7 \times 6 = 168$ 行；核心评估指标列零 NaN / Inf。
