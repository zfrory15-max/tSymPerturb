# tSymPerturb

[English](README.md) | **简体中文**

**tSymPerturb** 用于分析源症状状态或有向路径改变后，影响如何沿已拟合的纵向症状网络传播。框架包括六项功能：**虚拟敲除、虚拟敲降、虚拟剂量扰动、通信阻断、组合扰动和干预顺序优化**。七个效用维度描述候选靶点的不同价值，tVPPS 用于候选集合内的综合排序，不确定性单独报告。

## 框架总览

[![原稿图1：tSymPerturb 的源状态扰动、转移路径阻断与干预策略分析](docs/figures/figure1.png)](docs/figures/figure1.png)

*图 1 来自作者提供的原稿文件包，经作者授权原样使用。点击图片可查看原始分辨率。[图片来源](docs/figures/README.md)。*

源状态分析改变 T1 症状，通信阻断改变 T1→T2 转移路径，组合与顺序分析比较联合或有序策略。结果是给定模型下的干预假设，临床获益仍需独立验证。

七维效用综合满剂量状态响应、源节点阻断和预设配对价值；剂量曲线、边效应与时序轨迹单独输出。分析可从已拟合模型或匹配的两波观测数据开始。

## 仓库包含什么

| 资源 | 用途 |
|---|---|
| [Agent Skill](skills/tsymperturb/SKILL.md) | 辅助设计、执行、审计和解释时序虚拟扰动分析 |
| [NumPy 参考实现](skills/tsymperturb/scripts/tsymperturb.py) | 固定模型下的扰动、分布矩传播、七维效用与 JSON 命令行 |
| [可选拟合适配器](skills/tsymperturb/scripts/fit_tsymperturb.py) | 匹配两波数据的 Lasso 拟合与参与者层面的完整流程 bootstrap |
| [22 症状模型示例](examples/synthetic_model.json) / [配对观测示例](examples/synthetic_paired.json) | 包含 22 个症状、四个模块的构造示例 |
| [策略演示](examples/strategy_demo.py) | 多波传播、重复额外改善、有界结局和靶点纳入顺序 |
| [中文方法说明](docs/zh-CN/METHODS.md) / [英文方法说明](docs/METHODS.md) | 数学假设、实现约定与报告要求 |
| [逐式方法规范](skills/tsymperturb/references/methods.md) | 论文公式编号、估计对象和核验要点 |
| [独立验证报告](docs/VALIDATION_REPORT.zh-CN.md) | 数值 oracle、Monte Carlo、工作流检查与测试边界 |

仓库组织与说明方式参考 [SymPerturb 中文 README](https://github.com/zfrory15-max/SymPerturb/blob/main/README_zh-CN.md)。tSymPerturb 使用有向、已拟合的时序转移，不采用横断面协方差求逆或平衡态网络操作；两者的七维效用定义不能直接互换。

## 时序虚拟扰动模型

两波参考模型为：

```text
X2 = a + B X1 + ε
μ2 = a + B μ1
Σ2 = B Σ1 Bᵀ + Ψ
```

`B[j,i]` 表示 **T1 的症状 i 指向 T2 的症状 j**。列表示源症状的输出路径，行表示结局；对角线保留自回归持续效应。所有症状必须统一为数值越高、状态越差。

锚点 `c` 是预设的源症状参照状态，应与拟合模型处于相同量纲。若拟合时进行标准化，则用 `c_z = (c_raw − μ1_raw) / SD1_raw` 转换原始量表锚点。标准化后的 0 通常表示样本均值，不能直接当作无症状状态。

一般 location–scale 状态扰动分别改变目标的均值位置与残余变异。对目标集合 S：

```text
X1,S* = cS + Dμ(μ1,S − cS) + Dσ(X1,S − μ1,S)
R = μ2 − μ2* = B(μ1 − μ1*)
Σ2* = B Σ1* Bᵀ + Ψ
```

非目标源坐标不变；目标与非目标之间的协方差也必须同步缩放。联动参考映射使用 `Dμ = Dσ = diag(1 − dose)`。当前 API 也支持独立指定 [0,1] 内的残余尺度乘数，包括仅改变均值、仅改变尺度的查询；不实现任意方差放大的 location–scale 扩展。

## 六项扰动分析

### 1. 虚拟敲除（t-vKO）

将源症状 i 固定在锚点：`X1,i* = ci`。目标方差及其交叉协方差归零，满剂量响应为 `R_i = B[:,i](μ1,i − ci)`。分析保持 `B`、`a` 和 `Ψ` 不变，不对常数列重新标准化或拟合。前向传播允许扰动后的源协方差为奇异矩阵。

### 2. 虚拟敲降（t-vKD）

将源症状部分移向锚点：

```text
X1,i* = ci + (1 − d)(X1,i − ci),  0 < d < 1
```

均值与锚点的距离乘以 `1−d`，方差乘以 `(1−d)²`，目标与非目标的协方差乘以 `1−d`。前述一般位置–尺度接口还允许分别控制均值与残余变异。

### 3. 虚拟剂量扰动（t-vDP）

虚拟剂量扰动使用**剂量参数化的算子族** `T_i(d): X1,i → ci + (1−d)(X1,i−ci)`，其中 `d ∈ [0,1]`。d=0 时保持原状态，中间剂量对应敲降，d=1 对应敲除。在预设剂量网格上应用这一算子族，即得到剂量–响应曲线。论文将网格评估归为强度–响应分析程序，其底层状态映射仍是算子（Methods 公式 23–27）。

固定、无界线性模型满足 `R_i(d)=dR_i(1)`。非线性曲线需要额外的非线性机制，例如显式的结局边界。这里的剂量表示扰动强度，尚未对应经过校准的临床剂量。

### 4. 通信阻断

边级阻断使用 `B*[j,i] = (1−q)B[j,i]`；源节点阻断使用 `B*[:,i] = (1−q)B[:,i]`。参考设置包含自回归，`cross_only=True` 可保留 `B[i,i]`。

固定源参照状态和截距后，有符号效应为 `(B − B*) x_ref`。通信阻断效用则衡量无符号有限时域传播容量的变化：

```text
Q_H(B) = 1ᵀ [Σ(h=1…H) (γ |B|)^h] 1
```

先逐元素取绝对值，再做矩阵幂，并报告 `H`、`γ` 和 `q`。容量下降仍可能伴随症状恶化；源总体均值为零时，总体平均的有符号阻断效应也为零。

### 5. 组合扰动

同时改变多个源状态，将联合靶点与两个单靶点放在同一结局集上比较。配对分析排除两个靶点自身的结局，采用相同权重，并保留相对较优单靶点的有符号增量和可加性交互对比。

线性参考模型满足 `R_{i,k}=R_i+R_k`，源症状相关时也成立。因此，共同结局集上的线性效用具有零可加性交互对比；仅凭优于单靶点，不能认定存在协同作用。

### 6. 干预顺序优化

有实际多波转移时，一次改善按时间顺序传播为 `B_(t+h−1)…B_t δ_t`。重复的转移前额外改善满足 `Δ⁻_(t+1)=B_t(Δ⁻_t+u_t)`。其中 `u_t` 表示额外改善，重复敲除还需根据当时状态另行定义。将一个两波矩阵重复用于 `B^h` 则增加了平稳性假设。

靶点纳入求解器处理另一类决策问题：在固定长度下，哪个纳入顺序使可加折扣目标最大？当前支持不超过九个预设候选和非负成本。22 症状示例选取 tVPPS 前八名，搜索两个靶点的纳入顺序。该结果描述模型内的决策次序，不能据此推断生物学治疗先后。

完整定义和接口见[逐式方法说明](skills/tsymperturb/references/methods.md)与 [API 说明](skills/tsymperturb/references/input-output.md)。

## 七个效用维度与 tVPPS

参考状态效用在满剂量下计算。令 `Δji = Rj←i / s2,j`；分母使用**未扰动、模型隐含的随访标准差** `s2,j = sqrt((BΣ1Bᵀ+Ψ)[j,j])`，不能随扰动重新变化。正值表示预测改善，负值表示预测恶化。

| 效用维度 | 参考定义 | 关键约定 |
|---|---|---|
| **1. 下游效力** | `Σj wj Δji / Σj wj` | 包含所有 T2 结局，**包括靶点自身的随访结局**（公式 47） |
| **2. 跨症状溢出效应** | `Σ(j≠i) wj Δji / Σ(j≠i) wj` | 排除靶点结局，并在剩余结局上重新归一化权重（公式 48） |
| **3. 广度** | p−1 个非靶点结局中 `Δji ≥ τ` 的比例 | 不加权计数；参考 `τ=0.05 SD`（公式 49） |
| **4. 跨模块覆盖** | 其他模块中，模块内平均改善 `≥ τm` 的模块比例 | 模块内不加权平均，再对模块等权计数；参考 `τm=0.03 SD`（公式 50） |
| **5. 通信阻断价值** | `[Q_H(B)−Q_H(B_block_i)] / Q_H(B)` | 阻断整列输出；参考 `q=0.8`；H、γ 需显式设定（公式 35） |
| **6. 组合价值** | 对预设伙伴 k 平均 `max(0, Iik)` | `Iik = Gik^(−ik) − max(Gi^(−ik), Gk^(−ik))`；三个效用均排除 i、k，且使用相同权重（公式 38–39） |
| **7. 溢出占比** | `Σ(j≠i) wj max(Δji,0) / [Σj wj max(Δji,0)+10⁻¹²]` | 衡量正向加权获益的分布，不是有符号净获益（公式 51） |

每个维度在预设候选靶点集合内进行 0–100 min–max 标准化。常数维度赋 **50 分**；只有一个候选时也采用这一规则。tVPPS 是七个标准化维度的加权平均，要求权重非负、总权重大于 0。等权只是方法学参考选择，不代表经过验证的临床偏好。

应始终同时报告七维原始效用与综合分数。**Robustness 不是第八个效用维度。** 无界线性参考模型中的剂量效率和低剂量响应性与效力冗余，不作为额外的时序评分维度。较高的正向溢出占比仍可能伴随其他结局恶化，必须查看完整有符号响应。

## 安装与快速开始

需要 **Python 3.10+**。在仓库根目录运行：

```bash
python -m venv .venv
# macOS / Linux
source .venv/bin/activate
# Windows PowerShell 改用：
# .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python skills/tsymperturb/scripts/tsymperturb.py examples/synthetic_model.json --output outputs/synthetic_results.json
```

本仓库是可移植 Skill 脚本，不是可编辑安装的 Python 发行包；请使用上述命令，而非 `pip install -e .`。`--output` 为必填项，指定路径已存在时会覆盖，请按需选择新文件名。

运行策略示例与测试：

```bash
python examples/strategy_demo.py
python -m pip install -r requirements-dev.txt
# 同时安装独立数值审计所需的 SciPy
python -m pip install -r requirements-fitting.txt
python -m pytest -q
python docs/validation/run_audit.py
python docs/validation/audit_synthetic22.py
python docs/validation/mc_audit.py
```

### 安装为 Agent Skill

把整个 `skills/tsymperturb/` 目录连同脚本与 references 一起复制到代理的 skills 目录。macOS/Linux 下的本地 Codex 示例：

```bash
mkdir -p "${CODEX_HOME:-$HOME/.codex}/skills"
cp -R skills/tsymperturb "${CODEX_HOME:-$HOME/.codex}/skills/"
```

已有同名 Skill 时先检查并备份；其他兼容代理可能使用不同路径。可以这样提问：

> 使用 tSymPerturb 审计我的 T1→T2 拟合模型。核对矩阵方向、锚点与量纲，计算七个效用维度，保留不利的有符号效应，并报告 tVPPS、组合增量及顺序分析依赖的假设。

## 使用自己的模型或数据

### 路线 A：输入已拟合模型

参照 [synthetic_model.json](examples/synthetic_model.json) 准备 JSON。核心命令行接收 JSON，不直接接收 CSV。

| 输入类别 | 必要内容 |
|---|---|
| 身份与方向 | 顺序完全一致的源/结局标签、`orientation="outcome_by_source"`、`higher_is_worse=true` 及明确的量纲说明 |
| 拟合模型 | `B`、`intercept`、源均值 `mu1`、源协方差 `sigma1`、残差协方差 `psi` |
| 科学定义 | 已转换到拟合量纲的 `anchors`、每个症状的模块归属，且至少两个模块 |
| 评分方案 | 明确的候选集、非自身伙伴集、结局权重、七维效用权重、传播时域与折扣 |
| 可选查询 | 状态剂量、独立残余尺度乘数、阈值或阻断比例的覆盖设置 |

1. 拟合并验证纵向模型，记录波次间隔、参与者匹配和预处理方式。
2. 统一症状方向，用拟合源均值和 SD 转换原始量表锚点。
3. 核对矩阵方向和标签顺序。源为行、结局为列的导出需要显式转置；保留自回归对角线。
4. 预设候选、伙伴、模块、权重、阈值、时域与折扣。
5. 运行 CLI，先查看有符号响应和原始效用，再解释排序。

完整七维接口要求至少三个症状、每个比较结局集均保留正总权重，以及非零传播容量。对无效或未定义的输入，程序报错而不是编造分数。详细字段见[输入输出契约](skills/tsymperturb/references/input-output.md)。

### 路线 B：输入匹配的两波观测

```bash
python -m pip install -r requirements-fitting.txt
python skills/tsymperturb/scripts/fit_tsymperturb.py examples/synthetic_paired.json --n-boot 1000 --seed 42 --output outputs/bootstrap.json
python docs/validation/fit_audit.py
```

适配器默认对标准化后的每个结局拟合 Lasso（`alpha=0.03`，无截距）。每次参与者 bootstrap 均重新执行配对行重抽样、分波标准化、锚点转换、拟合、扰动、七维效用、标准化与排序。输入必须已正确按参与者匹配；当前不包含身份匹配、缺失数据修复或聚类 bootstrap。

构造示例包含 250 名配对参与者、22 个症状和四个模块。`--n-boot 10` 仅适合 smoke test，不能支持可靠的不确定性结论。必须检查失败重复：选择频率和百分位摘要以成功重复为分母，失败可能造成偏差。Wilson 区间反映有限 bootstrap 次数的 Monte Carlo 误差，不代表临床确定性或经验证的区间覆盖率。

惩罚回归残差不一定与源变量经验正交，因此 `BΣ1Bᵀ+Ψ` 可能不同于经验 T2 协方差。适配器保留源–残差协方差诊断，不能假定结构协方差约定已由数据验证。详见[拟合与不确定性说明](skills/tsymperturb/references/fitting-and-uncertainty.md)。

## 主要输出

| JSON 字段 | 应检查的内容 |
|---|---|
| `candidates`、`dimensions` | 所有得分矩阵的行列顺序 |
| `raw`、`normalized`、`tvpps`、`ranks` | 七维原始效用、标准化效用、综合分数及降序竞争排名 |
| `standardized_responses` | 候选×结局的满剂量有符号响应；不能忽略恶化 |
| `pairs` | 共同结局集上的两个单靶点和联合效用、有符号增量、可加性对比 |
| `baseline` | 未扰动的模型隐含均值、协方差及 SD |
| `config`、`scale`、`source_labels` | 可追溯的设置与坐标约定 |
| `dose_results` | 请求剂量下的源/结局均值、协方差和响应 |

核心 CLI 输出 JSON。有界正态期望、边级效应、多波传播和纳入顺序搜索通过 Python API 调用，具体见策略演示。

## 示例数据与验证

所有可运行示例统一使用**构造的 22 症状、四模块系统**。S01–S22 为示例标签，A/B/C/D 四组分别包含 6/6/5/5 个症状。模型输入、250 名参与者的配对观测及策略演示使用同一系统；六档示例剂量为 0、0.2、0.4、0.6、0.8 和 1。

```bash
# 按固定设置重新生成示例输入与元数据
python examples/generate_synthetic22.py
```

构造采用论文描述的 53 条非对角路径（47 正、6 负）和参数范围。谱半径约为 0.7675，论文报告值为 0.827。矩阵、症状/模块分配和观测均为新建示例，不是原始研究数据。[生成设置与来源记录](examples/synthetic22_metadata.json)使用种子 20261001。

[验证报告](docs/VALIDATION_REPORT.zh-CN.md)分别记录公式测试、构造数据的工作流检查及原始模型复现状态。参与者 bootstrap 反映给定数据集条件下的抽样变异，不能替代论文的独立数据集恢复实验。

### 原始模型验证

原始数值复现仍为 **NOT_RUN（未运行）**。所提供材料缺少完整的原始生成矩阵、有序症状/模块映射及机器可读图表数据。[strict22 验证入口](docs/validation/verify_manuscript22.py)要求作者提供原始输入，不以构造示例替代。

```bash
# 无原始输入：预期退出码 2（NOT_RUN）
python docs/validation/verify_manuscript22.py --output outputs/manuscript22_readiness.json

# 提供原始生成参数与验证设置后
python docs/validation/verify_manuscript22.py author_model.json --output outputs/manuscript22_checks.json

# 构造输入的校验器与算子回归测试
python -m pytest -q tests/test_manuscript22_contract.py
```

字段要求和实际执行状态见[输入指南](docs/validation/MANUSCRIPT22_INPUTS.md)、[就绪报告](docs/validation/manuscript22_readiness.json)与[执行清单](docs/validation/manuscript22_execution_manifest.json)。入口检查来源声明、维数、已报告参数约束和确定性恒等式。有限检查通过返回 **3 / PARTIAL_VERIFIED**；输入无效或检查失败返回 **1**。来源声明本身不能独立认证数组真实性。

该入口不执行论文的每靶点 250,000 次抽样、600 个独立数据集恢复实验、基准排序或四波/重复干预结果复现。这些研究层面的复现需要原始输入，不能由软件测试代替。

## 结果解释

报告 tVPPS 时，同时保留有符号响应和原始效用。tVPPS 是候选集合内的相对分数，不宜直接跨队列、网络或候选集合比较。Bootstrap 与敏感性分析描述既定分析下的不确定性，不能证明临床有效性。

有界正态期望是可选扩展，不能以单纯截断均值代替。当前实现未涵盖缺失处理、潜在混杂控制、序数/非线性模型拟合、临床验证或一般约束下的治疗排程。

## 引用、来源与许可证

方法依据：Zheng Zhu、Junwen Yu、Tiantian Hu、Zhongfang Yang、Jiaqing Wang，*tSymPerturb converts longitudinal symptom networks into time-indexed intervention strategies*，所提供稿件日期为 **2026 年 8 月 17 日**。算子与分析程序见 Methods 公式 13–45；七维效用与 tVPPS 见表 2 和公式 35、38–39、46–53；验证与不确定性见公式 54–58。

来源说明见 [SOURCE.md](SOURCE.md)，贡献与检查要求见 [CONTRIBUTING.md](CONTRIBUTING.md)。稿件投稿标识不作为公开 arXiv 编号或 DOI；请引用经确认的论文记录以及实际使用的仓库 commit/release，不应补入未经核实的 DOI 或正式版本。

本包尚未选择最终代码/内容再分发许可证，见 [LICENSING.md](LICENSING.md)。仓库仅收录获授权的框架图，不分发论文 PDF 或私有数据。正式声称论文级复现前，应与原始研究代码逐项核对、核实生成输入、独立复审分析并归档版本。
