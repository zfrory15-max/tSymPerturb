# tSymPerturb Skills

[English](README.md) | **简体中文**

**tSymPerturb** 将拟合后的纵向症状转移模型转化为可审计、带有明确时间索引的虚拟扰动假设。本仓库依据所提供的 tSymPerturb 论文，整理可移植的 Agent Skill、NumPy 参考实现、合成示例和中英文文档。

> **科学解释边界：** 输出是模型推导的假设，不能认定为已识别的因果治疗效应。两波模型不能识别反复治疗的先后顺序。自带示例用于核验程序行为，不是论文原始 22 节点模拟的复现。

## 仓库内容

- [Agent Skill](skills/tsymperturb/SKILL.md)：辅助设计、执行和审计 tSymPerturb 分析
- [Python 参考脚本](skills/tsymperturb/scripts/tsymperturb.py)：在已拟合的模型上执行扰动
- [四症状合成输入](examples/synthetic_model.json)：明确量纲、矩阵方向与评分设置
- [英文方法说明](docs/METHODS.md)与[中文方法说明](docs/zh-CN/METHODS.md)
- [中英文 GitHub 上传指南](docs/GITHUB_UPLOAD.md)

仓库呈现参考 [SymPerturb](https://github.com/zfrory15-max/SymPerturb) 的双语文档及 Skill 与参考程序并列结构。时间算子、效用定义和解释边界依据 tSymPerturb 论文编写，不能直接套用横断面 SymPerturb 的效用维度。

## 当前实现范围

| 接口 | 已支持功能 |
|---|---|
| JSON 命令行 | 七个原始效用及标准化分数、tVPPS 与排序、全剂量响应、有符号配对增量及可加性检查、可选状态剂量请求、基线矩与分析设定 |
| Python API | 锚点转换、位置–尺度来源扰动与协方差传播、有向边/来源节点阻断、参考状态的有符号预测变化、有限时域无符号路径容量、有界正态边际期望、多波传播与固定附加改善量递推、折扣时域效用，以及最多九个候选的可加靶点加入顺序精确搜索 |
| 可选 scikit-learn 适配器 | 配对两波 Lasso 估计与参与者层面的完整流程 bootstrap；记录每次重复的拟合诊断及失败情况 |
| 方法指导 | 敏感性分析与报告规范 |

**尚未实现：** 缺失数据处理、聚类重抽样、因果识别、临床效果验证及论文原始数值研究复现。有界期望和多波函数为独立 API；默认评分命令仍采用无界线性参考模型。

## 快速开始

在解压后的仓库根目录执行以下命令，需先安装 Python 3.10 或更新版本及 pip。

```bash
python -m venv .venv
```

激活环境：

```bash
# macOS / Linux
source .venv/bin/activate
# Windows PowerShell
.venv\Scripts\Activate.ps1
```

安装运行依赖并执行合成示例：

```bash
python -m pip install -r requirements.txt
python skills/tsymperturb/scripts/tsymperturb.py examples/synthetic_model.json --output outputs/synthetic_results.json
```

运行开发检查：

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
```

本仓库提供可移植 Skill 脚本，不是通过 editable install 安装的 Python 分发包。示例只有四个匿名合成症状，不包含参与者记录。

## 安装为 Skill

将完整的 `skills/tsymperturb/` 文件夹连同 scripts 与 references 复制到所用智能体的 skills 目录。macOS/Linux 上的本地 Codex 可使用：

```bash
mkdir -p "${CODEX_HOME:-$HOME/.codex}/skills"
cp -R skills/tsymperturb "${CODEX_HOME:-$HOME/.codex}/skills/"
```

若目标位置已有同名 Skill，请先检查并备份，再决定是否替换。仓库保留为可编辑源文件。其他兼容 Agent Skills 的工具可能采用不同安装路径。

示例提示词：

> 使用 tSymPerturb 审计我的 T1→T2 拟合模型。先检查矩阵方向和临床锚点，再计算七个效用维度；报告 tVPPS、有符号不利效应及限制解释的假设。

## 使用自己的拟合模型

1. 在核心评分脚本外拟合并验证 CLPN，保存症状顺序、量纲转换、截距、来源矩和残差协方差。
2. 统一症状方向，使较高值表示较差状态。在原始量表上定义锚点；必要时使用拟合时的 T1 均值和标准差转换。
3. 按[自带示例](examples/synthetic_model.json)准备 JSON。`B` 必须为**结局行、来源列**：`B[j,i]` 表示 `i(T1) → j(T2)`。来源与结局标签须按相同顺序排列，不能静默转置来源行为主的导出结果。
4. 预先确定候选集合、配对集合、模块、结局及效用权重、阈值、阻断比例、传播时域和衰减系数。示例参数仅作演示，并非恢复出的论文设定或已验证的临床默认值。
5. 将命令中的输入路径换为自己的 JSON。解释综合排序前，先检查原始响应、有害的有符号变化和配对增量。
6. 对独立参与者的完整配对两波数据，可使用下方可选适配器评估抽样不确定性；也可建立经过验证的外部重抽样流程，重复全部转换、拟合与评分。

完整字段与 Python API 见[输入输出规范](skills/tsymperturb/references/input-output.md)。评分脚本采用由 `B Σ1 Bᵀ + Ψ` 推导的固定基线**模型隐含** T2 标准差；与采用经验标准差的分析比较时，须说明这一约定。

命令行读取 JSON，不直接读取 CSV。如果外部估计器导出系数 CSV，先明确排序并转换格式。例如行是 `T2_A,T2_B`、列是 `T1_A,T1_B` 的表格已符合方向要求；若来源位于行，则需显式转置。保留对角线上的自身持续系数。

## 可选拟合与 bootstrap

独立适配器对每个标准化 T2 结局拟合 Lasso（默认 `alpha=0.03`、不拟合截距），并对配对参与者行重抽样。每次重复均重新完成各波标准化、锚点转换、拟合、七个效用、标准化和排序。

```bash
python -m pip install -r requirements-fitting.txt
python skills/tsymperturb/scripts/fit_tsymperturb.py examples/synthetic_paired.json --n-boot 1000 --seed 42 --output outputs/bootstrap.json
```

示例包含 60 个确定性生成的合成配对参与者和四个症状。快速检查可改为 `--n-boot 10`，但十次重复不足以形成可靠的不确定性摘要。输入前必须完成参与者身份配对与症状方向统一；适配器不按标识匹配，也不修复缺失值或执行聚类 bootstrap。

输出包含点估计、重复诊断、百分位摘要、进入前 1/3/5 名的频率及全部失败记录。频率分母为成功重复次数；失败可能使摘要偏倚，必须检查。Wilson 区间仅表示有限 bootstrap 次数带来的蒙特卡洛误差，不表示临床或因果确定性。

**拟合约定提醒：** Lasso 残差与来源变量未必经验正交，因此核心模型隐含协方差 `B Σ1 Bᵀ + Ψ` 可以与经验 T2 协方差不同。适配器记录来源–残差协方差诊断，应评估差异，不能假定拟合矩精确重建观测矩。完整约定见[拟合与不确定性说明](skills/tsymperturb/references/fitting-and-uncertainty.md)。参与者 bootstrap 与论文的独立数据集模拟研究不同。

## 策略 API 示例

```bash
python examples/strategy_demo.py
```

脚本向终端输出 JSON，展示平稳传播、反复附加改善、有符号抑制边阻断、假设有界结局和可加靶点加入顺序。示例将同一个合成矩阵重复用于三次转移，因此属于**依赖平稳性假设的模拟**，不是实际观测多波估计，也不推断经验治疗顺序。

## 输出如何解读

`outputs/synthetic_results.json` 包含：

- `candidates`、`dimensions`：得分数组的行列标签
- `raw`、`normalized`、`tvpps`、`ranks`：七个效用、候选集合内 0–100 分数、综合分及保留并列的排序
- `standardized_responses`：候选靶点 × 结局的全剂量改善
- `pairs`：共同结局集上的单靶点/配对效用、有符号增量与加性交互对比
- `baseline`、`source_labels`、`scale`：未扰动的后续矩与模型背景
- `config`、`dose_results`：分析设定及请求的状态扰动

七个维度依次为总下游效力、跨症状溢出、广度、跨模块覆盖、通信阻断价值、组合价值和溢出比例。常数维度取中性值 50。tVPPS 是候选集合内的相对分数，不能跨队列或不同候选集合直接比较。稳健性不是第八个效用维度。

## 解释结果前必须检查

- 在固定、无界的线性参考模型中，剂量响应线性，独立来源状态扰动可加
- 通信阻断路径容量无符号，需另看参考状态的有符号效应
- 精确敲除改变来源状态、保持拟合转移不变；不应重新拟合或标准化常数预测变量
- 观测多波传播使用实际波次转移矩阵；在观测区间外重复同一矩阵依赖平稳性假设
- 决策目标下的靶点加入顺序不能证明生物学治疗顺序

公式与报告要求见[中文方法说明](docs/zh-CN/METHODS.md)。

## 可重复性、引用与许可

另见[来源说明](SOURCE.md)、[贡献检查](CONTRIBUTING.md)和[许可状态](LICENSING.md)。

当前版本是依据论文规范构建的参考实现。数值单元测试核验公式与边界情形，不验证临床应用。在宣称达到发表级复现之前，仍需独立审查并与作者原始研究代码核对。

所提供论文题为 *tSymPerturb converts longitudinal symptom networks into time-indexed intervention strategies*，作者为 Zheng Zhu、Junwen Yu、Tiantian Hu、Zhongfang Yang 和 Jiaqing Wang。稿件中的投稿标识不作为公开 arXiv 编号或 DOI 使用。待正式信息可核实后，引用论文记录及实际仓库发布版本。

本文件包未替作者选择再分发许可证。公开发布前，请确定代码和文档的适用许可。包内不包含所提供的论文 PDF 或私人数据。发布前检查见[上传指南](docs/GITHUB_UPLOAD.md)。
