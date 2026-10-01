# Contributing / 贡献指南

## English

Before proposing a change, state the affected estimand and manuscript equation. Keep `B[outcome,source]`, source-scale anchors, fixed baseline outcome scaling and the distinction between state and transition operators explicit. Preserve signed pair increments and report uncertainty separately from the seven utility dimensions.

For changes to the operators, utilities or normalization, add tests covering the relevant algebraic identity and failure cases. For model fitting or resampling changes, test matched participant rows, per-replicate transformation/refitting, failed-replicate reporting and reproducibility. Do not add a causal or clinical efficacy claim based only on fitted CLPN coefficients.

From the repository root:

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
python skills/tsymperturb/scripts/tsymperturb.py examples/synthetic_model.json --output outputs/synthetic_results.json
```

If changing optional adapters, install their documented dependencies and run their tests and examples too. Keep the English and Chinese READMEs aligned. Update API/input contracts when behavior changes. Use synthetic fixtures; do not submit patient data, unpublished manuscripts, secrets or generated local environments. Discuss code/content licensing with the authors before distributing contributions under a new license.

## 中文

提出修改时，先说明受影响的估计对象和论文公式。保持 `B[结局,来源]`、来源量纲的锚点、固定基线结局标准化，以及状态算子与转移算子的区别。保留有符号配对增量，不将不确定性混入七个效用维度。

修改算子、效用或标准化时，应补充对应代数恒等式及失败情形的测试。修改拟合或重抽样流程时，应测试参与者配对、每次重复中的转换与重新拟合、失败重复的报告以及可重复性。不能仅凭 CLPN 系数增加因果或临床疗效声明。

在仓库根目录运行上方命令。修改可选适配器时，还需安装其文档指定的依赖并运行对应测试与示例。中英文 README 保持同步；行为改变时同步更新 API 与输入契约。仅使用合成测试材料，不提交患者数据、未发表论文、凭据或本地环境。以新许可证分发贡献前，应与作者确认许可。
