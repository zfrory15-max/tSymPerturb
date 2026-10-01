# GitHub upload and release guide / GitHub 上传与发布指南

[English README](../README.md) · [中文 README](../README.zh-CN.md)

## English

### Before uploading

1. Extract the archive. The repository root is the directory containing `README.md`, `skills/`, `examples/` and `tests/`. Upload its contents, not the ZIP as the only repository file and not an extra enclosing directory.
2. Review the license decision with the authors. This bundle does not grant a new open-source license. Confirm authorship, publication status and citation metadata; never substitute the manuscript submission identifier for a verified DOI/arXiv identifier.
3. Keep the manuscript PDF, participant data, credentials, local environments and generated personal outputs out of the public repository. The included synthetic model is a demonstration, not the original study data.
4. Run the tests and documented example. Resolve errors before release. Review changed files rather than relying on a test badge alone.

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
python skills/tsymperturb/scripts/tsymperturb.py examples/synthetic_model.json --output outputs/synthetic_results.json
```

### Recommended: upload with Git

Create an empty GitHub repository under your chosen account or organization. Replace `YOUR_ACCOUNT` and `YOUR_REPOSITORY` below with your actual destination. These commands are instructions for you; preparing this bundle does not publish it.

```bash
cd /path/to/extracted/tSymPerturb-skills
git init
git add .
git status --short
# Review every staged file before continuing.
git commit -m "Add tSymPerturb skill and reference implementation"
git branch -M main
git remote add origin https://github.com/YOUR_ACCOUNT/YOUR_REPOSITORY.git
git push -u origin main
```

`git add .` includes hidden files unless ignored. Inspect the staged list before committing. Do not use `git add *` as a substitute: shell wildcard expansion can omit dotfiles.

### If using GitHub's web upload

Use **Add file → Upload files**, then select the extracted repository contents. Preserve directory paths and verify that `.github/`, `.gitignore` and any other hidden configuration files are included. This bundle uses `skills/tsymperturb/` as its canonical skill folder; it does not require a duplicate `.agents/` copy. If you later add `.agents/`, preserve that hidden directory too. File managers often hide dot-prefixed folders; a Git-based upload is the safer option when the browser/file picker omits them.

After uploading, check that:

- the root page renders the README, with working English/Chinese links
- `skills/tsymperturb/SKILL.md`, its scripts and references are present
- the example JSON and tests are present
- `.github/workflows/` is present if included in the bundle; inspect the actual Actions result
- no private source manuscript, real participant records, secrets or local `.venv/` were uploaded

### Suggested About text

English, fewer than 350 characters:

> Portable tSymPerturb Agent Skill and reference code for time-indexed virtual perturbations in longitudinal symptom networks: temporal state and transition operators, seven-utility tVPPS, auditable synthetic examples, and bilingual documentation. Model-implied hypotheses, not causal treatment effects.

Suggested topics: `tsymperturb`, `symptom-networks`, `longitudinal-analysis`, `virtual-perturbation`, `agent-skills`, `python`.

### Before a formal release

Select the code/content license; verify citation metadata; reconcile the implementation with the original research code; archive the actual parameter files and figure source data needed for any claimed manuscript reproduction; choose a real version/tag; and record the environment and validation scope. Add a DOI only after an archive service actually assigns it. A public repository alone does not establish validation or a persistent release identifier.

## 中文

### 上传前

1. 解压文件包。包含 `README.md`、`skills/`、`examples/` 与 `tests/` 的目录即仓库根目录。应上传其中的内容，不能只把 ZIP 作为仓库唯一文件，也不要额外套一层文件夹。
2. 与作者确认许可选择。本包不自动授予新的开源许可。核实作者、发表状态和引用元数据，不能将投稿标识当作已经验证的 DOI/arXiv 编号。
3. 不要公开论文源 PDF、参与者数据、凭据、本地环境或个人分析输出。所附合成模型是演示数据，不是原始研究数据。
4. 运行测试及文档中的示例，解决错误后再发布。检查具体变更内容，不仅查看测试徽标。

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
python skills/tsymperturb/scripts/tsymperturb.py examples/synthetic_model.json --output outputs/synthetic_results.json
```

### 推荐使用 Git 上传

先在自己的账号或组织下新建空仓库，将上方 Git 命令中的 `YOUR_ACCOUNT` 与 `YOUR_REPOSITORY` 换为真实目标地址，再在解压后的根目录执行。命令仅供操作参考；生成文件包并不代表已经发布。

`git add .` 会纳入未被忽略的隐藏文件。提交前逐项检查暂存列表。不要以 `git add *` 替代，通配符可能遗漏点号开头的文件。

### 使用网页上传时

选择 **Add file → Upload files**，上传解压后的仓库内容。保留目录层级，检查 `.github/`、`.gitignore` 等隐藏配置文件是否完整。本包以 `skills/tsymperturb/` 为唯一 Skill 源目录，无需另建重复 `.agents/`；如果以后自行增加 `.agents/`，也应保留该隐藏目录。文件管理器可能隐藏点号目录，网页选择器无法纳入时建议改用 Git。

上传后检查：

- 根页面正确显示 README，中英文链接可互相跳转
- `skills/tsymperturb/SKILL.md`、脚本与 references 完整
- 示例 JSON 和测试文件完整
- 若包内包含 `.github/workflows/`，仓库中也应保留，并检查实际 Actions 执行结果
- 未公开私人论文源文件、真实参与者记录、密钥或 `.venv/`

### GitHub About 简介

中文，不超过 350 字符：

> 面向纵向症状网络的 tSymPerturb Agent Skill 与参考代码，支持带时间索引的状态扰动、转移阻断、七维 tVPPS、可审计合成示例及中英文文档。结果用于生成模型推导的干预靶点假设，不代表因果治疗效应。

建议 topics：`tsymperturb`、`symptom-networks`、`longitudinal-analysis`、`virtual-perturbation`、`agent-skills`、`python`。

### 正式发布前

确定代码与文档许可，核实引用信息，与原始研究代码核对；若声称复现论文，应归档所需真实生成参数及图表源数据。选择实际版本号与 tag，记录环境和验证范围。只有归档服务真实分配 DOI 后才加入 DOI。公开仓库本身不等于完成方法验证，也不自动产生持久发布标识。
