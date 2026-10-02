# Source and implementation status / 来源与实现状态

## English

The governing methodological source is the supplied manuscript:

**tSymPerturb converts longitudinal symptom networks into time-indexed intervention strategies**  
Zheng Zhu, Junwen Yu, Tiantian Hu, Zhongfang Yang and Jiaqing Wang.  
The supplied draft displays the date **17 August 2026** and a submission label. That label is not used as a public arXiv identifier or DOI. No verified publication DOI, release version or public manuscript URL is asserted here.

This repository was assembled as a new reference implementation from that specification. The bilingual presentation was informed by [the SymPerturb repository](https://github.com/zfrory15-max/SymPerturb), inspected on 30 September 2026. This does not imply that cross-sectional SymPerturb's algorithms or utility definitions apply unchanged to tSymPerturb.

The supplied PDF and extracted manuscript text are not redistributed in this repository. The author authorized inclusion of the unchanged manuscript Figure 1; its provenance and checksum are recorded in [docs/figures/README.md](docs/figures/README.md). The 22-symptom, four-module examples are newly constructed synthetic demonstration material with 250 paired observations and seed 20261001. Their parameters and generic labels are recorded in [examples/synthetic22_metadata.json](examples/synthetic22_metadata.json). The original 22-node parameter matrices, full generating configuration and figure source tables were not supplied as a machine-readable research bundle. Accordingly, this repository does not claim to reproduce the manuscript's reported numerical values.

Tests assess the implemented equations and workflows. Before a publication release, the authors should reconcile the reference implementation with their research code, review conventions and defaults, select the license, and verify citation/release metadata. Any later DOI or version should identify a real published record or archived release.

## 中文

本仓库采用的方法学依据为所提供论文：

**tSymPerturb converts longitudinal symptom networks into time-indexed intervention strategies**  
作者：Zheng Zhu、Junwen Yu、Tiantian Hu、Zhongfang Yang、Jiaqing Wang。  
所提供稿件显示日期为 **2026 年 8 月 17 日**，并带有投稿标识。本仓库不将该标识视为公开 arXiv 编号或 DOI，也不声明未经核实的正式发表 DOI、发布版本或公开论文链接。

本仓库依据论文规范重新编写参考实现。双语呈现方式参考了 [SymPerturb 仓库](https://github.com/zfrory15-max/SymPerturb)，查阅日期为 2026 年 9 月 30 日。横断面 SymPerturb 的算法和效用定义不能原样替代 tSymPerturb 的纵向定义。

仓库不再分发原始 PDF 或提取的全文。作者授权收录原稿图 1，图片原样保留，来源与校验值见 [docs/figures/README.md](docs/figures/README.md)。22 症状、四模块示例为新建合成演示数据，包含 250 名配对参与者，种子为 20261001；构造参数与通用标签见 [examples/synthetic22_metadata.json](examples/synthetic22_metadata.json)。原始 22 节点参数矩阵、完整生成配置与图表源数据未以机器可读研究包提供，因此本仓库不声称已复现论文所报告数值。

测试用于核验已实现公式与流程。正式发表前，作者应与研究代码核对，审查量纲约定及默认设定，确定许可证，并核实引用和发布元数据。后续加入的 DOI 或版本号应对应真实公开记录或归档版本。
