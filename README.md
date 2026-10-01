# tSymPerturb Skills

**English** | [简体中文](README.zh-CN.md)

**tSymPerturb** turns a fitted longitudinal symptom transition model into auditable, time-indexed virtual-perturbation hypotheses. This repository provides a portable Agent Skill, a NumPy reference implementation, synthetic examples and bilingual documentation based on the supplied tSymPerturb manuscript.

> **Scientific boundary:** outputs are model-implied hypotheses, not identified causal treatment effects. A two-wave model cannot identify repeated treatment order. The included example verifies software behavior; it does not reproduce the manuscript's original 22-node simulation.

## What is included

- [Portable Agent Skill](skills/tsymperturb/SKILL.md) for designing, running and auditing a tSymPerturb analysis
- [Reference Python script](skills/tsymperturb/scripts/tsymperturb.py) operating on an already fitted model
- [Four-symptom synthetic input](examples/synthetic_model.json) with explicit scale, direction and scoring settings
- [English method guide](docs/METHODS.md) and [Chinese method guide](docs/zh-CN/METHODS.md)
- [Bilingual GitHub upload guide](docs/GITHUB_UPLOAD.md)

The repository presentation follows the bilingual skill-plus-reference-code approach of [SymPerturb](https://github.com/zfrory15-max/SymPerturb). The temporal operators, utility definitions and interpretation rules follow the tSymPerturb manuscript; cross-sectional SymPerturb utilities are not interchangeable with them.

## Implemented scope

| Interface | Available functionality |
|---|---|
| JSON command line | Seven raw utilities, normalized utilities, tVPPS and ranks; full-dose response profiles; signed pair increments and additivity checks; optional state-dose requests; baseline moments and analysis settings |
| Python API | Anchor transformation; location–scale source perturbation and covariance propagation; edge/source-node blocking; reference-profile signed prediction changes; finite-horizon unsigned route capacity; bounded-normal marginal means; multi-wave propagation and fixed additional-improvement recursion; discounted horizon utility; exact additive target-acquisition search for at most nine candidates |
| Optional scikit-learn adapter | Matched two-wave Lasso estimation and participant-level complete-pipeline bootstrap; per-replicate fitting diagnostics and failure records |
| Method guidance | Sensitivity analysis and reporting requirements |

**Not implemented:** missing-data processing, clustered resampling, causal identification, clinical outcome validation, or reproduction of the paper's original numerical study. Bounded expectations and multi-wave functions are separate APIs; the default scoring CLI remains the unbounded linear reference analysis.

## Quick start

Run these commands from the extracted repository root. Python 3.10 or later and pip are required.

```bash
python -m venv .venv
```

Activate the environment:

```bash
# macOS / Linux
source .venv/bin/activate
# Windows PowerShell
.venv\Scripts\Activate.ps1
```

Install runtime dependencies and run the synthetic example:

```bash
python -m pip install -r requirements.txt
python skills/tsymperturb/scripts/tsymperturb.py examples/synthetic_model.json --output outputs/synthetic_results.json
```

For development checks:

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
```

This is a portable skill script, not an editable-install Python distribution. The example contains four anonymous synthetic symptoms and no participant records.

## Install as a skill

Copy the complete `skills/tsymperturb/` directory, including its scripts and references, into your agent's skills directory. For a local Codex installation on macOS/Linux:

```bash
mkdir -p "${CODEX_HOME:-$HOME/.codex}/skills"
cp -R skills/tsymperturb "${CODEX_HOME:-$HOME/.codex}/skills/"
```

If a `tsymperturb` skill already exists there, review and back it up before replacing it. Keep this repository as the editable source. Other skill-compatible agents may use a different installation path.

Example request:

> Use tSymPerturb to audit my fitted T1→T2 model. Check matrix orientation and clinical anchors, calculate the seven utilities, and report tVPPS with signed adverse effects and the assumptions that limit interpretation.

## Use your own fitted model

1. Fit and validate the CLPN outside the core scoring script. Retain symptom order, scale transformations, intercept, source moments and residual covariance.
2. Orient every symptom so higher means worse. Define anchors on the original scale and transform them using the fitted T1 mean and SD if needed.
3. Prepare JSON following [the included example](examples/synthetic_model.json). `B` must be **outcome-by-source**: `B[j,i]` is the path `i(T1) → j(T2)`. Source and outcome labels must have identical order. Do not silently transpose a source-by-outcome export.
4. Prespecify candidates, partners, modules, outcome/utility weights, thresholds, blocking fraction, propagation horizon and damping. Example choices are illustrative, not recovered paper settings or validated clinical defaults.
5. Run the command above with your JSON path. Inspect raw responses, harmful signed changes and pair increments before interpreting composite ranks.
6. For sampling uncertainty in independent, complete matched two-wave data, use the optional adapter below, or build a validated external resampling workflow that repeats all transformations, fitting and scoring.

See the [input/output contract](skills/tsymperturb/references/input-output.md) for required fields and Python API signatures. The scoring script uses fixed baseline **model-implied** T2 SDs derived from `B Σ1 Bᵀ + Ψ`; document this convention when comparing results with empirical-SD analyses.

The CLI accepts JSON, not CSV. If your external estimator exports a coefficient CSV, explicitly reorder and convert it first. For example, a table with rows `T2_A,T2_B` and columns `T1_A,T1_B` already has the required direction; a table with source rows requires an explicit transpose. Preserve the diagonal autoregressive coefficients.

## Optional fitting and bootstrap

The separate adapter fits standardized per-outcome Lasso models (`alpha=0.03`, no intercept by default) and resamples matched participant rows. It repeats wave-specific standardization, anchor transformation, fitting, all seven utilities, normalization and ranking in every replicate.

```bash
python -m pip install -r requirements-fitting.txt
python skills/tsymperturb/scripts/fit_tsymperturb.py examples/synthetic_paired.json --n-boot 1000 --seed 42 --output outputs/bootstrap.json
```

The example contains 60 deterministic synthetic matched participants and four symptoms. Use `--n-boot 10` for a short smoke test only; ten replicates do not support reliable uncertainty summaries. Input rows must already be matched by participant and symptoms oriented consistently; there is no identifier matching, missing-data repair or clustered bootstrap.

Outputs include point estimates, replicate diagnostics, percentile summaries, top-1/3/5 selection frequencies and all failed attempts. Frequencies use successful replicates as their denominator; failures can bias summaries and must be inspected. Wilson intervals describe finite-bootstrap Monte Carlo error, not clinical or causal certainty.

**Fitting caveat:** Lasso residuals need not be empirically orthogonal to the source variables. The core's model-implied covariance `B Σ1 Bᵀ + Ψ` may therefore differ from empirical T2 covariance. The adapter records source–residual covariance diagnostics; assess this discrepancy rather than assuming the fitted moments are exact. See [fitting and uncertainty](skills/tsymperturb/references/fitting-and-uncertainty.md) for conventions and limits. This participant bootstrap is distinct from the manuscript's independent-dataset simulation study.

## Strategy API example

```bash
python examples/strategy_demo.py
```

The script prints JSON illustrating stationary propagation, repeated additional improvements, signed inhibitory-edge blocking, hypothetical bounded outcomes and additive target-acquisition ordering. It repeats one synthetic matrix for three transitions, so these trajectories are **stationarity-assumption simulations**, not observed multi-wave estimates. No empirical treatment order is inferred.

## Interpret the output

`outputs/synthetic_results.json` contains:

- `candidates` and `dimensions`: row/column labels for score arrays
- `raw`, `normalized`, `tvpps`, `ranks`: seven utilities, 0–100 within-set scores, composite and competition ranks
- `standardized_responses`: candidate-by-outcome full-dose improvements
- `pairs`: shared-outcome single/pair utilities, signed increments and additive contrasts
- `baseline`, `source_labels`, `scale`: unperturbed follow-up moments and model context
- `config`, `dose_results`: analysis choices and requested state perturbations

The seven dimensions are downstream efficacy, cross-symptom spillover, breadth, cross-module reach, communication-block value, combination value and spillover fraction. Constant dimensions receive a neutral score of 50. tVPPS is relative to the candidate set; its numerical values are not directly comparable across different cohorts or candidate sets. Robustness is not an eighth utility.

## Checks that govern interpretation

- Under the fixed, unbounded linear reference model, dose response is linear and independent source-state perturbations are additive
- Communication-block route capacity is unsigned; inspect signed reference-state effects separately
- Exact knockout changes the source state while keeping the fitted transition fixed; do not refit or re-standardize a constant predictor
- Use observed wave-specific matrices for observed multi-wave propagation; repeating one matrix beyond its observed interval assumes stationarity
- Target-acquisition order under a decision objective does not establish biological treatment order

See the [method guide](docs/METHODS.md) for equations and reporting requirements.

## Reproducibility, citation and license

See [source status](SOURCE.md), [contribution checks](CONTRIBUTING.md) and [licensing decision](LICENSING.md).

This is a manuscript-derived reference implementation. Numerical unit tests assess stated equations and edge cases; they do not validate clinical use. Independent review and reconciliation with the authors' original research code remain necessary before claiming publication-level reproduction.

The supplied manuscript is titled *tSymPerturb converts longitudinal symptom networks into time-indexed intervention strategies*, by Zheng Zhu, Junwen Yu, Tiantian Hu, Zhongfang Yang and Jiaqing Wang. Its submission identifier is not treated here as a public arXiv identifier or DOI. Cite the confirmed manuscript record and the actual repository release when these become available.

No redistribution license is chosen by this bundle. The authors should select appropriate code/content licensing before public release. The supplied manuscript PDF and private data are not included. See the [upload guide](docs/GITHUB_UPLOAD.md) for the final checks.
