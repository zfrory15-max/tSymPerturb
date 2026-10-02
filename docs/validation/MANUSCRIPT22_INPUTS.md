# Original 22-symptom input gate

This checker is **not a complete manuscript reproduction runner**. The supplied manuscript and TeX archives do not include the original ordered labels, module map, model arrays or executable generator. Matching reported counts/ranges cannot establish that a model is the original model. Do not reverse-engineer a replacement from the published summaries.

Run from the package root:

```sh
python docs/validation/verify_manuscript22.py author_model.json --output report.json
```

Omitting the input, using a missing file, or omitting required fields returns `NOT_RUN` and exit code 2. Invalid data or a design/algebra discrepancy returns exit code 1. Even when design and algebra pass, the result is only `PARTIAL_VERIFIED`, exit code 3, with `replication_status: NOT_RUN`. There is deliberately no successful full-replication exit code.

## Input contract

The JSON must contain:

- `model`: the existing package model schema (`source_labels`, `outcome_labels`, `B`, `mu1`, `sigma1`, `psi`, `intercept`, `anchors`, `modules`, `orientation`, `higher_is_worse`, `scale`). Arrays must be the author-supplied **original generating-scale** model, not fitted standardised estimates. `B[j,i]` maps source i to outcome j; both label orders must be identical. The intercept is `(I-B)mu1`, and the original-scale anchors are zero
- `provenance`: `kind` must be `author_supplied_original`; `source_reference` must identify the supplied original source. This declaration is recorded but cannot independently authenticate original arrays
- `verification`: six distinct `dose_levels` between 0 and 1, and 22 nonnegative `outcome_weights` with positive total (and positive total after excluding every pair). Obtain the original values from the author; the checker does not invent them

Do not substitute `examples/synthetic_model.json`: it is a newly constructed 22-node illustrative demo, not an authenticated original system. Its stable spectral radius is approximately .767515, rather than the manuscript’s reported .827; matching selected counts and ranges does not establish original provenance. Do not substitute the synthetic 22-dimensional unit-test fixture either: its names and numbers are deliberately arbitrary and fail manuscript design checks.

## What is checked

- 22 nodes and four modules; 53 literal nonzero off-diagonal coefficients, 47 positive and six negative; no undisclosed edge threshold
- Published autoregressive, baseline-mean, baseline-SD and residual-SD range endpoints; half-unit rounding tolerance of 0.005 plus numeric slack
- Spectral radius within 0.0005 of the rounded 0.827, and less than 1; this is a summary consistency check, not matrix authentication
- Original zero anchors, diagonal residual covariance, stationary intercept, dimensions and covariance validity; all eleven names verifiable in supplied sources; seven named source-to-outcome bridges, without inventing their coefficients
- All 22 targets at six supplied doses, direct mean/covariance propagation, linear dose efficacy, all 231 two-target response and common-comparison-set additivity checks. Absolute implementation tolerance is `1e-10`; the paper's realised `2.78e-17` residuals are comparison values, not exact universal acceptance thresholds

The original factor-loading construction and full label/module identity cannot be validated from published summaries alone.

## Six framework components and limits

The dedicated unit tests use explicitly synthetic 22-dimensional arrays to exercise vKO, vKD, vDP, edge-level and node-centred communication blocking, combinations and sequence recursion. They are implementation tests, not original-result tests. In particular, sequence recursion checks the stated synthetic intervention policy; it does not reproduce the article's original top-eight sequence optimisation.

Communication utility reproduction additionally needs the original horizon and discount; combination utility needs original partner sets; population tVPPS needs the full scoring configuration. Four-wave and repeated-intervention replication additionally need the original discount, scale conventions and intervention policy. The 22-by-250,000 Monte Carlo run and 600 independent finite-sample datasets need the original generator, random-stream design and analysis conventions. These stages remain explicitly `NOT_RUN` in this checker. They require a subsequent reproduction workflow once the originals are available.

`manuscript22_spec.json` records the source design and historical result targets. `manuscript22_readiness.json` records the current missing-input run. `manuscript22_execution_manifest.json` distinguishes executed checker/unit tests from original-model work that has not run. Source statistics are never replaced by synthetic test outcomes.
