#!/usr/bin/env python3
"""Strict author-input gate + deterministic algebra audit, never full replication.

Usage: python docs/validation/verify_manuscript22.py [author_model.json] --output report.json
Exit 2: missing inputs (NOT_RUN); 1: invalid/failed; 3: partial checks passed,
full manuscript replication NOT_RUN. No path through this runner claims replication.
Input wraps the package model schema in `model`, with `provenance` containing
`kind: author_supplied_original` and `source_reference`, plus `verification`
containing six author-specified `dose_levels` and 22 `outcome_weights`.
Provenance is an author declaration, not independent proof of originality.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'skills/tsymperturb/scripts'))
from tsymperturb import Model, weights  # noqa: E402

MODEL_KEYS = ('source_labels', 'outcome_labels', 'B', 'mu1', 'sigma1', 'psi',
              'intercept', 'anchors', 'modules', 'orientation', 'higher_is_worse', 'scale')
BRIDGES = [('Insomnia', 'Fatigue'), ('Fatigue', 'Concentration difficulty'),
           ('Rumination', 'Anxiety'), ('Rumination', 'Insomnia'),
           ('Anxiety', 'Palpitations'), ('Pain', 'Insomnia'),
           ('Kinesiophobia', 'Activity avoidance')]
KNOWN_LABELS = {'Fatigue', 'Rumination', 'Anxiety', 'Pain', 'Anhedonia', 'Insomnia', 'Worry',
                'Concentration difficulty', 'Palpitations', 'Kinesiophobia', 'Activity avoidance'}
FULL_STAGES = ['250000_draws_per_target_monte_carlo', '600_dataset_finite_sample_recovery',
               'population_tvpps_and_benchmarks', 'four_wave_and_repeated_intervention']


def base_report():
    return {'schema_version': 1, 'status': 'NOT_RUN', 'replication_status': 'NOT_RUN',
            'scope': 'Input design consistency and deterministic package algebra only',
            'missing_inputs': [], 'checks': [], 'executed': [],
            'software': {'python': platform.python_version(), 'numpy': np.__version__},
            'runner_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'core_sha256': hashlib.sha256((ROOT/'skills/tsymperturb/scripts/tsymperturb.py').read_bytes()).hexdigest(),
            'not_run': FULL_STAGES.copy(),
            'disclaimer': 'No substitute 22-node reconstruction is an original-model replication.'}


def check_design(model):
    """Inspect declared original-scale arrays; summary consistency cannot prove identity."""
    checks = []
    def add(name, ok, actual, expected):
        checks.append(dict(name=name, passed=bool(ok), actual=actual, expected=expected))
    add('nodes', model.p == 22, model.p, 22)
    add('modules', len(set(model.modules)) == 4, len(set(model.modules)), 4)
    if model.p != 22:
        return checks
    add('known_source_labels', KNOWN_LABELS.issubset(model.labels),
        sorted(KNOWN_LABELS-set(model.labels)), 'all eleven source-verifiable labels present; actual lists missing labels')
    b = model.B
    off = b[~np.eye(model.p, dtype=bool)]
    # Literal non-zero counts: no hidden topology threshold is introduced.
    for name, actual, expected in [('nonzero_offdiagonal', int(np.count_nonzero(off)), 53),
                                    ('positive_offdiagonal', int(np.sum(off > 0)), 47),
                                    ('negative_offdiagonal', int(np.sum(off < 0)), 6)]:
        add(name, actual == expected, actual, expected)
    # Published ranges are rounded to two decimals; compare each endpoint within half a unit.
    for name, values, expected in [('autoregressive_range', np.diag(b), [0.43, 0.71]),
                                   ('baseline_mean_range', model.mu1, [1.44, 2.64]),
                                   ('baseline_sd_range', np.sqrt(np.diag(model.sigma1)), [0.63, 0.83]),
                                   ('residual_sd_range', np.sqrt(np.diag(model.psi)), [0.38, 0.52])]:
        actual = [float(values.min()), float(values.max())]
        add(name, np.all(np.abs(np.array(actual)-expected) <= 0.005+1e-12), actual,
            {'rounded_endpoints': expected, 'absolute_tolerance': 0.005})
    rho = float(np.max(np.abs(np.linalg.eigvals(b))))
    add('spectral_radius', abs(rho-0.827) <= 0.0005+1e-12 and rho < 1, rho,
        {'rounded_value': 0.827, 'absolute_tolerance': 0.0005, 'stable': True})
    add('independent_residuals', np.allclose(model.psi, np.diag(np.diag(model.psi)), atol=1e-12, rtol=0),
        float(np.max(np.abs(model.psi-np.diag(np.diag(model.psi))))), 'diagonal Psi (atol=1e-12)')
    delta = float(np.max(np.abs(model.intercept-(np.eye(model.p)-b)@model.mu1)))
    add('stationary_intercept', delta <= 1e-10, delta, 'a=(I-B)mu1; atol=1e-10')
    add('zero_original_scale_anchor', np.allclose(model.anchors, 0, atol=1e-12, rtol=0),
        model.anchors.tolist(), 'zero anchors on original generating scale')
    for source, outcome in BRIDGES:
        actual = None if source not in model.labels or outcome not in model.labels else float(b[model.index(outcome), model.index(source)])
        add('bridge:'+source+'->'+outcome, actual is not None and actual != 0, actual,
            'nonzero B[outcome,source]; magnitude and sign not inferred from the named bridge list')
    return checks


def check_algebra(model, dose_levels, outcome_weights):
    """Compare public state() outputs to separate direct mean propagation and identities."""
    d = np.asarray(dose_levels, dtype=float)
    if d.shape != (6,) or not np.all(np.isfinite(d)) or len(set(d)) != 6 or np.any((d < 0)|(d > 1)):
        raise ValueError('verification.dose_levels must contain six distinct finite values in [0,1]; do not guess manuscript values')
    w = weights(outcome_weights, model.p)
    oracle_error = dose_error = pair_error = covariance_error = 0.
    singles = []
    for i, label in enumerate(model.labels):
        full = model.state({label: 1})['standardized_response']
        singles.append(full)
        for dose in d:
            out = model.state({label: float(dose)})
            changed_mu = model.mu1.copy()
            changed_mu[i] = model.anchors[i]+(1-dose)*(model.mu1[i]-model.anchors[i])
            expected = (model.mu2-(model.intercept+model.B@changed_mu))/model.sd2
            oracle_error = max(oracle_error, float(np.max(np.abs(out['standardized_response']-expected))))
            dose_error = max(dose_error, abs(float(np.dot(w, out['standardized_response']-dose*full)/w.sum())))
            diagonal = np.eye(model.p); diagonal[i, i] = 1-dose
            expected_cov = model.B@diagonal@model.sigma1@diagonal@model.B.T+model.psi
            covariance_error = max(covariance_error, float(np.max(np.abs(out['outcome_covariance']-expected_cov))))
    pairs = 0
    for i, j in itertools.combinations(range(model.p), 2):
        out = model.state({model.labels[i]: 1, model.labels[j]: 1})['standardized_response']
        mask = np.ones(model.p, dtype=bool); mask[[i,j]] = False
        if w[mask].sum() <= 0:
            raise ValueError('every pair-excluded outcome set needs positive total weight')
        # Equation 36 is vector additivity; Eq37 uses the same excluded comparison set.
        pair_error = max(pair_error, float(np.max(np.abs(out-singles[i]-singles[j]))),
                         abs(float(np.dot(w[mask], (out-singles[i]-singles[j])[mask])/w[mask].sum())))
        pairs += 1
    actual = dict(targets=model.p, dose_checks=model.p*len(d), pairs=pairs,
                  direct_mean_max_error=oracle_error, dose_efficacy_max_error=dose_error,
                  pair_additivity_max_error=pair_error, linked_covariance_max_error=covariance_error)
    return {'passed': max(oracle_error,dose_error,pair_error,covariance_error) <= 1e-10,
            'absolute_tolerance': 1e-10, 'actual': actual,
            'published_comparison_only': {'dose_max_error': 2.78e-17, 'pair_max_error': 2.78e-17},
            'comparison_note': 'Published machine-precision residuals are observations, not exact equality targets.'}


def verify(payload):
    report = base_report()
    if not isinstance(payload, dict):
        report.update(status='INVALID_INPUT', error='top-level JSON must be an object'); return report, 1
    model_data = payload.get('model', {})
    if not isinstance(model_data, dict):
        report.update(status='INVALID_INPUT', error='model must be an object'); return report, 1
    report['missing_inputs'] = ['model.'+k for k in MODEL_KEYS if k not in model_data]
    for group, keys in [('provenance', ['kind', 'source_reference']),
                        ('verification', ['dose_levels', 'outcome_weights'])]:
        value = payload.get(group, {})
        if not isinstance(value, dict):
            report.update(status='INVALID_INPUT', error=group+' must be an object'); return report, 1
        report['missing_inputs'] += [group+'.'+k for k in keys if k not in value]
    if report['missing_inputs']:
        return report, 2
    provenance = payload['provenance']
    if provenance['kind'] != 'author_supplied_original' or not isinstance(provenance['source_reference'],str) or not provenance['source_reference'].strip():
        report.update(status='INVALID_INPUT', error='A cited author-supplied-original declaration is required. Synthetic checker fixtures are not manuscript inputs.')
        return report, 1
    report['provenance'] = provenance
    report['provenance_caveat'] = 'Declaration recorded; identity with original arrays is not independently established.'
    try:
        model = Model(**{k:model_data[k] for k in MODEL_KEYS})
        report['checks'] = check_design(model)
        report['executed'].append('model_validation_and_manuscript_summary_checks')
        if not all(c['passed'] for c in report['checks']):
            report['status'] = 'DESIGN_MISMATCH'; return report, 1
        report['algebra'] = check_algebra(model, **payload['verification'])
        report['executed'].append('22_target_six_dose_231_pair_deterministic_checks')
        report['status'] = 'PARTIAL_VERIFIED' if report['algebra']['passed'] else 'ALGEBRA_FAILED'
        return report, 3 if report['algebra']['passed'] else 1
    except (ValueError, TypeError, OverflowError, FloatingPointError, np.linalg.LinAlgError) as error:
        report.update(status='INVALID_INPUT', error=str(error)); return report, 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs='?', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.input is None or not args.input.is_file():
        report, code = base_report(), 2
        report['missing_inputs'] = ['author-supplied original model JSON' if args.input is None else str(args.input)]
    else:
        try:
            raw = args.input.read_bytes()
            report, code = verify(json.loads(raw))
            report['input_sha256'] = hashlib.sha256(raw).hexdigest()
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
            report, code = base_report(), 1
            report.update(status='INVALID_INPUT', error=str(error))
    content = json.dumps(report, indent=2, allow_nan=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(content)
    print(content, end='')
    return code

if __name__ == '__main__':
    raise SystemExit(main())
