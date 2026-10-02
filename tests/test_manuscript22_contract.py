"""Checker unit tests on SYNTHETIC 22-dimensional fixtures, not article replication.

Labels U00..U21, arbitrary four-group map, diagonal/cross coefficients, covariance
and dose grid below are invented solely to exercise dimensions and guardrails.
They deliberately fail original-design checks; never tune them to paper targets.
"""
import importlib.util
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from tsymperturb import Model, block, trajectories

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('verify_manuscript22', ROOT/'docs/validation/verify_manuscript22.py')
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


@pytest.fixture
def synthetic22():
    p = 22
    b = np.eye(p)*0.5
    for i in range(p-1):
        b[i+1,i] = 0.07 if i % 3 else -0.04
    mu = np.linspace(1,2,p)
    return Model(source_labels=[f'U{i:02d}' for i in range(p)],
                 outcome_labels=[f'U{i:02d}' for i in range(p)], B=b, mu1=mu,
                 sigma1=np.eye(p), psi=np.eye(p)*0.25,
                 intercept=(np.eye(p)-b)@mu, anchors=np.zeros(p),
                 modules=[f'fixture_group_{i%4}' for i in range(p)],
                 orientation='outcome_by_source', higher_is_worse=True,
                 scale='SYNTHETIC CHECKER FIXTURE; arbitrary units; not manuscript data')


def test_missing_original_is_not_run():
    result, code = audit.verify({})
    assert code == 2 and result['status'] == 'NOT_RUN'
    assert result['replication_status'] == 'NOT_RUN' and result['missing_inputs']
    assert not result['executed']


def test_missing_file_cli_nonzero(tmp_path):
    result = subprocess.run([sys.executable, str(ROOT/'docs/validation/verify_manuscript22.py'),
                             str(tmp_path/'absent.json')], capture_output=True, text=True)
    assert result.returncode == 2
    assert 'NOT_RUN' in result.stdout


def test_malformed_groups_rejected():
    result, code = audit.verify({'model': []})
    assert code == 1 and result['status'] == 'INVALID_INPUT'


def test_synthetic22_not_original_design(synthetic22):
    checks = audit.check_design(synthetic22)
    assert next(c for c in checks if c['name']=='nodes')['passed']
    assert next(c for c in checks if c['name']=='modules')['passed']
    assert not next(c for c in checks if c['name']=='nonzero_offdiagonal')['passed']
    assert not next(c for c in checks if c['name']=='spectral_radius')['passed']
    assert not all(c['passed'] for c in checks)


def test_six_doses_all_231_pairs_synthetic_only(synthetic22):
    # This grid is a UNIT-TEST choice; original paper's six levels are unavailable.
    result = audit.check_algebra(synthetic22, [0,.2,.4,.6,.8,1], np.ones(22))
    assert result['passed']
    assert result['actual']['dose_checks'] == 132
    assert result['actual']['pairs'] == 231


def test_invalid_dose_grid_rejected(synthetic22):
    with pytest.raises(ValueError, match='six distinct'):
        audit.check_algebra(synthetic22, [0,1], np.ones(22))


def test_vko_vkd_linked_moments_synthetic_only(synthetic22):
    m = synthetic22
    ko, kd = m.state({'U00':1}), m.state({'U00':.4})
    assert ko['source_mean'][0] == 0
    assert np.all(ko['source_covariance'][0] == 0)
    assert kd['source_mean'][0] == pytest.approx(.6*m.mu1[0])
    assert kd['source_covariance'][0,0] == pytest.approx(.36*m.sigma1[0,0])
    np.testing.assert_allclose(kd['response'], .4*ko['response'], atol=1e-12)


def test_communication_edge_and_node_synthetic_only(synthetic22):
    m = synthetic22
    edge = block(m, 'U00', .8, outcome='U01')
    expected = m.B.copy(); expected[1,0] *= .2
    np.testing.assert_allclose(edge, expected, atol=1e-12)
    node = block(m, 'U00', .8)
    np.testing.assert_allclose(node[:,0], .2*m.B[:,0], atol=1e-12)
    np.testing.assert_array_equal(node[:,1:], m.B[:,1:])
    cross = block(m, 'U00', .8, cross_only=True)
    assert cross[0,0] == m.B[0,0]
    np.testing.assert_allclose(cross[1:,0], .2*m.B[1:,0], atol=1e-12)


def test_sequence_recursion_synthetic_only(synthetic22):
    m = synthetic22
    delta = np.zeros(22); delta[0] = 1
    interventions = np.zeros((4,22)); interventions[1,1] = .3
    result = trajectories([m.B]*4, delta, interventions)
    manual = delta.copy()
    for t in range(4):
        manual = m.B@(manual+interventions[t])
        np.testing.assert_allclose(result[t],manual,atol=1e-12)
    # Algebra of the stated synthetic policy; no original sequence/rank replication.


@pytest.mark.parametrize('grid', [[0,0,.4,.6,.8,1], [0,.2,.4,.6,.8,float('nan')], [-.1,.2,.4,.6,.8,1]])
def test_invalid_six_doses_rejected(synthetic22, grid):
    with pytest.raises(ValueError, match='six distinct'):
        audit.check_algebra(synthetic22, grid, np.ones(22))


@pytest.mark.parametrize('kind,reference', [('synthetic_fixture','unit test'), ('author_supplied_original','')])
def test_missing_original_provenance_rejected(kind, reference):
    payload = {'model':{key:None for key in audit.MODEL_KEYS},
               'verification':{'dose_levels':[0,.2,.4,.6,.8,1],'outcome_weights':[1]*22},
               'provenance':{'kind':kind,'source_reference':reference}}
    report, code = audit.verify(payload)
    assert code == 1 and report['status'] == 'INVALID_INPUT'
    assert report['replication_status'] == 'NOT_RUN'


def test_null_weights_rejected(synthetic22):
    with pytest.raises(ValueError):
        audit.check_algebra(synthetic22, [0,.2,.4,.6,.8,1], None)


def test_invalid_json_cli_nonzero(tmp_path):
    path = tmp_path/'malformed.json'; path.write_text('{')
    result = subprocess.run([sys.executable, str(ROOT/'docs/validation/verify_manuscript22.py'), str(path)],
                            capture_output=True, text=True)
    assert result.returncode == 1
    assert 'INVALID_INPUT' in result.stdout and 'NOT_RUN' in result.stdout
