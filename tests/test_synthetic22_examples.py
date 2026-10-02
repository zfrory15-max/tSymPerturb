"""Public examples must use the reproducible constructed 22-symptom system."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from tsymperturb import Model
ROOT=Path(__file__).resolve().parents[1]

def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def test_generator_matches_all_checked_in_examples():
    gen=load_module('generator22',ROOT/'examples/generate_synthetic22.py')
    for filename,expected in zip(['synthetic_model.json','synthetic_paired.json','synthetic22_metadata.json'],gen.generate()):
        actual=json.loads((ROOT/'examples'/filename).read_text())
        # Linear algebra libraries may differ at the last bit across CI OSes.
        def compare(a,b):
            if isinstance(b,dict):
                assert set(a)==set(b)
                for key in b:compare(a[key],b[key])
            elif isinstance(b,list):
                assert len(a)==len(b)
                for left,right in zip(a,b):compare(left,right)
            elif isinstance(b,float):
                assert a==pytest.approx(b,rel=1e-12,abs=1e-12)
            else:assert a==b
        compare(actual,expected)


def test_constructed_model_dimensions_and_design():
    data=json.loads((ROOT/'examples/synthetic_model.json').read_text());m=Model(**data['model'])
    off=m.B.copy();np.fill_diagonal(off,0)
    assert m.p==22 and sorted(np.unique(m.modules,return_counts=True)[1])==[5,5,6,6]
    assert np.count_nonzero(off)==53 and np.sum(off>0)==47 and np.sum(off<0)==6
    assert np.diag(m.B).min()==pytest.approx(.43) and np.diag(m.B).max()==pytest.approx(.71)
    assert max(abs(np.linalg.eigvals(m.B)))<1
    np.testing.assert_allclose(m.mu2,m.mu1,atol=1e-12)
    assert data['provenance']['kind']=='synthetic_illustrative'
    assert data['provenance']['original_replication_status']=='NOT_RUN'
    assert len(data['provenance']['symptom_names'])==len(data['provenance']['module_map'])==22
    paired=json.loads((ROOT/'examples/synthetic_paired.json').read_text())
    assert np.shape(paired['t1'])==np.shape(paired['t2'])==(250,22)


def test_independent_synthetic22_equations():
    result=load_module('audit22',ROOT/'docs/validation/audit_synthetic22.py').audit()
    assert result['pairs']==231 and result['dose_checks']==132 and result['utilities']==7


def test_strategy_uses_22_and_explicit_top8_subset():
    run=subprocess.run([sys.executable,str(ROOT/'examples/strategy_demo.py')],capture_output=True,text=True,check=True)
    data=json.loads(run.stdout)
    assert data['dimensions']=={'symptoms':22,'modules':4,'transitions':4}
    assert np.shape(data['one_time_response_path_raw'])==(4,22)
    assert len(set(data['acquisition_candidate_subset']))==8
    for order in data['target_acquisition']['orders']:
        assert len(order)==2 and set(order)<=set(data['acquisition_candidate_subset'])
