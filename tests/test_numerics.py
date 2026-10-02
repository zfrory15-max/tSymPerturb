import copy
import json
from pathlib import Path
import numpy as np
import pytest
import tsymperturb as t

ROOT=Path(__file__).resolve().parents[1]
@pytest.fixture
def data():
    return json.loads((ROOT/'examples/synthetic_model.json').read_text())
@pytest.fixture
def model(data):
    return t.Model(**data['model'])

def test_baseline_and_knockout(model):
    zero=model.state({'S02':0})
    np.testing.assert_allclose(zero['source_covariance'],model.sigma1)
    np.testing.assert_allclose(zero['response'],0)
    one=model.state({'S02':1})
    np.testing.assert_allclose(one['response'],model.B[:,1]*(model.mu1[1]-model.anchors[1]))
    np.testing.assert_allclose(one['source_covariance'][1,:],0)
    np.testing.assert_allclose(one['standardized_response'],one['response']/model.sd2)

def test_partial_covariance(model):
    got=model.state({'S02':.4})
    assert got['source_covariance'][1,1] == pytest.approx(.6**2*model.sigma1[1,1])
    assert got['source_covariance'][0,1] == pytest.approx(.6*model.sigma1[0,1])
    assert got['source_covariance'][0,0] == model.sigma1[0,0]
    location=model.state({'S02':.4},{'S02':1})
    np.testing.assert_allclose(got['response'],location['response'])
    np.testing.assert_allclose(location['source_covariance'],model.sigma1)

def test_dose_linearity_and_additivity(model):
    for d in [0,.1,.5,.9,1]:
        np.testing.assert_allclose(model.state({'S01':d})['response'],d*model.state({'S01':1})['response'],atol=1e-15)
    np.testing.assert_allclose(model.state({'S01':.6,'S04':.8})['response'],
      model.state({'S01':.6})['response']+model.state({'S04':.8})['response'],atol=1e-15)

def test_block_orientation(model):
    b=t.block(model,'S01',.8)
    np.testing.assert_allclose(b[:,0],.2*model.B[:,0])
    np.testing.assert_allclose(b[:,1:],model.B[:,1:])
    bc=t.block(model,'S01',.8,cross_only=True)
    assert bc[0,0] == model.B[0,0]
    e=t.block(model,'S01',1,outcome='S02')
    expected=model.B.copy();expected[1,0]=0
    np.testing.assert_allclose(e,expected)
    profile=np.zeros(model.p);profile[0]=2
    expected_response=np.zeros(model.p);expected_response[1]=2*model.B[1,0]
    np.testing.assert_allclose(t.reference_block_response(model,e,profile),expected_response)
    np.testing.assert_allclose(t.reference_block_response(model,e,[0]*22),0)

def test_capacity_hand_calculation():
    b=np.array([[.5,-.2],[.1,.4]])
    a=.7*np.abs(b)
    assert t.propagation_capacity(b,3,.7) == pytest.approx((a+a@a+a@a@a).sum())

def test_seven_utilities_and_pairs(model,data):
    result=t.score(model,**data['scoring'])
    assert result['raw'].shape == (22,7)
    r=model.B[:,0]*(model.mu1[0]-model.anchors[0])/model.sd2
    w=np.array(data['scoring']['outcome_weights'])
    assert result['raw'][0,0] == pytest.approx((r*w).sum()/w.sum())
    assert result['raw'][0,1] == pytest.approx((r[1:]*w[1:]).sum()/w[1:].sum())
    pair=result['pairs'][0]
    r2=model.B[:,1]*(model.mu1[1]-model.anchors[1])/model.sd2
    gi=np.mean(r[2:]);gk=np.mean(r2[2:])
    assert pair['signed_increment'] == pytest.approx(gi+gk-max(gi,gk))
    assert max(abs(p['additive_contrast']) for p in result['pairs']) < 1e-14
    assert np.all((result['tvpps']>=0)&(result['tvpps']<=100))

def test_constant_dimensions_neutral(data):
    cfg=data['model'];cfg.update(B=(.5*np.eye(22)).tolist(),mu1=[1]*22,anchors=[0]*22,sigma1=np.eye(22).tolist(),psi=np.eye(22).tolist())
    m=t.Model(**cfg)
    data['scoring']['outcome_weights']=[1]*22
    r=t.score(m,**data['scoring'])
    np.testing.assert_allclose(r['normalized'],50)
    np.testing.assert_allclose(r['tvpps'],50)
    assert r['ranks']==[1]*22

def test_zero_capacity_rejected(data):
    data['model']['B']=np.zeros((22,22)).tolist()
    with pytest.raises(ValueError,match='Q_H'):
        t.score(t.Model(**data['model']),**data['scoring'])

def test_anchor_transform():
    np.testing.assert_allclose(t.transform_anchors([0,1],[2,3],[2,4]),[-1,-.5])

def test_noncommuting_multiwave_and_recursion():
    a=np.array([[1,1],[0,1.]])
    b=np.array([[1.,0],[1,1]])
    d=np.array([1.,2])
    r=t.trajectories([a,b],d)
    np.testing.assert_allclose(r[1],b@a@d)
    assert not np.allclose(r[1],a@b@d)
    u=np.array([[1,0],[0,2.]])
    np.testing.assert_allclose(t.trajectories([a,b],d,u)[1],b@(a@(d+u[0])+u[1]))
    assert t.horizon_utility(r,[1,2],.8)==pytest.approx((r[0]@[1,2]+.8*r[1]@[1,2])/3)

def test_acquisition_not_sequence():
    r=t.acquisition_order({'a':2,'b':1},{'a':0,'b':0},2,1)
    assert len(r['orders'])==2 and r['objective']==3
    r=t.acquisition_order({'a':2,'b':1},{'a':0,'b':0},2,.5)
    assert r['orders']==[['a','b']] and r['objective']==2.5

def test_bounded_gaussian():
    assert t.bounded_normal_mean(0,1,-1,1)==pytest.approx(0)
    np.testing.assert_allclose(t.bounded_normal_mean([-2,0,2],[0,0,0],-1,1),[-1,0,1])
    rng=np.random.default_rng(101)
    expected=np.clip(rng.normal(.6,.8,500000),0,1).mean()
    assert float(t.bounded_normal_mean(.6,.8,0,1))==pytest.approx(expected,abs=.002)

def test_monte_carlo_mean(model):
    rng=np.random.default_rng(17)
    x=rng.multivariate_normal(model.mu1,model.sigma1,100000)
    xp=x.copy();xp[:,0]=model.anchors[0]
    # Paired residuals cancel. This tests orientation and operator, not confidence intervals.
    empirical=((x-xp)@model.B.T).mean(axis=0)
    np.testing.assert_allclose(empirical,model.state({'S01':1})['response'],atol=.008)

@pytest.mark.parametrize('field,value',[
    ('orientation','source_by_outcome'),('higher_is_worse',False),
    ('outcome_labels',['S02','S01']+[f'S{i:02d}' for i in range(3,23)]),('anchors',[0,0]),
    ('B',[[1]]),('mu1',[float('nan')]+[1]*21),
    ('sigma1',np.diag([-1]+[1]*21).tolist()),('modules',['A']*22),
])
def test_invalid_models(data,field,value):
    data['model'][field]=value
    with pytest.raises(ValueError):t.Model(**data['model'])

@pytest.mark.parametrize('dose',[-.1,1.1,float('inf')])
def test_invalid_doses(model,dose):
    with pytest.raises(ValueError):model.state({'S01':dose})

def test_invalid_score_weights(model,data):
    data['scoring']['outcome_weights']=[1]+[0]*21
    with pytest.raises(ValueError,match='positive weight'):t.score(model,**data['scoring'])

def test_nonfinite_derived_moments(data):
    data['model']['B']=(np.eye(22)*1e200).tolist()
    with pytest.raises((ValueError,FloatingPointError)):
        t.Model(**data['model'])

def test_large_finite_utility_weights(model,data):
    base=t.score(model,**data['scoring'])
    data['scoring']['utility_weights']=[1e306]*7
    got=t.score(model,**data['scoring'])
    np.testing.assert_allclose(got['tvpps'],base['tvpps'])
