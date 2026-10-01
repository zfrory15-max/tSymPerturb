import numpy as np
import pytest
pytest.importorskip('sklearn')
import fit_tsymperturb as adapter
from tsymperturb import serializable


def data():
    rng = np.random.default_rng(831)
    x = rng.normal(size=(60, 4)) + [2,3,4,5]
    b = np.array([[.6,.2,0,0],[0,.7,.3,0],[.1,0,.6,.2],[.2,0,0,.8]])
    y = x@b.T+rng.normal(scale=.25, size=x.shape)
    labels = ['A','B','C','D']
    fitting = dict(labels=labels, modules=['M1','M1','M2','M2'],
                   raw_anchors=[0]*4, higher_is_worse=True)
    scoring = dict(candidates=labels, partners={s:[t for t in labels if t!=s] for s in labels},
                   outcome_weights=[1]*4, utility_weights=[1]*7, horizon=3,gamma=.7)
    return x,y,fitting,scoring


def test_fit_anchor_scale_orientation_and_residuals():
    x,y,f,s = data()
    model,d = adapter.fit_pairs(x,y,**f)
    np.testing.assert_allclose(model.anchors,-x.mean(0)/x.std(0))
    z1 = (x-x.mean(0))/x.std(0)
    z2 = (y-y.mean(0))/y.std(0)
    from sklearn.linear_model import Lasso
    expected = np.stack([Lasso(alpha=.03,fit_intercept=False,tol=1e-8,max_iter=10000).fit(z1,z2[:,j]).coef_ for j in range(4)])
    np.testing.assert_allclose(model.B,expected)
    residual = z2-z1@model.B.T
    np.testing.assert_allclose(model.psi, residual.T@residual/len(x), atol=1e-14)
    np.testing.assert_allclose(model.sigma1,z1.T@z1/len(x))


def test_reproducible_complete_pipeline_and_paired_rows(monkeypatch):
    x,y,f,s=data()
    original=adapter.fit_pairs
    calls=[]
    def spy(a,b,**kw):
        calls.append((a.copy(),b.copy()))
        return original(a,b,**kw)
    monkeypatch.setattr(adapter,'fit_pairs',spy)
    result=adapter.bootstrap_pairs(x,y,fitting=f,scoring=s,n_boot=5,seed=46)
    assert result['attempted']==result['succeeded']==5
    rng=np.random.default_rng(46)
    for b in range(5):
        ix=rng.integers(0,len(x),size=len(x))
        np.testing.assert_array_equal(calls[b+1][0],x[ix])
        np.testing.assert_array_equal(calls[b+1][1],y[ix])
        rep=result['replicates'][b]
        m,d=original(x[ix],y[ix],**f)
        expected=adapter.score(m,**s)
        for key in ('raw','normalized','tvpps','ranks'):
            np.testing.assert_allclose(rep['analysis'][key],expected[key])
        np.testing.assert_allclose(rep['fit']['transformed_anchors'],-x[ix].mean(0)/x[ix].std(0))
    second=adapter.bootstrap_pairs(x,y,fitting=f,scoring=s,n_boot=5,seed=46)
    assert serializable(result)==serializable(second)
    for row in result['summary']:
        assert row['top_k']['5']['probability']==1
        for item in row['top_k'].values():
            lo,hi=item['monte_carlo_wilson_interval']
            assert 0<=lo<=item['probability']<=hi<=1


def test_failures_counted_without_replacement(monkeypatch):
    x,y,f,s=data()
    original=adapter.fit_pairs
    count=0
    def fail_after_point(*args,**kw):
        nonlocal count
        count+=1
        if count>1: raise ValueError('Synthetic replicate failure')
        return original(*args,**kw)
    monkeypatch.setattr(adapter,'fit_pairs',fail_after_point)
    result=adapter.bootstrap_pairs(x,y,fitting=f,scoring=s,n_boot=3,seed=4)
    assert result['attempted']==result['failed']==3
    assert result['succeeded']==0 and result['summary']==[]
    assert result['status']=='no_successful_replicates'
    assert result['failure_reason_counts']=={'ValueError: Synthetic replicate failure':3}


def test_missing_constant_and_convergence_rejected():
    x,y,f,s=data()
    x[0,0]=np.nan
    with pytest.raises(ValueError,match='missing'): adapter.fit_pairs(x,y,**f)
    x,y,f,s=data();x[:,0]=1
    with pytest.raises(ValueError,match='Constant'): adapter.fit_pairs(x,y,**f)
    x,y,f,s=data()
    from sklearn.exceptions import ConvergenceWarning
    with pytest.raises(ConvergenceWarning): adapter.fit_pairs(x,y,**f,max_iter=1,tol=1e-16)


def test_top_k_all_exact_ties_included(monkeypatch):
    x,y,f,s=data()
    original=adapter.score
    def tied(*args,**kw):
        result=original(*args,**kw)
        result['tvpps']=np.full(4,50.)
        result['ranks']=[1]*4
        return result
    monkeypatch.setattr(adapter,'score',tied)
    result=adapter.bootstrap_pairs(x,y,fitting=f,scoring=s,n_boot=3,seed=2)
    assert sum(row['top_k']['1']['probability'] for row in result['summary'])==4
