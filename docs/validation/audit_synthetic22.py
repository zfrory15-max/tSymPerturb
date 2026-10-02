"""Independent equations on the public constructed 22-symptom example.

This validates implementation, not unavailable original manuscript results.
"""
import json
from pathlib import Path
import sys
import numpy as np
from numpy.testing import assert_allclose
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'skills/tsymperturb/scripts'))
import tsymperturb as t


def audit():
    data=json.loads((ROOT/'examples/synthetic_model.json').read_text())
    spec=data['model']; config=data['scoring']; m=t.Model(**spec)
    b=np.asarray(spec['B']); mu=np.asarray(spec['mu1']); anchor=np.asarray(spec['anchors'])
    sigma=np.asarray(spec['sigma1']); psi=np.asarray(spec['psi']); w=np.asarray(config['outcome_weights'])
    modules=np.asarray(spec['modules']); names=spec['source_labels']; p=len(names)
    sd=np.sqrt(np.diag(b@sigma@b.T+psi)); response=b*(mu-anchor)[None,:]; z=response/sd[:,None]
    def capacity(mat):
        return sum(np.linalg.matrix_power(config['gamma']*abs(mat),h).sum() for h in range(1,config['horizon']+1))
    raw=[]; max_pair=0.; max_dose=0.; count=0
    for i in range(p):
        other=np.arange(p)!=i; cb=b.copy(); cb[:,i]*=1-config['q']; increments=[]
        for k in range(p):
            if k==i:continue
            rest=(np.arange(p)!=i)&(np.arange(p)!=k)
            gi=np.average(z[rest,i],weights=w[rest]); gk=np.average(z[rest,k],weights=w[rest])
            joint=np.average(z[rest,i]+z[rest,k],weights=w[rest]); increments.append(max(0,joint-max(gi,gk)))
            if k>i:
                expected=response[:,i]+response[:,k]
                max_pair=max(max_pair,float(np.max(abs(m.state({names[i]:1,names[k]:1})['response']-expected))))
                count+=1
        pos=np.maximum(z[:,i],0)*w
        raw.append([np.average(z[:,i],weights=w),np.average(z[other,i],weights=w[other]),
                    np.mean(z[other,i]>=config['breadth_threshold']),
                    np.mean([np.mean(z[modules==g,i])>=config['module_threshold'] for g in set(modules) if g!=modules[i]]),
                    (capacity(b)-capacity(cb))/capacity(b),np.mean(increments),pos[other].sum()/(pos.sum()+1e-12)])
        for dose in data['provenance']['dose_grid']:
            state=m.state({names[i]:dose}); f=np.ones(p);f[i]=1-dose
            max_dose=max(max_dose,float(np.max(abs(state['response']-dose*response[:,i]))))
            expected_cov=b@(f[:,None]*sigma*f[None,:])@b.T+psi
            assert_allclose(state['outcome_covariance'],expected_cov,atol=1e-12)
    raw=np.asarray(raw); spread=np.ptp(raw,axis=0); norm=np.full_like(raw,50.)
    for j in range(7):
        if spread[j]!=0:norm[:,j]=100*(raw[:,j]-raw[:,j].min())/spread[j]
    result=t.score(m,**config)
    assert_allclose(result['raw'],raw,atol=1e-12)
    assert_allclose(result['normalized'],norm,atol=1e-12)
    assert_allclose(result['tvpps'],np.average(norm,axis=1,weights=config['utility_weights']),atol=1e-12)
    assert count==231 and max_pair<1e-12 and max_dose<1e-12
    return dict(status='PASS',scope='SYNTHETIC ILLUSTRATIVE ONLY; original replication NOT_RUN',
                symptoms=p,modules=len(set(modules)),utilities=7,pairs=count,dose_checks=p*6,
                max_response_pair_error=max_pair,max_response_dose_error=max_dose,
                spectral_radius=float(max(abs(np.linalg.eigvals(b)))),
                top5=[names[i] for i in np.argsort(-result['tvpps'],kind='stable')[:5]])

if __name__=='__main__':print(json.dumps(audit(),indent=2))
