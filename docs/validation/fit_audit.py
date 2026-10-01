"""Independent end-to-end check of optional fitting/bootstrap adapter."""
from pathlib import Path
import json,sys,tempfile,subprocess
import numpy as np
from numpy.testing import assert_allclose as eq
from sklearn.linear_model import Lasso
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'skills'/'tsymperturb'/'scripts'))
import fit_tsymperturb as f
rng=np.random.default_rng(20260930)
n=70;p=4
x=rng.normal(size=(n,p))*[.5,2.,1.,3.]+[2.,3.,4.,5.]
B=np.array([[.5,.1,0,0],[.2,.4,.1,0],[0,.1,.5,.1],[.1,0,.2,.4]])
y=x@B.T+rng.normal(size=(n,p))*.2+[2,0,-1,1]
names=['a','b','c','d']
fitting=dict(labels=names,modules=['A','A','B','B'],raw_anchors=[0,.5,1,0],higher_is_worse=True)
scoring=dict(candidates=names,partners={name:[o for o in names if o!=name] for name in names},outcome_weights=[1,2,3,4],utility_weights=[1]*7,horizon=3,gamma=.7)
r=f.bootstrap_pairs(x,y,fitting=fitting,scoring=scoring,n_boot=8,seed=91)
assert r['succeeded']==8 and r['failed']==0
bs=np.random.default_rng(91)
for replicate in r['replicates']:
 ids=bs.integers(0,n,n);xx=x[ids];yy=y[ids]
 zx=(xx-xx.mean(0))/xx.std(0);zy=(yy-yy.mean(0))/yy.std(0)
 bb=np.array([Lasso(alpha=.03,fit_intercept=False,max_iter=10000,tol=1e-8,selection='cyclic').fit(zx,zy[:,j]).coef_ for j in range(p)])
 eq(replicate['fit']['B'],bb,atol=1e-13)
 eq(replicate['fit']['transformed_anchors'],(np.array(fitting['raw_anchors'])-xx.mean(0))/xx.std(0),atol=1e-13)
 residual=zy-zx@bb.T;residual-=residual.mean(0)
 cov=bb@(zx.T@zx/n)@bb.T+residual.T@residual/n
 eq(replicate['fit']['normalization_sd'],np.sqrt(np.diag(cov)),atol=1e-13)
 raw=np.asarray(replicate['analysis']['raw']);sp=np.ptp(raw,axis=0)
 norm=np.empty_like(raw)
 for j in range(7):norm[:,j]=50 if sp[j]==0 else 100*(raw[:,j]-raw[:,j].min())/sp[j]
 eq(replicate['analysis']['normalized'],norm,atol=1e-12);eq(replicate['analysis']['tvpps'],norm.mean(1),atol=1e-12)
for i,s in enumerate(r['summary']):
 ranks=np.array([rep['analysis']['ranks'][i] for rep in r['replicates']])
 for k,item in s['top_k'].items():
  eq(item['probability'],np.mean(ranks<=int(k)))
  lo,hi=item['monte_carlo_wilson_interval'];assert 0<=lo<=item['probability']<=hi<=1
work=Path(tempfile.mkdtemp(prefix='tsymperturb-fit-audit-'))
(work/'input.json').write_text(json.dumps(dict(t1=x.tolist(),t2=y.tolist(),fitting=fitting,scoring=scoring)))
res=subprocess.run([sys.executable,f.__file__,str(work/'input.json'),'--output',str(work/'out.json'),'--n-boot','8','--seed','91'],cwd=work,capture_output=True,text=True)
assert res.returncode==0,res.stderr
cli=json.loads((work/'out.json').read_text());eq(cli['point_estimate']['tvpps'],r['point_estimate']['tvpps']);assert cli['succeeded']==8
print('PASS independent adapter: matched resampling, replicate standardization/anchor transform, direct Lasso coefficients, structural covariance convention, rerun normalization, rank probabilities, Wilson endpoints, isolated CLI (8/8 bootstrap successes)')
# Real degenerate resamples (not mocked) must be logged without replacement.
x3=np.array([[1.,2.,3.],[2.,4.,6.],[3.,6.,9.]])
n3=['a','b','c'];fit3=dict(labels=n3,modules=['A','B','B'],raw_anchors=[0,0,0],higher_is_worse=True)
sc3=dict(candidates=n3,partners={n:[o for o in n3 if o!=n] for n in n3},outcome_weights=[1]*3,utility_weights=[1]*7,horizon=2,gamma=.7)
bad=f.bootstrap_pairs(x3,x3,fitting=fit3,scoring=sc3,n_boot=30,seed=7)
assert bad['failed']>0 and bad['succeeded']>0 and bad['failed']+bad['succeeded']==30
assert bad['status']=='partial_failures'
assert len(bad['failures'])==bad['failed']
print('PASS real degenerate resamples logged:',bad['succeeded'],'successful,',bad['failed'],'failed; exactly 30 attempts')
