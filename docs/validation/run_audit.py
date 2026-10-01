import sys,json,subprocess,copy,warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'skills'/'tsymperturb'/'scripts'))
import tsymperturb as t
from oracle import *
from numpy.testing import assert_allclose as eq
names=['a','b','c','d']
spec=dict(source_labels=names,outcome_labels=names,B=B.tolist(),mu1=mu.tolist(),sigma1=Sigma.tolist(),psi=Psi.tolist(),intercept=[.2,-.1,.3,.4],anchors=anchor.tolist(),modules=modules.tolist(),orientation='outcome_by_source',higher_is_worse=True,scale='raw')
sc=dict(candidates=names,partners={n:[m for m in names if n!=m] for n in names},outcome_weights=w.tolist(),utility_weights=[1]*7,horizon=3,gamma=.7)
m=t.Model(**spec);out=t.score(m,**sc)
eq(out['raw'],raw,atol=1e-13);eq(out['normalized'],scores,atol=1e-12)
for i,n in enumerate(names):
 for dose in [0,.001,.2,.7,1]:
  r=m.state({n:dose});eq(r['response'],dose*R[:,i],atol=1e-14)
  f=np.ones(4);f[i]=1-dose
  eq(r['source_covariance'],np.diag(f)@Sigma@np.diag(f))
  eq(r['outcome_covariance'],B@np.diag(f)@Sigma@np.diag(f)@B.T+Psi)
  eq(r['standardized_response'],dose*Z[:,i],atol=1e-14)
  loc=m.state({n:dose},{n:1});eq(loc['response'],r['response']);eq(loc['source_covariance'],Sigma)
for i in range(4):
 for k in range(i+1,4): eq(m.state({names[i]:1,names[k]:1})['response'],R[:,i]+R[:,k])
eq(t.transform_anchors(anchor,mu,np.sqrt(np.diag(Sigma))),(anchor-mu)/np.sqrt(np.diag(Sigma)))
eq(t.block(m,'a',.8)[:,0],.2*B[:,0]);edge=t.block(m,'a',.8,outcome='c');target=B.copy();target[2,0]*=.2;eq(edge,target)
cross=t.block(m,'a',1,cross_only=True);eq(cross[0,0],B[0,0]);eq(cross[1:,0],0)
eq(t.reference_block_response(m,edge,mu),(B-edge)@mu)
eq(t.reference_block_response(m,edge,[0]*4),[0]*4)
b1=np.array([[1.,2],[0,1]]);b2=np.array([[1.,0],[3,1]]);d=np.array([1.,2.]);u=np.array([[.2,.3],[-.1,.4]])
eq(t.trajectories([b1,b2],d),[b1@d,b2@b1@d]);eq(t.trajectories([b1,b2],d,u),[b1@(d+u[0]),b2@(b1@(d+u[0])+u[1])])
assert t.acquisition_order({'a':2,'b':1},{'a':0,'b':0},2,.5)['orders']==[['a','b']]
assert len(t.acquisition_order({'a':2,'b':1},{'a':0,'b':0},2,1)['orders'])==2
# Numerical quadrature independent from implemented closed form
from scipy.integrate import quad
from scipy.stats import norm
for mm,ss,ll,uu in [(0,1,-1,2),(4,.3,0,3),(-2,.8,0,3),(1,0,0,3)]:
 expected=np.clip(mm,ll,uu) if ss==0 else ll*norm.cdf((ll-mm)/ss)+quad(lambda x:x*norm.pdf(x,loc=mm,scale=ss),ll,uu)[0]+uu*norm.sf((uu-mm)/ss)
 eq(t.bounded_normal_mean(mm,ss,ll,uu),expected,atol=1e-12)
# candidate-set normalization and singleton constants
one=t.score(m,**dict(sc,candidates=['a'],partners={'a':['b']}));eq(one['normalized'],np.full((1,7),50));eq(one['tvpps'],[50])
# Bad inputs must reject rather than silently compute
bad=[]
def reject(label,call):
 try:call()
 except (ValueError,TypeError,KeyError,FloatingPointError):return
 bad.append(label)
reject('orientation',lambda:t.Model(**dict(spec,orientation='source_by_outcome')))
reject('labels mismatch',lambda:t.Model(**dict(spec,outcome_labels=names[::-1])))
reject('non PSD',lambda:t.Model(**dict(spec,sigma1=-Sigma)))
reject('dose',lambda:m.state({'a':1.01}))
reject('partners same',lambda:t.score(m,**dict(sc,partners={**sc['partners'],'a':['a']})))
reject('weights zero comparison',lambda:t.score(m,**dict(sc,outcome_weights=[1,0,0,0])))
reject('Q zero',lambda:t.score(t.Model(**dict(spec,B=np.zeros((4,4)))),**sc))
reject('gamma endpoint',lambda:t.score(m,**dict(sc,gamma=1)))
reject('noninteger H',lambda:t.score(m,**dict(sc,horizon=2.5)))
# Isolated CLI forward task, no repository cwd dependencies
import tempfile
work=Path(tempfile.mkdtemp(prefix='tsymperturb-audit-'));inp=str(work/'input.json');dst=str(work/'output.json')
open(inp,'w').write(json.dumps(dict(model=spec,scoring=sc,dose_requests=[dict(doses={'a':.5})])))
r=subprocess.run([sys.executable,t.__file__,inp,'--output',dst],cwd='/tmp',capture_output=True,text=True);assert r.returncode==0,r.stderr
cli=json.load(open(dst));eq(cli['raw'],raw);eq(cli['dose_results'][0]['result']['response'],R[:,0]*.5)
print('PASS independent oracle: seven utilities, normalization, orientation, anchor transform, covariance, dose linearity, pair additivity, blocks, noncommuting transitions, interventions, acquisition, clipped Gaussian quadrature, isolated CLI')
print('INPUT FAILURES',bad);assert not bad
# Large finite inputs must not yield nonfinite public API result
with warnings.catch_warnings():
 warnings.simplefilter('ignore')
 try:
  giant=t.Model(**dict(spec,B=B*1e200))
  raise AssertionError('Overflow model was silently accepted')
 except (ValueError,FloatingPointError) as e:print('PASS overflow model rejected:',type(e).__name__)
# Unequal module sizes and threshold equality (zero source covariance; residual SD=1)
n5=['a','b','c','d','e'];bb=np.zeros((5,5));bb[:,0]=[.2,.03,.03,.03,.03];bb[:,1]=[.05,.2,.06,-.03,.06];mm=t.Model(source_labels=n5,outcome_labels=n5,B=bb,mu1=np.ones(5),sigma1=np.zeros((5,5)),psi=np.eye(5),intercept=np.zeros(5),anchors=np.zeros(5),modules=['A','B','C','C','C'],orientation='outcome_by_source',higher_is_worse=True,scale='raw')
rr=t.score(mm,candidates=['a','b'],partners={'a':['b'],'b':['a']},outcome_weights=[1,2,3,4,5],utility_weights=[1]*7,horizon=2,gamma=.5)
eq(rr['raw'][0,3],1);eq(rr['raw'][1,2],.75);eq(rr['raw'][1,3],1)
print('PASS unequal module sizes, unweighted module means, inclusive breadth/reach thresholds')
reject('outcome-weight overflow',lambda:t.score(m,**dict(sc,outcome_weights=[1e308]*4)))
reject('utility-weight overflow',lambda:t.score(m,**dict(sc,utility_weights=[1e308]*7)))
assert not bad,bad
print('PASS weight-sum overflow rejected')
