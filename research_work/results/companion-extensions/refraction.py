"""Axisymmetric positive-permittivity pilot, exact-third source retained."""
from common import *
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.interpolate import RegularGridInterpolator
import hashlib

class Grid:
 def __init__(self,gal,n=128,L=160.):
  self.gal=gal;self.n=n;self.L=L
  self.rf=L*np.linspace(0,1,n+1)**1.4;self.zf=.02*np.expm1(np.linspace(0,np.log1p(L/.02),n+1))
  self.r=(self.rf[:-1]+self.rf[1:])/2;self.z=(self.zf[:-1]+self.zf[1:])/2
  R,Z=np.meshgrid(self.r,self.z,indexing='ij');self.phi,_,_,self.rb,self.rd=gal.field(R,Z)
  ids=np.arange(n*n).reshape(n,n);dr=np.diff(self.rf);dz=np.diff(self.zf);ar=.5*np.diff(self.rf**2)
  self.ii=np.r_[ids[:-1,:].ravel(),ids[:,:-1].ravel()];self.jj=np.r_[ids[1:,:].ravel(),ids[:,1:].ravel()]
  self.area=np.r_[(self.rf[1:-1,None]*dz).ravel(),np.broadcast_to(ar[:,None],(n,n-1)).ravel()]
  self.hi=np.r_[np.broadcast_to(dr[:-1,None]/2,(n-1,n)).ravel(),np.broadcast_to(dz[None,:-1]/2,(n,n-1)).ravel()]
  self.hj=np.r_[np.broadcast_to(dr[1:,None]/2,(n-1,n)).ravel(),np.broadcast_to(dz[None,1:]/2,(n,n-1)).ravel()]
  self.c0=self.area/(self.hi+self.hj)
  self.bi=np.r_[ids[-1,:],ids[:,-1]]
  self.ba=np.r_[self.rf[-1]*dz,ar];self.bh=np.r_[np.full(n,dr[-1]/2),np.full(n,dz[-1]/2)]
  self.bp=np.r_[gal.field(np.full(n,L),self.z)[0],gal.field(self.r,np.full(n,L))[0]]
  self.brb=np.r_[gal.field(np.full(n,L),self.z)[3],gal.field(self.r,np.full(n,L))[3]]
  self.brd=np.r_[gal.field(np.full(n,L),self.z)[4],gal.field(self.r,np.full(n,L))[4]]
 def solve(self,einf,crit,driver,constant=False):
  density=self.rb if driver=='ordinary' else self.rd;bd=self.brb if driver=='ordinary' else self.brd
  eps=np.full(self.n*self.n,einf) if constant else (einf+(1-einf)*density/(density+crit)).ravel()
  be=np.full(len(self.bi),einf) if constant else einf+(1-einf)*bd/(bd+crit)
  c=self.area/(self.hi/eps[self.ii]+self.hj/eps[self.jj]);bc=self.ba/self.bh*be;bc0=self.ba/self.bh
  diag=np.bincount(np.r_[self.ii,self.jj,self.bi],weights=np.r_[c,c,bc],minlength=len(eps))
  A=coo_matrix((np.r_[diag,-c,-c],(np.r_[np.arange(len(eps)),self.ii,self.jj],np.r_[np.arange(len(eps)),self.jj,self.ii])),shape=(len(eps),len(eps))).tocsr()
  ph=self.phi.ravel();flux=(c-self.c0)*(ph[self.jj]-ph[self.ii]);bflux=(bc-bc0)*(self.bp-ph[self.bi])+bc*(1/einf-1)*self.bp
  rhs=np.bincount(np.r_[self.ii,self.jj,self.bi],weights=np.r_[flux,-flux,bflux],minlength=len(eps))
  delta=spsolve(A,rhs);res=float(np.linalg.norm(A@delta-rhs)/max(np.linalg.norm(rhs),1.))
  assert res<1e-8
  if constant:assert np.max(np.abs(delta-(1/einf-1)*ph))/np.max(np.abs(ph))<1e-9
  d=delta.reshape(self.n,self.n);dR,dZ=np.gradient(d,self.r,self.z,edge_order=2)
  ir=RegularGridInterpolator((self.r,self.z),dR,bounds_error=False,fill_value=None);iz=RegularGridInterpolator((self.r,self.z),dZ,bounds_error=False,fill_value=None)
  corr=ir(np.c_[Robs,np.full(len(Robs),self.z[0])]);g0=self.gal.field(Robs,0.)[1];v=np.sqrt(np.maximum(0,Robs*(g0+corr)))
  gz=self.gal.field(Rz,1.1)[2]+iz(np.c_[Rz,np.full(len(Rz),1.1)]);K=gz/(2*np.pi*G*1e6)
  return dict(epsilon_infinity=einf,rho_critical=crit,driver=driver,scores=scores(v),predicted_kms=v.tolist(),vertical_surface_equivalent=K.tolist(),vertical_provisional_RMS=float(np.sqrt(np.mean((K-Kobs)**2))),linear_residual=res,epsilon_min=float(eps.min()),epsilon_max=float(eps.max()))

runs=[];refinement=[]
for variant in ['I','II']:
 gal=Galaxy(variant);grid=Grid(gal)
 grid.solve(.8,1e6,'ordinary',True)
 control=grid.solve(1.,1e6,'ordinary');expected=np.sqrt(Robs*gal.field(Robs,0.)[1]);assert np.max(np.abs(expected-control['predicted_kms']))<1e-10
 local=[]
 for driver in ['ordinary','companion']:
  for e in [.7,.85,1.,1.15,1.3]:
   for rho in [1e5,1e6,1e7]:
    r=grid.solve(e,rho,driver);r['baryons']=variant;local.append(r)
  best=min([r for r in local if r['driver']==driver],key=lambda r:r['scores']['inner']['RMSE_kms'])
  for n,L,label in [(192,160.,'grid'),(192,240.,'domain')]:
   refined=Grid(gal,n,L).solve(best['epsilon_infinity'],best['rho_critical'],driver)
   refinement.append(dict(baryons=variant,driver=driver,kind=label,parameters=[best['epsilon_infinity'],best['rho_critical']],max_speed_change_kms=float(np.max(np.abs(np.array(best['predicted_kms'])-refined['predicted_kms']))),vertical_max_change=float(np.max(np.abs(np.array(best['vertical_surface_equivalent'])-refined['vertical_surface_equivalent']))),refined=refined))
  print(variant,driver,'control',control['scores'],'best',best['epsilon_infinity'],best['rho_critical'],best['scores'],flush=True)
 runs.extend(local)
result=dict(scope='Static axisymmetric response pilot; gas monopole, fixed source, no relativistic lens completion',input_sha256={str(p.relative_to(P.parents[2])):hashlib.sha256(p.read_bytes()).hexdigest() for p in [BASE/'third-radiation-retention-results.json',P.parent/'milky-way-capture/inputs.json',P/'common.py',P/'protocol.md']},runs=runs,refinement=refinement)
(P/'refraction-results.json').write_text(json.dumps(result,separators=(',',':'),allow_nan=False)+'\n',encoding='utf-8',newline='\n')
