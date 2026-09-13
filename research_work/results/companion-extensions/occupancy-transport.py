"""Spherical empty-site capture with an escaping release channel; stationary diagnostic."""
from pathlib import Path
import hashlib, io, json, zipfile
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import coo_matrix
from scipy.optimize import minimize_scalar

P=Path(__file__).resolve().parent; OLD=P.parent/'isotropic-galaxy-transfer'
BASE=P.parents[2]/'temporal_candidate_audit/data'
fitfile=OLD/'third-radiation-retention-results.json'; predfile=OLD/'third-radiation-retention-predictions.json'
localfile=P/'local-capacity-results.json'
reference=json.loads(fitfile.read_text()); model=reference['models']['attenuated']
saved=json.loads(predfile.read_text()); local={r['galaxy']:r for r in json.loads(localfile.read_text())['rows']}
for name,digest in reference['input_sha256'].items():
    assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==digest
def eta(x):
    q=np.cbrt(x); return q/(1+q)

class Rays:
    def __init__(self,nr,nmu,nline):
        self.x=np.r_[0.,np.geomspace(1e-5,1e3,nr-1)]
        self.z=np.log1p(self.x)
        mu,self.mw=leggauss(nmu); self.mu=mu; nodes,weights=leggauss(nline)
        x=self.x[:,None]; t=x*mu; b2=x*x*(1-mu*mu); B=np.sqrt(1+b2)
        upper=np.arctan(t/B); width=upper+np.pi/2
        alpha=-np.pi/2+width[:,:,None]*(nodes+1)/2
        radius=np.sqrt(b2[:,:,None]+B[:,:,None]**2*np.tan(alpha)**2)
        factor=width[:,:,None]/2*weights*np.cos(alpha)**2/B[:,:,None]**3
        self.A=self.operator(radius,factor)
        self.shape=(len(self.x),nmu)
        exact=t/(2*(1+b2)*(1+b2+t*t))+(np.arctan(t/B)+np.pi/2)/(2*B**3)
        assert np.max(abs((self.A@np.ones(len(self.x))).reshape(self.shape)-exact))<1e-12
        assert np.max(abs(self.field(np.zeros(len(self.x)),1)-1))<1e-12
        # Complete chords, for a separate global absorbed-power integration.
        nenergy=256
        en,ew=leggauss(nenergy); self.phi=(en+1)*np.pi/4; self.ew=ew*np.pi/4
        b=np.tan(self.phi); BB=np.sqrt(1+b*b)
        alpha=nodes*np.pi/2
        radius=np.sqrt(b[:,None]**2+BB[:,None]**2*np.tan(alpha)**2)
        factor=np.pi/2*weights*np.cos(alpha)**2/BB[:,None]**3
        self.chord=self.operator(radius,factor)
        self.b=b
    def operator(self,radius,factor):
        zz=np.log1p(radius)
        idx=np.clip(np.searchsorted(self.z,zz)-1,0,len(self.z)-2)
        frac=np.clip((zz-self.z[idx])/(self.z[idx+1]-self.z[idx]),0,1)
        rr=np.broadcast_to(np.arange(radius.size//radius.shape[-1]).reshape(radius.shape[:-1]+(1,)),radius.shape)
        return coo_matrix((np.r_[(factor*(1-frac)).ravel(),(factor*frac).ravel()],
            (np.r_[rr.ravel(),rr.ravel()],np.r_[idx.ravel(),(idx+1).ravel()])),
            shape=(radius.size//radius.shape[-1],len(self.z))).tocsr()
    def field(self,vacancy,ka):
        tau=ka*(self.A@vacancy).reshape(self.shape)
        return np.exp(-tau)@self.mw/2
    def solve(self,X,ka):
        low=self.field(np.ones(len(self.x)),ka); high=np.ones(len(self.x))
        for iteration in range(200):
            newlow=self.field(1-eta(X*low),ka); newhigh=self.field(1-eta(X*high),ka)
            assert min(newlow-low)>-1e-11 and max(newhigh-high)<1e-11
            low,high=newlow,newhigh
            if max(high-low)<1e-9:break
        else:raise AssertionError('Stationary brackets did not close')
        J=(low+high)/2; f=eta(X*J)
        # Total power / (c u_external a^2), via chords versus volume absorption.
        cross=2*np.pi*np.sum(self.ew*self.b/(np.cos(self.phi)**2)*(-np.expm1(-ka*(self.chord@(1-f)))))
        sample=np.interp(np.log1p(np.tan(self.phi)),self.z,J)
        volume=4*np.pi*ka*np.sum(self.ew*np.sin(self.phi)**2*(1-eta(X*sample))*sample)
        return J,dict(iterations=iteration+1,bracket_width=float(max(high-low)),
            absorbed_power_dimensionless=float(cross),volume_absorbed_power_dimensionless=float(volume),
            global_power_relative_difference=float(abs(volume/cross-1)))
    def masses(self,radii,a,X,J):
        nodes,weights=leggauss(128)
        x=radii[:,None]/a*(nodes+1)/2
        jj=np.interp(np.log1p(x),self.z,J)
        # Analytic old full-opacity field, same angular quadrature as this grid.
        mu=self.mu
        xx=x[:,:,None];t=xx*mu;b2=1+xx*xx*(1-mu*mu);b=np.sqrt(b2)
        ka=model['k0_per_kpc']*a
        tau=ka*(t/(2*b2*(b2+t*t))+(np.arctan(t/b)+np.pi/2)/(2*b**3))
        oldJ=np.exp(-tau)@self.mw/2
        radial=x*x/(1+x*x)**2
        new=np.sum(radial*eta(X*jj)*weights,axis=1)
        old=np.sum(radial*oldJ*eta(X)*weights,axis=1)
        return new/old

out=dict(scope='Local occupancy feeds back on spherical capture opacity; all release assigned to non-recaptured escaping energy. Static endpoint, not dynamics or a microscopic channel.',
    input_sha256=reference['input_sha256'],model_input_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [fitfile,predfile,localfile]},rows=[],scores=[])
grids=[Rays(193,40,80),Rays(385,80,160)]
refined_grid=None
cache=[]
with zipfile.ZipFile(BASE/'Rotmod_LTG.zip') as archive:
    for count,row in enumerate(saved):
        name=row['galaxy'];prior=local[name];X=prior['X'];a=prior['a_kpc']
        arr=np.atleast_2d(np.loadtxt(io.BytesIO(archive.read(name+'_rotmod.dat'))))
        vb=arr[:,3]*abs(arr[:,3])+.5*arr[:,4]*abs(arr[:,4])+.7*arr[:,5]*abs(arr[:,5])
        good=np.isfinite(arr).all(axis=1)&(arr[:,0]>0)&(arr[:,1]>0)&(arr[:,2]>0)&(vb>0)
        arr,vb=arr[good],vb[good]
        assert arr[:,0].tolist()==row['R_kpc'] and arr[:,1].tolist()==row['observed_kms']
        extra=np.array(row['predicted_kms'])**2-vb
        predictions=[]
        for grid in grids:
            J,checks=grid.solve(X,model['k0_per_kpc']*a)
            ratio=grid.masses(arr[:,0],a,X,J)
            predictions.append(np.sqrt(vb+extra*ratio))
        v=predictions[-1];error=float(max(abs(v-predictions[0])))
        final_grid=385
        if error>.05:
            if refined_grid is None:refined_grid=Rays(769,80,160)
            J,checks=refined_grid.solve(X,model['k0_per_kpc']*a)
            ratio=refined_grid.masses(arr[:,0],a,X,J)
            newv=np.sqrt(vb+extra*ratio);error=float(max(abs(newv-v)));v=newv;final_grid=769
        assert error<.1,(name,error)
        assert checks['global_power_relative_difference']<.005,(name,checks)
        assert min(v-np.array(prior['predicted_kms']))>-.02
        cache.append(dict(split=row['split'],vb=vb,extra=extra,ratio=ratio,obs=arr[:,1]))
        out['rows'].append(dict(galaxy=name,split=row['split'],X=X,a_kpc=a,
            R_kpc=arr[:,0].tolist(),predicted_kms=v.tolist(),extra_mass_ratio=ratio.tolist(),
            central_field=float(J[0]),central_occupancy=float(eta(X*J[0])),
            grid_speed_difference_kms=error,final_radial_grid_size=final_grid,**checks))
        if count%30==0:print('Solved',count+1,flush=True)
for branch in ['reference','feedback']:
    def loss(amplitude,split):
        return float(np.mean([np.mean((np.sqrt(r['vb']+amplitude*r['extra']*(1 if branch=='reference' else r['ratio']))-r['obs'])**2) for r in cache if r['split']==split]))
    fit=minimize_scalar(lambda a:loss(a,'train'),bounds=(0,4),method='bounded',options={'xatol':1e-10})
    assert fit.success and 1e-7<fit.x<4-1e-7
    for mode,amp in [('frozen',1),('training_C_refit',float(fit.x))]:
        scores={s:float(np.sqrt(loss(amp,s))) for s in ['train','validation','test']}
        if branch=='reference' and mode=='frozen':
            for s,value in scores.items():
                assert min(abs(value-model[k][s]['RMSE_kms']) for k in ['scores','finer_scores'])<1e-8
        out['scores'].append(dict(branch=branch,mode=mode,C_multiplier=amp,RMSE_kms=scores))
out['maximum_grid_speed_difference_kms']=max(r['grid_speed_difference_kms'] for r in out['rows'])
out['maximum_global_power_relative_difference']=max(r['global_power_relative_difference'] for r in out['rows'])
assert len(saved)==149 and sum(len(r['R_kpc']) for r in saved)==3150
(P/'occupancy-transport-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out['scores'],indent=2));print(out['maximum_grid_speed_difference_kms'],out['maximum_global_power_relative_difference'])
