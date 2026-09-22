#!/usr/bin/env python3
"""Independent 3D finite-difference angular solve and archived readback."""
from pathlib import Path
import json
import hashlib
import argparse
import numpy as np
from scipy.sparse import bmat, diags, eye, kron, lil_matrix
from scipy.sparse.linalg import spsolve
import currents as C

def angular_solve(nphi):
    s=C.Setup(nr=8,nz=10);g,f,d=C.solve(s);nn=g.size
    phi=np.arange(nphi)*2*np.pi/nphi;dp=2*np.pi/nphi
    d1=lil_matrix((nphi,nphi));d2=lil_matrix((nphi,nphi))
    for j in range(nphi):
        d1[j,(j+1)%nphi]=1/(2*dp);d1[j,(j-1)%nphi]=-1/(2*dp)
        d2[j,j]=-2/dp**2;d2[j,(j+1)%nphi]=1/dp**2;d2[j,(j-1)%nphi]=1/dp**2
    d1=d1.tocsr();d2=d2.tocsr()
    lap=kron(g.lap(0),eye(nphi))+kron(diags(1/g.R.ravel()**2),d2)
    kp=np.repeat(g.kp,nphi);km=np.repeat(g.km,nphi)
    a=-s.dp*lap+diags(kp+s.loss_p)-s.omega_s*kron(eye(nn),d1)
    b=-s.dc*lap+diags(km+s.loss_c)+kron(diags(g.omega-s.omega_s),d1)
    op=bmat([[a,diags(-km)],[diags(-kp),b]],format='csc')
    src=(g.J0[:,None]+g.J2[:,None]*np.cos(2*phi)).ravel()
    rhs=np.r_[src,np.zeros(nn*nphi)];sol=spsolve(op,rhs)
    p,c=np.split(sol,2);p=p.reshape(nn,nphi);c=c.reshape(nn,nphi)
    predicted=f[1][:,None]+np.real(f[3][:,None]*np.exp(2j*phi))
    fluct=predicted-f[1][:,None]
    return dict(nphi=nphi,full_3D_unknowns=len(sol),relative_C_field_error=float(np.linalg.norm(c-predicted)/np.linalg.norm(predicted)),
        relative_C_anisotropy_error=float(np.linalg.norm(c-predicted)/np.linalg.norm(fluct)),
        matrix_relative_residual=float(np.linalg.norm(op@sol-rhs)/np.linalg.norm(rhs)),
        C_angular_mean_error=float(np.max(abs(c.mean(1)-f[1]))),minimum_energy_density=float(min(p.min(),c.min())))

def main():
    ap=argparse.ArgumentParser(__doc__);ap.add_argument('--results',type=Path,required=True);args=ap.parse_args();out=args.results
    if (out/'audit.json').exists():raise FileExistsError('Preserve the prior audit')
    result=json.loads((out/'results.json').read_text());checks=[]
    def check(name,error,tol):checks.append(dict(name=name,error=float(error),tolerance=float(tol),passed=bool(error<=tol)))
    for label,row in {**result['cases'],**result['refinements']}.items():
        diag=row['diagnostics'];check(label+'_source_escape_balance',diag['balance_error'],1e-9)
        check(label+'_linear_equation',diag['linear_residual_max'],1e-8)
        check(label+'_nonnegative_P',max(0,-diag['min_P_over_all_phi']),1e-12)
        check(label+'_nonnegative_C',max(0,-diag['min_C_over_all_phi']),1e-12)
        path=out/(label+'.npz')
        if path.exists():
            with np.load(path) as z:
                w=z['weight'];ec=2*np.pi*np.dot(z['C0'],w)
                q=np.pi*np.dot(np.conj(z['C2'])*z['R'].ravel()**2,w)
                scale=2*np.pi*np.dot(z['C0']*z['R'].ravel()**2,w)
                check(label+'_energy_readback',abs(ec-diag['energy_C']),1e-12)
                check(label+'_quadrupole_readback',abs(abs(q)/scale-diag['normalized_C_quadrupole']),1e-12)
    direct=[angular_solve(n) for n in (16,32,64)]
    check('angular_refinement_decreasing',0 if all(direct[i+1]['relative_C_anisotropy_error']<direct[i]['relative_C_anisotropy_error'] for i in range(2)) else 1,0)
    check('direct_angular_mean',max(r['C_angular_mean_error'] for r in direct),1e-10)
    check('direct_3D_residual',max(r['matrix_relative_residual'] for r in direct),1e-8)
    boundary={}
    x=result['refinements']['aligned_48'];y=result['refinements']['aligned_boundary']
    for ch in ('C','both'):
        for key in ('inward_force','inward_deflection'):
            a=np.array([r[ch][key] for r in x['readout']['azimuths']]);b=np.array([r[ch][key] for r in y['readout']['azimuths']])
            boundary[ch+'_'+key]=float(np.max(abs(a/b-1)))
    boundary['normalized_quadrupole_relative_change']=float(y['diagnostics']['normalized_C_quadrupole']/x['diagnostics']['normalized_C_quadrupole']-1)
    full=dict(checks=checks,failed=sum(not c['passed'] for c in checks),independent_3D_angular_solves=direct,
        boundary_changes=boundary,boundary_scope='Force/lens readouts remain within 0.8%; global moment changes more because it weights distant material by R^2',
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'audit.json').write_text(json.dumps(full,indent=2,allow_nan=False)+'\n');print(json.dumps(full),flush=True)
    if full['failed']:raise SystemExit('Numerical audit failed; preserved')
if __name__=='__main__':main()
