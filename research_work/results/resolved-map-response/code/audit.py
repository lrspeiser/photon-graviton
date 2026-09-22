"""Independent saved-output readback and selected full replay for JR-8.

All 288 outputs are checked arithmetically. Three primary states are re-solved;
this is not a replay of every historical experiment in the repository.
"""
import argparse,json,hashlib,sys
from pathlib import Path
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import splu
from response import Mapped,Setting
from fits_images import read_images,channel

def main():
    ap=argparse.ArgumentParser(__doc__);ap.add_argument('--root',type=Path,required=True);args=ap.parse_args();root=args.root;out=root/'results/audit.json'
    if out.exists():raise FileExistsError('Use a fresh archive copy for a new audit')
    maps=json.loads((root/'results/maps/maps.json').read_text());rows=json.loads((root/'results/response/primary.json').read_text());ref=json.loads((root/'results/response/refinement.json').read_text());sup=json.loads((root/'results/supplement.json').read_text());checks=[]
    def chk(name,value,bound=1e-8):
        value=float(value);checks.append(dict(name=name,value=value,bound=bound,passed=bool(np.isfinite(value) and abs(value)<=bound)))
    chk('primary_count',len(rows)-288,0);chk('refinement_count',len(ref)-12,0)
    for j,r in enumerate(rows+ref):
        d=r['diagnostics'];chk(f'energy_{j}',.5*d['Ep']+.05*d['Ec']-1);chk(f'companion_balance_{j}',d['forward']-d['reverse']-.05*d['Ec']);chk(f'linear_{j}',d['linear_residual']);chk(f'negative_{j}',min(d['min_state'],0))
        if r['setting']['family']=='constant':chk(f'constant_Ec_{j}',d['Ec']-5/3);chk(f'constant_Ep_{j}',d['Ep']-11/6)
        for c in d['sensitivity'].get('finite_difference_checks',[]):chk(f'adjoint_{j}_{c["kind"]}_{c["sign"]}',c['relative_error'],1e-5)
    for rec in maps:
        ident=rec['id'];f=next((root/'inputs/maps').glob(f'manga-{ident}-*.fits.gz'));chk(f'input_hash_{ident}',int(hashlib.sha256(f.read_bytes()).hexdigest()!=rec['input_sha256']),0)
        h=read_images(f);old=np.load(root/'inputs/maps'/f'{ident}-maps.npz');ha=channel(h['EMLINE_GFLUX'][0],'Ha-6564')
        for key,a in [('light',h['SPX_MFLUX'][1]),('star_velocity',h['STELLAR_VEL'][1]),('Halpha',h['EMLINE_GFLUX'][1][ha])]:chk(f'FITS_archive_{ident}_{key}',np.max(abs(a-old[key])),0)
        M=Mapped(rec);g=M.g;s=Setting();A,rhs,gas,*_=M.system(s);lu=splu(A);y=lu.solve(rhs);archive=np.load(root/'results/response'/(ident+'-fields.npz'));chk(f'replay_{ident}',np.max(abs(y-np.r_[archive['P'],archive['C']]))/np.max(abs(y)),1e-10)
        # Rotate both projected inputs by exactly pi/2, not just one source.
        shift=g.m.nphi//4;gg=np.roll(gas.reshape(g.shape),shift,axis=2).ravel();JJ=np.roll(rhs[:g.n].reshape(g.shape),shift,axis=2).ravel();AA,_,*_=M.system(s,gg);yy=splu(AA).solve(np.r_[JJ,np.zeros(g.n)])
        expected=np.r_[np.roll(y[:g.n].reshape(g.shape),shift,axis=2).ravel(),np.roll(y[g.n:].reshape(g.shape),shift,axis=2).ravel()]
        chk(f'whole_rotation_{ident}',np.max(abs(yy-expected))/np.max(abs(expected)),1e-10)
        # Implicit Euler removal: exact DISCRETE source+escape ledger. The
        # timestep is not being claimed as a continuum rate/trajectory result.
        weights=np.r_[g.w,g.w];loss=np.r_[.5*g.w,.05*g.w];initial=float(weights@y);escaped=0.;step=.01;prop=splu(eye(len(y),format='csc')+step*A);history=[]
        for _ in range(50):
            y=prop.solve(y);escaped+=step*float(loss@y);history.append(float(weights@y));chk(f'removal_ledger_{ident}_{_}',(weights@y+escaped-initial)/initial,1e-9)
        chk(f'removal_monotone_{ident}',max(np.diff(history).max(),0),1e-10)
    for r in sup['adjoint_observables']:
        for c in r['adjoint_checks']:chk(f'{r["object"]}_{c["observable"]}_adj',c['relative_error'],1e-5)
        for j,c in enumerate(r['counterfactuals']):
            if c['status']!='solved':chk('failed_counterfactual',1,0);continue
            chk(f'counterfactual_mass_{r["object"]}_{j}',c['gas_inventory_error']);chk(f'counterfactual_force_{r["object"]}_{j}',c['force_relative_change'],1e-8)
    for r in sup['boundary_controls']:chk(f'expanded_balance_{r["object"]}_{r["phase"]}',r['diagnostics']['energy_error'])
    cluster=json.loads((root/'results/cluster-response.json').read_text());cf=root/'inputs/clusters/cl2-inputs-xcop-profiles.json';chk('cluster_input_hash',int(hashlib.sha256(cf.read_bytes()).hexdigest()!=cluster['input_sha256']),0);chk('cluster_derivative',cluster['derivative_absolute_error'],1e-7);chk('cluster_count',len(cluster['clusters'])-12,0)
    result=dict(check_count=len(checks),failed=sum(not c['passed'] for c in checks),scope=__doc__,checks=checks)
    out.write_text(json.dumps(result,indent=2)+'\n');print('AUDIT',result['check_count'],result['failed'])
    if result['failed']:sys.exit(1)
if __name__=='__main__':main()
