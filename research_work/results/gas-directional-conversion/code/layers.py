"""Post-screen toy attribution: fixed gas column, clumping and layer order.

Uses exactly the original JR-5 transfer rates; no new fit, no parameters selected
from observational residuals. No claim that real gas has these layers.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from core import transfer


def layers(parts,alpha=.3,beta=.3):
    widths=np.array([p[0] for p in parts]);gas=np.array([p[1] for p in parts])[:,None]
    edges=np.r_[0.,np.cumsum(widths)];f=transfer(edges,gas,alpha,beta)[:,0]
    state=np.array([.5,.5]);history=[]
    for width,rho in parts:
        a=alpha*rho/1e6;b=beta*rho/1e6*rho/1e7
        generator=np.array([[-a,b],[a,-b]])
        state=expm(generator*width)@state
        history.append(state.tolist())
    return dict(layers=[dict(width_kpc=w,rho_Msun_kpc3=r) for w,r in parts],
        column_Msun_pc2=float(np.sum(widths*gas[:,0])/1e6),
        squared_density_path_integral=float(np.sum(widths*gas[:,0]**2)),
        final_companion_fraction=float(f[-1]),final_photon_fraction=float(1-f[-1]),
        stationary_bound_source_ratio_to_gasfree=float(2*f[-1]),
        independent_matrix_exponential_error=float(np.max(abs(np.array(history)[:,1]-f))),
        energy_sum_error=float(np.max(abs(np.array(history).sum(1)-1))))


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    trials=dict(diffuse_then_dense=[(2.5,2e6),(.1,5e7)],dense_then_diffuse=[(.1,5e7),(2.5,2e6)],
       homogeneous=[(1.,1e7)],clumped_same_length=[(.2,5e7),(.8,0.)],spread_same_column=[(5.,2e6)])
    result=dict(scope='post-screen attribution toy; no added data or refit',alpha=.3,beta=.3,
        cases={k:layers(v) for k,v in trials.items()})
    result['order_ratio']=result['cases']['dense_then_diffuse']['final_companion_fraction']/result['cases']['diffuse_then_dense']['final_companion_fraction']
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
