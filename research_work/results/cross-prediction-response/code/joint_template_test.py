#!/usr/bin/env python3
"""Score a single model's motion AND lensing predictions without unit mistakes.

Input NPZ (allow_pickle=False) contains vectors observed, baseline, candidate,
block and a full covariance matrix. block entries must be motion or lensing.
Covariance must include the noise/systematic model appropriate to the scored
observables and preserve known cross-block correlations. Optional nuisance
marginalization: nuisance_jacobian B, nuisance_covariance S, with zero-mean
Gaussian nuisance prior; C_eff=C+B*S*B.T. This is a local linear approximation.

The shape and amplitude of candidate must be specified without fitting the
scored data for a predictive interpretation. A=0 is baseline, A=1 is candidate.
The fitted amplitudes are DIAGNOSTICS, not permission to add independently
adjustable matter/light couplings to the physical model. No p-values or
physical-validity claims are produced. Full nonlinear posterior/image-model
analyses remain necessary. No galaxy simulation is performed by this script.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor, cho_solve


def evaluate(y, base, candidate, covariance, labels, B=None, S=None) -> dict:
    y,base,candidate=(np.asarray(x,dtype=float) for x in (y,base,candidate))
    labels=np.asarray(labels).astype(str)
    n=y.size
    if y.ndim!=1 or base.shape!=y.shape or candidate.shape!=y.shape or labels.shape!=y.shape:
        raise ValueError('Observed, baseline, candidate and block must be equal-length vectors')
    if not set(labels)=={'motion','lensing'}:raise ValueError('Both motion and lensing blocks required')
    C=np.asarray(covariance,dtype=float)
    if C.shape!=(n,n) or not np.isfinite(C).all():raise ValueError('Invalid covariance')
    if not all(np.isfinite(x).all() for x in (y,base,candidate)):raise ValueError('Nonfinite vectors')
    if not np.allclose(C,C.T,rtol=1e-10,atol=1e-12):raise ValueError('Covariance must be symmetric')
    if (B is None)!=(S is None):raise ValueError('Provide both nuisance matrices or neither')
    if B is not None:
        B,S=np.asarray(B,float),np.asarray(S,float)
        if B.shape[0]!=n or S.shape!=(B.shape[1],B.shape[1]):raise ValueError('Nuisance shape mismatch')
        if not np.isfinite(B).all() or not np.isfinite(S).all() or not np.allclose(S,S.T):
            raise ValueError('Invalid nuisance matrices')
        if np.min(np.linalg.eigvalsh(S)) < -1e-12:raise ValueError('Prior covariance must be positive semidefinite')
        C=C+B@S@B.T
    factor=cho_factor(C,lower=True,check_finite=True)
    inv=lambda a:cho_solve(factor,a)
    residual=y-base;template=candidate-base
    info=float(template@inv(template))
    chi0=float(residual@inv(residual))
    chi1=float((residual-template)@inv(residual-template))
    result=dict(status='Prediction diagnostic, not full theory validation',n_observations=n,
      n_motion=int(np.sum(labels=='motion')),n_lensing=int(np.sum(labels=='lensing')),
      baseline_chi2=chi0,fixed_candidate_chi2=chi1,
      fixed_candidate_chi2_improvement=chi0-chi1,
      nuisance_marginalized=B is not None,
      expected_squared_SNR_for_candidate_vs_baseline=info,
      conditional_on_covariance_and_local_forward_templates=True)
    if not np.isfinite(info) or info<=1e-20:
        result['amplitude_identifiable']=False
        return result
    A=float(template@inv(residual)/info)
    result.update(amplitude_identifiable=True,common_diagnostic_amplitude=A,
      common_amplitude_standard_error=float(1/np.sqrt(info)),
      fitted_common_amplitude_chi2=float((residual-A*template)@inv(residual-A*template)))
    X=np.column_stack([np.where(labels==block,template,0.) for block in ('motion','lensing')])
    F=X.T@inv(X)
    if np.linalg.matrix_rank(F,tol=max(np.max(abs(F)),1.)*1e-12)<2:
        result['separate_diagnostic_amplitudes_identifiable']=False
        return result
    acov=np.linalg.inv(F);avec=acov@X.T@inv(residual);diff=np.array([1.,-1.])
    ds=float(np.sqrt(diff@acov@diff))
    result.update(separate_diagnostic_amplitudes_identifiable=True,
      diagnostic_amplitude_motion=float(avec[0]),diagnostic_amplitude_lensing=float(avec[1]),
      diagnostic_amplitude_covariance=acov.tolist(),
      motion_minus_lensing_amplitude=float(diff@avec),
      amplitude_difference_standard_error=ds,
      amplitude_difference_in_standard_errors=float(diff@avec/ds))
    return result


def self_test() -> dict:
    # Synthetic algebra control, explicitly not galaxy data.
    t=np.array([1.,2.,-.5,1.2,2.1,-.6]);base=np.arange(6,dtype=float)
    C=np.diag([.5,.4,.3,.2,.7,.9]);C+=.03*np.ones((6,6))
    lab=np.array(['motion']*3+['lensing']*3)
    result=evaluate(base+1.7*t,base,base+t,C,lab)
    assert abs(result['common_diagnostic_amplitude']-1.7)<1e-12
    assert abs(result['diagnostic_amplitude_motion']-1.7)<1e-12
    assert abs(result['diagnostic_amplitude_lensing']-1.7)<1e-12
    return {'status':'synthetic algebra control only','passed':True,'results':result}


def main():
    ap=argparse.ArgumentParser(__doc__)
    ap.add_argument('--input',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--self-test',action='store_true')
    args=ap.parse_args()
    if args.output.exists():raise FileExistsError('Choose a fresh output filename')
    if args.self_test:result=self_test()
    else:
        if args.input is None:ap.error('--input required unless --self-test')
        with np.load(args.input,allow_pickle=False) as d:
            result=evaluate(d['observed'],d['baseline'],d['candidate'],d['covariance'],d['block'],
              d['nuisance_jacobian'] if 'nuisance_jacobian' in d else None,
              d['nuisance_covariance'] if 'nuisance_covariance' in d else None)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))

if __name__=='__main__':main()
