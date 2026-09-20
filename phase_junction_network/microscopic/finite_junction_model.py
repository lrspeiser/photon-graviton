#!/usr/bin/env python3
"""Finite-state checks for the Phase Junction microscopic prototype."""
from __future__ import annotations
import argparse, itertools, json
import numpy as np


def spin_ops(S: float):
    m=np.arange(-S,S+1,dtype=float); d=len(m)
    E=np.diag(m).astype(complex); R=np.zeros((d,d),complex)
    for j,x in enumerate(m[:-1]): R[j+1,j]=np.sqrt(S*(S+1)-x*(x+1))
    R/=np.max(np.abs(R))
    return E,R,m


def kron_all(xs):
    out=xs[0]
    for x in xs[1:]: out=np.kron(out,x)
    return out


def embed(op,link,d):
    I=np.eye(d, dtype=complex)
    return kron_all([op if i==link else I for i in range(4)])


def plaquette(S: float,U=1.0,K=0.4):
    E,R,m=spin_ops(S); d=len(m)
    es=[embed(E,i,d) for i in range(4)]
    rs=[embed(R,i,d) for i in range(4)]; ls=[x.conj().T for x in rs]
    G=[es[0]+es[3], es[1]-es[0], -es[1]-es[2], es[2]-es[3]]
    W=rs[0]@rs[1]@ls[2]@ls[3]
    HE=(U/2)*sum(x@x for x in es); H=HE-(K/2)*(W+W.conj().T)
    comm=max(float(np.max(np.abs(H@g-g@H))) for g in G)
    neutral=[]
    for i,s in enumerate(itertools.product(m,repeat=4)):
        e0,e1,e2,e3=s
        if max(abs(e0+e3),abs(e1-e0),abs(-e1-e2),abs(e2-e3))<1e-12: neutral.append(i)
    Hn=H[np.ix_(neutral,neutral)]; En=HE[np.ix_(neutral,neutral)]
    eener=np.unique(np.round(np.diag(En).real,14))
    return {
      "spin":S,"local_link_dimension":d,"full_hilbert_dimension":d**4,
      "neutral_sector_dimension":len(neutral),"gauge_commutator_max_abs":comm,
      "neutral_electric_energies":[float(x) for x in eener],
      "neutral_hamiltonian_eigenvalues":[float(x) for x in np.linalg.eigvalsh(Hn)],
      "electric_stiffness_nontrivial":bool(len(eener)>1)}


def inv(a,p): return pow(int(a)%p,-1,p)


def rank_mod(A,p):
    A=np.asarray(A,dtype=np.int64)%p; r=0
    for c in range(A.shape[1]):
        q=next((i for i in range(r,A.shape[0]) if A[i,c]%p),None)
        if q is None: continue
        A[[r,q]]=A[[q,r]]; A[r]=(A[r]*inv(A[r,c],p))%p
        for i in range(A.shape[0]):
            if i!=r and A[i,c]%p: A[i]=(A[i]-A[i,c]*A[r])%p
        r+=1
        if r==A.shape[0]: break
    return r


def idx(x,y,z,L): return (x*L+y)*L+z


def derivative(L,axis,p):
    N=L**3; D=np.zeros((N,N),dtype=np.int64); half=inv(2,p)
    for x,y,z in itertools.product(range(L),repeat=3):
        q=[x,y,z]; plus=q.copy(); minus=q.copy()
        plus[axis]=(plus[axis]+1)%L; minus[axis]=(minus[axis]-1)%L
        row=idx(x,y,z,L); D[row,idx(*plus,L)]=half; D[row,idx(*minus,L)]=-half
    return D%p


def gravity_code(L,p):
    N=L**3; dx,dy,dz=[derivative(L,a,p) for a in range(3)]
    d2=[x@x%p for x in (dx,dy,dz)]; half=inv(2,p)
    C=np.zeros((N,6*N),dtype=np.int64)
    C[:,0*N:1*N]=-(d2[1]+d2[2]); C[:,1*N:2*N]=-(d2[0]+d2[2]); C[:,2*N:3*N]=-(d2[0]+d2[1])
    C[:,3*N:4*N]=2*(dx@dy); C[:,4*N:5*N]=2*(dx@dz); C[:,5*N:6*N]=2*(dy@dz); C%=p
    G=np.zeros((3*N,6*N),dtype=np.int64)
    G[0*N:1*N,0*N:1*N]=dx; G[0*N:1*N,3*N:4*N]=half*dy; G[0*N:1*N,4*N:5*N]=half*dz
    G[1*N:2*N,1*N:2*N]=dy; G[1*N:2*N,3*N:4*N]=half*dx; G[1*N:2*N,5*N:6*N]=half*dz
    G[2*N:3*N,2*N:3*N]=dz; G[2*N:3*N,4*N:5*N]=half*dx; G[2*N:3*N,5*N:6*N]=half*dy; G%=p
    rc,rg=rank_mod(C,p),rank_mod(G,p); logical=6*N-rc-rg
    return {"lattice_size":L,"prime":p,"sites":N,"component_qudits":6*N,
      "scalar_constraint_rank":rc,"vector_constraint_rank":rg,
      "exact_commutator_nonzero_entries":int(np.count_nonzero((C@G.T)%p)),
      "logical_qudits":logical,"expected_logical_qudits":2*N+4,
      "local_modes_per_site_after_global_modes":float((logical-4)/N)}


def oscillator(d):
    a=np.zeros((d,d),complex)
    for n in range(1,d): a[n-1,n]=np.sqrt(n)
    q=(a+a.conj().T)/np.sqrt(2); mom=(a-a.conj().T)/(1j*np.sqrt(2))
    C=(q@mom-mom@q)/1j; I=np.eye(d); P=I[:,:d-1]
    return {"dimension":d,"commutator_over_i_eigenvalues":[float(x) for x in np.linalg.eigvalsh(C).real],
      "full_space_operator_error":float(np.linalg.norm(C-I,2)),
      "low_subspace_operator_error":float(np.linalg.norm(P.conj().T@(C-I)@P,2))}


def collective(S):
    return {"spin":S,"local_dimension":int(2*S+1),
      "one_excitation_relative_error":float(1/S),
      "two_excitation_relative_error":float(min(1,2/S))}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--output"); args=ap.parse_args()
    em=[plaquette(S) for S in (0.5,1.0,1.5,2.0)]
    grav=[gravity_code(3,5),gravity_code(5,7)]
    oscillators=[oscillator(d) for d in (3,4,5)]
    report={
      "status":"finite-state kinematic prototype; dynamics and impedance closure remain open",
      "electromagnetic_quantum_link":em,
      "gravity_discrete_constraint_code":grav,
      "canonical_commutator":{
        "finite_dimension_no_go":"No finite matrices Q,P satisfy [Q,P]=iI on the full local Hilbert space because tr([Q,P])=0 while tr(iI)=i*d.",
        "truncated_oscillators":oscillators,
        "collective_spin_low_excitation_errors":[collective(S) for S in (1.,2.,4.,8.,16.)]},
      "decisions":{
        "minimal_useful_electromagnetic_link_spin":1.0,
        "minimal_useful_electromagnetic_link_dimension":3,
        "gravity_exact_finite_option":"Exact GF(p) Weyl constraints leave two local logical qudits per site plus four periodic global modes.",
        "gravity_continuous_option":"Continuous canonical frame symmetry can only be emergent below a truncation or collective-spin boundary.",
        "not_yet_solved":["linear helicity-2 finite dynamics","nonlinear closure","Z_g/Z_A","chiral matter and mass hierarchy","vacuum-volume term"]}}
    assert not em[0]["electric_stiffness_nontrivial"] and em[1]["electric_stiffness_nontrivial"]
    assert max(x["gauge_commutator_max_abs"] for x in em)<1e-12
    assert all(x["exact_commutator_nonzero_entries"]==0 and x["logical_qudits"]==x["expected_logical_qudits"] for x in grav)
    assert max(x["low_subspace_operator_error"] for x in oscillators)<1e-12
    text=json.dumps(report,indent=2,sort_keys=True); print(text)
    if args.output:
        with open(args.output,"w",encoding="utf-8") as f: f.write(text+"\n")

if __name__=="__main__": main()
