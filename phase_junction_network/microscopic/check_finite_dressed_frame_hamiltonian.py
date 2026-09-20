#!/usr/bin/env python3
"""Finite dressed-frame regulator: exact GF(p) quotient, clock spectra, and lock test."""
from __future__ import annotations
import argparse, json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from check_constrained_local_reduction import analyze_mode, lattice_momentum

COMP=("xx","yy","zz","xy","xz","yz")
DEFAULT={5:(4,5,6,7,8),7:(4,5,6,7),11:(3,4,5,6)}

def inv(a,p): return pow(int(a)%p,-1,p)

def rref(A,p):
    A=np.asarray(A,dtype=np.int64).copy()%p; piv=[]; r=0
    for c in range(A.shape[1]):
        q=next((i for i in range(r,A.shape[0]) if A[i,c]%p),None)
        if q is None: continue
        A[[r,q]]=A[[q,r]]; A[r]=A[r]*inv(A[r,c],p)%p
        for i in range(A.shape[0]):
            if i!=r and A[i,c]%p: A[i]=(A[i]-A[i,c]*A[r])%p
        piv.append(c); r+=1
        if r==A.shape[0]: break
    return A,piv

def rank(A,p): return len(rref(A,p)[1])

def null(A,p):
    A=np.asarray(A,dtype=np.int64)%p; R,piv=rref(A,p); free=[i for i in range(A.shape[1]) if i not in piv]; out=[]
    for f in free:
        v=np.zeros(A.shape[1],dtype=np.int64); v[f]=1
        for i,c in enumerate(piv): v[c]=-R[i,f]%p
        out.append(v)
    return np.asarray(out,dtype=np.int64) if out else np.zeros((0,A.shape[1]),dtype=np.int64)

def minv(A,p):
    n=len(A); R,piv=rref(np.c_[np.asarray(A,dtype=np.int64)%p,np.eye(n,dtype=np.int64)],p)
    if piv[:n]!=list(range(n)): raise ValueError("singular modular matrix")
    return R[:,n:]%p

def quotient(kernel,stab,p):
    base=np.asarray(stab,dtype=np.int64)%p; r0=rank(base,p); out=[]
    for v in sorted(kernel,key=lambda x:(np.count_nonzero(x),tuple(int(y) for y in x))):
        if rank(np.vstack([base,*out,v]),p)>r0+len(out): out.append(v%p)
    return np.asarray(out,dtype=np.int64)

def center(A,p): return (((np.asarray(A,dtype=np.int64)%p)+p//2)%p)-p//2

def logical_block(p):
    C=np.array([[1,1,0,0,0,0]],dtype=np.int64)%p
    G=np.zeros((3,6),dtype=np.int64); G[0,4]=G[1,5]=inv(2,p); G[2,2]=1
    Z=quotient(null(G,p),C,p); X0=quotient(null(C,p),G,p); X=minv(Z@X0.T%p,p).T@X0%p
    return {"prime":p,"component_order":list(COMP),"scalar_Z_stabilizer":C.tolist(),"vector_X_stabilizers":G.tolist(),
      "logical_qudits":len(Z),"expected_logical_qudits":6-rank(C,p)-rank(G,p),
      "logical_Z_exponents_centered":center(Z,p).tolist(),"logical_X_exponents_centered":center(X,p).tolist(),
      "logical_pairing_mod_p":(Z@X.T%p).tolist(),"maximum_Z_vector_commutator_mod_p":int(np.max((G@Z.T)%p)),
      "maximum_X_scalar_commutator_mod_p":int(np.max((C@X.T)%p)),"uses_transverse_traceless_projector":False}

def derivative(L,p):
    D=np.zeros((L,L),dtype=np.int64); h=inv(2,p)
    for x in range(L): D[x,(x+1)%L]=h; D[x,(x-1)%L]=-h
    return D%p

def chain_count(L,p):
    D=derivative(L,p); D2=D@D%p; C=np.zeros((L,6*L),dtype=np.int64); C[:,:L]=-D2; C[:,L:2*L]=-D2
    G=np.zeros((3*L,6*L),dtype=np.int64); h=inv(2,p); G[:L,4*L:5*L]=h*D; G[L:2*L,5*L:6*L]=h*D; G[2*L:3*L,2*L:3*L]=D; C%=p; G%=p
    k=6*L-rank(C,p)-rank(G,p)
    return {"lattice_size":L,"prime":p,"scalar_constraint_rank":rank(C,p),"vector_constraint_rank":rank(G,p),
      "scalar_vector_commutator_nonzero_entries":int(np.count_nonzero(C@G.T%p)),
      "logical_qudits_including_global_sector":k,"expected_logical_qudits_including_global_sector":2*L+4,
      "propagating_logical_qudits_after_k0_separation":2*L,"global_sector_qudits_separated":4}

def digits(states,p,L):
    powers=p**np.arange(L-1,dtype=np.int64); r=((states[:,None]//powers)%p).astype(np.int16); q=np.zeros((len(states),L),dtype=np.int16); q[:,1:]=r
    return powers,r,q

def encode(r,p,powers): return np.sum((r%p)*powers,axis=1,dtype=np.int64)

def clock_H(p,L,scale):
    dim=p**(L-1); s=np.arange(dim,dtype=np.int64); powers,r,q=digits(s,p,L); d=(np.roll(q,-1,axis=1)-q)%p
    rows=[s]; cols=[s]; data=[scale*(np.sum(1-np.cos(2*np.pi*d/p),axis=1)+L)]; a=-.5*scale
    for j in range(L-1):
        rp=r.copy(); rp[:,j]=(rp[:,j]+1)%p; rm=r.copy(); rm[:,j]=(rm[:,j]-1)%p
        rows += [s,s]; cols += [encode(rp,p,powers),encode(rm,p,powers)]; data += [np.full(dim,a),np.full(dim,a)]
    rows += [s,s]; cols += [encode(r-1,p,powers),encode(r+1,p,powers)]; data += [np.full(dim,a),np.full(dim,a)]
    return coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),shape=(dim,dim)).tocsr(),r,powers

def translate_perm(r,p,powers):
    q=np.zeros((len(r),r.shape[1]+1),dtype=np.int16); q[:,1:]=r; qp=np.roll(q,1,axis=1); qp=(qp-qp[:,[0]])%p
    return encode(qp[:,1:],p,powers)

def resolve(E,V,perm,L,tol=1e-7):
    ip=np.empty_like(perm); ip[perm]=np.arange(len(perm)); out=[]; i=0
    while i<len(E):
        j=i+1
        while j<len(E) and abs(E[j]-E[i])<tol: j+=1
        W=V[:,i:j]; phase,U=np.linalg.eig(W.conj().T@W[ip,:]); WV=W@U
        for n,z in enumerate(phase):
            k=int(np.rint((np.angle(z)%(2*np.pi))*L/(2*np.pi)))%L; v=WV[:,n]; v/=np.linalg.norm(v); out.append((float(E[i]),k,float(abs(z)),v))
        i=j
    return sorted(out,key=lambda x:x[0])

def boundary(v,r,p):
    L=r.shape[1]+1; q=np.zeros((len(r),L),dtype=np.int16); q[:,1:]=r; d=center(np.roll(q,-1,axis=1)-q,p); cut=np.abs(d)==p//2; pr=np.abs(v)**2
    return {"mean_fraction_of_bonds_on_compact_cut":float(pr@np.mean(cut,axis=1)),"probability_any_bond_on_compact_cut":float(pr@np.any(cut,axis=1))}

def lock_gap(p):
    w=np.exp(2j*np.pi/p); Z=np.diag([w**x for x in range(p)]); X=np.zeros((p,p),complex)
    for x in range(p): X[(x+1)%p,x]=1
    H=p/(4*np.pi)*(4*np.eye(p)-X-X.conj().T-Z-Z.conj().T); e=np.linalg.eigvalsh(H); return float(e[1]-e[0])

def clock_row(p,L,scale,k):
    H,r,powers=clock_H(p,L,scale); E,V=eigsh(H,k=min(k,H.shape[0]-2),which="SA",tol=1e-11,maxiter=30000); ix=np.argsort(E); E=E[ix]; V=V[:,ix]
    states=resolve(E,V,translate_perm(r,p,powers),L); g=min(states,key=lambda x:x[0]); plus=next(x for x in states if x[1]==1 and x[0]-g[0]>1e-9); minus=next(x for x in states if x[1]==L-1 and abs(x[0]-plus[0])<1e-7)
    gap=plus[0]-g[0]; kh=2*np.sin(np.pi/L); lg=lock_gap(p)
    return {"prime":p,"lattice_size":L,"one_polarization_hilbert_dimension_after_k0_fix":H.shape[0],"ground_energy":g[0],
      "positive_momentum_index":plus[1],"negative_momentum_index":minus[1],"translation_eigenvalue_modulus_positive":plus[2],
      "single_branch_gap":gap,"two_polarization_first_gap_degeneracy":2,"lattice_momentum":kh,"gap_over_lattice_momentum":gap/kh,
      "connection_lock_gap":lg,"lock_gap_over_tensor_gap":lg/gap,"ground_boundary":boundary(g[3],r,p),"tensor_boundary":boundary(plus[3],r,p)}

def fit(rows):
    k=np.array([x["lattice_momentum"] for x in rows]); g=np.array([x["single_branch_gap"] for x in rows]); power,lp=np.polyfit(np.log(k),np.log(g),1); A=np.c_[k,k**3]; c,d=np.linalg.lstsq(A,g,rcond=None)[0]; pred=A@np.array([c,d]); ratio=g/k
    return {"gap_proportional_to_lattice_momentum_power":float(power),"log_fit_prefactor":float(np.exp(lp)),"linear_speed_c_g":float(c),"cubic_cutoff_coefficient":float(d),
      "maximum_relative_linear_plus_cubic_fit_residual":float(np.max(abs(pred-g)/g)),"gap_over_momentum_coefficient_of_variation":float(np.std(ratio)/np.mean(ratio)),"target_power":1.0}

def op_embed(A,i,n,p):
    out=np.array([[1.]],complex); I=np.eye(p)
    for j in range(n): out=np.kron(out,A if i==j else I)
    return out

def lock_demo(p=3):
    A=np.array([[1,1],[1,2]])%p; nc,nf=A.shape; n=nc+nf; X,Z=clock_ops(p); xs=[op_embed(X,i,n,p) for i in range(n)]; zs=[op_embed(Z,i,n,p) for i in range(n)]; I=np.eye(p**n); H=np.zeros_like(I,dtype=complex); s=p/(4*np.pi)
    for a in range(nc):
        R=zs[nf+a].copy()
        for i in range(nf): R=R@np.linalg.matrix_power(zs[i],int(-A[a,i])%p)
        H += s*(4*I-xs[nf+a]-xs[nf+a].conj().T-R-R.conj().T)
    xb=[]
    for i in range(nf):
        u=xs[i].copy()
        for a in range(nc): u=u@np.linalg.matrix_power(xs[nf+a],int(A[a,i])%p)
        xb.append(u)
    F=sum((2*I-u-u.conj().T for u in xb),np.zeros_like(H)); b=zs[1]@zs[0].conj().T; F=.125*(F+2*I-b-b.conj().T)
    e,v=np.linalg.eigh(H); mask=np.isclose(e,e[0],rtol=1e-10,atol=1e-10); P=v[:,mask]@v[:,mask].conj().T
    return {"prime":p,"hilbert_dimension":p**n,"ground_band_dimension":int(mask.sum()),"expected_ground_band_dimension":p**nf,
      "lock_gap":float(e[mask.sum()]-e[0]),"maximum_lock_frame_commutator":float(np.max(abs(H@F-F@H))),
      "ground_band_leakage_operator_norm":float(np.linalg.norm((I-P)@F@P,2))}

def clock_ops(p):
    w=np.exp(2j*np.pi/p); Z=np.diag([w**x for x in range(p)]); X=np.zeros((p,p),complex)
    for x in range(p): X[(x+1)%p,x]=1
    return X,Z

def symbol_check():
    sizes=(12,16,24,32,48,64); axial=[]; split=err=0.
    for L in sizes:
        k=lattice_momentum((1,0,0),L); r=analyze_mode(k); target=np.linalg.norm(k); split=max(split,abs(r.frequencies[1]-r.frequencies[0])/target); err=max(err,max(abs(x-target) for x in r.frequencies)/target); axial.append({"lattice_size":L,"frequencies":r.frequencies,"target_lattice_momentum":float(target)})
    an=[]
    for L in sizes:
        a=analyze_mode(lattice_momentum((3,0,0),L)).frequencies[0]; b=analyze_mode(lattice_momentum((2,2,1),L)).frequencies[0]; an.append({"lattice_size":L,"relative_anisotropy":abs(a-b)/(.5*(a+b))})
    slope,inter=np.polyfit(np.log([x["lattice_size"] for x in an]),np.log([x["relative_anisotropy"] for x in an]),1)
    return {"method":"full scalar/vector constrained reduction of the local Fierz-Pauli symbol; no TT projector","tested_lattice_sizes":list(sizes),"maximum_relative_polarization_split":split,"maximum_relative_frequency_error":err,"anisotropy_proportional_to_L_power":float(slope),"anisotropy_prefactor":float(np.exp(inter)),"uses_transverse_traceless_projector":False}

def parse_sets(s): return {int(x.split(':')[0]):tuple(map(int,x.split(':')[1].split(','))) for x in s.split(';')}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--size-sets",default="5:4,5,6,7,8;7:4,5,6,7;11:3,4,5,6"); ap.add_argument("--frame-scale",type=float,default=.25); ap.add_argument("--eigenpairs",type=int,default=10); ap.add_argument("--output"); a=ap.parse_args(); sets=parse_sets(a.size_sets)
    blocks=[logical_block(p) for p in sorted(sets)]; counts=[chain_count(3,5),chain_count(5,7)]; rows=[clock_row(p,L,a.frame_scale,a.eigenpairs) for p,Ls in sets.items() for L in Ls]; fits={str(p):fit([x for x in rows if x["prime"]==p]) for p in sorted(sets)}; reps=[max((x for x in rows if x["prime"]==p),key=lambda x:x["lattice_size"]) for p in sorted(sets)]; lock=lock_demo(); symbol=symbol_check()
    public_rows=[{k:x[k] for k in ("prime","lattice_size","one_polarization_hilbert_dimension_after_k0_fix","positive_momentum_index","negative_momentum_index","single_branch_gap","two_polarization_first_gap_degeneracy","lattice_momentum","gap_over_lattice_momentum","connection_lock_gap","lock_gap_over_tensor_gap")} for x in rows]
    public_reps=[{"prime":x["prime"],"lattice_size":x["lattice_size"],"ground_boundary":x["ground_boundary"],"tensor_boundary":x["tensor_boundary"]} for x in reps]
    report={"status":"finite dressed linear frame Hamiltonian/positive transfer regulator verified on the exact first-class quotient","scope":"Issue #2 linear finite regulator; coefficient derivation, nonlinear closure, and universality remain open",
      "finite_parent_construction":{"frame_coordinate":"Zbar_i = Z_hi","frame_shift":"Xbar_i = X_hi * product_a X_Ca^(A_ai)","first_class_treatment":"exact GF(p) scalar/vector stabilizer quotient or equivalent local group-averaged transfer step","one_polarization_hamiltonian":"lambda_g/2 sum_z [(2-X-Xdagger)+(2-Znext Zdagger-Z Znextdagger)]","polarizations":2,"frame_scale_lambda_g":a.frame_scale,"uses_transverse_traceless_projector":False},
      "finite_logical_reduction":blocks,"real_space_constraint_counts":counts,"finite_clock_spectrum_rows":public_rows,"dispersion_fits_by_prime":fits,"boundary_representatives":public_reps,"dressed_connection_lock":lock,"complete_three_dimensional_harmonic_symbol":symbol,"global_sector_policy":"propagating test excludes k=0; four periodic global qudits remain separate",
      "acceptance":{"declared_finite_positive_hamiltonian_or_transfer_matrix":True,"exact_finite_first_class_quotient":True,"two_low_energy_tensor_branches":True,"propagating_connection_scalar_or_vector_state":False,"linear_gap_with_cubic_cutoff_correction":True,"connection_lock_leakage_below_tolerance":True,"uses_nonlocal_TT_projector":False,"finite_dimension_boundary_error_measured":True},
      "limitations_owned_elsewhere":["#3 shared stiffnesses","#5 nonlinear gravity and volume/vacuum","#6 finite Coulomb phase/common cone","#8 interacting continuum and regulator universality"]}
    assert all(x["logical_qudits"]==x["expected_logical_qudits"]==2 and x["logical_pairing_mod_p"]==[[1,0],[0,1]] and x["maximum_Z_vector_commutator_mod_p"]==x["maximum_X_scalar_commutator_mod_p"]==0 for x in blocks)
    assert all(x["scalar_vector_commutator_nonzero_entries"]==0 and x["logical_qudits_including_global_sector"]==x["expected_logical_qudits_including_global_sector"] for x in counts)
    assert all(.95<x["gap_proportional_to_lattice_momentum_power"]<1.12 and x["maximum_relative_linear_plus_cubic_fit_residual"]<5e-3 and x["gap_over_momentum_coefficient_of_variation"]<.03 for x in fits.values())
    assert min(x["lock_gap_over_tensor_gap"] for x in rows)>2 and all(x["positive_momentum_index"]==1 and x["negative_momentum_index"]==x["lattice_size"]-1 for x in rows)
    bv=[x["tensor_boundary"]["mean_fraction_of_bonds_on_compact_cut"] for x in reps]; assert all(b<a for a,b in zip(bv,bv[1:])) and bv[-1]<1e-3
    assert lock["maximum_lock_frame_commutator"]<1e-12 and lock["ground_band_leakage_operator_norm"]<1e-12
    assert symbol["maximum_relative_polarization_split"]<1e-12 and symbol["maximum_relative_frequency_error"]<1e-12 and abs(symbol["anisotropy_proportional_to_L_power"]+2)<.1
    text=json.dumps(report,indent=2,sort_keys=True); print(text)
    if a.output: open(a.output,"w",encoding="utf-8").write(text+"\n")
if __name__=="__main__": main()
