#!/usr/bin/env python3
"""Verify a local first-order frame/connection representation of Fierz-Pauli stiffness."""
from __future__ import annotations
import argparse, json
import numpy as np


def sym_basis():
    out=[]
    for i in range(3):
        x=np.zeros((3,3)); x[i,i]=1.; out.append(x)
    for i,j in ((0,1),(0,2),(1,2)):
        x=np.zeros((3,3)); x[i,j]=x[j,i]=1/np.sqrt(2); out.append(x)
    return np.asarray(out)


def conn_basis():
    out=[]
    for a in range(3):
        for b in range(3):
            for c in range(b,3):
                x=np.zeros((3,3,3)); x[a,b,c]=1.
                if b!=c: x[a,c,b]=1.
                out.append(x)
    return np.asarray(out)

SB=sym_basis(); CB=conn_basis()


def hmat(v): return np.tensordot(v,SB,axes=1)

def hvec(x): return np.asarray([np.sum(b*x) for b in SB])

def cvec(x): return np.asarray([np.sum(b*x)/np.sum(b*b) for b in CB])


def q_connection(C):
    value=0.
    for i in range(3):
        for a in range(3):
            for b in range(3):
                value += C[a,b,i]*C[b,a,i] - C[a,i,i]*C[b,a,b]
    return float(value)


def quadratic_matrix(basis,value):
    n=len(basis); M=np.zeros((n,n))
    for i in range(n):
        M[i,i]=value(basis[i])
        for j in range(i):
            M[i,j]=M[j,i]=0.5*(value(basis[i]+basis[j])-value(basis[i])-value(basis[j]))
    return M

J=quadratic_matrix(CB,q_connection)


def gamma(h,k):
    G=np.zeros((3,3,3))
    for a in range(3):
        for b in range(3):
            for c in range(3):
                G[a,b,c]=0.5*(k[b]*h[a,c]+k[c]*h[a,b]-k[a]*h[b,c])
    return G


def gamma_map(k):
    return np.column_stack([cvec(gamma(b,k)) for b in SB])


def fp_value(h,k):
    k2=float(k@k); tr=float(np.trace(h)); kh=k@h
    return float(k2*np.sum(h*h)-2*kh@kh+2*(kh@k)*tr-k2*tr*tr)


def fp_matrix(k): return quadratic_matrix(SB,lambda h:fp_value(h,k))


def gauge_matrix(k):
    return np.column_stack([hvec(np.outer(k,np.eye(3)[i])+np.outer(np.eye(3)[i],k)) for i in range(3)])


def null_space(A,tol=1e-11):
    _,s,vh=np.linalg.svd(A,full_matrices=True); rank=int(np.sum(s>tol))
    return vh[rank:].T


def tt_basis(k):
    constraints=np.zeros((4,6))
    for j,b in enumerate(SB): constraints[:3,j]=b@k; constraints[3,j]=np.trace(b)
    return null_space(constraints)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--samples",type=int,default=500); ap.add_argument("--seed",type=int,default=43); ap.add_argument("--output"); a=ap.parse_args()
    rng=np.random.default_rng(a.seed)
    max_schur=max_gauge=max_tt=max_identity=0.; failures=0
    representative=None
    for _ in range(a.samples):
        k=rng.normal(size=3); k2=float(k@k)
        A=gamma_map(k); M=fp_matrix(k); schur=-4*A.T@J@A
        max_schur=max(max_schur,float(np.linalg.norm(schur-M,ord=2))/max(k2,1e-30))
        max_gauge=max(max_gauge,float(np.linalg.norm(M@gauge_matrix(k),ord=2))/max(k2**1.5,1e-30))
        T=tt_basis(k); eig=np.linalg.eigvalsh(T.T@M@T)
        max_tt=max(max_tt,float(np.max(np.abs(eig-k2)))/k2)
        h=rng.normal(size=6); H=hmat(h); G=gamma(H,k)
        lhs=fp_value(H,k); rhs=-4*q_connection(G)
        max_identity=max(max_identity,abs(lhs-rhs)/max(abs(lhs),k2*np.dot(h,h),1e-30))
        if T.shape!=(6,2) or np.any(eig<=0) or np.linalg.matrix_rank(J)!=18: failures+=1
        if representative is None: representative={"k_squared":k2,"tt_eigenvalues":[float(x) for x in eig]}
    report={
      "status":"local first-order classical architecture verified; finite quantum implementation and nonlinear closure remain open",
      "samples":a.samples,"failures":failures,"connection_components":18,"connection_quadratic_rank":int(np.linalg.matrix_rank(J)),
      "connection_quadratic_eigenvalue_range":[float(np.min(np.linalg.eigvalsh(J))),float(np.max(np.linalg.eigvalsh(J)))],
      "max_relative_schur_complement_error":max_schur,"max_relative_gamma_identity_error":max_identity,
      "max_relative_gauge_null_error":max_gauge,"max_relative_tt_eigenvalue_error":max_tt,
      "representative":representative,
      "identity":"V_FP(h,k) = -4 Q(Gamma[h]); V_aux(h,C)=4(Q(C)-2B(C,Gamma[h])) has stationary C=Gamma[h] and Schur complement V_FP.",
      "decision":"Advance the frame-plus-independent-connection architecture. It is local with one derivative in the mixed term and escapes the k^3 invariant-square obstruction, but it has not yet supplied an exact finite continuous gauge symmetry or a nonlinear quantum phase."}
    assert failures==0 and max_schur<1e-12 and max_identity<1e-12 and max_gauge<1e-12 and max_tt<1e-12
    text=json.dumps(report,indent=2,sort_keys=True); print(text)
    if a.output:
        with open(a.output,"w",encoding="utf-8") as f: f.write(text+"\n")

if __name__=="__main__": main()
