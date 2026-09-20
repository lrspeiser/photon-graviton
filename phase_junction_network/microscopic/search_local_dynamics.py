#!/usr/bin/env python3
"""Find the minimum derivative order of exact local gravity invariants."""
from __future__ import annotations
import argparse, json
import numpy as np

COMP=("xx","yy","zz","xy","xz","yz")


def rank_mod(A,p):
    A=np.asarray(A,dtype=np.int64)%p; r=0
    for c in range(A.shape[1]):
        q=next((i for i in range(r,A.shape[0]) if A[i,c]%p),None)
        if q is None: continue
        A[[r,q]]=A[[q,r]]; A[r]=(A[r]*pow(int(A[r,c]),-1,p))%p
        for i in range(A.shape[0]):
            if i!=r and A[i,c]%p: A[i]=(A[i]-A[i,c]*A[r])%p
        r+=1
        if r==A.shape[0]: break
    return r


def monomials(n):
    return [(x,y,n-x-y) for x in range(n+1) for y in range(n-x+1)]


def add(a,b): return tuple(a[i]+b[i] for i in range(3))


def coordinate_matrix(order):
    src=monomials(order); dst=monomials(order+1); where={x:i for i,x in enumerate(dst)}
    M=np.zeros((3*len(dst),6*len(src)),dtype=np.int64)
    trans={
      "xx":((0,0,2),),"yy":((1,1,2),),"zz":((2,2,2),),
      "xy":((0,1,1),(1,0,1)),"xz":((0,2,1),(2,0,1)),"yz":((1,2,1),(2,1,1))}
    unit=((1,0,0),(0,1,0),(0,0,1))
    for ci,c in enumerate(COMP):
        for mi,e in enumerate(src):
            col=ci*len(src)+mi
            for gauge,axis,coef in trans[c]:
                M[gauge*len(dst)+where[add(e,unit[axis])],col]+=coef
    return M


def momentum_matrix(order):
    src=monomials(order); dst=monomials(order+2); where={x:i for i,x in enumerate(dst)}
    M=np.zeros((len(dst),6*len(src)),dtype=np.int64)
    direction={
      "xx":(((0,2,0),-1),((0,0,2),-1)),
      "yy":(((2,0,0),-1),((0,0,2),-1)),
      "zz":(((2,0,0),-1),((0,2,0),-1)),
      "xy":(((1,1,0),2),),"xz":(((1,0,1),2),),"yz":(((0,1,1),2),)}
    for ci,c in enumerate(COMP):
        for mi,e in enumerate(src):
            col=ci*len(src)+mi
            for g,coef in direction[c]: M[where[add(e,g)],col]+=coef
    return M


def analyze(sector,order,M,primes):
    ranks={str(p):rank_mod(M,p) for p in primes}
    if len(set(ranks.values()))!=1: raise RuntimeError(f"prime-dependent rank: {sector} order {order}")
    r=next(iter(ranks.values()))
    return {"sector":sector,"derivative_order":order,"unknown_coefficients":M.shape[1],
      "equation_count":M.shape[0],"ranks_by_prime":ranks,"nullity":M.shape[1]-r}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--max-order",type=int,default=4); ap.add_argument("--output"); a=ap.parse_args()
    primes=(1_000_003,1_000_033)
    cr=[analyze("coordinate_under_vector_gauge",n,coordinate_matrix(n),primes) for n in range(a.max_order+1)]
    mr=[analyze("momentum_under_scalar_gauge",n,momentum_matrix(n),primes) for n in range(a.max_order+1)]
    cmin=next(x["derivative_order"] for x in cr if x["nullity"]); mmin=next(x["derivative_order"] for x in mr if x["nullity"])
    power=cmin+mmin
    report={
      "coordinate_invariants":cr,"momentum_invariants":mr,
      "minimum_coordinate_derivative_order":cmin,"minimum_momentum_derivative_order":mmin,
      "manifest_local_square_hamiltonian":{
        "coordinate_stiffness_power":2*cmin,"momentum_stiffness_power":2*mmin,
        "predicted_omega_power":power,"prediction":f"omega proportional to k^{power}"},
      "interpretation":"Squares of the lowest exact local invariants give a k^3 branch. Linear dispersion needs a globally invariant two-derivative bilinear, an area/frame factor linear in curvature, auxiliary variables, or emergent rather than exact microscopic continuous symmetry."}
    assert (cmin,mmin,power)==(2,1,3)
    text=json.dumps(report,indent=2,sort_keys=True); print(text)
    if a.output:
        with open(a.output,"w",encoding="utf-8") as f: f.write(text+"\n")

if __name__=="__main__": main()
