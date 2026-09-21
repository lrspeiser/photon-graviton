"""Area/reference move and geometry-complete finite spatial matter candidate.

New architecture is explicit: a coherent loop/reference channel and a coalesced
site-block move decomposition. Existing x, Clifford and slab parameters are
inherited, not re-fit. This module does not assert gravitational constraints.
"""
from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import block_diag

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'microscopic'))
import chiral_matter_model as slab
from all_sector_move_common import X_ROOT

GAMMA = np.asarray(slab.SPATIAL_GAMMAS)
BETA = slab.GAMMA5
SPIN = np.asarray([np.kron(sigma, slab.I2)/2 for sigma in (slab.SX, slab.SY, slab.SZ)])


def completed(move: np.ndarray, x: float = X_ROOT) -> np.ndarray:
    return -x*(move+move.conj().T)+x*x*(move.conj().T@move+move@move.conj().T)


def area_move(area: np.ndarray, holonomy: np.ndarray) -> np.ndarray:
    """The difference is a microscopic reference-path amplitude, not a fit."""
    return area @ (holonomy-np.eye(len(holonomy))) / (2j)


def area_energy(area: np.ndarray, holonomy: np.ndarray, x: float = X_ROOT) -> float:
    return float(np.trace(completed(area_move(area, holonomy), x)).real / len(area))


def geometry(inverse_frames: np.ndarray):
    frames = np.asarray(inverse_frames, dtype=float)
    if frames.ndim != 3 or frames.shape[1:] != (3,3):
        raise ValueError('inverse_frames must have shape (sites,3,3)')
    det = np.linalg.det(frames)
    if not np.all(np.isfinite(frames)) or np.any(det <= 0):
        raise ValueError('oriented invertible finite frames required')
    volume = 1/det
    gamma = np.einsum('nai,abc->nibc', frames, GAMMA)
    metric = np.einsum('nai,naj->nij', frames, frames)
    area = np.empty((len(frames),3,3,4,4), dtype=complex)
    for i in range(3):
        for j in range(3):
            area[:,i,j] = volume[:,None,None]*(gamma[:,i]@gamma[:,j]-gamma[:,j]@gamma[:,i])/(2j)
    return volume, gamma, metric, area


def lattice(L: int):
    if L < 3:
        raise ValueError('L>=3 avoids ambiguous coalescing of periodic directed edges')
    sites = list(np.ndindex(L,L,L))
    lookup = {s:n for n,s in enumerate(sites)}
    neighbors = np.asarray([[lookup[tuple((v+(a==i))%L for a,v in enumerate(s))]
                              for i in range(3)] for s in sites])
    return sites, neighbors


def matter_matrix(L: int, width: int, charge: int, volume: np.ndarray,
                  gamma: np.ndarray, metric: np.ndarray, spin_links: np.ndarray,
                  gauge_phases: np.ndarray | None = None,
                  antiperiodic: bool = True) -> dict:
    """Canonical H=nu^-1/2 [Dirac flux + beta covariant Wilson] nu^-1/2.

    D_i=T_i-I; Wilson = 1/2 sum D_i^dagger (nu g^ij) D_j.
    Dirac = 1/2 sum {nu Gamma^i,(T_i-T_i^dagger)/(2i)}.
    All fifth-direction, mass and Wilson operators are present. Independent
    spin links allow a torsionful background; compatibility is not imposed.
    """
    sites, neighbors = lattice(L)
    N = len(sites)
    if width < 1 or volume.shape != (N,) or np.any(volume <= 0):
        raise ValueError('invalid width or volume')
    block = width*4
    dim = N*block
    phase = np.zeros((N,3)) if gauge_phases is None else np.asarray(gauge_phases)
    nu = np.diag(np.repeat(volume,block))
    invsqrt = np.diag(np.repeat(volume**(-.5),block))
    beta = np.kron(np.eye(N*width), BETA)
    shifts=[]
    for i in range(3):
        T=np.zeros((dim,dim), dtype=complex)
        for a,s in enumerate(sites):
            b=neighbors[a,i]
            bc=-1 if antiperiodic and s[i]==L-1 else 1
            T[a*block:(a+1)*block,b*block:(b+1)*block] = bc*np.exp(1j*charge*phase[a,i])*np.kron(np.eye(width),spin_links[a,i])
        shifts.append(T)
    identity=np.eye(dim)
    gradients=[t-identity for t in shifts]
    kinetic=np.zeros((dim,dim),dtype=complex)
    for i,t in enumerate(shifts):
        flux=block_diag(*[np.kron(np.eye(width), volume[a]*gamma[a,i]) for a in range(N)])
        K=(t-t.conj().T)/(2j)
        kinetic += (flux@K+K@flux)/2
    wilson=np.zeros_like(kinetic)
    for i in range(3):
        for j in range(3):
            weight=np.diag(np.repeat(volume*metric[:,i,j],block))
            wilson += gradients[i].conj().T@weight@gradients[j]/2
    kinetic=invsqrt@kinetic@invsqrt
    wilson=beta@invsqrt@wilson@invsqrt
    wall=slab.slab_hamiltonian((0,0,0),slab.wilson_mass(charge),width)
    onsite=np.kron(np.eye(N),wall)
    target=onsite+kinetic+wilson
    error=float(np.max(abs(target-target.conj().T)))
    if error>1.e-10:
        raise ArithmeticError(f'non-Hermitian target: {error}')
    # Normal-ordered ONE-BODY component of the common-rule Fock lift:
    # one on-site move -O/2 and one move -K_ab per unordered spatial pair.
    # The exact quartic partner is defined/tested by full_fock_completion;
    # it must not be silently discarded in the interacting theory.
    diagonal_completion=np.zeros_like(target)
    for a in range(N):
        a_slice=slice(a*block,(a+1)*block)
        O=target[a_slice,a_slice]
        diagonal_completion[a_slice,a_slice] += .5*O@O
        for b in range(a+1,N):
            b_slice=slice(b*block,(b+1)*block)
            K=target[a_slice,b_slice]
            diagonal_completion[a_slice,a_slice] += K@K.conj().T
            diagonal_completion[b_slice,b_slice] += K.conj().T@K
    lifted=X_ROOT*target+X_ROOT**2*diagonal_completion
    return {'target':target,'completed':lifted,'kinetic':kinetic,'wilson':wilson,
            'mass_wall':onsite,'completion':diagonal_completion,
            'measure':nu,'shifts':shifts}


def plaquette_energy(L: int, area: np.ndarray, spin_links: np.ndarray) -> float:
    _,neighbors=lattice(L)
    total=0.
    for a in range(len(neighbors)):
        for i in range(3):
            for j in range(i+1,3):
                ai,aj=neighbors[a,i],neighbors[a,j]
                W=spin_links[a,i]@spin_links[ai,j]@spin_links[aj,i].conj().T@spin_links[a,j].conj().T
                total+=area_energy(area[a,i,j],W)
    return total


def lapse_source(H:np.ndarray, lapse:np.ndarray, block:int)->np.ndarray:
    """Declared symmetric energy-density assignment; includes every H term."""
    N=np.diag(np.repeat(lapse,block))
    return (N@H+H@N)/2


def transform_sources(gamma, area, links, neighbors, spin_rotations):
    S=spin_rotations
    gamma_new=S[:,None]@gamma@S[:,None].conj().swapaxes(-2,-1)
    area_new=S[:,None,None]@area@S[:,None,None].conj().swapaxes(-2,-1)
    link_new=np.empty_like(links)
    for i in range(3):
        link_new[:,i]=S@links[:,i]@S[neighbors[:,i]].conj().swapaxes(-2,-1)
    return gamma_new,area_new,link_new


def bilinear_fock(A: np.ndarray) -> np.ndarray:
    """Second-quantize c^dagger A c with exact CAR signs; small controls only."""
    n=len(A)
    if A.shape!=(n,n) or n>10:
        raise ValueError('Fock control supports square matrices up to ten modes')
    out=np.zeros((1<<n,1<<n),dtype=complex)
    entries=np.argwhere(abs(A)>1.e-15)
    for mask in range(1<<n):
        for i,j in entries:
            i,j=int(i),int(j)
            if not(mask>>j)&1:
                continue
            removed=mask^(1<<j)
            sign=(-1)**((mask&((1<<j)-1)).bit_count())
            if (removed>>i)&1:
                continue
            sign*=(-1)**((removed&((1<<i)-1)).bit_count())
            out[removed|(1<<i),mask]+=sign*A[i,j]
    return out


def full_fock_completion(target: np.ndarray, block: int, x: float=X_ROOT) -> dict:
    """Apply the rule AFTER second quantization, retaining quartic interactions.

    For bilinears Ahat,Bhat, Ahat Bhat = c†ABc +
    sum A_ij B_kl c_i† c_k† c_l c_j. Thus single-particle completion alone is
    only the normal-ordered quadratic component, not the full matter theory.
    """
    n=len(target)
    if n%block or n>10:
        raise ValueError('small complete-site Fock controls only')
    onebody=np.zeros_like(target,dtype=complex)
    full=np.zeros((1<<n,1<<n),dtype=complex)
    for a in range(n//block):
        sa=slice(a*block,(a+1)*block)
        move=np.zeros_like(target,dtype=complex)
        move[sa,sa]=-target[sa,sa]/2
        onebody+=completed(move,x); full+=completed(bilinear_fock(move),x)
        for b in range(a+1,n//block):
            sb=slice(b*block,(b+1)*block)
            move=np.zeros_like(target,dtype=complex)
            move[sa,sb]=-target[sa,sb]
            onebody+=completed(move,x); full+=completed(bilinear_fock(move),x)
    quadratic=bilinear_fock(onebody)
    return {'full':full,'onebody':onebody,'quadratic_fock':quadratic,'quartic_fock':full-quadratic}
