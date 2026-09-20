"""CC-2S linear current response at toy and illustrative galactic speeds."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import subprocess
import numpy as np
from scipy.integrate import quad
from run_audit import HERE,ROOT,save,digest

G=.001;KAPPA=.5;R=.7;MASS=6.


def ring(n,speed,sense):
    theta=np.arange(n)*2*np.pi/n
    positions=R*np.column_stack((np.cos(theta),np.sin(theta),np.zeros(n)))
    velocities=sense*speed*np.column_stack((-np.sin(theta),np.cos(theta),np.zeros(n)))
    energy=MASS/n/np.sqrt(1-speed**2)
    return positions,velocities,energy


def field(point,positions,velocities,energy,eta):
    d=point-positions;distance=np.linalg.norm(d,axis=1)
    u=-G**2*energy*np.sum(1/distance)/(4*np.pi)
    grad=G**2*energy*np.sum(d/distance[:,None]**3,axis=0)/(4*np.pi)
    beta=-(KAPPA*eta)**2*energy*np.sum(velocities/distance[:,None],axis=0)/(4*np.pi)
    curl=(KAPPA*eta)**2*energy*np.sum(np.cross(d,velocities)/distance[:,None]**3,axis=0)/(4*np.pi)
    return u,grad,beta,curl


def evaluate(n,speed,sense,ratio,radius_ratio):
    positions,velocities,energy=ring(n,speed,sense)
    eta=G*ratio;r=R*radius_ratio
    u,grad,beta,curl=field(np.array([r,0,0]),positions,velocities,energy,eta)
    acceleration=-np.cross([0,speed,0],curl)
    star_ratio=-acceleration[0]/grad[0]
    offset=r-positions[:,1]
    scalar_bend=-G**2*energy*np.sum(1/offset)/(2*np.pi)
    vector_bend=-(KAPPA*eta)**2*energy*np.sum(velocities[:,0]/offset)/(2*np.pi)
    maxu=abs(u);maxbeta=np.linalg.norm(beta)
    for x in np.linspace(-20*r,20*r,81):
        local=field(np.array([x,r,0]),positions,velocities,energy,eta)
        maxu=max(maxu,abs(local[0]));maxbeta=max(maxbeta,np.linalg.norm(local[2]))
    return dict(n=n,speed=speed,sense=sense,coupling_ratio=ratio,radius_ratio=radius_ratio,
                scalar_star_acceleration=grad[0],directional_star_acceleration=-acceleration[0],star_fraction=star_ratio,
                scalar_photon_bend=scalar_bend,directional_photon_bend=vector_bend,photon_fraction=vector_bend/scalar_bend,
                max_abs_U=maxu,max_abs_beta=maxbeta,weak_amplitude_flag=bool(max(maxu,maxbeta)>.01),
                inferred_ratio_for_unit_star_effect=ratio/np.sqrt(abs(star_ratio)),
                inferred_photon_fraction_at_unit_star_effect=vector_bend/scalar_bend/abs(star_ratio))


def main():
    out=HERE/'slow-motion-v1';out.mkdir(exist_ok=False)
    files=[HERE/name for name in ('slow-motion-protocol.md','slow_motion.py','run_audit.py')]
    save(out/'manifest.json',dict(git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                                hashes={p.relative_to(ROOT).as_posix():digest(p) for p in files}))
    rows=[];refinements=[]
    for speed in (.2,.02,.002,200/299792.458):
        for sense in (-1,1):
            for ratio in (1,10,100,1000):
                for radius in (2,4,8):
                    row=evaluate(2048,speed,sense,ratio,radius)
                    fine=evaluate(4096,speed,sense,ratio,radius)
                    keys=('scalar_star_acceleration','directional_star_acceleration','scalar_photon_bend','directional_photon_bend')
                    error=max(abs(row[k]-fine[k])/max(abs(fine[k]),1e-30) for k in keys)
                    rows.append(row);refinements.append(dict(index=len(rows)-1,relative_error=error,passed=bool(error<=1e-6)))
    quadratures=[]
    speed=200/299792.458
    for sense in (-1,1):
        for ratio in (1,1000):
            pos,vel,energy=ring(2048,speed,sense);impact=2*R;eta=G*ratio
            def integrand(z,component):
                value=field(np.array([z,impact,0]),pos,vel,energy,eta)
                return -value[1][1] if component==0 else value[3][2]
            actual=[quad(lambda z:integrand(z,j),-np.inf,np.inf,epsabs=1e-20,epsrel=1e-10,limit=200)[0] for j in (0,1)]
            row=next(r for r in rows if r['speed']==speed and r['sense']==sense and r['coupling_ratio']==ratio and r['radius_ratio']==2)
            expected=[row['scalar_photon_bend'],row['directional_photon_bend']]
            error=max(abs(a-b)/abs(b) for a,b in zip(actual,expected))
            quadratures.append(dict(sense=sense,ratio=ratio,actual=actual,reference=expected,relative_error=error,passed=bool(error<=1e-7)))
    summary=dict(cases=len(rows),refinement_passes=sum(r['passed'] for r in refinements),
                 max_refinement_error=max(r['relative_error'] for r in refinements),
                 quadrature_passes=sum(r['passed'] for r in quadratures),max_quadrature_error=max(r['relative_error'] for r in quadratures),
                 weak_amplitude_flag_count=sum(r['weak_amplitude_flag'] for r in rows),
                 illustrative=[r for r in rows if r['speed']==speed and r['sense']==1 and r['coupling_ratio']==1])
    save(out/'cases.json',rows);save(out/'refinements.json',refinements);save(out/'quadratures.json',quadratures);save(out/'summary.json',summary)
    save(out/'hashes.json',{p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file()})
    print((out/'summary.json').read_text(encoding='utf8'))


if __name__=='__main__':main()
