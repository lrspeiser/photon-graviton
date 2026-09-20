"""Verify RB-1 archive, physical units and active coefficient against SR-1."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from decimal import Decimal as D
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from spatial_model import local

ROOT=Path(__file__).resolve().parent


def main():
    result=json.loads((ROOT/'radiation-budget-v1/results.json').read_text())
    manifest=json.loads((ROOT/'radiation-budget-v1/manifest.json').read_text())
    checks=[]
    def check(name,value):checks.append(dict(name=name,passed=bool(value)))
    for name,digest in manifest['sources'].items():
        check('current '+name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
        blob=subprocess.check_output(['git','show',manifest['source_commit']+':research_work/experiments/common_cone/'+name],cwd=ROOT)
        check('committed '+name,hashlib.sha256(blob).hexdigest()==digest)
    check('counts',len(result['budgets'])==54 and len(result['profiles'])==36)
    c=D('299792458');G=D('6.67430e-11');year=D('31557600');target=D('1e10')*year
    for i,row in enumerate(result['budgets']):
        v=D(str(row['speed_m_per_s']));mass=D(str(row['source_kg']));eps=D(str(row['available_fraction']));active=D(1+row['s'])
        lum=v*v*c*c*c/(G*active);life=eps*mass*c*c/lum
        required=v*v*c*target/(eps*G*mass)
        check('budget '+str(i),abs(row['luminosity_watt']/float(lum)-1)<1e-12 and abs(row['lifetime_years']/float(life/year)-1)<1e-12 and abs(row['required_active_coefficient']/float(required)-1)<1e-12 and row['target_lifetime_passed']==(life>=target))
    for i,row in enumerate(result['profiles']):
        v=D(str(row['speed_m_per_s']));r=D(str(row['radius_m']));r0=D(str(row['emission_radius_m']));active=D(1+row['s'])
        mass=v*v*(r-r0)/(G*active);velocity=(v*v*(1-r0/r)).sqrt()
        check('profile '+str(i),abs(row['wave_equivalent_kg']/float(mass)-1)<1e-12 and abs(row['circular_speed_m_per_s']/float(velocity)-1)<1e-12 and row['passed'])
    # Direct field-energy derivative from the implemented Hamiltonian, with
    # particle energy canceled and plane-wave kinetic/gradient equipartition.
    derivative_errors=[]
    for s in (0,1):
        for amplitude in (.1,1.,3.):
            pi=amplitude*np.array([.2,.3,.4,.5]);direction=np.array([1.,2.,3.])/np.sqrt(14)
            gradient=direction[:,None]*pi;f=np.zeros(4);p=np.array([.7,.3,.2]);eps=1e-5
            def field(value):
                f[0]=value
                full=local(f,pi,gradient,p,1,s,0,eta=0,omega=0)[0]
                particle=local(f,np.zeros(4),np.zeros((3,4)),p,1,s,0,eta=0,omega=0)[0]
                return full-particle
            numerical=(field(eps)-field(-eps))/(2*eps)
            expected=.08*(1+s)*np.dot(pi,pi)
            error=abs(numerical-expected)/max(1,abs(expected));derivative_errors.append(error)
            check(f'active coefficient s={s}, amplitude={amplitude}',error<2e-8)
    check('physical rejection recorded',result['persistent_budget_pass_count']==0 and not any(r['target_lifetime_passed'] for r in result['budgets']))
    output=dict(passed=all(x['passed'] for x in checks),count=len(checks),max_active_source_derivative_error=max(derivative_errors),checks=checks)
    (ROOT/'radiation-budget-audit.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in output.items() if k!='checks'}))
    if not output['passed']:raise RuntimeError('RB-1 audit failed')


if __name__=='__main__':main()
