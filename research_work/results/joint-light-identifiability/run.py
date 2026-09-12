"""Conditional joint-light invariance and local DES metadata audit."""
from pathlib import Path
import hashlib,json
import numpy as np
from astropy.io import fits
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def main():
    audit=json.loads((HERE.parent/'timing-foundation/input-audit.json').read_text(encoding='utf-8'))
    relevant=[v for v in audit['input_manifest'] if v['path'].endswith('_HEAD.FITS.gz') or v['path'].endswith('README.md') or v['path'].endswith('_DES.README')]
    for v in relevant:assert hashlib.sha256((ROOT/v['cache_relative_path']).read_bytes()).hexdigest()==v['sha256']
    head=next(v for v in relevant if v['path'].endswith('_HEAD.FITS.gz'))
    with fits.open(ROOT/head['cache_relative_path']) as h:columns=list(h[1].columns.names)
    z=np.array([.001,.01,.1,.5,1.2]);S=1+z;D=np.array([10.,30.,100.,500.,1000.])
    t=np.linspace(-300,300,1001)[:,None]
    def predict(b,ew,el):
        width=30*S**(b+ew)
        peak=2*S**(el-1-b)/(4*np.pi*D**2)
        flux=peak[None,:]*np.exp(-.5*(t/width[None,:])**2)
        fluence=peak*width*np.sqrt(2*np.pi)
        return width,peak,flux,fluence
    b,ew,el=1.,.12,-.07;reference=predict(b,ew,el);checks=[]
    for delta in [-2.,-1.,-.3,.2,1.,2.]:
        changed=predict(b+delta,ew-delta,el+delta)
        errors=[float(np.max(abs(a-v))/np.max(abs(a))) for a,v in zip(reference,changed)]
        assert max(errors)<1e-12
        checks.append({'delta':delta,'relative_max_errors_width_peak_lightcurve_fluence':errors})
    jacobian=np.array([[1.,1.,0.],[-1.,0.,1.]])
    null=np.array([1.,-1.,1.]);assert np.array_equal(jacobian@null,np.zeros(2))
    ranks={'width_and_flux':int(np.linalg.matrix_rank(jacobian)),
       'plus_fluence':int(np.linalg.matrix_rank(np.vstack([jacobian,[0,1,1]]))),
       'plus_independent_duration_evolution_constraint':int(np.linalg.matrix_rank(np.vstack([jacobian,[0,1,0]]))),
       'plus_independent_peak_luminosity_evolution_constraint':int(np.linalg.matrix_rank(np.vstack([jacobian,[0,0,1]])))}
    assert ranks=={'width_and_flux':2,'plus_fluence':2,'plus_independent_duration_evolution_constraint':3,'plus_independent_peak_luminosity_evolution_constraint':3}
    out={'scope':'Algebraic identifiability in a declared static-geometric, bolometric model with power-law source evolution. No data fit or new causal mechanism.', 'checks':checks,'jacobian':jacobian.tolist(),'null_direction':null.tolist(),'ranks':ranks,'head_columns':columns,'verified_inputs':relevant,'real_flux_used':False,'new_independent_distances_acquired':False}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'ranks':ranks,'maximum_invariance_error':max(max(v['relative_max_errors_width_peak_lightcurve_fluence']) for v in checks),'head_column_count':len(columns)}))
if __name__=='__main__':main()
