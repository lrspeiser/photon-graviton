"""Test uniform baryon-attached loading against the existing effective source."""
from pathlib import Path
import ast
import hashlib
import importlib.util
import json
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CODE=HERE.parent/'full-bar-source-audit/run.py'
spec=importlib.util.spec_from_file_location('source_audit',CODE)
source=importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)
from field import bar, nuclei
SOURCE=HERE.parent/'full-bar-source-audit/results.json'
DISK=HERE.parent/'conservative-field-completion/run.py'
FIELD=HERE.parent/'bar-field-foundation/field.py'
BAR=HERE.parent/'bar-field-foundation/published_bar.py'

def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()

def main():
    files=[Path(__file__),CODE,SOURCE,DISK,FIELD,BAR]
    hashes={str(p.relative_to(ROOT)):digest(p) for p in files}
    data=json.loads(SOURCE.read_text())
    for path,expected in data['source_hashes'].items():
        assert digest(ROOT/path)==expected
    rows=data['rows']
    xyz=np.array([r['xyz_kpc'] for r in rows])
    R=np.hypot(xyz[:,0],xyz[:,1]);z=xyz[:,2]
    # Read the declared disk parameters from the actual cached-field builder.
    tree=ast.parse(DISK.read_text())
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='Disk')
    init=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
    assignment=next(n for n in init.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='params' for t in n.targets))
    params=ast.literal_eval(assignment.value)
    components={'bar':bar(xyz),'nuclei':nuclei(xyz)}
    for i,(sigma,rd,height,hole,kind) in enumerate(params):
        surface=sigma*np.exp(-hole/R-R/rd)
        vertical=(np.exp(-abs(z)/height)/(2*height) if kind=='exp' else
                  np.exp(-2*np.logaddexp(z/(2*height),-z/(2*height))+2*np.log(2))/(4*height))
        components[f'disk_{i}']=surface*vertical
    components['softened_centre']=3*4.1e6*.001**2/(4*np.pi*(np.sum(xyz*xyz,axis=1)+.001**2)**2.5)
    baryon=sum(components.values())
    assert np.isfinite(baryon).all() and (baryon>0).all()
    targets={
        'finer':np.array([r['rho_finer_Msun_kpc3'] for r in rows]),
        'fine':np.array([r['rho_fine_Msun_kpc3'] for r in rows]),
        'requested_divQ':np.array([r['rho_from_requested_divQ_steps'][-1] for r in rows]),
    }
    summaries={}
    for label,target in targets.items():
        assert (target>0).all()
        q=target/baryon
        lo,hi=float(q.min()),float(q.max())
        eta=2*lo*hi/(lo+hi)
        error=(hi-lo)/(hi+lo)
        actual=np.max(abs(eta/q-1))
        assert abs(actual-error)<1e-12
        # Independent scalar optimizer of the convex maximum-relative-error loss.
        from scipy.optimize import minimize_scalar
        opt=minimize_scalar(lambda value:np.max(abs(value/q-1)),bounds=(lo,hi),method='bounded',options={'xatol':1e-12})
        assert opt.success and abs(opt.fun-error)<1e-7
        summaries[label]=dict(required_loading_minimum=lo,required_loading_maximum=hi,
            required_loading_dynamic_range=hi/lo,minimax_relative_error=error,
            minimax_loading=eta,minimax_multiplicative_factor=float(np.sqrt(hi/lo)),
            minimum_location={k:rows[int(q.argmin())][k] for k in ['R_kpc','z_kpc','phi_rad']},
            maximum_location={k:rows[int(q.argmax())][k] for k in ['R_kpc','z_kpc','phi_rad']})
    result_rows=[]
    for i,r in enumerate(rows):
        result_rows.append(dict(R_kpc=r['R_kpc'],z_kpc=r['z_kpc'],phi_rad=r['phi_rad'],
            declared_baryon_density_Msun_kpc3=float(baryon[i]),
            equivalent_extra_density_Msun_kpc3=float(targets['finer'][i]),
            required_extra_per_baryon=float(targets['finer'][i]/baryon[i]),
            components_Msun_kpc3={name:float(values[i]) for name,values in components.items()}))
    for path in files:
        assert digest(path)==hashes[str(path.relative_to(ROOT))]
    out=dict(scope='Shape test of constant loading on declared baryonic matter, not an observational fit or a rejection of all deposits.',
        postulate='rho_extra(x)=eta*rho_b(x), one nonnegative position-independent eta',
        assumptions=['Deposits remain attached to ordinary matter with identical loading per baryonic mass.',
                     'Extra source obeys Newtonian Poisson sourcing; no additional response kernel.',
                     'Uniform accumulated capture per baryonic mass, or equal saturated capacity per baryonic mass.'],
        samples=len(rows),summaries=summaries,rows=result_rows,disk_parameters=params,
        source_hashes=hashes,holdouts_opened=False,energy_supply_budget_evaluated=False,
        nonuniform_capture_ruled_out=False,modified_gravity_response_ruled_out=False)
    (HERE/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(summaries,indent=2))

if __name__=='__main__':
    main()
