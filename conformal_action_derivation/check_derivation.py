"""Numerical checks of algebraic identities, not new observational validation."""
from pathlib import Path
import csv,json,math
P=Path(__file__).resolve().parent
kappa=.000077315 # inverse million light-years, fixed from previous work
pc=149597870700.*648000/math.pi
mpc_in_mly=pc/(299792458.*31557600.)
kap=kappa*mpc_in_mly
Rc=4547.068376184085 # Mpc, previous exploratory fit; not fitted here
rows=[]
for z in [.1,.5,1.,2.]:
    S=1+z;ne=1/S;no=1.
    # q=-dt^2/n^2+spatial metric. Local matter clock rate per t = B/n.
    # Emitted photon coordinate frequency = B_e/n_e in units of nu_atom=1.
    # Conserved spatial momentum fixes received coordinate frequency.
    def measured_redshift(Be,Bo):
        emitted=Be/ne
        received=emitted*ne/no
        laboratory=Bo/no
        return laboratory/received-1
    null=measured_redshift(1.,1.)
    recovered=measured_redshift(ne,no)
    assert abs(null)<1e-14
    assert abs(recovered-z)<1e-14
    R=math.log1p(z)/kap
    D=Rc*math.sinh(R/Rc)
    DLstatic=S*D;DLmetric=S*D
    DAstatic=D;DAmetric=D/S
    assert abs(DLmetric/DAmetric-S*S)<1e-12
    assert abs(DLstatic/DAstatic-S)<1e-12
    rows.append(dict(z=z,universal_static_redshift=null,
        conformal_material_expansion_redshift=recovered,
        current_geometric_distance_Mpc=R,
        static_transport_DL_Mpc=DLstatic,metric_coasting_DL_Mpc=DLmetric,
        static_transport_DA_Mpc=DAstatic,metric_coasting_DA_Mpc=DAmetric))
with (P/'distance_comparison.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(P/'checks.json').write_text(json.dumps(dict(status='Algebra verified numerically; not a new fit or successful nonexpanding model',
    kappa_per_Mly=kappa,previous_negative_radius_Mpc=Rc,rows=rows),indent=2)+'\n')
print(json.dumps(rows,indent=2))
