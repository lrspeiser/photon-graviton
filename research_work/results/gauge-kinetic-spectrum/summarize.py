"""Check carrier calculations and plot the resonant outgoing spectrum."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
base=json.loads((HERE/'results.json').read_text())
res=json.loads((HERE/'resonant-results.json').read_text())
coarse,fine=res['cases']
checks={
    'resonant_scalar_fraction_absolute_refinement_change':abs(coarse['final_scalar_fraction']-fine['final_scalar_fraction']),
    'resonant_mean_ratio_absolute_refinement_change':abs(coarse['mean_wavenumber_ratio']-fine['mean_wavenumber_ratio']),
    'resonant_outgoing_power_absolute_refinement_change':abs(coarse['outgoing_field_power_ratio']-fine['outgoing_field_power_ratio']),
    'resonant_scalar_fraction':fine['final_scalar_fraction'],
    'resonant_energy_weighted_mean_frequency_ratio':fine['mean_wavenumber_ratio'],
    'resonant_wave_action_proxy_ratio':fine['outgoing_wave_action_proxy_ratio'],
    'resonant_power_over_action_ratio':fine['outgoing_field_power_ratio']/fine['outgoing_wave_action_proxy_ratio'],
    'resonant_scalar_plus_windowed_outgoing_power_minus_one':fine['final_scalar_fraction']+fine['outgoing_field_power_ratio']-1,
    'max_energy_drift_all_cases':max(r['max_total_energy_relative_drift'] for r in base['cases']+res['cases']),
}
assert checks['resonant_scalar_fraction_absolute_refinement_change']<1e-8
assert checks['resonant_mean_ratio_absolute_refinement_change']<1e-8
assert fine['field_weighted_kinetic_coefficient_departure']<1e-8
assert abs(checks['resonant_scalar_plus_windowed_outgoing_power_minus_one'])<1e-6
assert checks['max_energy_drift_all_cases']<1e-6
sp=res['spectrum']; k=np.array(sp['k'])
p0=np.array(sp['initial_normalized_power']); p1=np.array(sp['final_normalized_power'])
assert abs(p0.sum()-1)<1e-12 and abs(p1.sum()-1)<1e-12
recomputed=float((k*p1).sum()/(k*p0).sum())
assert abs(recomputed-fine['mean_wavenumber_ratio'])<1e-12
checks['spectral_mean_independently_recomputed']=recomputed
(HERE/'verification.json').write_text(json.dumps(checks,indent=2)+'\n',newline='\n')

dk=k[1]-k[0]
fraction=1-fine['final_scalar_fraction']
target=np.interp(k/fraction,k,p0,left=0,right=0)/fraction
target/=target.sum()
fig,axes=plt.subplots(1,2,figsize=(11,4.4))
axes[0].plot(k,p0/dk,label='Freely propagated input',color='#687787')
axes[0].plot(k,p1/dk,label='Calculated outgoing light',color='#067a87')
axes[0].plot(k,target/dk,'--',label='Uniform redshift if loss alone set it',color='#bf5328')
axes[0].set(xlim=(.2,2.1),xlabel='Wavenumber (dimensionless)',ylabel='Normalized field-power density')
axes[0].legend(fontsize=8)
axes[0].set_title('The surviving spectrum is not uniformly redshifted')
axes[1].bar(['Initial light','Outgoing light','Scalar field'],[100,100*fine['outgoing_field_power_ratio'],100*fine['final_scalar_fraction']],color=['#687787','#067a87','#9363ac'])
axes[1].set(ylabel='Percent of initial total energy',ylim=(0,112))
axes[1].set_title('Energy transfers, but that does not set frequency')
for i,v in enumerate([100,100*fine['outgoing_field_power_ratio'],100*fine['final_scalar_fraction']]):
    axes[1].text(i,v+2,f'{v:.2f}%',ha='center')
fig.suptitle('Closed photon–scalar collision: resonant dimensionless example',fontsize=13)
fig.text(.5,.01,'Outgoing energy uses symmetric packets and a localized free-field diagnostic; agreement with the full energy ledger is within 0.00002 percentage points.',ha='center',fontsize=8)
fig.tight_layout(rect=(0,.04,1,.95))
fig.savefig(HERE/'spectrum-comparison.png',dpi=170)
plt.close(fig)
print(json.dumps(checks,indent=2))
