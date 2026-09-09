"""Photon-count accounting and a scoped FIRAS residual comparison."""
from pathlib import Path
import hashlib,json,os
import numpy as np
from scipy.constants import h,k,c
from scipy.integrate import quad
from scipy.optimize import least_squares,minimize_scalar

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'thermal-conversion'

def planck(nu,temperature):
    return 2*h*nu**3/c**2/np.expm1(h*nu/(k*temperature))/1e-20 # MJy/sr

def main():
    protocol=json.loads((HERE/'protocol.json').read_text(encoding='utf-8'))
    path=ROOT/protocol['input'];data=np.loadtxt(path)
    assert data.ndim==2 and data.shape[1]==5 and np.all(data[:,3]>0)
    nu=data[:,0]*100*c;observed=data[:,2];sigma=data[:,3];galaxy=data[:,4]
    reference=planck(nu,protocol['reference_temperature_K'])
    weighted_galaxy=galaxy/sigma
    def model(temperature,amplitude):return 1000*(amplitude*planck(nu,temperature)-reference)
    def profile(temperature,amplitude,full=False):
        difference=(observed-model(temperature,amplitude))/sigma
        g=float(np.dot(weighted_galaxy,difference)/np.dot(weighted_galaxy,weighted_galaxy))
        residual=difference-g*weighted_galaxy
        return (float(np.dot(residual,residual)),g,residual) if full else float(np.dot(residual,residual))
    rows=[]
    for q in protocol['energy_survival_factors']:
        amplitude=q**(-3)
        fit=minimize_scalar(lambda temp:profile(temp,amplitude),bounds=(.5,4),method='bounded',options={'xatol':1e-12})
        assert fit.success and .5001<fit.x<3.9999
        chi,g,residual=profile(fit.x,amplitude,True)
        # Cross-check optimizer against an independent nonlinear residual solver.
        repeat=least_squares(lambda x:(model(x[0],amplitude)+x[1]*galaxy-observed)/sigma,
                             [fit.x,g],bounds=([.5,-np.inf],[4,np.inf]),x_scale='jac',
                             ftol=1e-12,xtol=1e-12,gtol=1e-12)
        chi2=float(np.dot(repeat.fun,repeat.fun))
        assert repeat.success and abs(chi2-chi)<1e-5*max(1,chi)
        rows.append({'photon_energy_survival_q':q,'required_blackbody_amplitude':amplitude,
                     'profile_color_temperature_K':float(repeat.x[0]),'Galaxy_template_coefficient':float(repeat.x[1]),
                     'diagonal_chi_square':chi2,'degrees_of_freedom_if_diagonal_model_were_complete':int(len(nu)-2),
                     'maximum_absolute_residual_over_published_sigma':float(np.max(abs(repeat.fun))),
                     'weighted_residuals':repeat.fun.tolist()})
    free=least_squares(lambda x:(model(x[1],np.exp(x[0]))+x[2]*galaxy-observed)/sigma,
                       [0.,2.725,0.],bounds=([-5,.5,-np.inf],[5,4,np.inf]),x_scale='jac',
                       ftol=1e-12,xtol=1e-12,gtol=1e-12)
    assert free.success
    free_chi=float(np.dot(free.fun,free.fun));assert free_chi<=rows[0]['diagonal_chi_square']+1e-5
    # Integrate dimensionless photon-number and energy spectra independently.
    def occupation(x):return 0. if x>700 else 1/np.expm1(x)
    initial_n=quad(lambda x:x*x*occupation(x),0,np.inf)[0]
    initial_u=quad(lambda x:x**3*occupation(x),0,np.inf)[0]
    moments=[]
    for q in [.5,.9,.99]:
        new_n=quad(lambda x:q**-3*x*x*occupation(x/q),0,np.inf)[0]
        new_u=quad(lambda x:q**-3*x**3*occupation(x/q),0,np.inf)[0]
        assert abs(new_n/initial_n-1)<1e-9 and abs(new_u/initial_u-q)<1e-9
        # Uniform photon survival q^3 restores a standard Planck spectrum at qT.
        moments.append({'q':q,'number_ratio':new_n/initial_n,'energy_ratio':new_u/initial_u,
                        'photon_survival_required_for_Planck_restoration':q**3,
                        'total_energy_ratio_after_shift_and_removal':q**4})
    result={'scope':protocol['scope'],'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'protocol_sha256':hashlib.sha256((HERE/'protocol.json').read_bytes()).hexdigest(),
            'number_of_channels':len(nu),'columns_used':[1,3,4,5],
            'number_and_energy_checks':moments,'fixed_conversion_profiles':rows,
            'free_amplitude_comparator':{'amplitude':float(np.exp(free.x[0])),
                                         'color_temperature_K':float(free.x[1]),
                                         'Galaxy_template_coefficient':float(free.x[2]),
                                         'diagonal_chi_square':free_chi},
            'conditional_conclusion':'Number-preserving redshift at constant c in a fixed volume does not preserve standard blackbody normalization. This tests an initially thermal, unreplenished homogeneous bath only.',
            'unresolved':['Full FIRAS covariance, calibration and foreground likelihood',
                          'Companion law at microwave frequencies and actual radiation source history',
                          'Photon emission/absorption, escape, thermalization and receiving-energy ledger',
                          'Background origin and angular correlations without imported cosmology']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'thermal-conversion-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'channels':len(nu),'profiles':[{k:v for k,v in row.items() if k!='weighted_residuals'} for row in rows],
                      'free_comparator':result['free_amplitude_comparator']},indent=2))

if __name__=='__main__':main()
