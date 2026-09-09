"""One-zone energy transport: loss, escape and capture are separate processes.

This is a conditional rate model, not a microscopic or spatial halo solution.
Escaped energy remains explicitly accounted for outside the chosen volume.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.linalg import expm


@dataclass(frozen=True)
class TransportRates:
    photon_conversion: float
    companion_capture: float
    photon_escape: float = 0.
    companion_escape: float = 0.
    companion_energy_loss: float = 0.
    deposit_energy_loss: float = 0.

    def validate(self):
        if any(not np.isfinite(v) or v < 0 for v in asdict(self).values()):
            raise ValueError('Rates must be finite and nonnegative.')


def operator(rates, source):
    rates.validate()
    h,g,ep,ec,lc,ld = asdict(rates).values()
    # State: photon, companion, deposit, receiving sector, escaped, fuel, 1.
    a=np.zeros((7,7))
    a[0,0],a[0,6]=-(h+ep),source
    a[1,0],a[1,1]=h,-(g+ec+lc)
    a[2,1],a[2,2]=g,-ld
    a[3,1],a[3,2]=lc,ld
    a[4,0],a[4,1]=ep,ec
    a[5,6]=-source
    return a


def solve(times,rates,initial_photons=1.,fuel=0.,source=0.,duration=0.):
    times=np.asarray(times,dtype=float)
    if times.ndim!=1 or not len(times) or not np.isfinite(times).all() or np.any(times<0) or np.any(np.diff(times)<0):
        raise ValueError('Times must be finite, nonnegative and ordered.')
    if any(not np.isfinite(x) or x<0 for x in [initial_photons,fuel,source,duration]):
        raise ValueError('Invalid initial/source parameters.')
    stop=min(duration,fuel/source) if source else 0.
    on,off=operator(rates,source),operator(rates,0.)
    initial=np.array([initial_photons,0.,0.,0.,0.,fuel,1.])
    at_stop=expm(on*stop)@initial
    return np.array([(expm(on*t)@initial if t<=stop else expm(off*(t-stop))@at_stop)[:6] for t in times])


def permanent_capture_fraction(rates):
    """Infinite-time impulse yield for permanent deposits and finite fixed rates."""
    rates.validate()
    if rates.deposit_energy_loss:
        raise ValueError('Permanent yield requires zero deposit leakage.')
    h,g,ep,ec,lc,_=asdict(rates).values()
    if not h or not g:
        return 0.
    return h/(h+ep)*g/(g+ec+lc)
