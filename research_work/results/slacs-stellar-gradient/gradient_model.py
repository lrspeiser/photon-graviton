"""Positive, bounded radial mass/light sensitivity; observed tracer light stays fixed."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from model import ComponentModel, G

class GradientModel(ComponentModel):
    def __init__(self, *args, eta=0., **kwargs):
        super().__init__(*args, **kwargs)
        if eta == 0:
            return
        r = self.r
        # Fixed angular pivot from the original published Re, not a fitted scale.
        pivot = 1.8153*self.a
        shape = 1 + eta/(1+(r/pivot)**2)
        assert np.all(shape > 0)
        density = self.nu*shape
        slope = np.log(density[1]/density[0])/np.log(r[1]/r[0])
        initial = 4*np.pi*density[0]*r[0]**3/(3+slope)
        enclosed = initial+4*np.pi*cumulative_trapezoid(r*r*density, r, initial=0)
        fraction = enclosed/enclosed[-1]
        assert np.all(np.diff(fraction) >= 0)
        self.mass_interp = PchipInterpolator(np.log(r), fraction, extrapolate=False)
        self.inner_power = 3+slope
        gb = G*1e11*fraction/r**2
        gct = self.A*self.astar*(G*1e11*self.mass_fraction(self.rt)/self.rt**2/self.astar)**self.p
        self.forces = np.array([gb, np.where(r <= self.rt, self.A*self.astar*(gb/self.astar)**self.p,
                                            gct*(self.rt/r)**2)])
        # self.nu, aperture light weights, and self.den remain the measured tracer.
