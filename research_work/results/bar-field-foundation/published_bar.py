"""Bar density adapted from AGAMA, Mattia Sormani and Eugene Vasiliev.
Source commit f302756b8af2b763db58e278e30478517dc8eea3,
py/example_mw_bar_potential.py. See AGAMA-LICENSE.txt.
Only density functions are retained; no halo or disk is imported.
Density units Msun/kpc^3, Cartesian input in kpc. Known published model.
"""
import numpy
def makeXBar(**params):
    densityNorm = params['densityNorm']
    x0   = params['x0']
    y0   = params['y0']
    z0   = params['z0']
    xc   = params['xc']
    yc   = params['yc']
    c    = params['c']
    alpha= params['alpha']
    cpar = params['cpar']
    cperp= params['cperp']
    m    = params['m']
    n    = params['n']
    outerCutoffRadius = params['outerCutoffRadius']
    def density(xyz):
        r  = numpy.sum(xyz**2, axis=1)**0.5
        a  = ( ( (abs(xyz[:,0]) / x0)**cperp + (abs(xyz[:,1]) / y0)**cperp )**(cpar/cperp) +
            (abs(xyz[:,2]) / z0)**cpar )**(1/cpar)
        ap = ( ((xyz[:,0] + c * xyz[:,2]) / xc)**2 + (xyz[:,1] / yc)**2 )**(0.5)
        am = ( ((xyz[:,0] - c * xyz[:,2]) / xc)**2 + (xyz[:,1] / yc)**2 )**(0.5)
        return (densityNorm / numpy.cosh(a**m) * numpy.exp( -(r/outerCutoffRadius)**2) *
            (1 + alpha * (numpy.exp(-ap**n) + numpy.exp(-am**n) ) ) )
    return density

def makeLongBar(**params):
    densityNorm = params['densityNorm']
    x0   = params['x0']
    y0   = params['y0']
    cpar = params['cpar']
    cperp= params['cperp']
    scaleHeight = params['scaleHeight']
    innerCutoffRadius   = params['innerCutoffRadius']
    outerCutoffRadius   = params['outerCutoffRadius']
    innerCutoffStrength = params['innerCutoffStrength']
    outerCutoffStrength = params['outerCutoffStrength']
    def density(xyz):
        R = (xyz[:,0]**2 + xyz[:,1]**2)**0.5
        a = ( (abs(xyz[:,0]) / x0)**cperp + (abs(xyz[:,1]) / y0)**cperp )**(1/cperp)
        return densityNorm / numpy.cosh(xyz[:,2] / scaleHeight)**2 * numpy.exp(-a**cpar
            -(R/outerCutoffRadius)**outerCutoffStrength - (innerCutoffRadius/R)**innerCutoffStrength)
    return density

def makeBarDensity():
    params = numpy.array(
      # short/thick bar
      [ 3.16273226e+09, 4.90209137e-01, 3.92017253e-01, 2.29482096e-01,
        1.99110223e+00, 2.23179266e+00, 8.73227940e-01, 4.36983774e+00,
        6.25670015e-01, 1.34152138e+00, 1.94025114e+00, 7.50504078e-01,
        4.68875471e-01] +
      # long bar 1
      [ 4.95381575e+08, 5.36363324e+00, 9.58522229e-01, 6.10542494e-01,
        9.69645220e-01, 3.05125124e+00, 3.19043585e+00, 5.58255674e-01,
        1.67310332e+01, 3.19575493e+00] +
      # long bar 2
      [ 1.74304936e+13, 4.77961423e-01, 2.66853061e-01, 2.51516920e-01,
        1.87882599e+00, 9.80136710e-01, 2.20415408e+00, 7.60708626e+00,
       -2.72907665e+01, 1.62966434e+00]
    )
    ind=0
    densityXBar = makeXBar(
         densityNorm=params[ind+0],
         x0=params[ind+1],
         y0=params[ind+2],
         z0=params[ind+3],
         cpar=params[ind+4],
         cperp=params[ind+5],
         m=params[ind+6],
         outerCutoffRadius=params[ind+7],
         alpha=params[ind+8],
         c=params[ind+9],
         n=params[ind+10],
         xc=params[ind+11],
         yc=params[ind+12])
    ind+=13
    densityLongBar1 = makeLongBar(
        densityNorm=params[ind+0],
        x0=params[ind+1],
        y0=params[ind+2],
        scaleHeight=params[ind+3],
        cperp=params[ind+4],
        cpar=params[ind+5],
        outerCutoffRadius=params[ind+6],
        innerCutoffRadius=params[ind+7],
        outerCutoffStrength=params[ind+8],
        innerCutoffStrength=params[ind+9] )
    ind+=10
    densityLongBar2 = makeLongBar(
        densityNorm=params[ind+0],
        x0=params[ind+1],
        y0=params[ind+2],
        scaleHeight=params[ind+3],
        cperp=params[ind+4],
        cpar=params[ind+5],
        outerCutoffRadius=params[ind+6],
        innerCutoffRadius=params[ind+7],
        outerCutoffStrength=params[ind+8],
        innerCutoffStrength=params[ind+9] )
    ind+=10
    assert len(params)==ind, 'invalid number of parameters'
    return lambda x: densityXBar(x) + densityLongBar1(x) + densityLongBar2(x)
