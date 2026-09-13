# Angular escape channels and accumulated direction changes

13 September 2026. Local geometric diagnostic; not a radiative-transfer or image fit.

## Finding

Using less illuminated outgoing directions can favor the proposed loading interaction, but it also changes photon directions. The simplest isolated-source geometry does not force an unacceptable direction spread in every case: a very distant star occupies a tiny angular patch, so a small deflection can reach outside it. This does not establish a working cosmological mechanism. The surrounding sky, scattering kernel and filling of outgoing modes must be calculated.

## Geometry and formula provenance

**Hypothesis:** one uniformly bright circular source occupies an angular disk of radius alpha. The same disk is illuminated at both closely spaced frequencies. Inside it each paired mode has occupation n; outside it the relevant outgoing modes initially have zero occupation. Every interaction changes direction by a fixed angle theta with uniformly random azimuth. We use flat angular-plane geometry, appropriate to the small apparent source sizes tested.

Let d=theta/alpha. The fraction f of outgoing directions falling outside the bright disk, averaged over uniform incoming directions, is one minus the overlap area of two translated equal disks:

\[
f(d)=1-\frac{2\arccos(d/2)-d\sqrt{1-d^2/4}}\pi,\quad 0\le d\le2,
\]

and f=1 for d>=2. This is known circle-intersection geometry, derived here by integrating intersecting chords. The physical choice of this source and scattering kernel is ours.

Averaging the previously proposed bosonic rates over inside and outside channels gives

\[
a=\Gamma n[1+(1-f)n],\quad b=\Gamma(1-f)n(1+n),
\]
\[
r=\frac{(1-f)(1+n)}{1+(1-f)n},\quad
f=\frac{(1-r)(1+n)}{1+(1-r)n}.
\]

These equations assume equal forward/reverse matrix weights and equal occupations within the illuminated patches. They describe this paired-channel averaging, not a full angular collision operator. At low occupation, r=0.9 requires f approximately 0.1, which corresponds to theta=0.157242 alpha. At high occupation a larger empty fraction is needed. Exact thermal frequency differences are neglected for this diagnostic.

## Repeated-kick comparison

For independent isotropic kicks of identical size, the directional correlation after N kicks is (cos theta)^N. In the small-spread limit the directional RMS is approximately sqrt(N) theta. This is random-walk mathematics; it is not a new redshift law or the angular width of an observed image.

For continuity with earlier calculations, use a 2 eV initial photon, fixed transfer delta=1e-8 eV and the energy-loss fraction stored in quantum-bridge-results.json. The net number of forward energy increments needed is N=1,520,448.73. We use that mean-count scale as an optimistic kick-count diagnostic. Reverse events add traffic; correlated deflections, selection of rays reaching the observer and varying scattering locations require a separate calculation. For a Poisson kick count the correlation is exp[N(cos theta-1)], rather than the fixed-count expression; the leading small-angle RMS is the same.

Assume a Sun-sized emitting disk of radius 6.957e8 m, dilute modes and r=0.9, holding its apparent size fixed throughout this diagnostic:

| Distance from source to interaction region | Required kick size (arcseconds) | Small-angle accumulated directional RMS (arcseconds) |
|---:|---:|---:|
| 1 parsec | 0.00073125 | 0.90168 |
| 1,000 parsecs | 7.3125e-7 | 0.00090168 |
| 1,000,000 parsecs | 7.3125e-10 | 9.0168e-7 |

At one astronomical unit the required kick is about 150.83 arcseconds and the cumulative spread is not small; the exact fixed-kick directional correlation is about 0.666. We do not report a small-angle RMS there. A chosen one-arcsecond RMS comparison corresponds to about 0.902 parsecs in this fixed-geometry example. One arcsecond is an illustrative benchmark, not a measured universal tolerance.

## What the favorable numbers do not establish

- Source distance here means source-to-interaction distance, not source-to-observer distance. Along a real ray the angular source size changes; interactions near the emitter may contribute disproportionate deflection.
- A galaxy contains many emitters, gas and diffuse backgrounds. Its angular occupation pattern is not one isolated uniformly bright disk with empty surroundings.
- Outgoing channels accumulate photons unless transport empties them. The previous source/outlet calculation must be connected to these actual modes before f can be treated as a maintained empty fraction.
- A physical interaction must produce the extremely narrow angular kernel without losing its required rate, and account for momentum recoil. The earlier spatial-response calculation already showed that narrowing a kernel can strongly suppress its total rate.
- Random directional spread is not an observer's point-spread function: image position depends on the locations of kicks and which rays reach the telescope. Reverse kicks and correlations must also be included.

Consequently this result keeps angular mode separation open as a candidate, but supplies no fitted redshift or lensing prediction. A one-dimensional empty-channel model cannot claim image compatibility from these numbers alone. Polarization-separated or companion-only channels remain distinct alternatives requiring their own observational consequences. The exact one-third reference is unchanged.

## Verification

Thirty-six configurations cover three source occupations, three target rate ratios and four source distances. Nine independently integrated overlap areas verify the analytic circle formula and inversion to tolerance 1e-8. The JSON retains one minus directional correlation using expm1 so tiny spreads are not silently rounded away. Inputs and executable are local; no observational data are fitted in this experiment.

Files: `angular-outlet.py`, `angular-outlet-results.json`. Related: [source/outlet](source-outlet-report.md), [spatial response](spatial-response-report.md).
