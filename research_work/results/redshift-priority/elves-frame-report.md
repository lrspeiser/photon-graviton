# Direction-only frame audit before evaluation

The ELVES source describes querying SIMBAD for redshift information. [SIMBAD's documentation](https://simbad.cds.unistra.fr/Pages/guide/ch15.htx) distinguishes radial velocities and redshifts and describes solar-relative velocity information; its object-page help labels the frame heliocentric. The ELVES machine-readable table simply labels its column radial velocity. This is insufficient to certify every original measurement's convention, transformation and epoch. We therefore compute a conditional frame diagnostic, not a new correction to the data.

## Declared external input and meaning

The diagnostic adopts the [Planck 2018 results I](https://doi.org/10.1051/0004-6361/201833880) solar dipole speed 369.82 km/s and Galactic apex (264.021, 48.253) degrees. The transformation uses a kinematic interpretation of that dipole and local special relativity. It does not adopt the paper's expanding cosmology, galaxy peculiar-velocity model, dark matter or Big-Bang explanation. Whether the fictional theory retains this kinematic interpretation is a separate assumption that must remain visible.

Established Lorentz boost and spectral-factor definitions, applied to the declared frames; not a new formula:

F_solar_to_CMB = 1/[gamma (1 - beta dot n_solar)],

1 + z_CMB = F_solar_to_CMB (1 + z_solar),

gamma = 1/sqrt(1 - beta^2).

n_solar points toward the source in the solar comparison frame; the received photon propagates in direction -n_solar. This convention matters: using a CMB-frame direction would require accounting for aberration instead of inserting it unchanged. Astropy 7.1.1 converts the adopted Galactic apex to ICRS; the tabulated J2000 positions are treated as solar/barycentric astrometric directions for this diagnostic. Exact solar-versus-barycentric and catalog epoch details remain part of the eventual precision audit.

The sign can be checked directly: looking toward solar motion, an observer receives a blueshift relative to the CMB frame, so the CMB redshift factor is larger. Toward the opposite direction it is smaller. Transversely the factor is 1/gamma with the stated solar-frame direction. The code verifies all three limits and, for each of 29 directions, independently boosts a photon four-vector and reverses the boost.

## Result without reading candidate outcomes

All 29 originally staged targets are retained with their current decisions; 26 are still pending. Across those pending directions, c times (F-1) ranges from -224.924 to +367.107 km/s. This is a convenient zero-redshift equivalent, not an inferred velocity or the actual correction for a particular galaxy. The actual redshift change also multiplies the unknown factor (1+z_solar).

The size is material relative to the current roughly 450 km/s exploratory model scatter. It does not explain that scatter: the existing CF4 data already carry a CMB-frame label, and applying this factor again would be a double correction. Nor does it measure any galaxy's independent motion. No candidate distance or velocity values were read by the script; no spectral factor was applied to a catalog outcome and no model parameter was fitted.

Before evaluation, recover the source measurement and conversion convention for each usable ELVES velocity, determine whether a heliocentric-to-CMB transform is actually needed, and quantify uncertainty in the adopted frame assumptions. Approximate velocity-to-redshift conversion must not silently replace the exact convention already documented in observation-model.md. Continue to retain motion, endpoint and model uncertainties separately; the accurately measured solar dipole does not supply their missing distributions.

## Two unresolved group aliases

Further HyperLEDA requests for PGC 5808772 and PGC 4231240 timed out. CDS Sesame resolves AGC 740112 to SDSS J104955.41+230406.7 and resolves LV J1017+2922 to its named entry with object-type code PofG, but the queried PGC identifiers return no match. The saved metadata adds independent positions and source names; it does not explicitly certify either PGC alias or settle host/group identity. Both candidates remain withheld pending confirmation, and the source's object classification also needs interpretation before treating the latter as an independent galaxy. A failed lookup is not evidence that a historical overlap is absent.

The updated frame work changes the next action: source-level spectral provenance must precede any numerical new-sample evaluation. The final physical/common model and credible uncertainty still need to be frozen; 26 pending identities are not 26 certified fresh test objects.
