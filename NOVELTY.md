# What is borrowed, what may be new

**25 September 2026, round 19.** A first novelty audit, requested by the review of 25 September 2026. For every central
claim of the hot-companion work it names the closest published work found, what that work already establishes, and what,
if anything, is left for us. It is a targeted search, not an exhaustive priority review: each reference below was located
and checked (authors, journal, year, main result) on 25 September 2026, but only abstracts and summaries were read, not
the full papers, and the comparison is not yet at the level of equations. The automatic formula guard
(`research_work/tools/formula_guard.py`) checks whether a formula is secretly MOND, Newton or dark matter; it says nothing
about priority.

**Labels.** *Borrowed*: the ingredient is established physics, used as is. *Known class*: the kind of effect is
published; our version differs in a stated way that still needs an equation-level comparison. *New combination*: the
parts are known, their combination was not found. *Not found*: no precedent turned up in this search (which is not the
same as new).

## 1. The law's architecture

| our claim | closest published work | what it establishes | status for us |
|---|---|---|---|
| The field equation ∇²Φ = −∇·h, with h built from the Newtonian field | M. Milgrom, *Quasi-linear formulation of MOND*, MNRAS 403, 886 (2010) | a modified-potential theory in which a nonlinear algebraic step on g_N feeds a linear Poisson equation; derivable from an action, conserving momentum | **borrowed** |
| In the cold limit the law is g = g_N ν(g_N/a) | Milgrom's MOND (1983 onward) | the whole family of MOND interpolating functions and the external-field effect | **borrowed form**; what we add is an interpretation (a = 2ℓ/u, the switch from the release factor), which is only as good as the mechanism below |
| Hot matter weighs more: k = 3σ²/u² | R. C. Tolman (1930) and E. T. Whittaker, *On Gauss' theorem and the concept of mass in general relativity*, Proc. R. Soc. A 149, 384 (1935): the active gravitational mass density is ρ + 3p/c² | for a gas with velocity dispersion σ per axis, p = ρσ², so its active mass is ρ(1 + 3σ²/c²): **the same form as our heat weight, with c where we have u** | **known form.** Our differences: the speed is the companion's u = 169 km/s, not c, so the effect is large; and colliding gas gets k = 0, while in general relativity a gas's pressure counts whether or not its particles collide. Both differences are ours to justify |
| Collisions switch the heat off | R. H. Dicke, *The effect of collisions upon the Doppler width of spectral lines*, Phys. Rev. 89, 472 (1953) | frequent velocity-changing collisions average away motion's effect on emission (Dicke narrowing) | **borrowed physics**; applying it to the companion's emission is ours, and needs the response times and the survival of the internal state under collisions |
| A companion streaming from all matter, whose flux gives the pull | N. Fatio (1690), G.-L. Le Sage (1748): gravity from the shadowing of streams of particles; historical objections summarized, e.g., in the Wikipedia article *Le Sage's theory of gravitation* and in Edwards (ed.), *Pushing Gravity* (2002) | a streaming medium can mediate attraction; the classic objections are **drag** on moving bodies and **heating** by the absorbed stream | **known class, with known objections.** Round 18 met the drag objection head on (waves carried by the stream push every emitter) and the absorbing stream moves it into heating (33–48% of the wave power absorbed); both must be answered quantitatively |

## 2. The machinery inside matter

| our claim | closest published work | what it establishes | status for us |
|---|---|---|---|
| Bodies oscillating in step, coupled through a medium, attract | C. A. Bjerknes (1870s–1880s) and V. Bjerknes, *Fields of Force* (Columbia University Press, 1906): pulsating spheres in a fluid attract when in phase and repel when out of phase, proposed as a hydrodynamic analogy of gravitation; L. A. Crum, *Bjerknes forces on bubbles in a stationary sound field*, J. Acoust. Soc. Am. 57, 1363 (1975) | in-phase emitters attract through the wave they share (secondary Bjerknes force) | **known class.** Ours: the phase is selected by the bodies' own dynamics (an inverted, self-sustained oscillator locks a quarter cycle ahead), not imposed by a driving field |
| An active body is pulled toward the source of the wave it amplifies | A. Mizrahi and Y. Fainman, *Negative radiation pressure on gain medium structures*, Opt. Lett. 35, 3405 (2010); D. L. Gao, R. Shi, Y. Huang, W. H. Ni, L. Gao, *Fano-enhanced pulling and pushing optical force on active plasmonic nanoparticles*, Phys. Rev. A 96, 043826 (2017); J. Chen, J. Ng, Z. Lin, C. T. Chan, *Optical pulling force*, Nat. Photonics 5, 531 (2011) | light amplification can pull a gain body toward the light source; negative (pulling) forces from interference | **known class.** Our round-15 result that a sub-threshold amplifier is pulled in proportion to the wave's energy agrees with this class. What we have not found: an **above-threshold, self-sustained** oscillator locked to a weak passing wave and pulled in proportion to the wave's **height**, the property the law's square root needs. A directed search on forces on injection-locked or lasing particles is owed |
| Gravity-like 1/r² forces mediated by radiation | J. Luis-Hita, M. I. Marqués, R. Delgado-Buscalioni, N. de Sousa, L. S. Froufe-Pérez, F. Scheffold, J. J. Sáenz, *Light induced inverse-square law interactions between nanoparticles: "mock gravity" at the nanoscale*, Phys. Rev. Lett. 123, 143201 (2019); M. M. Burns, J.-M. Fournier, J. A. Golovchenko, *Optical binding*, Phys. Rev. Lett. 63, 1233 (1989) | resonant particles in an isotropic random light field attract with an inverse-square law in near and far field; light binds particles into organized structures | **known class.** Ours is sourced by the matter's own emission, not an external illumination, and aims at a 1/r pull growing with the emitters' motion |
| Collective, self-sustained, narrow-line emitters | D. Meiser, J. Ye, D. R. Carlson, M. J. Holland, *Prospects for a millihertz-linewidth laser*, Phys. Rev. Lett. 102, 163601 (2009) (the superradiant laser) | inverted emitters sharing one mode emit collectively and phase-lock; the Bloch equations we use | **borrowed** |
| Waves in a moving medium, one-way where the flow outruns them, with negative-energy waves | W. G. Unruh, *Experimental black-hole evaporation?*, Phys. Rev. Lett. 46, 1351 (1981); C. Barceló, S. Liberati, M. Visser, *Analogue gravity*, Living Rev. Relativ. 14, 3 (2011); M. V. Nezlin, *Negative-energy waves and the anomalous Doppler effect*, Sov. Phys. Usp. 19, 946 (1976) | the acoustic metric of a flowing fluid, horizons where the flow exceeds the wave speed, negative-energy waves and their instabilities | **borrowed** (round 18's carried-wave medium is this analogue-gravity equation in one dimension) |
| A pressure-driven wind is subsonic inside its source | R. A. Chevalier and A. W. Clegg, *Wind from a starburst galaxy nucleus*, Nature 317, 44 (1985) | the analytic wind from uniform mass and energy injection | **borrowed** (used to rule out a pressure-driven companion as a one-way medium) |
| One-way (non-reciprocal) coupling between emitters | C. W. Gardiner, Phys. Rev. Lett. 70, 2269 (1993), and H. J. Carmichael, Phys. Rev. Lett. 70, 2273 (1993) (cascaded open systems); A. Metelmann and A. A. Clerk, *Nonreciprocal photon transmission and amplification via reservoir engineering*, Phys. Rev. X 5, 021025 (2015); K. Fang et al., Nat. Phys. 13, 465 (2017); P. Lodahl et al., *Chiral quantum optics*, Nature 541, 473 (2017); N. Hatano and D. R. Nelson, Phys. Rev. Lett. 77, 570 (1996); S. Yao and Z. Wang, Phys. Rev. Lett. 121, 086803 (2018) | one-way coupling through a unidirectional reservoir (cascaded systems); directionality from matching a coherent coupling with its dissipative counterpart through a shared reservoir, over a wide bandwidth, in isolators and amplifiers; its experimental realization; non-reciprocal many-body chains and the non-Hermitian skin effect | **known class.** Round 17's imposed rule is a cascaded coupling; round 18's absorbing stream is a spatially distributed directional loss. **Most useful lead of this audit:** Metelmann and Clerk's construction makes a coupling one-way by balancing its coherent part against a dissipative part through a reservoir, and stays passive. The companion's stream could be that reservoir; this is the next medium to build (round 19, §29.2) |

## 3. The astrophysical programme

| our claim | closest published work | status |
|---|---|---|
| Lensing follows the galaxies in colliding clusters without dark matter | the dark-matter interpretation of the Bullet Cluster (Clowe et al. 2006) and of the 72-collision stack (Harvey et al. 2015); MOND's difficulty with them is published | our mechanism (hot galaxies carry their companion, cold gas does not) **not found** elsewhere in this search |
| Ellipticals lens more than spirals because their stars move randomly | KiDS-1000 (Brouwer et al. 2021) measured the difference; dark-matter and MOND readings are published | a velocity-dispersion-weighted source term **not found** (modified-gravity models in the search, e.g. MOG, scale with mass and size, not with the source's velocity dispersion) |
| Solar-System release, wide binaries, Cassini | the MOND external-field literature and the Cassini constraints (Hees et al. 2014; Park, Hees, Famaey, Desmond, Durakovic 2026) | our release length is a fitted amendment, not a derived mechanism |

## 4. A novelty statement that could survive

> Building on established wave-mediated forces between in-phase emitters (Bjerknes), gain-induced pulling forces and
> directional (reservoir-engineered) coupling, we show that inverted, self-sustained emitters in a shared medium are
> pulled in proportion to the wave's height, that random motion raises their output while collisions suppress it, and
> that a directional medium lets warm sources pull harder, as the square root of their extra output, without the drag
> of a moving medium; the effect predicts that lensing at fixed visible mass grows with the stars' independently
> measured velocity dispersion.

The first two clauses need the directed search on injection-locked particles; the third needs the passive directional
medium actually built (round 18 reached 73–84% of the square root of the model's own glow with an absorbing stream,
about half of the law's heat gain; round 19 found no heat gain on average in the medium the data allow, see
STEP-BACK-AUDIT.md); the prediction needs a frozen law.

## 5. Searches still owed

1. Forces on injection-locked, lasing or otherwise self-sustained particles in a weak external field (the "height"
   scaling).
2. Velocity-dependent or temperature-dependent gravitating mass beyond Tolman–Whittaker: laboratory tests of the
   temperature dependence of weight, and astrophysical proposals with a velocity-dispersion-weighted source.
3. Streaming-medium and preferred-frame models of gravity (aether and "flowing space" pictures) and their Solar-System
   preferred-frame limits, since the companion's stream defines a local frame.
4. Radiation-mediated many-body synchronization with forces (optical matter, driven-dissipative crystals).
5. Follow the references and citations of every paper above.
