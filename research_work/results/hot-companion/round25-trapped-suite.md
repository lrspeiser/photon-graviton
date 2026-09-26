# Round 25: the trapped companion on the full suite (registered before it was run)

Registered 26 September 2026, before either candidate below was run on any test of the suite. The results go in
README section 36.9; this file is not edited after the run except to add a pointer to them and to list any
deviation from it.

## What is tested

The self-energy repair (the owner's, README §36.8) gives a hot-glow packet of free energy ε0 the energy ε = n ε0 in the
aligned case, with

    n = 1 - 2 (sqrt(Q + T) - sqrt(T)) / sqrt(a),        n < 0  exactly where  Q > a/4 + sqrt(a T),

Q the ordinary-source term (|g_N| in a round system) and T the hot brightness where the packet is. In a static field a
packet keeps its energy, so glow born where n < 0 can never reach a place where n > 0 (consequence A, trapping). The
candidate adds A to the adopted law (round 12) and nothing else. The fill-up of the trapped regions to n = 0
(consequence B) is not included: it moves SPARC's constants (a and g_d, §36.8) and needs its own refit.

## The rule, for every test

1. **Glow born where n < 0 is heard only where n < 0.** A place where n >= 0 hears only the glow born where n >= 0.
   A place where n < 0 hears all of it, as in the law (its S and g_hot unchanged).
2. **Systems the models resolve** (SPARC's bulges, the Milky Way's bulge and stellar halo, the six SLACS lenses): n is
   evaluated with the law's own fields at each place (Q = |g_N| of all the visible matter, T = the law's S there). Not
   iterated: a place just outside a zone, where the missing trapped glow would push n below 0, is counted outside;
   the size of that band is reported.
3. **Galaxies the models do not resolve** (the KiDS and Mistele lenses; the galaxies of the 12 X-COP clusters and of
   the five collision models) are given measured sizes: early types (KiDS red and bulge samples, Mistele's early
   types, cluster galaxies) a Hernquist profile with the circularised half-light radius, late types (blue, disc,
   Mistele's late types) an exponential disk with the exponential scale length, each the median of SDSS DR17 galaxies
   of the same stellar mass and split (`code/lens_sizes_sdss_v25.py` → `data/lens_sizes_sdss_v25.json`: the same
   selection as the heat measured in round 12; for example red at 10^10.55–10^10.65, R_e = 2.74 kpc; blue, R_d = 3.19
   kpc; red at 10^11.0–10^11.1, R_e = 4.61 kpc). Inside such a galaxy Q(r) = G M(<r)/r² (its baryons following its
   stars) and T(r) = its own glow (its heat k times the two-way shell sum of its stars) plus, in a cluster, the
   cluster's brightness there. The galaxy's trapped share is the fraction of its stars where n < 0; the rest of its
   glow leaves it.
   * Isolated lenses (KiDS, Mistele): the measured heat k (round 12) becomes k × (escaping share), at each sample's
     mass (KiDS: 10^10.6, the heat's 10^10.3–10.9 mean; Mistele: each mass bin), with the static law's mass and size
     factors. Their graded radii (44–300 kpc) lie outside every zone.
   * Cluster galaxies: one typical galaxy, M* = 10^11.05 Msun (the SDSS bin 10^11.0–11.1; the stellar-mass-weighted
     typical cluster galaxy is near M*), R_e = 4.61 kpc, taken in the surveys' units (no static conversion). The
     cluster's brightness at each place is the glow that escaped its galaxies, solved together with the escaping
     share (iterated from the law's brightness to the first self-consistent state). Places in the cluster (X-COP's
     radii, the collision maps' cells, the pre-collision clusters' Jeans profiles) lie outside the galaxies' zones and
     hear only the escaped glow.
4. **Which heat a cluster galaxy's own glow carries inside it.** The law gives every cluster star the cluster's heat,
   k = 3σ²/u² with σ ≈ 600–1,100 km/s, but inside one galaxy its stars move relative to one another only at the
   galaxy's internal speeds. The two readings are both run:
   * `trapped_own_heat_r12` (primary): the galaxy's own glow, inside it, carries the heat of its internal motions,
     measured by SDSS for red galaxies at 10^11.05 (σ_e 178 km/s; k = 3.16 at u = 169.4 km/s, scaling as 1/u²); the
     cluster's heat applies to the glow between galaxies (the law's cluster brightness).
   * `trapped_cluster_heat_r12`: the law's formula taken literally inside the galaxy, the own glow at the cluster's
     heat. T is then about k_cl ≈ 50–100 times the galaxy's own pull, n stays positive (n < 0 would need Q > about
     k_cl a), and the clusters and collisions should reduce to the law; only the galaxies change.
5. **The dwarfs** hear the Galaxy's glow born outside its zone (rule 2's Milky Way zone, per shell of the bulge and
   the stellar halo); their own n is positive (0.25–0.87, §36.8), so their own glow is unchanged. The precision tests
   carry no heat and change only through u (the Double Pulsar's emission a u / 2).
6. **Refit:** u on the 12 X-COP clusters (the suite's refit, `"refit": ["u"]`); a and g_d held (the repair leaves the
   static response unchanged).

Not included, and stated as limits: the cluster's own pull at a galaxy is left out of the galaxy's Q (it changes Q by
a few percent on average); one typical cluster galaxy stands for the whole luminosity function (the brightest cluster
galaxy would trap more of its glow in the cores); the Milky Way's spherical stand-ins for its hot components take the
escaping share of each shell from the real (flattened) densities on that shell.

## What was expected before running

§36.8 used the cluster's brightness alone for T inside a cluster galaxy (reading 4, primary, without the galaxy's own
glow) and found u ≈ 111–115 km/s, X-COP rms 0.203 and worst radius 0.15, and, at u = 113, KiDS all / red / blue
+0.075 / +0.068 / +0.055 dex with gaps 0.18–0.20. With the galaxy's own glow in T, less is trapped, so the refitted u
should come out higher than 113 in the primary reading and near 169 in the second.

## What would count

Against the round-12 baseline (full tier; 77 graded checks: 59 pass, 11 close, 7 fail; Cassini's Q2 is marked worse
for every law since its target changed in round 18 and is not counted):
* **Supported** (a replacement worth pursuing): at least 59 passes and at most 7 fails, and X-COP no worse.
* **Mixed**: within two passes of the baseline either way, with trades.
* **Not supported**: fewer than 57 passes or more than 9 fails.

Every change of grade is traced to the part of the rule that causes it. Implementation: `code/trapping_v25.py` (the
switch `companion_trapping` in `regression/law_config.py`, installed by `regression/common.apply_distances`; 'none'
must reproduce the baseline to the last digit).

## Results (added after the run)

README §36.9; runs in `run-trapped-v25/suite/`, table `run-trapped-v25/suite_table_v25.json`. Against the round-12
baseline (59 / 11 / 7): `trapped_own_heat_r12` (u refitted to 146.0 km/s) 57 / 12 / 8, inside the "mixed" band by its
numbers but with no grade gained (SPARC's bulge-dominated galaxies, the Milky Way at 15–27 kpc and Mistele's ellipticals
lose a grade); `trapped_cluster_heat_r12` (u 169.4 km/s, clusters and collisions as the law) 57 / 9 / 11, not supported.

Deviations from the text above: none in the rules. Implementation details not fixed above: each hook reads the switch
when it runs (so the suite's import order is unchanged); X-COP's brightness at its emitting shells is computed on 160
radii and interpolated; the lens sizes are interpolated in log M* between the SDSS bins' centres; the maps' iteration
stops when no cell's escaping share moves by more than 10⁻⁴ (four steps), the spherical ones at 10⁻⁶. One exploratory
run was added after looking (`candidates/trapped_own_heat_refit_au_r12.json`: a refitted too), reported as such.
