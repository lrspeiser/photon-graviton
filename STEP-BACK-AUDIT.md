# Stepping back: what did not go our way, and which doors are still open

*25 September 2026, after round 19. Written at the owner's request: "Make a list of the situations we would prefer to
have a different outcome and what we tried so far that was not ideal and why. I want to make sure we didn't close the
door or get into a corner with our derivation work."*

**How this was done.** Eight separate read-throughs covered the whole record:
* the hot-companion log (`research_work/results/hot-companion/README.md` §1–§29, rounds 1–19);
* `THEORY.md`, `PAPER-ROADMAP.md`, `NOVELTY.md`, `RULES.md` and `BLOG.md`;
* the earlier photon–companion era, 8–22 September (`CURRENT-STATUS.md`, `CHANGELOG.md`,
  `research_work/results/RESEARCH-CHECKPOINT.md`).

Each read-through listed three things: every result we would have wanted to come out differently, every idea we set
aside, and whether each "ruled out" holds in general or only inside one particular model. The numbers quoted here
were re-checked against the result files. **This audit changes nothing in the law, its constants or the locked
forecasts.** It corrects a few statements in the write-ups (§5).

---

## The short answer

1. **The law is not cornered.** It still meets the tests it was built on: galaxies, clusters, colliding clusters,
   galaxy lensing and the Solar System, with four shared constants plus the distance scale. Its standing misses (§1)
   each have several remedies that were never tried. Several of those misses trace back to our own conventions
   (distances, star masses, how cluster stars' speeds are computed), not to physics.
2. **The derivation work is in a corner, but we made it ourselves, and there is a marked way out.**
   * Since round 10, every small-scale model has made the pull the same way. Matter behaves like a laser medium: it
     keeps a fixed rhythm with the companion wave, a quarter of a beat ahead, and is pulled toward the source of the
     wave it amplifies.
   * Every setback with warm matter since then comes from needing that shared rhythm. That includes the diluted pull
     (round 10), warm sources falling out of tune (rounds 16–17), the need for a one-way wave (rounds 17–18), and
     round 19's "no extra pull from heat".
   * Round 13 wrote down the way out ("a receiver that responds to the wave's energy rather than its phase";
     "a pull drawn from the companion's energy flow rather than its phase"). It was never taken.
3. **Almost every "ruled out" from the models is narrower than it sounds.** The force experiments never had:
   * moving matter: "heat" was a frozen random mixing inside each piece, with no Doppler shifts (round 15's
     glow-only calculation did move its pieces, but no run has evolved hot motion, the forces, the wave and the recoil
     together; narrowed after review, §8);
   * cluster-level heat: heat weights of at most 16, against about 25–150 for galaxies in clusters;
   * large sources: never more than about 100 pieces, or more than 4.5 wavelengths across.
   
   Only a handful of verdicts are general (§3.1), and none of them blocks the law.
4. **Three decisions need the owner's call** (§4.3–4.5):
   * the field equation's form, which is borrowed from QUMOND, a MOND formulation that `RULES.md` names;
   * a dwarf-galaxy fix shelved only for lack of a derivation;
   * the distance law, and the no-expanding-universe rule that is not yet written into `RULES.md`.
5. **The audit found a few slips in our own write-ups** (§5). They are corrected alongside this document, and none
   changes the law or a forecast.

---

## 1. Results we would prefer were different: the data

The current law's standing (regression suite, full tier): **59 pass, 11 close, 7 fail.** The seven failures are five
faint dwarf galaxies, the Bullet Cluster's smaller half, and galaxy lensing by red galaxies.

### 1.1 Five faint dwarf galaxies move about three times too slowly

| | Draco | Ursa Minor | Sextans | Crater II | Antlia 2 | (Carina, close) |
|---|---|---|---|---|---|---|
| measured spread of star speeds (km/s) | 9.1 ± 1.2 | 9.5 ± 1.2 | 7.9 ± 1.3 | 2.7 ± 0.3 | 5.7 ± 1.1 | 6.6 ± 1.2 |
| the law | 2.8 | 3.4 | 2.1 | 1.0 | 1.1 | 3.5 |

* **Tried:**
  * Heavier stars, a mass-to-light ratio up to 3: Draco reaches 3.5 and Ursa Minor 4.2, not enough.
  * Weakening how much the Milky Way's pull holds back each dwarf's companion: over the ten dwarfs the misfit
    (χ²) falls from 135 to 60 as the hold goes to zero. MOND with the same stars gets 119.
  * A candidate with no hold (`no_hold`) scored 64 pass / 9 close / 3 fail in round 9 and broke nothing. It was
    shelved because "there is no derivation yet", and because it conflicts with the law's single square root for
    all pulls.
* **Why not enough:** Draco and Ursa Minor stay at 4.1–4.2 km/s whatever the hold.
* **Resting on:**
  * round, settled dwarfs with no tides;
  * the Galaxy's pull applied uniformly and at full strength;
  * our own speed estimator, which gives about 1.7× less than the one used in the MOND literature, so part of the
    miss may be method.
* **Door: open.** Remedies never tried: tides and dwarfs out of equilibrium; a no-hold rule with a derivation
  (history-dependent screening, §2.7, is one candidate source); the estimator.

### 1.2 The Bullet Cluster's smaller half has about half its measured lensing mass

* **Wanted** 2.47–2.85 × 10¹⁴ suns inside 301 kpc; **got** 1.38 × 10¹⁴.
* **Tried:**
  * the companion's memory of its source (round 4);
  * heat from streaming: negligible once memory is included;
  * tidal shaking: this is ordered motion, not heat;
  * a larger lost population (1:3): disfavoured at 2.5–3σ by the galaxies' speeds, and by the light (0.46 observed
    against 1.55 needed);
  * moving to our own distances: 0.94 → 1.29 × 10¹⁴;
  * heat gained in the crossing: +12–18%, about a fifth of the gap, tried on the Bullet only.
* **Resting on:**
  * Clowe's star masses at a mass-to-light ratio of 2, which that paper gives as an upper limit and allows from 0.5
    to 3;
  * a pre-collision mass ratio of 1:8 taken from dark-matter merger simulations;
  * target masses from lens models computed in standard distances, then converted.
* **Door: open.** Remedies never tried: the crossing heat with its full history, diffuse stripped stars, the star
  masses, a mass ratio not borrowed from dark-matter simulations.

### 1.3 Galaxy lensing (KiDS) sits 16% above the law, and red galaxies fail

* **Got:** all lenses +0.065 dex (close); red +0.078 (fail); disks +0.066 (close).
* **Tried:**
  * hot gas haloes (round 2; undone when gas stopped counting as hot);
  * the lenses' heat measured from SDSS spectra;
  * star-mass gradients (at most 0.09 dex);
  * a heat exponent of 1.75: this fixes KiDS but breaks the ellipticals' lensing speeds (rms z 2.6 → 5.2);
  * a slower companion: this fixes the red lenses but breaks the same test.
* **What the audit found:** in the survey's own distances the level is **+0.013 (red +0.029, a pass).** Most of the
  offset comes from converting into our static distance law: at the lenses' depth our brightness distance makes
  the star masses × 0.87 (§4.5).
* **The proposed fix:** gas around the lenses of about the stars' mass. It has not been measured.
* **Door: open.** It is tied to the distance law and to unmeasured gas.

### 1.4 Strong lenses need unusually heavy stars

* Six SLACS lenses need star masses 0.43 dex above the standard (Chabrier) value: 1.4–1.9 × Salpeter in our
  distances. Spectra suggest about Salpeter.
* In standard distances the need is 1.05–1.35 × Salpeter. The distance convention alone moves it by 0.13–0.16 dex.
* The studies we compare against derived their star masses with dark-matter haloes.
* **Door: open.** Remedies: star masses from spectra; the distance law.

### 1.5 The Milky Way runs 3–6% slow

* With McMillan's (2017) matter, which was fitted together with a dark halo, the Sun's orbital speed comes out
  209 km/s against 229–234.
* Round 19 refitted the Galaxy with its matter held to independent measurements:
  * the shape follows the 2026 Cepheid curve;
  * the level is 3–6% low (a common scale of 0.935–0.966 against 1 ± 0.02);
  * the inner points want a heavy bulge (2.1–2.3 × 10¹⁰ suns);
  * the vertical pull comes out 77–78 against 70 ± 5.
* **Tried:**
  * a later switch-off (g_d × 1.5): the Sun reaches 217.6, but the strong-lens and vertical-pull checks break;
  * more disk: this exceeds the star counts.
* **Reading:** the level matches the law's own offset in the 149 SPARC galaxies at this pull (+0.026 dex). It is a
  property of the law's switch-over region, not of the Milky Way.
* **Door: open.** The switch-off's shape is free: the galaxies cannot tell the shapes tried apart (differences of
  0.9 km/s or less).

### 1.6 Clusters: a 25% typical miss, and a steady trend with radius

* The law misses the 12 X-COP clusters' masses by 25%. Dark-matter haloes fitted to each cluster miss by 11%.
* The miss follows a steady trend with radius: the law predicts about 10% too little near the centre and about 20%
  too much at R500, in every refit so far.
* **A soft spot we made:** the stars' speeds in clusters are computed from the X-ray-measured gravity, which is the
  very thing being predicted. Using the law's own gravity instead raises the miss from 25% to 39% (round 3).
* **Star masses may be too low:** the X-COP star paper's own error budget allows up to 60% more stars, plus up to
  50% from light between the galaxies. More stars would mean a faster companion u, and u (169 km/s) is set by
  these clusters alone. It varies from 148 to 177 km/s between halves of the sample.
* **Door: open.** Remedies: self-consistent stars; the star-mass uncertainty; weak-lensing masses of clusters (planned
  on 22 September and never done).

### 1.7 Galaxy speeds in three colliding clusters come out low

| | MACS J0025 | Abell 520 | El Gordo (NW) |
|---|---|---|---|
| galaxy speeds, law vs measured | 669 vs 835 ± 59 km/s | rms z 2.7 | 1,014 vs 1,290 ± 134 km/s |

Abell 520's P6 clump also lenses 2.5 vs 4.3 × 10¹³. All three rest on star-mass conventions (a mass-to-light ratio
of 2; an inferred basis for MACS J0025), each of which can move the stars by up to 2×. Infrared star masses were
planned in rounds 8–14 and never done.

### 1.8 The locked wide-binary forecast now sits above the law's own prediction

* **The lock:** made in round 10 with round 9's constants. It says 4% more pull than Newton at 7,000 AU and 9% at
  20,000 AU.
* **The drift:** refits since then have lowered the law's own prediction: 1.076 at 20,000 AU in round 12, and
  1.057–1.068 with the refitted Galaxy in round 19. That is below the lock's "supported" window (1.08–1.10) but
  inside its "not refuted" range (1.02–1.2).
* **The consequence:** if Gaia finds about 6%, the old lock will read "not supported" while the current law fits.
  That is the price of locking. The amendments are labelled as such and must always be quoted beside the lock.

### 1.9 Galaxies: good, not better than good

* 149 galaxies: a typical miss of 15.9 km/s. MOND gets 16.1; dark matter gets 7.5 with 298 adjustable numbers.
* A trend of the outer residuals with surface density (+0.055 ± 0.025) was noted in round 16 and not followed up.

*Cassini is no longer a problem: within 0.2–0.4σ with the refitted Galaxy (round 19).*

---

## 2. Results we would prefer were different: the theory and the models

### 2.1 The pull needs laser-like matter (round 10, confirmed in round 15)

* **The rule for the mechanism our models use:** in a travelling wave far from its source, a body that takes its
  pull from the wave's recoil is pulled toward the source only if it feeds the wave. *(Narrowed after review, §8: this
  is not a general law of physics. Conservative forces between bodies that share a field, such as Casimir forces, and
  pulling forces from redirected scattered light work differently.)*
  * Ordinary absorbing matter is pushed.
  * An amplifier below its threshold is pulled only as the wave's intensity, which gives a Newton-like 1/r² pull.
  * Only an inverted, self-sustained emitter locked a quarter of a beat ahead is pulled in proportion to the wave's
    amplitude, as the law needs.
* **Not ideal:**
  * that kind of matter is an assumption, "the weakest link" (round 9);
  * a locked pull needs a shared rhythm, and that is the root of §2.2–2.5.

### 2.2 Many independent sources dilute the pull (round 10)

* **The problem:** sources with independent rhythms give a weaker pull when they surround the body: to 0.28–0.65 of
  the law's value inside clusters.
* **The data reject that dilution:**
  * galaxies: 19.6 against 15.9 km/s;
  * clusters: 0.38 against 0.23 in log units, i.e. 46% against 26% typical miss.
* **What the data select** behaves like "a single coherent flow": the strength comes from all the companion present,
  and the direction from its net flow.
* **Still missing:** no model of ours has produced that rule for warm matter. Simple incoherent schemes (independent
  phases, or a saturating amplifier) give the diluted form too.

### 2.3 Scrambled, warm sources weaken the chorus (round 13)

* **The first negative:** with a fixed rhythm offset, free random motion halves the pull (0.36–0.49 of rest at the
  higher speeds). The text reads: "excluded: the heat term's mechanism as the blog stated it … for any one fixed
  offset and for free amplitudes."
* **Then a nonlinear rule, an offset that depends on the wave's strength (rounds 13–14), gave the right pattern:**
  * motion raises the pull ×1.3–2.0;
  * collisions suppress the rise;
  * rotation gives none;
  * the pull follows √(cold + released power).
  
  Its problems: the released power grows more slowly than σ², and bodies lose step beyond a certain speed.
* **It was dropped, not refuted.** Round 15 switched to the store model, and the rule was never re-run with the
  later tools that book every watt exactly.

### 2.4 Warm sources fall out of tune (rounds 16–17)

* **The model:** "one kind of matter", where every piece both sends and receives and energy is booked exactly.
* **What happens:**
  * a cold source pulls;
  * warm sources lose their shared beat;
  * receivers lose step (0.55 → 0.16–0.27);
  * at heat weight 8 the net force is a push.
* **Internal designs:** twelve were tried and none protected the beat. Reciprocity suggests why, but only for designs
  that are linear, passive and unchanging and that send and receive through one channel. The two best designs were
  never run in the full model.
* **The glow:** it grows only about half as fast as the law's 1 + k (×4.7–5.0 at k = 8, against ×9).

### 2.5 A one-way wave fixes it, but only when imposed by hand (rounds 17–19)

* **By hand:** with coupling allowed only outward, warm sources keep one beat and pull 1.98× the cold value at heat
  weight 8. But the imposed rule is not a passive medium: it could create energy.
* **Physical versions (round 18):**
  * Waves carried by the stream are one-way, but they push every emitter with a drag as large as the law's own
    pull. Ruled out.
  * A stream that absorbs counter-moving waves has no push. It gave "73–84% of the square root the law needs". That
    was **measured against the model's own glow**; against the law's heat gain the raw pull ratio at k = 8 is about
    1.6 against 3, about half (correction, §5).
* **Round 19:**
  * **The exact medium:** it is passive only if its absorbers move at the wave's own speed. The energy bill needs the
    companion's stream to move at least twice as fast, so the absorbers cannot simply ride the stream (correction,
    §5).
  * **The ray description** overstated how one-way the medium is when absorption is strong.
  * **What the clusters allow:** only weak absorption. A wave must survive about 300 kpc inward. X-COP alone would
    accept 33 kpc; what breaks the galaxy-lensing tests is the slower companion the clusters then ask for.
  * **The heat gain:** in that medium, averaged over three arrangements, there was none: 0.90–0.97 of the cold pull
    at k = 2, where the law needs 1.73, and a push at k = 8.
  * **The scatter** between arrangements is very large (correction, §5). With the plain two-way wave, one arrangement
    reached 2.2× and 1.8× at k = 2 (the law's √3 = 1.73), and at k = 8 one reached 1.7× and 3.5× while the others
    pushed.
  * **Size:** the models' pull depends strongly on source size at fixed mass (2–15×), which real galaxies do not show.

### 2.6 No field theory yet

* As written, the heat term breaks Newton's third law.
* Round 2's repair, a QUMOND-type action, pushed hot cluster stars outward by 32–73% of gravity; data exclude that
  (round 3).
* It was replaced by a postulate: the companion carries momentum. That leaves a small self-force on lopsided
  systems, 12% of the Bullet's mean pull.
* There is no action and no relativistic version. The field equation's form is borrowed from QUMOND, and how light
  responds is assumed to follow Einstein's theory.

### 2.7 The switch-off near stars is assumed

* The form exp(−|g_N|/g_d) was a choice. The galaxies cannot tell it from other shapes: 15.86–16.08 against
  15.87 km/s.
* No matter model has produced a switch-off; the laser-like receivers only saturate at a ceiling.
* The release length L is bounded by Cassini, not derived.
* The "blockers" derivation supplied to the project gives both from one lifetime and predicts a step in wide binaries
  (round 19). It is not yet found in the matter model.

### 2.8 Inertia, and why everything falls alike

Never started. The models define no inertia. Round 17 found that equal falling would need every kind of matter to
share the same radiator share, and to have an inertia in proportion to its quiet store.

### 2.9 Energy

* If the heat term is paid for by continuous extra emission, that power cannot come from the random motion itself:
  it would drain the motion in 85 million years. It must then come from a store inside matter, and that store is
  postulated. *(Narrowed after review, §8: this bill belongs to continuous-emission mechanisms; a static,
  conservative attraction does not need a running power supply to hold a settled arrangement together.)*
* The Sun would lose 2.3 parts in 10¹⁵ of its mass a year, about 30 times below current precision.

### 2.10 Gas not counting as hot

* This was derived in a toy model: collisions scramble the rhythm (Dicke narrowing). It rests on an unestablished
  assumption about cluster plasma: ions re-scattered every thousand years or so, against a mean free path of about
  10 kpc.
* The clusters themselves barely care: 0.217–0.228 whatever weight the gas gets. The choice rests on the colliding
  clusters.

### 2.11 One speed doing three jobs

* The companion speed u sets three things: the heat weight k = 3σ²/u²; the constant a = 2ℓ/u; and how fast the
  companion travels, which governs its memory.
* Round 10 noted these are three distinct speeds in the model. Round 16 made them one by construction.
* u has only ever been fitted to cluster X-ray masses, and those fits use stars' speeds taken from the same X-ray
  masses (§1.6). It has never been measured independently.

---

## 3. Doors: really closed, closed only inside our models, and never opened

### 3.1 Closed by data or by a general law of physics (these should stay closed)

* **Rivals and missing pieces, rejected by data:**
  * plain Newton with visible matter;
  * any version without heat, on clusters (misses of ×2.9 or worse);
  * MOND on clusters.
* **Near stars:**
  * an extra pull that is never switched off near stars: at Saturn it would be 6 × 10⁻⁸ m/s², which planetary
    tracking rules out;
  * an instant switch-on, which misses Cassini by about 8σ.
* **How contributions add up** (these tests are what make the law distinct):
  * cold matter adding up without cancelling: galaxies, 18.9 against 15.9 km/s;
  * heat cancelling like arrows: X-COP, 0.288 against 0.221;
  * a diluted pull from independent waves: galaxies and X-COP.
* **Collisions:** lensing sitting on the gas in colliding clusters.
* **General physics, within the mechanisms they apply to** (§8):
  * for continuous-emission mechanisms, heat powered by motion (the 85-million-year drain);
  * for travelling-wave recoil, a pull without feeding;
  * passive or below-threshold receivers giving the amplitude term (they respond to intensity, a Newton-like pull);
  * crests moving at light speed, ruled out by the planets and by the galaxies.
* **Media:**
  * hot matter heard strictly one way, with outer shells unheard: 11 suite failures even after refitting;
  * waves swept past emitters (drag);
  * point-like absorbers (drag).
* **Distant cluster stars:** older, heavier stars there. Their colours show younger stars.
* **From the earlier era:**
  * redshift from energy loss alone, which fails supernova time-stretching (DES: b = 1.003 ± 0.011, where energy
    loss predicts 0);
  * mass-energy funded by light, short by thousands of times;
  * linear and steady "memory" responses;
  * radiated 1/r mediators, which would spend a galaxy's rest energy in about 70,000 years;
  * companion populations that carry their own conserved amount (now also excluded by the no-dark-matter rule).

### 3.2 Closed only inside our simplified models (all still open)

| verdict | what it rests on |
|---|---|
| "Warm matter cannot pull harder" (rounds 13, 16, 17, 19) | In the force experiments, heat as a frozen random mixing: no moving pieces, no Doppler shifts (round 15 moved pieces only to measure their glow). Heat weights of at most 16; cluster galaxies have about 25–150. 24–100 pieces. Sources 1–4.5 wavelengths across. One scalar wave carrying both cold and hot glow. Senders and receivers identical. A pull that works only through a locked rhythm. Three arrangements, with very large scatter between them (§2.5) |
| "No internal structure keeps warm matter in tune" (round 17) | Designs that are linear, passive, unchanging and single-channel, in a reduced model; the two best designs never run in full |
| "Free amplitudes do not escape" (round 13) | One seed, at rest, one sparse cloud |
| The strength-dependent rhythm rule (rounds 13–14) | Dropped, not refuted |
| "Strong one-way absorption is excluded" (round 19) | One constant absorption rate. One radial stream per cluster. Hot and cold glow in one channel. The companion speed fitted to X-COP only. Never tested: absorption strong inside galaxies but weak between them; the colliding clusters' maps |
| "Gas does not carry heat" | A toy model plus an unestablished plasma assumption |
| "The companion cannot be a medium" (round 4) | Shown only for a fixed medium, not for one carried along with matter |
| "An action with reaction forces is excluded" (round 3) | Only that one action, with stars hot |
| The dwarfs' no-hold rule (round 9) | Shelved for lack of a derivation, not excluded by data |
| "Separate systems cannot pull separately" | Holds only given the heat term's single square root |
| Our distance law's shape; the other ("metric") geometry | Judged only under the current constants and star masses. The supernovae prefer an extra term (Δχ² ≈ 101) that was never derived or adopted |
| Heat exponent 2 rather than 1.75 | A trade-off between two lensing tests; it sits behind a switch |
| The switch-off's shape | The galaxies cannot tell the shapes apart |
| "Sources smaller than the wavelength as a way out" | Never actually reached; the smallest source tried was one wavelength |
| From the earlier era: the "gravity cliff" and the "depth law" | Dropped by a change of direction, not by a negative result. The depth law met X-COP's needs to within 8% (median) without being fitted to clusters, using a depth variable close to our heat weight. Many earlier campaigns also stopped at numerical-accuracy gates, not at physics |

### 3.3 Never tried

1. **A pull set by the companion's energy and flow rather than its rhythm.** This was round 13's own recommendation.
   * The data select a specific structure: the strength comes from all the companion present, the direction comes
     from its net flow, and warm contributions add without cancelling.
   * A body that senses the companion's total energy as one number, and its net flow as a direction (roughly, the
     way a sailor feels both how strong the wind is and where it comes from), would give exactly that structure. It
     would need no shared beat at all.
   * Neither our locked-rhythm models nor simple incoherent receivers can give that structure (§2.2), so it needs a
     genuinely new receiver or medium. **This is the most promising unopened door.**
   * *Refined after review (§8):* receivers that respond to intensity were tested (round 15's passive and
     below-threshold amplifiers give a Newton-like, intensity pull). What was never built is a nonlinear response
     that turns the companion's total energy and flow into the law. And "energy instead of rhythm" is not enough on
     its own: a body whose energy depends only on the local strength of the companion cannot give both the square
     root of mass and the inverse distance. The target is a force from the energy and stress of the coupled
     matter–medium state.
2. **Real motion in the force experiments:** moving pieces with Doppler shifts at cluster-level heat, with the
   forces, the wave and the recoil evolved together (round 15 moved pieces only to measure their glow).
3. **Large, dilute sources:** many wavelengths across, with hundreds or thousands of pieces, measured for the pull.
   So far only their glow has been measured.
4. **Separate channels for the hot and the cold glow.** For example, the cold glow guided along gravity's lines and
   the hot glow spreading in all directions, as the law's two terms already suggest.
5. **A heat term produced at the receiver**, rather than heard from the sources.
6. **A vector or tensor companion.** Every model so far used a scalar wave.
7. **Receivers that differ from sources**, and internal parts that are nonlinear, change in time or are
   non-reciprocal.
8. **The strength-dependent rhythm rule of rounds 13–14, redone with the exact energy booking.**
9. **A field theory for the whole law, with an action**, and light bending derived from the same coupling.
10. **Self-consistent cluster stars:** their speeds taken from the law's own gravity.
11. **Tides for the dwarfs**, and the blocker history applied to the dwarfs.
12. **Tests of the heat term dropped along the way:**
    * weak-lensing masses of clusters;
    * relaxed against disturbed clusters;
    * galaxy groups;
    * Abell 1689;
    * measuring u independently of the clusters.
13. **The distance law's shape, and a mechanism for redshift and supernova time-stretching in a universe that does not
    expand.** This has been parked since the earlier era: the static law assumes both without deriving them.

---

## 4. Corners we may have put ourselves in, and the way out of each

1. **The rhythm-locked pull, since round 10.**
   * Way out: model a pull that responds to the companion's energy and flow (§3.3, item 1).
   * Keep the rhythm models as one possible microscopic picture, not the only one.
2. **The paper's working scope, "motion-enhanced attraction in an active streaming medium".**
   * The models no longer support the "motion-enhanced" part.
   * Until the heat term has a mechanism, the bounded claim is the cold law and the collision rule (as the blog's §9
     now says). The data's case for the heat term stands on its own, and the frozen lensing test (blog §8,
     prediction 13) checks it directly.
3. **The field equation borrowed from QUMOND (owner's call).**
   * The conflict: `RULES.md` §1 lists "AQUAL or QUMOND field equations" among what we do not adopt as our law. The
     law uses that architecture, credited as borrowed, and without the action that justified it (§2.6).
   * The cost: its structure for an outside pull is where the law's most MOND-like problems come from (the dwarfs,
     and originally Cassini).
   * The choice: accept it as borrowed machinery, or make deriving our own field equation a priority. Doing that would
     also supply the missing action and light's response.
4. **The single-square-root rule vetoing the dwarf fix (owner's call).**
   * `no_hold` was the best-scoring candidate (dwarfs χ² 60 against 135) and was shelved for lack of a derivation.
   * The release length L was adopted without one.
   * One standard should apply to both. The blocker picture might supply the derivation.
5. **Our distance law (owner's call on where the rule is written).**
   * What stays: the no-expanding-universe rule is the owner's.
   * What is only a choice: the specific formula (1 + z = e^{αD}, fixed geometry, (1 + z)² dimming). It is one
     choice among several and was never derived, and the supernovae disfavour its shape (Δχ² ≈ 101).
   * What it costs: it creates most of the KiDS offset and part of the strong lenses' heavy-star need, and every
     constant moves with its scale (a from 5.0 to 7.2 × 10⁻¹¹ m/s² over its tested range).
   * Way out: treat its shape as open, and re-grade the lensing tests under the alternatives.
   * The rule itself is not written in `RULES.md`. It comes from the owner's round-9 instruction, and the blog lists it
     as a fourth rule.
6. **One u for three jobs, fitted only to cluster X-ray masses** (§2.11, §1.6).
   * Way out: self-consistent cluster stars.
   * Measure u separately: from travel (collision timing, memory) and from galaxy lensing.
7. **The locked wide-binary forecast** (§1.8).
   * Keep it locked, as intended.
   * Always quote the current law's amended prediction beside it.

---

## 5. Corrections made alongside this audit

None changes the law, a constant, a forecast or a figure.

1. **The exact absorbing medium (round 19).**
   * What was shown: it is passive for absorbers that move at the wave's own speed (the calculation's units set
     u = c).
   * The problem: the energy bill needs the companion's stream to move at least twice as fast as the wave's crests.
     Absorbers carried by that stream would amplify forward waves within about 60° of it, the "pumping" that §29.2
     used to rule out point absorbers.
   * The correction: "a physically consistent absorbing stream exists" holds only if the absorbers are a separate
     component moving at the wave's speed, not the stream itself. Recorded in the log (§29.2 and §30), the theory
     notes, the blog and the page.
2. **Round 18's "73–84% of the square root the law needs".** It was measured against the model's own glow. Against
   the law's heat gain, the raw pull ratio at k = 8 is about 1.6 against 3, about half. Corrected in the blog and the
   page, and noted in the log.
3. **The frozen prediction's fingerprint.** 01219170… is the SHA-256 of the prediction text as the script wrote it,
   without the file's final newline; `sha256sum` of the file gives c0b10567…. The content is unchanged. Noted in the
   log (§29.7).
4. **Round 19's "no heat gain".**
   * The quoted errors treat two distances × three arrangements as six independent values, which they are not.
   * With the two-way wave, the arrangements scatter from a push to more than the law's gain (§2.5).
   * The honest reading is "not reproduced on average, with very large scatter; three arrangements cannot settle it".
     Noted in the log.
5. **Unit labels (rounds 17–18).** The wavelength is 1 in those models, not 2π.
   * Round 17's compact source has a radius of 0.6 wavelengths.
   * Round 18's balls have mean spacings of about 1.3 wavelengths (radius 3) and 13 wavelengths (radius 30).
   * The ray-form passivity note compares against a wavelength of 1.
   * Corrected in the log; no conclusion changes.
6. **Mistele et al.'s four ratios (round 11)** are bins of stellar mass, each averaged over 50–300 kpc. They are not
   radii. Noted in the log.
7. **Round 17's one-way runs in the full model** booked energy but did not check momentum balance. Noted in the log.

---

## 6. The earlier era (8–22 September), briefly

* **Closed by data or physics:**
  * redshift from energy loss alone, which fails supernova time-stretching;
  * a time field that slows clocks without expansion: protecting atomic clocks forces rulers to shrink, and
    laboratory cavities rule out the version with fixed atoms by 9–47×;
  * mass-energy funded by light, short by a median of about 5,400× even converting all starlight for 10 billion
    years;
  * companion populations that need 10²–10⁵ times the cosmic mean density or form outer envelopes;
  * linear and steady-memory responses;
  * radiated 1/r mediators.
* **Closed only by a change of direction** (worth a second look):
  * the "gravity cliff" factor;
  * the empirical "depth law";
  * the hint that the galaxy constant is about cα/10 (α being the redshift rate);
  * the lagging response of rotating disks;
  * the lensing–dynamics cross-prediction protocol (JR-9);
  * the Milky Way test stars, never opened.
* **Many campaigns stopped at numerical-accuracy gates** (for example 9.2% against a 5% gate), not at physics.
* **Every law that fitted cold disks collapsed to the deep-MOND form.** The baryonic Tully–Fisher relation forces
  that, and our law's cold limit meets the same relation.

---

## 7. What to do next, and why this is good news

1. **Open the unopened door:** a receiver or medium whose pull follows the companion's energy and flow, not its
   rhythm. If it works, the heat term, the direction rule and the absence of dilution come together naturally.
2. **Put real motion into the matter models**, and push to cluster-level heat and larger sources.
3. **Make the cluster stars' speeds self-consistent**, and weigh clusters by weak lensing.
4. **Reopen the distance law's shape** and its redshift mechanism, and re-grade the lensing tests under it.
5. **The owner's decisions** (§4.3–4.5).
6. **The frozen lensing test** checks the heat term in data directly, whatever its mechanism turns out to be.

**Why this is good news.** Nothing found here rules out the law.

* The recurring obstacle in the models comes from one modelling choice made in round 10. The record itself already
  points to the alternative.
* Several of the data misses trace to conventions we can revisit cheaply.
* The doors that are genuinely closed are the ones that should stay closed: Newton, MOND and dark matter as answers;
  an instant switch-on near stars; dilution. They are also what makes the law distinctive.

---

## 8. Refinements after review (25 September 2026)

A review of this audit (supplied by the owner, with a new calculation; `research_work/results/hot-companion/
energy-shift-v20/`, taken up in the log's §31) narrowed four statements. The corrections are made above.

1. **"A wave pulls a body only if the body feeds it"** holds for the travelling-wave recoil mechanism our models use,
   not for every force. Conservative forces between bodies sharing a field (Casimir forces) and pulling forces from
   redirected scattered momentum (optical pulling, Chen et al. 2011) are different.
2. **The 85-million-year drain** applies when the heat term is continuous extra emission paid for by motion. It is not
   a cost of every static attraction.
3. **"Our simulations never had moving matter"** was too broad. Round 15's `wave_dark_v15.py` moves its pieces
   (free, colliding, rotating, boosted) and recomputes their couplings, but it measures only the released glow. The
   missing experiment is hot motion with the forces, the wave and the recoil all evolved together.
4. **"An energy-sensitive receiver was never tried"** was too broad. Round 15's `receivers_v15.py` has passive and
   below-threshold receivers whose forces follow the wave's intensity. What is missing is the specific nonlinear,
   many-direction response that would turn the companion's energy and flow into the law.

The review also sharpened the way out: not "energy instead of rhythm" but "a force from the energy and stress of a
coupled matter–medium state". Its new calculation shows attraction through correlations between pieces with no
rhythm on any single piece (with up to eight pieces and any number of them excited), and it proves two limits. A body
whose energy depends only on the local companion strength cannot give both of the law's exponents, and a fixed blend
of short-range forces cannot fall more slowly than 1/r².
