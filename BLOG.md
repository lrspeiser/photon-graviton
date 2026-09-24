# Gravity that streams

*Hot-companion gravity: one law for rotating galaxies, bending light, galaxy clusters and
colliding clusters, with no dark matter and no expanding universe.*

**Rewritten from scratch on 23 September 2026 (rev 12); updated that day and the next (revs 13–24).**
* Rev 13 added the companion's memory (§3.10).
* **Rev 14** adds §4, the law piece by piece: where each part may come from, and why it works
  so widely. It also brings the Bullet Cluster's galaxy speeds and strong-lensing masses into
  line (§6.3), and corrects how rev 13 read one lensing measurement.
* **Rev 15** tests rev 14 against new data:
  * a star count from the Legacy Survey confirms the main cluster's outskirts;
  * the stars the "bigger smaller half" needs don't show up, so its lensing strength is open
    again (§6.3);
  * the wide-binary prediction, 19% extra pull, is new (§6.7).
* **Rev 16** answers a direct question: have we checked the law against the Milky Way, and
  against the stars and lenses the field treats as the standard tests? Now we have, side by
  side with MOND and dark matter, with every measured value traced to its paper:
  * the Milky Way (§6.8);
  * the standard lensing examples (§6.9);
  * the precision tests of gravity (§6.10);
  * one coverage table of everything (§7).

  The results:
  * **Agree:**
    * the Milky Way from 15 to 27 kpc;
    * its mass out to 200 kpc and its escape speed;
    * the pull above the disk;
    * S2, pulsars and the planets;
    * the absolute lensing of both spirals and ellipticals, where our law beats MOND.
  * **Fall short:**
    * the Sun's own orbital speed (8% slow);
    * six faint dwarf galaxies;
    * Cassini's measurement of the Galaxy's field inside the Solar System.

  For each shortfall we found what would have to change, and each fix is testable.
* **Rev 17** does two things you asked for:
  * **A regression suite (§6.12)**: every test the law has faced, 89 checks (76 of them
    graded against published measurements), rerun in one command. It shows what any change fixes and what it
    breaks. On the unchanged law it reproduces every number published so far.
  * **Three more colliding clusters (§6.11)**, run with the Bullet's machinery and published
    inputs only:
    * **MACS J0025.4−1222** passes on lensing masses, peak positions and galaxy speeds.
    * **Abell 520**: its disputed "dark core" comes out of ordinary gas and the surrounding
      galaxies' heat, with a lensing mass between the two teams' measurements.
    * **El Gordo** passes if its stars weigh about twice the colour-based estimate, which is
      within that estimate's stated uncertainty.

  The suite's first use tested rev 16's three proposed fixes:
  * the Cassini fix costs nothing elsewhere;
  * the dwarf fix must depend on speed;
  * the Sun's speed is better fixed in the Milky Way's disk than in the law.
* **Rev 18** answers your follow-up: the Cassini fix was called free, but never adopted.
  * **Adopted (§6.10).** The companion now needs about 700 years of travel, 0.15 parsec, to
    break free of its source. On the full regression suite Cassini passes (4.6 against
    (3 ± 3) × 10⁻²⁷ s⁻²) and nothing else gets worse. The law now scores 62 pass, 8 close and
    6 fail (from 61, 8, 7). One prediction changes: wide binaries get 9% more pull at 20,000 AU,
    not 19%.
  * **What it opens (§6.13).** The fix also removes what blocked the dwarf-galaxy fix. With the
    Milky Way's hold on the dwarfs taken away, 6 of the 10 dwarfs agree instead of 4. That beats
    MOND with the same stars. We have not yet found the physical reason, so this second change is
    not adopted.
  * **A plan for everything still open (§9)**, item by item: what is wrong, what would have to
    change, and the next test.
  * **An audit for borrowed assumptions (§9, item 1).** Our tests of the far clusters had
    quietly used the expanding universe's distances, and published star masses cap star ages at
    the Big Bang's timeline. With the project's own distances, El Gordo is 58% larger and its gas
    68% heavier. Its gap has to be recomputed before we trust it.
  * **Ten proposals for finishing the theory (§9.1)**, each evaluated: where we stand, how we
    would work on it, and the first step.
* **Rev 19** takes the next steps toward first principles, and takes the far clusters out of the
  expanding universe's distances.
  * **Where the pull comes from (§4.14).** In a model where the only interaction is local, matter
    that feeds the passing companion a quarter beat ahead of it is pulled toward the companion's
    source. The pull comes out exactly in proportion to the companion's strength, as the law
    assumes. A 3D simulation of the wave confirms it to within 4%.
  * **How the companion adds up (§4.14).** We tested the three natural ways in which the
    companion of many pieces of matter could combine. Only the one our law has used since round 3
    fits both the 149 galaxies and the 12 clusters. The other two miss the galaxies by 19–20 km/s
    (ours: 15.9), and one also misses the clusters by 47% (ours: 25%). So three properties of the
    companion are now set by the data rather than by our choice.
  * **Our own distances (§6.14).** The far colliding clusters and the six strong lenses were
    recomputed with the project's own distance law, and the suite now grades them that way. The
    same light now means 20–29% less star mass, and five far-cluster grades slip: the suite reads
    58 pass, 11 close, 7 fail. Stars 1.4 times heavier would restore them all, and more (63 pass).
    Stars older than the Big Bang's timeline allows would weigh about that much, so that is the
    next thing to compute. *(Rev 20: neither is needed; see below.)*
  * **A first locked forecast (§6.7).** Our prediction for wide binary stars was recorded in the
    repository before Gaia's next data release.
* **Rev 20** takes rev 19's three next steps.
  * **Why orderly matter adds up like Newton (§4.15).** If the companion loses no energy on its
    way out, has no whirlpools, and matter feels it as one stream, its flow must follow Newton's
    field lines exactly. That is the rule our law has used since round 3, now derived instead of
    picked by testing. The alternative would carry more energy out of a galaxy than the galaxy puts
    in. The galaxies and clusters demand the rule to within a few percent.
  * **A mistake in our inputs, found and fixed (§6.15).** The 12 clusters that set the companion's
    speed list their stars inside a circle on the sky, which includes stars in front and behind.
    We had used those as the stars inside a sphere, so we had 1.3 to 2 times too many. Corrected,
    and with these clusters also put into our own distances, the companion's speed is **163 km/s,
    not 197**, and the clusters fit slightly better. Rev 20 adopts the refitted constants.
  * **Everything in our own distances (§6.15).** The Bullet Cluster, galaxy lensing and the
    calibrating clusters now use the project's distance law, as the far clusters have since rev 19.
    With the companion's speed measured the same way, the far clusters come back: El Gordo's lensing
    agrees with its published stars, and rev 19's "1.4 times the stars" is no longer needed. Rev
    19's idea that the far clusters' stars are older than the Big Bang allows is withdrawn: the
    observations show they are younger.
  * **Two sharper tests (§6.15).** In our own distances, lensing around spiral galaxies is about 20%
    stronger than our law predicts. That depends on one constant of our distance law, the rate at
    which light loses energy, and a 13% lower rate would fit. So galaxy lensing now measures our
    distance law. *(Rev 21: not so. With the galaxies refitted alongside, the rate barely moves it;
    §6.16.)* Ellipticals' extra lensing over spirals comes out 0.23 dex against 0.15 measured,
    which depends on the companion's speed. In MACS J0025 the north-western lensing peak slides onto
    the gas unless the collision is younger than about 350 million years.
  * The regression suite: 57 pass, 9 close, 10 fail (rev 19: 58, 11, 7).
* **Rev 21** works through rev 20's next steps one at a time, updating this page as each lands
  (§6.16).
  * **MACS J0025's age, from its own shock fronts.** Its two radio relics, which trace the shock waves
    the collision launched, sit right beside its two groups of galaxies. So the shocks have not yet
    outrun the galaxies: the collision is young, 0.1–0.36 billion years old in our distances, not the
    0.5 we had assumed. At that age its north-western lensing peak sits on its galaxies, as observed. The
    suite reads 58 pass, 9 close, 9 fail.
  * **The lensing galaxies' star speeds, measured.** The comparison with galaxy lensing had assumed
    one star speed for every elliptical-like lens and no heat at all for spiral-like ones. SDSS spectra
    of 119,000 galaxies of the same masses show the red lenses' stars are a little slower than assumed,
    and the blue ones carry real heat in their central bulges. With the measured speeds, the extra
    lensing of ellipticals over spirals comes out 0.13–0.16 dex against 0.15 measured: **it matches**.
    The suite reads 59 pass, 10 close, 8 fail.
  * **Our distance law's one number, fitted to everything at once.** Our law turns a redshift into a
    distance with one rate. Fitted jointly to 1,365 supernovae and to the 81 galaxies whose distances
    come from their redshifts, it comes out 5% lower than before (a Hubble-constant-like 70.9 instead
    of 74.6), and the companion's speed moves from 163 to 169 km/s. Galaxy lensing turns out not to
    depend on this rate, as rev 20 had thought. Its remaining 16% excess is about the lenses' mass:
    gas around them, of the amount the KiDS team itself estimates but leaves out of the visible mass,
    accounts for all of it. The suite reads 59 pass, 11 close, 7 fail.
  * **A working model of the companion.** Of five simple rules for how the companion's energy could
    move, only one, following ordinary gravity's field lines outward, keeps all the energy, travels as
    one stream and forms no whirlpools, and a short argument shows why it must. That turns rev 20's
    last assumption into a physical property: the companion is a wave guided by gravity's own field
    lines, as some waves in the Sun's hot gas are guided by its magnetic field.
* **Rev 22** does what you asked next: it joins our two working pieces, the pull (§4.14) and the scrambling
  (§3.5–3.6), in one experiment where only the motion changes (§6.17).
  * **With one fixed timing, heat weakens the pull.** Free random motion halves it, and collisions protect it. The
    reason is exact: a body is pulled by the recoil of what it feeds the passing wave, and bodies in step all feed or
    all absorb together, so a pull needs a loud chorus that motion can only break up. §3.5's "scrambled contributions
    don't cancel, so they pull harder" is ruled out as the mechanism, though the data still demand the heat term.
  * **Another solution works in part:** emitters that feed a quiet wave but absorb a loud one. A cold, dense ball then
    hushes itself yet still pulls. Free random motion strengthens the pull up to twice, frequent collisions hold the
    gain back, orderly rotation gives none, and the law's square root appears by itself. It works while the motion is
    slower than about a quarter of the re-timing speed; beyond that the test bodies lose step.
  * **The data on "doubling":** galaxy lensing pins how steeply the heat term must grow. Doubling the speed spread
    multiplies the extra pull by 1.8–2.0; the law has 2.0.
* **Rev 23** works through rev 22's next steps one at a time, updating this page as each lands (§6.18).
  * **A gentler heat term, tested on everything.** Growing as σ^1.75 instead of σ² erases galaxy lensing's 16%
    excess, but makes the lensing of massive ellipticals worse: 60 pass, 10 close, 7 fail against 59, 11, 7. It
    trades one tension for another, so the law keeps σ². Along the way we fixed a slip in the test suite's refit.
  * **Where "feed quiet waves, absorb loud ones" comes from.** A fixed supply of power plus a loss that grows with the
    wave's loudness sets each piece of matter's timing, like a generator or a motor on the grid. With it, the heat
    pattern comes back with no memory; and the same loss holds the companion back in a shape the 149 galaxies accept
    as well as our law's exponential release.
  * **How far the release goes.** The darker cold matter is, the lower the speed at which heat takes over, and test
    bodies keep step up to about the re-timing speed, so a deep enough hold covers the whole heat range. What doesn't
    match yet is the steepness: the release grows about 1.6 times per doubling of speed, not the 4 times our law needs.
* **Rev 24** takes up an independent calculation you sent, "motion-opened radiation from ordinary matter's internal
  oscillations", and joins it to our pull (§6.19).
  * **The idea:** each piece of matter holds a quiet internal vibration (its store of energy) and a few loud ones, and
    relative motion mixes the quiet into the loud, so the leak grows as σ². We reran it: every number matches.
  * **The wave alone gives half of it.** In a cloud of particles that talk only through the wave, spinning or moving
    as a whole releases nothing, and collisions hold the release back, as our law needs. But the release grows too
    gently (1.3–1.8 times per doubling of speed), for a clear reason: a cloud's quietness has no clean edge. The
    quiet store has to be inside each piece of matter. That also explains round 14's shortfall.
  * **Which bodies a wave pulls.** An ordinary self-sustained oscillator, the calculation's included, falls a quarter
    beat behind a passing wave and is pushed. An *inverted* one, holding energy it is ready to give, falls a quarter
    beat ahead by itself and is pulled in proportion to the wave's height. That is our round-10 rule, derived for the
    first time.
  * **The whole chain in one experiment:** the pull grows as √(cold + released glow), our law's heat term, with no
    square root put in. Doubling the speed spread doubles the extra pull where heat dominates (×2.48, 2.20, 1.97).
    Collisions hold it back; stopped, it vanishes.
  * **The energy bill:** being pulled at our law's strength costs a kilogram twice its cold output if the companion is
    slow, as ours is. A light-speed companion would cost 3,500 times more, which the planets' orbits already rule
    out.
* Every number below is computed from public data by a script in this repository, named
  where the number appears (§12).
* The earlier notebook (revisions 1–11), with all its retracted and retired claims left
  visible, is archived in
  [research_work/blog-archive/BLOG-notebook-rev11.md](research_work/blog-archive/BLOG-notebook-rev11.md).
* The technical record is in
  [research_work/results/hot-companion/README.md](research_work/results/hot-companion/README.md).

---

## In one minute

**The puzzle.** Galaxies and clusters of galaxies pull harder than the matter we can see
should allow:
* Stars at a galaxy's edge orbit too fast.
* Light passing a galaxy bends too much.
* Clusters hold their hot gas with 4 to 27 times the pull their visible matter provides.

When two clusters collide, the extra pull follows the galaxies, not the gas, even though the
gas outweighs the galaxies ten to fifty times. The usual answer is invisible "dark matter".

**Our answer, in plain words.**
1. **Every piece of ordinary matter feeds a faint, streaming companion to its gravity.**
2. **Orderly matter feeds it in step.** Contributions from opposite sides cancel, just as
   Newton's pulls do. Far from a galaxy, what survives falls off slowly enough to keep the
   outer stars orbiting fast. *Rev 20: this now follows from energy conservation (§4.15). Rev 21:
   and from a companion guided along gravity's field lines, which a working model confirms (§6.16).*
3. **Matter whose parts move randomly and freely, like stars in elliptical galaxies or
   galaxies in clusters, feeds it out of step.** Nothing cancels, so it pulls much harder.
   That is why clusters need so much extra pull.
   *Rev 19: points 2 and 3 are no longer just assumptions. Of three natural ways the companion
   could add up, only this one fits both galaxies and clusters (§4.14). Rev 22: the data still demand point 3, but
   "out of step" is not why. In a simulation, scrambling alone weakens the pull. What works is matter that holds its
   companion back when cold and releases it when moving freely (§6.17). Rev 23: that follows from a power balance,
   a fixed supply against a loss that grows with the companion's loudness (§6.18). Rev 24: with a quiet store inside
   each piece of matter that motion opens, the whole chain works in one experiment: the pull grows as the square root
   of cold plus released output, and doubling the speed spread doubles the extra pull (§6.19).*
4. **Gas does not count as "hot" in this sense.** Its particles collide so often that the
   companion sees them as sitting still. That is why, in a collision, the extra pull stays
   with the galaxies and leaves the gas behind.
5. **Strong gravity holds the companion back, and breaking free takes about 850 years of
   travel,** so the Solar System feels nothing. *Rev 23: the hold may be the same loss that sets matter's timing in
   loud waves; the 149 galaxies accept its shape as well as ours (§6.18).*
6. **The companion is slow, and it keeps moving the way its source moved.** After two clusters
   collide, each is still wrapped in the companion it had before, travelling on with its
   galaxies. That is why the smaller cluster in the Bullet Cluster still bends light like the
   whole cluster it used to be.

**The scoreboard.**

| Test | Ours | MOND | Dark matter |
|---|---|---|---|
| 149 galaxies' rotation speeds (typical miss) | **15.9 km/s** | 16.1 km/s | 7.5 km/s |
| 12 galaxy clusters' masses (typical miss) | **25%** | 2.9× | 11% |
| Bullet Cluster: lensing on the galaxies, not the gas | **yes** | no | yes |
| Bullet Cluster: main half's lensing mass, galaxy speeds and star count | **agree** | no | yes (fitted) |
| Bullet Cluster: smaller half's lensing mass | about half of what is measured (open) | no | yes (fitted) |
| 72 colliding clusters: lensing stays with the galaxies | **yes** | no | yes |
| MACS J0025.4−1222: lensing masses, peaks on the galaxies, galaxy speeds | **lensing and both peaks agree** at the age its shock fronts give (§6.16); speeds 2.9 error bars low | no | yes (fitted) |
| Abell 520's "dark core" (lensing without galaxies) | **3.5 × 10¹³ suns from gas and heat; measured 3.13–3.69** (our distances) | no | not expected |
| El Gordo: lensing mass and galaxy speeds | **lensing agrees with its published stars** (our distances, §6.15); NW galaxy speeds 22% low | no | yes (fitted) |
| Ellipticals bend light more than spirals (0.15 ± 0.04 dex, median over bins) | **yes: 0.13–0.16**, with the lenses' star speeds measured (§6.16) | no difference | yes, with tuned haloes |
| Solar System | **no extra pull** | small extra pull | no extra pull |
| Wide binary stars (data disputed) | **8% extra pull at 20,000 AU, 3% at 7,000 AU** (§6.7) | 43% | none |
| Milky Way rotation 15–27 kpc (4 Gaia studies) | **within 1–6%** | within 3% | 6–11% too fast |
| Milky Way rotation at the Sun (229–234 km/s) | 209 km/s, 9% slow | 223 | 234 |
| Milky Way mass inside 100 and 200 kpc | **agrees** | 40% high at 200 kpc | agrees |
| Galaxy lensing (KiDS): spirals and ellipticals | **agrees with both** in the standard distances; in ours the difference between them matches, and all lenses together sit 16% above our law, which gas around them as heavy as their stars would close (§6.16) | too high for spirals, too low for ellipticals | yes, with tuned haloes |
| Ten small galaxies around the Milky Way | 4 agree, 6 too slow; 6 agree without the Milky Way's hold (reason still sought, §6.13) | the same 4 | fitted one by one |
| Cassini: the Galaxy's field inside the Solar System | **passes** since rev 18 (§6.10) | about 10× too big | passes |
| Adjustable numbers | **4 in total**: 3 fitted to galaxies and clusters, 1 release length set by Cassini (the distance law's rate is set by supernovae) | 1 | 2 per galaxy or cluster (~320) |

**What is still open** (§9 has the plan, item by item):
* **six faint dwarf galaxies**, 1.5–5 times too slow. The Milky Way's pull dilutes their own
  companion. Without that dilution 6 of the 10 agree (§6.13); we are looking for the reason it
  would be weaker;
* **the Bullet Cluster's smaller half.** Its lensing mass is twice what our law gives it. The
  leading candidate is heat its galaxies picked up while crossing the main cluster;
* **the Sun's orbital speed**, 8% slow. A more compact Milky Way disk is the lever to test;
* **galaxy lensing in our own distances (§6.15–6.16):** the difference between ellipticals and
  spirals matches; all lenses together sit 16% above our law. Gas around them weighing as much as
  their stars, the KiDS team's own middle estimate, would close it; that gas is to be weighed. A gentler heat
  term would also close it, at a cost to massive ellipticals (§6.18);
* **Abell 1689**, the next cluster to model;
* **the heat term's mechanism (§6.17–6.19):** since rev 24 the whole chain runs in one experiment, given a quiet
  internal store in each piece of matter that motion opens; what that store is physically, and the number it sets
  for the companion's speed, are still to be found;
* a full relativistic version of the law.

Settled along the way: the main Bullet cluster's galaxy speeds (§6.3, rev 14–15), Cassini's
test (§6.10, rev 18), and the far clusters' star masses (§6.15, rev 20).

---

## 1. The problem, as plainly as we can put it

Four measurements have to be explained together:

1. **Rotation.** In a spiral galaxy, stars far from the centre should orbit slowly, the way
   Neptune orbits the Sun more slowly than Mercury. They don't; the speeds stay flat. Across
   149 well-measured galaxies (the SPARC catalogue, 3,150 speed measurements), Newton's law
   applied to the visible stars and gas misses the measured speeds by 45.6 km/s on average.
2. **Light bending (lensing).** Galaxies and clusters bend the light of objects behind them
   more than their visible matter can.
3. **Clusters.** In 12 well-measured galaxy clusters (the X-COP sample), the pull needed to
   hold the hot gas in place is up to 27 times what the visible matter supplies near the
   centre, and 4 times at the edge. Newton misses these masses by a factor of about 9.
4. **Collisions.** In the Bullet Cluster, two clusters passed through each other about 150
   million years ago. Their galaxies sailed through. Their gas, which outweighs the galaxies
   ten to fifty times, collided and was left behind in the middle. The lensing, which traces
   the pull, sits on the galaxies. A stack of 72 such collisions says the same.

There are three standard responses:
* **Dark matter.** Add invisible mass, tuned galaxy by galaxy, until the sums work. It fits
  well, but it needs about two adjustable numbers per object, and the particle has never been
  detected.
* **MOND.** Change the law of gravity at low accelerations. It predicts galaxy rotation well
  with one constant, but it fails in clusters by a factor of about 3, and it cannot put the
  collision's lensing on the galaxies.
* **Newton with visible matter only.** It fails everywhere outside the Solar System.

## 2. The ground rules we work under

The owner of this project has set four rules (recorded in [RULES.md](RULES.md)):
* no dark matter;
* no MOND or anything derived from it as an input;
* no plain Newtonian gravity as the answer;
* no expanding universe.

Every candidate formula goes through an automatic check,
[research_work/tools/formula_guard.py](research_work/tools/formula_guard.py). It asks three
questions:
* **Is it secretly MOND?** For an isolated cold mass, does it reduce to a function of the
  Newtonian pull alone? If so, how close is it to each published MOND formula? Across real
  galaxies, how much of its prediction depends on anything else?
* **Is it secretly Newton?** Is its boost just a constant, which is Newton with rescaled masses?
* **Is it secretly dark matter?** Does it need extra, adjustable, invisible mass?

We also work in a fixed order: first principles, then a toy model, then real data. When a
test fails, the first question is what would have to change to make it pass.

## 3. The idea, built from first principles

### 3.1 Matter feeds a companion to its gravity

Suppose every kilogram of ordinary matter steadily feeds a tiny amount of energy, ℓ watts
per kilogram, into a companion to its gravity. The energy streams outward at a speed u. This
is the project's founding intuition, that gravity is produced and carried rather than simply
attached. It is the only new ingredient; everything else follows from how such a stream
adds up.

### 3.2 Orderly matter adds up the way Newton does

Energy flowing away from a source spreads over a sphere, so its flux falls as 1/distance².
That is the same geometry as Newton's pull. Add up the flows from many orderly sources
arrow by arrow, as vectors, and the total flow points exactly along Newton's pull:

```
F_coherent = (ℓ / 4πG) · g_N
```

Here g_N is the ordinary Newtonian pull of all the visible matter. This is Gauss's geometry
and holds exactly; our script checks it to 2 parts in 10¹⁵. Opposite contributions cancel,
just as Newton's pulls do.

*Rev 19 tested this against the alternatives: if orderly matter's opposing flows did not cancel,
the inner parts of galaxies would spin much too fast (§4.14). Rev 20 derives it from energy
conservation (§4.15).*

### 3.3 The companion pulls with its strength

Give the companion an energy density A²/8πG, and let its pull equal its amplitude A. This
has the same form as the usual energy density of a gravitational field, g²/8πG. Energy balance through a sphere
around a mass M then says that what streams out equals what is fed in:

```
4πr² · u · A²/8πG = ℓ M   ⟹   A = √(G M a) / r,   with   a = 2ℓ / u
```

This pull falls as 1/r, not 1/r². A star orbiting at distance r then has
v² = r·A = √(G M a), which is **the same at every radius: a flat rotation curve.** It also
gives the observed rule **v⁴ = G M a**, which links a galaxy's mass to its rotation speed.

*Why the pull should equal the companion's strength is derived, as a first step, in §4.14
(rev 19).*

The constant a is not a new constant of nature here. It is twice the power per kilogram
divided by the companion's speed. With the fitted values, ℓ = 5.3 × 10⁻⁶ watts per
kilogram. Over 13 billion years, the age of the oldest stars, that is two and a half thousandths
of one percent of a kilogram's rest energy.

### 3.4 Strong gravity holds the companion back

Near the Sun the pull is billions of times stronger than in a galaxy's outskirts. If the
companion pulled there, planetary orbits would show it; they don't. So the companion must
be held back, "attached", where ordinary gravity is strong, and released where it is weak:

```
released fraction = exp( −|g_N| / g_d )
```

Here g_d = 2.03 × 10⁻¹⁰ m/s², measured on the galaxies. At every planet the released
fraction is zero to the precision of a computer. The same factor also improves the galaxy
fits.

**Breaking free takes time (rev 18).** The companion is not released the instant it leaves the
strong-pull zone. It has to travel about 0.15 parsec, 30,000 times the Earth–Sun distance or
about 900 years at its speed, to break free:

```
released share builds up as  1 − exp( −distance from its source / L ),   L = 0.15 pc
```

Cassini's radio tracking of Saturn requires this (§6.10). Even the smallest galaxies are a
thousand times larger than L, so nothing there changes.

### 3.5 Random motion scrambles the companion, and scrambled contributions don't cancel

A moving source Doppler-shifts whatever it emits. Take a population whose parts move
randomly, with a spread of speeds σ. Each part shifts its companion differently, so the
contributions drift out of step. Out-of-step contributions add as plain totals, not as
arrows, so they no longer cancel. We call the plain total the "scalar sum":

```
S = G ∫ k ρ / d²,   heat weight  k = 3σ² / u²
```

The heat weight k comes from the same Doppler effect. A moving source's received power goes
up with the cube of its Doppler factor. Averaged over random directions, that gives an extra
3σ²/u².

We checked the scrambling numerically. With emitters all in step, the time-averaged strength
is the arrow sum: 160.50, against a predicted 160.49. With emitters moving randomly at
0.3 u, it is the plain total: 210.55, against 209.20.

Why this matters is geometry. Inside a big, spread-out cloud, Newton's pulls from all sides
largely cancel, but the plain total does not. In a cluster, that is exactly the extra pull
needed: large near the centre and smaller toward the edge.

*Rev 22: joined with the pull of §4.14 in one experiment, scrambling alone makes the pull weaker, not stronger
(§6.17). The heat term stands on the data, but its mechanism is different: cold matter holds its companion back,
and free random motion releases it.*

### 3.6 Collisions switch the scrambling off

**This is the idea that solved colliding clusters.** Laboratory physics has a well-known
exception to Doppler scrambling. If an emitter changes direction many times before it has
moved one wavelength, its Doppler shifts average away. Its emission then looks like that of
an emitter at rest. The effect is called Dicke narrowing (R. H. Dicke, 1953). Its cousin,
the Mössbauer effect, is why atoms locked in a crystal emit perfectly sharp lines.

Gas particles collide, and in clusters they also spiral around magnetic fields, constantly.
Stars never collide; galaxies pass through each other. So two things at the same
temperature behave oppositely:

| | Moves freely? | Companion |
|---|---|---|
| Stars in an elliptical galaxy; galaxies in a cluster | yes | scrambled: **hot** |
| Gas and plasma, including the Sun's interior | no, constant collisions | in step: **cold** |

Our simulation shows it:
* An emitter that collides 100 times per companion wavelength stays in step about 70 times
  longer than a free one (102 time units against 1.4).
* A population of colliding emitters adds up almost like cold matter: 1.2 times the arrow
  sum, where free movers at the same temperature reach 11 times it.

So **only freely moving matter, stars and galaxies, carries the heat term.** Gas still counts
in full as ordinary mass.

### 3.7 The companion pulls along its net flow

Out-of-step waves don't add their strengths as arrows, but their energy *flows* still add
as arrows. So the net flow of the companion points along Newton's pull plus the heat-weighted
flow from the free-moving matter. The extra pull goes that way:

```
direction of the extra pull = (g_N + g_hot) / (|g_N| + |g_hot|)
g_hot = the Newtonian pull of the free-moving matter, weighted by its heat k
```

In a round, settled galaxy or cluster, both flows point to the centre, and nothing changes.
In a collision, the gas dominates Newton's pull but the galaxies dominate the heat-weighted
flow. The pull then bends toward the galaxies, and so does the lensing.

Where the two flows point in opposite directions, they partly cancel and the extra pull fades
smoothly to zero. An earlier version used the flow's direction alone, which flipped abruptly
there and would have planted thin sheets of *negative* lensing. The version above does not.

### 3.8 The whole law

```
g_N  = Newton's pull of all ordinary matter (stars, gas, everything)
S    = G ∫ k ρ_free / d²                   the scrambled total of the free-moving matter
g_hot = G ∫ k ρ_free (x′ − x) / |x′ − x|³   its heat-weighted flow
k    = 3σ²/u² for stars and galaxies,  0 for gas and plasma

extra pull  = exp(−|g_N|/g_d) · √( a (|g_N| + S) )   along   (g_N + g_hot)/(|g_N| + |g_hot|)
h = g_N + extra pull ;   ∇²Φ = −∇·h ;   the pull on matter and on light is  g = −∇Φ

around each source the companion builds up over L = 0.15 pc:  extra → (1 − e^(−r/L)) · extra
```

The last line turns the rule into a proper field that conserves energy for anything moving
through it. Mathematically it has the same shape as Milgrom's QUMOND field equation, which
we use for its structure, not its physics (§11).

**Four constants, all universal, none set per object:**

| Constant | Value | Measured on |
|---|---|---|
| a | 6.30 × 10⁻¹¹ m/s² | 149 galaxies |
| g_d (release level) | 2.03 × 10⁻¹⁰ m/s² | 149 galaxies |
| u (companion speed) | 169 km/s | 12 clusters |
| L (release length, rev 18) | 0.15 pc (at least) | Cassini, in the Solar System |

*Rev 21: measured with every data set in the project's own distances, including the 81 galaxies
whose distances come from their redshifts, and with the distance law's rate fitted to supernovae
(§6.16). Rev 20 had 6.55 × 10⁻¹¹, 2.11 × 10⁻¹⁰ and 163 km/s; revs 12–19, 6.56 × 10⁻¹¹, 2.26 × 10⁻¹⁰
and 197 km/s.*

**For round, settled objects the law is simply**

```
g = g_N + exp(−|g_N|/g_d) · √( a (|g_N| + S) )
```

That is what we fit to galaxies and clusters. For matter that has recently changed its motion,
such as colliding clusters, the companion's parts are summed where the companion now is
(§3.10).

§4 takes the law apart piece by piece: where each part may come from, and what each explains.

### 3.9 Momentum: the companion carries it

A field that streams can carry momentum, as light does. A lamp does not recoil when its light
later pushes on something far away.

We tried the alternative, making ordinary matter alone balance every push and pull. It fails:
* the stars, only about 5% of a cluster's mass, would have to absorb the reaction to steering
  the gas;
* they would be pushed outward with a third to three quarters of the force of gravity;
* cluster galaxies would then orbit as if clusters held a fraction of their mass, and they
  don't.

So we take the streaming option: gas and galaxies feel the same pull, as observed.

The cost is a small self-push on lopsided systems. For the Bullet Cluster it is 12% of the
average pull. That is 13% of the momentum the companion's stream can carry, so the books
balance.

### 3.10 The companion is slow, and it remembers

Two things follow from what we already have. Nothing new is added.

**It is slow.** It streams at u = 169 km/s, measured on clusters. At that speed it takes about
580 million years to travel 100 kpc. So most of the companion around a galaxy or a cluster
today was given off hundreds of millions of years ago.

**It keeps moving the way its source moved,** like a ball thrown from a moving train keeps
the train's speed. It has to. Galaxies move through space at hundreds of km/s, faster than the
companion itself. If the companion moved at u through a fixed background, a galaxy's gravity
would depend on how fast the galaxy happened to be moving, and galaxies would not all obey one
rotation law. They do.

For anything settled or moving steadily, this changes nothing: the companion travels with its
source. It matters only where matter has suddenly changed its motion, as in a collision:
* **The galaxies sail through,** so they keep their motion, and the old companion travels on
  with them.
* **The gas is stopped.** Around it, a fresh companion has had time to build up only out to
  u × (time since it was stopped): about 30 kpc for the Bullet Cluster.
* **So each cluster is still wrapped in the companion of the whole settled cluster it was
  before the collision,** gas, galaxies and all, centred on its galaxies.

In the law of §3.8 this affects only the companion's parts: the scrambled total S, its flow,
and the direction of the extra pull. They are summed over the companion where it is now, that
is, over the matter where it was when it gave that companion off. Newton's ordinary pull and
the release factor use the matter where it is now.

**An earlier version of this idea failed, and now it doesn't.** In round 2 the companion
streamed at 874 km/s. The stopped gas rebuilt its companion within about 300 million years,
and the lensing drifted back to the gas, which 72 observed collisions rule out. Round 3's refit,
with the stars carrying the heat, lowered u to 197 km/s for a separate reason (163 since rev
20). The fresh companion around stopped gas now grows 5.4 times more slowly, and it never catches
up with gas separating from its galaxies at 1,000 km/s or more.

## 4. The law, piece by piece: what each part is, where it may come from, and why it keeps working

The whole law fits in six lines (§3.8). This section takes it apart. For each piece we give:
* **What it is:** the maths, and what it does.
* **Where it may come from:** our best guess at the root cause. Where a guess goes beyond what
  we have tested, we say so plainly. These guesses are how we plan the next steps.
* **What it explains:** the observations that depend on it.
* **How to test the guess.**

At the end, a table shows which pieces each observation needs, and why so few pieces cover so
much.

```
g_N   = Newton's pull of all ordinary matter
S     = G ∫ k ρ_free / d²                    g_hot = G ∫ k ρ_free (x′ − x) / |x′ − x|³
k     = 3σ²/u²  for stars and galaxies,  0 for gas and plasma
extra = exp(−|g_N|/g_d) · √( a (|g_N| + S) )   along   (g_N + g_hot) / (|g_N| + |g_hot|)
h = g_N + extra ;   ∇²Φ = −∇·h ;   pull on matter and light  g = −∇Φ ;   a = 2ℓ/u
```

### 4.1 g_N: ordinary gravity, and the shape of every flow

**What it is.** Newton's pull of all the visible matter: stars, gas, everything. It falls off
as 1/distance².

**Where it may come from.** Ordinary gravity; Einstein's theory in weak fields. In our picture
it has a second meaning. Every kilogram feeds the companion at the same rate, and a flow that
spreads over spheres adds up exactly like Newton's pulls (§3.2). So the companion's orderly
flow has Newton's shape, and **ordinary gravity and the companion share one source: mass.**
That is why g_N turns up in every other part of the law: inside the square root, in the
release factor, and in the direction.

**What it explains.**
* The Solar System and the inner parts of galaxies.
* **The extra pull can be predicted from the visible matter alone.** In 149 galaxies the measured
  pull tracks the visible matter's pull point by point (Figure 1), and a galaxy's visible mass
  fixes its rotation speed. Dark matter has to be tuned galaxy by galaxy to reproduce that
  tightness. Here it is automatic, because both pulls come from the same matter.

### 4.2 ℓ: every kilogram feeds the companion, 5.3 × 10⁻⁶ watts per kilogram

**What it is.** A steady power fed by each kilogram of matter into the companion. It is tiny.
Over 13 billion years, the age of the oldest stars, it adds up to two and a half thousandths of
one percent of a kilogram's rest energy. For the Sun it is 1.1 × 10²⁵ watts, 2.8% of the Sun's
light.

**Where it may come from.** This is the project's founding intuition: energy converts into
gravity. Our working guess is that matter continuously sheds a sliver of its energy into a
gravitational companion, the way a warm object radiates.

There is a numerical clue. The project's redshift work (no expansion) has light losing energy
at a steady rate, c·α ≈ 70–75 km/s per megaparsec. Multiplied by the speed of light, that rate
is an acceleration, c·α ≈ 7 × 10⁻¹⁰ m/s². Our two acceleration constants sit close to it:
* a ≈ c·α / 11;
* g_d ≈ c·α / 3.4.

If light and matter both exchange energy with the same background, these numbers would be
linked. A similar near-match has long been noticed for MOND's constant, so this could be a
coincidence. If it is not, it is the strongest hint of where the feed comes from.

**What it explains.** It sets the strength of everything the companion does, through a = 2ℓ/u.

**How to test it.** The Sun's extra mass loss is predicted at 1.9 × 10⁻¹⁵ of its mass per year,
below what planetary tracking can see today. A theory that links ℓ to the redshift rate would
turn "a ≈ c·α / 11" into an exact prediction.

### 4.3 u: the companion's speed, 169 km/s

**What it is.** How fast the companion streams away from the matter that fed it, measured
relative to that matter.

**Where it may come from.** A slow wave. Sound moves at a speed set by the stiffness and weight
of the air; a wave on a guitar string, by its tension and weight. A companion at 0.07% of the
speed of light suggests a heavy, sluggish mode of the gravitational field. It is not an ordinary
gravitational wave: those travel at the speed of light, as the 2017 neutron-star merger showed,
and they are a separate thing.

Two clues, which may be coincidences:
* **u is the rotation speed of a big galaxy.** Our law gives a flat rotation speed of exactly u to
  a galaxy of 1.0 × 10¹¹ suns (from v⁴ = G a M), somewhat more than the Milky Way's visible mass.
* **u²/a ≈ 15 kpc**, the size of a large galaxy's disk.

Both hint that galaxies may be sized by the companion.

**What it explains.** This is the key to the whole law. Through the heat weight k = 3σ²/u², u is
the dividing line between cold and hot:

| System | Random speed of its stars or galaxies, σ | Heat weight k |
|---|---:|---:|
| Disk of a spiral galaxy | 10–30 km/s | 0.01–0.09 |
| Elliptical galaxy | 150–250 km/s | 2.4–6.5 |
| Galaxies in a cluster | 500–1,200 km/s | 26–150 |

So with a single number:
* disks come out cold, and follow the MOND-like cold limit;
* ellipticals come out warm, and bend light more than spirals;
* clusters come out hot, and get the large boosts they need.

u also sets the companion's memory after a collision: 580 million years per 100 kpc (§3.10).

**How to test it.**
* In a collision, lensing returns around stopped gas only as a sphere growing about 170 kpc per
  billion years.
* The Bullet subcluster's galaxies should move at the speeds they had before the collision.

### 4.4 a = 2ℓ/u and the square root: how strong the companion pulls

**What it is.** Energy balance. Give the companion an energy density A²/8πG, the same form as the
energy stored in an ordinary gravitational field, and let it stream outward at u. Through any
sphere, the energy flowing out must equal what the matter inside feeds in:

```
4πr² · u · A²/8πG = ℓ M    ⟹    A = √(G M a) / r = √(a · g_N),    a = 2ℓ/u
```

**Why a square root.** Energy goes as the square of a field's strength, so the strength is the
square root of the energy. The energy flow thins as 1/r² over ever-bigger spheres, so the
strength thins only as 1/r. With random motion, the energy flow also carries the scalar sum S,
which is why the law reads √(a(|g_N| + S)).

**Where it may come from.** The guess is that the companion is a gravitational field in its own
right: its energy has gravity's form, and its pull equals its strength, just as for ordinary
gravity. That is the most economical choice: nothing new is invented, only the familiar form
reused.

**What it explains.**
* **Flat rotation curves.** A star at distance r orbits at v² = r·A = √(G M a), the same at every
  distance.
* **The rule v⁴ = G M a.** It links a galaxy's visible mass to its rotation speed with a power of
  exactly 4, as observed across galaxies.
* **Where the flat part starts.** The companion equals Newton's pull where g_N ≈ a,
  6.5 × 10⁻¹¹ m/s², in the outskirts of galaxies.

**How to test it.** The power of exactly 4 in v⁴ ∝ M, and the same switch-over pull a in every
galaxy.

### 4.5 The release factor exp(−|g_N|/g_d): held back where gravity is strong

**What it is.** The share of the companion that is free to pull. It is almost zero where
ordinary gravity is strong and almost one where it is weak. The level is g_d = 2.03 × 10⁻¹⁰ m/s².

**Where it may come from.** Two guesses.
1. **Escape over a barrier, like evaporation.** Suppose the companion must climb a fixed height ℓ_d
   against the local pull to break free of its source, and its energies are spread like those of
   a warm gas, about u²/2 per kilogram on average. The share with enough energy to climb is
   exp(−2 g ℓ_d / u²). That is exactly our form, with g_d = u²/2ℓ_d, so ℓ_d ≈ 2.0 kpc. The same
   maths describes molecules evaporating from a liquid and stars leaking out of a star cluster.
2. **Screening.** In a plasma, electric charges are hidden beyond a short distance, and the
   hiding follows an exponential. Exponential cut-offs are the signature of screening.

Guess 1 ties g_d to u and to a length of a few kiloparsecs, about the size of a galaxy's core.
If that length is also the companion's wavelength, it would join this piece to the collision
rule (§4.7).

**The release takes time (rev 18).** Any escape takes time: molecules leave a liquid at a rate,
not all at once. So the released share builds up as the companion travels away from its source,
as 1 − e^(−r/L), and L = 0.15 pc is about 900 years of travel.
* **What it explains:** Cassini's measurement of the Galaxy's field inside the Solar System
  (§6.10).
* **What it costs:** one more number, bounded rather than fitted. It must be at least 0.15 pc for
  Cassini, and at most about 10 pc so that the smallest dwarf galaxies are untouched.
* **How to test it:**
  * wide binary stars: 3.5% more pull than Newton at 7,000 AU and 8% at 20,000 AU at the minimum
    length, less for a longer one;
  * a longer length would also let the dwarf galaxies' fix through (§6.13).

**What it explains.**
* **The Solar System is silent.** At the Earth, g_N is 28 million times g_d, so the released
  share is e^(−28,000,000): zero.
* **Galaxy centres behave like Newton.** Where g_N ≈ 10⁻⁹ m/s², only about 1% is released.
* **The smooth bend of rotation curves** between the Newtonian centre and the flat outskirts. The
  galaxy fits improved when this factor was added.
* **Clusters are mostly released**: 62–95% for g_N between 10⁻¹⁰ and 10⁻¹¹ m/s².

**How to test it.**
* **Wide binary stars near the Sun.** They sit in the Galaxy's pull, where our law has only about
  half of the companion released. It predicted 19% more pull than Newton beyond about 7,000 AU,
  against MOND's 43%. With the release length (rev 18) it is 9% at 20,000 AU, 8% with rev 20's
  constants (§6.7).
* **Compact galaxies.** If ℓ_d is the companion's wavelength, stars on orbits smaller than a few
  kiloparsecs should partly lose their heat.

### 4.6 The heat weight k = 3σ²/u² and the scalar sum S: random motion doesn't cancel

**What it is.** When the sources move randomly, each source's companion drifts out of step with
the others. The parts that drift add as plain totals, which never cancel, weighted by k.

**How it is derived: the Doppler effect.** A source moving toward you at speed v gets three
boosts, each by the Doppler factor D ≈ 1 + v/u:
* each wave carries more energy;
* the waves arrive more often;
* they crowd together ahead of the source.

So the received power goes as D³ ≈ 1 + 3v/u + 3v²/u². Averaged over random directions, the middle
term cancels (as many sources come toward you as go away), and **3σ²/u² is left over.** That is k.

The steady part of each source's companion keeps its timing, because on average the source is
not going anywhere, so it adds as arrows. The fluctuating Doppler excess wanders out of step, so
it adds as a plain total. Our toy simulation confirms both halves: 160.50 against 160.49
predicted for sources in step, and 210.55 against 209.20 for sources moving randomly (§3.5).

**Where it may come from.** Nothing new is needed: it is the difference between a laser and a
lamp. Light from a laser adds wave by wave, so it can cancel. Light from a lamp adds by
brightness, so it never cancels. Random motion turns part of the companion from "laser" into
"lamp".

*Rev 22: not quite. In a simulation that joins this with the pull (§6.17), a "lamp" pulls less, because a body
can't keep in step with it. What works is closer to a laser held just below its threshold: cold matter holds its
companion back, and free random motion releases it, with the energy coming from matter's own supply. Galaxy lensing
confirms the σ² (§6.17).*

**Why it matters so much in clusters.** Inside a big, spread-out cloud, Newton's pulls from all
sides cancel toward the centre, but plain totals don't. So near a cluster's centre S stays large
while g_N shrinks toward zero. That is exactly where clusters need 4 to 27 times the visible
matter's pull.

**What it explains.**
* **Cluster masses:** a 25% typical miss, against MOND's factor of 2.9.
* **Ellipticals bend light 0.17–0.27 dex more than spirals** with the same stars.
* **Galaxies with big, hot bulges.**
* **Why MOND works for disk galaxies (k ≈ 0) and fails for clusters (k ≈ 30–160).**

**How to test it.**
* At equal visible mass, systems whose stars move randomly pull harder than those whose stars
  circle in step.
* The elliptical–spiral lensing gap should stay flat beyond 100 kpc.

### 4.7 k = 0 for gas: collisions keep the companion in step

**What it is.** Gas and plasma particles collide, and in clusters also spiral around magnetic
fields, so often that their Doppler shifts average away. They feed the companion as if they
stood still. They still count in full in g_N, as mass.

**Where it comes from.** This is established physics applied to the companion: Dicke narrowing of
spectral lines (Dicke 1953), and its cousin the Mössbauer effect. In our toy, an emitter that
changes direction 100 times per companion wavelength stays in step about 70 times longer than a
free one.

**What it needs.** The companion's wavelength has to be longer than the distance a gas particle
travels between changes of direction, and shorter than the stars' orbits. For cluster gas this
means scattering by magnetic fields at least every 0.1–1 million years.

**What it explains.**
* **The Bullet Cluster.** The lensing is not on the gas, even though the gas outweighs the
  galaxies.
* **The Sun** has no heat term, so its extra mass loss is tiny.
* **Cluster cores.** Gas at 100 million degrees adds no heat; the galaxies do.
* **Ellipticals.** Any hot gas halo counts only as mass.

**How to test it.** In colliding clusters, lensing should never peak on gas alone. Abell 520 is a
test: one team reports a lensing peak on gas with few galaxies there, and another finds none.
Our law sides with none.

### 4.8 The direction: along the net flow of companion energy

**What it is.** direction = (g_N + g_hot)/(|g_N| + |g_hot|).
* Where the orderly flow and the heat-weighted flow agree, nothing changes.
* Where they oppose, they partly cancel, and the extra pull fades smoothly instead of flipping.

**Where it may come from.** Energy flows add as arrows even when the waves are out of step. In a
room lit by many lamps, the light is incoherent, yet its net flow still points away from the
brighter side. The pull follows the flow of companion energy.

**What it explains.**
* **Round galaxies and clusters are unchanged.**
* **In a collision, the pull bends toward the hot galaxies** and away from the heavy, cold gas.
  Drop this piece and the Bullet's lensing lands on the gas (§6.3).

**How to test it.** Lensing peaks in other collisions should be skewed from the galaxies toward
the gas only by the gas's own mass.

### 4.9 The field equation ∇²Φ = −∇·h: one landscape for light and matter

**What it is.** The raw pull h is turned into the slope of a single landscape Φ, a potential.
Any part of h that swirls in circles is dropped.

**Where it may come from.** Our guess is that the companion reshapes the same landscape that
ordinary gravity shapes, so everything that feels gravity rolls on it, light included. Where
the companion's flow converges (∇·h < 0), the landscape dips as if extra mass were there. That
"as if" mass is what lensing maps show, and what dark-matter analyses call dark matter.

**The swirl.** The dropped part is where "gravity that swirls", one of the project's founding
ideas, could live. It would carry spin but no pull. We have not needed it yet; it is an open
door.

**What it explains.**
* **Light and matter feel the same pull.** The six SLACS lenses agree to −0.017 ± 0.024 dex.
* **Orbits conserve energy.**
* **One potential covers rotation curves, cluster gas and lensing.**

**How to test it.** In relaxed clusters, the galaxies' orbits and the lensing should give the
same mass.

### 4.10 Momentum: the companion carries it

**What it is.** Gas and galaxies feel the same pull, and any recoil is carried off by the
companion's stream.

**Where it comes from.** Streaming fields carry momentum. Light can push a sail, and the lamp
does not recoil when its light arrives somewhere else.

**What it explains.** The galaxies in clusters orbit as the lensing says they should. If matter
alone had to balance the books, cluster stars would be pushed outward with 30–70% of the force
of gravity, which is ruled out (§3.9).

**How to test it.** Lopsided systems get a small self-push that must fit within the stream's
budget. For the Bullet Cluster it is 12% of the average pull, 13% of the budget.

### 4.11 Travel: the companion keeps its source's motion

**What it is.** Companion given off by a moving source keeps that source's velocity. S, g_hot and
the flow are summed over the companion where it actually is now.

**Where it comes from.** It is required, not added: otherwise a galaxy's gravity would depend on
how fast the galaxy moves through space (§3.10). Think of a ball thrown from a moving train.

**What it explains.**
* **The Bullet Cluster's smaller half** still bends light like the bigger cluster it used to be.
  Its measured lensing is stronger still, about that of a cluster a third the size of the main
  one, yet the galaxies that would imply are not there. Heat its galaxies picked up in the
  crossing, carried by this same rule, is the leading candidate (§6.3).
* **72 collisions:** the lensing does not follow the gas.
* **Settled systems:** no change at all.

**How to test it.**
* Lensing around stopped gas comes back at about 170 kpc per billion years.
* Older collisions show extra lensing around the smaller clump.

### 4.12 Why so few pieces cover so much

Each piece is a switch set by something measurable about the matter itself. None is a number
tuned object by object:

| Switch | Set by | Separates |
|---|---|---|
| Release factor | how strong the pull is | the Solar System from galaxies |
| Release length | how far the companion has travelled from its source | the Sun's surroundings (Cassini, wide binaries) from everything larger |
| Heat weight | how randomly the matter moves | disks, ellipticals and clusters |
| Collision rule | whether its particles collide | gas from stars |
| Direction | which way its energy flows | round systems from collisions |
| Travel | how recently it changed its motion | settled systems from collisions |

Each kind of system flips a different combination of switches. The constants never change: three
fitted to galaxies and clusters, and the release length from the Solar System.

**Which pieces each observation needs.** ● means essential; ○ means involved.

| Observation | Newton | √(a ·) companion | Release | Heat k, S | Gas cold | Direction | One landscape | Travel |
|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| Solar System and the Sun | ● | | ● | | ● | | | |
| Cassini's test; wide binary stars | ● | ● | ● (and the release length) | | ● | | ● | |
| Flat rotation curves, v⁴ = G M a | ● | ● | ○ | | | | ● | |
| The bend in rotation curves | ● | ● | ● | | | | ● | |
| Galaxies with big bulges | ● | ● | ● | ● | | | ● | |
| Ellipticals lens more than spirals | ● | ● | ○ | ● | ● | | ● | |
| Strong lenses: light = matter | ● | ● | ○ | ● | | | ● | |
| Cluster masses | ● | ● | ○ | ● | ● | | ● | |
| Bullet Cluster: lensing on galaxies | ● | ● | ○ | ● | ● | ● | ● | ○ |
| Bullet Cluster: strengths | ● | ● | ○ | ● | ● | ● | ● | ● |
| 72 colliding clusters | ● | ● | ○ | ● | ● | ● | ● | ● |

**MOND uses only the first three columns.** That is why it works where only those matter
(disk galaxies) and fails where the rest do (clusters and collisions).

### 4.13 Clues to a deeper origin

These are hypotheses, recorded because they point to what a microscopic theory must produce.

1. **a ≈ c·α/11 and g_d ≈ c·α/3.4,** where c·α is the redshift rate of light turned into an
   acceleration (§4.2). If the companion and the redshift share a cause, the ratios 11, 3.4 and
   λ = g_d/a = 3.22 are numbers that theory must derive. (Galaxy lensing now also measures α,
   §6.15.)
2. **u⁴ = G a × (8 × 10¹⁰ suns).** The companion's speed is the rotation speed of a galaxy a
   little heavier than the Milky Way's visible matter.
3. **u²/(2g_d) ≈ 2.0 kpc,** a detachment length, if the release factor is an escape over a
   barrier (§4.5).
4. **The Sun's companion carries 2.8% as much power as its light.**
5. **Breaking free takes about 900 years (0.15 pc at u).** If the release is an escape, this is
   the escape time, and theory must give it along with ℓ_d.

What would pin these down:
* a field theory that gives the companion's speed (§9);
* a measurement of its wavelength, from which systems count as hot and which as cold.

### 4.14 Where the pull comes from, and how the companion adds up: first derivations (rev 19)

Three postulates of the law say how the companion behaves:
* it pulls with its strength (§3.3);
* orderly matter's contributions cancel like Newton's pulls (§3.2);
* scrambled contributions add as plain totals (§3.5).

Rev 19 takes the first steps toward deriving them rather than assuming them, as the project's
owner proposed (§9.1, proposal 1). It also tests the alternatives on the data.

**The model.**
* The companion is a wave.
* Every piece of matter is a small, self-powered emitter of that wave, like a clock that ticks
  by itself.
* Its ticking falls into step with the wave passing it.
* The only interaction is local: each emitter is pushed by the slope of the wave where it sits.

**The result is like pushing a swing.**
* **A quarter beat ahead:** an emitter running a quarter beat ahead of the passing wave adds
  energy to it, the way a well-timed push adds to a swing. Energy added to a wave travelling
  away from its source carries momentum away with it, so the emitter recoils the other way,
  toward the source. **The pull is exactly proportional to the passing wave's strength, at every
  distance.** That is the postulate of §3.3, derived rather than assumed.
* **In step with the wave:** the pull falls as 1/distance², like Newton's.
* **A quarter beat behind:** the emitter takes energy from the wave, like catching a swing, and
  is pushed away.

We checked this with a full 3D simulation of the wave equation, which knows nothing of the
formula. From a quarter of a wavelength out to three wavelengths, the pull matches it to within
4%.

**What the microscopic model must then contain** (derived, not assumed):
1. **Matter must feed the companion, not absorb it,** like the medium of a laser. Absorbing
   matter would be pushed away.
2. **Each tiny emitter must lock at its own slightly different pitch.** If a body's emitters all
   locked together, they would emit in chorus. A galaxy's extra pull would then grow with its mass
   in the wrong way (v² ∝ M instead of the observed v⁴ ∝ M).
3. **The companion's wave crests must move much more slowly than its energy, at most 1/20 as
   fast.** The pull is paid for by the energy each emitter adds. If crests moved as fast as the
   energy, galaxies would amplify their own companion by about as much again, and v⁴ = G M a
   would not be as tight as it is.

**The pull, read locally.** The same calculation gives a general rule. Wherever an emitter sits,
its pull equals the local strength of the companion times how quickly the wave's crests sweep
past it.
* For a single wave moving away from its source, the crests sweep past at the full rate, and the
  pull is the full pull.
* Where companion arrives from several directions at once, the crests partly stand still, as in
  the standing wave on a guitar string, and the pull weakens.
* For example, add a second wave with half the strength running the opposite way. The pull drops
  to 81% of the first wave's alone, although there is now a quarter more companion.

**Three ways the companion could add up, tested.** When the companion from many pieces of matter
meets, three things could happen:
1. **the waves pass through each other**, like the light of many candles;
2. **everything merges into one stream** in the net direction, and nothing cancels;
3. **what our law has assumed since round 3:** orderly matter's opposing flows cancel, as
   Newton's pulls do; scrambled companion adds in full; and all of it pulls as one stream, in the
   net direction.

We worked out what each rule predicts and fitted each to the data, with its constants
readjusted:

| How the companion adds up | 149 galaxies' rotation (typical miss) | 12 clusters' masses (typical miss) | The Milky Way |
|---|---|---|---|
| 1. The waves pass through each other | 19.6 km/s | 47% | speeds rise outward, to 249 km/s at 30 kpc |
| 2. One stream, nothing cancels | 19.0 km/s | 25% | speeds fall steeply, to 143 km/s at 30 kpc |
| 3. **Our law** | **15.9 km/s** | **25%** | 211 km/s at the Sun, 185 at 30 kpc |

MOND scores 16.1 km/s on the same galaxies. The Milky Way numbers are quick local estimates.
*(Rev 19's numbers, with the constants and conventions of the time; rev 20's rerun in its own
conventions gives the same verdicts, §4.15.)*

* **The galaxies rule out rule 2.** Inside a disk, a star's neighbours pull from all sides. Their
  plain total is typically four times the net pull. If that total counted, the inner parts of
  galaxies would spin much too fast.
* **The clusters rule out rule 1.** In a cluster, the scrambled companion of the galaxies arrives
  from all sides. If it weakened itself the way waves passing through each other do, the extra
  pull would be 35–70% too weak, most of all near the centre.
* **Only our law passes both.**

**So three properties of the companion are now set by the data, not by our choice:**
1. **Orderly matter's opposing flows cancel,** leaving only the net flow.
2. **Scrambled companion adds its full strength.**
3. **Matter feels it as one stream whose crests sweep past at the full rate.** Waves that pass
   through each other fail this test (rule 1). So would a companion that behaves like a
   superfluid: its flows slow down when they meet head-on instead of cancelling, and the pull
   weakens in just the same way.

**Where this points.** Waves that always travel at full speed and annihilate when they meet
head-on are known in "active" media: chemical waves in some reactions, or waves of activity in heart
tissue, where two colliding fronts simply go out. In those media the medium sets the waves'
strength, whereas the companion's strength must fall off away from its source. A model that does
both is the next step, and the same galaxies, clusters and Milky Way will grade it. (Superfluid
dark matter theories are dark-matter theories, which the rules exclude, and we borrow nothing
from them.)

**One more consequence.** The companion's speed enters the law twice: in how far it travels
(memory, reach) and in the heat weight k = 3σ²/u². In this model these are two different speeds:
the speed of its energy, and a locking rate divided by its wavenumber. Measuring the travel speed
on its own (§9.1, proposal 7) therefore tests something new.

Scripts: `code/first_principles_v10.py` (the pull), `code/combination_rules_v10.py` (the three
rules).

### 4.15 Why orderly matter adds up like Newton: energy bookkeeping (rev 20)

§4.14 found, by testing, that orderly matter's companion must add up the way Newton's pulls do.
Rev 20 shows why, from energy conservation.

**The argument, in four steps.**
1. **No energy is lost on the way.** In a steady state, every watt that matter feeds into the
   companion flows outward. Through any closed surface, the power flowing out equals the power fed
   in by the matter inside.
2. **Matter feels the full pull only from one stream.** §4.14 showed that the pull is full only
   where the companion passes as a single stream whose crests sweep by at the full rate. The
   stream's energy is then its energy flow divided by its speed, and the pull goes as the square
   root of that flow.
3. **No whirlpools.** Suppose the flow also has no swirl.
4. **Then there is exactly one possible flow.** A swirl-free flow fixed by where its sources are is
   Newton's field-line pattern, Gauss's geometry. So the stream's strength follows Newton's pull,
   and the extra pull goes as the square root of Newton's pull. That is the rule the law has used
   since round 3, which §4.14 had only picked out by testing.

**The alternative breaks energy conservation.** §4.14's rule 2, "everything merges into one stream
and nothing cancels", would carry more power out of a region than the matter inside it puts in:

| Where | power carried out ÷ power fed in |
|---|---|
| inside a uniform ball, at a quarter / half / all / twice its radius | 11.8 / 5.5 / 1.5 / 1.06 |
| the Milky Way model, spheres of 4 / 8.2 / 12 / 20 / 30 kpc | 1.70 / 1.36 / 1.22 / 1.11 / 1.06 |

So it fails twice: on the galaxies (§4.14), and on energy conservation. Waves that pass through
each other (rule 1) do conserve energy, but they weaken the pull. Our law's rule is the only one
that is both lossless and undiluted.

**How tightly the data demand it.** We built one formula with three dials that contains all the
rules, and turned each dial with the constants refitted at every setting (400 random resamplings
of the galaxies and clusters for the spread):

| Dial | What it controls | 149 galaxies (typical miss) | 12 clusters (typical miss, rms) | Verdict |
|---|---|---|---|---|
| w, from 0 to 1 | how completely orderly flows cancel | 18.9 → **15.9** km/s | 0.212 → 0.221 | the galaxies demand full cancellation (at 95% the miss is already 16.1) |
| γ, from 0 to 1 | how much the pull is diluted | **15.9** → 27.2 km/s | **0.221** → 0.505 | both demand no dilution (γ ≤ 0.1) |
| h, from 0 to 1 | whether heat adds as plain totals | no change | 0.288 → **0.221** | the clusters demand plain totals, in every resampling |

The clusters lean very slightly to less cancellation (0.212 against 0.221), but only by shrinking
the constant a 7.5-fold and slowing the companion to 60 km/s. The galaxies rule that out.

**What is still assumed.** The "no whirlpools" condition. A companion that scatters off itself,
the way heat spreads by diffusion, would settle into such a flow on its own. Building that as a
working model of waves is the next step (§9, item 9). *Rev 21: built. Scattering keeps the flow
free of whirlpools but not one stream; a companion guided along gravity's field lines is both, and
needs no assumption about whirlpools at all (§6.16).*

Script: `code/companion_flow_v11.py`.

## 5. MOND, as the cold limit of our law, not an input

Switch off the heat (S = 0), as for a cold disk galaxy, and the law becomes

```
g = g_N · ν(g_N/a),   ν(y) = 1 + e^(−y/λ) / √y,   λ = g_d/a = 3.22
```

That is a MOND-type formula. **But it came out of the companion mechanism; it was not put
in:**

| MOND's ingredient | In MOND | Here |
|---|---|---|
| Its constant a₀ | a new constant of nature | 2 × power per kilogram ÷ companion speed |
| The switch between Newton and the deep regime | chosen by hand | set by the release factor. Not identical to any published choice; the nearest (McGaugh, Lelli & Schombert's exponential) is 0.028 dex (7%) away |
| Why galaxies obey it | assumed | disks are cold |
| Why clusters don't | needs extra matter | cluster galaxies are hot |
| Why collisions don't | needs extra matter | galaxies carry the heat, gas doesn't |

This is how Newton's law comes out of Einstein's in weak gravity: a special case of a larger
law. The formula check still flags the cold limit as MOND-like, as it should. Every result we
claim as ours rests on where the law is not MOND.

## 6. Tests against real data

### 6.1 Rotation curves: 149 galaxies, 3,150 measured speeds

Average miss per galaxy (km/s; lower is better). The galaxies were split in advance into
training, validation and test sets.

| | All | Training | Validation | Test | Adjustable numbers |
|---|---:|---:|---:|---:|---|
| **Ours** | **15.94** | 16.11 | **19.11** | **12.46** | 3 in total |
| MOND | 16.13 | 15.85 | 19.80 | 13.52 | 1 |
| Newton, visible matter only | 45.58 | 45.94 | 48.64 | 41.68 | 0 |
| Dark-matter halo fitted to each galaxy | 7.52 | | | | 298 |

* Ours beats MOND overall and on the held-back galaxies; MOND is slightly better on the
  training set.
* **Galaxies with big bulges** were the risk. Their bulges are hot, and the new heat rule
  makes hot stars count much more. Those 25 galaxies came out better: 29.7 km/s against
  MOND's 30.4.
* Dark matter fits more tightly because it has a hundred times more adjustable numbers.

*Rev 20's constants; revs 12–19 gave 15.85, 15.95, 19.21 and 12.40 (and 29.2 for the bulges).*

Script: `code/run_v3.py`; the regression suite for rev 20's numbers.

### 6.2 Galaxy clusters: 12 clusters with measured gas, temperatures and stars

The quantity is the typical miss in the mass needed to hold the gas, at six radii from the
core to the edge.

| | Typical miss | Adjustable numbers |
|---|---:|---|
| **Ours** | **25%** (rms 0.221 in ln M) | same 3 constants |
| Ours, clusters held out of the fit (924 splits) | 27% (0.236) | |
| Ours, stars' speeds computed by our law from visible matter alone (rev 12) | 39% (0.329) | |
| MOND | ×2.9 (1.062) | 1 |
| MOND with its constant refitted on clusters | ×1.4 (0.338); needs a constant 9.7× its galaxy value | 1 |
| Newton | ×9.1 (2.211) | 0 |
| Dark matter (NFW fit per cluster) | 11% (0.101) | 24 |

**Why the stars can do it.** The stars are concentrated in the middle, where clusters need
the most extra pull:

| Radius (fraction of the cluster's size, R500) | 0.02 | 0.05 | 0.1 | 1 |
|---|---:|---:|---:|---:|
| Mass in stars ÷ mass in gas | 1.4 to 7.5 | 0.35 to 1.25 | 0.16 to 0.33 | 0.03 to 0.06 |

Their random speeds, roughly 300 to 1,200 km/s, make each kilogram count dozens of times
over.

The companion speed comes out at u = 163 km/s. Across the 924 ways of splitting the clusters
in half, it stays between 148 and 178 km/s in 90% of them. *(Rev 21: 169 km/s, with the distance
law's rate fitted to supernovae and the galaxies, §6.16.)*

*Rev 20: the table uses the stars correctly counted inside spheres, and the clusters in our own
distances (§6.15). Revs 12–19 used the release's star masses inside circles on the sky, which made
the stars 1.3–2 times too heavy, and found u = 197 km/s (178–218) and a typical miss of 0.227
(held out 0.244). The MOND and dark-matter rows are rev 12's.*

Script: `code/run_v3.py`, with stellar mass profiles from the X-COP release (Ghizzardi et al.
2021); rev 20: `code/xcop_static_v11.py`.

### 6.3 Colliding clusters

**The Bullet Cluster, on its published numbers.** Three kinds of measurement exist.

1. **Clowe et al. (2006)** measured the gas, the stars and the lensing strength κ inside
   100-kiloparsec circles at four places:

   | Place | Gas (10¹² suns) | Stars (10¹² suns) | Measured lensing κ |
   |---|---:|---:|---:|
   | Main cluster's galaxies | 5.5 | 0.54 | at least 0.36 |
   | Main cluster's gas | 6.6 | 0.23 | 0.05 ± 0.06 (extra) |
   | Subcluster's galaxies | 2.7 | 0.58 | at least 0.20 |
   | Subcluster's gas | 5.8 | 0.12 | 0.02 ± 0.06 (extra) |

   Their paper states that its method underestimates κ in cluster cores, so the two galaxy
   values are **lower bounds**.
2. **Strong lensing** gives the total projected mass inside 250 kpc of each brightest galaxy
   (10¹⁴ suns):
   * Bradač et al. (2006): main 2.8 ± 0.2, subcluster 2.3 ± 0.2;
   * Paraficz et al. (2016): main 2.5 ± 0.1, subcluster 2.0 ± 0.2.
3. **Galaxy speeds** (Barrena et al. 2002):
   * main cluster: 1,249 ± 100 km/s from 71 galaxies;
   * subcluster: 212 ± 60 km/s from a tight group of 7.

**How we modelled it.**
* We built the gas and the stars to match Clowe et al.'s eight masses.
* We solved our field equation in 3D and projected the result exactly as a lensing map is
  made.
* We applied the same analysis the observers used.
* The stars' speeds were computed by our own law from the visible matter.
* **The companion's memory (rev 13, §3.10).** Around each cluster, the companion is that of the
  settled cluster it was before the collision.
* **The main cluster's outskirts (rev 14).** Rev 13's model stopped the main cluster's galaxies
  at 1.5 Mpc and had half the stars typical of cluster outskirts. The 7 X-COP clusters with
  measured star profiles have stars at 3.5–7.3% of the gas mass at their edge. Rev 14 lets the
  galaxies continue and adds the outer stars that the measured galaxy speed requires.
* *Rev 20: the Bullet is now computed in our own distances with the adopted constants (§6.15).
  The pattern and the main half's mass hold; the smaller half is still about half its measured
  mass. The tables below are the record of revs 12–15.*

**Where the lensing sits.** This is the pattern test.

| | Main galaxies | Subcluster galaxies | Main gas | Subcluster gas | Where the peaks are |
|---|---:|---:|---:|---:|---|
| Measured | at least 0.36 | at least 0.20 | 0.05 ± 0.06 | 0.02 ± 0.06 | on the galaxies |
| **Our law (rev 14)** | **0.61** | **0.49** | **0.08** | **0.10** | **on the galaxies** |
| Our law, without memory (rev 12) | 0.51 | 0.07 | 0.04 | 0.04 | on the galaxies |
| Old direction rule (pull along ordinary gravity only) | 0.35 | 0.13 | **0.21** | 0.07 | **on the gas** |
| Gas counted as hot (round-1 rule) | 0.46 | 0.21 | **0.34** | 0.08 | **on the gas** |
| No heat at all (MOND-like) | 0.14 | 0.07 | 0.06 | 0.04 | **on the gas** |

**How strong it is.** These are the calibrated numbers: masses in 10¹⁴ suns, speeds in km/s.

| | Main mass inside 250 kpc | Subcluster mass inside 250 kpc | Main galaxies' speed |
|---|---:|---:|---:|
| Measured | 2.5 ± 0.1 to 2.8 ± 0.2 | 2.0 ± 0.2 to 2.3 ± 0.2 | 1,249 ± 100 |
| Rev 13 (memory; lighter stars; subcluster 1/8 of the main) | 1.66 | 0.92 | 981 |
| Rev 14, galaxies not cut off at 1.5 Mpc | | | 1,150 |
| Rev 14, plus the outer stars the speed needs; subcluster 1/8 | 2.54 | 1.04 | 1,249 |
| Rev 14, subcluster 1/4 of the main before the collision | 2.49 | 1.52 | 1,249 |
| Rev 14, subcluster 1/3 of the main before the collision (not borne out in rev 15) | 2.42 | 1.93 | 1,249 |

**What this shows:**
* **The pattern is reproduced.** Both lensing peaks sit on the galaxies, 31 and 17 kpc from
  the two brightest, and the gas regions carry little extra lensing, as measured. **This is the
  pattern usually called impossible without dark matter.** Both ideas of §3.6–3.7 are needed:
  drop either one and the lensing lands on the gas.
* **The main cluster's galaxy speed is settled.**
  * Half of rev 13's gap was a modelling artefact. Our model's galaxies stopped at 1.5 Mpc,
    and real clusters don't end there. Letting them continue gives 1,150 km/s.
  * The rest needs a normal amount of stars in the outskirts: about 5% of the gas mass at the
    cluster's edge, inside the 3.5–7.3% of the X-COP clusters. Orbits that are partly radial,
    common in cluster outskirts, would need fewer.
  * In our law the galaxy speed works like a star count. The heat term makes the speed grow
    with the hot stellar mass.
* **The main cluster's lensing mass agrees too.** The same outer stars bring it to
  2.4–2.5 × 10¹⁴, against 2.5 and 2.8 measured. The speed and the lensing, two independent
  measurements, now point to the same cluster.
* **A correction to rev 13.** Rev 13 compared our 100 kpc κ values with Clowe et al.'s as if
  they were measurements. They are floors. So:
  * rev 13's "all four within about 1σ" should read "above both floors";
  * the lighter star masses it used were never required.

  The strong-lensing masses are the calibrated test, and they are what we now compare with.
* **Rev 14's answer for the smaller half: it was a bigger cluster than we assumed.** (Rev 15's
  checks, below, do not bear this out.)
  * At 1/8 of the main before the collision, our law gives it half its measured lensing mass.
  * It matches if it held **about a third of the main cluster's visible matter** (1.93 × 10¹⁴
    against 2.0–2.3).
  * Its galaxies then moved at about 800 km/s before the collision. Barrena et al. estimated
    about 700 km/s from its hot gas.
  * The lensing masses themselves are in the ratio 0.8 inside 250 kpc. Dark-matter
    reconstructions use 1/6 to 1/10 in total mass, so this is where our picture and theirs
    differ, and it can be checked.
* **Rev 14's predictions.**
  * About 7 × 10¹² suns of stars came with the smaller half. They now surround it, mostly
    counted as main-cluster galaxies, and carry its 616 km/s offset in line-of-sight speed. A
    larger spectroscopic sample can separate them.
  * About 10¹⁴ suns of its gas was stripped into the main cluster.
  * The main cluster's outskirts hold stars at about 5% of the gas mass.
* **The costs.** The gas residuals rise to 0.08 and 0.10 (within 1.3σ). The main peak moves to
  31 kpc from its galaxies, because the bigger subcluster's broad field tilts the combined map.

**Rev 15: the checks.** We tested rev 14 against two new sets of data. (Scripts:
`code/bullet_light_v6.py`, `code/bullet_members_v6.py`.)
* **A star count confirms the main cluster's outskirts.** We counted 1,652 cluster galaxies in
  the Legacy Survey, a public sky survey, chosen by colour-based distance and with the field
  level subtracted.
  * Their light keeps rising out to 2.5–3 Mpc, and 23–28% of it lies beyond 1.5 Mpc.
  * Rev 13's model had no stars there. Rev 14's, inferred from the galaxy speed alone, has 13%
    there, and its total out to 3 Mpc is within 15% of what we count.
* **The smaller half's lost stars are not there.** If it had been a third of the main, about
  7 × 10¹² suns of its galaxies would still travel with it at +616 km/s.
  * Among the 71 measured galaxy speeds outside its core, 12–15 such galaxies are expected; the
    data prefer none (about 3σ for a compact group, 1.4σ if they spread wide and fast).
  * The starlight within 250 kpc of the smaller half is 0.46 of the main cluster's central
    light. Its core alone gives 0.32; the lost galaxies at 1:3 would bring it to 1.55.
* **So the smaller half's lensing strength is open again.** In our law its lensing comes mainly
  from its fast-moving stars; more gas alone does not replace them. What would have to change,
  in order of promise:
  1. **Heat picked up in the crossing.** While its galaxies crossed the main cluster at
     3,000–4,500 km/s, their speed relative to the stars around them was several times their
     own random speed, so they fed the companion about ten times harder. By the memory rule
     (§3.10) that companion travels with them. We have not yet modelled this history.
  2. **Diffuse starlight** stripped from its galaxies and moving with it, which a galaxy count
     cannot see and deep images can.
  3. **Lost galaxies spread wider and moving faster** than assumed, which the speeds allow at
     about 1.4σ. About 450–1,900 galaxy speeds would settle it.

**Two effects that build up later.**
* **Streams passing through each other.** Two galaxy streams crossing at 3,000 km/s do look
  "hot" to each other. Switched on instantly, that would lift the smaller cluster from 0.07
  to 0.11, but it would also smear its peak. The companion is slow, though: 150 million years
  after the pass this new heat has spread only 30–60 kpc. It builds up over the next billion
  years.
* **Tidal shaking.** The main cluster's pull kicked the smaller cluster's stars by 170–600
  km/s. A fast pass kicks every star in one place the same way, which is an orderly squeeze
  rather than random motion. It turns into heat only as the orbits mix, over half a billion to
  a billion years.
* Both predict stronger lensing around the smaller clump in older collisions.

**72 colliding clusters.** Harvey et al. (2015) stacked 72 pieces of colliding clusters. The
lensing sits within 5.8 ± 8.2 kpc of the galaxies; as a fraction of the way to the gas,
−0.04 ± 0.07.
* We ran 20 pieces of simulated collision under our law with the memory, at different stages
  and mass ratios, plus a reference with the gas left on its galaxies.
* **Moving the gas 40–300 kpc away moved the lensing by −2 to +10 kpc (median 3 kpc).** As a
  fraction of the way to the gas, that is 0.03 (from −0.01 to 0.12). We agree within 1σ.
* The fresh companion around the gas grows at the companion's speed (169 km/s since rev 21;
  the stack gives β = 0.019 with it), so it never catches up with gas moving away at about
  1,000 km/s. That is exactly where round 2's version failed.

Scripts: `code/bullet_main_v5.py` (rev 14), `code/bullet_v4.py`, `code/collisions_v4.py`,
`code/bullet_speeds_v4.py`, `code/stream_tidal_v4.py` (rev 13), `code/bullet_v3.py`,
`code/collisions_v3.py` (rev 12).

### 6.4 Lensing by galaxies

**Ellipticals versus spirals (KiDS survey).** Brouwer et al. (2021) found that elliptical
galaxies bend light more than spirals with the same stars. The mean difference is 0.17 dex
(splitting by shape) and 0.27 dex (splitting by colour), at more than 5.7σ. A law that treats
all matter alike cannot produce that.
* In ours, an elliptical's stars move randomly at 150–200 km/s and never collide, so they are
  hot. A spiral's stars circle in step.
* **Our prediction: 0.17 to 0.27 dex, roughly constant beyond 100 kpc.**
* No hidden gas is needed. If ellipticals do have large hot gas haloes, the gas adds only its
  mass: 0.02–0.07 dex more for haloes of 0.3–1 times the stars. So we predict those haloes
  are modest.

Script: `code/kids_v3.py`.

**Six strong lenses (SLACS).** Each lens has a measured ring size and measured star speeds.
With standard distances and published star masses:
* **The stars need to be 1.05–1.35 times heavier than the "Salpeter" standard assumption**
  (1.37–1.85 times in our own distances with rev 20's constants, §6.15).
  Spectra of giant ellipticals already suggest star populations this heavy. Measuring these
  six lenses' star populations directly is the clean test.
* The mass needed to bend the light and the mass needed to move the stars agree to
  −0.017 ± 0.024 dex. **Light and matter feel the same pull.**

Script: `code/lenses_t35.py --constants ../run-v3/results.json`.

### 6.5 The Solar System and the Sun

* **Planets.** At every planet, the release factor switches the extra pull off entirely.
* **The Sun.** Its interior is colliding plasma, so it has no heat term. It feeds the
  companion at ℓ = 5.3 × 10⁻⁶ W/kg: 1.1 × 10²⁵ W, or 2.8% of its light. That means an
  extra mass loss of 1.9 × 10⁻¹⁵ of its mass per year, below what planetary tracking can see.

### 6.6 The formula check

* **Not MOND.** For the same Newtonian pull, random speeds of 0–1,000 km/s change the
  prediction by a factor of 10 (1.02 dex). Across the 3,150 galaxy measurements, 0.029 dex
  of the prediction depends on more than the local Newtonian pull.
* **Not Newton.** The boost varies from object to object.
* **Not dark matter.** Nothing invisible is added, and no number is set per object.
* **The cold limit is MOND-like**, and the check reports it (§5).

### 6.7 Wide binary stars: a prediction

Pairs of stars 1,000–30,000 AU apart pull on each other weakly. Beyond about 5,000 AU the pull
falls below 10⁻¹⁰ m/s², where galaxies stop following Newton. The Gaia data are disputed:
* Chae (2023, 2024) reports about 1.4 times Newton's pull at the lowest accelerations;
* Banik et al. (2024) report Newton and rule MOND out.

**Our law.** A binary is cold, and it sits in the Galaxy's own pull, 1.58 × 10⁻¹⁰ m/s² near the
Sun (the visible matter's part of a 230 km/s rotation). That holds back half of the companion.
We solved our field equation exactly for a pair in that pull and averaged over orientations.
With the Galaxy removed, the calculation reproduces the law's isolated form to 3 decimals.

| Separation (thousand AU) | 2 | 3 | 5 | 7 | 10 | 30 |
|---|---:|---:|---:|---:|---:|---:|
| **Our law: pull ÷ Newton's** | 1.000 | 1.003 | 1.08 | 1.17 | 1.19 | 1.19 |
| MOND (simple function) | 1.05 | 1.11 | 1.26 | 1.38 | 1.43 | 1.43 |
| Newton | 1 | 1 | 1 | 1 | 1 | 1 |

These are for a pair of 1.5 suns; pairs of 1–2 suns reach the same 1.19.

* **19% more pull beyond about 7,000 AU, so orbits about 9% faster**, less than half of MOND's
  effect. It switches on later, because the release factor keeps close pairs Newtonian.
* **If the Galaxy's pull near the Sun is 20% weaker or stronger**, the plateau moves to 1.26 or
  1.14.
* **Our law sits between the two published answers and agrees with neither.** That makes wide
  binaries a clean test once the analyses converge.
  * If Newton wins, the extra pull near the Sun would have to follow the Galaxy's own heat flow
    (§3.7) rather than each pair's pull, which we have not yet computed.
  * If Chae wins, more of the companion must be released near the Sun.

**Rev 16 update.** Cassini's radio tracking of Saturn (§6.10) now points to the first case. The
fix it needs, a companion that takes time to be released, lowers the wide-binary prediction to
**about 4% extra pull at 7,000 AU and 9% at 20,000 AU** (rev 17 correction: rev 16 said 1–5%,
which holds at 7,000 AU, or at 20,000 AU only for a release length of 100,000 AU).

**Rev 18 update.** The release is now part of the law. Our prediction is **4% more pull than
Newton at 7,000 AU and 9% at 20,000 AU** (for a pair weighing one sun in all, with the shortest
release length Cassini allows; less if the release is longer). If the Milky Way's hold on separate systems turns
out weaker (§6.13), it rises to 19–36%. So the binaries will measure both.

**Rev 19: the forecast is locked.** Before Gaia's fourth data release, our prediction was
committed to the repository with a fingerprint of its numbers
(`forecasts/wide_binaries_gaia_dr4_v10.json`). Pull ÷ Newton for pairs weighing about one sun:

| Separation (AU) | 3,000 | 5,000 | 7,000 | 10,000 | 20,000 | 30,000 | 50,000 |
|---|---|---|---|---|---|---|---|
| **Our law** | 1.002 | 1.02 | 1.04 | 1.05 | 1.09 | 1.12 | 1.16 |
| The no-hold route for the dwarfs (§6.13, not adopted) | 1.000 | 1.005 | 1.015 | 1.04 | 1.19 | 1.43 | 2.15 |

The two differ most beyond 20,000 AU, so the data will separate them.

**Rev 20: the adopted constants move the prediction slightly, inside the locked window.** With
u = 163 km/s the pull is 1.035 × Newton's at 7,000 AU and 1.083 at 20,000 AU. The locked forecast
counts 1.03–1.05 and 1.08–1.10 as support, so the test stands as locked.

**Rev 21: the constants move again (§6.16).** The pull is now 1.032 × Newton's at 7,000 AU and 1.076
at 20,000 AU. The second sits just under the window's edge, well inside the forecast's allowance
for the analysis errors; what would refute it (below 1.02 or above 1.2) is far away. The forecast
stands as locked.

Script: `code/wide_binaries_v6.py`.

### 6.8 The Milky Way, measured star by star

Our own galaxy is the one measured in the most detail, so it is the hardest test. We built its
visible matter from published maps (McMillan 2017: a central bulge, two stellar disks, two gas
disks, plus the thin halo of old stars) and ran our law on it with nothing adjusted.

**How the law sees the Milky Way.**
* **Inside about 3 kpc** the pull is stronger than g_d. The companion is held back and gravity is
  nearly Newton's.
* **At the Sun** the ordinary pull is 1.3 × 10⁻¹⁰ m/s², just under g_d. About half the companion
  is released, so the Sun sits right in the switch-over.
* **Beyond about 15 kpc** it is fully released. The extra pull then fades slowly, keeping the
  rotation high.
* **Heat** comes only from the bulge and the halo stars. They are a small part of the Galaxy, so
  the heat adds just a few percent.

**What we compared against:**
* the rotation speed from 5 to 27 kpc, measured from hundreds of thousands of stars by the
  Gaia satellite, in four independent analyses (Eilers et al. 2019, Zhou et al. 2023, Ou et
  al. 2024, Jiao et al. 2023);
* the pull towards the disk 1.1 kpc above the Sun (Bovy & Rix 2013; Holmberg & Flynn 2004);
* the mass inside 20, 50, 100 and 200 kpc, from star streams, globular clusters and satellites
  (five studies);
* the escape speed at the Sun (seven studies);
* the inner Galaxy's share of pull from its own stars, from star counts and microlensing;
* ten small galaxies in orbit around it.

| Measurement | Measured | Our law | MOND | Dark-matter model |
|---|---|---|---|---|
| Rotation at 20 kpc | 200–207 km/s | **198** | 206 | 223 |
| Rotation at 25–27 kpc | 173–201 km/s | **191–193** | 201–202 | 220–221 |
| Rotation at the Sun | 229–234 km/s | 211 | 223 | 234 |
| Pull above the disk (as a surface density) | 68 ± 4 and 74 ± 6 M☉/pc² | **73** (with the matter counted near the Sun) | 84 | 74 |
| Mass inside 100 kpc | 6.1–7.3 × 10¹¹ suns | **6.7** | 8.0 | 8.5 |
| Mass inside 200 kpc | 11.0 +2.7/−2.2 × 10¹¹ | **12.7** | 15.6 | 12.8 |
| Mass inside 20 kpc | 1.9–2.1 × 10¹¹ | 1.7 | 1.8 | 2.2 |
| Escape speed at the Sun | 445–580 km/s (seven studies) | **509–525** | 546–564 | 560–570 |
| Share of the inner pull from stars | 0.88 ± 0.07 | 0.97–0.99 | 0.91–0.93 | 0.88 (built in) |

The dark-matter model here is McMillan's (2017), which was fitted before the newest Gaia data.
Our column uses rev 20's constants; the Milky Way is nearby, so they change it by at most 1%.

**What it means.**
* **Far out, our law runs right along the Gaia measurements**, including the gentle decline the
  newest Gaia data show. The standard dark-matter model of the Galaxy stays at about 220 km/s
  there, too fast by 6–11%.
* **Near the Sun our law is 8% too slow.** This is not special to the Milky Way. At the same
  strength of pull, the typical SPARC galaxy sits 0.03 dex (7%) above our law too. The switch,
  as fitted, holds the companion back a little too much at this strength of pull.
  * **What would fix it: a gentler hold.** Setting g_d 1.5–2 times higher releases more of the
    companion at the Sun and speeds it up to 218–221 km/s.
  * The 149 SPARC galaxies are just as happy with that: their typical miss goes from 15.85 to
    15.80 km/s.
  * Adding 25% more stars to the disk would also work, but the stars counted near the Sun
    don't allow it.
* **The pull towards the disk is right** when we use the stars and gas actually counted near
  the Sun: 74 against 68–74. MOND gives 84.
  * Rotation and vertical pull together say the extra pull near the Sun should lean a little
    more towards the Galaxy's centre. Our law already bends the extra pull towards the flow
    from hot, randomly moving stars. In the Milky Way those stars are only a tenth of the
    matter, so the lean is small. This is the "direction-dependent gravity" idea, now with a
    number attached.
* **The Galaxy's total mass out to 100 and 200 kpc agrees with our law.** MOND's is 10–40% too
  high there, because its extra pull is stronger at the lowest pulls.

**Ten small galaxies around the Milky Way.** Dwarf galaxies hold a few hundred thousand to
twenty million stars, 76–260 kpc away. We predicted how fast their stars should move, taking
into account the Milky Way's own pull on them:

| Dwarf | Measured (km/s) | Our law | MOND |
|---|---|---:|---:|
| Fornax | 11.7 ± 0.9 | **12.4** | 13.6 |
| Leo I | 9.2 ± 1.4 | **9.3** | 10.2 |
| Leo II | 6.6 ± 0.7 | 5.3 | 5.9 |
| Sculptor | 9.2 ± 1.4 | 6.7 | 7.5 |
| Carina | 6.6 ± 1.2 | 3.5 | 4.0 |
| Sextans | 7.9 ± 1.3 | 2.1 | 2.4 |
| Draco | 9.1 ± 1.2 | 2.8 | 3.2 |
| Ursa Minor | 9.5 ± 1.2 | 3.5 | 3.9 |
| Crater II | 2.7 ± 0.3 | 1.0 | 1.2 |
| Antlia 2 | 5.7 ± 1.1 | 1.1 | 1.3 |

* **The four bright or distant ones agree** within one or two error bars.
* **The six faint ones near the Galaxy, or very spread out, come out 1.5–5 times too slow.**
  MOND, calculated the same careful way, has the same problem.
  * Published MOND numbers for these dwarfs used a rougher formula that gives about 1.7 times
    higher speeds.
* **Most of the shortfall is the Milky Way's hold.** In our law, as in MOND, a dwarf sitting in
  the Galaxy's field gets much less extra pull.
  * With that hold removed, Crater II comes out at 3.4 km/s (measured 2.7) and Antlia 2 at 4.2
    (measured 5.7).
  * Draco and Ursa Minor would still need about twice their stars. They are hard for MOND too,
    and are thought to be disturbed by the Galaxy's tides.
* **A natural reason the hold could be weaker.** A dwarf moves past the Galaxy at 100–300 km/s,
  about the companion's own speed. Its companion and the Galaxy's may then not add in step.
  That is our next calculation.

Scripts: `code/milky_way_v7.py`, `code/mw_dwarfs_v7.py`, `code/transition_check_v7.py`; the
solver is `code/mw_model.py`.

### 6.9 Lensing: the standard examples

* **Light bending at the Sun:** exactly Einstein's 1.75 arcseconds. The companion is switched
  off there.
* **Microlensing towards the bulge**, when a star passes in front of another. The lensing star's
  pull at its "Einstein ring" (about 3 AU across) is a million times g_d, so every event is
  standard. The rate simply counts stars, as observed.
* **The Einstein Cross**, a quasar split into four by a galaxy's bulge. The ring is 0.65 kpc
  from the centre, where the pull is 23 times g_d. Our law says the lensing mass inside it is
  just the stars. Measured: dark matter at most 15–20%, consistent with none.
* **Six strong lenses (SLACS):** light and matter agree (§6.4).
* **Galaxy lensing out to 3 Mpc (KiDS-1000; Brouwer et al. 2021).** Around about 259,000
  isolated galaxies, the lensing measures the average pull. We put our law through exactly the
  paper's steps, using its public data. In our law spirals are cold and ellipticals are hot:

| Lenses | Our law | MOND |
|---|---|---|
| All | **+0.02 dex** | +0.05 dex |
| Spirals (blue) | **−0.005 dex** | −0.07 dex |
| Ellipticals (red) | **+0.02 dex** | +0.11 dex |

*Median of log(measured ÷ predicted) where the isolation is reliable; zero is a perfect match.
These use the paper's standard distances. In our own distances (rev 20, §6.15) they read +0.06,
+0.08 and +0.04: spirals then lens 20% more than our law predicts, by an amount our distance law's
scale decides.*

* **Our law matches both spirals and ellipticals.** MOND's single curve is too high for spirals
  and too low for ellipticals.
  * We had worried that our law's weaker deep-regime pull would under-predict galaxy lensing.
    It doesn't: spirals sit exactly where the cold law puts them.
* **Flat lensing speeds to about a million light years.** Mistele et al. (2024) find the circular
  speeds from lensing stay flat to about 1 Mpc. So do ours.
  * **A prediction:** our law's stay flat to the companion's reach, 1.7–2.2 Mpc, and then fall.
  * For ellipticals our speeds match theirs (184 / 221 / 270 against 197 ± 12 / 211 ± 8 /
    266 ± 5 km/s).
  * For the lightest spirals we are 20–30% lower than their numbers, but not lower than
    Brouwer et al.'s conversion of the same data. So the two conversions need reconciling first.
* **Colliding clusters beyond the Bullet** (Abell 520, MACS J0025, El Gordo): now modelled, in
  §6.11. The dense cluster Abell 1689 is next; its published measurements are collected.

Script: `code/lensing_census_v7.py`.

### 6.10 The precision tests: planets, pulsars, the Galactic Centre, and Cassini

* **Planets.** The extra pull is switched off by a factor 10^(−12,590) at Neptune and
  10^(−7.6×10⁷) at Mercury. Even at Sedna's far point, 937 AU out, it is 10⁻¹⁴ of Newton's.
* **The star S2 around the Galaxy's central black hole.** Its orbit shows Einstein's
  precession, 1.10 ± 0.19 (GRAVITY 2020). Our law gives exactly Einstein's value: the pull there
  is 30 million to 8 billion times g_d.
* **Binary pulsars.** Einstein's orbital decay, exactly. Our one new effect, the companion's
  tiny mass loss, is 0.3% (Hulse–Taylor) and 1.3% (Double Pulsar) of the measurement error.
* **Cassini's test of the Galaxy's field inside the Solar System: failed by the law as it stood;
  passed since rev 18.**
  * **What it is.** In laws like ours and MOND, the Galaxy's pull slightly reshapes the Sun's own
    extra pull thousands of AU out. That reshaping reaches the planets as a tiny stretching
    force, called Q2.
  * **What was measured.** Nine years of radio tracking of the Cassini probe at Saturn found
    Q2 = (3 ± 3) × 10⁻²⁷ per second squared (Hees et al. 2014). Their work excluded MOND with its
    usual forms.
  * **What the law gave (rev 16).** Computed exactly, **2.4–3.1 × 10⁻²⁶**, about 10 times too
    big. MOND gives the same.
  * **What had to change: the release takes time.** The companion leaving the Sun breaks free
    gradually as it travels. The released share builds up as 1 − e^(−distance/L) with
    L = 0.15 parsec, 30,000 times the Earth–Sun distance, about 840 years at 169 km/s.
    * **Q2 drops to 4.6 × 10⁻²⁷** (4.5 with rev 20's constants, 4.4 with rev 21's), inside
      Cassini's measurement.
    * Every galaxy, dwarf and cluster is thousands of times larger than that, so none of them
      changes. The full regression suite confirms it: Cassini moves from fail to pass, none of
      the other 75 graded checks changes its grade, and the only other number that moves is the
      wide-binary prediction (§6.12).
    * Wide binaries get 4% more pull at 7,000 AU and 9% at 20,000 AU, instead of 19% (3.5% and
      8% with rev 20's constants, §6.7).
  * **Adopted in rev 18.** This is "gravity that builds up over time", one of the project's
    founding ideas, now with a measured length.
    * It costs one number, the release length: at least 0.15 parsec (Cassini), at most about
      10 parsecs (the smallest dwarf galaxies). We use the minimum.
    * Checked against the dwarfs directly: even at 1 parsec, at most about 0.5% of a dwarf's
      companion is still held back at its half-light radius (`code/release_hold_v9.py`).

Script: `code/strong_field_v7.py`.

### 6.11 Three more colliding clusters

*Rev 19: this section keeps round 8's results, computed with the expanding universe's
distances. §6.14 redoes all three with the project's own. Rev 20 (§6.15): with the companion's
speed measured in our own distances too, they no longer need heavier stars, and El Gordo's lensing
agrees with its published ones.*

The Bullet Cluster is the famous collision, but it is one object. Round 8 runs our law on the
three other collisions that astronomers have mapped best, with nothing adjusted:
* the gas and star masses come from the papers;
* the machinery is the Bullet's (§3.10): the old companion rides with each cluster's galaxies,
  and a fresh one grows around the stopped gas at the companion's speed (197 km/s in rev 17;
  163 since rev 20).

Every input and its source is listed in the script (`collisions_v8.py`).

**MACS J0025.4−1222: a second Bullet, and a clean pass.** Two near-equal clusters crossed in the
plane of the sky about half a billion years ago at about 2,000 km/s. Their gas stuck in the
middle; their galaxies flew on, 370 and 170 kpc beyond it.

| | Measured (Bradač et al. 2008) | Ours |
|---|---|---|
| Lensing mass within 300 kpc of the SE galaxies | 2.5 (+1.0 / −1.7) × 10¹⁴ suns | **2.0** |
| Lensing mass within 300 kpc of the NW galaxies | 2.6 (+0.5 / −1.4) × 10¹⁴ | **1.6** |
| Where the lensing peaks sit | on the galaxies, > 4σ from the gas | **15 and 43 kpc from the galaxies** |
| Galaxies' speed spread | 835 ± 59 km/s | **770** |

All four agree within the measurement errors, and the answer does not depend on the uncertain
age of the collision (0.26–1 billion years give the same numbers).

**Abell 520: the "dark core", explained by ordinary gas.** This "train wreck" has five or six
clumps in a line. In its middle sits a lensing clump, P3, on top of the stopped gas but with
very few galaxies. Two teams measured it differently: 3.35 ± 0.34 (Jee et al. 2014) and
2.84 ± 0.64 (Clowe et al. 2012, unsmoothed) × 10¹³ suns inside 150 kpc. Taken at face value it
seemed to need dark matter that had separated from the galaxies, which dark-matter models do not
expect.

| Mass inside 150 kpc (10¹³ suns) | P1 | P2 | P3 (dark core) | P4 | P5 | P6 |
|---|---|---|---|---|---|---|
| Jee et al. 2014 | 2.10 | 4.05 | 3.35 | 4.23 | 2.93 | – |
| Clowe et al. 2012 | 2.81 | 4.16 | 2.84 | 5.59 | 3.17 | 3.68 |
| **Ours** | **2.56** | **4.84** | **3.06** | **3.30** | **2.83** | 2.10 |

* **Our law puts a lensing clump on the dark core, with a mass between the two teams' values.**
  * There is 0.76 × 10¹³ suns of ordinary matter there, mostly gas.
  * The companion multiplies it about four times, as it does for any ordinary matter.
  * The heat of the hot galaxies all around adds to it, since heat never cancels.
  * Nothing dark is needed, and no galaxies either.
* Five of the six clumps agree within the errors, counting the gap between the two teams as
  part of the uncertainty. P6 is 2.3 error bars low.
* The whole cluster inside 710 kpc: 4.3 against 5.0 ± 0.55 × 10¹⁴ suns (Mahdavi et al. 2007).
* The galaxies' speed spreads come out about a third below those measured near each clump
  (from only 6–9 galaxies each): 430–590 against 580–810 km/s.

**El Gordo: a very distant giant, where the stars' mass decides everything.** El Gordo is
seen as it was 7 billion years ago (z = 0.87). It is one of the most massive clusters known:
two halves 750 kpc apart that passed through each other about half a billion years ago at
2,400 km/s.

| | Measured | Ours, published star masses | Ours, stars × 2 |
|---|---|---|---|
| Lensing mass within 0.5 / 1 / 1.5 Mpc (10¹⁴ suns) | 6.7 / 15.8 / 22.6 (Kim et al. 2021) | 4.6 / 11.3 / 17.0 | **7.0 / 16.9 / 24.8** |
| NW / SE galaxies' speed spread | 1,290 ± 134 / 1,089 ± 200 km/s | 944 / 839 | **1,187 / 1,047** |

* With the published star masses, the lensing mass and the galaxy speeds both come out 25–30%
  low.
* **Twice the published star masses recovers both at once**, within 10% for the lensing and
  within one error bar for the speeds.
* The published values come from colour fitting, which its authors say is good only to a factor
  of two. At the star-to-gas ratio of the X-COP clusters that calibrated u, El Gordo's stars would
  be about 1.5 times the published values.
* The same doubling would overshoot MACS J0025. So it is not a rule for every cluster. It says
  that El Gordo's stars are the number to measure better.
* **Rev 18 caution.** These El Gordo numbers, like the published masses they are compared with,
  use the expanding universe's distances. With the project's own distance law, El Gordo is
  1.58 times larger, its gas 1.68 times heavier, its stars 0.71 times as heavy and its lensing
  mass 1.60 times higher. The comparison has to be redone before "twice the stars" can be
  trusted (§9, item 1). *Redone in rev 19 (§6.14): against the masses measured inside fixed
  circles, El Gordo needs 1.6–2 times its published stars.*
* Our SE lensing peak sits 60 kpc from the SE galaxies, on the side facing the centre. That
  puts it 160 kpc from the cool gas core, against about 100 kpc measured. This depends on where
  the bulk of the gas sits, which no paper maps, so the suite tracks it without grading it.

### 6.12 Every test, every time: the regression suite

A change that fixes one measurement can quietly break another. So every test the law has faced
now lives in one program: the regression suite.
* It runs the law through **89 checks**. 76 are graded against published measurements:
  galaxies, clusters, lensing, the Milky Way, its dwarfs, the Solar System, the Bullet Cluster,
  the 72-collision stack, and the three collisions above. The other 13 are numbers it tracks,
  such as the law's constants.
* Each check is graded the same way:
  * **pass**: within two error bars;
  * **close**: within three;
  * **fail**: further out.
* Every run is compared with a saved **baseline**, and the report lists what changed:
  * what a change **fixes**;
  * what it **breaks**;
  * what moved without changing its grade.
* The quick run takes about a minute; everything, including the collisions, takes about 20
  minutes.
* On an unchanged law it reproduces every number published so far to the last digit.

**The law today (rev 21): 59 pass, 11 close, 7 fail** of the 77 graded checks (rev 20: 57, 9, 10 of 76;
rev 19: 58, 11, 7; rev 18: 62, 8, 6; rev 17: 61, 8, 7). §6.15 and §6.16 have what moved and why.
* **Fail (7):**
  * five faint dwarf galaxies: Draco, Ursa Minor, Sextans, Crater II and Antlia 2;
  * the lensing mass of the Bullet's smaller half;
  * galaxy lensing around red, elliptical-like lenses (KiDS), in our distances.
* **Close (11):**
  * the dwarf Carina;
  * the Sun's orbital speed;
  * the Milky Way's mass inside 50 kpc;
  * galaxy lensing of all lenses together and around disks (KiDS), and the lensing speeds of spirals
    and of ellipticals (Mistele et al.), all in our distances;
  * MACS J0025's galaxy speeds;
  * Abell 520's clump P6 and its galaxy speeds;
  * El Gordo's NW galaxies' speeds.

**The first thing we did with it: test the fixes proposed in rev 16.**

| Change | Fixes | Breaks |
|---|---|---|
| Release over 30,000 AU (for Cassini) | Cassini | nothing (adopted in rev 18) |
| Weaker hold by outside galaxies, 10% (for the dwarfs) | Cassini, Carina; every dwarf improves | wide binaries would pull 2.5× Newton's, more than any analysis allows |
| Later switch-off, g_d × 1.25 | nothing yet (the Sun 211 → 215 km/s) | the pull above the disk, just |
| Later switch-off, g_d × 1.5 | the Sun's speed (218 km/s) | strong lenses and the pull above the disk move to "close" |

What we learned, in a few minutes of computer time:
* **The Cassini fix is free.** Rev 18 adopts it (§6.10).
* **The dwarf fix works, but the hold must weaken only for systems moving past the Galaxy's
  companion.** Dwarfs move past it at 100–300 km/s; a pair of stars moves with it. So the hold
  should depend on speed, not be a single number.
  * Rev 18 found a simpler way out: the adopted release length, if about half a parsec or
    longer, protects the binaries whatever the hold (§6.13).
* **The Sun's speed is better fixed in the Milky Way's matter than in the law.** A later
  switch-off helps the Sun but lifts the pull above the disk too. A more compact disk, as in
  Bovy & Rix (2013), gives 217 km/s at the Sun with less matter near the Sun, not more.

### 6.13 What the Cassini fix opens: the faint dwarf galaxies

Six of the ten dwarf galaxies around the Milky Way come out too slow, five of them by a factor
of 2.7–5. Rev 18 found why. It also found that the Cassini fix removes what stood in the way of
fixing them.

**Why they are slow.** Inside a dwarf, our law's extra pull depends on the square root of the
total pull there, and that includes the Milky Way's. For Draco the Milky Way's pull at the
half-light radius is twice the dwarf's own; for Sextans 9 times, for Crater II 54 and for
Antlia 2 82. Under a square root a big steady pull swamps a small one, like a whisper added to a
shout. The dwarf's own companion counts for little, so its stars move slowly. MOND has the same
effect; in our law it is stronger.

Taking parts of the Milky Way out of each dwarf's formula shows which part is to blame (χ² over
the ten dwarfs; lower is better):

| The Milky Way in each dwarf's formula | χ² |
|---|---|
| all of it (the law as it stands) | 135 |
| its pull, without its heat | 133 |
| its heat, without its pull | 80 |
| neither | **60** |
| *MOND with the same stars, for comparison* | *119* |

**What stood in the way.** Anything that weakens the Milky Way's hold on a dwarf also weakens
it on the Sun and on pairs of stars. Without the hold, a pair 20,000 AU apart would pull three
times harder than Newton's law allows, and no one has seen that.

**What the Cassini fix changes.** The companion now builds up over a short distance around each
star. If that distance is half a parsec or more, still far too small to touch any galaxy, pairs
of stars stay inside the measured range even with no hold:

| Release length | Pair at 20,000 AU, full hold | Pair at 20,000 AU, no hold |
|---|---|---|
| none (the law before rev 18) | 1.19 | 2.97 |
| 0.15 parsec (adopted) | 1.09 | 1.96 |
| 0.5 parsec | 1.03 | 1.36 |
| 1 parsec | 1.02 | 1.19 |

*Pull ÷ Newton's. The two published analyses find 1.0 and 1.4–1.5.*

**The dwarfs with no hold** (km/s):

| Dwarf | Measured | Full hold | No hold |
|---|---|---|---|
| Draco | 9.1 ± 1.2 | 2.8 | 4.1 |
| Ursa Minor | 9.5 ± 1.2 | 3.5 | 4.2 |
| Sextans | 7.9 ± 1.3 | 2.1 | 4.5 (close) |
| Carina | 6.6 ± 1.2 | 3.5 (close) | **4.4** ✓ |
| Crater II | 2.7 ± 0.3 | 1.0 | 3.4 (close) |
| Antlia 2 | 5.7 ± 1.1 | 1.1 | **4.2** ✓ |
| Fornax, Leo I, Leo II, Sculptor | 6.6–11.7 | agree | agree |

* **6 of the 10 agree instead of 4, and two more are close.** χ² falls from 135 to 60, half of
  MOND's 119 with the same stars.
* **The regression suite** runs this as the candidate `no_hold` (no hold, release over about
  1 parsec): 64 pass, 9 close and 3 fail in all, from 62, 8 and 6. Carina and Antlia 2 move to pass, Sextans and Crater II to close, and no check drops a grade. Cassini's Q2 falls to zero, still a pass, and wide binaries sit at 1.19.
* **Why it is not adopted yet: we have not found the reason.**
  * Our law already says how the companions of galaxies moving at random inside a cluster
    combine: their contributions add up inside one square root.
  * So "each system pulls on its own" is ruled out: it would change every cluster.
  * The real reason has to weaken the Milky Way's hold on a dwarf while leaving galaxies and
    clusters alone.
  * The ten dwarfs already narrow it down. Crater II wants a little hold: 2.05 km/s at 5% of the
    hold, 3.4 with none, against 2.7 ± 0.3. Antlia 2 wants none.
* **Draco and Ursa Minor** stay about half as fast as measured, whatever the hold. They are also
  the hardest cases for MOND. Both are among the closest of the ten, and a dwarf without dark
  matter is easily stirred by the Milky Way's tides, so that is the next thing to check.

Script: `code/release_hold_v9.py`.

### 6.14 The far clusters and the strong lenses in our own distances (rev 19)

Rev 18's audit found that our tests of distant objects had borrowed the expanding universe's
distances (§9, item 1). Rev 19 takes them out.

**What changes.** The project's own distance law is static. Redshift comes from light slowly
losing energy on its way: 1 + z = e^(αD). Using it instead of the expanding universe's law to turn
what is measured (angles, brightness, X-rays, lensing) into sizes and masses gives, at the same
angle and brightness:

| | z | sizes | star masses (from light) | gas masses (from X-rays) | measured lensing masses |
|---|---|---|---|---|---|
| Abell 520 | 0.20 | ×1.08 | ×0.80 | ×1.00 | ×1.10 |
| MACS J0025 | 0.59 | ×1.36 | ×0.73 | ×1.36 | ×1.38 |
| El Gordo | 0.87 | ×1.58 | ×0.71 | ×1.68 | ×1.55 |

Each lensing mass is converted with its own paper's assumptions, both about the universe and
about how far away its background galaxies are. El Gordo is now compared with the lensing mass
measured inside fixed circles (Kim et al. 2021), not with their fit of two dark-matter-shaped
profiles.

**Results** (ours ÷ measured for lensing):

| | Expanding-universe distances (rev 17) | Our distances (rev 19) |
|---|---|---|
| MACS J0025: lensing inside 300 kpc, SE / NW | 0.78 / 0.61 | 0.64 / 0.53, both within the errors |
| MACS J0025: galaxies' speed spread (835 ± 59 km/s) | 770 | 707, 2.2 error bars low |
| Abell 520: six clumps inside 150 kpc | 5 of 6 within the errors | 4 of 6 (P4 2.1 and P6 2.6 error bars low) |
| Abell 520: the galaxy-poor "dark core" P3 | 3.06 against 2.84–3.35 × 10¹³ suns | 2.99 against 3.13–3.69: still agrees |
| Abell 520: inside 710 kpc | 0.87 | 0.76 |
| El Gordo: lensing inside 0.5 / 1 / 1.5 Mpc | 0.69 / 0.72 / 0.75 | 0.75 / 0.75 / 0.86 |
| El Gordo: galaxies' speed spread, NW / SE (1,290 / 1,089 km/s) | 944 / 839 | 917 / 822 |

**The regression suite now grades these in our distances.** Its tally moves from 62 pass,
8 close, 6 fail to **58 pass, 11 close, 7 fail**. All five changes are far-cluster checks, and
all move the same way.

**Why: the stars.** In our distances the same starlight means 20–29% less star mass, while the
measured lensing masses rise by 10–55%. In our law the freely moving stars and galaxies carry the
heat term, which makes most of a cluster's extra pull. So lighter stars mean less lensing and
slower galaxies.

**How heavy would the stars have to be?** We reran all three clusters with every star mass
raised by one common factor:

| Star masses × | MACS J0025 | Abell 520 | El Gordo |
|---|---|---|---|
| 1 (published, in our distances) | galaxy speeds 2.2 error bars low | P4, P6, the 710-kpc mass and the speeds low | lensing 0.75–0.86; speeds low |
| **1.4** | **everything agrees** | everything but P6 and the speeds (both close) | **lensing agrees**; SE speeds agree, NW close |
| 1.8 | everything agrees | P2 now too heavy; the speeds still close | everything agrees |

* With 1.4 times the stars in all three, the suite would read **63 pass, 7 close, 6 fail**,
  better than before the switch.
* A factor of 1.4 roughly undoes the conversion, which made the stars 0.71–0.80 times as heavy.

*Rev 20: both the factor 1.4 and the "older stars" reading below are withdrawn. The calibrating
clusters' stars had been overcounted, and the companion's speed measured in the standard
distances; measured consistently, the far clusters need no heavier stars. And distant cluster stars
are observed to be younger, not older (§6.15).*

**Is that plausible? (Rev 19's reasoning, kept for the record.)** The published star masses
carry assumptions of their own:
* El Gordo's come from colour fits that did not allow stars older than 7 billion years, the Big
  Bang's age at that distance;
* MACS J0025's use one fixed ratio of mass to light;
* Abell 520's use a ratio of 2 that we took from the Bullet Cluster, on the low side for old
  stars.

Older stars weigh more for the same light. In standard models of starlight, stars 10–13 billion
years old, like those in nearby clusters, weigh roughly 1.3–1.7 times more for the same light
than stars 6–7 billion years old. Without the Big Bang's timeline there is no reason for the
distant clusters' stars to be young. Refitting their star masses from the photometry with no age
cap (§9, item 1) will tell.

**The six strong lenses.** Light and matter still agree: −0.012 ± 0.023 dex with our distances
(−0.017 ± 0.024 with the standard ones). Their stars must be heavier, though: 1.44–1.95 times
Salpeter's rule instead of 1.05–1.35. Our distances again give less star mass for the same light
(1.25 times less; rev 19 said 1.48, the energy-loss-only figure, while the lens code itself used the
right one) and more lensing mass (1.08–1.17 times). This is prediction 1 (§8), now in our own
distances.

Scripts: `code/collisions_v10.py`, `code/collisions_star_sweep_v10.py`; the lenses in
`regression/t_lensing.py`.

### 6.15 Everything in our own distances, and the calibrating clusters' stars corrected (rev 20)

Rev 19 moved the far clusters and the strong lenses into the project's own distances. Rev 20 moves
everything else that depends on distance: the Bullet Cluster, galaxy lensing (KiDS), and the 12
clusters that set the companion's speed. On the way, a check of where every cluster's star masses
come from found a mistake in our own inputs, dating from round 1.

**The mistake: the calibrating clusters' stars were counted in projection.**
* The X-COP release gives each cluster's star mass inside a circle on the sky, that is, inside a
  long cylinder through the cluster, which includes stars in front of it and behind it. The
  release's paper says so, and converts to the mass inside a sphere by multiplying by 0.75.
* We had used the cylinder masses as sphere masses since round 1. That made the stars 1.3 times
  too heavy at the clusters' edges and up to about twice too heavy in their middles.
* Corrected with the paper's own galaxy profile, our sphere masses match its published table to
  within 1–13%. The star-to-gas ratio at the edge becomes 2.7–5.6%, not 3.5–7.4%.

**The companion's speed, remeasured.** The clusters' stars carry the heat term, whose weight is
3σ²/u². Fewer stars at the same speeds must be balanced by a slower companion. With the clusters
also converted into our distances (radii ×0.97–1.00, gas ×0.89–0.92, stars ×0.84–0.86):

| X-COP clusters | companion speed u | typical miss (rms) |
|---|---|---|
| as used in rounds 1–19 (standard distances, stars in projection) | 197 km/s | 0.227 |
| stars corrected | 179 km/s | 0.221 |
| **stars corrected, our distances (adopted)** | **163 km/s** | **0.221** |

Refitting all three constants together gives a = 6.55 × 10⁻¹¹ m/s², g_d = 2.11 × 10⁻¹⁰ m/s² and
u = 163 km/s, with the 149 galaxies at 15.9 km/s. **Rev 20 adopts these.** It is the same law, with
its constants now measured in the project's own distances on correctly counted stars. The numbers
that follow from u move with it: the companion's reach is 2.2 Mpc (was 2.6), breaking free takes
about 900 years (was 700), and each kilogram feeds the companion 5.3 × 10⁻⁶ W (was 6.5).

**The far clusters come back.** With the constants measured consistently, and every star mass on
the calibrating clusters' basis (z-scores; within ±2 passes):

| | rev 19 | rev 20 |
|---|---|---|
| MACS J0025: lensing inside 300 kpc, SE / NW; galaxy speeds | −0.53 / −0.87; 707 km/s (−2.2) | −0.65 / −0.97; 665 km/s (−2.9) |
| Abell 520: clumps P4 and P6; inside 710 kpc; galaxy speeds (rms) | −2.1, −2.6; −2.2; 3.2 | **−1.0**, −2.3; **−1.0**; 2.7 |
| El Gordo: lensing inside 0.5 / 1 / 1.5 Mpc; NW / SE galaxy speeds | −2.1 / −2.1 / −1.2; 917 / 822 km/s | **−1.1 / −1.2** / −0.3; **1,002 / 892** km/s |

* **Of the 17 graded checks of these three clusters, 12 pass, 4 are close and 1 fails** (rev 19:
  9, 7, 1).
* **El Gordo's lensing now agrees with its published star masses.** Rev 19's "1.4 times the stars"
  came mostly from mixing conventions: the far clusters were in our distances, while the companion's
  speed had been measured in the standard ones, on overcounted stars.
* **"Older stars" is withdrawn.** Rev 19 suggested that the far clusters' stars might be older, and
  so heavier, than the Big Bang's timeline allows. Observations say the opposite, without using any
  distance law: at fixed speed spread, cluster ellipticals at z = 0.83 are 0.24 magnitudes bluer
  than in the nearby Coma cluster, and a standard age indicator of their starlight (the 4000 Å
  break) weakens steadily with distance. Distant cluster stars are younger and lighter for their
  light. With u measured consistently, they don't need to be heavier.

**MACS J0025, the exception.** Its star masses turn out to be on a heavier basis (Salpeter's)
than the calibrating clusters' (Chabrier's), so on the same basis they are 1.78 times lighter than
published. Two things follow:
* its galaxy speeds come out 665 km/s against 835 ± 59; with the published masses they would be 804;
* its north-western lensing peak moves from the galaxies onto the gas. The galaxies there are only
  231 kpc from the gas, and the lensing map between them is a nearly flat ridge, highest at the gas
  by 9%. It depends on the collision's age and the star masses (distance of the NW peak from its
  galaxies, kpc):

  | Time since the clusters crossed | stars on the calibrating basis | × 1.33 | as published (× 1.78) |
  |---|---|---|---|
  | 0.26 billion years | **73** | 60 | 55 |
  | 0.35 billion years | **71** | 59 | 54 |
  | 0.5 billion years (the suite) | **209, on the gas** | 58 | 54 |

  Bradač et al. put the crossing "a few 10⁸ years" ago. A younger collision, or a third more star
  mass, puts the peak back on the galaxies. The suite keeps the 0.5 billion years fixed in rev 17,
  so this check now fails. It is a knife edge, and a clean test of the collision's age. *(Rev 21: the
  collision's own shock fronts date it at 0.1–0.36 billion years, and the peak passes; §6.16.)*

**The Bullet Cluster in our own distances** (sizes ×1.14, gas ×1.08, stars ×0.78, lensing ×1.18):

| | rev 19: ours | measured (standard distances) | rev 20: ours | measured, in our distances |
|---|---|---|---|---|
| lensing strength on the main / smaller half's galaxies | 0.68 / 0.14 | at least 0.36 / 0.20 | 0.72 / **0.26** | at least 0.36 / 0.20 |
| main half: stars in its outskirts (10¹² suns) | 6.2 | 3.9–6.7 | 3.3 | 3.0–5.2 |
| main half: lensing mass inside 250 kpc (286 in ours), 10¹⁴ suns | 2.40 | 2.5–2.8 | 2.83 | 2.94–3.29 |
| smaller half: the same | 0.94 | 2.0–2.3 | 1.29 | 2.35–2.70 |

The pattern holds: both peaks sit on the galaxies, 14 and 21 kpc from them. The smaller half gains
37% but is still about half its measured mass.

**The six strong lenses.** Light and matter still agree: −0.032 ± 0.023 dex. Their stars must be
1.37–1.85 times Salpeter's rule. *A correction to §6.14:* our distances make their stars 1.25 times
lighter for the same light, not 1.48; the lens code itself had it right.

**Galaxy lensing in our own distances: two new tests.** In our distances the KiDS points move: the
visible matter's pull ×0.64 and the measured lensing pull ×0.93 at the lenses' typical distance.
Median of log(measured ÷ predicted), zero being a perfect match:

| | all | spirals (blue) | ellipticals (red) | disks | bulges | extra lensing of ellipticals over spirals (measured 0.15 ± 0.04 this way) |
|---|---|---|---|---|---|---|
| standard distances (rev 19) | +0.02 | −0.005 | +0.02 | +0.03 | −0.02 | 0.18 |
| **our distances (rev 20)** | **+0.06** | **+0.08** | +0.04 | **+0.11** | −0.01 | **0.23** |

Two separate things have moved:
1. **The level, set by our distance law's scale.** The lensing pull inferred from the same data
   scales with α, the rate at which light loses energy in our distance law, while the visible
   matter's pull does not depend on it. At α 10% lower the spirals read +0.03 and the whole sample
   +0.02; the whole sample centres at 13% lower, an "H0-like" 65 km/s per Mpc instead of 74.6. So
   **galaxy lensing now measures our distance law.** α also sets SPARC's distances for about half of
   its galaxies, so the fair test refits them together; that is next. *(Rev 21: done, and this
   reading was wrong. Refitted together, the rate barely moves the lensing level, which measures the
   lenses' mass instead; §6.16.)*
2. **The ellipticals' extra lensing, set by the companion's speed.** An elliptical's stars carry
   heat 3σ²/u². With u = 163 km/s it is 1.5 times what it was at 197, and the ellipticals' extra
   lensing rises from 0.18 to 0.23 dex against 0.15 measured. *(Rev 21: the real cause was two assumed
   star speeds; measured, the extra lensing matches, §6.16.)* Clusters now want a slower companion
   than single ellipticals do. The first thing to check is the ellipticals' star masses: spectra
   say giant ellipticals make stars heavy for their light (about 1.85 times Chabrier's rule), and
   that would change the clusters' calibration and KiDS's ellipticals together.

*The extra lensing is the median over the same bins; Brouwer et al.'s mean differences, 0.17 and
0.27 dex (§6.4), are a different average of the same data.* Mistele et al.'s lensing speeds tell
the same story: spirals 17–40% above ours, ellipticals within 9%.

**Two ways to turn angles into sizes.** The project's distance law allows two geometries. We ran
every test in both. The galaxy lenses lean slightly (0.01 dex) to the second; the clusters clearly
prefer the one we use (El Gordo's lensing −1.1 against −2.4 error bars; the Bullet's smaller half
0.26 against 0.12). We keep it.

**The suite, with the adopted constants: 57 pass, 9 close, 10 fail** (rev 19: 58, 11, 7).
* Improved: Abell 520's clump P4, its mass inside 710 kpc and its galaxy speeds; El Gordo's lensing
  inside 0.5 and 1 Mpc.
* Worse: galaxy lensing around spirals and disks (fail), all lenses and the ellipticals' extra
  lensing (close), Mistele's spirals (fail), MACS J0025's north-western peak (fail).

The balance is honest rather than flattering: what got worse is exactly where our own distances
now make a sharper test, of the distance law's scale and of the companion's speed.

Scripts: `code/xcop_static_v11.py`, `code/kids_static_v11.py`, `code/bullet_static_v11.py`,
`code/distance_variants_v11.py`, `code/macs_peak_scan_v11.py`; the audit of every star mass is
`literature/star_mass_audit_v11.md`.

### 6.16 Round 12, step by step (rev 21)

Rev 20 ended with four next steps. This section adds each result as it lands.

**Step 1: MACS J0025's age, read from its own shock fronts.** When two clusters collide, the collision
launches shock waves through their gas. Radio telescopes see where those shocks have lit up particles:
"radio relics". Riseley et al. (2017) found two relics in MACS J0025, one on each side, lying across
the direction of the collision.
* We measured where their radio contours are centred: 150–160 kpc north-west of the centre and
  315–325 kpc south-east of it, in the paper's units. The two groups of galaxies sit 165 and 360 kpc
  out. **The shocks have not yet outrun the galaxies, which only happens early in a collision.**
* A shock moves at most about 1.9 times the speed of sound in the gas (the relic's radio colour limits
  it), and that speed is about 1,300 km/s. So the shocks took at least 0.06–0.16 billion years to get
  where they are.
* The simpler clock, separation over speed, gives 0.26 billion years: 540 kpc between the galaxy groups
  at about 2,000 km/s.
* In our distances sizes are 1.36 times larger, so the clocks read **0.08–0.36 billion years**. The
  0.5 we had used came from averaging with a clock based on the galaxies' star formation (0.5–1 billion
  years), which dates when star formation stopped. That can begin before the clusters reach each other.

| Age used | Fresh companion around the gas | North-western peak from its galaxies |
|---|---|---|
| **0.3 billion years (now used)** | 50 kpc | **72 kpc: on the galaxies** |
| 0.4 billion years | 67 kpc | 70 kpc: on the galaxies |
| 0.45 billion years | 75 kpc | 207 kpc: on the gas |
| 0.5 billion years (rev 17–20) | 83 kpc | 209 kpc: on the gas |

**The peak stays on the galaxies for any age up to 0.4 billion years**, beyond every dynamical clock, and
nothing else about MACS J0025 depends on the age. The suite now uses 0.3 billion years: **58 pass,
9 close, 9 fail.** Script: `code/macs_timing_v12.py`.

**Step 2: how hot are the galaxies that lens? Measured, not assumed.** In our law a galaxy's extra
lensing depends on how fast its stars move at random. For galaxy lensing (KiDS) we had used two round
numbers: 160 km/s for every star of an elliptical-like (red) lens, and a cold disk, with no heat at all,
for every spiral-like (blue) lens. We measured instead, from SDSS spectra of 119,000 galaxies with the
same stellar masses, split the same way the KiDS team split theirs.
* **Red galaxies** of the lenses' typical mass have stars moving at 140–160 km/s, a little cooler than
  assumed.
* **Blue galaxies are not cold.** At these masses they have central bulges whose stars move at
  100–120 km/s and hold 10–35% of their light.

| Heat weight (0 = cold disk) | red lenses | blue lenses | bulge-dominated | disk-dominated |
|---|---|---|---|---|
| assumed until rev 20 | 2.9 | 0.1 | 2.9 | 0.1 |
| **measured** | **2.1** | **0.55** | **2.2** | **0.30** |

What that does (median log of measured ÷ predicted lensing, zero is perfect; our distances):

| | all | spirals | ellipticals | disks | bulges | ellipticals over spirals (measured 0.15) | the same by shape (measured 0.15) |
|---|---|---|---|---|---|---|---|
| assumed heat | +0.06 | +0.08 | +0.04 | +0.11 | −0.01 | 0.23 | 0.23 |
| **measured heat** | +0.07 | **+0.02** | +0.08 | +0.08 | +0.04 | **0.13** | **0.16** |

* **The difference between ellipticals and spirals now matches, both ways of splitting them.** It was
  never a sign that the companion's speed was wrong. Two assumed numbers caused it.
* Lensing speeds around spirals (Mistele et al.) improve from 4.4 to 3.0 error bars. Those around
  ellipticals move from 1.0 to 2.5, as they join the spirals at a common level.
* **What is left is one common level:** all the lens samples sit 0.02–0.08 dex above our law. Step 3
  asks whether our distance law's scale sets it.
* **Heavy stars (the IMF), checked:** spectra show giant ellipticals' stars are heavy for their light
  mainly in their centres. Averaged over whole galaxies that adds at most 0.09 dex, and only above
  2 × 10¹¹ suns (Domínguez Sánchez et al. 2019). That is too small to matter for the lenses or for
  the clusters that calibrate u. It trims the strong lenses' extra star need from 0.45 to about 0.36
  dex.

The suite now uses the measured heat and also grades the split by shape: **59 pass, 10 close, 8 fail**
(of 77). Scripts: `code/lens_heat_sdss_v12.py`, `code/kids_heat_v12.py`.

**Step 3: our distance law's one number, fitted to everything at once.** Our law turns a galaxy's
redshift into its distance with one number: the rate at which light loses energy on its way. So far
it came from 164 nearby groups of galaxies, and corresponded to a Hubble-constant-like 74.6 km/s per
Mpc. Rev 20 suggested that galaxy lensing wanted it 13% lower. We fitted it to four kinds of data at
once, refitting the law's three constants at every trial value.
* **Supernovae.** 1,365 exploding stars, their brightness calibrated on Cepheid stars in nearby
  galaxies, with the Pantheon+ team's full error bookkeeping. They prefer a rate 4.5% lower (± 1.3%).
* **The galaxies' rotation.** 81 of our 149 galaxies have distances that come from their redshifts;
  the other 68 have distances measured directly. Put in our own distances, the 81 agree best with the
  68 at a rate 13% lower (8–18%), or 6% lower by a second way of scoring.
* **Clusters and galaxy lensing hardly care.** When the rate moves, the law's constants move with it
  and cancel most of the change.

| Data | The rate, compared with before | Hubble-constant-like value |
|---|---|---|
| 1,365 supernovae | 4.5% lower (± 1.3%) | 71.2 |
| 81 galaxies placed by their redshifts | 13% lower (8–18%) | 65 |
| **together (adopted)** | **5% lower (± 1.3%)** | **70.9** |

**We adopt it.** All our data now sit in the same distances, including those 81 galaxies, the last
set that did not. The constants move a little: a from 6.55 to 6.30 × 10⁻¹¹ m/s², g_d from 2.11 to
2.03 × 10⁻¹⁰ m/s², and the companion's speed from 163 to 169 km/s. The galaxies fit slightly better
(a typical miss of 15.87 km/s against 15.94).

**A correction.** Rev 20 said galaxy lensing measures this rate. It does not. That held only with the
constants kept fixed. Refitted, the lensing level moves by just 0.02 dex across a 25% range of the rate.

**What the lensing level does measure.** We changed one thing at a time (median log of measured ÷
predicted, all lenses together; zero is perfect):

| Change | All lenses |
|---|---|
| none (as adopted) | +0.065 |
| the other way of turning angles into sizes | +0.058 |
| lenses 5% farther away, judged by their brightness | +0.042 |
| lenses 5% larger on the sky | +0.066 |
| gas around each lens weighing half as much as its stars | +0.026 |
| **gas around each lens weighing as much as its stars** | **−0.002** |

* The level follows how much mass the lenses have for their light. Geometry hardly matters. The
  supernovae now pin the lenses' distances to about 1.4%, so what is left is mass that the count of
  visible matter leaves out.
* KiDS counts stars and cold gas. Galaxies are also wrapped in thin, hot gas that is hard to weigh.
  The KiDS team's own middle estimate is as much again as the stars, within 100 kpc. **With that much
  gas, our law matches all the lenses together.** In our law every bit of ordinary matter pulls, so
  this is a prediction: X-ray and ultraviolet measurements of the gas around such galaxies can check
  it. X-ray haloes are seen mostly around ellipticals, which is where our law needs the most gas (red
  lenses +0.078, blue +0.011).

**The suite, with the rate fitted: 59 pass, 11 close, 7 fail** (step 2: 59, 10, 8). No check got a
worse grade. Lensing around disks moved from fail to close, and 14 other checks moved closer to their
measurements: the galaxies, most of the KiDS samples, the strong lenses, and the lensing masses of El
Gordo and Abell 520 among them. Seven moved slightly away without changing grade: six Milky Way
numbers (the Sun's speed 210.7 → 209.2 km/s, against 229–234) and the colour gap (0.132 → 0.126,
against 0.153). The locked wide-binary forecast moves from 1.083 to 1.076 at 20,000 AU (§6.7).

Scripts: `code/sn_scale_v12.py`, `code/distance_scale_v12.py`, `code/kids_level_v12.py`.

**Step 4: a working model of the companion.** Rev 20 showed that energy conservation fixes how
orderly matter's companion adds up, provided its flow has no whirlpools (§4.15). That last condition
was an assumption. We built a small computer model: a flat sheet on which matter gives off the
companion's energy, and five candidate rules for how that energy moves:
1. **fly straight**, in every direction, passing through each other;
2. **bounce off each other at random**, the way heat spreads;
3. **destroy each other** when they meet head-on;
4. **turn to follow their neighbours**, like a flock of birds;
5. **follow the lines of ordinary gravity** outward, the way certain waves in the Sun's hot gas run
   along its magnetic field lines.

We graded each on the three requirements, with the matter laid out four ways: a disk, two equal
clumps, an unequal pair and three clumps.

| Rule | Energy kept | One stream (1 = fully) | Off Newton's pattern | Whirlpools |
|---|---|---|---|---|
| fly straight | all | 0.83–0.87 | 0.01–0.02 dex | none |
| bounce at random | all | 0.12–0.13 | 0.06 dex | none |
| destroy head-on | 76–81% | 0.83–0.87 | 0.10–0.13 dex | none |
| follow neighbours | all | 0.98 | 0.35–0.75 dex | **they form on their own** |
| **follow gravity's lines** | **all** | **1.00** | **0.01–0.02 dex** | **none** |

* **Flying straight** keeps every watt and follows Newton's pattern, but where streams from different
  places cross, the energy is not one stream. That weakens the pull, which the galaxies rule out.
* **Bouncing at random** keeps every watt, but its energy barely flows at all (one stream: 0.12).
* **Head-on destruction** throws energy away and still leaves crossing streams.
* **Following neighbours** makes one stream, but whirlpools appear by themselves, as they do in flocks
  of birds and schools of fish, and the flow strays far from Newton's pattern.
* **Following gravity's lines** meets all three. Its small misses are the grid's: they halve each time
  the grid is made twice as fine (0.038, 0.020,
  0.011 dex), while the neighbours' whirlpools stay the same
  (0.59 and 0.58 on our whirlpool scale, where Newton's own pattern scores 0.001).

**Why following gravity's lines works, exactly.** Picture a bundle of gravity's field lines, like a
bundle of drinking straws, starting at the centre of a clump of matter. As the bundle passes through
matter, Newton's pull through it grows in proportion to the matter it passes. If the companion flows
along the same lines, its flow through the bundle grows in proportion to the same matter, and none of
it leaks out of the sides. Both start at zero at the same point. So they stay in the same proportion
everywhere: the companion's flow follows Newton's pattern exactly, with no whirlpools, and nothing had
to be assumed. The model shows it; the argument proves it.

**What this means.** The last assumption in rev 20's derivation becomes a physical property: **the
companion is a wave guided by gravity's own field lines.** Nature has waves like this. Alfvén waves in
the Sun's corona carry their energy along magnetic field lines, whichever way their crests face. Two
more things fit:
* random motion scrambles a wave's direction, so hot matter's companion is not guided and adds up
  without cancelling, which is what the clusters need (§4.15);
* an Alfvén wave is carried along with the gas it lives in, and a guided companion carried along with
  its matter would have exactly the memory the colliding clusters need (§3.10). That is the next
  model: sources that move.

Script: `code/companion_toy_v12.py`.

### 6.17 Round 13: the pull and the scrambling in one experiment (rev 22)

**What you asked.** Two pieces of the theory worked, but only separately:
* **the pull** (§4.14): matter that feeds a passing wave a quarter beat ahead of it is pulled toward the wave's source;
* **the scrambling** (§3.5–3.6): random motion knocks the companion out of step, and frequent collisions keep it in
  step.

You asked for one experiment that joins them, with the same matter and the same rules, changing only how the matter
moves: orderly, randomly and freely, or randomly with frequent collisions. Small test bodies would feel only the wave
where they sit, with nothing of our law built in. Success would be:
* free random motion pulls *harder*;
* collisions take the extra pull away;
* doubling the speed spread roughly doubles the extra pull, without that being programmed in.

This was proposal 3 of §9.1. Midway, when the first answer came out backwards, you asked us to look for another
solution. We found one, and it works in part.

**The set-up.** Picture 100 tiny singers in a ball about two companion wavelengths across. Each one:
* sings the companion wave;
* listens to the wave around it;
* keeps adjusting its timing to stay a fixed amount ahead of, or behind, what it hears.

Thirty-two more singers of the same kind sit six wavelengths away as test bodies. We measure:
* how hard, on average over long times, the test bodies are pulled toward the ball;
* how loudly the ball sings compared with 100 singers who ignore one another.

Speeds are measured against how fast a singer can re-time itself. At speed 1, a moving singer drifts out of step
exactly as fast as it can correct. That also turns the speeds into our law's heat weight: speed 1 is k ≈ 1, speed 3
is k ≈ 9 (a giant elliptical's stars), and speed 10 is k ≈ 100 (a rich cluster's galaxies). The model passed its
check: one singer and one test body give exactly the pull that §4.14 derived, to five digits.

**First try, with the rule from rev 19.** Every singer runs a quarter beat ahead of the wave it hears.
* **At rest** the singers fall into step and sing 4.4 times louder than independent singers would. The test bodies
  are pulled twice as hard as by independent singers.
* **Moving freely at random,** the chorus breaks up, and the pull drops to about half (0.36–0.49 of the pull at rest).
* **Colliding,** the chorus survives and so does the pull (0.87–0.96). Collisions protect exactly as §3.6 says, and by
  the amount the Dicke formula behind §3.6 predicts.
* **Rotating in an orderly way,** the chorus survives, but its pattern sweeps past the test bodies, which can't follow
  it.

So the scrambling half works, but heat weakens the pull instead of strengthening it.

**Why, exactly.** Two facts decide it, and neither can be tuned away.
1. **A body is pulled by the recoil of what it adds to the passing wave:** the pull equals the power it feeds the wave,
   divided by the wave's speed. It is like a rowing boat pushing water backwards: no feeding, no pull.
2. **When singers are in step, every one of them, test bodies included, sits at the same timing relative to the wave
   it hears.** So either all of them feed or all of them absorb. If all feed, the ball sings louder than independent
   singers: a chorus. If all absorb, it sings more quietly: a hush. A hushed ball's test bodies absorb too, so they are
   pushed.

We checked this at 16 different timings. At every one, the test bodies were pulled exactly when the ball sang louder
than independent singers. With one fixed timing, heat can only break a chorus up, so it can only weaken the pull.

**What the opposite timing showed.** With singers a quarter beat *behind* the wave:
* the ball at rest is hushed, at a fifth of independent singers' volume (a denser ball: a 35th);
* free random motion breaks the hush, making the ball up to four times louder;
* collisions keep the hush.

That is exactly the pattern our heat term needs, and the extra power comes out of the singers' own supply. It has to:
if it came out of their motion, a cluster's galaxies would lose their speed in 85 million years, and clusters have
stayed hot for billions. But the test bodies of a hushed ball are pushed. So we needed one rule that hushes a ball
from inside and still makes distant bodies feed.

**Another solution: timing that depends on how loud the wave is.** The rule: every singer runs a quarter beat ahead
of a *quiet* wave, feeding it, and falls a quarter beat behind a *loud* one, absorbing it. Everyone follows the same
rule.
* Inside a dense ball, the singers' waves overlap and are loud, so the ball hushes itself.
* Far away the wave is quiet, so the test bodies feed it and are pulled.

It is how an amplifier with a small built-in loss behaves: a faint signal gets boosted, but a loud one uses up the
gain, and the loss soaks it up.

**It works, in part: the first time the three pieces behave as the law needs, in one simulation.** In our densest
ball, 87 times quieter than independent singers when cold, the pull on the test bodies compared with the same ball
at rest was (four runs each; speeds as fractions of the re-timing speed):

| speed | 1/32 | 1/16 | 1/8 | 1/4 | 1/2 | 1 |
|---|---|---|---|---|---|---|
| **moving freely at random** | ×1.32 | ×1.60 | ×1.83 | **×1.99** | ×1.91 | ×1.66 |
| with frequent collisions | ×1.07 | ×1.07 | ×1.18 | ×1.39 | ×1.91 | ×1.99 |
| rotating in step | ×1.01 | ×0.98 | ×0.96 | ×0.90 | ×0.81 | ×0.55 |

* **Free random motion strengthens the pull, up to twice.** The motion breaks the hush and the ball gets louder
  (seven times louder at a quarter of the re-timing speed), and the test bodies feed the louder wave.
* **Frequent collisions hold the gain back,** at the slower speeds to a small fraction, because they keep the hush.
  That is what gas needs. In a second, less dense ball, collisions removed the gain entirely (×0.90–0.95, against
  ×1.10–1.35 moving freely).
* **Rotating in step gives no gain.**
* **Our law's square root appears by itself.** The wave reaching the test bodies grows exactly as the square root
  of the cold output plus the released output. That is the form √(|g_N| + S) our law uses (§3.8), and nobody put it
  in.
* **Doubling:** going from 1/32 to 1/16 of the re-timing speed multiplies the released output by 3.3 and the extra
  pull by 1.9. The law has 4 and 2. After that the gains shrink.
* **Where it stops:** once the motion is faster than about a quarter of the re-timing speed, the test bodies can't
  keep in step with the flickering wave, and their pull falls back even though the ball keeps getting louder.
* **Why that limit may matter less than it looks:** with a hush this deep, heat already dominates at 1/32 of the
  re-timing speed. So the conversion from speed to heat weight given above is too pessimistic. The deeper the cold
  hush, the more of the heat range lies at speeds where bodies keep in step.

**What the data say about "doubling".** We asked the real data how steeply the heat weight must grow with the speed
spread σ, instead of assuming σ². The clusters can't tell: their galaxies all move at similar speeds, and the
companion's speed adjusts to compensate. Galaxy lensing can, because its lens galaxies' stars move about six times more
slowly. The difference in lensing between ellipticals and spirals (0.15 dex measured) is matched only if **doubling σ
multiplies the extra pull by 1.8–2.0**. Our law has 2.0. Whatever the mechanism, it must produce that.

| heat weight grows as | σ¹ | σ^1.5 | σ^1.75 | **σ² (the law)** | σ^2.5 | σ³ |
|---|---|---|---|---|---|---|
| companion's speed the clusters then want | 33 km/s | 98 | 134 | **169** | 237 | 296 |
| clusters' typical miss | 25% | 25% | 25% | **25%** | 25% | 25% |
| ellipticals' extra lensing (colour split; 0.153 ± 0.04 measured) | 0.149 | 0.160 | 0.148 | **0.126** | 0.077 | 0.039 |
| the same, by galaxy shape (0.154 measured) | 0.201 | 0.210 | 0.188 | **0.155** | 0.089 | 0.043 |
| all lenses' common excess over the law | −0.22 | −0.05 | **+0.01** | +0.065 | +0.14 | +0.18 |

One side lead: at σ^1.75 the lenses' common 16% excess (§6.16) almost disappears (+0.01 dex), while the clusters fit
just as well. That is worth a full regression-suite run.

**What it means.**
* **The pull and the scrambling now live in one experiment.** Both behave as derived: the pull from feeding a passing
  wave, and collisions protecting the timing.
* **The simple story for heat, "scrambled contributions don't cancel, so they pull harder" (§3.5, §4.6), is ruled
  out:** scrambling weakens a chorus. The heat term itself still stands on the data (clusters, the elliptical–spiral
  lensing gap, colliding clusters), and the data demand its σ² scaling.
* **Its likely mechanism is now in view.** Cold matter holds its companion back in a hush, motion releases it, and
  collisions keep it held. The energy comes from matter's own supply. The one rule we found, "feed quiet waves, absorb
  loud ones", gives the right signs in the same experiment.
* **Still to do:**
  * find the physics behind that rule;
  * let bodies follow a wave that changes quickly (our test bodies lose the extra pull once the motion outpaces their
    re-timing);
  * make the release grow as steeply as the data say.
  * Also: whether the same rule is what our law calls "strong gravity holds the companion back" (§3.4). That hold
    works the same way, strong waves held back and weak ones free, so it is the first thing to check.

Scripts: `code/coherent_force_v13.py`, `code/joint_checks_v13.py`, `code/strength_offset_v13.py`,
`code/heat_exponent_v13.py`.

### 6.18 Round 14, step by step (rev 23)

Rev 22 ended with a list of next steps (§9, item 9). This section adds each result as it lands.

**Step 1: a gentler heat term, tested on everything.** Round 13 found that galaxy lensing allows the heat term to
grow a little more gently with the stars' speed spread, as σ^1.75 instead of σ², and that the gentler version would
erase the 16% by which all lenses sit above our law. We ran the whole regression suite that way.
* To do it fairly, we first made the exponent a single switch: every place the code computes the heat weight now
  goes through one function. With the switch at 2, the suite reproduces its saved results to the last digit.
* We also found and fixed a slip in the suite. When asked to refit the companion's speed, it used an old version of
  the cluster data (stars counted in projection, standard distances) and got 192 km/s instead of 169. Fixed, the
  refit reproduces the adopted law exactly.

The result, with the two constants refitted to the gentler term (the companion's speed becomes 132 km/s):

| | now (σ²) | gentler (σ^1.75) |
|---|---|---|
| overall | 59 pass, 11 close, 7 fail | 60 pass, 10 close, 7 fail |
| all lenses' common level (KiDS) | 16% above (close) | **3% (pass)** |
| the red lenses, failing since rev 20 | 20% above (fail) | **4% (pass)** |
| massive ellipticals: SLACS lensing against star speeds | −0.028 dex (pass) | −0.054 (close) |
| ellipticals' lensing speeds (Mistele et al.) | 2.6 (close) | 5.2 (fail) |
| bulge-dominated galaxies' rotation speeds | 29.4 km/s (pass) | 30.5 (close; MOND 30.35) |

It trades one tension for another. The ordinary ellipticals in Brouwer et al.'s analysis of the KiDS survey want a
little more heat; the massive ones (SLACS), and Mistele et al.'s analysis of the same survey, want a little less. We
keep σ² for now. The switch is ready for when the two KiDS analyses are reconciled. Scripts:
`regression/candidates/heat_p175.json`, `code/law.py`.

**Step 2: where "feed quiet waves, absorb loud ones" comes from.** Round 13's rule worked, but we had put it in by
hand. We tried two physical versions of a piece of matter and dropped each into the same experiment.

*A singer that also absorbs.* Every piece of matter keeps singing as before (a quarter beat ahead of the wave
around it), and also soaks up a little of any wave passing through it, as most materials do. A lone body then
gains from a quiet wave in proportion to its strength, but loses to a loud one in proportion to the strength
squared. So it feeds quiet waves and absorbs loud ones, without any switch.
* **A cold, dense ball goes dark and still pulls:** 47 times quieter than independent singers, with its test bodies
  pulled. The absorbing part swallows the ball's own chorus.
* **But heat doesn't release it.** When we made the absorbing part tune itself like a real resonance, random motion
  released nothing (×1.04–1.07 up to the re-timing speed). An absorber keeps absorbing whatever arrives, so motion
  can't free what it holds. A simpler version did show a big gain, up to four times the pull, but that came from
  motion shaking a jammed ball into step. It also had a memory: stopped again, two of four balls stayed bright.

*A singer with a power balance.* Give each piece of matter a fixed supply of power (our law's 5.3 × 10⁻⁶ watts per
kilogram) and a loss that grows with the loudness of the wave around it (the same absorption as above). Whatever it
doesn't lose, it must feed into the wave. That fixes its timing, the way the power balance fixes whether a machine
on the electrical grid runs as a generator or a motor:
* in a **quiet** wave its supply beats its loss: it runs ahead and feeds the wave, and is pulled;
* in a **loud** wave its loss beats its supply: it falls behind and draws from the wave, and is pushed.

That is round 13's rule, derived from two physical numbers. In our densest ball (four runs each):

| speed | 1/32 | 1/16 | 1/8 | 1/4 | 1/2 | 1 |
|---|---|---|---|---|---|---|
| **moving freely at random** | ×1.30 | ×1.51 | **×1.81** | ×1.71 | ×1.65 | ×1.69 |
| with frequent collisions | ×0.98 | ×1.01 | ×1.08 | ×1.31 | ×1.68 | ×1.75 |
| rotating in step | ×0.98 | ×0.95 | ×0.93 | ×0.88 | ×0.80 | ×0.52 |

(pull on the test bodies compared with the same ball at rest)
* The heat pattern comes back: free random motion strengthens the pull, frequent collisions hold the gain back, and
  rotation gives none.
* **No memory:** stopped again, all four balls went back to their dark, cold state. The gain is a true response to
  the motion.
* Beyond about an eighth of the re-timing speed the test bodies lose step again, as in round 13 (step 3 looks at
  this).

**And the release factor may be the same thing.** Our law holds the companion back where ordinary gravity is strong,
using an exponential we fitted, not derived (§3.4). Both singers above hold back in loud waves too, with shapes of
their own. We let the 149 galaxies judge (typical miss in rotation speed, constants refitted):

| how the companion is held back | typical miss |
|---|---|
| the law's exponential | 15.87 km/s |
| the singer that absorbs | **15.86 km/s** |
| the singer with a power balance | 16.08 km/s (and the best fit by the measure the fits use) |
| variants that push in strong gravity, or let the pull come back there | 16.19–16.78 km/s |

The galaxies can't tell our exponential from the absorbers' hold. So "strong gravity holds the companion back" and
"loud waves are absorbed" may be one effect: the same loss that sets each singer's timing.

Scripts: `code/singer_v14.py`, `code/release_shape_v14.py`.

**Step 3: how far the release goes, and whether bodies keep up.** We ran the power-balance singers in three balls of
increasing density, 33, 76 and 107 times quieter than independent singers when cold, at speeds from 1/64 to twice the
re-timing speed (three runs each).

| free random motion | 1/32 | 1/16 | 1/8 | 1/4 | 1/2 | 1 | 2 |
|---|---|---|---|---|---|---|---|
| released ÷ held, densest ball | 0.8 | 3.3 | 5.1 | 7.7 | 9.3 | 14 | 24 |
| pull ÷ at rest, densest ball | ×1.52 | **×1.95** | ×1.80 | ×1.68 | ×1.68 | ×1.39 | ×0.76 |
| pull ÷ at rest, middle ball | ×1.22 | ×1.43 | ×1.71 | ×1.68 | ×1.63 | **×1.73** | ×0.82 |

* **The darker the cold ball, the earlier heat takes over.** The released output passes the held output at 1/8, 1/16
  and 1/24 of the re-timing speed: about four times the ball's cold darkness.
* **Bodies keep up further than we feared.** Their pull stays raised up to about the full re-timing speed, because
  the released wave grows faster than their step slips; only at twice that speed do they lose it. So the range of
  speeds where heat dominates and bodies keep step widens as matter gets darker: 8-fold, 16-fold, 24-fold. Real
  matter has vastly more singers per wavelength than our 100, so it can be far darker still.
* **The square root holds everywhere:** the wave reaching the bodies is always the square root of cold plus released
  output.
* **What doesn't match yet is the steepness.** After a steep start, the released output grows only about 1.6 times
  per doubling of speed. Our law, and the lensing data, want about 4 times. The hold empties too gently. Finding
  what makes it empty as σ² (a bigger ball, a loss that grows faster with loudness, or matter of several kinds) is
  the next task.

Script: `code/release_depth_v14.py`.

### 6.19 Round 15: an independent idea for where the heat term comes from, joined to the pull (rev 24)

You sent us an independent calculation with a new idea for the heat term and asked us to build toward a solid proof
from first principles. Here is the idea, and what happened when we joined it to our pull.

**The idea: matter has a quiet store, and motion opens it.**
* Picture each piece of matter holding a **quiet internal vibration**. It carries the piece's internal energy, but
  its parts swing against each other, so almost nothing leaks out.
* It also has three **loud vibrations** that radiate freely, one for each direction.
* When the piece's parts move relative to each other, their rhythms drift slightly apart. That mixes a little of the
  quiet vibration into the loud ones.
* The mixing grows in proportion to the speed, and the leaked energy with its square, σ². That is exactly the form
  of our law's heat term.
* Collisions keep changing the direction of motion, so the mixing never builds up.

We reran the calculation here and every number matches; four of its five result files are identical to the last
byte:
* each doubling of speed multiplies the extra glow by 4.00, 3.99, 3.98 and 3.92, until the store starts to run low;
* collisions leave 50.7%, 9.2% and 4.8% of it, where the simple formula says 50%, 9.1% and 4.8%.

Two things were missing, and the calculation's own README says so. It did not show a pull, and its pieces were not
yet sending real waves to anything. That was our job.

**Step 1: does the wave alone give this?** We built a small cloud of 40 ordinary vibrating particles. They talk to
each other only through the wave, and we assumed nothing about quiet stores or drifting rhythms.
* **Spinning the cloud, or moving it as a whole, releases exactly nothing** (less than one part in a trillion). The
  wave cares only about the distances between particles, which is just what the idea assumed.
* **Random motion releases energy, and collisions hold it back.** With 1, 10 and 100 direction changes per unit
  time, 38%, 7% and 0.4% of the extra glow survives. That is our law's "gas is cold" rule, from the wave itself.
* **But the glow grows only 1.3–1.8 times per doubling of speed, not 4 times.** That is the same shortfall as in round
  14.
  * The reason: a cloud's quietness comes in every shade. Its patterns of vibration leak at rates spread across ten
    powers of ten, with no clean line between quiet and loud.
  * The σ² needs one quiet store that is sharply separated from the loud channels, a million times quieter in the
    calculation.
* **That is a useful answer.** The quiet store has to be something inside each piece of matter, not a result of how
  the pieces are arranged. It also explains why round 14 fell short: its quietness came from the arrangement.

**Step 2: which bodies does a passing wave pull?** We put six kinds of body, one at a time, in a steady wave from a
distant source, and let each respond by its own equations:

| body | falls into step | pulled or pushed | pull grows with the wave's |
|---|---|---|---|
| an ordinary oscillator, just driven by the wave | — | pushed | energy |
| an "inverted" one, below the point where it would ring by itself (like a laser amplifier) | — | pulled | energy |
| our round-10 rule, put in by hand | a quarter beat ahead | pulled | height |
| **the calculation's own oscillator** (a self-sustained ringer with a power supply) | **a quarter beat behind** | **pushed** | height |
| **an inverted, self-sustained ringer** | **a quarter beat ahead, by itself** | **pulled** | **height** |

* An ordinary self-sustained ringer, like a clock or an ordinary laser, falls into step a quarter beat *behind* a
  passing wave. It soaks up energy and is pushed away. The calculation's oscillator does this.
* An **inverted** ringer holds energy it is ready to give, like the excited atoms in a laser. It falls into step a
  quarter beat *ahead*, feeds the wave and is pulled, in proportion to the wave's height. That is exactly our
  round-10 rule, which until now we had to put in by hand.
* It comes out of standard equations for such emitters, the ones used for a real laboratory device, the
  "superradiant laser". Matter holds an enormous store, its rest energy, so treating it as inverted toward the
  companion is natural.
* In very strong waves its pull stops growing, because it can only feed as fast as it is refilled. That is a
  ceiling, not the switch-off our law uses near the Sun, so our release factor stays a separate question.

**Step 3: the whole chain in one experiment.**
* 100 pieces of matter, each with its quiet store and its three loud vibrations, sit in a ball one wavelength across.
  Their waves go out for real.
* 64 test bodies sit six wavelengths away. Each feels only the wave where it is. Three kinds share every wave: our
  round-10 bodies, the inverted ringers (whose timing nobody sets), and amplifiers.
* We measured force, timing, energy and momentum each on its own. Nowhere does the code turn power into force.

With k the heat weight (how much motion adds to the cold glow) and the pull compared with the same ball at rest:

| speed spread | k | glow released ÷ cold glow | pull, round-10 bodies | pull, inverted ringers | pull, amplifiers |
|---|---|---|---|---|---|
| 1 | 0.19 | ×1.19 | ×1.13 | ×1.15 | ×1.27 |
| 4 | 3 | ×4.1 | ×2.26 | ×2.41 | ×4.9 |
| 8 | 12 | ×13.2 | ×4.1 | ×4.4 | ×16 |
| 16 | 48 | ×48 | ×7.8 | ×8.2 | ×59 |
| 32 | 192 | ×162 (the fuel starts to run low) | ×14.3 | ×14.8 | ×202 |

* **The pull grows as the square root of cold plus released glow.** That is exactly our law's heat term,
  √(|g_N| + S), and no square root was put in anywhere.
* **Doubling the speed spread doubles the extra pull where heat dominates:** ×2.48, 2.20 and 1.97 on successive
  doublings. This is the decisive outcome we set ourselves in round 13, now met.
* **The inverted ringers do the same with no timing put in.** Amplifiers follow the glow itself, not its square root.
* **Collisions hold the heat back twice.**
  * At the source, the extra glow falls in proportion to 1/(1 + collision rate), to within 1%.
  * At the test bodies, the colliding pieces' waves flicker. Bodies that can't re-time fast enough lose step and gain
    almost nothing. Either way, colliding matter adds almost no heat, as our law says of gas.
* **No memory.** Stopped, the extra pull disappears; only the spent fuel is gone.
* **Momentum checks out.** Every body's pull equals the power it feeds the wave divided by the wave's speed, within
  1.1%, in every run.

**What being pulled costs.** A body pulled by feeding a wave pays for it, and the price depends on how fast the wave
travels.
* Our law has the companion streaming at u = 169 km/s. Being pulled at our law's a then costs a kilogram 1.1 × 10⁻⁵
  watts. That is exactly twice the power it gives off when cold, and about 4 parts in 10¹⁵ of its mass a year.
* If the companion travelled at light speed, the bill would be 3,500 times larger. The Sun would lose 5 parts in 10¹²
  of its mass each year to the Galaxy's companion. The planets' orbits would show that: their limit is about 1 part
  in 10¹³.
* **So energy and momentum require a slow companion,** which is what the clusters told us independently when they
  set u.

**What the constants say.** Two of the idea's numbers now have values or bounds:
* If a piece's quiet store is its rest energy, the store leaks at 3 × 10⁻²³ of itself per second, fixed by our law's
  5.3 × 10⁻⁶ watts per kilogram.
* Our law's two rules pin down how long the loud vibrations ring:
  * stars carry heat, so the loud vibrations must ring out faster than stars change direction, about every 200
    million years;
  * cluster gas doesn't, so they must ring longer than its ions take to circle, about 100 seconds.

  So each loud vibration rings for between about 100 seconds and 200 million years.
* If the drift in rhythm is an ordinary Doppler shift of the companion, its wavelength lies between about 60 parsecs
  and 500 million parsecs.

**Where the proof stands:**
* **Now derived, in explicit equations:**
  * moving matter releases companion power as σ², collisions hold it back, and spinning releases none;
  * the power travels as waves and adds up;
  * inverted matter falls into step ahead of the wave and is pulled by its height;
  * so the pull grows as √(cold + released), and doubling σ doubles the extra pull.
* **Still assumed:**
  * that each piece of matter has a quiet internal store with a clean gap to its loud channels. The wave alone
    gives the symmetries and the collision rule, but not the gap;
  * the size of the rhythm drift per unit speed, which with the ringing time sets u = 169 km/s.
* **Still separate:** the release factor that silences the companion near the Sun.

**Next:**
* find a physical quiet store with a gap (an internal pair of vibrations that relative motion pulls apart, as an
  electric field mixes a quiet atomic state into a bright one) and derive u from it;
* one kind of matter that is both source and receiver, in one cloud;
* whether a receiver drained by its neighbours' loud waves can switch its pull off, which would be the release factor.

Scripts: `code/wave_dark_v15.py`, `code/receivers_v15.py`, `code/reservoir_force_v15.py`; the independent
calculation and its reproduction are in `independent-r15/`.

## 7. How this compares

| | Ours | MOND | Dark matter |
|---|---|---|---|
| Rotation speeds, 149 galaxies | **15.9 km/s** | 16.1 km/s | 7.5 km/s (298 numbers) |
| Cluster masses, 12 clusters | **25%** | ×2.9 | 11% (24 numbers) |
| Bullet Cluster: lensing on galaxies | **yes** | no | yes |
| Bullet Cluster: main half | **lensing mass 2.97 against 3.09–3.46 (our distances); galaxy speed 1,249 km/s; star count agrees** | no | yes (fitted) |
| Bullet Cluster: smaller half's lensing mass | about half of the measured 2.47–2.85 (our distances; open) | no | yes (fitted) |
| Wide binary stars (data disputed) | **8% extra pull at 20,000 AU, 3% at 7,000 AU** | 43% | none |
| Collisions: lensing stays with galaxies | **yes** | no | yes |
| MACS J0025.4−1222 (a second Bullet) | **lensing masses and both peaks agree** at the age its shock fronts give; speeds 2.9σ low (our distances) | no | yes, fitted |
| Abell 520's galaxy-poor lensing clump | **from its gas and the galaxies' heat** | no | a puzzle |
| El Gordo (z = 0.87) | **lensing agrees with its published stars** (our distances); NW speeds 22% low | no | yes, fitted |
| Ellipticals lens more than spirals | **yes, from their stars**: 0.13–0.16 dex against 0.15 measured, with the lenses' star speeds measured (our distances) | no | yes, via tuned haloes |
| Strong lenses: light and stars agree | **yes** (−0.028 dex with our own distances) | | yes |
| Galaxy lensing (KiDS), spirals / ellipticals | **−0.005 / +0.02 dex** (standard distances); +0.01 / +0.08 in ours, with measured star speeds | −0.07 / +0.11 dex | yes, tuned |
| Milky Way rotation 15–27 kpc | **within 1–6%** | within 3% | 6–11% fast (McMillan) |
| Milky Way rotation at the Sun | 9% slow | 3% slow | right |
| Milky Way pull above the disk | **73 against 68–74** | 84 | 74 |
| Milky Way mass inside 100 / 200 kpc | **6.5 / 12.4 against 6.1–7.3 / 11.0** × 10¹¹ | 8.0 / 15.6 | 8.5 / 12.8 |
| Escape speed at the Sun | **509–525 against 445–580** | 546–564 | 560–570 |
| Ten Milky Way dwarfs | 4 agree, 6 too slow (6 agree without the Milky Way's hold, §6.13) | the same | fitted |
| Planets, S2, pulsars, light bending | **Einstein's, exactly** | small effects | Einstein's |
| Cassini Q2 | **passes** (release over 0.15 pc) | about 10× too big | passes |
| Solar System | **silent** | small effect | silent |
| Explains *why* | disks are cold; cluster galaxies are hot; gas collides | no | no |
| Adjustable numbers | **4**: 3 fitted, 1 release length set by Cassini | 1 | ~320 |

Dark matter fits individual objects more tightly because it is tuned object by object. Ours
fits everything with four shared numbers, three fitted to galaxies and clusters and one set by
the Solar System. It also says why each kind of system behaves as it does.

**Every standard test, in one list (rev 21).** ✓ = agrees within the errors; ~ = close (within
about two error bars); ✗ = a shortfall, with the fix we are testing.

| Test | Measured | Ours | |
|---|---|---|---|
| 149 galaxies' rotation (SPARC) | | 15.9 km/s typical miss | ✓ |
| 12 clusters' masses (X-COP), our distances, stars in spheres | | 25% typical miss | ✓ |
| Milky Way rotation, 15–27 kpc | 173–217 km/s | 2–12 km/s low | ✓ |
| Milky Way rotation at the Sun | 229–234 km/s | 209 (217 with a more compact disk, rev 17's constants) | ~ the Milky Way's disk shape is the lever to test (§9) |
| Pull above the Milky Way's disk | 68 ± 4, 74 ± 6 | 73 | ✓ |
| Milky Way mass inside 100 / 200 kpc | 6.1–7.3 / 11.0 × 10¹¹ | 6.5 / 12.4 | ✓ |
| Milky Way mass inside 20–50 kpc | 1.9–4.5 × 10¹¹ | 13–21% low | ~ (same cause as the Sun's speed) |
| Escape speed at the Sun | 445–580 km/s | 509–525 | ✓ |
| Inner Galaxy: share of pull from stars | 0.88 ± 0.07 | 0.97–0.99 | ~ |
| Dwarfs Fornax, Leo I, Leo II, Sculptor | 6.6–11.7 km/s | within 1–2 error bars | ✓ |
| Six faint or spread-out dwarfs | 2.7–9.5 km/s | 1.5–5× too slow | ✗ with no hold by the Milky Way, 2 agree and 2 are close; Draco and Ursa Minor stay 2× slow (§6.13) |
| Planets, Sedna, S2, pulsars, light bending | Einstein's | Einstein's, exactly | ✓ |
| Cassini Q2 | (3 ± 3) × 10⁻²⁷ s⁻² | 4.4 × 10⁻²⁷ (release over 0.15 pc, rev 18) | ✓ |
| Wide binaries | disputed | 8% at 20,000 AU, 3% at 7,000 AU | open |
| Microlensing, Einstein Cross | stars only | stars only | ✓ |
| Six strong lenses (SLACS), our own distances | | light = matter (−0.028 dex) | ✓ |
| Galaxy lensing (KiDS), spirals and ellipticals, our own distances | | the ellipticals' extra lensing matches (0.13–0.16 against 0.15); all lenses together 0.065 dex above ours, closed by gas around them as heavy as their stars | ~ the gas to be weighed (§6.16) |
| Flat lensing speeds to 1 Mpc | | flat to 2.25 Mpc | ✓ |
| Bullet Cluster: lensing on galaxies, main half | | yes | ✓ |
| Bullet Cluster: smaller half's lensing mass | 2.0–2.3 × 10¹⁴ | about half | open |
| 72 collisions: lensing with the galaxies | | yes | ✓ |
| MACS J0025.4−1222: lensing inside 300 kpc, SE / NW (our distances) | 3.6 / 3.8 × 10¹⁴ | 2.1 / 1.9 | ✓ |
| MACS J0025.4−1222: lensing peaks | on the galaxies | both on their galaxies at the age its shock fronts give (0.1–0.36 billion years; on the gas beyond 0.45) | ✓ |
| MACS J0025.4−1222: galaxies' speed spread | 835 ± 59 km/s | 669 (stars on the calibrating basis); about 800 with the published masses | ~ stars to measure |
| Abell 520: lensing of six clumps inside 150 kpc (our distances) | 2.3–6.2 × 10¹³ | 5 of 6 within the errors (P6 2.3σ low) | ✓ |
| Abell 520: the galaxy-poor "dark core" P3 (our distances) | 3.29–3.89 × 10¹³ | 3.72 | ✓ |
| El Gordo: lensing inside 1 Mpc (our distances) | 24.3 × 10¹⁴ | 21.4 (published stars) | ✓ |
| Abell 1689 | collected | not yet run | next |
| The regression suite (§6.12) | 77 graded checks | 59 pass, 11 close, 7 fail (everything in our distances) | |

## 8. Predictions anyone can check

1. **Heavy stars in the six SLACS lenses:** 1.37–1.85 times Salpeter with the project's own
   distances (1.05–1.35 with the standard ones), measurable from their spectra. Dark-matter
   models expect about 1.0.
2. **The elliptical/spiral lensing gap is flat beyond 100 kpc.** A hot-gas-halo explanation
   makes it grow with radius.
3. **In relaxed clusters, galaxy orbits and gas agree.** Both feel the same pull.
4. **In every collision, lensing stays with the galaxies**, at every stage.
5. **The Bullet Cluster's smaller half carries a companion heated during its crossing,** so its
   lensing is strong although few galaxies travel with it. Deep images should also show
   whatever diffuse starlight travels with it.
6. **The Bullet Cluster's main outskirts** hold stars at about 5% of the gas mass at the
   cluster's edge. Confirmed by the star count in rev 15.
7. **Wide binary stars:** 4% more pull than Newton at 7,000 AU and 9% at 20,000 AU, a fifth of
   MOND's 43% (with the release length adopted in rev 18, §6.10); 3.5% and 8% with rev 20's
   constants, 3.2% and 7.6% with rev 21's, inside the locked window within its errors. If the Milky Way's hold turns out weaker (§6.13), 19–36%.
   **Locked in the repository before Gaia's next data release (rev 19, §6.7).**
8. **Lensing in a collision shows the clusters as they were before it.** Around gas that has
   been stopped, lensing comes back only inside a sphere growing at 163 km/s, about 170 kpc per
   billion years.
9. **Older collisions** (half a billion to a billion years after the pass) show extra lensing
   around the smaller clump, as the stream heat and tidal shaking build up.
10. **The extra pull tracks random motion.** At equal visible mass, systems whose stars move
   randomly and freely pull harder than those whose stars circle in step.
11. **The Sun's extra mass loss** is 1.9 × 10⁻¹⁵ per year.
12. **The Milky Way keeps slowing down beyond 25 kpc:** about 189 km/s at 30 kpc and 179 at
   50 kpc, not flat.
13. **Lensing speeds around isolated galaxies stay flat to about 1.7–2.2 Mpc and then fall**,
   where the companion has not yet reached.
14. **The Einstein Cross and similar compact lenses need only their stars** inside the ring.
15. **Ellipticals and spirals follow two different lensing curves**, set by how randomly their
   stars move (KiDS already agrees).
16. **The far clusters' published star masses are about right** (rev 20, §6.15). Rev 19 predicted
   1.4–2 times heavier, older stars; with the companion's speed measured in our own distances
   that is no longer needed, and distant cluster stars are observed to be younger. El Gordo's
   lensing agrees with its published stars. Near-infrared imaging and spectra can still check them.
17. **Abell 520's galaxy-poor clump carries about 3 × 10¹³ suns inside 150 kpc**, of which a
   quarter is visible gas. Its lensing should follow the gas as that gas settles.
18. **MACS J0025's clusters crossed less than about 400 million years ago.** Otherwise, in our
   law, its north-western lensing peak would sit on the gas, where none is seen (§6.15). *Rev 21: its
   radio relics and its separation date the crossing at 0.1–0.36 billion years (§6.16).* Deeper X-ray
   data that measure the shocks' Mach numbers would pin it down. Its star
   masses set its galaxy speeds: 665 km/s with the stars on the calibrating clusters' basis, 804
   with the published masses, against 835 ± 59.
19. **Isolated galaxies are wrapped in gas weighing about as much as their stars within 100 kpc,
   more around ellipticals than around spirals.** With that gas our law matches KiDS's lensing in our
   own distances (§6.16). X-ray and ultraviolet measurements of the gas around such galaxies can check
   it. *(Rev 20's prediction 19, that galaxy lensing measures our distance law's rate, was withdrawn
   in rev 21: refitted together, the rate barely moves it.)*

## 9. What is still open, and why we are optimistic

**Settled so far:** the main Bullet cluster's galaxy speeds and lensing mass (rev 14), its star
count (rev 15), Cassini's test (rev 18), and the far clusters' star masses (rev 20).

For each open item below: what is wrong, what would have to change, and the next test. They are
ordered by how directly they can move the scoreboard, and each ends in a regression-suite run,
so we will know at once what a change fixes and what it costs. §9.1 adds ten longer-range
proposals for finishing the theory.

1. **Take out what the far clusters borrowed from the expanding universe (new in rev 18).**
   * *What is wrong:* our tests of the far clusters convert angles into kiloparsecs, and light
     into star masses, with the standard expanding-universe distance law. The project's own
     static law gives different numbers. In it, redshift comes from light slowly losing energy,
     1 + z = e^(αD). The difference grows with distance. At a fixed angle, static ÷ expanding:

     | | z | sizes | star masses | gas masses | lensing masses |
     |---|---|---|---|---|---|
     | Abell 520 | 0.20 | ×1.08 | ×0.80 | ×1.00 | ×1.10 |
     | Bullet Cluster | 0.30 | ×1.14 | ×0.78 | ×1.08 | ×1.18 |
     | MACS J0025 | 0.59 | ×1.36 | ×0.73 | ×1.36 | ×1.39 |
     | El Gordo | 0.87 | ×1.58 | ×0.71 | ×1.68 | ×1.60 |

     *Lensing for typical background-galaxy distances.*

     The six strong lenses (§6.4) were computed both ways. Light and matter agree under either
     law: −0.012 ± 0.023 dex with ours, −0.017 ± 0.024 with the standard one. Since rev 19 the
     suite grades ours.
   * *Two more borrowed pieces:*
     * Star-mass fits for distant galaxies usually forbid stars older than the Big Bang allows at
       that distance: 6.3 billion years at El Gordo's, 7.9 at MACS J0025's. Older stars weigh
       more for the same light. That may be why the far clusters now need 1.4–2 times their
       published stars (§6.14).
     * Some published lensing masses are fits of a dark-matter-shaped profile (NFW), such as El
       Gordo's. Masses measured inside apertures are cleaner.
   * *What would have to change:* nothing in the law, only its inputs. Whether El Gordo's gap
     closes or widens has to be computed: its gas rises and its stars fall.
   * *Done in rev 19 (§6.14):* MACS J0025, Abell 520, El Gordo and the six strong lenses, with each
     paper's own assumptions about the universe and its background galaxies, and El Gordo's
     aperture masses instead of its dark-matter-shaped fit. The suite now grades them this way:
     58 pass, 11 close, 7 fail. Every grade that slipped is a far cluster whose stars came out
     lighter.
   * *Done in rev 20 (§6.15):* the Bullet Cluster, KiDS, Mistele's lensing speeds and X-COP in our
     distances, and the constants refitted there. X-COP's stars turned out to have been counted in
     projection since round 1; corrected, u = 163 km/s. The far clusters then recover without
     heavier stars. The age cap turned out not to matter: distant cluster stars are younger.
   * *Next:* what our distances sharpened, galaxy lensing (new item 11).

     The nearby tests (SPARC, X-COP, the Milky Way) change by only a few percent, mostly because
     our distance law's scale corresponds to a Hubble-constant-like 74.6 rather than 70 (70.9
     since rev 21, fitted to supernovae).
   * *One more to watch:* the companion's reach, 2.2 Mpc, assumes matter has been emitting for
     13 billion years, the age of the oldest stars. Without a Big Bang, matter may be older than
     its stars. No current test reaches that far, but prediction 13 depends on it, so the reach
     should be measured rather than assumed.
2. **Six faint dwarf galaxies (§6.13).**
   * *What is wrong:* Draco, Ursa Minor, Sextans, Crater II and Antlia 2 come out 2.7–5 times too
     slow; Carina is close.
   * *What would have to change:* the Milky Way's pull must dilute a dwarf's own companion much
     less than our law says. Without the dilution, 6 of the 10 agree and χ² falls from 135 to 60,
     against MOND's 119. The Cassini fix already keeps the binaries safe.
   * *Next:*
     * find the physical reason (§9.1, proposal 5). It must leave galaxies, clusters and the
       heat rule alone, and match the dwarf-by-dwarf pattern: a little hold for Crater II, none
       for Antlia 2;
     * turn it into a suite candidate.
   * *Also next:* Draco and Ursa Minor stay about half as fast as measured whatever the hold.
     * Compute how far the Milky Way's tides reach into them under our law, and compare with how
       far their stars extend.
     * Look in Gaia and spectroscopic surveys for the stretching and velocity trends that tides
       would leave.
3. **The Bullet Cluster's smaller half (§6.3).**
   * *What is wrong:* its lensing mass inside 250 kpc is 0.94 × 10¹⁴ suns against a measured
     2.0–2.3 (rev 20, our distances: 1.29 inside 286 kpc against 2.35–2.70). Rev 14's idea, a
     bigger cluster before the crash, is not borne out by its galaxies and starlight.
   * *What would have to change:* more pull around it, from something travelling with it.
   * *Next:*
     * build the heat its galaxies picked up while crossing the main cluster (§6.3) as a suite
       candidate, and in time a full transport calculation (§9.1, proposal 6);
     * redo it with our own distances (item 1): *done in rev 20*;
     * deep images for diffuse starlight moving with it, and more galaxy speeds around it.
4. **The Sun's orbital speed and the Milky Way's middle (§6.8).**
   * *What is wrong:* 211 km/s against 229–234 at the Sun, and the mass inside 50 kpc is 20% low.
     Both are "close".
   * *What would have to change:* more pull at 8–50 kpc.
     * Changing the law does it, but costs the strong lenses and the pull above the disk (§6.12).
     * The Milky Way's own matter is less certain than the law. Bovy & Rix's more compact disk
       already gives 217 km/s, which would pass.
   * *Next:* add the disk's length and the bar's mass as Milky Way options in the suite. Grade
     them against the Sun's speed, the pull above the disk, the 15–27 kpc curve and the masses
     together.
5. **Abell 1689.**
   * *What it is:* a heavy, settled cluster that bends light strongly. MOND needs extra unseen
     mass there.
   * *Next:* model it from its published gas and stars, with our distances, and add it to the
     suite. Grade it against its lensing mass profile and its Einstein ring.
6. **How much cluster stars weigh (§6.11, §6.15).**
   * *Rev 20:* with the companion's speed measured consistently, the far clusters no longer need
     heavier stars. Each cluster's star masses are now on the calibrating clusters' basis.
     MACS J0025 is the exception: its galaxy speeds prefer its published masses, which are on a
     heavier basis. *Rev 21:* its north-western lensing peak is settled by its age, which its shock
     fronts put at 0.1–0.36 billion years (§6.16).
   * *Next:*
     * MACS J0025's star masses from infrared light;
     * star masses from infrared light and spectra, with declared assumptions (§9.1, proposal 9);
     * Abell 520 modelled as two main clusters before the crash, not one group per clump.
7. **Lensing speeds of the lightest spirals (§6.9).**
   * *What is wrong:* ours are 20–30% below Mistele et al.'s in the standard distances, 17–40% in
     ours, but agree with Brouwer et al.'s conversion of the same KiDS data in the standard ones.
   * *Next:* grade the law against the published lensing profiles directly, without either
     conversion (and item 11).
8. **Wide binary stars (§6.7).**
   * *The prediction:* 4% at 7,000 AU and 9% at 20,000 AU (3.5% and 8% with rev 20's constants,
     3.2% and 7.6% with rev 21's);
     19–36% if item 2 goes the way of no hold. The two published analyses disagree (0 and about
     40%).
   * *Next:* Gaia's fourth data release. The binaries will measure the release length and decide
     between the two routes for the dwarfs.
9. **The theory**, the subject of most of §9.1:
   * *Done in rev 19 (§4.14):* the pull, derived from a local interaction, with the three
     conditions it imposes; and the way the companion adds up, tested against two alternatives.
     Only the law's own rule fits both galaxies and clusters.
   * *Done in rev 20 (§4.15):* energy conservation derives how orderly matter adds up, given one
     stream and no whirlpools.
   * *Done in rev 21 (§6.16):* a working model of five local rules. Only a companion guided along
     gravity's field lines keeps every watt, travels as one stream and forms no whirlpools, and a
     short argument shows it must follow Newton's pattern exactly.
   * *Done in rev 22 (§6.17):* the pull and the scrambling in one experiment. With one fixed timing, heat weakens
     the pull, for an exact reason. Emitters that feed quiet waves but absorb loud ones give the right pattern at
     slow speeds: a pull up to twice as strong with free random motion, held back by collisions, none from rotation.
   * *Done in rev 23 (§6.18):* the rule derived from a power balance (a fixed supply, a loss that grows with
     loudness); its hold fits the galaxies as well as the release factor; the darker cold matter is, the wider the
     range where heat dominates and bodies keep step.
   * *Done in rev 24 (§6.19):* an independent calculation's quiet internal store, opened by motion, gives the σ²; in
     one experiment with real waves, the pull grows as √(cold + released) and doubling σ doubles the extra pull. An
     inverted, self-sustained body falls into step a quarter beat ahead by itself (round 10's rule, derived). A cloud's
     own quietness gives the symmetries and the collision rule but not the σ², because it has no clean edge.
   * *Next:* a physical quiet store with a gap, and u from it; one kind of matter as both source and receiver;
     whether a drained receiver can switch its pull off (the release factor); derive why the hold stops at zero
     instead of pushing; then the guided companion with moving sources, for the collisions' memory;
   * a full field theory for the companion, with its travel, its release length and a
     relativistic form, so that lensing is derived;
   * what the release length is;
   * u and g_d derived rather than fitted.
10. **Cosmology.** The microwave background and the growth of large-scale structure are outside
   the law's tested scope. Under the project's no-expansion rule they need their own treatment.
11. **Galaxy lensing in our own distances (new in rev 20, §6.15).**
   * *Rev 21:* the ellipticals' extra lensing matches once the lenses' star speeds are measured
     (0.13–0.16 against 0.15; §6.16). What remains is a common level, +0.065 dex for all lenses
     together. It is not our distance law's rate: fitted jointly with supernovae and the galaxies,
     the rate barely moves it. It measures the lenses' mass: gas around them weighing as much as
     their stars, the KiDS team's own middle estimate, would close it.
   * *What was wrong (rev 20):* around spirals and disks the lensing is 20–28% stronger than our law
     predicts (+0.08 and +0.11 dex), and ellipticals' extra lensing over spirals is 0.23 dex
     against 0.15 ± 0.04.
   * *What rev 20 thought would have to change:* for the first, the distance law's scale (rev 21:
     not so); for the second, the ellipticals' heat (rev 21: measured, it matches).
   * *Done in rev 21:* α refitted jointly on supernovae, SPARC's distances, X-COP and KiDS; the
     ellipticals' star speeds measured; bottom-heavy stars checked (too small a change at these
     masses).
   * *Next:* weigh the gas around galaxies like KiDS's lenses, from X-ray stacks (ellipticals) and
     ultraviolet absorption (spirals), and predict the lensing level from it rather than fit it;
     check what such gas does to the outermost points of the 149 rotation curves.

Every item is a concrete calculation or measurement, and the suite now tells us within minutes
whether a change helps.

### 9.1 Ten proposals for finishing the theory, and how we would work on them

Rev 18 adds ten proposals from the project's owner. Each describes what a *completed* solution
would look like; **none of them is a result yet.** For each: where we stand, how we would work
on it, what to watch out for, and the first step.

The standing rules apply throughout: nothing borrowed from dark matter, an expanding universe or
a Big Bang, even indirectly. Item 1 above shows three places where such assumptions had crept
into our inputs.

| # | Proposal | Where we stand | First step | How hard |
|---|---|---|---|---|
| 1 | Derive the whole law from one microscopic interaction | **rev 19: the pull derived and simulated; of three ways to add up, only the law's fits galaxies and clusters (§4.14). Rev 20: that rule derived from energy conservation (§4.15). Rev 21: a companion guided along gravity's field lines gives it with no further assumption (§6.16)** | the guided companion with moving sources: does it keep its memory in a collision? | hard; the foundation |
| 2 | Predict unseen data with a locked model | **rev 19: the wide-binary forecast locked (§6.7)** | lensing profiles by speed spread for a new galaxy sample | doable now |
| 3 | Show that random motion makes a steady pull and collisions remove it | **rev 22: done in one experiment (§6.17). One fixed timing gives the opposite, for an exact reason; a timing that depends on the wave's loudness gives the pattern at slow speeds (×2 free, less with collisions, none rotating). Rev 23: that timing derived from a power balance, with no memory; bodies keep step up to the re-timing speed (§6.18). Rev 24: with a quiet internal store opened by motion, the full pattern: pull ∝ √(cold + released), ×2 per doubling of σ where heat dominates, collisions suppress as 1/(1 + ν/γ) (§6.19)** | a physical quiet store with a gap, and u from it | medium; the test bench exists |
| 4 | Light and matter from one coupling | lensing assumed to follow the pull; six lenses agree | a light-like mode in the same toy | hard; after 1 |
| 5 | One mechanism for Cassini, dwarfs and wide binaries | release adopted; the dwarfs need a weaker hold, reason unknown | derive the release time and the hold from one transition | hard; biggest payoff |
| 6 | Evolve the companion through collisions | memory imposed; the Bullet's smaller half at half weight | a time-dependent transport solver | medium to hard |
| 7 | Measure u independently | u = 169 km/s from the clusters (rev 21: every data set in our distances); the lensing gap agrees with it once the lenses' star speeds are measured (rev 21) | the lensing signature of the growing fresh companion | medium; needs data |
| 8 | Derive the constants | fitted; clues such as a ≈ cα/11 | the toy's stiffness and inertia give u | hard; after 1 |
| 9 | Stellar masses without using gravity | **rev 20: every star mass on one basis; the calibrating clusters' corrected; no age-cap correction needed (§6.15)** | infrared light and spectra; the ellipticals' IMF | doable now |
| 10 | Detect the energy cost of making gravity | the Sun loses 1.9 × 10⁻¹⁵ of its mass a year to the companion. **Rev 24: being pulled costs power too, force × the companion's speed: 2ℓ per kilogram at the law's a with a slow companion; a light-speed one is excluded by the planets (§6.19)** | precision ranging | medium; the test is future |

**1. One interaction, the whole law.**
* *Where we stand:* the law is assembled from nine postulates, and two have toy models: the
  Dicke toy (§3.6) and the 3D field equation (§4.9).
  * *Rev 19 (§4.14):* the weakest link below is now derived, and checked by a 3D simulation. The
    way the companion adds up is pinned down by the data: of three natural rules, only the law's
    own fits both galaxies and clusters.
* *The weakest link* is postulate 4, "the companion pulls with its amplitude". A field's steady
  push normally scales with its energy, the amplitude squared, and that would give Newton's
  1/r² rather than the 1/r that flat rotation curves need.
* *How we would work on it:* build a toy medium. Sources feed a wave field with a bound
  (near-source) state and a free (travelling) state, and test bodies move according to their
  interaction energy with the field. The most promising route to a pull proportional to the
  amplitude is interference. A body's own companion and a passing one share an interaction
  energy proportional to the product of their amplitudes, so it is linear in the passing one.
  Then check whether the release factor, the heat term, the direction rule and the field
  equation all come out of the same toy.
* *Watch out for:* borrowing a ready-made relativistic theory. The known ones in this area
  (TeVeS, AeST) are built on MOND, which the rules exclude.
* *First step:* a one-dimensional version with two sources and a test body. Measure the
  time-averaged force against distance and source strength. *Done in 3D in rev 19. Rev 20:
  energy conservation fixes how orderly matter adds up, given one stream and no whirlpools
  (§4.15). Rev 21: a companion guided along gravity's field lines supplies both, in a working model
  and by a short proof (§6.16). Next: the same with moving sources.*

**2. A locked prediction.**
* *Where we stand:* the suite re-checks known data; passing it again is not a new prediction.
  §8 lists 18 predictions.
* *How:* freeze the law, its constants, the distance law and the analysis code in a tagged
  commit, then publish forecasts before looking at the data. Candidates:
  * Gaia's fourth data release, for wide binaries: 4% extra pull at 7,000 AU and 9% at
    20,000 AU;
  * lensing around galaxies split by their stars' speed spread (KiDS-Legacy, Euclid): full
    profiles per bin, from each galaxy's own light and spectrum;
  * El Gordo's star mass from infrared light (proposal 9).
* *Watch out for:* the conversions. New samples come with distances and star masses computed
  for an expanding universe, and they must be redone with ours first.
* *First step:* commit the wide-binary forecast file now; its timestamp is the record. *Done in
  rev 19 (§6.7).*

**3. Random motion makes a steady pull, and collisions take it away.**
* *Where we stand:* the Dicke toy showed that randomly moving emitters add as plain totals
  (210.55 against 209.20 predicted), and that emitters colliding 100 times per wavelength stay in
  step (1.2 times the arrow sum, against 11 times for free ones). It did not measure the pull on
  test bodies, and k = 3σ²/u² came from averaging the Doppler boost by hand.
* *How:* add test bodies to the toy.
  * Measure the time-averaged force for populations with the same mass profile and different
    speed spreads. Doubling σ should give four times k, and twice the extra pull where heat
    dominates.
  * Raise the collision rate at fixed energy and watch the extra pull fade smoothly, with no
    "gas" switch.
  * Split one population into arbitrary catalogue groups and check that nothing changes.
* *Watch out for:* a "force" that is only jitter. Average over long times and check its sign.
* *First step:* reuse the toy's emitters with proposal 1's interaction. *Done in rev 22 (§6.17), with 135 runs plus
  checks. With one fixed timing, free random motion halves the pull (the reason is exact). With a timing that depends
  on the wave's loudness, free random motion strengthens the pull up to ×2, collisions hold that back, and rotation
  adds nothing, at speeds up to a quarter of the re-timing rate. Next: the rule's physics, and bodies that keep in
  step at higher speeds.*

**4. Light and matter from one coupling.**
* *Where we stand:* we assume light feels the same landscape Φ as matter, with Einstein's
  factor 2 for bending. The six SLACS lenses agree to 0.017 dex and KiDS agrees, but that
  supports an assumption rather than a derivation.
* *How:* give proposal 1's medium a fast, light-like mode, and derive how it bends and how slow
  bodies fall through the same disturbance. Newton's own calculation for light gives half the
  measured bending; the factor 2 comes from gravity affecting clocks and distances equally, and
  our medium would have to produce both.
* *Watch out for:* importing Einstein's geometry as an assumption.
* *First step:* after proposal 1, trace rays through the toy's field.

**5. One mechanism for Cassini, the dwarfs and wide binaries.**
* *Where we stand (rev 18):*
  * the release length is adopted: Cassini passes, and binaries are at 4% and 9% (3.5% and 8%
    with rev 20's constants);
  * the dwarfs need a much weaker hold by the Milky Way;
  * a release length of 0.5 pc or more already protects the binaries;
  * there is no mechanism for the weaker hold yet.
* *How:* derive both from one bound-to-free transition.
  * Its relaxation time gives the release, R = 1 − e^(−t/τ), with τ calculated rather than
    chosen.
  * The same transition should say how two companions lose overlap when they move relative to
    each other.
  * The test: co-moving binaries keep the hold, dwarfs moving at 100–300 km/s lose most of it,
    and the parts of one rotating galaxy stay in step, so the heat rule survives.
* *Watch out for:* your completion criterion is the right one. The dwarf deficits must close,
  not shrink. Draco and Ursa Minor need something beyond the hold, and tides come first
  (item 2 above).
* *First step:* in the 1D toy, measure the release time of a bound disturbance and the overlap
  between two sources in relative motion.

**6. Evolve the companion through collisions.**
* *Where we stand:* the memory rule is imposed. The old companion rides with the galaxies and a
  fresh sphere grows around stopped gas. It explains the 72-collision stack, MACS J0025 and
  Abell 520's gas-centred clump, but gives the Bullet's smaller half half its measured lensing
  mass.
* *How:* a transport equation for the companion: emission, travel with its source plus
  streaming at u, release, and heat from streams passing through each other.
  * Solve it in time for two clusters built from their observed gas and galaxies: the Bullet
    first, then the others with the same code.
  * Derive when a lensing peak can sit on gas, from density, surrounding heat and time since the
    gas stopped.
* *Watch out for:* starting conditions or timings taken from dark-matter simulations. Use
  observed gas and galaxies, our own dynamics and our own distances.
* *First step:* a 2D transport solver, checked against today's memory rule in the steady case.

**7. Measure the companion's speed independently.**
* *Where we stand:* u = 169 km/s (rev 21, every data set in our own distances, the stars correctly
  counted) comes only from the 12 X-COP clusters, through the heat term. *(Rev 20 read KiDS's
  ellipticals as preferring about 200 km/s; rev 21 traced that to assumed star speeds, §6.16.)*
* *How:*
  * Derive how the fresh companion around stopped gas shows in projected lensing as it grows:
    25 kpc after 150 million years, 83 kpc after 500 million, 166 kpc after a billion.
  * Pick collisions whose timing is known from their shock fronts (X-ray and radio), not from
    lensing or dark-matter simulations.
  * Fit u across several systems.
* *Watch out for:* the other clock, the reach. We took 13 billion years of emission, the age of
  the oldest stars. Without a Big Bang, matter may be older, so the outer lensing transition of
  isolated galaxies (prediction 13) measures u times the emission time rather than u alone.
* *First step:* the projected lensing profile of a growing fresh sphere around a β-model gas
  cloud.

**8. Derive the constants.**
* *Where we stand:* a, g_d and u are fitted and L is bounded, with g_d/a = 3.22.
  * A clue: a ≈ cα/11 and g_d ≈ cα/3.4, where α is the rate at which light loses energy in the
    project's static universe. That is a local process, not expansion. Galaxy lensing now
    measures α too (§6.15).
* *How:* in the toy medium, stiffness κ and inertia χ set u² = κ/χ. The source coupling sets ℓ.
  The binding sets both the escape length L_d (g_d = u²/2L_d, L_d ≈ 2.0 kpc) and the escape time
  (L = 0.15 pc of travel), two different processes, a barrier height and a rate. If light loses
  its energy to the same medium, α and ℓ may share a cause, which would predict a/cα and g_d/cα.
* *Watch out for:* working numbers backwards from the fits and calling them derived.
* *First step:* after proposal 1.

**9. Stellar masses without using gravity.**
* *Where we stand (rev 20, §6.15):* every cluster's star masses traced to their source and put on
  one basis (`literature/star_mass_audit_v11.md`). The calibrating clusters' had been counted in
  projection; corrected. With the companion's speed measured consistently, the far clusters need
  no heavier stars. The SLACS lenses need 1.37–1.85 times Salpeter.
* *The age cap, checked:* distant cluster galaxies are observed to be bluer, with weaker 4000 Å
  breaks, and lighter for their light, measured without any distance law. Their stars are
  younger, so the cap does not hide mass.
* *How:*
  * Infrared photometry and spectra, with declared assumptions.
  * The IMF: spectra say giant ellipticals' stars are heavy for their light (about 1.85 times
    Chabrier's rule). Test it in the clusters and in KiDS together (§9, item 11).
  * Lock, and predict lensing and galaxy speeds.
* *Watch out for:* this is where borrowed assumptions hide most, because star-mass codes
  build in a timeline.
* *First step:* the ellipticals' IMF in X-COP and KiDS together. Within reach now.

**10. The energy cost of making gravity.**
* *Where we stand:* each kilogram feeds ℓ = au/2 = 5.32 × 10⁻⁶ W into its companion. Paid from
  mass (E = Mc²), that is a loss of 1.9 × 10⁻¹⁵ of the Sun's mass a year.
* *A consequence we can already draw:* the energy cannot come from heat.
  * The Earth's companion power would be 3.2 × 10¹⁹ W, 680,000 times the heat flowing out of its
    interior.
  * A typical white dwarf's would be about 16 times its light.
  * Both would cool almost at once, and they do not. So the companion must be paid for by mass
    itself.
* *How:* check energy and momentum balance in the toy (proposal 1), including the recoil of
  lopsided sources. Then look for the mass loss. It acts like gravity weakening by
  1.9 × 10⁻¹⁵ a year, about 35 times below the precision of today's lunar laser ranging.
* *First step:* the energy balance in the toy; then follow planetary and lunar ranging as it
  improves.

## 10. How we got here

* **Revisions 1–9.** A "streaming gravity" law plus a light-versus-matter "slip". The law
  turned out to be MOND's, so under the new rules it was retired (revision 10). The slip later
  proved an artefact of measuring against MOND.
* **Revision 10 (round 1).** Hot-companion gravity, with the gas as the hot component. It fit
  galaxies and clusters (u = 874 km/s), but put the Bullet's lensing on the gas.
* **Revision 11 (round 2).**
  * MOND derived as the cold limit.
  * A field equation, and an action to conserve momentum.
  * Hot gas haloes to explain the ellipticals.
  * The Bullet partly solved by "companion memory", but that failed the 72-collision stack.
* **Revision 12 (round 3).** The collision rule and the flow direction. Stars now carry the
  heat (u = 197 km/s), and the Bullet's pattern is reproduced.
* **Revision 13 (round 4).** The companion's memory, now safe because the companion is slow.
  The 72-collision stack still holds. We also corrected rev 12's galaxy-speed check for the
  main cluster, which showed a gap to close.
* **Revision 14 (round 5).**
  * §4, the law piece by piece.
  * The main cluster's galaxy speed matched, half by removing a modelling artefact and half with
    a normal amount of stars in its outskirts.
  * Both strong-lensing masses within 1–2σ, with the smaller half a third of the main before the
    crash.
* **Revision 15 (round 6).**
  * A Legacy Survey star count confirms the main cluster's outskirts.
  * The galaxy speeds and starlight don't show the stars a bigger smaller half would need, so
    its lensing strength is open again.
  * The wide-binary prediction.
* **Revision 16 (round 7).**
  * The full check against the Milky Way and the standard stars and lenses, with every
    measured value traced to its paper.
  * Passes: the outer Milky Way, its mass and escape speed, the pull above the disk, S2,
    pulsars and planets, and galaxy lensing of spirals and ellipticals (better than MOND).
  * Shortfalls, each with a candidate fix: the Sun's speed, Cassini's Q2, six faint dwarfs.
* **Revision 17 (round 8).** A regression suite: every test in one command, graded the same
  way, compared with a saved baseline. Three more collisions: MACS J0025 passes; Abell 520's dark
  core comes from gas and heat; El Gordo's stars are the number to measure.
* **Revision 18 (round 9).**
  * The Cassini fix adopted, after the full suite showed it costs nothing: 62 pass, 8 close,
    6 fail.
  * It also unblocks the dwarf fix: with no hold by the Milky Way, 6 of 10 dwarfs agree
    (χ² 60 against MOND's 119). The physical reason is still to be found.
  * A plan for every open item.
* **Revision 19 (round 10).**
  * The pull derived from a local interaction and checked by a 3D simulation.
  * The way the companion adds up tested against two alternatives: only the law's own rule fits
    both galaxies and clusters, which fixes three properties of the companion.
  * The far clusters and the strong lenses in the project's own distances, now graded that way:
    58 pass, 11 close, 7 fail. The slipped grades all trace to lighter stars.
  * A first locked forecast (wide binaries).
* **Revision 20 (round 11).**
  * How orderly matter adds up, derived from energy conservation; the data demand it to within a
    few percent.
  * The calibrating clusters' stars corrected (counted in projection since round 1), and every
    test in the project's own distances. The companion's speed, measured again: 163 km/s. The far
    clusters recover without heavier stars.
  * Galaxy lensing becomes a test of our distance law's scale and of the companion's speed.
    Suite: 57 pass, 9 close, 10 fail.
* **Revision 21 (round 12).** Rev 20's next steps, one at a time:
  * MACS J0025's age, read from its shock fronts: young, and its north-western lensing peak sits on
    its galaxies;
  * the lensing galaxies' star speeds, measured from SDSS spectra: the ellipticals' extra lensing
    matches;
  * the distance law's rate, fitted jointly to supernovae and the galaxies: 5% lower, with every data
    set in our own distances. Galaxy lensing's remaining level points to the gas around the lenses.
    Suite: 59 pass, 11 close, 7 fail;
  * a working model of the companion: only a companion guided along gravity's field lines keeps every
    watt, travels as one stream and forms no whirlpools.
* **Revision 22 (round 13).** The pull and the scrambling in one experiment:
  * with one fixed timing, heat weakens the pull, and the reason is exact;
  * emitters that feed quiet waves but absorb loud ones give the right pattern at slow speeds;
  * galaxy lensing pins the heat term's σ² (doubling σ gives 1.8–2.0 times the extra pull).
* **Revision 23 (round 14).** Rev 22's next steps, one at a time (§6.18):
  * a gentler heat term (σ^1.75) on the full suite: it closes galaxy lensing's common level but costs the massive
    ellipticals, so σ² stays;
  * "feed quiet waves, absorb loud ones" derived from a power balance (a fixed supply, a loss growing with loudness);
    the same loss holds the companion back in a shape the galaxies accept as well as the release factor;
  * how far the release goes: the darker cold matter is, the earlier heat takes over, and bodies keep step up to about
    the re-timing speed; the release's steepness is the open question.
* **Revision 24 (round 15, this page).** An independent calculation joined to the pull (§6.19):
  * its quiet internal store, opened by motion, reproduced to the last digit;
  * the wave alone gives its symmetries and the collision rule, not its σ²: a cloud's quietness has no clean edge;
  * inverted, self-sustained matter falls into step a quarter beat ahead by itself: round 10's rule, derived;
  * in one experiment with real waves, the pull grows as √(cold + released) and doubling σ doubles the extra pull;
  * being pulled costs force × the companion's speed, which requires a slow companion.

**Superseded along the way, kept on the record:**
* round 2's hot-gas-halo explanation of the ellipticals (now it is their stars);
* round 2's reaction force on hot matter (now the companion carries momentum);
* revs 12–21's mechanism for the heat term, "scrambled contributions don't cancel, so they pull
  harder": joined with the pull in one experiment, scrambling alone weakens it (rev 22, §6.17);
* rev 23's guess that a larger ball or a steeper loss would make the release grow as σ²: a cloud's own quietness has
  no clean edge and grows too gently however it is arranged; the σ² needs a quiet store inside each piece of matter
  (rev 24, §6.19);
* round 1's companion speed (874 km/s; 197 in revs 12–19; 163 since rev 20, measured in our own
  distances with the clusters' stars correctly counted);
* round 2's solar mass-loss figure (1.4 × 10⁻¹⁴ per year; 2.3 × 10⁻¹⁵ in revs 12–19; now
  1.9 × 10⁻¹⁵);
* rev 12's "the main cluster's galaxy speeds check out" (a single radius; the average over the
  surveyed region was 20% low in rev 13's model, §6.3);
* rev 13's "all four Bullet lensing strengths within about 1σ": Clowe et al.'s κ values are
  lower bounds, as their paper states. The calibrated test is the strong-lensing mass, and by
  that measure rev 13's smaller half was half as heavy as observed;
* rev 13's lighter star masses (mass-to-light 1–1.5): never required;
* rev 14's "the smaller half was a third of the main before the crash": the galaxies and
  starlight that would need are not seen (rev 15);
* rev 16 revisits two earlier statements:
  * "the Solar System is silent" still holds for the planets' own pull. But Cassini's test of
    the Galaxy's field is failed as the law stands (§6.10);
  * the wide-binary 19% becomes 1–5% if the Cassini fix is adopted (rev 17: 4% at 7,000 AU,
    9% at 20,000 AU; rev 18 adopts it);
* rev 16–17's "the dwarf fix must depend on speed": not forced any more, since the release
  length protects the binaries (§6.13);
* rev 19's "the far clusters need about 1.4 times their stars" and "their stars may be older than
  the Big Bang allows": the first came from mixing conventions, and distant cluster stars are
  observed to be younger (rev 20, §6.15);
* rev 19's "our distances make the strong lenses' stars 1.48 times lighter": 1.25 (the lens code
  itself was right);
* rounds 1–10's X-COP star masses, used as masses inside spheres although the release gives them
  inside circles on the sky (rev 20);
* earlier retractions:
  * a cluster claim (revision 5, retracted in revision 7);
  * a misattributed group-lensing figure (revision 8).

## 11. What is borrowed and what is ours

**Borrowed, and credited:**
* Newton's and Einstein's gravity in strong fields.
* Gauss's flux geometry.
* Dicke narrowing (Dicke 1953) and the Mössbauer effect, as known physics that we apply to
  the companion.
* The mathematical form of Milgrom's QUMOND field equation.
* The data:
  * SPARC (Lelli, McGaugh & Schombert);
  * X-COP (Eckert, Ettori, Ghirardini, Ghizzardi and collaborators);
  * Clowe et al. 2006;
  * Barrena et al. 2002;
  * Harvey et al. 2015;
  * KiDS-1000 (Brouwer et al. 2021; data release) and Mistele et al. 2024;
  * SLACS (Bolton, Auger, Treu and collaborators);
  * the Milky Way: McMillan 2017; Gaia rotation curves by Eilers et al. 2019, Zhou et al. 2023,
    Ou et al. 2024 and Jiao et al. 2023; Bovy & Rix 2013; Holmberg & Flynn 2004; McKee et al.
    2015; Posti & Helmi 2019; Watkins et al. 2019; Vasiliev 2019; Deason et al. 2019, 2021;
    Correa Magnus & Vasiliev 2022; Piffl et al. 2014; Monari et al. 2018; Koppelman & Helmi
    2021; Necib & Lin 2022; Prudil et al. 2022; Roche et al. 2024; Wegg, Gerhard & Portail 2016;
  * dwarf galaxies: McConnachie 2012; Caldwell et al. 2017; Torrealba et al. 2016, 2019;
  * precision tests: GRAVITY Collaboration 2020, 2024; Weisberg & Huang 2016; Kramer et al.
    2021; Bertotti et al. 2003; Hees et al. 2014;
  * lenses: van de Ven et al. 2010; Trott et al. 2010; Mróz et al. 2019;
  * star masses and ages (rev 20): Ghizzardi et al. 2021 and van der Burg et al. 2015 (X-COP's
    stars and their deprojection); Bradač et al. 2008 and Drory et al. 2004 (MACS J0025);
    Menanteau et al. 2012 and Hilton et al. 2013 (El Gordo); Holden et al. 2010, Moresco et al.
    2012, van Dokkum & van der Marel 2007 and Saglia et al. 2010 (distant stars are younger);
    Treu et al. 2010, Auger et al. 2010, Cappellari et al. 2013 and Conroy & van Dokkum 2012 (the
    IMF of giant ellipticals).

**Ours, as far as we have found:**
* the companion mechanism;
* heat as loss of step;
* collisions switching it off;
* the pull along the companion's net flow;
* the companion's memory after collisions;
* the companion breaking free over a release length;
* MOND's constant and switch derived from it.

A full literature search is still owed before any claim of priority.

## 12. Reproduce it

All scripts are in `research_work/results/hot-companion/code/`. Each writes to a fresh
output folder.

```
python dicke_toy.py      --output-dir ../run-dicke          # collisions switch scrambling off
python run_v3.py         --output-dir ../run-v3             # galaxies, clusters, constants, guard
python bullet_v3.py      --output-dir ../run-bullet-v3      # the Bullet Cluster
python collisions_v3.py  --output-dir ../run-collisions-v3  # 20 simulated collision pieces
python bullet_v4.py      --output-dir ../run-bullet-v4      # the Bullet with the companion's memory
python collisions_v4.py  --output-dir ../run-collisions-v4  # the collision stack with memory
python bullet_speeds_v4.py --output-dir ../run-bullet-speeds-v4  # galaxy speeds, as measured
python stream_tidal_v4.py  --output-dir ../run-stream-tidal-v4   # stream heat and tidal shaking
python bullet_main_v5.py   --output-dir ../run-bullet-main-v5    # main cluster's speeds and lensing masses
python bullet_main_v5.py   --output-dir ../run-bullet-main-v5-wide --dx 25   # the same, 4.8 Mpc box
python bullet_members_v6.py --output-dir ../run-bullet-members-v6  # lost galaxies in 78 velocities
python bullet_light_v6.py   --output-dir ../run-bullet-light-v6    # Legacy Survey star count
python wide_binaries_v6.py  --output-dir ../run-wide-binaries-v6   # wide binary stars near the Sun
python milky_way_v7.py      --output-dir ../run-milky-way-v7        # the Milky Way against Gaia and more
python mw_dwarfs_v7.py      --output-dir ../run-mw-dwarfs-v7        # ten dwarf galaxies
python strong_field_v7.py   --output-dir ../run-strong-field-v7     # planets, S2, pulsars, Cassini
python lensing_census_v7.py --output-dir ../run-lensing-census-v7   # KiDS, Mistele, Einstein Cross
python transition_check_v7.py --output-dir ../run-transition-check-v7  # a gentler hold (larger g_d), on SPARC
python kids_v3.py        --output-dir ../run-kids-v3        # ellipticals vs spirals
python lenses_t35.py     --output-dir ../run-lenses-v3 --constants ../run-v3/results.json
python derive_mond.py    --output-dir ../run-derive         # MOND as the cold limit
python field_equation.py --output-dir ../run-field          # the 3D field equation
python collisions_v8.py  --output-dir ../run-collisions-v8  # MACS J0025, Abell 520, El Gordo
python release_hold_v9.py --output-dir ../run-release-hold-v9  # release length and the Milky Way's hold
python first_principles_v10.py --output-dir ../run-first-principles-v10  # the pull of a locked emitter
python combination_rules_v10.py --output-dir ../run-combination-rules-v10  # how the companion adds up
python collisions_v10.py --output-dir ../run-collisions-v10    # the far collisions in our own distances
python collisions_star_sweep_v10.py --output ../run-collisions-v10/star_sweep_v10.json  # how heavy their stars must be
python forecast_wide_binaries_v10.py --output ../forecasts/wide_binaries_gaia_dr4_v10.json  # the locked forecast
python companion_flow_v11.py --output-dir ../run-companion-flow-v11    # energy bookkeeping and the three dials
python xcop_static_v11.py   --output-dir ../run-xcop-static-v11        # X-COP deprojected, in our distances; the constants
python kids_static_v11.py   --output-dir ../run-kids-static-v11        # galaxy lensing in our distances; the distance scale
python bullet_static_v11.py --output-dir ../run-bullet-static-v11      # the Bullet Cluster in our distances
python distance_variants_v11.py --output-dir ../run-distance-variants-v11  # the two static geometries
python macs_peak_scan_v11.py --output-dir ../run-collisions-v10        # MACS J0025's NW peak by age and star basis
python macs_timing_v12.py   --output-dir ../run-collisions-v10         # MACS J0025's age from its shock fronts
python lens_heat_sdss_v12.py --output ../data/lens_heat_sdss_v12.json  # the lenses' star speeds, from SDSS
python kids_heat_v12.py     --output-dir ../run-kids-heat-v12          # galaxy lensing with the measured heat
python sn_scale_v12.py      --output ../run-distance-scale-v12/sn_scale_v12.json      # the distance rate from supernovae
python distance_scale_v12.py --output-dir ../run-distance-scale-v12    # the rate fitted jointly (--adopt 0.95: the constants)
python kids_level_v12.py    --output ../run-distance-scale-v12/kids_level_v12.json    # what sets the lensing level
python companion_toy_v12.py --output ../run-companion-toy-v12/companion_toy_v12.json  # five local rules for the companion
python companion_toy_v12.py --convergence --output ../run-companion-toy-v12/convergence_v12.json  # the same at three resolutions
python coherent_force_v13.py --output ../run-coherent-force-v13/coherent_force_v13.json  # the pull and the scrambling in one experiment
python joint_checks_v13.py  --output ../run-coherent-force-v13/joint_checks_v13.json    # 16 timings, doubling the speed, free amplitudes
python strength_offset_v13.py --output ../run-coherent-force-v13/strength_offset_v13.json  # feed quiet waves, absorb loud ones
python heat_exponent_v13.py --output ../run-coherent-force-v13/heat_exponent_v13.json  # how steeply the heat term must grow
python singer_v14.py        --output ../run-coherent-force-v13/singer_v14.json          # a singer that absorbs; a singer with a power balance
python release_shape_v14.py --output ../run-coherent-force-v13/release_shape_v14.json   # the release factor against the absorbers' hold, on SPARC
python release_depth_v14.py --output ../run-coherent-force-v13/release_depth_v14.json   # how far the release goes, in three densities
python wave_dark_v15.py     --output ../run-reservoir-force-v15/wave_dark_v15.json      # dark states of wave-coupled emitters, moved
python receivers_v15.py     --output ../run-reservoir-force-v15/receivers_v15.json      # which bodies a passing wave pulls
python reservoir_force_v15.py --output ../run-reservoir-force-v15/reservoir_force_v15.json  # the quiet store, real waves, the pull
```

The independent calculation of rev 24 reruns from its own folder (about a minute; see `independent-r15/REPRODUCTION.md`):

```
cd ../independent-r15/first_principles_test
python finite_reservoir_test.py
python test_dark_bright.py --out results --part static    # also: collisions, controls, no-sink
```

The regression suite runs everything at once and compares with the saved baseline:

```
cd ../regression
python run_suite.py                        # the adopted law (rev 21), quick tier, about a minute
python run_suite.py --tier full            # everything, about 20 minutes
python run_suite.py --law no_hold          # a candidate change, scored against the baseline
python run_suite.py --law round11          # rev 20's constants and distance rate
python run_suite.py --law round9           # revs 18-19's constants
python run_suite.py --law round3           # the law before rev 18
python run_suite.py --law heat_p175 --tier full   # rev 23: the heat weight as sigma^1.75, a and u refitted
```
