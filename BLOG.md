# Gravity that streams

*One proposed law for spinning galaxies, bending light, galaxy clusters and colliding clusters: a few stated
assumptions, worked out step by step and tested on public data, with no dark matter and no expanding universe.*

24 September 2026, updated 25 September 2026 with round 19 (the companion's medium tested against the data, the
Milky Way refitted with its matter held to independent measurements, a frozen prediction) and a step back over the
whole record ([STEP-BACK-AUDIT.md](STEP-BACK-AUDIT.md), §9), round 20 (a second way to derive the pull, tested in
small models, and the clusters checked without their X-ray input), and round 21 (three proposals from the project's
owner checked against the data: a faster glow from colliding clusters, a refined distance law, and a first
relativistic version, §5.3 and §9), and round 22 (the faster glow taken as far as the model allows: how it must spread,
where its energy comes from, and what the companion's medium must do, §5.3 and §9). This is a fresh write-up
of where the project stands. The full working notebook, with every
step, revision and correction along the way, is archived in
[research_work/blog-archive/BLOG-notebook-rev26.md](research_work/blog-archive/BLOG-notebook-rev26.md), and the
technical record is in [research_work/results/hot-companion/README.md](research_work/results/hot-companion/README.md).
Every number below is computed from public data by a script in this repository (§10).

**Contents**
1. [The puzzle](#1-the-puzzle)
2. [How we work](#2-how-we-work)
3. [The idea, step by step](#3-the-idea-step-by-step)
4. [How the physics works inside matter](#4-how-the-physics-works-inside-matter)
5. [The data](#5-the-data)
6. [Why stars and light lens the way we see](#6-why-stars-and-light-lens-the-way-we-see)
7. [How it compares](#7-how-it-compares)
8. [Predictions anyone can check](#8-predictions-anyone-can-check)
9. [What is still open, and why we are optimistic](#9-what-is-still-open-and-why-we-are-optimistic)
10. [Reproduce it](#10-reproduce-it)

---

## In one minute

**The puzzle.** Galaxies and clusters of galaxies pull harder than the matter we can see should allow. Stars at a
galaxy's edge orbit too fast, light passing a galaxy bends too much, and clusters hold their hot gas with 4 to 27
times the pull their visible matter provides. When two clusters collide, the extra pull follows the galaxies, not the
gas, even though the gas outweighs the galaxies ten to fifty times. The usual answer is invisible "dark matter".

**Our answer.** Every piece of ordinary matter feeds a faint companion to its gravity, which streams outward at
169 km/s and pulls in proportion to its strength. Five plain facts about how that companion behaves explain the four
puzzles together:

1. **Orderly matter feeds it in step,** so opposite contributions cancel the way Newton's pulls do. What survives falls
   off slowly enough to keep a galaxy's outer stars orbiting fast.
2. **Matter whose parts move randomly and freely** (stars in elliptical galaxies, galaxies in clusters) **feeds it much
   harder.** That is why clusters need so much extra pull, and why ellipticals bend light more than spirals.
3. **Gas doesn't count as moving randomly,** because its particles collide far too often. So in a collision the extra
   pull stays with the galaxies and leaves the gas behind.
4. **Strong gravity holds the companion back,** so the Solar System feels nothing.
5. **The companion is slow and keeps moving the way its source moved,** so after two clusters collide, each is still
   wrapped in the companion it had before.

**The scoreboard,** with four constants shared by everything (three fitted to galaxies and clusters, one set by the
Solar System):

| Test | Our law | MOND | Dark matter |
|---|---|---|---|
| 149 galaxies' rotation speeds (typical miss) | **15.9 km/s** | 16.1 km/s | 7.5 km/s, with 298 adjustable numbers |
| 12 galaxy clusters' masses (typical miss) | **25%** | ×2.9 | 11%, with 24 adjustable numbers |
| Colliding clusters: lensing sits on the galaxies, not the gas | **yes** (Bullet Cluster, 72-collision stack, three more) | no | yes |
| Ellipticals bend light more than spirals (measured 0.15 ± 0.04 dex) | **yes: 0.13–0.16**, from their stars' random motion | no difference | yes, with tuned haloes |
| Light and matter feel the same pull (six strong lenses) | **yes** (−0.03 ± 0.02 dex), with light's response assumed as in Einstein's theory | | yes |
| Milky Way rotation from 15 to 27 kpc | **within 1–6%**; with the Galaxy's matter held to its own measurements (2026 Cepheids included) the curves' shape is followed, 3–6% slow overall (§5.4) | within 3% | 6–11% too fast |
| Solar System, planets, pulsars | **no measurable extra pull** (switched off where gravity is strong) | small effects | Einstein's |
| Cassini's limit on the Galaxy's distortion of the Sun's field (2026 re-analysis) | **within 0.2–0.4σ** with the Galaxy's pull at the Sun from the refitted Galaxy (1.5σ with a rough estimate of it) | common form about 20× too big | Einstein's |

These are typical misses and simple checks, not a full statistical comparison: each model here uses its own treatment
of the uncertain inputs every model needs (star masses, gas, distances), and "four constants" does not count those
inputs. A proper comparison, with shared inputs, their uncertainties and fair baselines, is still to be done (§9).

---

## 1. The puzzle

Four measurements have to be explained together.

1. **Rotation.** In a spiral galaxy, stars far from the centre should orbit slowly, the way Neptune orbits the Sun
   more slowly than Mercury. They don't: the speeds stay flat. Across 149 well-measured galaxies (the SPARC catalogue,
   3,150 measured speeds), Newton's law applied to the visible stars and gas misses by 45.6 km/s on average.
2. **Light bending (lensing).** Galaxies and clusters bend the light of objects behind them more than their visible
   matter can.
3. **Clusters.** In 12 well-measured galaxy clusters (the X-COP sample), the pull needed to hold the hot gas in place is
   up to 27 times what the visible matter supplies near the centre, and 4 times at the edge.
4. **Collisions.** In the Bullet Cluster, two clusters passed through each other about 150 million years ago. Their
   galaxies sailed through; their gas, which outweighs the galaxies ten to fifty times, collided and was left behind in
   the middle. The lensing, which traces the pull, sits on the galaxies. A stack of 72 such collisions says the same.

The three standard responses each fail somewhere:
* **Dark matter** adds invisible mass, tuned object by object until the sums work. It fits well, but needs about two
  adjustable numbers per galaxy or cluster, and the particle has never been found.
* **MOND** changes the law of gravity at low accelerations. It predicts galaxy rotation well with one constant, but
  it misses cluster masses by a factor of about 3 and cannot put a collision's lensing on the galaxies.
* **Newton with visible matter only** fails everywhere outside the Solar System.

## 2. How we work

The project runs under four rules ([RULES.md](RULES.md)): no dark matter; no MOND or anything derived from it as an
input; no plain Newtonian gravity as the answer; and no expanding universe. Distances are the project's own, from a
static universe in which light loses energy slowly as it travels (its rate is fitted to 1,365 supernovae). A
refinement of that law, and what it requires, was tested in round 21 (§9).

Every candidate formula goes through an automatic check,
[research_work/tools/formula_guard.py](research_work/tools/formula_guard.py), that asks whether it is secretly MOND
(a function of the Newtonian pull alone), secretly Newton (a constant boost) or secretly dark matter (extra, adjustable,
invisible mass).

The order of work is always the same: first principles, then a small model that can be solved exactly, then real
data. A regression suite reruns every test the law has faced (77 graded checks) in one command, so any change shows at
once what it fixes and what it breaks. It is an engineering tool: a "pass" means within two standard deviations of one
measurement, the checks are not all independent, and the tally is not a measure of how likely the law is to be right.

## 3. The idea, step by step

### 3.1 Matter feeds a companion to its gravity

Suppose every kilogram of ordinary matter steadily feeds a tiny amount of energy, ℓ watts per kilogram, into a
companion to its gravity, and that this energy streams outward at a speed u. That is the one new ingredient. The value
that fits the data is ℓ = 5.3 × 10⁻⁶ watts per kilogram: over 13 billion years a kilogram gives off two and a half
thousandths of one percent of its rest energy this way. For the Sun it is 2.8% of its light, a mass loss of
2 parts in 10¹⁵ a year, too little to notice.

### 3.2 Orderly matter adds up the way Newton does

Energy flowing away from a source spreads over ever-larger spheres, so its flux falls as 1/distance², the same
geometry as Newton's pull. Add up the flows from many sources that feed in step, arrow by arrow, and the total points
exactly along Newton's pull g_N of all the visible matter (this is Gauss's geometry; our script checks it to 2 parts in
10¹⁵). Opposite contributions cancel. We also derived this from energy conservation, and a small working model shows
that a companion guided along gravity's field lines keeps every watt and forms no whirlpools.

### 3.3 The companion pulls with its strength

Give the companion an energy density A²/8πG, the familiar form of the energy stored in a gravitational field, and let
its pull equal its strength A. Energy balance through a sphere around a mass M says that what streams out equals what
is fed in:

```
4πr² · u · A²/8πG = ℓ M    ⟹    A = √(G M a) / r = √(a · g_N),    a = 2ℓ/u
```

This pull falls as 1/r, not 1/r². A star orbiting at distance r then has v² = r·A = √(G M a), **the same at every
radius: a flat rotation curve.** It also gives the observed rule v⁴ = G M a, linking a galaxy's visible mass to its
rotation speed. The constant a = 6.30 × 10⁻¹¹ m/s² is not new: it is twice the power per kilogram divided by the
companion's speed. §4 shows where "the pull equals the strength" comes from inside matter.

### 3.4 Strong gravity holds the companion back

Near the Sun the pull is billions of times stronger than in a galaxy's outskirts, and planetary orbits show no
companion at all. So the companion must be held back where ordinary gravity is strong and released where it is weak:

```
released share = exp( −|g_N| / g_d ),    g_d = 2.03 × 10⁻¹⁰ m/s²
```

At the Earth the pull is 28 million times g_d, so the released share is exactly zero. In a galaxy's centre about 1% is
released; in its outskirts and in clusters, most of it. Breaking free also takes time: the companion must travel about
0.15 parsec (about 900 years at its speed) from its source before it is released, which Cassini's radio tracking of
Saturn requires.

### 3.5 Random motion makes more companion

The companion streams out of all matter at u = 169 km/s. A piece of matter sitting still therefore still feels its
companion rushing past at 169 km/s, and that relative motion is what opens its internal store: it is the ordinary
("cold") glow. A piece moving at random with speed spread σ feels the stream go by faster, on average by
u² + 3σ² instead of u², the 3 coming from the three directions it can move in. So its glow grows by the **heat
weight**

```
k = 3σ²/u²
```

and the extra adds up as a plain total S, because the random parts don't cancel as arrows do. Inside a big,
spread-out cloud, Newton's pulls from all sides largely cancel toward the centre, but S does not. In a cluster that is
exactly the extra pull needed: large near the centre, smaller toward the edge. With one speed, u:

| System | Random speed of its stars or galaxies, σ | Heat weight k |
|---|---:|---:|
| Disk of a spiral galaxy | 10–30 km/s | 0.01–0.09 (cold) |
| Elliptical galaxy | 150–250 km/s | 2.4–6.5 (warm) |
| Galaxies in a cluster | 500–1,200 km/s | 26–150 (hot) |

A spinning disk does not count as hot, because the companion moves along with the matter that made it (§3.8).

### 3.6 Collisions switch the heat off

Laboratory physics has a well-known exception to motion's effects on emission: if an emitter changes direction many
times before its internal response can follow, the effect averages away (Dicke narrowing, 1953; its cousin, the
Mössbauer effect, is why atoms locked in a crystal emit perfectly sharp lines). Gas particles collide constantly, and in
clusters they also spiral around magnetic fields. Stars never collide; galaxies pass through each other. So:

| | Moves freely? | Companion |
|---|---|---|
| Stars in an elliptical galaxy; galaxies in a cluster | yes | **hot**: extra glow |
| Gas and plasma, including the Sun's interior | no, constant collisions | **cold**: no extra glow |

Gas still counts in full as ordinary mass.

### 3.7 The companion pulls along its net flow

The extra pull points along the net flow of companion energy: Newton's direction plus the flow from the hot matter,

```
direction = (g_N + g_hot) / (|g_N| + |g_hot|),    g_hot = the Newtonian pull of the free-moving matter, weighted by k
```

In a round, settled galaxy or cluster both point to the centre and nothing changes. In a collision the cold gas
dominates Newton's pull but the hot galaxies dominate the companion's flow, so the pull, and the lensing, bend toward
the galaxies.

### 3.8 The companion carries momentum, is slow, and remembers

A streaming field carries momentum, as light does: a lamp does not recoil when its light later pushes on something far
away. So gas and galaxies feel the same pull, as observed, and the stream takes up any recoil.

At 169 km/s the companion takes about 580 million years to cross 100 kpc, so most of the companion around a galaxy
today was given off long ago. And it keeps moving the way its source moved, like a ball thrown from a moving train.
It has to: galaxies move through space at hundreds of km/s, and a galaxy's gravity cannot depend on how fast it
happens to be moving. For settled systems this changes nothing. After a collision it matters a lot: the galaxies sail
through and keep their motion, so the old companion travels on with them, while around the stopped gas a fresh
companion has had time to grow only about 30 kpc.

### 3.9 The whole law

```
g_N   = Newton's pull of all ordinary matter (stars, gas, everything)
S     = G ∫ k ρ_free / d²                    the plain total from the free-moving matter
g_hot = G ∫ k ρ_free (x′ − x) / |x′ − x|³    its heat-weighted flow
k     = 3σ²/u²  for stars and galaxies,   0 for gas and plasma

extra pull = exp(−|g_N|/g_d) · √( a (|g_N| + S) )   along   (g_N + g_hot) / (|g_N| + |g_hot|)
h = g_N + extra pull;    ∇²Φ = −∇·h;    the pull on matter and on light is  g = −∇Φ
```

For round, settled systems it is simply **g = g_N + exp(−|g_N|/g_d) · √(a (|g_N| + S))**. The field equation in the
last line turns the rule into one landscape Φ that everything rolls on, light included (§6).

| Constant | Value | Measured on |
|---|---|---|
| a = 2ℓ/u | 6.30 × 10⁻¹¹ m/s² | 149 galaxies |
| g_d, the release level | 2.03 × 10⁻¹⁰ m/s² | 149 galaxies |
| u, the companion's speed | 169 km/s | 12 clusters |
| L, the release length | 0.15 pc | Cassini, in the Solar System |

None is set object by object. Each piece of the law is a switch set by something measurable about the matter itself:

| Switch | Set by | Separates |
|---|---|---|
| Release factor | how strong the pull is | the Solar System from galaxies |
| Heat weight | how randomly the matter moves | disks, ellipticals and clusters |
| Collision rule | whether its particles collide | gas from stars |
| Direction | which way the companion's energy flows | round systems from collisions |
| Memory | how recently the matter changed its motion | settled systems from collisions |

**MOND comes out as a special case.** Switch the heat off (cold disk galaxies) and the law becomes
g = g_N · ν(g_N/a) with ν(y) = 1 + e^(−y/3.22)/√y, a MOND-type formula. But here MOND's constant is 2ℓ/u rather than
a new constant of nature, its switch between Newton and the flat regime comes from the release factor, and the law says
why disks obey it (they are cold) and why clusters and collisions don't (their galaxies are hot, their gas is not).
Everything we claim as ours lies where the law is not MOND.

### 3.10 What is assumed, what is derived, and what is still open

The law rests on a few assumptions. Some consequences follow from them exactly, some have been derived in small
working models (§4), and some pieces are still assumptions. Keeping them apart:

| Piece of the law | Status |
|---|---|
| Matter feeds a companion at ℓ watts per kilogram, which streams outward at u | **assumed** (ℓ and u fitted to galaxies and clusters) |
| Its energy density is A²/8πG and it pulls with its strength A | **assumed**; §4.1 shows a kind of matter (inverted, self-sustained) that a wave pulls in proportion to its height |
| Flat rotation curves, v⁴ = G M a, and contributions adding up along Newton's direction | **derived** from the two lines above (energy balance and Gauss's geometry); the working models now reproduce the square root of mass and the fall with distance, with uncertainties (§4.4) |
| Strong gravity holds the companion back, exp(−\|g_N\|/g_d), released over L = 0.15 pc | **assumed** (g_d fitted to galaxies, L set by Cassini); a candidate process now gives both factors from one lifetime (§5.5), not yet found in the matter model |
| The heat weight k = 3σ²/u², with the same u | **derived** in a working model (§4.2), given a quiet store inside matter that the flowing companion opens; the store itself is assumed |
| Collisions switch the heat off | **derived** in the same working model; the physics is borrowed (Dicke narrowing) |
| The pull points along the companion's net flow; the companion remembers its source's motion | **assumed** (motivated, not derived) |
| The field equation ∇²Φ = −∇·h | **borrowed** form (Milgrom's QUMOND) |
| Light responds to Φ as in Einstein's theory | **assumed**; to be derived in a relativistic version (a first proposed action was checked in round 21: what it gets right and four fixes, §9) |
| Warm matter pulls harder, by the square root of its extra glow | **required by the data** (clusters, the lensing of ellipticals, collisions); shown in the working models only with strongly one-way waves, which the clusters now rule out (§4.4). The mechanism is **open** again |
| What a piece's inertia is, and why all matter falls alike | **open**: the working models do not yet say |
| MOND as the cold limit | **derived** from the law |

So "from first principles" means here: a few stated assumptions, followed exactly wherever we can, and tested at each
step. Several links are still assumptions, and the working models that support others are simplified.

## 4. How the physics works inside matter

Section 3 is the law. This section is about the machinery underneath it: what kind of matter, and what kind of wave,
make a companion that pulls with its strength, glows more when it moves, and holds its beat. Here we work with small
models solved exactly, with every watt of energy accounted for, and we say plainly what is derived and what is still
assumed.

### 4.1 Matter that holds energy ready to give is pulled by a passing wave

Picture each piece of matter as holding a quiet internal oscillation, like the atoms of a laser that hold energy ready
to give. Such a piece is "inverted": it sustains itself, and when a companion wave passes, it falls into step with it
by itself. The question is how.

![Which bodies a passing wave pulls](blog-figures/which-bodies-a-wave-pulls.png)

*One body in a steady wave, from its own equations. An ordinary self-sustained oscillator falls a quarter beat behind
the wave, takes energy from it, and is pushed. An inverted one falls a quarter beat ahead, feeds energy into the wave,
and is pulled, in proportion to the wave's height. Script `code/receivers_v15.py`.*

* **The sign comes from the body's own oscillation.** An inverted self-sustained body locks a quarter beat ahead of a
  passing wave by itself, feeds energy into it, and is pulled toward its source. The pull equals the power it feeds
  divided by the wave's speed, to four decimal places.
* **It is pulled in proportion to the wave's height, not its energy**: ten times the height gives 9.3 times the pull
  in weak waves. The height of the companion is the square root of its power, so this is where the law's square root
  comes from.
* A body that only amplifies, below its own threshold, is pulled in proportion to the wave's energy instead, which
  falls as 1/r² like Newton's pull. The law's flat rotation curves need self-sustained, inverted matter.

### 4.2 Where the heat comes from: the companion streaming past

Each piece of matter holds a quiet store that the companion flowing past opens a little, and random motion opens more.

![Motion measured against the flow](blog-figures/heat-from-the-flowing-companion.png)

*Pieces whose quiet store opens to their speed relative to the companion streaming past. The glow relative to rest lies
on the law's 1 + 3σ²/u² (dashed): 1.70, 3.88, 12.7, 48.6 and 198 at σ/u = 0.5, 1, 2, 4, 8, against 1.75, 4, 13, 49
and 193. Script `code/stream_store_v16.py`.*

* The heat weight comes out exactly, factor 3 and all, with the **same** 169 km/s that sets how fast the companion
  travels. Our law has always used one speed for both jobs; here that is automatic.
* **Collisions hold it back** as the law needs: at collision rates of 1, 3, 10 and 30 (in units of the store's
  response rate) the extra glow is 5.90, 2.96, 1.07 and 0.38, against 6, 3, 1.09 and 0.39 expected.
* A new prediction falls out: matter falling inward, against the outflow, glows more.

### 4.3 One kind of matter, every piece both sender and receiver

The real test is to build sources and receivers out of the same matter, let every piece talk to every other through
the complete wave (near field and far field), and account for every watt. A ball of 48 pieces is the source; more
pieces of exactly the same matter sit farther out as receivers. Nothing sets anyone's timing, and nothing converts
power into force by hand.

![One kind of matter: a warming source](blog-figures/one-kind-of-matter.png)

*As the source warms, the wave reaching distant receivers grows as the square root of the source's glow (dashed), as
the law needs. With the wave travelling both ways, though, the pull on the receivers does not follow (§4.4). Script
`code/one_matter_v16.py`.*

* **The books balance:** energy to about 1 part in a trillion, and momentum as far as it can be measured.
* **A cold source pulls** distant pieces of the same matter: they settle a quarter beat ahead of its wave and are
  pulled, reaching 0.92 on a scale where 1 is perfect step.
* **A warming source glows more, but by less than the heat rule's 1 + k:** ×1.2, 2.1, 5.0 and 8.1 of its cold glow at
  heat weights 0.5, 2, 8 and 16, where 1 + k gives 1.5, 3, 9 and 17. Two things hold it back. Each piece's extra glow
  draws on its own supply: pieces far apart reach only 76–87% of 1 + k (×2.6, 7.0 and 11.8 at k = 2, 8 and 16), as
  their inner reserve falls by up to a quarter. And in a dense source, whose pieces sit closer than a wavelength, the
  shared wave holds back part of the rest (45–60% of 1 + k). In the law this would read as a larger effective u for
  densely packed matter; whether that is the same for all real systems is an open question (§9). The wave reaching
  distant matter grows as the square root of whatever glow there is (×1.01, 1.31, 2.24, 2.81): the square root the law
  needs is physically there in the wave.
* **An exact rule about matter:** a piece's own loud inner vibrations, stirred by the same wave, kick back and take a
  share 2f of its pull, where f is the share of their energy they send into the companion. So matter's loud vibrations
  must mostly ring inside and only whisper into the companion.

### 4.4 Keeping warm matter in tune

In the version of §4.3, a warm source's pieces fall out of tune with each other, distant matter can't follow its wave,
and the pull doesn't grow with the heat. We found exactly why, and a way out.

![Which part of the warm pushes does the damage](blog-figures/what-does-the-damage.png)

*Distant matter keeping step with a source at a heat weight of 8. Script `code/rhythm_protect_v17.py`,
`code/one_way_v17.py`.*

* **What does the damage.** Warm pieces push on each other's beats in two ways: one that always lets a group settle
  into a common rhythm, and a sideways one that can keep it churning forever. Switch off only the sideways half and
  distant matter keeps step at 0.88 with a warm source, better than the 0.63 of a cold one, because the warm source's
  wave is stronger.
* **None of the twelve inner structures we tried removes it.** At equal glow the pushes came out the same size in all
  of them, within 4%. A basic rule of waves (reciprocity) suggests why, for any structure that sends and receives
  through the same channel: whatever lets a piece send its glow into the wave lets it receive its neighbours' glow just
  as strongly. What does help: motion must enter the equations the way a velocity does (it flips sign when time runs
  backwards, like the Coriolis force), which is a correction from first principles.
* **A wave that only travels outward removes the damage completely.** The companion streams outward from the matter
  that makes it. If its crests are carried outward by that stream faster than they can move against it, a piece hears
  only matter nearer the centre, and never its own echo.

![Keeping step at a distance](blog-figures/keeping-step-at-a-distance.png)

*With a two-way wave (purple) distant matter loses a warm source; with a wave that only travels outward (green) it keeps
step at every distance, better than with a cold source (dashed). Script `code/one_way_v17.py`.*

* With the outward-only wave, a warm source keeps one beat (70 times tighter than with the two-way wave) and distant
  matter keeps step at 0.98, 0.96, 0.85 and 0.81 at four distances. What matters is that the wave has no way back: any
  one-way order works, and an outflow supplies the natural one, from the centre outward.
* **In the full simulation, warm matter then pulls distant matter harder than cold matter**, as the law needs. The net
  pull on receivers of the same matter at four distances (units of 10⁻⁴; a minus sign is a push):

  | Source | r = 6 | 9 | 13.5 | 20 |
  |---|---:|---:|---:|---:|
  | Cold, outward-only wave | +1.61 | +1.14 | +0.41 | +0.53 |
  | **Warm (heat weight 8), outward-only wave** | **+2.36** | **+2.12** | **+1.75** | **+1.07** |
  | Warm but colliding, outward-only wave | +1.53 | +0.99 | +0.38 | +0.50 |
  | Warm, wave travelling both ways | −0.89 | −0.55 | −0.03 | −0.09 |
  | Cold, stream that absorbs inward waves (round 18, κ = 5) | +1.36 | +0.76 | +0.56 | +0.47 |
  | Warm, stream that absorbs inward waves (round 18, κ = 5) | +1.72 | +1.38 | +1.23 | +0.60 |

  Measured against the square root of the wave's extra strength, the warm source's extra pull comes out 1.06 on
  average, close to the 1 the law's square root requires; distance by distance it ranges from 0.6 to 1.8, because the
  cold source that serves as the reference keeps step unevenly. Colliding sources pull like cold ones, and a warm source brought
  to rest goes back to about the cold pull. With the wave travelling both ways, the same warm source pushes the nearest
  matter away instead: the outward-only wave is what turns the push into a pull.
* **The glow stays in proportion to mass, and the pull grows about as its square root, as the law needs.** In the full
  simulation the glow per piece is the same within 8% for sources of 24, 48 and 96 pieces, and four times the mass
  gives ×1.7 the pull for a cold source and ×2.1 for a warm one, where the law's √(G M a)/r gives ×2 (two arrangements
  each; a proper fit of the exponent, with its uncertainty, over a wider range is still to be done). With the rule
  imposed in these tests every piece hears every piece nearer the centre, which for a round source is Newton's rule
  that only the mass inside a radius pulls there; round 18 found that the physical one-way media hear less (below).
* **Status, and what round 18 found.** The one-way wave was imposed in these tests, by deleting every inward coupling
  by hand. Round 18 gave the wave its own local equations, first in one dimension, solved exactly and checked against a
  direct simulation of the medium (forces and powers agree within 1–4%):
  * **Waves carried by the companion's own stream are one-way whenever the stream outruns them** (the stream's speed u
    above the waves' speed c relative to it). This is now derived, not imposed, and for a stream that moves at u
    everywhere, as ours does, it holds right into the centre of a source. With both of the stream's waves included, the
    pull's energy bill fits the galaxies' limit (§4.5) when c is at least 71% of u.
  * **But a stream flowing past a body pushes on anything that emits into it,** like wind on a sail: each piece is
    pushed downstream by about its emitted power divided by the stream's speed. For the companion that push is as
    large as the law's own pull in the outskirts of galaxies, so this version fails as it stands.
  * **A version with no push:** waves that travel at their own speed through matter's frame, with the stream absorbing
    the part that moves against it. Emission is then symmetric (no push), the coupling is one-way (a wave loses a
    factor e^(−κ) for every unit of distance it travels inward), and the energy bill is the static one (crests at most
    u/2). In the full three-dimensional simulation (table above) it removes the push, colliding sources pull like cold
    ones, and a warm source pulls distant matter harder than a cold one at every distance, with every watt booked
    (the stream absorbs 33–48% of what the pieces give the wave) and the medium passive, unlike the rule imposed by
    hand, which could in principle create energy. How much harder: 73–84% of the square root of the model's own extra
    glow; but the model's glow grows only about half as fast as the law's, so against the law's heat gain the pull
    rises about 1.6 times where the law needs 3 (corrected 25 September; an earlier version of this page said "73–84%
    of the square-root growth the law needs").
  * **What it does not yet do** is keep a warm source on one beat. In it a piece hears only the inner matter on its
    own side of the source, because a wave that must cross the centre is absorbed; for strong absorption a distant
    receiver hears mainly the near half of a source. The imposed rule let every piece hear everything nearer the
    centre, which is what gave one beat and Newton's rule. So those two results belonged to the rule, not to the
    physical media found so far, and the law's counting of all matter (its |g_N| and S) now has to be reconciled
    with a medium that hears less: the review's hot-shell benchmark, in general form.
* **What round 19 found.** Three things, in order: the absorbing stream works as a real wave medium; the clusters say
  how one-way it may be; and with that allowed amount, the working models give the law's cold part but not its heat
  part.
  * **The absorbing stream as an exact wave medium.** Solved exactly (for a uniform stream, the geometry a distant
    receiver sees), it needs absorbers that move along the stream at the wave's own speed. A wave moving with the
    stream then doesn't oscillate at all as they see it, and nothing passive can absorb a wave that doesn't oscillate,
    so the medium is one-way with no tuning: the Doppler effect alone does it. The absorbers must be at least about a
    wavelength across; point-like ones would drag every emitter downstream, which rules them out. With that, the medium
    can never create energy, at any absorption strength, and it pulls a lone emitter only slightly toward the source's
    centre (at the weak absorption the data allow, below, by 1.5% of its emitted power divided by the wave speed,
    within what the 149 galaxies allow). Round 18's
    simple description, a loss of e^(−κ) per unit distance travelled inward, is accurate when the absorption is weak
    (an absorption length of a wavelength or more). A passive medium cannot be strongly one-way over only a few
    wavelengths, so round 18's strong-absorption runs overstated how one-way it was. *A catch found by the step back
    (§9):* those absorbers must move at the wave's speed, but to pay the pull's energy bill the companion's own stream
    must move at least twice as fast. Absorbers carried by the stream would amplify some waves instead of absorbing
    them, so they would have to be something separate from the stream, and what they could be is open.
  * **What the clusters allow (the hot-shell benchmark).** In a cluster, most of the law's heat term at a given radius
    comes from hot galaxies farther out: 69% at a tenth of the cluster's radius, 28% at its edge. A strictly one-way
    medium would never let that inner point hear them. Tested on the 12 X-COP clusters, the law, which hears all hot
    matter, misses their masses by a typical 25%. Hearing only the matter nearer the centre misses by 50%, still 32%
    with the companion's speed refitted, and then fails 11 of the regression suite's 49 checks against the law's 6.
    A stream that absorbs inward waves is tolerated only if a wave survives about 300 kpc or more before it is
    absorbed. That is longer than a galaxy and comparable to a cluster's radius. With stronger absorption the clusters
    can be rescued only by slowing the companion (to 141 km/s at 333 kpc, 96 km/s at 10 kpc), which makes hot stars
    louder in every galaxy and spoils the lensing of ellipticals and the strong lenses. For a round source there is
    also an exact result: a stream that absorbs every inward wave hears exactly half of every inner shell of matter,
    so the only real difference from the law is the outer shells.

![What the clusters allow](blog-figures/what-the-medium-hears.png)

*The typical miss in 12 clusters' masses against how far an inward-travelling wave survives before the stream absorbs
it. Orange: the constants as fitted. Teal: with the companion's speed refitted to the clusters (numbers: that speed in
km/s; the law's is 169), which the rest of the suite rejects below about 300 kpc. Dashed: a strictly one-way medium.
Script `code/hot_shell_v19.py`.*

  * **The law the working models produce.** With that weak absorption (one absorption length per wavelength, the most
    the clusters allow for sources a few wavelengths across), and with the plain two-way wave for comparison, we ran
    84 full, energy-balanced simulations: sources of 24 to 96 pieces, cold, warm and colliding, packed at different
    densities, with receivers of the same matter at two distances and three arrangements of each. For cold sources
    the pull grows as mass to the power 0.58 ± 0.13 and falls with distance as the power 1.11 ± 0.22 (0.80 ± 0.29
    and 0.93 ± 0.20 with the two-way wave), against the law's 0.5 and 1. Colliding sources pull 0.78–0.90 times as
    hard as cold ones, where the law says 1. **But warm sources do not pull harder:** 0.90–0.97 times the cold pull
    where the law needs 1.73 at the milder heat, and at the stronger heat they push the nearer receivers away. These
    are averages over three arrangements that scatter widely: with the plain two-way wave one arrangement did show the
    law's full gain (2.2 and 1.8 times the cold pull at the two distances) while the others pulled less or pushed, so
    three arrangements cannot settle it. And the models' pull depends on a source's size at fixed mass (sources with a radius of one or two wavelengths pull 2–15
    times less than those of three or four), which real galaxies don't show. Every watt is booked in every run.
  * **What this means.** The law's cold part now comes out of the working models with honest error bars, and so does
    the collision rule. The heat term is still required by the data, which make their case on their own (the
    clusters' masses, the extra lensing of ellipticals, the pull following galaxies in collisions), but its mechanism
    is open again: the strong one-way absorption that made warm sources pull harder in rounds 17–18 is the kind the
    clusters exclude. The next mechanism has two clear targets: the hot glow must reach inward, and a warm source must
    keep its beat without strong one-way absorption.

### 4.5 What being pulled costs

A body pulled by feeding a wave pays for it: the price is the force times the speed of the wave's crests. The energy it
feeds flows on outward and strengthens the companion further out. The 149 galaxies allow that as long as the crests
move at no more than about half the companion's travel speed (85 km/s). Being pulled then costs a kilogram at most
what it gives off when cold: 2 parts in 10¹⁵ of its mass a year. Light-speed crests are ruled out twice over, by the
planets and by the galaxies.

### 4.6 A second route: attraction from the energy of a shared state (round 20)

After the step back (§9), a review of the audit proposed a different kind of mechanism. Every piece of matter would no
longer need to keep one rhythm with the companion. Instead, bodies would attract because where they sit changes the
energy of the combined state of matter and companion, much as two marbles on a stretched sheet roll together because
that lowers the sheet's energy. The review came with a new calculation, which we reproduced exactly.

* **What it already shows.** In a small model of eight pieces sharing one medium, attraction comes from correlations
  between the pieces, with no rhythm on any single piece, whether one piece or several are excited. When every piece
  is excited the pull vanishes, so more stored energy does not automatically mean more pull.
* **A medium that stiffens when energized reaches less far, not farther.** We gave the shared medium a stiffening term
  and energized it (exact statistics on a 32 × 32 × 32 lattice). Its reach fell from 3.2 lattice spacings to 0.9 as
  the energy rose, and the pull fell off faster. Energizing only a region around a source did the same. With the
  medium tuned so that its reach is unlimited, a weak source pulls like Newton (1/r²) and a strong one saturates. The
  law's square root of mass with a 1/r fall never appears.
* **A probe outside a cluster.** For clusters of 3 to 13 pieces and a separate probe, solved exactly (1,260 cases):
  * the probe is always attracted;
  * the pull fades with distance at least as fast as the medium's own short range, and much faster when only a few
    pieces are excited;
  * how it grows with the cluster's size depends on how many pieces are excited, from about in proportion to the size
    down to about its cube root, never the law's square root.

![How fast the pull fades](blog-figures/how-fast-the-pull-fades.png)

*The pull on a probe against its distance from a source cluster of nine pieces, relative to its value at distance
2.5. The law needs 1/r (teal) and Newton falls as 1/r² (dashed); the shared-state model falls far faster, with half
its pieces excited (purple) or one (orange). Script `code/finite_population_probe_v20.py`.*

* **What it means.** A medium that has settled into equilibrium cannot give the law's long reach and its square root
  of mass: whatever we did, it only changed how far the medium reaches. That points somewhere specific. The companion
  is never settled, because matter feeds it energy that flows continuously outward. A medium carrying a steady
  outward flow of energy is the next thing to build, with the same careful bookkeeping.
* **The target is known.** Read the companion as a field whose energy density equals that of the outflow. The balance
  of energy flowing outward is then exactly the equation Bekenstein and Milgrom wrote for MOND in 1984 (in its deep
  regime, with our constant a). So any route of this kind ends up with the MOND family's equation for cold matter, as
  the galaxy data demand. What can be new is where it comes from, the heat term, collisions and the switch-off near
  stars.

## 5. The data

### 5.1 Galaxies: 3,150 measured speeds, three constants

![3,150 measured accelerations in 149 galaxies](blog-figures/rotation-curves-rar.png)

*Every point is a measured orbit in one of 149 galaxies (SPARC), turned into an acceleration and plotted against what
the visible matter alone would give. Newton would put every point on the dashed line. The green curve is our law for
cold disks; hot bulges add to it. Script `code/run_v3.py`.*

| Typical miss per galaxy (km/s) | All 149 | Held-back test galaxies | Adjustable numbers |
|---|---:|---:|---|
| **Our law** | **15.9** | **12.5** | 3 in total |
| MOND | 16.1 | 13.5 | 1 |
| Newton, visible matter only | 45.6 | 41.7 | 0 |
| A dark-matter halo fitted to each galaxy | 7.5 | | 298 |

* The measured pull tracks the visible matter's pull point by point, and a galaxy's visible mass fixes its rotation
  speed. Dark matter has to be tuned galaxy by galaxy to reproduce that tightness; here it is automatic, because the
  companion and ordinary gravity share one source, mass.
* The 25 galaxies with big, hot bulges were the risk, because hot stars count much more. They came out better than
  MOND (29.7 against 30.4 km/s).

### 5.2 Galaxy clusters: the stars sit where the pull is needed

![Measured over predicted cluster mass](blog-figures/clusters.png)

*Twelve clusters with measured gas, temperatures and stars (X-COP). At six radii, the mass needed to hold the gas over
the mass each law predicts; a perfect law sits on ×1. Chart from `code/run_v3.py`, drawn with an earlier fit of
the constants; the current fit, in the project's own distances (`code/xcop_static_v11.py`), moves our points by at most
4%.*

| | Typical miss | Adjustable numbers |
|---|---:|---|
| **Our law** | **25%** | the same 3 constants |
| Our law, clusters held out of the fit | 27% | |
| MOND | ×2.9 | 1 |
| Newton | ×9.1 | 0 |
| Dark matter (a halo fitted to each cluster) | 11% | 24 |

Why the stars can do it: they are concentrated in the middle, exactly where clusters need the most extra pull (near the
very centre there is 1.4 to 7.5 times more mass in stars than in gas), and their random speeds of 300 to 1,200 km/s
make each kilogram count dozens of times over.

**A check without the X-ray input (round 20).** The stars' random speeds used above are worked out from the pull
measured with X-rays, which is the very thing the law is asked to predict.
* With the stars' speeds taken from the law's own pull instead, the typical miss grows from 25% to 40%, or 36% with the
  companion's speed refitted to 186 km/s.
* In the outskirts the law's own pull comes out 26–51% stronger than the X-ray pull. X-COP's own analysis allows only
  about 6–10% of extra support there (inferred using the standard cosmology's share of ordinary matter).

So the cluster fit leans on the X-ray input. Either the heat term is too strong in cluster outskirts, or the outer
galaxies move on stretched orbits or are still falling in; both would lower their heat weight. Measured speed spreads
of the cluster galaxies and lensing masses are the independent test.

### 5.3 Colliding clusters: the pull follows the galaxies

![The Bullet Cluster, our law against its cold limit](blog-figures/bullet-cluster.png)

*The Bullet Cluster, solved in 3D and projected exactly as a lensing map is made. Left: our law puts both lensing peaks
on the galaxies (dots) and little on the gas (diamonds), as measured. Right: the same law with the heat switched off
(MOND-like) puts the lensing on the gas. Chart from `code/bullet_main_v5.py` and `code/bullet_v3.py`; the current
constants in the project's own distances (`code/bullet_static_v11.py`) keep the same pattern.*

* **The pattern usually called impossible without dark matter is reproduced:** both lensing peaks sit on the galaxies,
  and the gas regions carry little extra lensing. It takes both ideas of §3.6–3.7: drop either and the lensing lands
  on the gas.
* **The main cluster matches** its measured galaxy speed (1,249 km/s), its star count (a Legacy Survey count of 1,652
  galaxies confirms its outskirts) and its lensing mass (3.1–3.2 against 3.09–3.46 × 10¹⁴ suns inside 250 kpc,
  including the heat its galaxies picked up crossing the other cluster).
* **72 collisions stacked** (Harvey et al. 2015): the lensing sits a fraction −0.04 ± 0.07 of the way from the galaxies
  to the gas; ours, 0.03. The fresh companion around stopped gas grows at 169 km/s, so it never catches up with gas
  moving away at about 1,000 km/s.
* **The smaller half is still short, and round 21 tested a way to close it.** With the law alone it has 1.4 × 10¹⁴
  suns inside 250 kpc against 2.5–2.9 measured; the heat its galaxies picked up crossing the other cluster raises that
  to 1.6. The project's owner suggested that the heat made in a collision spreads faster than the settled companion,
  at its own speed. At 600 km/s both halves land inside their measured masses (2.5 and 3.4), with both lensing peaks
  still on the galaxies (chart below). The catch is energy: a glow that spreads 3.5 times faster with the same power is
  3.5 times thinner, and with that counted the smaller half reaches only 1.8. So the idea works if the collision puts
  about 3.3 times more power into this fast glow than the law's heat rule gives: a definite target for the physics to
  explain.

![The Bullet Cluster against the speed of the collision's glow](blog-figures/bullet-fast-glow.png)

*The Bullet Cluster's two halves against the speed at which the collision's heat spreads. Shaded: the measured lensing
masses inside 250 kpc. Teal: as proposed; both halves are in range from about 570 to 820 km/s. Orange: with the energy
counted, a faster glow is thinner and the smaller half stays near 1.7–1.8. Solid: the smaller half; dashed: the main
cluster. Script `code/hot_mode_speed_v21.py`.*

Round 22 checked whether the fast glow could simply be glow left behind in the other cluster, like a boat's wake, which
would need no new speed and no extra energy. It cannot: left behind, the glow adds almost nothing to the smaller half
(at most 1.7) and drags its lensing up to 190 kpc onto the gas, which is not what is seen. So the glow must spread
evenly around the stars that made it: a genuinely faster kind of wave in the companion's medium.

![Three more collisions](blog-figures/three-more-collisions.png)

*Three more collisions, with published inputs only (lines: our lensing; circles: galaxies; diamonds: gas). Script
`code/collisions_v10.py`.*

* **MACS J0025.4−1222**, a second Bullet: both lensing peaks sit on their galaxies at the age its shock fronts give
  (0.1–0.36 billion years since the crossing), and its lensing masses agree within their large errors.
* **Abell 520**, a "train wreck" whose middle clump has lensing but few galaxies, long called a dark core: ours gives
  3.7 × 10¹³ suns there from its gas and the surrounding galaxies' heat, against 3.3–3.9 measured, and five of its six
  clumps agree.
* **El Gordo**, seen 7 billion years ago: its lensing inside 1 Mpc comes out 21.4 against 24.3 × 10¹⁴ suns with its
  published stars.

The heat of crossing had never been added to these clusters. Round 21 added it to the two with a clean two-body
history, and it moves both toward their measured masses: El Gordo's lensing inside 1 Mpc from 21.4 to 22.1 × 10¹⁴ suns
with the law's own speed, or 23.7 with the faster glow as proposed (24.3 measured); MACS J0025's two halves from 2.1
and 1.9 to 2.2–2.6 and 2.0–2.5 (3.6 and 3.8 measured, with large errors), with its north-western lensing peak moving
from 79 to 18–44 kpc from its galaxies.

### 5.4 Our own galaxy

![The Milky Way's circular speed](blog-figures/milky-way.png)

*The Milky Way's circular speed from 2 to 28 kpc, against four Gaia measurements. Script `code/milky_way_v7.py`, drawn
with an earlier fit of the constants; the current fit lowers our curve by about 2 km/s at the Sun.*

* **Agrees:** the rotation from 15 to 27 kpc (within 1–6%), the pull above the disk (73 against 68–74 in the usual
  units), the Galaxy's mass inside 100 and 200 kpc (6.5 and 12.4 against 6.1–7.3 and 11.0 × 10¹¹ suns), and the
  escape speed at the Sun (509–525 against 445–580 km/s).
* **Refitted in round 19, with the matter held to independent measurements.** The figure above uses McMillan's (2017)
  matter, which was fitted together with a dark halo: its stars at the Sun (45.8 suns per square parsec) exceed the
  local census (33.4 ± 3) by four standard deviations, and its disk is long. Held instead to the census, the disk's
  measured lengths, the gas and the bulge, and fitted to the newest curve (903 Cepheids, Feng et al. 2026) and the
  pull above the disk together, the law follows the measured curves' shape, with every piece of matter within its
  measured range except a bulge on the heavy side. What remains is an overall level: the law runs 3–6% slow against
  all four Gaia curves. That is about the offset it has for the typical galaxy at this pull (3% in speed in SPARC), so
  the Milky Way is no longer an outlier. Newton with the same matter is 27% slow.

![The Milky Way against 903 Cepheids](blog-figures/milky-way-2026.png)

*The Milky Way from 6.6 to 17.6 kpc: Cepheids measured in 2026 (points) against our law with the Galaxy's matter held to
independent measurements (solid), the same raised by 5.4% (dashed), and Newton with the same matter (orange). Script
`code/mw_joint_v19.py`.*

![Ten Milky Way dwarf galaxies](blog-figures/dwarf-galaxies.png)

*Measured speed spreads of the stars in ten small galaxies around the Milky Way, against predictions. Script
`code/mw_dwarfs_v7.py`.*

* Four agree (Fornax, Leo I, Leo II, Sculptor). Six faint or spread-out ones come out 1.5 to 5 times too slow, because
  the Milky Way's own strong pull holds their companion back. Without that hold, six of the ten agree; the reason the
  hold would be weaker is still being sought.
* **Rerun on today's law without the hold** (round 20; recorded as a comparison, not adopted):
  * the whole suite scores 61 pass, 12 close and 4 fail, against 59, 11 and 7;
  * Carina and Antlia 2 agree, and Sextans and Crater II come close;
  * Draco and Ursa Minor stay below half their measured speeds;
  * wide binaries would then show 18% extra pull at 20,000 AU.

### 5.5 The Solar System and wide binary stars

* **Planets, the star S2 around the Galaxy's central black hole, and binary pulsars:** the release factor switches the
  law's extra pull off where gravity is strong (at the Earth it is e to the power −28 million), so these orbits are as
  in Einstein's theory, which we assume holds in strong fields. The companion's other possible effects there are
  computed only in part (its emission changes the Double Pulsar's orbit by 1% of the measurement error).
* **Cassini's radio tracking** of Saturn limits the Galaxy's distortion of the Sun's field. The 2026 re-analysis (Park,
  Hees, Famaey, Desmond & Durakovic) gives (1.6 ± 1.8) × 10⁻²⁷ s⁻². With the Galaxy's pull and heat at the Sun taken
  from the refitted Galaxy (§5.4), ours is 2.0–2.4 × 10⁻²⁷, within half a standard deviation (with the rougher estimate
  used before, 4.4 × 10⁻²⁷, 1.5 standard deviations). MOND in its widely used "simple" form gives 3 × 10⁻²⁶, about
  twenty times the measured value (forms that switch faster between its two regimes can pass).
* **Wide binary stars** (pairs 5,000–30,000 AU apart, where the mutual pull is weak): the forecast locked in the
  repository in round 10, ahead of Gaia's next release, is 4% more pull than Newton at 7,000 AU and 9% at 20,000 AU,
  against MOND's 43%; it stays as it was. Amended with the refitted Galaxy (labelled as an amendment, not a new
  forecast): 2.4–2.8% and 5.7–6.8%. The Gaia data are disputed.
* **Where the switch-off could come from (round 19).** A derivation supplied to the project proposes one process for
  both of the law's switch-offs: the companion carries temporary "blockers" that stop it coupling to matter; strong
  pull makes more of them, each lasts the same time, and emission adds one. The chance of holding none gives the
  screening exactly, and the extra blocker gives the release length, both from one lifetime (868 years for 0.15 pc).
  We checked the derivation and followed it to the Sun, and it changes the picture near stars. The Sun's companion is
  launched in the Sun's enormous pull and carries out a huge load of blockers, which take many lifetimes to clear. So
  the release becomes a sudden switch far out, not a gradual one. Cassini then only sets a minimum lifetime (about 6
  to 60 years, depending on how emission loads the companion), instead of a value to be tuned. Wide binaries become
  the decisive test: in this version the extra pull is **absent below a certain separation and then jumps to about
  12%**, where the adopted law rises gradually. With its own 868-year lifetime, the Sun's companion stays closed out to
  several light-years and wide binaries are exactly Newtonian. It is recorded as a named candidate, not adopted.

![Wide binaries: a gradual rise or a step](blog-figures/wide-binaries-step.png)

*Extra pull between two Sun-like stars, against their separation. Solid: the adopted law. Dashed, and purple with ten
blocker sites: the blocker switch with the same Cassini standing. Orange: the switch with an 868-year lifetime (none at
all). Far out, every version reaches the same 12% the Galaxy's own pull sets. Script `code/blocker_release_v19.py`.*

## 6. Why stars and light lens the way we see

Lensing is the bending of light by gravity: a galaxy or cluster in front magnifies and distorts the galaxies behind it,
and the distortion measures the pull. Four observations about it, and how the law explains each.

**1. Light and matter feel the same pull.** In our law the companion reshapes the same landscape Φ that ordinary
gravity shapes (∇²Φ = −∇·h), and everything that feels gravity rolls on it, light included. Where the companion's
flow converges, the landscape dips as if extra mass were there. That "as if" mass is what lensing maps show, and what
dark-matter analyses call dark matter. So the mass needed to bend light and the mass needed to move stars must agree,
and in six strong lenses (SLACS), with measured ring sizes and star speeds, they do: to −0.03 ± 0.02 dex in the
project's own distances. We take light to respond to the landscape the way it does in Einstein's theory, twice as
strongly as a slow particle; deriving that factor from the companion itself is part of the relativistic version still
to be built.

**2. Galaxies bend light far beyond their stars, and the bending stays strong far out.** The companion's pull falls as
1/r instead of 1/r², so around an isolated galaxy the lensing signal stays flat far beyond the visible galaxy. In our
law it stays flat to about 2 Mpc, the distance the companion has travelled in the age of the oldest stars, and then
falls: a prediction for the next lensing surveys.

**3. Ellipticals bend light more than spirals with the same stars.** This is the observation that most clearly
separates our law from MOND.

![Weak lensing around 259,000 isolated galaxies](blog-figures/galaxy-lensing-kids.png)

*Weak lensing around about 259,000 isolated galaxies (KiDS-1000). For each acceleration from visible matter, the
observed one, separately for spirals (cold) and ellipticals (hot). Chart from `code/lensing_census_v7.py`; the numbers
below, with the lenses' star speeds measured and in the project's own distances, from `code/kids_heat_v12.py`.*

* An elliptical's stars move randomly at 150–250 km/s and never collide, so in our law they are warm (heat weight
  2.4–6.5) and feed the companion harder. A spiral's stars circle in step and are cold.
* With the lenses' star speeds measured, the ellipticals' extra lensing comes out 0.13–0.16 dex, against the measured
  0.15 ± 0.04. MOND, which treats all matter alike, gives no difference; dark matter has to tune separate haloes.
* All the lenses together sit 16% above our law in the project's own distances. Gas around them weighing about as much
  as their stars, the survey team's own middle estimate, would close that; X-ray and ultraviolet measurements can weigh
  it.

**4. In strong lenses the stars must be somewhat heavy.** Inside a strong lens's ring the stars dominate, and our law
needs them 1.4–1.9 times heavier than the standard "Salpeter" assumption (in the project's own distances). Spectra of
giant ellipticals already suggest star populations this heavy, and measuring these six lenses' stars directly is a
clean test. Compact lenses like the Einstein Cross need only their stars, as observed.

**And in collisions, the lensing follows the companion's flow** (§5.3). Soon after a crossing the hot galaxies carry
their old companion with them, while the colliding gas is cold and its fresh companion has barely begun to grow, so the
light bends around the galaxies rather than the gas, as the Bullet Cluster and 72 other collisions show. The rule is
conditional, not "always on the galaxies": where dense gas sits among hot galaxies whose flow converges on it, as in
Abell 520's galaxy-poor clump, the lensing can sit on the gas, and around gas that has stopped for long enough the
fresh companion grows back.

## 7. How it compares

| | Our law | MOND | Dark matter |
|---|---|---|---|
| Rotation speeds, 149 galaxies | **15.9 km/s** | 16.1 km/s | 7.5 km/s (298 numbers) |
| Cluster masses, 12 clusters | **25%** | ×2.9 | 11% (24 numbers) |
| Collisions: lensing stays with the galaxies | **yes** | no | yes |
| Abell 520's galaxy-poor lensing clump | **from its gas and the galaxies' heat** | no | a puzzle |
| Ellipticals lens more than spirals | **yes, from their stars** | no | yes, via tuned haloes |
| Strong lenses: light and stars agree | **yes**, with light's response assumed | | yes |
| Milky Way rotation 15–27 kpc | **within 1–6%** | within 3% | 6–11% fast |
| Milky Way rotation at the Sun | 9% slow with McMillan's matter; 3–6% with the matter held to its measurements | 3% slow | right |
| Milky Way mass inside 100 / 200 kpc | **agrees** | 30–40% high | agrees |
| Ten Milky Way dwarfs | 4 agree, 6 too slow | the same | fitted |
| Solar System, pulsars | **no extra pull** | small effects | Einstein's |
| Cassini, 2026 re-analysis | **within 0.2–0.4σ** (refitted Galaxy) | common form about 20× too big | passes |
| Wide binary stars (data disputed) | **4% / 9% at 7,000 / 20,000 AU** (locked; 2–3% / 6–7% amended) | 43% | none |
| Fitted constants (not counting the uncertain inputs every model needs) | **4** | 1 | about 320 (two per halo) |

Dark matter fits individual objects more tightly because it is tuned object by object; our law uses four shared
constants and ties the differences between systems to things measured about their matter (how randomly it moves,
whether it collides). But this table is not a statistical comparison: the three columns do not treat the uncertain
inputs (star masses, gas, distances, inclinations) the same way, the misses are not independent, and some data helped
choose the law's form. On the regression suite's 77 graded checks the law scores 59 pass, 11 close and 7 fail,
everything in the project's own distances; that tally tracks regressions, not probability. A likelihood comparison
with shared inputs and a frozen law is on the list (§9).

## 8. Predictions anyone can check

1. **Heavy stars in the six SLACS lenses:** 1.4–1.9 times Salpeter (project distances), measurable from their spectra.
   Dark-matter models expect about 1.
2. **The elliptical/spiral lensing gap stays flat beyond 100 kpc.** A hot-gas-halo explanation makes it grow with
   radius.
3. **Lensing around isolated galaxies stays flat to about 2 Mpc and then falls,** where the companion has not yet
   reached.
4. **Wide binary stars:** 4% more pull than Newton at 7,000 AU and 9% at 20,000 AU, locked ahead of Gaia's next release
   (amended with the refitted Galaxy: 2.4–2.8% and 5.7–6.8%). If the switch-off works by blockers (§5.5), the extra
   pull instead jumps from zero to about 12% at one separation, or is absent altogether.
5. **After a collision, lensing stays with the galaxies while the gas moves away from them faster than a fresh companion
   grows around it** (169 km/s). Around stopped gas it comes back only inside a sphere growing at 169 km/s, about
   170 kpc per billion years, so in old collisions whose gas has long stopped, some lensing should return to the gas.
6. **Older collisions show extra lensing around the smaller clump,** as the heat from crossing and tidal shaking builds
   up.
7. **At equal visible mass, systems whose stars move randomly and freely pull harder** than those whose stars circle in
   step.
8. **Matter falling inward against the companion's outflow glows more;** matter flowing outward with it glows less.
9. **The Milky Way keeps slowing down beyond 25 kpc:** about 189 km/s at 30 kpc and 179 at 50 kpc, not flat.
10. **Isolated galaxies are wrapped in gas weighing about as much as their stars within 100 kpc,** more around
    ellipticals than spirals.
11. **Abell 520's galaxy-poor clump carries about 3 × 10¹³ suns inside 150 kpc,** a quarter of it visible gas.
12. **The Sun loses 2 parts in 10¹⁵ of its mass a year** to its companion.
13. **Frozen on 25 September 2026: lensing grows with the stars' speed at fixed visible mass.** For isolated,
    bulge-dominated galaxies of one stellar mass (4 × 10¹⁰ suns), the pull measured by weak lensing should rise with
    the central speed spread of their stars: relative to 200 km/s, −0.18, −0.09, +0.08 and +0.14 dex at 100, 150, 250
    and 300 km/s, so galaxies at 250 km/s lens about 46% more than those at 150. MOND predicts no difference; dark
    matter predicts whatever difference the halo masses of such galaxies happen to have. The numbers and their
    fingerprint are in `research_work/results/hot-companion/run-frozen-prediction-v19/`.

## 9. What is still open, and why we are optimistic

**A step back (25 September 2026).** Before going further we went through the whole record, every round from the
first idea to round 19, and listed every result we would prefer were different, what we tried, and whether each
"ruled out" holds in general or only inside one particular model. The full list is in
[STEP-BACK-AUDIT.md](STEP-BACK-AUDIT.md). The main points:

* **The law is not cornered.** Its misses (five faint dwarf galaxies, the Bullet Cluster's smaller half, lensing by
  red galaxies, heavy stars in strong lenses, the Milky Way's 3–6%) each have remedies not yet tried, and several
  come from our own conventions rather than physics. For example, galaxy lensing sits 16% above the law in our
  distances but only 3% in the survey's own, so the exact form of our distance law matters; and the speeds of the
  stars in clusters are worked out from the X-ray-measured pull that the law is asked to predict.
* **The models are in a corner we made ourselves, and the record marks the way out.** Since round 10, every small
  model makes the pull the same way: matter behaves like a laser medium that keeps a fixed rhythm with the companion
  wave. Every setback with warm matter since then (§4.4) comes from needing that shared rhythm. Round 13 already
  named the alternative, a pull that follows the companion's energy and flow rather than its rhythm, and it was
  never tried. It is the one that fits what the data select: the pull's strength set by all the companion present,
  its direction by the companion's net flow, warm contributions adding without cancelling.
* **Most "ruled out"s from the models are narrow.** The force experiments never had moving matter (heat was a
  frozen random mixing; an earlier calculation of the glow alone did move its pieces), never reached the heat of
  cluster galaxies, and never had more than about a hundred pieces.
* **What stays closed is what should:** Newton, MOND and dark matter as answers; an extra pull that is not switched
  off near stars; a pull diluted by sources on all sides.
* **Seven small slips in our write-ups were corrected** (listed in the audit); none changes the law or a forecast.

**Round 20: a second route, and a stricter cluster check.** Following a review of the step back, a second way to
derive the pull now has its own line of work: attraction from the energy of the shared state of matter and companion,
not every piece keeping one rhythm (§4.6).
* **The first results:** the attraction is robust, but a settled medium only changes how far it reaches. The law's long
  reach must come from the companion's steady outward flow, which is the next calculation.
* **On the data side:** without their X-ray input the clusters fit worse (40% instead of 25%, §5.2). That makes
  independent measurements of the cluster galaxies' speeds and lensing masses the priority, alongside the frozen
  lensing test.

**Round 21: three proposals checked.** The project's owner sent three ideas for the hardest open problems; each was run
against the data.
* **A faster glow from collisions** (§5.3). Reproduced: at 600 km/s both halves of the Bullet Cluster land in their
  measured masses, and the same heat of crossing brings MACS J0025 and El Gordo closer to theirs. The condition: the
  collision must put about 3.3 times more power into this glow than the law's heat rule gives, which the physics now
  has to explain.
* **A refined distance law.** A factor that grows smoothly with distance, whose simplest value (η = 1/2) may be
  derivable, fits 1,365 supernovae as well as the standard expanding-universe model does, and in the regression suite
  it removes most of the 16% offset in galaxy lensing (62 checks pass, against 59). Two conditions came out. The light
  must arrive as the same photons, each with less energy: the alternative, the same energy spread over more photons,
  is rejected by the supernovae. And the stretching of distant supernovae's light curves, exact to 1% (White et al.
  2024), requires something in intergalactic space that changes slowly over time, about 7.5 parts in 100 billion per
  year, without reaching laboratories: if nothing changed with time, pulses sent one second apart would arrive one
  second apart. It stays a registered comparison until that process is identified.
* **A relativistic version.** The proposed action has the right skeleton (the amount of companion sets the pull's
  size, its flow the direction), and four concrete fixes came out: a sign, a finite response speed, keeping the
  companion's fields independent of the matter's, and putting the extra pull into space-time itself, so that
  gravitational waves feel it as light does. GW170817's waves and light arrived within 1.7 seconds; if only light felt
  the extra pull, the Milky Way's alone would have separated them by about three and a half years.

**Round 22: the fast glow, taken as far as the model allows.** Three results sharpen what the collision glow must be.
* **It must spread evenly** around the stars that made it: a glow left behind in the other cluster adds no mass and
  moves the lensing onto the gas (§5.3).
* **Its energy must come from inside matter.** Even the law's present heat of crossing carries 1.8 times the smaller
  cluster's energy of motion, and the fast glow 5.8 times, so it cannot be paid for by slowing the collision down. It
  needs about 0.05% of the stars' mass-energy, a quarter of what the law already asks of stars in cluster galaxies
  over 10 billion years.
* **The companion's medium must have two regimes.** A collision shakes each star at the rate at which it crosses the
  other cluster's slow wave crests: hundreds to thousands of times the companion's own rate. For the glow to travel 3.5
  times faster there, the medium's wave speed must rise gently at such high rates; a medium of the simplest kind would
  make the glow hundreds of times faster and far too thin to matter. That is now a concrete design target for the
  microscopic model.
* The evenly spreading glow also passes the 72-collision stack: the lensing moves toward the gas by between 0.00 and
  −0.02 of the gas's lag, against −0.04 ± 0.07 measured.

An independent review (25 September 2026) set out what a paper would need. In order:

1. **The companion as a flowing medium.** Done in round 19 for a uniform stream (§4.4): the absorbing stream is a
   proper wave medium (absorbers moving at the wave's speed, at least a wavelength across; the step back found they
   cannot simply be the companion's own stream, which must move faster), it never creates energy and it barely pulls
   on itself. Still to do: the round geometry of a real source exactly, and what sets the absorption.
2. **The law's numbers from the working models.** Done in round 19 (§4.4). The clusters allow only weak one-way
   absorption (a wave must survive about 300 kpc or more inward): the hot glow has to reach inward, as the law's
   total S says. With that, the models give the cold law (mass to the power 0.58 ± 0.13, distance to −1.11 ± 0.22,
   against 0.5 and −1) and the collision rule, but no extra pull from heat on average, and a dependence on a source's size that
   real galaxies don't show. The heat term's mechanism is the open problem now: something that keeps warm matter in
   step, lets the hot glow reach inward, and doesn't depend on size.
3. **Inertia, and why all matter falls alike.** The working models do not yet define a piece's inertia, so equal
   acceleration for different kinds of matter is a requirement, not yet a result.
4. **The release factor and light's response,** derived rather than assumed, and one prescription for the Solar
   System, wide binaries and dwarf galaxies. With the refitted Galaxy, Cassini's 2026 value is met at the adopted
   0.15 pc (§5.4, §5.5); a supplied derivation gives screening and release from one lifetime and predicts a step in
   wide binaries (§5.5), still to be found in the matter model.
5. **Collisions as a calculation from before the crossing,** with the companion's emission, transport and the heat of
   crossing followed in time. The Bullet Cluster's smaller half still has only about 60% of its measured lensing mass;
   a faster glow from the collision closes the gap if the collision puts about 3.3 times more power into it (round
   21, §5.3); round 22 found it must spread evenly around the stars and draw on matter's internal store.
6. **A proper statistical comparison:** a frozen law, a full list of fitted and measured inputs with their
   uncertainties, likelihoods, fair baselines for MOND and dark matter, and at least one test chosen in advance. The
   test is now frozen (§8, prediction 13): lensing against the stars' measured speeds at fixed visible mass.
7. **The remaining misses,** each explained or stated as a limit: six faint dwarf galaxies too slow, all galaxy lenses
   16% above the law (the gas around them to be weighed), and strong lenses needing heavy stars (to be measured). The
   Milky Way is now down to the law's own typical 3–6% (§5.4).
8. **A full literature search** before any claim of priority, and a paper with one bounded claim; the review suggests
   "motion-enhanced attraction in an active streaming medium" first, with the astronomy as motivation. A first check,
   claim by claim, is in [NOVELTY.md](NOVELTY.md). Since the models no longer show the heat gain in the medium the data
   allow, the bounded claim for now is the cold law and the collision rule from an active medium.
9. **A frozen, reproducible release:** every figure regenerated with one version of the law, with inputs, seeds and
   commands archived.

Why optimistic: the four constants have held across galaxies, clusters, collisions, lenses and the Solar System, and
every piece of the law is tied to something measured about matter. Round 19 made the Milky Way ordinary (3–6%, like
the typical galaxy), put Cassini within half a standard deviation, froze a new test that neither MOND nor dark matter
predicts in the same way, and turned the cold law and the collision rule into results of an explicit, energy-balanced
model with error bars. The one link that resists, heat making the pull stronger, is now boxed in precisely: the
data need it and say how it must work (the hot glow reaching inward, weak one-way absorption at most, no dependence on
a source's size), and the models say which mechanisms don't do it. That is a much smaller search than before, and the
frozen lensing test will check the heat term directly in data, whatever its mechanism turns out to be. Round 21
added three concrete leads: a way to close the Bullet Cluster's gap with a stated energy requirement, a distance
factor that fits the supernovae as well as the standard model and fixes galaxy lensing, and a relativistic skeleton
with a short list of fixes. Each open item is a concrete calculation or measurement.

**What is borrowed and what is ours.** Borrowed and credited: Newton's and Einstein's gravity in strong fields; Gauss's
flux geometry; Dicke narrowing and the Mössbauer effect as known physics; the Bloch equations of inverted, self-sustained
emitters (as in the "superradiant laser"); the mathematical form of Milgrom's QUMOND field equation; and the data (SPARC,
X-COP, Clowe et al. 2006, Barrena et al. 2002, Harvey et al. 2015, KiDS-1000, SLACS, Gaia-based Milky Way studies,
dwarf-galaxy catalogues, and the precision tests of gravity; full list in the archived notebook). A first literature
check, claim by claim, is in [NOVELTY.md](NOVELTY.md): it finds, among others, that the heat weight has the form of
Tolman and Whittaker's active mass (with u in place of the speed of light), that in-phase emitters attracting is
Bjerknes', and that one-way coupling through a reservoir is established physics. Ours, as far as we
have found: the companion mechanism, heat as extra glow opened by motion against the companion's flow, collisions
switching it off, the pull along the companion's net flow, its memory after collisions, its release length, MOND's
constant and switch derived from it, and the stream that absorbs inward waves as an exact, passive wave medium. A full literature
search is still owed before any claim of priority.

## 10. Reproduce it

All scripts are in `research_work/results/hot-companion/code/`, and each writes to a fresh output folder. The main
ones:

```
python run_v3.py            --output-dir ../run-v3                 # galaxies, clusters, the constants, the guard
python xcop_static_v11.py   --output-dir ../run-xcop-static-v11    # clusters in the project's own distances
python bullet_static_v11.py --output-dir ../run-bullet-static-v11  # the Bullet Cluster
python collisions_v10.py    --output-dir ../run-collisions-v10     # MACS J0025, Abell 520, El Gordo
python milky_way_v7.py      --output-dir ../run-milky-way-v7       # the Milky Way against Gaia
python mw_dwarfs_v7.py      --output-dir ../run-mw-dwarfs-v7       # ten dwarf galaxies
python strong_field_v7.py   --output-dir ../run-strong-field-v7    # planets, S2, pulsars, Cassini
python kids_heat_v12.py     --output-dir ../run-kids-heat-v12      # galaxy lensing with the measured heat
python derive_mond.py       --output-dir ../run-derive             # MOND as the cold limit
python receivers_v15.py     --output ../run-reservoir-force-v15/receivers_v15.json      # which bodies a wave pulls
python stream_store_v16.py  --output ../run-stream-store-v16/stream_store_v16.json        # heat from the flowing companion
python one_matter_v17.py    --set oneway_distance --output ../run-one-matter-v17/oneway_distance.json  # one kind of matter
python hot_mode_speed_v21.py --grid coarse --output ../run-hot-mode-v21/hot_mode_speed_coarse.json    # the Bullet, faster glow
python distance_eta_v21.py  --output ../run-distance-eta-v21/distance_eta_v21.json         # supernovae and the path factor
python crossing_frame_v22.py --mu 1,0.85,0 --output ../run-crossing-frame-v22/frame.json    # the collision glow's frame
```

The regression suite runs everything at once and compares with the saved baseline:

```
cd research_work/results/hot-companion/regression
python run_suite.py                  # the adopted law, quick tier, about a minute
python run_suite.py --tier full      # everything, about 20 minutes
```

The complete list of scripts and run folders is in the archived notebook's §12 and in the technical record.
