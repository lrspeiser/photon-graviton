# Gravity that streams

*Hot-companion gravity: one law for rotating galaxies, bending light, galaxy clusters and
colliding clusters, with no dark matter and no expanding universe.*

**Rewritten from scratch on 23 September 2026 (rev 12); updated the same day (revs 13–15).**
* Rev 13 added the companion's memory (§3.10).
* **Rev 14** adds §4, the law piece by piece: where each part may come from, and why it works
  so widely. It also brings the Bullet Cluster's galaxy speeds and strong-lensing masses into
  line (§6.3), and corrects how rev 13 read one lensing measurement.
* **Rev 15** tests rev 14 against new data:
  * a star count from the Legacy Survey confirms the main cluster's outskirts;
  * the stars the "bigger smaller half" needs don't show up, so its lensing strength is open
    again (§6.3);
  * the wide-binary prediction, 19% extra pull, is new (§6.7).
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
   outer stars orbiting fast.
3. **Matter whose parts move randomly and freely, like stars in elliptical galaxies or
   galaxies in clusters, feeds it out of step.** Nothing cancels, so it pulls much harder.
   That is why clusters need so much extra pull.
4. **Gas does not count as "hot" in this sense.** Its particles collide so often that the
   companion sees them as sitting still. That is why, in a collision, the extra pull stays
   with the galaxies and leaves the gas behind.
5. **Strong gravity holds the companion back,** so the Solar System feels nothing.
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
| Ellipticals bend light 0.17–0.27 dex more than spirals | **0.17–0.27** | no difference | yes, with tuned haloes |
| Solar System | **no extra pull** | small extra pull | no extra pull |
| Wide binary stars (data disputed) | **19% extra pull predicted** | 43% | none |
| Adjustable numbers | **3 in total** | 1 | 2 per galaxy or cluster (~320) |

**What is still open:**
* **the Bullet Cluster's smaller half.** Its lensing mass is twice what our law gives it.
  Rev 14's explanation, that it was a bigger cluster before the crash, is not borne out by the
  galaxies and starlight around it (rev 15). The leading candidate now is heat its galaxies
  picked up while crossing the main cluster, carried along with them;
* a full relativistic version of the law.

Both are in §9, with what could close them. Rev 13's open item, the main cluster's galaxy
speeds, is settled: half of the gap was a modelling artefact, and a normal amount of stars in
the cluster's outskirts closes the rest (§6.3).

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

The constant a is not a new constant of nature here. It is twice the power per kilogram
divided by the companion's speed. With the fitted values, ℓ = 6.5 × 10⁻⁶ watts per
kilogram. Over 13.8 billion years that is three thousandths of one percent of a kilogram's
rest energy.

### 3.4 Strong gravity holds the companion back

Near the Sun the pull is billions of times stronger than in a galaxy's outskirts. If the
companion pulled there, planetary orbits would show it; they don't. So the companion must
be held back, "attached", where ordinary gravity is strong, and released where it is weak:

```
released fraction = exp( −|g_N| / g_d )
```

Here g_d = 2.26 × 10⁻¹⁰ m/s², measured on the galaxies. At every planet the released
fraction is zero to the precision of a computer. The same factor also improves the galaxy
fits.

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
```

The last line turns the rule into a proper field that conserves energy for anything moving
through it. Mathematically it has the same shape as Milgrom's QUMOND field equation, which
we use for its structure, not its physics (§11).

**Three constants, all universal, none set per object:**

| Constant | Value | Measured on |
|---|---|---|
| a | 6.56 × 10⁻¹¹ m/s² | 149 galaxies |
| g_d (release level) | 2.26 × 10⁻¹⁰ m/s² | 149 galaxies |
| u (companion speed) | 197 km/s | 12 clusters |

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

**It is slow.** It streams at u = 197 km/s, measured on clusters. At that speed it takes about
500 million years to travel 100 kpc. So most of the companion around a galaxy or a cluster
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
with the stars carrying the heat, lowered u to 197 km/s for a separate reason. The fresh
companion around stopped gas now grows 4.4 times more slowly, and it never catches up with gas
separating from its galaxies at 1,000 km/s or more.

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

### 4.2 ℓ: every kilogram feeds the companion, 6.5 × 10⁻⁶ watts per kilogram

**What it is.** A steady power fed by each kilogram of matter into the companion. It is tiny.
Over the age of the universe it adds up to three thousandths of one percent of a kilogram's
rest energy. For the Sun it is 1.3 × 10²⁵ watts, 3.4% of the Sun's light.

**Where it may come from.** This is the project's founding intuition: energy converts into
gravity. Our working guess is that matter continuously sheds a sliver of its energy into a
gravitational companion, the way a warm object radiates.

There is a numerical clue. The project's redshift work (no expansion) has light losing energy
at a steady rate, c·α ≈ 70–75 km/s per megaparsec. Multiplied by the speed of light, that rate
is an acceleration, c·α ≈ 7 × 10⁻¹⁰ m/s². Our two acceleration constants sit close to it:
* a ≈ c·α / 10;
* g_d ≈ c·α / 3.

If light and matter both exchange energy with the same background, these numbers would be
linked. A similar near-match has long been noticed for MOND's constant, so this could be a
coincidence. If it is not, it is the strongest hint of where the feed comes from.

**What it explains.** It sets the strength of everything the companion does, through a = 2ℓ/u.

**How to test it.** The Sun's extra mass loss is predicted at 2.3 × 10⁻¹⁵ of its mass per year,
below what planetary tracking can see today. A theory that links ℓ to the redshift rate would
turn "a ≈ c·α / 10" into an exact prediction.

### 4.3 u: the companion's speed, 197 km/s

**What it is.** How fast the companion streams away from the matter that fed it, measured
relative to that matter.

**Where it may come from.** A slow wave. Sound moves at a speed set by the stiffness and weight
of the air; a wave on a guitar string, by its tension and weight. A companion at 0.07% of the
speed of light suggests a heavy, sluggish mode of the gravitational field. It is not an ordinary
gravitational wave: those travel at the speed of light, as the 2017 neutron-star merger showed,
and they are a separate thing.

Two clues, which may be coincidences:
* **u is the rotation speed of a big galaxy.** Our law gives a flat rotation speed of exactly u to
  a galaxy of 1.7 × 10¹¹ suns (from v⁴ = G a M), which is a Milky-Way-sized galaxy.
* **u²/a ≈ 19 kpc**, the size of a large galaxy's disk.

Both hint that galaxies may be sized by the companion.

**What it explains.** This is the key to the whole law. Through the heat weight k = 3σ²/u², u is
the dividing line between cold and hot:

| System | Random speed of its stars or galaxies, σ | Heat weight k |
|---|---:|---:|
| Disk of a spiral galaxy | 10–30 km/s | 0.008–0.07 |
| Elliptical galaxy | 150–250 km/s | 1.7–4.8 |
| Galaxies in a cluster | 500–1,200 km/s | 19–110 |

So with a single number:
* disks come out cold, and follow the MOND-like cold limit;
* ellipticals come out warm, and bend light more than spirals;
* clusters come out hot, and get the large boosts they need.

u also sets the companion's memory after a collision: 500 million years per 100 kpc (§3.10).

**How to test it.**
* In a collision, lensing returns around stopped gas only as a sphere growing about 200 kpc per
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
  6.6 × 10⁻¹¹ m/s², in the outskirts of galaxies.

**How to test it.** The power of exactly 4 in v⁴ ∝ M, and the same switch-over pull a in every
galaxy.

### 4.5 The release factor exp(−|g_N|/g_d): held back where gravity is strong

**What it is.** The share of the companion that is free to pull. It is almost zero where
ordinary gravity is strong and almost one where it is weak. The level is g_d = 2.26 × 10⁻¹⁰ m/s².

**Where it may come from.** Two guesses.
1. **Escape over a barrier, like evaporation.** Suppose the companion must climb a fixed height ℓ_d
   against the local pull to break free of its source, and its energies are spread like those of
   a warm gas, about u²/2 per kilogram on average. The share with enough energy to climb is
   exp(−2 g ℓ_d / u²). That is exactly our form, with g_d = u²/2ℓ_d, so ℓ_d ≈ 2.8 kpc. The same
   maths describes molecules evaporating from a liquid and stars leaking out of a star cluster.
2. **Screening.** In a plasma, electric charges are hidden beyond a short distance, and the
   hiding follows an exponential. Exponential cut-offs are the signature of screening.

Guess 1 ties g_d to u and to a length of a few kiloparsecs, about the size of a galaxy's core.
If that length is also the companion's wavelength, it would join this piece to the collision
rule (§4.7).

**What it explains.**
* **The Solar System is silent.** At the Earth, g_N is 26 million times g_d, so the released
  share is e^(−26,000,000): zero.
* **Galaxy centres behave like Newton.** Where g_N ≈ 10⁻⁹ m/s², only about 1% is released.
* **The smooth bend of rotation curves** between the Newtonian centre and the flat outskirts. The
  galaxy fits improved when this factor was added.
* **Clusters are mostly released**: 64–96% for g_N between 10⁻¹⁰ and 10⁻¹¹ m/s².

**How to test it.**
* **Wide binary stars near the Sun.** They sit in the Galaxy's pull, where our law has only about
  half of the companion released. It predicts 19% more pull than Newton beyond about 7,000 AU,
  against MOND's 43% (§6.7).
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

**Why it matters so much in clusters.** Inside a big, spread-out cloud, Newton's pulls from all
sides cancel toward the centre, but plain totals don't. So near a cluster's centre S stays large
while g_N shrinks toward zero. That is exactly where clusters need 4 to 27 times the visible
matter's pull.

**What it explains.**
* **Cluster masses:** a 25% typical miss, against MOND's factor of 2.9.
* **Ellipticals bend light 0.17–0.27 dex more than spirals** with the same stars.
* **Galaxies with big, hot bulges.**
* **Why MOND works for disk galaxies (k ≈ 0) and fails for clusters (k ≈ 20–100).**

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
* Lensing around stopped gas comes back at about 200 kpc per billion years.
* Older collisions show extra lensing around the smaller clump.

### 4.12 Why so few pieces cover so much

Each piece is a switch set by something measurable about the matter itself. None is a number
tuned object by object:

| Switch | Set by | Separates |
|---|---|---|
| Release factor | how strong the pull is | the Solar System from galaxies |
| Heat weight | how randomly the matter moves | disks, ellipticals and clusters |
| Collision rule | whether its particles collide | gas from stars |
| Direction | which way its energy flows | round systems from collisions |
| Travel | how recently it changed its motion | settled systems from collisions |

Each kind of system flips a different combination of switches. The three constants never change.

**Which pieces each observation needs.** ● means essential; ○ means involved.

| Observation | Newton | √(a ·) companion | Release | Heat k, S | Gas cold | Direction | One landscape | Travel |
|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| Solar System and the Sun | ● | | ● | | ● | | | |
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

1. **a ≈ c·α/10 and g_d ≈ c·α/3,** where c·α is the redshift rate of light turned into an
   acceleration (§4.2). If the companion and the redshift share a cause, the ratios 10, 3 and
   λ = g_d/a = 3.45 are numbers that theory must derive.
2. **u⁴ = G a × (1.7 × 10¹¹ suns).** The companion's speed is the rotation speed of a
   Milky-Way-sized galaxy.
3. **u²/(2g_d) ≈ 2.8 kpc,** a detachment length, if the release factor is an escape over a
   barrier (§4.5).
4. **The Sun's companion carries 3.4% as much power as its light.**

What would pin these down:
* a field theory that gives the companion's speed (§9);
* a measurement of its wavelength, from which systems count as hot and which as cold.

## 5. MOND, as the cold limit of our law, not an input

Switch off the heat (S = 0), as for a cold disk galaxy, and the law becomes

```
g = g_N · ν(g_N/a),   ν(y) = 1 + e^(−y/λ) / √y,   λ = g_d/a = 3.45
```

That is a MOND-type formula. **But it came out of the companion mechanism; it was not put
in:**

| MOND's ingredient | In MOND | Here |
|---|---|---|
| Its constant a₀ | a new constant of nature | 2 × power per kilogram ÷ companion speed |
| The switch between Newton and the deep regime | chosen by hand | set by the release factor. Not identical to any published choice; the nearest (McGaugh, Lelli & Schombert's exponential) is 0.030 dex (7%) away |
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
| **Ours** | **15.85** | 15.95 | **19.21** | **12.40** | 3 in total |
| MOND | 16.13 | 15.85 | 19.80 | 13.52 | 1 |
| Newton, visible matter only | 45.58 | 45.94 | 48.64 | 41.68 | 0 |
| Dark-matter halo fitted to each galaxy | 7.52 | | | | 298 |

* Ours beats MOND overall and on the held-back galaxies; MOND is slightly better on the
  training set.
* **Galaxies with big bulges** were the risk. Their bulges are hot, and the new heat rule
  makes hot stars count much more. Those 25 galaxies came out better: 29.2 km/s against
  MOND's 30.4.
* Dark matter fits more tightly because it has a hundred times more adjustable numbers.

Script: `code/run_v3.py`.

### 6.2 Galaxy clusters: 12 clusters with measured gas, temperatures and stars

The quantity is the typical miss in the mass needed to hold the gas, at six radii from the
core to the edge.

| | Typical miss | Adjustable numbers |
|---|---:|---|
| **Ours** | **25%** (rms 0.227 in ln M) | same 3 constants |
| Ours, clusters held out of the fit (924 splits) | 28% (0.244) | |
| Ours, stars' speeds computed by our law from visible matter alone | 39% (0.329) | |
| MOND | ×2.9 (1.062) | 1 |
| MOND with its constant refitted on clusters | ×1.4 (0.338); needs a constant 9.7× its galaxy value | 1 |
| Newton | ×9.1 (2.211) | 0 |
| Dark matter (NFW fit per cluster) | 11% (0.101) | 24 |

**Why the stars can do it.** The stars are concentrated in the middle, where clusters need
the most extra pull:

| Radius (fraction of the cluster's size, R500) | 0.02 | 0.05 | 0.1 | 1 |
|---|---:|---:|---:|---:|
| Mass in stars ÷ mass in gas | 2 to 11 | about 1 | 0.3 to 0.5 | 0.02 to 0.07 |

Their random speeds, roughly 300 to 1,200 km/s, make each kilogram count dozens of times
over.

The companion speed comes out at u = 197 km/s. Across the 924 ways of splitting the clusters
in half, it stays between 178 and 218 km/s in 90% of them.

Script: `code/run_v3.py`, with stellar mass profiles from the X-COP release (Ghizzardi et al.
2021).

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
* The fresh companion around the gas grows at 197 km/s, so it never catches up with gas
  moving away at about 1,000 km/s. That is exactly where round 2's version failed.

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
* **The stars need to be 1.05–1.35 times heavier than the "Salpeter" standard assumption.**
  Spectra of giant ellipticals already suggest star populations this heavy. Measuring these
  six lenses' star populations directly is the clean test.
* The mass needed to bend the light and the mass needed to move the stars agree to
  −0.017 ± 0.024 dex. **Light and matter feel the same pull.**

Script: `code/lenses_t35.py --constants ../run-v3/results.json`.

### 6.5 The Solar System and the Sun

* **Planets.** At every planet, the release factor switches the extra pull off entirely.
* **The Sun.** Its interior is colliding plasma, so it has no heat term. It feeds the
  companion at ℓ = 6.5 × 10⁻⁶ W/kg: 1.3 × 10²⁵ W, or 3.4% of its light. That means an
  extra mass loss of 2.3 × 10⁻¹⁵ of its mass per year, below what planetary tracking can see.

### 6.6 The formula check

* **Not MOND.** For the same Newtonian pull, random speeds of 0–1,000 km/s change the
  prediction by a factor of 8.6 (0.93 dex). Across the 3,150 galaxy measurements, 0.022 dex
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

Script: `code/wide_binaries_v6.py`.

## 7. How this compares

| | Ours | MOND | Dark matter |
|---|---|---|---|
| Rotation speeds, 149 galaxies | **15.9 km/s** | 16.1 km/s | 7.5 km/s (298 numbers) |
| Cluster masses, 12 clusters | **25%** | ×2.9 | 11% (24 numbers) |
| Bullet Cluster: lensing on galaxies | **yes** | no | yes |
| Bullet Cluster: main half | **lensing mass 2.4–2.5 against 2.5–2.8; galaxy speed 1,249 km/s; star count agrees** | no | yes (fitted) |
| Bullet Cluster: smaller half's lensing mass | about half of the measured 2.0–2.3 (open) | no | yes (fitted) |
| Wide binary stars (data disputed) | **19% extra pull beyond 7,000 AU** | 43% | none |
| Collisions: lensing stays with galaxies | **yes** | no | yes |
| Ellipticals lens more than spirals | **yes, 0.17–0.27 dex, from their stars** | no | yes, via tuned haloes |
| Strong lenses: light and stars agree | **yes** | | yes |
| Solar System | **silent** | small effect | silent |
| Explains *why* | disks are cold; cluster galaxies are hot; gas collides | no | no |
| Adjustable numbers | **3** | 1 | ~320 |

Dark matter fits individual objects more tightly because it is tuned object by object. Ours
fits everything with three shared numbers, and says why each kind of system behaves as it
does.

## 8. Predictions anyone can check

1. **Heavy stars in the six SLACS lenses:** 1.05–1.35 times Salpeter, measurable from their
   spectra. Dark-matter models expect about 1.0.
2. **The elliptical/spiral lensing gap is flat beyond 100 kpc.** A hot-gas-halo explanation
   makes it grow with radius.
3. **In relaxed clusters, galaxy orbits and gas agree.** Both feel the same pull.
4. **In every collision, lensing stays with the galaxies**, at every stage.
5. **The Bullet Cluster's smaller half carries a companion heated during its crossing,** so its
   lensing is strong although few galaxies travel with it. Deep images should also show
   whatever diffuse starlight travels with it.
6. **The Bullet Cluster's main outskirts** hold stars at about 5% of the gas mass at the
   cluster's edge. Confirmed by the star count in rev 15.
7. **Wide binary stars:** 19% more pull than Newton beyond about 7,000 AU, less than half of
   MOND's 43%.
8. **Lensing in a collision shows the clusters as they were before it.** Around gas that has
   been stopped, lensing comes back only inside a sphere growing at 197 km/s, about 200 kpc per
   billion years.
9. **Older collisions** (half a billion to a billion years after the pass) show extra lensing
   around the smaller clump, as the stream heat and tidal shaking build up.
10. **The extra pull tracks random motion.** At equal visible mass, systems whose stars move
   randomly and freely pull harder than those whose stars circle in step.
11. **The Sun's extra mass loss** is 2.3 × 10⁻¹⁵ per year.

## 9. What is still open, and why we are optimistic

**Settled:** the main Bullet cluster's galaxy speeds and lensing mass (rev 14), and its star
count (rev 15).

1. **The Bullet Cluster's smaller half (§6.3).** Its measured lensing mass, 2.0–2.3 × 10¹⁴ suns,
   is twice what our law gives it. Rev 14's explanation, a bigger cluster before the crash, is
   not borne out: the galaxies and starlight it would need are not there. Next:
   * **model the heat its galaxies picked up while crossing the main cluster**, carried with
     them by the memory rule;
   * look for diffuse starlight moving with it;
   * measure a few hundred more galaxy speeds around it.
2. **A full field theory for the companion.** It needs its momentum, its travel time (the
   memory of §3.10), and a relativistic form, so that lensing is derived rather than assumed
   equal to the pull on matter. The data support that equality: the six lenses agree to
   0.017 dex.
3. **The microscopic numbers.** We want to derive u and g_d, and pin down the companion's
   wavelength. That wavelength must lie between the gas's collision length and the stars'
   orbits. For cluster gas this requires particles to be scattered, mainly by magnetic
   fields, at least every 0.1–1 million years. That is thought likely but is not measured.
4. **Wide binary stars (§6.7).** Our prediction is now computed: 19% more pull than Newton
   beyond about 7,000 AU. The two published analyses find about 40% and none. Whichever holds
   up will tell us how much of the companion is released near the Sun.
5. **Cosmology.** The cosmic microwave background and the growth of large-scale structure are
   outside this law's scope as tested so far. Under the project's no-expansion rule they need
   their own treatment.

Each of these is a concrete calculation or measurement, not a wall.

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
* **Revision 15 (round 6, this page).**
  * A Legacy Survey star count confirms the main cluster's outskirts.
  * The galaxy speeds and starlight don't show the stars a bigger smaller half would need, so
    its lensing strength is open again.
  * The wide-binary prediction.

**Superseded along the way, kept on the record:**
* round 2's hot-gas-halo explanation of the ellipticals (now it is their stars);
* round 2's reaction force on hot matter (now the companion carries momentum);
* round 1's companion speed (874 km/s, now 197);
* round 2's solar mass-loss figure (1.4 × 10⁻¹⁴ per year, now 2.3 × 10⁻¹⁵);
* rev 12's "the main cluster's galaxy speeds check out" (a single radius; the average over the
  surveyed region was 20% low in rev 13's model, §6.3);
* rev 13's "all four Bullet lensing strengths within about 1σ": Clowe et al.'s κ values are
  lower bounds, as their paper states. The calibrated test is the strong-lensing mass, and by
  that measure rev 13's smaller half was half as heavy as observed;
* rev 13's lighter star masses (mass-to-light 1–1.5): never required;
* rev 14's "the smaller half was a third of the main before the crash": the galaxies and
  starlight that would need are not seen (rev 15);
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
  * KiDS-1000 (Brouwer et al. 2021);
  * SLACS (Bolton, Auger, Treu and collaborators).

**Ours, as far as we have found:**
* the companion mechanism;
* heat as loss of step;
* collisions switching it off;
* the pull along the companion's net flow;
* the companion's memory after collisions;
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
python kids_v3.py        --output-dir ../run-kids-v3        # ellipticals vs spirals
python lenses_t35.py     --output-dir ../run-lenses-v3 --constants ../run-v3/results.json
python derive_mond.py    --output-dir ../run-derive         # MOND as the cold limit
python field_equation.py --output-dir ../run-field          # the 3D field equation
```
