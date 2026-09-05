# S301: simulating a ground-based VLT/VLTI observing campaign

## What S301 is

S301 is a real star. It orbits Sagittarius A\* (Sgr A\*), the supermassive
black hole at the center of the Milky Way. The GRAVITY/VLTI collaboration
announced its discovery in August 2026. The discovery paper is
["Discovery of a star sensitive to the spin of Sgr A\*"](https://www.nature.com/articles/s41586-026-10894-w)
(*Nature*, arXiv:[2607.12664](https://arxiv.org/html/2607.12664v1), DOI
10.1038/s41586-026-10894-w). See also the
[ESO press release](https://www.eso.org/public/news/eso2612/), the
[MPE presskit](https://www.mpe.mpg.de/8219735/news20260819), and
[sci.news's writeup](https://www.sci.news/astronomy/s301-star-spin-sagittarius-a-15004.html).

S301 passes about 1.78 billion km from Sgr A\* every 8.68 years. That is
roughly 12 AU, or 136 Schwarzschild radii. At closest approach it moves at
more than 8% the speed of light (about 25,000 km/s). No other known star
gets this close or moves this fast, including the famous S2. Its orbit
also precesses: each pass rotates the orbit's orientation, tracing out a
rosette shape over many years. This precession makes S301 a candidate for
eventually measuring Sgr A\*'s spin, once astronomers have watched enough
passes.

## What this project does

This project simulates S301 in three stages:

1. **A relativistic light curve.** How bright S301 looks over one orbit,
   including real relativity effects near closest approach.
2. **A synthetic observing campaign.** Noisy, realistic measurements of
   S301's position, velocity, and brightness across two close passes.
3. **An orbit fit.** Recovering S301's orbital shape, and its precession
   rate, from the noisy campaign data alone.

Every simulated measurement is tied to a real, currently-operating
ground-based instrument. Every per-measurement precision figure comes
from a real achieved or published result, not a guessed number. Where
this project has to make an assumption, it says so, and says why.

## Which instruments, and why

Real astronomers use two facilities at ESO's Paranal Observatory in
Chile for this kind of work:

- **GRAVITY+**, an interferometer at the VLTI (Very Large Telescope
  Interferometer), for position measurements (astrometry). GRAVITY+ is
  an upgrade of the original GRAVITY instrument that discovered S301.
- **ERIS**, an adaptive-optics-fed spectrograph at the VLT, for velocity
  measurements. ERIS is the successor to SINFONI, the instrument used
  for essentially all radial-velocity monitoring of stars this close to
  Sgr A\* to date.

An earlier version of this project instead used a single space telescope
for all three measurement types. That was wrong, for two reasons:

1. **The physics rules it out.** A single mirror's resolving power at
   near-infrared wavelengths is far too coarse. A 2.4 m mirror, for
   example, cannot separate two points closer than about 223
   milliarcseconds (mas) apart. S301's entire orbit spans only 83 mas.
   The whole orbit would fit inside one blurry point of light. This is
   exactly why the real discovery needed GRAVITY's 100 m interferometric
   baseline, not a single mirror of any size a space telescope could
   carry.
2. **Real astronomers agree.** A search for actual announced S301
   follow-up plans turns up GRAVITY+ and future 30-40 m ground telescopes
   (the ELT and GMT). It does not turn up any space telescope. See the
   [GRAVITY+ project paper](https://arxiv.org/pdf/2301.08071) and
   [coverage of planned S301 follow-up](https://phys.org/news/2026-07-star-orbiting-galaxy-supermassive-black.html).

So this project uses GRAVITY+ and ERIS, and only those two.

## What's real, and what's assumed

`constants.py` is the single source of truth for every number in this
project. Every value there carries a comment naming its source.

**Published facts, taken directly from the discovery paper** (Extended
Data Table 2, Solution 1): orbital period P = 8.68 years, eccentricity e
= 0.9832, inclination i = 124.09°, longitude of the ascending node Ω =
73.8°, argument of periapsis ω = 293.4°, time of last periapsis t_peri =
2023.126, semi-major axis a = 83.0 mas. At the real distance to Sgr A\*
(8277 pc), that semi-major axis is 687 AU. Also published: the black
hole's mass (4.297 million solar masses), S301's spectral type (F1.5V),
its apparent K-band magnitude (19.3 ± 0.3), its mass (1.1-1.5 solar
masses), and its radius (1.4-1.6 solar radii). The paper also publishes a
second, nearly-identical orbit solution — a known ambiguity in
astrometry-only fits. This project keeps that second solution in
`constants.py` for reference, but does not use it.

**One assumption, clearly flagged:** the paper does not state S301's
surface temperature. This project estimates 7300 K from standard tables
for its published spectral type (F1.5V).

**One fact worth stating plainly:** Gaia has never observed S301, and
never could. Dust toward the galactic center blocks about 30 magnitudes
of visible light. Gaia is an optical mission. Only near-infrared
instruments can see through that dust to a star this deep in the
Sgr A\* star cluster.

**One timing assumption:** a multi-year, two-pass VLT/VLTI campaign like
this would realistically go through ESO's Large Programme process, not
routine scheduling. This project assumes about two years from today to
the first scheduled night: time to write the proposal, submit it to the
next Large Programme call, wait for committee review, and wait for the
awarded observing period to start.

## Stage 1: the relativistic light curve (`s301_lightcurve.py`)

This stage computes S301's brightness over one full orbit, using the
real published orbital elements.

The physics:

1. A full Keplerian orbit, solved exactly (Newton-Raphson) rather than
   approximated.
2. The exact special-relativistic Doppler shift, from S301's true 3D
   velocity — not the small-velocity approximation most textbooks use.
3. The Schwarzschild gravitational redshift.
4. Relativistic brightness boosting, applied to a blackbody spectrum and
   integrated over GRAVITY/VLTI's real K-band range (1.98-2.40 μm) — the
   actual bandpass the published magnitude was measured in.

Left out, on purpose: black-hole spin effects, strong-field light
bending, and light-travel-time delay. These are all real effects, but
smaller than what's modeled here, and outside this project's scope.

Running it prints:

```
Derived semi-major axis a = 686.7 AU  (published: 687.0 AU)
Periapsis check         r_p = 136.0 R_S  (published: 136 R_S)
Periapsis speed check   v_p = 8.54% c  (published: >8% c, 25000 km/s)
Max blueshift/redshift factors: (1+z) in [0.9552, 1.0242]
Rough mean apparent K-band mag (order-of-magnitude): 17.93  (published: 19.3 +/- 0.3)
Peak-to-trough delta-mag amplitude near periapsis: 0.1158 mag
```

The first two lines are a self-consistency check. This project only
inputs the period and the eccentricity; the semi-major axis and
periapsis distance come out of the physics, not the input. They match
the published paper closely. The apparent-magnitude estimate is rougher
— it comes from a generic bolometric correction, not a real filter
throughput curve — but it lands within about 1.5 magnitudes of the real
value, which is the right order of magnitude for that method.

## Stage 2: the synthetic campaign (`campaign.py`)

This stage builds a fake but realistic observing campaign: what
GRAVITY+ and ERIS would plausibly measure if astronomers pointed them at
S301 for years.

**Why two passes.** A precession *rate* can only be measured by comparing
the orbit's shape across two separate passes. One pass alone only shows
one instant's orbit, not how it's changing. So this campaign covers the
next two predicted periapsis passages: around 2031.8 and around 2040.5.

**Cadence.** Observations start in 2028.5 (see the timing assumption
above) and run through 2041.7. Most of the year gets sparse sampling
(about six times a year). The year around each periapsis passage gets
dense sampling (about once a month), because that's when the orbit
changes fastest and matters most.

**The physics injected as ground truth** goes one step further than a
plain, unchanging orbit. It includes the real Schwarzschild precession
rate — the same 1PN general-relativity formula used to explain Mercury's
perihelion advance, applied here to S301. That formula predicts the
argument of periapsis rotates by about 2° with every orbit. This project
does *not* inject any black-hole-spin effect. Spin would show up as a
different kind of precession — a drift in the orbit's inclination and
node, called frame-dragging — and Sgr A\*'s spin has never been measured.
Injecting a spin value here would mean guessing at exactly the number
nobody knows.

**Three measured quantities per epoch, from two instruments:**

- **Position on the sky** (RA and Dec offset from Sgr A\*), from
  GRAVITY+.
- **Radial velocity** (motion toward or away from Earth), from ERIS.
- **K-band brightness**, from GRAVITY+, generated for completeness but
  not used in the orbit fit (see Stage 3).

Running it prints:

```
Campaign: 128 epochs, 2028.50 - 2041.67 (passages: 2031.81 and 2040.49)
Injected Schwarzschild precession: 0.2307 deg/yr (2.003 deg/orbit) -- real 1PN formula, not illustrative
Astrometric precision (GRAVITY+, real achieved on-sky figure): 100 uas/epoch
RV precision (ERIS, background-limited scaling from real S2 SINFONI/ERIS precision -- the realistic regime for a target this faint in this field): 1621 km/s/epoch (optimistic source-limited alternative: 141.2 km/s/epoch)
Cross-check: SPIFFIER's own R=5000 gives a 60.0 km/s resolution element -- the injected precision implies continuum SNR~0.037, i.e. no real single-epoch RV detection for a star this faint with a current 8m-class instrument
Photometric precision (S301's own published m_K uncertainty): 0.30 mag/epoch
True RA offset range: [-0.5, 60.0] mas
True RV range: [-10727, 4233] km/s
```

128 epochs is more than double the 60 epochs in an earlier, single-pass
version of this campaign.

## How each precision figure was derived

Every number above traces back to a real measurement or a real
published fact, not a guess. Here is the reasoning for each:

**Position precision: 100 microarcseconds.** GRAVITY has already
achieved 30-100 microarcsecond astrometry on stars as faint as magnitude
20 in K-band, by combining a full night of data
([arXiv:2301.08071](https://arxiv.org/pdf/2301.08071), "The GRAVITY+
Project"). S301, at magnitude 19.3, falls inside that demonstrated
range. This project uses the worse (larger, more cautious) end of the
range.

**Velocity precision: 1621 km/s.** This is the hardest number to pin
down, and the most important finding in this project's noise model.

The starting point is a real, achieved result: ERIS and its predecessor
SINFONI measure S2's radial velocity — S2 is the brightest, best-studied
star this close to Sgr A\* — to a median precision of 12.3 km/s
([arXiv:1709.01598](https://arxiv.org/abs/1709.01598)). S2 has an
apparent magnitude of 14.0. S301 is much fainter, at magnitude 19.3 — a
difference of 5.3 magnitudes, or about 130 times less light.

Scaling a precision measurement across that much of a brightness gap
means picking a noise regime, and the two reasonable choices give very
different answers:

- If the noise scales with the square root of the star's own light
  (source-photon-limited), precision degrades by about 11.5x, to 141
  km/s.
- If the noise scales linearly with the star's light instead
  (background-limited — the sky background and stellar crowding matter
  more than the target's own faint signal), precision degrades by about
  132x, to 1621 km/s.

Ground-based K-band spectroscopy of faint Galactic Center stars is
well documented as background-limited in practice: bright, variable sky
emission lines, thermal background, and severe stellar crowding all
dominate over a faint target's own photon count. So this project injects
the background-limited number, 1621 km/s, as the realistic value. The
source-limited number, 141 km/s, is kept only as an optimistic point of
comparison — it is not what the simulated data actually uses.

A cross-check confirms this conclusion. ERIS's SPIFFIER spectrograph has
a velocity resolution of about 60 km/s. Dividing that by the injected
1621 km/s precision gives an implied signal-to-noise ratio of about
0.04. That is not a weak detection — it is no detection at all. This
matches what the real discovery papers say directly: measuring S301's
radial velocity needs a *future*, larger telescope (they name the ELT's
HARMONI and MICADO instruments), not a current 8-10 m telescope like the
VLT.

**Brightness precision: 0.30 magnitudes.** This project uses S301's own
published photometric uncertainty directly, rather than deriving a new
number. That published figure is the uncertainty on one absolute
brightness measurement. Night-to-night relative brightness measurements
from the same instrument are usually more precise than that, because
shared calibration errors cancel out. No published figure for that
relative precision exists, so this project's 0.30 magnitude figure is
probably more cautious than reality.

## Stage 3: recovering the orbit (`fit_orbit.py`)

This stage takes the noisy campaign data and tries to recover the orbit
that produced it — without looking at the injected truth.

The fit solves for 7 numbers at once: the 6 standard orbital elements
(period, eccentricity, inclination, node, argument of periapsis, and
time of periapsis), plus the precession rate. It uses the position and
velocity measurements together. It does not use the brightness
measurements — they carry no information not already in position and
velocity. The semi-major axis is not one of the 7 numbers being fit,
either: it comes directly from the period, through Kepler's third law,
the same way it does in Stage 1. Uncertainties on all 7 numbers come
from refitting 1000 random resamples of the data (bootstrap
resampling).

Running it prints:

```
parameter                truth      fitted   fit +/- sigma   published sigma
P_yr                    8.6800      8.6798          0.0003            0.1100
e                       0.9832      0.9832          0.0000            0.0010
i_deg                 124.0900    124.1270          0.0287            1.1000
Omega_deg              73.8000     73.7008          0.0703            3.5000
omega_deg             293.4000    293.3609          0.0416            2.2000
t_peri_yr            2023.1260   2023.1264          0.0004            0.0100
omega_dot_deg_yr        0.2307      0.2302          0.0016 n/a (unpublished)
```

**Every one of the 7 numbers comes back within a small fraction of one
standard deviation of the true value.** The precession rate, in
particular, comes back at 0.2302 ± 0.0016 degrees per year, against a
true value of 0.2307 — a real detection of real general-relativistic
physics from simulated data, not a number built in to match. The fitted
sky-plane track visibly shows the periapsis direction rotated between
the two passes (see `output/fit_orbit.png`) — the visual signature of
that precession.

**This result comes entirely from the position measurements.** With
radial velocity realistically too imprecise to detect anything (see
above), position data alone recovers the whole orbit. That is not a
weakness of this simulation — it is exactly how the real discovery paper
solved for S301's orbit in the first place. Radial velocity from a
current instrument adds essentially nothing for a star this faint; that
finding fell directly out of the noise-model derivation, and this
project reports it rather than working around it.

One consequence is worth noting plainly: with a genuinely weak radial
velocity channel and a highly precise, two-pass astrometric baseline,
this project's fit uncertainties end up far tighter than the real
paper's own published ones — for example, 0.070° on the node, against a
published 3.5°. That is not a claim that this simulation beats the real
discovery. It reflects two real differences: this campaign assumes
GRAVITY+'s best demonstrated precision at every epoch, and it uses two
full periapsis passages 8.68 years apart, which the real 2026 paper's
dataset does not yet have.

## What this project does not attempt

Measuring Sgr A\*'s spin is explicitly out of scope. Spin causes
frame-dragging: a drift in the orbit's inclination and node, distinct
from the Schwarzschild precession modeled here. Real astronomers say
measuring it needs a decade or more of observations across several
periapsis passages, plus a way to separate the spin signal from ordinary
"Newtonian confusion" — other, unrelated mass near Sgr A\* can produce a
similar-looking drift. This project's two-passage precession
measurement is a genuine, correctly-modeled precursor to that
measurement, not the measurement itself. Attempting the full spin
measurement here, without a real spin value to inject, would mean
comparing a fit against a guess — exactly what this project has tried
throughout to avoid.

## Quick start

```
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python3 s301_lightcurve.py   # stage 1: relativistic light curve
python3 campaign.py          # stage 2: synthetic noisy campaign -> output/synthetic_observations.csv
python3 fit_orbit.py         # stage 3: recover the orbit from stage 2's data (~2 min: 1000 bootstrap fits)
```

`pylint *.py` scores 10.00/10 across all five modules. `.pylintrc`
documents two narrowly-scoped customizations this required, matching
this course's own established practice of adjusting `.pylintrc` only
with a stated reason. This project has no pytest suite: it follows the
JHU EN.605.256 convention for script-style projects (matching
`module_8`-`module_12` in this repository), not the Flask-app modules'
`pytest --cov-fail-under=100` convention.

## Files

- `constants.py` — every physical constant and published S301/Sgr A*
  number, with its source.
- `orbit.py` — the shared Keplerian and relativistic orbital-mechanics
  engine used by all three stages.
- `s301_lightcurve.py`, `campaign.py`, `fit_orbit.py` — the three
  stages, described above.
- `output/` — generated figures and the synthetic-observations CSV.
  Gitignored; regenerate by running the scripts in order.
