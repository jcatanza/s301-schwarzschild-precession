# S301 precession feasibility study

Simulation code and manuscript for *Catching Up with the Fastest Star in
the Galaxy: A Two-Passage GRAVITY+ Campaign to Detect S301's Schwarzschild
Precession* (J. Catanzarite, 2026; `article.tex`).

S301 is the fastest known star in the Galaxy: it orbits Sagittarius A\*
every 8.68 yr on an e = 0.983 orbit that brings it within 136
Schwarzschild radii of the black hole at more than 8% of the speed of light.
It was discovered by the GRAVITY/VLTI collaboration ("Discovery of a star
sensitive to the spin of Sgr A\*", Nature 2026,
[arXiv:2607.12664](https://arxiv.org/abs/2607.12664)); the companion paper
"S301 and friends" ([arXiv:2607.24931](https://arxiv.org/abs/2607.24931))
develops the multi-star spin-forecast methodology. This project asks a
narrower, prior question: with the instrument that exists today, can a
two-passage GRAVITY+ astrometric campaign (2028.5-2041.7, 131 epochs)
measure S301's Schwarzschild (1PN) apsidal precession -- 0.2307 deg/yr,
about 2.0 deg per orbit -- and what does the systematic error budget for
that measurement look like? Everything here is a synthetic-data study: the
published Solution 1 orbit is injected as truth, realistic noise is added,
and the orbit is recovered blind.

## What it produces

- `results/<script>.json` -- every number each analysis script computes
  that the paper quotes (bootstrap uncertainties, systematic biases,
  cadence comparisons, ...).
- `numbers.tex` -- LaTeX macros generated from those JSONs by
  `make_numbers.py`; `article.tex` `\input`s it.
- `output/` -- the manuscript's figures (`s301_lightcurve.png`,
  `fit_orbit.png`, `s301_field_context.png`) and the rosette animation
  `s301_precession.mp4`, plus untracked scratch products (synthetic CSVs).
- `article.pdf` -- compiled with TeX Live in Docker.
- `arxiv/` -- the flattened, ready-to-upload submission package.

## Quick start

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

pytest -q                  # ~30 s: Kepler solver, physics values, pipeline smoke test
./run_all.sh --skip-slow   # ~15 min: every result except the cadence search and movie
./run_all.sh               # ~1 h: everything, then article.pdf, then arxiv/
```

`run_all.sh` activates `.venv/` if present, prints elapsed time per step,
and skips the PDF and `arxiv/` stages with a notice if `docker` is not
installed. `precession_movie.py` needs `ffmpeg`.

## Scripts and what they feed

| Script | Manuscript | Approx. runtime |
|---|---|---|
| `s301_lightcurve.py` | Sec. 3.1 light curve; Fig. 3; validates a, r_p, v_p against the published solution | seconds |
| `campaign.py` | Sec. 3.2 campaign design; Table 1 context; writes `output/synthetic_observations.csv` | ~2 min (Fisher design search) |
| `fit_orbit.py` | Sec. 3.4 / Sec. 4; **Table 3** (7-parameter fit + 1000-resample bootstrap); Fig. 4 | ~3 min |
| `optimal_design.py` | Sec. 3.2; **Table 2** (cadence designs, heuristic vs. Fisher-optimal) | ~25 min (slow) |
| `cadence_alternatives.py` | Sec. 3.2 / 5; Monte Carlo cross-check of the bootstrap | ~15 min (slow) |
| `photometry_checks.py` | Sec. 3.3 photometric cross-check | ~2 min |
| `rv_channel_test.py` | Sec. 3.2: adding an RV channel changes nothing | ~3 min |
| `precision_sensitivity.py` | Sec. 4: result at 207 uas (achieved, fiducial) vs. 100 uas (forecast) | ~3 min |
| `correlated_noise_test.py` | Sec. 5: per-run common-mode and slowly varying reference systematics | ~5 min |
| `reference_frame_error.py` | Sec. 4.2 reference-frame tie; Table 4 | ~2 min |
| `extended_mass_error.py` | Sec. 4.1 extended-mass (Newtonian) confusion; Table 4 | ~1 min |
| `confusion_error.py` | Sec. 4.3 source confusion / crowding; Table 4 | ~1 min |
| `lensing_error.py` | Sec. 4 weak-lensing astrometric shift; Table 4 | ~1 min |
| `roemer_delay.py` | Sec. 4.5 Roemer delay: bias if ignored, cost of including | ~3 min |
| `spin_contamination.py` | Sec. 4.4 Lense-Thirring contamination bound | seconds |
| `mass_distance_test.py` | Sec. 5: fit with M and R0 free under GRAVITY 2022 priors | ~2 min |
| `solar_conjunction.py` | Sec. 3.2: passage-1 conjunction blackout | ~1 min |
| `precession_movie.py` | Fig. 2 caption; `output/s301_precession.mp4` (arXiv ancillary) | ~3 min (slow) |
| `make_context_figure.py` | Fig. 1b (`s301_field_context.png`); skips itself if the ESO source image is absent | seconds |
| `make_numbers.py` | `numbers.tex` from `results/*.json` | seconds |

Shared modules: `orbit.py` (Kepler solve, relativistic state, Roemer
emission-time solve, precession rates), `photometry.py` (band-integrated
K-band flux model), `constants.py` (physical constants and the published
elements; `TRUTH` = Solution 1), `results_io.py` (results registry).
Runtimes are for a laptop core; the bootstrap-heavy steps scale with
`fit_orbit.N_BOOTSTRAP`.

## Numbers: `results/*.json` -> `numbers.tex` -> `article.tex`

Every script ends by calling `results_io.write_results(name, {...})`, which
writes `results/<name>.json` with each value and its display format.
`make_numbers.py` turns all of those into `\newcommand` macros
(e.g. `results/fit_orbit.json` key `sigma_omega_dot_deg_yr` becomes
`\fitorbitSigmaOmegaDotDegYr`), and `article.tex` uses the macros instead
of typed numbers. The point is drift-proofing: a referee-requested change
to the model (the Roemer delay becoming default-on, the photometric model
changing from a bolometric proxy to a band-integrated one) propagates to
every quoted value on the next `run_all.sh`, and the paper cannot quote a
number the code did not produce. `results/` and `numbers.tex` are tracked
so the record is auditable in the repository history.

## Reproducibility notes

- **Seeds.** `campaign.py` and `fit_orbit.py` use `RNG_SEED = 301`
  (`numpy.random.default_rng`), for the injected noise and the bootstrap
  resampling respectively; reruns reproduce `results/*.json` exactly under
  the pinned `requirements.txt`.
- **PDF.** `article.pdf` and `arxiv/` are built with
  `texlive/texlive:latest` in Docker (pdflatex, bibtex, pdflatex x2), so
  the TeX environment does not depend on the host.
- **Astrometric precision.** The campaign adopts 207 uas per epoch, the
  precision *achieved* on S301 in the discovery data. The discovery
  paper's own GRAVITY+ forecast of 100 uas is the optimistic case;
  `precision_sensitivity.py` reruns the headline fit at both so the paper
  can state the result either way.
- **Correlated systematics.** `correlated_noise_test.py` adds
  per-observing-run common-mode offsets and a slowly varying reference
  term on top of the white noise and measures the precession precision
  from independent realizations; `error_budget.py` folds the excess into
  the bottom line.
- **Roemer delay** is part of both the truth and the fitted model
  (`orbit.orbit_state_observed`, default on), as in the discovery paper's
  own fit; `roemer_delay.py` is the only place it is switched off, to
  measure the bias of ignoring it.
- **Season filter.** Sgr A\* is observable from Paranal only roughly
  March-September; every synthetic epoch respects that. Passage 1's
  periapsis (2031-10-22) falls just outside the season, and the campaign
  reports that gap rather than patching around it.
- The four cadence designs in Table 2 are deterministic given the seed;
  the "naive" and "robust" Fisher designs are kept as documented failures.

## Deliberately not modeled

- Black-hole spin (Kerr metric, Lense-Thirring nodal precession).
  `spin_contamination.py` only bounds how an aligned-spin apsidal term
  would contaminate the fitted rate; the multi-star confusion-subtraction
  method of arXiv:2607.24931 is out of scope.
- Strong-field ray tracing. `lensing_error.py` quantifies the weak
  point-lens astrometric deflection as a systematic; it is not applied in
  the model.
- A radial-velocity channel. ERIS/SPIFFIER precision scaled from the real
  S2 result is ~1600 km/s at m_K = 19.3, so RV was never proposed;
  `rv_channel_test.py` confirms it would not tighten the fit.
- Stellar effective temperature: T_eff = 7300 K is a spectral-type table
  value (not published), affecting only the photometric channel.
- Weather, instrument downtime, and any sub-season scheduling detail
  beyond the visibility window.

## Tests

`pytest -q` runs `tests/` in well under 90 s with no network or Docker:
the Kepler solver on its documented divergent case and a 2000-pair random
sweep, the physics values the paper leans on (precession rate,
Lense-Thirring rate, light-curve amplitude, band-integrated vs. proxy
photometry), and a campaign -> single-fit smoke test including the Roemer
emission-time bound. `pylint --rcfile=.pylintrc *.py` is expected to
report 10.00.

## Licenses

- Code: MIT (`LICENSE`).
- Manuscript, figures and movie: CC BY 4.0 (`LICENSE-MANUSCRIPT.md`),
  with the ESO and EHT images credited as stated there.

## How to cite

Placeholder until the preprint is posted:

```bibtex
@article{catanzarite2026s301,
  author  = {Catanzarite, J.},
  title   = {Catching Up with the Fastest Star in the Galaxy: A Two-Passage
             GRAVITY+ Campaign to Detect S301's Schwarzschild Precession},
  year    = {2026},
  journal = {arXiv e-prints},
  note    = {arXiv:XXXX.XXXXX}
}
```
