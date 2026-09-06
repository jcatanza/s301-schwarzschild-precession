# S301 feasibility study -- project note

Simulation code and manuscript for a feasibility study of measuring the
Schwarzschild apsidal precession of S301 (the fastest known star orbiting
Sgr A*, discovered by GRAVITY/VLTI; arXiv:2607.12664, companion
arXiv:2607.24931) with a synthetic two-passage GRAVITY+/VLTI astrometric
campaign, plus a systematic error budget. Single author: J. Catanzarite.
`README.md` is the project overview; read it first.

## Ground rules

- **`article.tex` is the canonical manuscript.** `article.pdf` is its
  compiled form (Docker TeX Live; see `run_all.sh`). There is no other
  copy of the text: `README.md` is a project README, not the paper, and
  `arxiv/` is a generated flattened package -- never edit it by hand.
- **The paper never quotes a number the code did not produce.** Every
  analysis script writes `results/<script>.json` through
  `results_io.write_results`; `make_numbers.py` turns all of them into
  `numbers.tex` macros that `article.tex` `\input`s. To change a number,
  change the code and re-run; do not hand-edit `numbers.tex` or type a
  number into `article.tex`.
- **`run_all.sh`** runs the whole pipeline in dependency order, then the
  PDF, then regenerates `arxiv/`. `./run_all.sh --skip-slow` skips
  `optimal_design.py`, `cadence_alternatives.py` and `precession_movie.py`.
- **Tests:** `pytest -q` (`tests/`; < 90 s, no network, no Docker). Run
  them and `pylint --rcfile=.pylintrc *.py` (expected 10.00) before
  committing code.
- **Fixed seeds:** `RNG_SEED = 301` in `campaign.py` and `fit_orbit.py`;
  reruns reproduce `results/*.json` exactly.
- `constants.py` is the single source of truth for physical constants and
  the published S301 elements (`TRUTH` = Solution 1). `orbit.py` is a
  standalone physics engine that deliberately does not import it.

## Scripts (see README.md for the table mapping each to the manuscript)

- Core pipeline: `s301_lightcurve.py` -> `campaign.py` -> `fit_orbit.py`
  -> `optimal_design.py`; shared modules `orbit.py`, `photometry.py`,
  `constants.py`, `results_io.py`.
- Robustness tests: `cadence_alternatives.py`, `photometry_checks.py`,
  `rv_channel_test.py`, `precision_sensitivity.py`,
  `mass_distance_test.py`, `solar_conjunction.py`.
- Systematics (error budget): `reference_frame_error.py`,
  `extended_mass_error.py`, `confusion_error.py`, `lensing_error.py`,
  `roemer_delay.py`, `spin_contamination.py`.
- Figures: `precession_movie.py` (mp4, ffmpeg), `make_context_figure.py`
  (needs the ESO source image in `output/`; skips itself if absent).
- Manuscript: `make_numbers.py`, `article.tex`, `references.bib`,
  `run_all.sh`, `arxiv/`.

## Conventions that trip people up

- Orbit frame: +Z toward the observer; `v_los` positive = receding.
- `omega_deg` in the fit is the argument of periapsis at the fitted
  `t_peri` (`orbit_state_precessing`'s `t_ref`); the semi-major axis is
  derived from (GM, P), never fit.
- The Roemer delay is on by default (`orbit.orbit_state_observed`,
  `light_time=True`) in both truth and model; `roemer_delay.py` is the
  only place `light_time=False` is used, to quantify the bias.
- There is no radial-velocity channel in the campaign (feasibility
  rejection documented in `campaign.py`); `rv_channel_test.py` verifies
  it would not help.
- Sgr A* is only observable from Paranal roughly March-September;
  `campaign.py` filters every epoch to that season, and passage 1's
  periapsis (2031-10-22) falls just outside it.

## Editorial pass

For a writing/structure review of `article.tex`, use the user-level skill
`/manuscript-review` (`~/.claude/skills/manuscript-review/SKILL.md`).
