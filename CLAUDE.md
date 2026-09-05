# S301 Observing Campaign — Project Context

## What this project is
An article + simulation code for a ground-based GRAVITY+/VLTI observing
campaign of S301, the fastest known star in the Milky Way, discovered
orbiting Sagittarius A* (announced August 2026, Nature). The project has
three parts:
1. A relativistic light curve for S301 over one orbit
2. A synthetic observing campaign (noisy, realistic measurements tied to
   the real GRAVITY+/VLTI instrument; NO radial-velocity channel — see
   "RV is dropped" below)
3. An orbit fit recovering the orbital shape and precession from the
   synthetic data alone

## Where things stand
- `README.md` is the article itself — substantially written, current as
  of 2026-09-05.
- Simulation code: `campaign.py`, `orbit.py`, `fit_orbit.py`,
  `optimal_design.py`, `s301_lightcurve.py`, `constants.py`. `pylint`
  scores 10.00/10 across all six.
- Generated figures: `output/fit_orbit.png` (2 panels: sky-track,
  photometry — no RV panels), `output/s301_lightcurve.png`.
- Reference PDFs identified and cited: arXiv:2607.12664 (discovery
  paper) and arXiv:2607.24931 ("S301 and friends", Piran et al.,
  spin-forecast companion paper) — both cited in README.

## What's outstanding

**Complete and verified (2026-09-05 session):**
- Two-passage Schwarzschild-precession extension: `omega_dot` machinery
  in orbit.py/campaign.py/fit_orbit.py, fit recovers it within ~1sigma
  of the real 1PN prediction.
- Every plot labels its instrument in the title/legend; the sky-track
  path is epoch-colored (LineCollection) so the real precession
  self-crossing near apoapsis doesn't read as a rendering artifact.
- Seasonal-visibility bug fixed: `campaign.py` filters every epoch to
  the real Mar-Sep Paranal visibility window. Passage 1's periapsis
  (2031-10-22) falls just outside it, so its dense cluster anchors to
  the nearest visible date (`_season_anchor`) instead — a real,
  reported scheduling gap.
- **Cadence, a real methodological journey** (`optimal_design.py`): a
  naive Fisher/D-optimal design (linearized at the truth) predicted
  huge improvements but catastrophically failed a real bootstrap test
  (piled epochs at 2 extrapolation-extreme clusters, skipped periapsis
  entirely, 4-600x worse in practice — a nonlinear "cycle-slip"
  aliasing pitfall). A pseudo-Bayesian fix (average Fisher over a
  parameter prior) still skipped periapsis (Fisher info can't see
  structural non-identifiability). The winning design adds a hard
  periapsis-coverage floor (~30% of budget, domain knowledge, not
  statistically derived) plus Fisher-optimal allocation of the rest.
  This is `campaign.py`'s actual default `build_epoch_grid()` now
  (`build_epoch_grid_heuristic()` kept only as the baseline it beats).
  Real bootstrap result: beats the old uniform heuristic on every
  parameter, most importantly `omega_dot` (0.0013 vs. 0.0043 deg/yr).
- **RV dropped from the actual campaign entirely.** ERIS/SINFONI-scaled
  RV precision (~1621 km/s vs. a ~15,000 km/s signal, SNR~0.04) means no
  real Time Allocation Committee would grant ERIS time for this
  measurement — so it was never proposed, not simulated-then-discarded.
  `campaign.compute_rv_precision()` stays as the documented feasibility-
  rejection calculation; no RV column exists in the CSV, fit_orbit.py's
  residuals, or optimal_design.py's Fisher computation. Verified directly
  (not just argued): refitting with vs. without a simulated RV channel
  changed no parameter's uncertainty outside bootstrap noise.
  `fit_orbit.py`'s figure is now 2 panels (sky-track, photometry), down
  from 4.
- Checked (WebSearch, 2026-09-05) whether Roman (ruled out: confusion
  floor + coarser resolution, neither fixed by longer integration) or
  JWST/NIRSpec (real resolution is decent and no OH airglow, but real
  Galactic-Center-specific operational problems — saturation, MSA
  leakage, guide-star misidentification — and no citable achieved
  precision for a comparable faint point source) could do better. Real
  future fix is the ELT's HARMONI, but ESO's latest schedule has science
  first light slipping to December 2030 (later than this campaign's
  2028.5 start) — noted in README.

**Still open / deferred:**
- RV integration-time vs. cadence trade-off — moot now that RV isn't
  proposed at all.
- Spin-forecast (Lense-Thirring) extension: `orbit.lense_thirring_
  apsidal_rate` and campaign.py's `extra_omega_dot` hook already exist
  and were validated against the real paper's own quoted number
  (0.114 deg/orbit), but `spin_forecast.py` itself was never created.
  The prior session flagged major overlap with the discovery paper's own
  mock spin forecast and the companion paper's multi-star Newtonian-
  confusion methodology — needs an honest framing discussion with the
  user before building, not just an implementation.
- User requested (2026-09-05) a full ApJ-style journal-article draft
  summarizing the project — in progress/next.

Next session: check with the user for direction among the above.

## History
The prior session ended due to a platform-side error unrelated to the
article's content — not a code or writing problem. A full transcript of
that session is preserved at
/mnt/c/Users/jcata/OneDrive/Documents/astrophysics/S301-project/s301_session_FULL.jsonl
if a specific past decision or piece of reasoning ever needs to be checked —
read it only if something specific requires it, not as a first move.
