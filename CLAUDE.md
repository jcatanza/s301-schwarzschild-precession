# S301 Observing Campaign — Project Context

## What this project is
An article + simulation code for a ground-based VLT/VLTI observing campaign
of S301, the fastest known star in the Milky Way, discovered orbiting
Sagittarius A* (announced August 2026, Nature). The project has three parts:
1. A relativistic light curve for S301 over one orbit
2. A synthetic observing campaign (noisy, realistic measurements tied to
   real instruments — GRAVITY+ for astrometry, ERIS for radial velocity)
3. An orbit fit recovering the orbital shape and precession from the
   synthetic data alone

## Where things stand
- `README.md` (358 lines) is the article itself — substantially written.
- Simulation code: `campaign.py`, `orbit.py`, `fit_orbit.py`,
  `s301_lightcurve.py`, `constants.py`.
- Generated figures: `output/fit_orbit.png`, `output/s301_lightcurve.png`.
- Reference PDF in the project's OneDrive folder — confirm which paper
  this is before citing it further.

## What's outstanding

**Complete and verified:**
- Two-passage Schwarzschild-precession extension (see
  /home/jcatanz/.claude/plans/virtual-knitting-stearns.md): `omega_dot`
  machinery in orbit.py/campaign.py/fit_orbit.py, fit recovers it to
  within ~1sigma of the real 1PN prediction.
- Reference PDFs identified and cited: arXiv:2607.12664 (discovery paper)
  and arXiv:2607.24931 ("S301 and friends", Piran et al., spin-forecast
  companion paper) — both now cited in README.
- Every plot labels its instrument (GRAVITY+/VLTI vs. ERIS/VLT) in the
  title/legend; the sky-track path is epoch-colored (LineCollection) so
  the real precession self-crossing near apoapsis doesn't read as a
  rendering artifact.
- Seasonal-visibility bug fixed: `campaign.py` now filters every epoch to
  the real Mar-Sep Paranal visibility window (`_in_visibility_season`).
  Passage 1's periapsis (2031-10-22) falls just outside that window, so
  its dense cluster is anchored to the nearest visible date
  (`_season_anchor`, ~2031-09-30) instead of periapsis itself — a real,
  reported scheduling gap, not patched around.
- Cadence redesigned to the user's original spec: ~90% of epochs within
  +/-20 days of a periapsis (131 total, 90.1% verified at runtime), not
  spread evenly. This resolves the earlier "cadence allocation" open
  question. Real, reported trade-off: P_yr/e got *more* precise, but the
  angular elements (i, Omega, omega) and omega_dot got *less* precise
  than the old (evenly-spread-per-year) design — physically sensible
  (periapsis-only sampling gives timing leverage but less angular
  diversity to pin down 3D orientation). All 7 parameters still recover
  within ~1.5sigma of truth. See README's Stage 2/3 sections.

**Still open:**
- RV integration-time vs. cadence trade-off — never quantitatively
  tested (see README/transcript). Low priority: unlikely to change the
  "RV is undetectable with current instruments" finding.
- User asked (2026-09-05) whether a Fisher-information/D-optimal design
  would be more principled than the +/-20-day heuristic, given the
  orbit's extreme eccentricity (e=0.9832) and the just-observed
  precision trade-off across parameters. Discussed conceptually; not yet
  implemented — ask the user before building it (real effort: needs a
  per-parameter sensitivity Jacobian across a time grid plus a greedy or
  convex-optimal epoch-selection routine).
- Spin-forecast (Lense-Thirring) extension: `orbit.lense_thirring_
  apsidal_rate` and campaign.py's `extra_omega_dot` hook already exist
  and were validated against the real paper's own quoted number
  (0.114 deg/orbit), but `spin_forecast.py` itself was never created.
  The prior session flagged major overlap with the discovery paper's own
  mock spin forecast and the companion paper's multi-star Newtonian-
  confusion methodology — needs an honest framing discussion with the
  user before building, not just an implementation.

Next session: check with the user for direction among the above.

## History
The prior session ended due to a platform-side error unrelated to the
article's content — not a code or writing problem. A full transcript of
that session is preserved at
/mnt/c/Users/jcata/OneDrive/Documents/astrophysics/S301-project/s301_session_FULL.jsonl
if a specific past decision or piece of reasoning ever needs to be checked —
read it only if something specific requires it, not as a first move.
