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
The two-passage Schwarzschild-precession extension described in
/home/jcatanz/.claude/plans/virtual-knitting-stearns.md is **complete and
verified** (re-confirmed 2026-09-05): `orbit.py`/`campaign.py`/`fit_orbit.py`
all have the additive `omega_dot` machinery, `campaign.py` produces 128
epochs across both passages, and `fit_orbit.py`'s bootstrap fit recovers
`omega_dot` = 0.2302 +/- 0.0016 deg/yr against a true 0.2307 — matching the
numbers already written into README.md. `output/fit_orbit.png` shows the
expected two rotated loops. No known outstanding work on this extension.

Done: initial git commit made (`59a4c7c`). The OneDrive folder's two
reference PDFs were identified — arXiv:2607.12664 ("Discovery of a star
sensitive to the spin of Sgr A*", Nature, already cited) and arXiv:2607.24931
("S301 and friends: Measuring the spin of Sgr A*", Piran et al.) — and the
latter is now cited in README's "What this project does not attempt"
section (commit `07e8a69`).

All plots now label their instrument (GRAVITY+/VLTI for astrometry+
photometry, ERIS/VLT for RV) in the title/legend, and the sky-track path
in `fit_orbit.py`'s `_plot_sky_track` is epoch-colored (LineCollection)
rather than flat gray, so the real precession-induced self-crossing near
apoapsis reads as "two different times" instead of a stray artifact line.
Also fixed a real inconsistency: fit_orbit.py's docstring wrongly
attributed photometry to ERIS instead of GRAVITY+ (matches campaign.py/
README, which were already correct).

Campaign-parameter optimization (a separate, NOT-yet-done task, found by
re-reading the preserved session transcript): the prior session diagnosed
but never implemented a real bug — `campaign.py` samples epochs uniformly
year-round, ignoring that Sgr A* is only observable from Paranal ~March-
September — plus two unexplored trade-offs (cadence allocation given
e=0.9832's slow apoapsis/fast periapsis asymmetry; RV integration-time vs.
cadence). The session ended by asking the user whether to fix the seasonal
bug now or fold it into the queued-up spin-forecast (Lense-Thirring)
extension; unanswered when the crash happened. Ask the user before
proceeding on this.

Next session: check with the user for new direction (e.g. the campaign-
optimization fix above, the spin-forecast extension, or further article
polish).

## History
The prior session ended due to a platform-side error unrelated to the
article's content — not a code or writing problem. A full transcript of
that session is preserved at
/mnt/c/Users/jcata/OneDrive/Documents/astrophysics/S301-project/s301_session_FULL.jsonl
if a specific past decision or piece of reasoning ever needs to be checked —
read it only if something specific requires it, not as a first move.
