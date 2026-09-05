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

Not yet done: the repo has no git commits (fresh `.git`, no history) —
confirm with the user whether/when to make an initial commit. Also open:
CLAUDE.md's earlier note to confirm which paper the OneDrive reference PDF
is before citing it further in the article.

## History
The prior session ended due to a platform-side error unrelated to the
article's content — not a code or writing problem. A full transcript of
that session is preserved at
/mnt/c/Users/jcata/OneDrive/Documents/astrophysics/S301-project/s301_session_FULL.jsonl
if a specific past decision or piece of reasoning ever needs to be checked —
read it only if something specific requires it, not as a first move.
