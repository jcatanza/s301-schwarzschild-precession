"""
Quantifies the one gap in this project's error budget that is different
IN KIND from the three already-quantified systematics (extended-mass
confusion, astrometric reference-frame tie, source confusion/crowding):
Lense-Thirring (frame-dragging) contamination of the apsidal precession
rate itself.

The other three systematics are nuisances with a known or boundable
magnitude that can be modeled and subtracted (or shown subdominant).
Black-hole spin is not: for a spin aligned with S301's orbital angular
momentum, frame-dragging produces a PURELY APSIDAL shift (see
orbit.lense_thirring_apsidal_rate's docstring) -- i.e. it adds directly,
indistinguishably in form, to the same omega_dot this project's entire
pipeline fits for. Since Sgr A*'s spin is not robustly measured, there
is no analogous GRAVITY-quoted upper bound to inject the way
extended_mass_error.py used the real 3000 Msun limit. This script
instead shows the SCALE of the possible contamination across the range
of real, published (if contested) spin estimates -- not a bound, a
plausibility check.

Real, quoted spin values used (no recalled/fabricated numbers):

  - chi=1.0 (maximal, aligned): the S301 discovery paper's OWN fiducial
    benchmark (GRAVITY Collaboration 2026, arXiv:2607.12664): "We
    optimistically adopt chi=1 and an orientation approximately aligned
    with the stellar orbital angular momentum" for their own spin-
    detectability forecast. orbit.lense_thirring_apsidal_rate was
    already validated in this project against their quoted value at
    this exact benchmark (0.114 deg/orbit).
  - chi=0.90 +/- 0.06: an independent, real (if debated -- spin
    estimates for Sgr A* vary by method) measurement via the radio/
    X-ray outflow method: Daly et al. (2024, MNRAS 527, 428,
    arXiv:2310.12108), "New Black Hole Spin Values for Sagittarius A*
    Obtained with the Outflow Method."
  - chi=0.5: a conservative, sub-maximal illustrative case, spanning
    the range of estimates found across the spin literature (values as
    low as ~0.5 have also been reported by other methods/epochs).

Method: orbit.lense_thirring_apsidal_rate (already in this project's
physics engine, unchanged) evaluated at S301's real orbital parameters
for xi_deg=0 (spin aligned with the orbital plane -- the maximal-
contamination, purely-apsidal case; a misaligned spin would instead
produce some nodal precession, distinguishable in principle, but that
general case is explicitly not implemented here -- see the function's
own docstring), compared against the Schwarzschild signal and this
project's fitted omega_dot precision (both the photon-noise-only value
from Table 3 and the frame-tie-marginalized value from
instrumental_error.py).
"""

import numpy as np

import constants as k
import orbit

PERIOD_SEC = k.TRUTH["P_yr"] * k.year
SMA_M = orbit.semi_major_axis_from_period(k.GM_BH, PERIOD_SEC)
ECC = k.TRUTH["e"]

SCHWARZSCHILD_SIGNAL_DEG_YR = 0.2307
SIGMA_PHOTON_NOISE_ONLY = 0.0011       # Table 3
SIGMA_FRAME_TIE_MARGINALIZED = 0.0028  # instrumental_error.py, central Plewa et al. case

SPIN_CASES = (
    ("conservative sub-maximal", 0.5),
    ("Daly et al. (2024) outflow-method measurement", 0.90),
    ("discovery paper's own fiducial (chi=1, aligned)", 1.0),
)


def lt_rate_deg_yr(spin_chi, xi_deg=0.0):
    """Lense-Thirring apsidal precession rate (deg/yr) for the given
    dimensionless spin and orbit-spin misalignment angle."""
    rate = orbit.lense_thirring_apsidal_rate(k.GM_BH, SMA_M, ECC, PERIOD_SEC, spin_chi, xi_deg)
    return np.degrees(rate) * k.year


def main():
    """Print the Lense-Thirring contamination table for each real spin
    case, and the closing note on why this is not a subtractable term."""
    print("Lense-Thirring apsidal contamination of S301's measured omega_dot")
    print("(xi=0: spin aligned with orbital angular momentum -- maximal, purely-apsidal case)\n")
    print(f"{'case':<45}{'chi':>6}{'deg/yr':>12}{'% of signal':>14}"
          f"{'sigma (noise)':>16}{'sigma (frame-tie)':>20}")
    for label, chi in SPIN_CASES:
        rate = lt_rate_deg_yr(chi)
        pct_signal = rate / SCHWARZSCHILD_SIGNAL_DEG_YR * 100
        sigma_noise = rate / SIGMA_PHOTON_NOISE_ONLY
        sigma_frame = rate / SIGMA_FRAME_TIE_MARGINALIZED
        print(f"{label:<45}{chi:>6.2f}{rate:>12.5f}{pct_signal:>13.2f}%"
              f"{sigma_noise:>16.1f}{sigma_frame:>20.1f}")

    print("\nUnlike the three bounded systematics elsewhere in this project's error budget, "
          "this is not a subtractable nuisance: Sgr A*'s spin is not robustly measured, so "
          "there is no analogous upper limit to inject and marginalize over. A measured "
          "omega_dot from S301 alone cannot distinguish 'pure GR monopole' from 'GR monopole "
          "plus aligned frame-dragging' -- exactly the degeneracy Piran et al. (2026)'s multi-star "
          "method exists to break, and exactly why this project's Scope paragraph excludes a "
          "spin measurement rather than claiming this pipeline delivers a clean point-mass test.")


if __name__ == "__main__":
    main()
