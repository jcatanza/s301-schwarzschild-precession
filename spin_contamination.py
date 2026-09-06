"""
Quantifies the one gap in this project's error budget that is different
IN KIND from the three already-quantified systematics (extended-mass
confusion, astrometric reference-frame systematics, source confusion/
crowding): Lense-Thirring (frame-dragging) contamination of the apsidal
precession rate itself.

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
own docstring), compared against the Schwarzschild signal
(fit_orbit.TRUE_OMEGA_DOT_DEG_YR) and this project's fitted omega_dot
precision: the photon-noise-only value from results/fit_orbit.json
(Table 3) and the reference-systematic-marginalized value from
reference_frame_error.py (results/reference_frame_error.json, key
joint_sigma_fid) when that file exists, falling back to a documented
placeholder otherwise.

Results go to results/spin_contamination.json.
"""

import os

import numpy as np

import constants as k
import fit_orbit
import orbit
import results_io

PERIOD_SEC = k.TRUTH["P_yr"] * k.year
SMA_M = orbit.semi_major_axis_from_period(k.GM_BH, PERIOD_SEC)
ECC = k.TRUTH["e"]

SCHWARZSCHILD_SIGNAL_DEG_YR = fit_orbit.TRUE_OMEGA_DOT_DEG_YR
# Used only if reference_frame_error.py has not yet written its results
# file; the real reference-marginalized sigma is read from there.
SIGMA_REF_FALLBACK = 0.0011

SPIN_CASES = (
    ("chi50", "conservative sub-maximal", 0.5),
    ("chi90", "Daly et al. (2024) outflow-method measurement", 0.90),
    ("chi100", "discovery paper's own fiducial (chi=1, aligned)", 1.0),
)
CHI_VALUES_STR = "0.5, 0.90, 1.0"   # the SPIN_CASES values as the manuscript quotes them


def lt_rate_deg_yr(spin_chi, xi_deg=0.0):
    """Lense-Thirring apsidal precession rate (deg/yr) for the given
    dimensionless spin and orbit-spin misalignment angle."""
    rate = orbit.lense_thirring_apsidal_rate(k.GM_BH, SMA_M, ECC, PERIOD_SEC, spin_chi, xi_deg)
    return np.degrees(rate) * k.year


def load_sigmas():
    """(photon-noise-only sigma, reference-systematic-marginalized sigma)
    on omega_dot in deg/yr. The first comes from results/fit_orbit.json;
    the second from results/reference_frame_error.json if it exists (it
    is produced by a separate script and may not have been run yet), else
    SIGMA_REF_FALLBACK with a printed note."""
    sigma_noise = results_io.read_results("fit_orbit")["sigma_omega_dot_deg_yr"]
    path = os.path.join(results_io.RESULTS_DIR, "reference_frame_error.json")
    sigma_ref = None
    if os.path.exists(path):
        sigma_ref = results_io.read_results("reference_frame_error").get("joint_sigma_fid")
    if sigma_ref is None:
        print(f"NOTE: {path} (key joint_sigma_fid) not available -- reference_frame_error.py has not "
              f"been run; using the fallback reference-marginalized sigma {SIGMA_REF_FALLBACK} deg/yr.\n")
        sigma_ref = SIGMA_REF_FALLBACK
    return sigma_noise, sigma_ref


def main():
    """Print the Lense-Thirring contamination table for each real spin
    case, write the results, and close with why this is not a
    subtractable term."""
    sigma_noise, sigma_ref = load_sigmas()
    print("Lense-Thirring apsidal contamination of S301's measured omega_dot")
    print("(xi=0: spin aligned with orbital angular momentum -- maximal, purely-apsidal case)")
    print(f"Schwarzschild signal {SCHWARZSCHILD_SIGNAL_DEG_YR:.4f} deg/yr; sigma(omega_dot) = "
          f"{sigma_noise:.4f} (photon noise only), {sigma_ref:.4f} (reference systematics marginalized)\n")
    print(f"{'case':<45}{'chi':>6}{'deg/yr':>12}{'% of signal':>14}"
          f"{'sigma (noise)':>16}{'sigma (ref-marg.)':>20}")
    results = {"chi_values": CHI_VALUES_STR, "sigma_ref_used": (sigma_ref, ".4f")}
    for tag, label, chi in SPIN_CASES:
        rate = lt_rate_deg_yr(chi)
        pct_signal = rate / SCHWARZSCHILD_SIGNAL_DEG_YR * 100
        nsigma_noise = rate / sigma_noise
        nsigma_ref = rate / sigma_ref
        print(f"{label:<45}{chi:>6.2f}{rate:>12.5f}{pct_signal:>13.2f}%"
              f"{nsigma_noise:>16.1f}{nsigma_ref:>20.1f}")
        results.update({
            f"lt_rate_{tag}": (rate, ".5f"),
            f"lt_pct_{tag}": (pct_signal, ".2f"),
            f"lt_nsigma_noise_{tag}": (nsigma_noise, ".1f"),
            f"lt_nsigma_ref_{tag}": (nsigma_ref, ".1f"),
        })
    results_io.write_results("spin_contamination", results)

    print("\nUnlike the three bounded systematics elsewhere in this project's error budget, "
          "this is not a subtractable nuisance: Sgr A*'s spin is not robustly measured, so "
          "there is no analogous upper limit to inject and marginalize over. A measured "
          "omega_dot from S301 alone cannot distinguish 'pure GR monopole' from 'GR monopole "
          "plus aligned frame-dragging' -- exactly the degeneracy Piran et al. (2026)'s multi-star "
          "method exists to break, and exactly why this project's Scope paragraph excludes a "
          "spin measurement rather than claiming this pipeline delivers a clean point-mass test.")


if __name__ == "__main__":
    main()
