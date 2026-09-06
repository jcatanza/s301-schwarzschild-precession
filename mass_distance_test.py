"""
Does holding Sgr A*'s mass and distance fixed over-constrain the fit?

fit_orbit.py derives the semi-major axis from (GM, P) with GM and R0
fixed at the discovery paper's values, so the angular scale of the orbit
is pinned by the period alone. Real S-star fits float both. This script
refits with M and R0 as free parameters under Gaussian priors from
GRAVITY Collaboration (2022, A&A 657, L12): M = (4.297 +/- 0.012) x 10^6
Msun and R0 = 8277 +/- 9 pc statistical, with ~40,000 Msun and ~30 pc
systematic, combined in quadrature here (1.0% and 0.38%).

Astrometry alone constrains only the combination GM / R0^3 (the angular
semi-major axis is a/R0 with a^3 = GM P^2 / 4 pi^2), the classic
mass-distance degeneracy that radial velocities normally break; the
priors are what make M and R0 separately identifiable in this fit. The
parameters are the dimensionless ratios m = M/M_pub and d = R0/R0_pub.

Reports each parameter's bootstrap precision with M and R0 free versus
Table 3's fixed-scale values, and whether omega_dot moves. Results go to
results/mass_distance_test.json.
"""

import numpy as np
from scipy.optimize import least_squares

import constants as k
import fit_orbit
import orbit
import results_io

M_STAT_MSUN, M_SYS_MSUN = 0.012e6, 0.040e6
R0_STAT_PC, R0_SYS_PC = 9.0, 30.0
SIGMA_M_FRAC = np.hypot(M_STAT_MSUN, M_SYS_MSUN) / 4.297e6
SIGMA_R0_FRAC = np.hypot(R0_STAT_PC, R0_SYS_PC) / 8277.0
N_BOOT = 300


def model_observables_md(params9, epochs_yr):
    """RA/Dec (mas) for 7 orbital parameters plus m = M/M_pub, d = R0/R0_pub."""
    orbit_params, m_ratio, d_ratio = params9[:7], params9[7], params9[8]
    state = fit_orbit.model_state(orbit_params, epochs_yr, grav_param=k.GM_BH * m_ratio)
    return orbit.sky_offset_mas(state, k.D_OBS * d_ratio)


def residuals_md(params9, epochs, ra_obs, dec_obs, sigma_ra, sigma_dec):
    """Astrometric residuals plus Gaussian priors on m and d."""
    ra_model, dec_model = model_observables_md(params9, epochs)
    return np.concatenate([
        (ra_model - ra_obs) / sigma_ra,
        (dec_model - dec_obs) / sigma_dec,
        [(params9[7] - 1.0) / SIGMA_M_FRAC, (params9[8] - 1.0) / SIGMA_R0_FRAC],
    ])


def fit_md(epochs, ra_obs, dec_obs, sigma_ra, sigma_dec, x0_9):
    """Levenberg-Marquardt fit of the 9-parameter (orbit + M + R0) model."""
    lo = np.concatenate([fit_orbit.BOUNDS_LO, [0.8, 0.9]])
    hi = np.concatenate([fit_orbit.BOUNDS_HI, [1.2, 1.1]])
    return least_squares(residuals_md, x0=x0_9, bounds=(lo, hi),
                         args=(epochs, ra_obs, dec_obs, sigma_ra, sigma_dec)).x


def main():
    """Fit with M and R0 free (priors), bootstrap, compare to Table 3."""
    data = fit_orbit.load_observations()
    base = results_io.read_results("fit_orbit")
    x0_9 = np.concatenate([fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION, [1.0, 1.0]])
    best = fit_md(data["epoch_yr"], data["ra_offset_mas"], data["dec_offset_mas"],
                  data["sigma_ra_mas"], data["sigma_dec_mas"], x0_9)

    rng = np.random.default_rng(fit_orbit.RNG_SEED)
    n = len(data["epoch_yr"])
    samples = np.zeros((N_BOOT, 9))
    for k_iter in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        samples[k_iter] = fit_md(data["epoch_yr"][idx], data["ra_offset_mas"][idx],
                                 data["dec_offset_mas"][idx], data["sigma_ra_mas"][idx],
                                 data["sigma_dec_mas"][idx], best)
    sig = samples.std(axis=0)

    print(f"Priors: sigma(M)/M = {SIGMA_M_FRAC * 100:.2f}%, sigma(R0)/R0 = {SIGMA_R0_FRAC * 100:.2f}%")
    print(f"{'parameter':<18}{'fixed-scale sigma':>20}{'M,R0 free sigma':>18}{'ratio':>8}")
    results = {"sigma_m_pct": (SIGMA_M_FRAC * 100, ".2f"), "sigma_r0_pct": (SIGMA_R0_FRAC * 100, ".2f"),
               "n_boot": N_BOOT}
    for j, name in enumerate(fit_orbit.PARAM_NAMES):
        fixed = base[f"sigma_{name}"]
        print(f"{name:<18}{fixed:>20.5f}{sig[j]:>18.5f}{sig[j] / fixed:>8.2f}")
        results[f"sigma_{name}"] = (sig[j], ".4f")
        results[f"ratio_{name}"] = (sig[j] / fixed, ".2f")
    idx_od = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")
    shift = best[idx_od] - base["fit_omega_dot_deg_yr"]
    print(f"\nomega_dot: {best[idx_od]:.4f} +/- {sig[idx_od]:.4f} deg/yr (shift {shift:+.5f} vs Table 3)")
    print(f"recovered m = {best[7]:.4f} +/- {sig[7]:.4f}, d = {best[8]:.4f} +/- {sig[8]:.4f}")
    results.update({
        "fit_omega_dot": (best[idx_od], ".4f"),
        "omega_dot_shift": (shift, ".5f"),
        "omega_dot_shift_nsigma": (abs(shift) / base["sigma_omega_dot_deg_yr"], ".2f"),
        "m_ratio": (best[7], ".4f"), "m_ratio_sigma": (sig[7], ".4f"),
        "d_ratio": (best[8], ".4f"), "d_ratio_sigma": (sig[8], ".4f"),
    })
    results_io.write_results("mass_distance_test", results)


if __name__ == "__main__":
    main()
