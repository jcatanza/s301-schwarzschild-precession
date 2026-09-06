"""
Combine the error-budget terms into one honest bottom-line precision.

Two nuisance sets cost precision when marginalized -- the Sgr A*
reference zero-point offset (reference_frame_error.py) and the mass and
distance priors (mass_distance_test.py). Multiplying their separate
inflation factors assumes they are independent, which need not hold, so
this script fits all of them AT ONCE: an 11-parameter model (7 orbital
parameters + x0, y0 offset with the fiducial prior + m = M/M_pub and
d = R0/R0_pub with the GRAVITY 2022 priors) and bootstraps it. That
sigma is the marginalized random+nuisance precision.

The bounded, non-fitted systematics (extended mass at the conservative
3-sigma bound, lensing, confusion) are then added in quadrature as bias
magnitudes. Roemer delay contributes nothing: it is part of the model.
Spin (Lense-Thirring) is NOT included: it has no bound to inject and is
reported separately as the reason the test is spin-agnostic, not
spin-free.

Results go to results/error_budget.json.
"""

import numpy as np
from scipy.optimize import least_squares

import constants as k
import fit_orbit
import mass_distance_test as md
import orbit
import reference_frame_error as rf
import results_io

N_BOOT = 300


def model_11(params, epochs):
    """RA/Dec (mas) for orbit + reference offset + mass/distance ratios."""
    orbit_params = params[:7]
    x0_mas, y0_mas, m_ratio, d_ratio = params[7:]
    state = fit_orbit.model_state(orbit_params, epochs, grav_param=k.GM_BH * m_ratio)
    ra_mas, dec_mas = orbit.sky_offset_mas(state, k.D_OBS * d_ratio)
    return ra_mas + x0_mas, dec_mas + y0_mas


def residuals_11(params, epochs, ra_obs, dec_obs, sigma_ra, sigma_dec, prior_offset_mas):
    """Astrometric residuals plus Gaussian priors on all four nuisances."""
    ra_model, dec_model = model_11(params, epochs)
    return np.concatenate([
        (ra_model - ra_obs) / sigma_ra,
        (dec_model - dec_obs) / sigma_dec,
        params[7:9] / prior_offset_mas,
        [(params[9] - 1.0) / md.SIGMA_M_FRAC, (params[10] - 1.0) / md.SIGMA_R0_FRAC],
    ])


def fit_11(epochs, ra_obs, dec_obs, sigma_ra, sigma_dec, x0_11, prior_offset_mas):
    """Levenberg-Marquardt fit of the 11-parameter model."""
    lo = np.concatenate([fit_orbit.BOUNDS_LO, [-rf.OFFSET_BOUNDS_MAS] * 2, [0.8, 0.9]])
    hi = np.concatenate([fit_orbit.BOUNDS_HI, [rf.OFFSET_BOUNDS_MAS] * 2, [1.2, 1.1]])
    return least_squares(residuals_11, x0=x0_11, bounds=(lo, hi),
                         args=(epochs, ra_obs, dec_obs, sigma_ra, sigma_dec, prior_offset_mas)).x


def main():
    """Combined nuisance fit, then quadrature with the bounded systematics."""
    data = fit_orbit.load_observations()
    base = results_io.read_results("fit_orbit")
    sigma_base = base["sigma_omega_dot_deg_yr"]
    idx_od = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")
    prior_offset_mas = rf.OFFSET_PRIOR_UAS / 1000.0

    x0_11 = np.concatenate([fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION, [0, 0, 1, 1]])
    best = fit_11(data["epoch_yr"], data["ra_offset_mas"], data["dec_offset_mas"],
                  data["sigma_ra_mas"], data["sigma_dec_mas"], x0_11, prior_offset_mas)
    rng = np.random.default_rng(fit_orbit.RNG_SEED)
    n = len(data["epoch_yr"])
    samples = np.zeros((N_BOOT, 11))
    for k_iter in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        samples[k_iter] = fit_11(data["epoch_yr"][idx], data["ra_offset_mas"][idx],
                                 data["dec_offset_mas"][idx], data["sigma_ra_mas"][idx],
                                 data["sigma_dec_mas"][idx], best, prior_offset_mas)
    sig = samples.std(axis=0)
    sigma_combined = sig[idx_od]
    print(f"11-parameter (orbit + offset + M + R0) bootstrap: omega_dot = {best[idx_od]:.4f} "
          f"+/- {sigma_combined:.5f} deg/yr ({sigma_combined / sigma_base:.2f}x Table 3)")

    extmass = results_io.read_results("extended_mass_error")
    lensing = results_io.read_results("lensing_error")
    confusion = results_io.read_results("confusion_error")
    bias_extmass = abs(extmass["bias_rate_three"])
    bias_lensing = abs(lensing["omega_dot_bias"])
    # Confusion enters as extra per-epoch scatter; express its effect as the
    # extra variance it adds to the random term.
    extra_confusion = sigma_combined * np.sqrt(max(confusion["sigma_ratio"] ** 2 - 1.0, 0.0))
    bottom_line = np.sqrt(sigma_combined ** 2 + bias_extmass ** 2 + bias_lensing ** 2 + extra_confusion ** 2)
    print(f"Bounded systematics in quadrature: extended mass {bias_extmass:.5f}, lensing {bias_lensing:.6f}, "
          f"confusion {extra_confusion:.5f} -> bottom line {bottom_line:.5f} deg/yr "
          f"({bottom_line / sigma_base:.2f}x Table 3)")

    results = {
        "n_boot": N_BOOT,
        "combined_omega_dot": (best[idx_od], ".4f"),
        "combined_sigma": (sigma_combined, ".4f"),
        "combined_ratio": (sigma_combined / sigma_base, ".2f"),
        "bias_extmass": (bias_extmass, ".5f"),
        "bias_lensing": (bias_lensing, ".6f"),
        "extra_confusion": (extra_confusion, ".5f"),
        "bottom_line": (bottom_line, ".4f"),
        "bottom_line_ratio": (bottom_line / sigma_base, ".2f"),
        "bottom_line_pct_of_signal": (bottom_line / fit_orbit.TRUE_OMEGA_DOT_DEG_YR * 100, ".2f"),
        "offset_x_uas": (best[7] * 1000, ".1f"), "offset_y_uas": (best[8] * 1000, ".1f"),
        "m_ratio": (best[9], ".4f"), "d_ratio": (best[10], ".4f"),
    }
    for j, name in enumerate(fit_orbit.PARAM_NAMES):
        results[f"combined_sigma_{name}"] = (sig[j], ".5f")
        results[f"combined_ratio_{name}"] = (sig[j] / base[f"sigma_{name}"], ".2f")
    results_io.write_results("error_budget", results)


if __name__ == "__main__":
    main()
