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
magnitudes, and so is the EXCESS scatter from correlated systematics
(correlated_noise_test.py: per-run common-mode offsets and a slowly
varying reference term), scaled by the same factor as the marginalized
term. Roemer delay contributes nothing: it is part of the model.
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

# "Unknown unknowns": real GRAVITY orbit fits scatter more than their formal
# errors predict. The S2 Schwarzschild-precession fit has reduced chi^2 =
# 1.5 (GRAVITY Collaboration 2020, A&A 636, L5, arXiv:2004.07187) and the
# four-star fit of GRAVITY Collaboration 2022 (A&A 657, L12,
# arXiv:2112.07478) has chi_r^2 = 2.17. Precision scales linearly with the
# per-epoch error, so inflating every uncertainty by sqrt(chi_r^2) is the
# empirical allowance for everything this budget has not named. The S301
# discovery fit itself has chi^2 = 26 for 34 degrees of freedom, i.e. no
# excess at the 207 uas floor, so the inflation is a conservative envelope.
CHI2R_EMPIRICAL_LOW = 1.5
CHI2R_EMPIRICAL_HIGH = 2.17


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


def bounded_terms(sigma_combined):
    """The non-fitted terms added in quadrature, as sigma-equivalents in
    deg/yr: the 3-sigma extended-mass bias, the lensing bias, and the
    EXCESS scatter from confusion and from correlated systematics. The two
    excesses are inflation factors measured on the random term, applied
    here to the marginalized term."""
    extmass = results_io.read_results("extended_mass_error")
    lensing = results_io.read_results("lensing_error")
    confusion = results_io.read_results("confusion_error")
    corr = results_io.read_results("correlated_noise_test")

    def excess(ratio):
        return sigma_combined * np.sqrt(max(ratio ** 2 - 1.0, 0.0))

    return {
        "extmass": abs(extmass["bias_rate_three"]),
        "lensing": abs(lensing["omega_dot_bias"]),
        "confusion": excess(confusion["sigma_ratio"]),
        "correlated": excess(corr["fid_ratio_omega_dot"]),
    }


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

    terms = bounded_terms(sigma_combined)
    bottom_line_white = np.sqrt(sigma_combined ** 2 + terms["extmass"] ** 2 + terms["lensing"] ** 2
                                + terms["confusion"] ** 2)
    bottom_line = np.sqrt(bottom_line_white ** 2 + terms["correlated"] ** 2)
    print(f"Bounded systematics in quadrature: extended mass {terms['extmass']:.5f}, "
          f"lensing {terms['lensing']:.6f}, confusion {terms['confusion']:.5f} -> white-noise "
          f"bottom line {bottom_line_white:.5f} deg/yr ({bottom_line_white / sigma_base:.2f}x Table 3)")
    print(f"Correlated systematics excess {terms['correlated']:.5f} -> bottom line {bottom_line:.5f} deg/yr "
          f"({bottom_line / sigma_base:.2f}x Table 3, "
          f"{fit_orbit.TRUE_OMEGA_DOT_DEG_YR / bottom_line:.0f} sigma detection)")
    inflation = np.sqrt(CHI2R_EMPIRICAL_HIGH)
    bottom_line_inflated = bottom_line * inflation
    # The unbounded spin term, expressed against the headline precision so
    # the manuscript compares like with like.
    spin = results_io.read_results("spin_contamination")
    spin_nsigma_headline = {tag: abs(spin[f"lt_rate_chi{tag}"]) / bottom_line_inflated
                            for tag in ("90", "100")}
    print(f"Empirical sqrt(chi_r^2) inflation x{inflation:.2f} (chi_r^2 = {CHI2R_EMPIRICAL_HIGH}) -> "
          f"{bottom_line_inflated:.5f} deg/yr, "
          f"{fit_orbit.TRUE_OMEGA_DOT_DEG_YR / bottom_line_inflated:.0f} sigma detection")

    results = {
        "n_boot": N_BOOT,
        "combined_omega_dot": (best[idx_od], ".4f"),
        "combined_sigma": (sigma_combined, ".4f"),
        "combined_ratio": (sigma_combined / sigma_base, ".2f"),
        "bias_extmass": (terms["extmass"], ".5f"),
        "bias_lensing": (terms["lensing"], ".6f"),
        "extra_confusion": (terms["confusion"], ".5f"),
        "extra_correlated": (terms["correlated"], ".5f"),
        "bottom_line_white": (bottom_line_white, ".4f"),
        "bottom_line_white_ratio": (bottom_line_white / sigma_base, ".2f"),
        "bottom_line": (bottom_line, ".4f"),
        "bottom_line_ratio": (bottom_line / sigma_base, ".2f"),
        "bottom_line_pct_of_signal": (bottom_line / fit_orbit.TRUE_OMEGA_DOT_DEG_YR * 100, ".2f"),
        "bottom_line_detection_nsigma": (fit_orbit.TRUE_OMEGA_DOT_DEG_YR / bottom_line, ".0f"),
        "correlated_share_pct": (terms["correlated"] ** 2 / bottom_line ** 2 * 100, ".0f"),
        "chi2r_low": CHI2R_EMPIRICAL_LOW, "chi2r_high": CHI2R_EMPIRICAL_HIGH,
        "inflation_low": (np.sqrt(CHI2R_EMPIRICAL_LOW), ".2f"),
        "inflation_high": (inflation, ".2f"),
        "bottom_line_inflated": (bottom_line_inflated, ".4f"),
        "bottom_line_inflated_ratio": (bottom_line_inflated / sigma_base, ".2f"),
        "bottom_line_inflated_pct_of_signal":
            (bottom_line_inflated / fit_orbit.TRUE_OMEGA_DOT_DEG_YR * 100, ".2f"),
        "bottom_line_inflated_detection_nsigma":
            (fit_orbit.TRUE_OMEGA_DOT_DEG_YR / bottom_line_inflated, ".0f"),
        "spin_nsigma_headline_low": (spin_nsigma_headline["90"], ".1f"),
        "spin_nsigma_headline_high": (spin_nsigma_headline["100"], ".1f"),
        "offset_x_uas": (best[7] * 1000, ".1f"), "offset_y_uas": (best[8] * 1000, ".1f"),
        "m_ratio": (best[9], ".4f"), "d_ratio": (best[10], ".4f"),
    }
    for j, name in enumerate(fit_orbit.PARAM_NAMES):
        results[f"combined_sigma_{name}"] = (sig[j], ".5f")
        results[f"combined_ratio_{name}"] = (sig[j] / base[f"sigma_{name}"], ".2f")
    results_io.write_results("error_budget", results)


if __name__ == "__main__":
    main()
