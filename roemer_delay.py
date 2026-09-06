"""
Quantifies the Roemer (light-travel-time) delay: what ignoring it would
cost, and what including it costs.

Since the pipeline upgrade that made orbit.orbit_state_observed the
default in campaign.py and fit_orbit.py, the headline result (Table 3)
already includes the Roemer delay on both the truth and model sides, as
the S301 discovery paper's own orbit fit does. This script therefore
answers two questions the manuscript's error budget needs:

  1. BIAS if ignored: fit the (light-time-correct) noise-free synthetic
     data with a light-time-UNAWARE model (fit_orbit.fit_once with
     light_time=False). The per-parameter offsets are the pure systematic
     a naive analysis would incur, isolated from noise.
  2. PRECISION cost of the correct treatment: bootstrap the light-time-
     aware fit on noisy data and compare its omega_dot sigma with the
     naive model's bootstrap on the same data. Unlike the reference-frame
     zero point (reference_frame_error.py), the light-time correction
     adds no free parameters, so it should be free in precision as well
     as exact in bias -- checked here rather than asserted.

Both results are written to results/roemer_delay.json.
"""

import numpy as np

import campaign
import constants as k
import fit_orbit
import orbit
import results_io

N_BOOT = 300


def light_time_days(epochs_yr):
    """Emission-minus-arrival time (days) at each epoch for the truth orbit."""
    elements, omega_dot = campaign.truth_elements()
    t_obs = epochs_yr * k.year
    t_emit = orbit.solve_emission_time(t_obs, elements, omega_dot, elements.t_peri, k.GM_BH)
    return (t_emit - t_obs) / 86400.0


def main():
    """Report the light-time scale, the naive-model bias, and the
    bootstrap precision of naive vs. light-time-aware fits."""
    data = fit_orbit.load_observations()
    epochs = data["epoch_yr"]
    elements, _ = campaign.truth_elements()
    max_swing_days = 2 * elements.sma / k.c / 86400.0
    lt_days = light_time_days(epochs)
    print(f"Maximum possible light-time swing (2a/c): {max_swing_days:.2f} days")
    print(f"Actual correction across campaign epochs: [{lt_days.min():.3f}, {lt_days.max():.3f}] days "
          f"(RMS {np.sqrt(np.mean(lt_days ** 2)):.3f} days)")

    baseline = results_io.read_results("fit_orbit")
    sigma_ref = np.array([baseline[f"sigma_{name}"] for name in fit_orbit.PARAM_NAMES])
    idx_od = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")

    # 1. Noise-free bias of a light-time-unaware model on light-time-correct data.
    ra_true, dec_true, _, _ = campaign.true_observables(epochs)
    sig = np.full_like(epochs, k.GRAVITY_PLUS_ASTROMETRY_MAS)
    x0 = fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION
    naive_fit = fit_orbit.fit_once(epochs, ra_true, dec_true, sig, sig, x0=x0, light_time=False)
    aware_fit = fit_orbit.fit_once(epochs, ra_true, dec_true, sig, sig, x0=naive_fit, light_time=True)
    bias = naive_fit - fit_orbit.TRUTH_VECTOR
    print("\n=== Light-time-unaware model on noise-free light-time-correct data ===")
    for name, b_val, s_val in zip(fit_orbit.PARAM_NAMES, bias, sigma_ref):
        print(f"  {name:<18} bias={b_val:+.5f}  ({abs(b_val) / s_val:.1f} sigma of Table 3 precision)")
    aware_resid = aware_fit[idx_od] - fit_orbit.TRUTH_VECTOR[idx_od]
    print(f"  light-time-aware model on the same data: omega_dot residual {aware_resid:+.2e} deg/yr")

    # 2. Bootstrap precision, naive vs. aware, on the real noisy dataset.
    x0_data = fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION
    best_aware = fit_orbit.fit_once(epochs, data["ra_offset_mas"], data["dec_offset_mas"],
                                    data["sigma_ra_mas"], data["sigma_dec_mas"], x0=x0_data)
    best_naive = fit_orbit.fit_once(epochs, data["ra_offset_mas"], data["dec_offset_mas"],
                                    data["sigma_ra_mas"], data["sigma_dec_mas"], x0=x0_data,
                                    light_time=False)
    sig_aware = fit_orbit.bootstrap_samples(data, best_aware, n_resamples=N_BOOT).std(axis=0)
    sig_naive = fit_orbit.bootstrap_samples(data, best_naive, n_resamples=N_BOOT,
                                            light_time=False).std(axis=0)
    print(f"\nBootstrap sigma(omega_dot), n={N_BOOT}: aware {sig_aware[idx_od]:.5f}, "
          f"naive {sig_naive[idx_od]:.5f} deg/yr (ratio {sig_aware[idx_od] / sig_naive[idx_od]:.3f})")

    results = {
        "max_swing_days": (max_swing_days, ".1f"),
        "lt_min_days": (lt_days.min(), ".2f"), "lt_max_days": (lt_days.max(), ".2f"),
        "lt_rms_days": (np.sqrt(np.mean(lt_days ** 2)), ".2f"),
        "naive_omega_dot_bias": (bias[idx_od], ".4f"),
        "naive_omega_dot_bias_pct": (abs(bias[idx_od]) / fit_orbit.TRUE_OMEGA_DOT_DEG_YR * 100, ".2f"),
        "naive_omega_dot_bias_nsigma": (abs(bias[idx_od]) / sigma_ref[idx_od], ".1f"),
        "naive_i_bias": (bias[2], ".2f"), "naive_Omega_bias": (bias[3], ".2f"),
        "naive_omega_bias": (bias[4], ".2f"),
        "naive_angle_nsigma_min": (min(abs(bias[j]) / sigma_ref[j] for j in (2, 3, 4)), ".0f"),
        "naive_angle_nsigma_max": (max(abs(bias[j]) / sigma_ref[j] for j in (2, 3, 4)), ".0f"),
        "naive_t_peri_bias_days": (bias[5] * 365.25, ".2f"),
        "aware_omega_dot_residual": (aware_resid, ".1e"),
        "sigma_aware": (sig_aware[idx_od], ".4f"),
        "sigma_naive_model": (sig_naive[idx_od], ".4f"),
        "sigma_ratio": (sig_aware[idx_od] / sig_naive[idx_od], ".2f"),
        "n_boot": N_BOOT,
    }
    results_io.write_results("roemer_delay", results)


if __name__ == "__main__":
    main()
