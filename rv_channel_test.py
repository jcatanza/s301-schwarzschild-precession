"""
What would a radial-velocity channel buy this campaign? Two precisions.

campaign.py proposes no RV channel (see its module docstring: at S301's
m_K=19.3, ERIS/VLT's background-limited precision scales to ~1621 km/s,
SNR~0.04 against the signal, so no TAC would grant the time). This
script verifies the consequence rather than asserting it: it adds a
simulated RV channel to fit_orbit.py's astrometric fit and measures the
change in every parameter's bootstrap sigma, at two precisions:

  eris     constants.SIGMA_RV_S301_KMS_LINEAR_SCALING (~1621 km/s), the
           realistic figure for the instrument that does today's
           S-cluster RV monitoring;
  harmoni  constants.SIGMA_RV_S301_KMS_HARMONI_ELT (~21 km/s), the same
           background-limited regime scaled to the ELT's collecting area
           and HARMONI's higher K-band resolving power -- indicative of
           what a future instrument could add, not a proposal-ready
           sensitivity (see the caveats in constants.py).

RV truth comes from campaign.true_observables (km/s, positive receding),
noise-injected with np.random.default_rng(campaign.RNG_SEED). Residuals
are the astrometric ones plus (rv_model - rv_obs)/sigma_rv with rv_model
from fit_orbit.model_observables. Each fit is warm-started from the
astrometry-only best fit (fit_orbit's perturbed start can fall outside
the basin once an extra channel is added). Uncertainties: 400 epoch
resamples with replacement, all columns together, as in
fit_orbit.bootstrap_samples; reported as ratios to the astrometry-only
sigmas in results/fit_orbit.json (negative pct_change = improvement).

Results go to results/rv_channel_test.json.
"""

import numpy as np
from scipy.optimize import least_squares

import campaign
import constants as k
import fit_orbit
import results_io

N_BOOT = 400
RV_CASES = (
    ("eris", "ERIS/VLT (background-limited scaling from S2)", k.SIGMA_RV_S301_KMS_LINEAR_SCALING),
    ("harmoni", "HARMONI/ELT (same regime, R and D^2 scaled)", k.SIGMA_RV_S301_KMS_HARMONI_ELT),
)


def residuals_with_rv(params, data):
    """Sigma-normalized (model - data) residuals: the astrometric channels
    exactly as fit_orbit.residuals, plus the RV channel."""
    ra_model, dec_model, rv_model = fit_orbit.model_observables(params, data["epoch_yr"])
    return np.concatenate([
        (ra_model - data["ra_offset_mas"]) / data["sigma_ra_mas"],
        (dec_model - data["dec_offset_mas"]) / data["sigma_dec_mas"],
        (rv_model - data["rv_kms"]) / data["sigma_rv_kms"],
    ])


def fit_with_rv(data, x0):
    """One Levenberg-Marquardt fit of the 7 parameters to astrometry + RV."""
    result = least_squares(residuals_with_rv, x0=x0,
                           bounds=(fit_orbit.BOUNDS_LO, fit_orbit.BOUNDS_HI), args=(data,))
    return result.x


def build_rv_dataset(obs, rv_true, sigma_rv):
    """The campaign CSV columns plus a noise-injected RV column at
    precision sigma_rv (km/s), as a dict of equal-length arrays so a
    bootstrap draw can index every column together."""
    rng = np.random.default_rng(campaign.RNG_SEED)
    data = {name: obs[name] for name in obs.dtype.names}
    data["rv_kms"] = rv_true + rng.normal(0.0, sigma_rv, size=rv_true.shape)
    data["sigma_rv_kms"] = np.full_like(rv_true, sigma_rv)
    return data


def bootstrap_with_rv(data, best_fit, n_boot=N_BOOT, seed=fit_orbit.RNG_SEED):
    """(n_boot, 7) matrix of astrometry+RV fits over epoch resamples drawn
    with replacement (all columns together), each refit from best_fit."""
    n = len(data["epoch_yr"])
    draws = np.random.default_rng(seed).integers(0, n, size=(n_boot, n))
    return np.array([
        fit_with_rv({name: arr[idx] for name, arr in data.items()}, best_fit) for idx in draws
    ])


def main():
    """Fit astrometry + RV at both precisions, bootstrap, compare to the
    astrometry-only sigmas, and write the results."""
    obs = fit_orbit.load_observations()
    base = results_io.read_results("fit_orbit")
    idx_od = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")
    astro_fit = fit_orbit.fit_once(
        obs["epoch_yr"], obs["ra_offset_mas"], obs["dec_offset_mas"],
        obs["sigma_ra_mas"], obs["sigma_dec_mas"],
        x0=fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION,
    )
    rv_true = campaign.true_observables(obs["epoch_yr"])[2]
    rv_max = float(np.abs(rv_true).max())
    results = {
        "n_boot": N_BOOT,
        "eris_rv_sigma": (k.SIGMA_RV_S301_KMS_LINEAR_SCALING, ".0f"),
        "harmoni_rv_sigma": (k.SIGMA_RV_S301_KMS_HARMONI_ELT, ".1f"),
        "harmoni_snr": (rv_max / k.SIGMA_RV_S301_KMS_HARMONI_ELT, ".0f"),
    }
    print(f"RV signal over the campaign: |rv| up to {rv_max:.0f} km/s; "
          f"ERIS sigma {k.SIGMA_RV_S301_KMS_LINEAR_SCALING:.0f} km/s (peak SNR "
          f"{rv_max / k.SIGMA_RV_S301_KMS_LINEAR_SCALING:.1f}), HARMONI sigma "
          f"{k.SIGMA_RV_S301_KMS_HARMONI_ELT:.1f} km/s (peak SNR "
          f"{rv_max / k.SIGMA_RV_S301_KMS_HARMONI_ELT:.0f})")

    for tag, label, sigma_rv in RV_CASES:
        data = build_rv_dataset(obs, rv_true, sigma_rv)
        best = fit_with_rv(data, astro_fit)
        print(f"\n=== {label}: sigma_rv = {sigma_rv:.1f} km/s, {N_BOOT} bootstrap resamples ===")
        samples = bootstrap_with_rv(data, best)
        sig = samples.std(axis=0)
        print(f"{'parameter':<18}{'astrometry-only':>18}{'+ RV':>14}{'ratio':>8}")
        for j, name in enumerate(fit_orbit.PARAM_NAMES):
            ratio = sig[j] / base[f"sigma_{name}"]
            print(f"{name:<18}{base[f'sigma_{name}']:>18.5f}{sig[j]:>14.5f}{ratio:>8.2f}")
            results[f"{tag}_sigma_{name}"] = (sig[j], ".5f")
            results[f"{tag}_ratio_{name}"] = (ratio, ".2f")
        pct_change = (sig[idx_od] / base["sigma_omega_dot_deg_yr"] - 1) * 100
        print(f"omega_dot: {best[idx_od]:.4f} +/- {sig[idx_od]:.4f} deg/yr; sigma change "
              f"{pct_change:+.1f}% vs. astrometry-only (negative = improvement)")
        results[f"{tag}_pct_change_omega_dot"] = (pct_change, ".1f")

    results_io.write_results("rv_channel_test", results)


if __name__ == "__main__":
    main()
