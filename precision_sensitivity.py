"""
Sensitivity of the headline result to the assumed per-epoch astrometric
precision.

The campaign adopts 207 uas per epoch, the precision ACHIEVED on S301 in
the discovery data (their pixel-grid statistical floor, Methods,
"Astrometric errors"). The S301 discovery paper's own GRAVITY+ forecast
is 100 uas ("as expected for the final performance of GRAVITY+"); it is
the optimistic case here. Both are run on the identical epoch grid and
noise seed so the manuscript can state the result at either. The result
keys keep the names "forecast" (100 uas) and "achieved" (207 uas).
Results go to results/precision_sensitivity.json.
"""

import numpy as np

import campaign
import fit_orbit
import results_io

SIGMA_CASES_UAS = {"forecast": 100.0, "achieved": 207.0}
N_BOOT = 400


def main():
    """Fit and bootstrap at each per-epoch precision."""
    epochs = campaign.build_epoch_grid()
    ra_true, dec_true, _, _ = campaign.true_observables(epochs)
    idx_od = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")
    results = {"n_boot": N_BOOT}
    for tag, sigma_uas in SIGMA_CASES_UAS.items():
        sigma_mas = sigma_uas / 1000.0
        rng = np.random.default_rng(campaign.RNG_SEED)
        ra_obs = ra_true + rng.normal(0, sigma_mas, epochs.shape)
        dec_obs = dec_true + rng.normal(0, sigma_mas, epochs.shape)
        sig = np.full_like(epochs, sigma_mas)
        data = {"epoch_yr": epochs, "ra_offset_mas": ra_obs, "dec_offset_mas": dec_obs,
                "sigma_ra_mas": sig, "sigma_dec_mas": sig}
        best = fit_orbit.fit_once(epochs, ra_obs, dec_obs, sig, sig,
                                  x0=fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION)
        boot = fit_orbit.bootstrap_samples(data, best, n_resamples=N_BOOT).std(axis=0)
        print(f"[{tag}: {sigma_uas:.0f} uas] omega_dot = {best[idx_od]:.4f} +/- {boot[idx_od]:.4f} deg/yr "
              f"(detection {best[idx_od] / boot[idx_od]:.0f} sigma); "
              f"sigma(i, Omega, omega) = {boot[2]:.4f}, {boot[3]:.4f}, {boot[4]:.4f} deg")
        results[f"{tag}_uas"] = (sigma_uas, ".0f")
        results[f"{tag}_omega_dot"] = (best[idx_od], ".4f")
        results[f"{tag}_omega_dot_sigma"] = (boot[idx_od], ".4f")
        results[f"{tag}_detection_nsigma"] = (best[idx_od] / boot[idx_od], ".0f")
        for j, name in enumerate(fit_orbit.PARAM_NAMES):
            results[f"{tag}_sigma_{name}"] = (boot[j], ".4f")
    results_io.write_results("precision_sensitivity", results)


if __name__ == "__main__":
    main()
