"""
The two photometry tests the manuscript quotes, now backed by code.

campaign.py simulates a K-band delta-magnitude channel (dmag_K, 0.30 mag
per epoch) from photometry.py's band-integrated relativistic model, but
fit_orbit.py deliberately fits astrometry alone. Two questions follow,
and this script answers both from the same synthetic CSV:

  (a) JOINT FIT -- does adding the photometric channel to the fit
      tighten any orbital element? The photometric model at a trial
      parameter vector is the same one that generated the data:
      orbit.redshift_factor on fit_orbit.model_state, through
      photometry.dmag_from_one_plus_z, median-normalized over the
      campaign epochs. Residuals are the astrometric ones (exactly as
      fit_orbit.residuals) plus (dmag_model - dmag_K)/sigma_dmag. The
      fit is warm-started from the astrometry-only best fit: once extra
      channels are added, fit_orbit's standard perturbed start can land
      outside the basin. Uncertainties come from 500 epoch-resample
      bootstraps (all columns resampled together, as in
      fit_orbit.bootstrap_samples) and are reported as ratios to the
      astrometry-only sigmas in results/fit_orbit.json.

  (b) INDEPENDENT VALIDATION -- taking the astrometry-only best fit as
      given, predict dmag at the campaign epochs and compare to the
      noisy dmag_K with a chi-square: the photometry then acts as a
      check of the light-curve physics that was never used to obtain
      the orbit. With no photometric free parameters the expectation is
      chi2 ~ n +/- sqrt(2 n).

  (c) Context numbers: the truth model's peak-to-trough amplitude over
      the campaign epochs and the per-epoch photometric sigma.

A note on the bootstrap zero point: the data's dmag_K were normalized to
the median over the full 131-epoch campaign. joint_residuals therefore
evaluates the model on the full epoch set (one orbit-state evaluation,
same cost as the astrometry-only fit) and then selects the resampled
residuals by index, so a bootstrap draw with duplicated epochs keeps the
zero point the data were generated with rather than re-centring the
model on the resample's own median.

Every number the manuscript quotes from here goes to
results/photometry_checks.json (results_io.py / make_numbers.py).
"""

import numpy as np
from scipy.optimize import least_squares

import campaign
import constants as k
import fit_orbit
import orbit
import photometry
import results_io

N_BOOT = 500


def joint_model(params, epochs):
    """RA (mas), Dec (mas), and K-band delta-mag at `epochs` for a
    7-parameter vector, from ONE orbit-state evaluation. The delta-mag is
    median-normalized over `epochs`, exactly as campaign.true_observables
    generates the data."""
    state = fit_orbit.model_state(params, epochs)
    ra_mas, dec_mas = orbit.sky_offset_mas(state, k.D_OBS)
    dmag = photometry.dmag_from_one_plus_z(orbit.redshift_factor(state, k.Rs))
    return ra_mas, dec_mas, dmag


def joint_residuals(params, data, idx):
    """Sigma-normalized (model - data) residuals, astrometric channels
    then photometric, for the epochs selected by `idx` (the identity for
    the headline fit; a with-replacement draw for a bootstrap resample).
    See the module docstring for why the model is evaluated on the full
    epoch set and indexed afterwards."""
    ra_model, dec_model, dmag_model = joint_model(params, data["epoch_yr"])
    return np.concatenate([
        (ra_model[idx] - data["ra_offset_mas"][idx]) / data["sigma_ra_mas"][idx],
        (dec_model[idx] - data["dec_offset_mas"][idx]) / data["sigma_dec_mas"][idx],
        (dmag_model[idx] - data["dmag_K"][idx]) / data["sigma_dmag"][idx],
    ])


def fit_joint(data, idx, x0):
    """One Levenberg-Marquardt fit of the 7 parameters to the astrometric
    + photometric channels of the epochs in `idx`, from x0."""
    result = least_squares(joint_residuals, x0=x0,
                           bounds=(fit_orbit.BOUNDS_LO, fit_orbit.BOUNDS_HI), args=(data, idx))
    return result.x


def astrometry_only_fit(data):
    """fit_orbit.py's headline astrometry-only fit (same perturbed start),
    used both as the warm start for the joint fit and as the fixed orbit
    for the chi-square validation."""
    x0 = fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION
    return fit_orbit.fit_once(
        data["epoch_yr"], data["ra_offset_mas"], data["dec_offset_mas"],
        data["sigma_ra_mas"], data["sigma_dec_mas"], x0=x0,
    )


def bootstrap_joint(data, best_fit, n_boot=N_BOOT, seed=fit_orbit.RNG_SEED):
    """(n_boot, 7) matrix of joint-fit parameters over epoch resamples
    drawn with replacement, each refit from best_fit."""
    rng = np.random.default_rng(seed)
    n = len(data["epoch_yr"])
    samples = np.zeros((n_boot, len(best_fit)))
    for k_iter in range(n_boot):
        idx = rng.integers(0, n, size=n)
        samples[k_iter] = fit_joint(data, idx, best_fit)
    return samples


def chi_square_validation(data, astro_fit):
    """Chi-square of the astrometry-only orbit's predicted dmag against the
    noisy dmag_K, with no photometric free parameters (n degrees of
    freedom)."""
    _, _, dmag_pred = joint_model(astro_fit, data["epoch_yr"])
    resid = (dmag_pred - data["dmag_K"]) / data["sigma_dmag"]
    n = len(resid)
    chi2 = float(np.sum(resid ** 2))
    print("\n=== Independent validation: astrometry-only orbit -> predicted dmag vs. dmag_K ===")
    print(f"  chi2 = {chi2:.1f} for n = {n} epochs (chi2/n = {chi2 / n:.2f}); "
          f"expected {n} +/- {np.sqrt(2 * n):.1f}, i.e. {(chi2 - n) / np.sqrt(2 * n):+.2f} sigma from "
          f"expectation; RMS residual {np.sqrt(np.mean(resid ** 2)):.2f} sigma_dmag")
    return {
        "chi2": (chi2, ".1f"),
        "n_epochs": n,
        "chi2_reduced": (chi2 / n, ".2f"),
        "rms_resid_sigma": (np.sqrt(np.mean(resid ** 2)), ".2f"),
        "expected_chi2_sigma": (np.sqrt(2 * n), ".1f"),
        "chi2_nsigma_from_expected": ((chi2 - n) / np.sqrt(2 * n), ".2f"),
    }


def main():
    """Joint astrometry+photometry fit with bootstrap, then the
    chi-square validation, then the context numbers; write results."""
    data = fit_orbit.load_observations()
    base = results_io.read_results("fit_orbit")
    all_idx = np.arange(len(data["epoch_yr"]))

    astro_fit = astrometry_only_fit(data)
    joint_fit = fit_joint(data, all_idx, astro_fit)
    print(f"Joint fit (astrometry + photometry), {N_BOOT} bootstrap resamples...")
    samples = bootstrap_joint(data, joint_fit)
    joint_sigma = samples.std(axis=0)

    print(f"\n{'parameter':<18}{'astrometry-only':>18}{'joint fit':>18}{'sigma (astro)':>16}"
          f"{'sigma (joint)':>16}{'ratio':>8}")
    results = {"joint_n_boot": N_BOOT}
    for j, name in enumerate(fit_orbit.PARAM_NAMES):
        ratio = joint_sigma[j] / base[f"sigma_{name}"]
        print(f"{name:<18}{astro_fit[j]:>18.5f}{joint_fit[j]:>18.5f}{base[f'sigma_{name}']:>16.5f}"
              f"{joint_sigma[j]:>16.5f}{ratio:>8.2f}")
        results[f"joint_sigma_{name}"] = (joint_sigma[j], ".5f")
        results[f"joint_ratio_{name}"] = (ratio, ".2f")

    results.update(chi_square_validation(data, astro_fit))

    dmag_true = campaign.true_observables(data["epoch_yr"])[3]
    amplitude = float(dmag_true.max() - dmag_true.min())
    sigma_dmag = float(np.mean(data["sigma_dmag"]))
    print(f"\nTruth-model dmag amplitude over the campaign epochs: {amplitude:.3f} mag "
          f"(vs. {sigma_dmag:.2f} mag per-epoch noise)")
    results["model_amplitude_at_epochs"] = (amplitude, ".3f")
    results["sigma_dmag"] = (sigma_dmag, ".2f")
    results_io.write_results("photometry_checks", results)


if __name__ == "__main__":
    main()
