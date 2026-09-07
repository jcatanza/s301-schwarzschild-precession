"""
Correlated (non-white) astrometric systematics and their effect on the
precession precision.

Every other precision in this project assumes white noise: each epoch an
independent Gaussian draw at GRAVITY_PLUS_ASTROMETRY_MAS, and each
systematic either a static bias or per-epoch scatter. Real interferometric
astrometry also carries errors that are COHERENT across epochs, which do
not average down as sqrt(N) and which an epoch bootstrap cannot see. Two
such terms are modelled here, each anchored to a published GRAVITY figure:

  1. Per-observing-run COMMON-MODE offset. GRAVITY Collaboration (2022,
     A&A 657, L12, arXiv:2112.07478) state that their calibration
     "adds a systematic uncertainty of 60 uas, divided by the square root
     of the number of available calibrations". We take one calibration
     per observing run: every epoch within a run (epochs separated by
     less than RUN_GAP_DAYS) shares one N(0, COMMON_MODE_UAS) offset per
     coordinate, drawn independently run to run. A dense periapsis
     cluster is therefore ONE run, the conservative reading.

  2. Slowly varying REFERENCE term. The near-infrared photocentre of
     Sgr A* is the accretion flow; GRAVITY Collaboration (2018, A&A 618,
     L10, arXiv:1810.12641) find the centre of the flare orbits
     "consistent with the centroid of the flare orbit, to within the
     +/-50 uas uncertainties" of the mass centroid from S2. Slowly moving
     confusing neighbours act the same way. We model the term as an
     Ornstein-Uhlenbeck (red-noise) process with stationary sigma
     RED_NOISE_UAS per coordinate and correlation time RED_NOISE_TAU_YR
     (one observing season), so that epochs within a season share most of
     it and epochs in different years do not.

Both are added to the white noise; the fit is told only the white sigma,
as a real analysis would be. Uncertainties come from N_REALIZATIONS
independent noise realizations of the whole campaign (the only estimator
that sees correlated noise), compared with the white-only realization
scatter, and the epoch bootstrap is also run on one correlated realization
to measure how much it UNDER-estimates. Amplitudes are scanned so the
reader can rescale. Results go to results/correlated_noise_test.json;
error_budget.py folds the fiducial excess into the bottom line.
"""

import numpy as np

import constants as k
import fit_orbit
import results_io

N_REALIZATIONS = 40
SEED_BASE = 5000
N_BOOT_NAIVE = 300

RUN_GAP_DAYS = 30.0
COMMON_MODE_UAS = 60.0                  # GRAVITY 2022 calibration systematic
COMMON_MODE_SCAN_UAS = (0.0, 30.0, 60.0, 100.0)
RED_NOISE_UAS = 50.0                    # GRAVITY 2018 flare-centre vs mass-centroid bound
RED_NOISE_TAU_YR = 0.5                  # one observing season
RED_NOISE_SCAN_UAS = (0.0, 50.0, 100.0)

OMEGA_DOT_INDEX = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")
SHORT_NAME = {"omega_dot_deg_yr": "omega_dot"}


def observing_runs(epochs_yr, gap_days=RUN_GAP_DAYS):
    """Integer run label per epoch: a new run starts wherever consecutive
    (sorted) epochs are more than gap_days apart."""
    order = np.argsort(epochs_yr)
    sorted_epochs = epochs_yr[order]
    new_run = np.concatenate([[True], np.diff(sorted_epochs) * 365.25 > gap_days])
    labels_sorted = np.cumsum(new_run) - 1
    labels = np.empty_like(labels_sorted)
    labels[order] = labels_sorted
    return labels


def ou_process(epochs_yr, sigma, tau_yr, rng):
    """Exact sample of a stationary Ornstein-Uhlenbeck process with
    standard deviation `sigma` and correlation time `tau_yr` at the given
    (unsorted) epochs."""
    order = np.argsort(epochs_yr)
    sorted_epochs = epochs_yr[order]
    values = np.empty_like(sorted_epochs)
    values[0] = rng.normal(0.0, sigma)
    for j in range(1, len(sorted_epochs)):
        decay = np.exp(-(sorted_epochs[j] - sorted_epochs[j - 1]) / tau_yr)
        values[j] = values[j - 1] * decay + rng.normal(0.0, sigma * np.sqrt(1.0 - decay ** 2))
    out = np.empty_like(values)
    out[order] = values
    return out


def simulate(epochs, ra_true, dec_true, rng, common_uas, red_uas):
    """White noise plus per-run common-mode offsets plus a red-noise
    reference term. The returned per-epoch sigma is the WHITE sigma only,
    which is what the fit is told."""
    sigma = k.GRAVITY_PLUS_ASTROMETRY_MAS
    ra_obs = ra_true + rng.normal(0.0, sigma, size=epochs.shape)
    dec_obs = dec_true + rng.normal(0.0, sigma, size=epochs.shape)
    if common_uas > 0:
        runs = observing_runs(epochs)
        n_runs = runs.max() + 1
        ra_obs = ra_obs + rng.normal(0.0, common_uas / 1000.0, size=n_runs)[runs]
        dec_obs = dec_obs + rng.normal(0.0, common_uas / 1000.0, size=n_runs)[runs]
    if red_uas > 0:
        ra_obs = ra_obs + ou_process(epochs, red_uas / 1000.0, RED_NOISE_TAU_YR, rng)
        dec_obs = dec_obs + ou_process(epochs, red_uas / 1000.0, RED_NOISE_TAU_YR, rng)
    return {
        "epoch_yr": epochs, "ra_offset_mas": ra_obs, "dec_offset_mas": dec_obs,
        "sigma_ra_mas": np.full_like(epochs, sigma), "sigma_dec_mas": np.full_like(epochs, sigma),
    }


def fit_data(data):
    """fit_orbit.fit_once from the standard perturbed-truth warm start."""
    return fit_orbit.fit_once(
        data["epoch_yr"], data["ra_offset_mas"], data["dec_offset_mas"],
        data["sigma_ra_mas"], data["sigma_dec_mas"],
        x0=fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION,
    )


def realization_fits(epochs, ra_true, dec_true, common_uas, red_uas, n=N_REALIZATIONS):
    """Best-fit parameter vectors for n independent noise realizations."""
    fits = np.empty((n, len(fit_orbit.PARAM_NAMES)))
    for i in range(n):
        rng = np.random.default_rng(SEED_BASE + i)
        fits[i] = fit_data(simulate(epochs, ra_true, dec_true, rng, common_uas, red_uas))
    return fits


def summarize(label, fits, white_sigma):
    """Print and return the realization scatter and mean bias of omega_dot."""
    sigma = fits.std(axis=0, ddof=1)
    bias = fits.mean(axis=0) - fit_orbit.TRUTH_VECTOR
    ratio = sigma[OMEGA_DOT_INDEX] / white_sigma
    print(f"[{label:<24}] sigma(omega_dot) = {sigma[OMEGA_DOT_INDEX]:.5f} deg/yr "
          f"({ratio:.2f}x white-noise realizations), mean bias {bias[OMEGA_DOT_INDEX]:+.5f}")
    return sigma, bias


def main():
    """White-only baseline, fiducial correlated case, amplitude scans, and
    the epoch bootstrap's under-estimate on one correlated realization."""
    epochs = fit_orbit.load_observations()["epoch_yr"]
    ra_true, dec_true, _ = fit_orbit.model_observables(fit_orbit.TRUTH_VECTOR, epochs)
    table3_sigma = results_io.read_results("fit_orbit")["sigma_omega_dot_deg_yr"]
    n_runs = observing_runs(epochs).max() + 1
    print(f"{len(epochs)} epochs in {n_runs} observing runs (gap > {RUN_GAP_DAYS:.0f} d); "
          f"white sigma {k.GRAVITY_PLUS_ASTROMETRY_MAS * 1000:.0f} uas; "
          f"{N_REALIZATIONS} realizations per case")

    white_fits = realization_fits(epochs, ra_true, dec_true, 0.0, 0.0)
    white_sigma_vec = white_fits.std(axis=0, ddof=1)
    white_sigma = white_sigma_vec[OMEGA_DOT_INDEX]
    print(f"[white only              ] sigma(omega_dot) = {white_sigma:.5f} deg/yr "
          f"({white_sigma / table3_sigma:.2f}x the Table 3 bootstrap)")

    results = {
        "n_realizations": N_REALIZATIONS, "n_runs": int(n_runs), "run_gap_days": RUN_GAP_DAYS,
        "common_mode_uas": COMMON_MODE_UAS, "red_noise_uas": RED_NOISE_UAS,
        "red_tau_yr": RED_NOISE_TAU_YR,
        "white_sigma_omega_dot": (white_sigma, ".5f"),
        "white_over_boot_omega_dot": (white_sigma / table3_sigma, ".2f"),
    }

    fid_fits = realization_fits(epochs, ra_true, dec_true, COMMON_MODE_UAS, RED_NOISE_UAS)
    fid_sigma, fid_bias = summarize("fiducial: 60 + 50 uas", fid_fits, white_sigma)
    excess = np.sqrt(max(fid_sigma[OMEGA_DOT_INDEX] ** 2 - white_sigma ** 2, 0.0))
    results.update({
        "fid_sigma_omega_dot": (fid_sigma[OMEGA_DOT_INDEX], ".5f"),
        "fid_ratio_omega_dot": (fid_sigma[OMEGA_DOT_INDEX] / white_sigma, ".2f"),
        "fid_ratio_vs_boot_omega_dot": (fid_sigma[OMEGA_DOT_INDEX] / table3_sigma, ".2f"),
        "fid_bias_omega_dot": (fid_bias[OMEGA_DOT_INDEX], ".5f"),
        "fid_bias_nsigma": (abs(fid_bias[OMEGA_DOT_INDEX]) / fid_sigma[OMEGA_DOT_INDEX], ".2f"),
        "fid_excess_sigma_omega_dot": (excess, ".5f"),
        "fid_pct_of_signal": (fid_sigma[OMEGA_DOT_INDEX] / fit_orbit.TRUE_OMEGA_DOT_DEG_YR * 100, ".2f"),
        "fid_detection_nsigma": (fit_orbit.TRUE_OMEGA_DOT_DEG_YR / fid_sigma[OMEGA_DOT_INDEX], ".0f"),
    })
    for name, sig, base in zip(fit_orbit.PARAM_NAMES, fid_sigma, white_sigma_vec):
        short = SHORT_NAME.get(name, name)
        results[f"fid_sigma_{short}"] = (sig, ".5f")
        results[f"fid_ratio_{short}"] = (sig / base, ".2f")

    for amp in COMMON_MODE_SCAN_UAS:
        fits = realization_fits(epochs, ra_true, dec_true, amp, 0.0)
        sig, _ = summarize(f"common-mode {amp:.0f} uas only", fits, white_sigma)
        results[f"common{amp:.0f}_ratio"] = (sig[OMEGA_DOT_INDEX] / white_sigma, ".2f")
        results[f"common{amp:.0f}_sigma"] = (sig[OMEGA_DOT_INDEX], ".5f")
    for amp in RED_NOISE_SCAN_UAS:
        fits = realization_fits(epochs, ra_true, dec_true, 0.0, amp)
        sig, _ = summarize(f"red noise {amp:.0f} uas only", fits, white_sigma)
        results[f"red{amp:.0f}_ratio"] = (sig[OMEGA_DOT_INDEX] / white_sigma, ".2f")
        results[f"red{amp:.0f}_sigma"] = (sig[OMEGA_DOT_INDEX], ".5f")

    # What the epoch bootstrap would have reported on ONE fiducial
    # correlated realization: the estimator a real analysis would use.
    one = simulate(epochs, ra_true, dec_true, np.random.default_rng(SEED_BASE),
                   COMMON_MODE_UAS, RED_NOISE_UAS)
    best = fit_data(one)
    boot_sigma = fit_orbit.bootstrap_samples(one, best, n_resamples=N_BOOT_NAIVE)[:, OMEGA_DOT_INDEX].std()
    print(f"[epoch bootstrap, 1 realization] sigma(omega_dot) = {boot_sigma:.5f} deg/yr = "
          f"{boot_sigma / fid_sigma[OMEGA_DOT_INDEX]:.2f}x the realization scatter")
    results["fid_boot_sigma_omega_dot"] = (boot_sigma, ".5f")
    results["fid_boot_over_mc"] = (boot_sigma / fid_sigma[OMEGA_DOT_INDEX], ".2f")
    results_io.write_results("correlated_noise_test", results)


if __name__ == "__main__":
    main()
