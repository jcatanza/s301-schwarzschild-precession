"""
Astrometric reference systematic for a GRAVITY+-only campaign, replacing
an earlier analysis that injected the Plewa et al. (2015) NACO-to-radio
frame tie (zero point plus 13-year linear drift).

Why the NACO frame tie does not apply here: GRAVITY Collaboration (2020,
A&A 636, L5) define the x0, y0, vx0, vy0 reference-frame parameters for
the NACO adaptive-optics imaging frame ("relative to the AO spectroscopic
or imaging frames"), and they limited the S2 precession detection because
that analysis combined decades of NACO data. GRAVITY's own positions are
referenced differently: "GRAVITY positions are directly referring to
Sgr A*, since it is visible in each exposure and since it is one of the
point sources in the multi-source model" (GRAVITY Collaboration 2022,
A&A 657, L12); the S301 discovery paper likewise states "We do not allow
for any coordinate system offsets, as our data are interferometrically
referenced directly to the near-infrared counterpart of Sgr A*."

What CAN go wrong for interferometric referencing, and what this script
quantifies:

  1. A constant zero-point OFFSET between Sgr A*'s near-infrared
     photocentre (the accretion flow / flares) and the mass centroid the
     orbit is actually bound to. This is bounded physically by the size
     of the emitting region: a few Schwarzschild radii, with 1 R_S =
     10 uas for Sgr A* (discovery paper) and the EHT ring diameter of
     51.8 uas (EHT Collaboration 2022). We adopt a 50 uas prior (fiducial)
     and 100 uas (conservative). Unlike a NACO drift, this offset is
     common to every epoch, so it is absorbed exactly if the fit carries
     two extra parameters (x0, y0) -- but it biases a fit that assumes the
     focus sits at the origin, which the 7-parameter headline model does.
  2. Per-epoch reference JITTER: the counterpart's photocentre wanders
     (flares orbit at a few R_S) and its per-epoch position has its own
     noise. GRAVITY's quoted 65-100 uas S2-minus-Sgr A* precision already
     includes this, so it is largely inside the campaign's 207 uas floor;
     we scan an additional 0 / 50 / 100 uas in quadrature as a sensitivity.
     Systematics that are CORRELATED across epochs (per-run common-mode
     offsets, slow photocentre drift) are the subject of
     correlated_noise_test.py, not this script.

Outputs (results/reference_frame_error.json): naive-model bias from an
ignored offset, the 9-parameter joint fit's residual and bootstrap
precision (the honest cost of the correct treatment), and the jitter scan.
"""

import numpy as np
from scipy.optimize import least_squares

import constants as k
import fit_orbit
import results_io

OFFSET_PRIOR_UAS = 50.0          # fiducial: a few R_S (1 R_S = 10 uas); EHT ring diameter 51.8 uas
OFFSET_PRIOR_UAS_CONSERVATIVE = 100.0
JITTER_SCAN_UAS = (0.0, 50.0, 100.0)
N_BOOT = 200
OFFSET_BOUNDS_MAS = 2.0


def model_with_offset(params9, epochs):
    """7-parameter orbit plus a constant (x0, y0) zero-point offset (mas)."""
    ra_orbit, dec_orbit, _ = fit_orbit.model_observables(params9[:7], epochs)
    return ra_orbit + params9[7], dec_orbit + params9[8]


def residuals_with_offset(params9, epochs, ra_obs, dec_obs, sigma_ra, sigma_dec, prior_sigma_mas):
    """Astrometric residuals plus Gaussian-prior residuals on (x0, y0)."""
    ra_model, dec_model = model_with_offset(params9, epochs)
    return np.concatenate([
        (ra_model - ra_obs) / sigma_ra,
        (dec_model - dec_obs) / sigma_dec,
        params9[7:] / prior_sigma_mas,
    ])


def fit_with_offset(epochs, ra_obs, dec_obs, sigma_ra, sigma_dec, x0_9, prior_sigma_mas):
    """Levenberg-Marquardt fit of the 9-parameter joint model."""
    lo = np.concatenate([fit_orbit.BOUNDS_LO, [-OFFSET_BOUNDS_MAS, -OFFSET_BOUNDS_MAS]])
    hi = np.concatenate([fit_orbit.BOUNDS_HI, [OFFSET_BOUNDS_MAS, OFFSET_BOUNDS_MAS]])
    result = least_squares(residuals_with_offset, x0=x0_9, bounds=(lo, hi),
                           args=(epochs, ra_obs, dec_obs, sigma_ra, sigma_dec, prior_sigma_mas))
    return result.x


def bootstrap_with_offset(data, best_9, prior_sigma_mas, n_resamples=N_BOOT):
    """Epoch-resampling bootstrap of the joint fit (same convention as fit_orbit)."""
    rng = np.random.default_rng(fit_orbit.RNG_SEED)
    n = len(data["epoch_yr"])
    samples = np.zeros((n_resamples, 9))
    for k_iter in range(n_resamples):
        idx = rng.integers(0, n, size=n)
        samples[k_iter] = fit_with_offset(
            data["epoch_yr"][idx], data["ra_offset_mas"][idx], data["dec_offset_mas"][idx],
            data["sigma_ra_mas"][idx], data["sigma_dec_mas"][idx], best_9, prior_sigma_mas)
    return samples


def inject_offset(data, x0_mas, y0_mas):
    """Copy of the dataset with a constant reference offset added."""
    shifted = {name: data[name].copy() for name in data.dtype.names}
    shifted["ra_offset_mas"] = data["ra_offset_mas"] + x0_mas
    shifted["dec_offset_mas"] = data["dec_offset_mas"] + y0_mas
    return shifted


def _fit7(data, x0_guess):
    """Baseline 7-parameter fit of a dataset dict (or structured array)."""
    return fit_orbit.fit_once(data["epoch_yr"], data["ra_offset_mas"], data["dec_offset_mas"],
                              data["sigma_ra_mas"], data["sigma_dec_mas"], x0=x0_guess)


def offset_case(data, baseline_fit, prior_uas, sigma_base, tag):
    """Inject a one-prior-sigma constant offset; report the 7-parameter
    model's bias and the 9-parameter joint fit's recovery and precision."""
    idx_od = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")
    x0_guess = fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION
    prior_mas = prior_uas / 1000.0
    component_mas = prior_mas / np.sqrt(2.0)      # |offset| = one prior sigma
    shifted = inject_offset(data, component_mas, component_mas)
    bias = _fit7(shifted, x0_guess) - baseline_fit
    print(f"\n=== Constant reference offset |d| = {prior_uas:.0f} uas, ignored by the 7-parameter model ===")
    for name, b_val in zip(fit_orbit.PARAM_NAMES, bias):
        print(f"  {name:<18} bias = {b_val:+.5f}")
    nsig = abs(bias[idx_od]) / sigma_base
    pct = abs(bias[idx_od]) / fit_orbit.TRUE_OMEGA_DOT_DEG_YR * 100
    print(f"  omega_dot bias: {bias[idx_od]:+.5f} deg/yr = {nsig:.2f} sigma, {pct:.2f}% of signal")

    joint = fit_with_offset(shifted["epoch_yr"], shifted["ra_offset_mas"], shifted["dec_offset_mas"],
                            shifted["sigma_ra_mas"], shifted["sigma_dec_mas"],
                            np.concatenate([x0_guess, [0.0, 0.0]]), prior_mas)
    resid = joint[idx_od] - baseline_fit[idx_od]
    sig_joint = bootstrap_with_offset(shifted, joint, prior_mas)[:, idx_od].std()
    print(f"  joint 9-parameter fit (prior {prior_uas:.0f} uas): recovered offset "
          f"({joint[7] * 1000:.1f}, {joint[8] * 1000:.1f}) uas vs injected "
          f"({component_mas * 1000:.1f}, {component_mas * 1000:.1f}); omega_dot residual {resid:+.5f} "
          f"({abs(resid) / fit_orbit.TRUE_OMEGA_DOT_DEG_YR * 100:.3f}%); bootstrap sigma "
          f"{sig_joint:.5f} deg/yr ({sig_joint / sigma_base:.2f}x Table 3)")
    return {
        f"naive_bias_{tag}": (bias[idx_od], ".4f"),
        f"naive_bias_nsigma_{tag}": (nsig, ".1f"),
        f"naive_bias_pct_{tag}": (pct, ".2f"),
        f"joint_residual_pct_{tag}": (abs(resid) / fit_orbit.TRUE_OMEGA_DOT_DEG_YR * 100, ".3f"),
        f"joint_sigma_{tag}": (sig_joint, ".4f"),
        f"joint_sigma_ratio_{tag}": (sig_joint / sigma_base, ".2f"),
        f"joint_offset_x_{tag}": (joint[7] * 1000, ".1f"),
        f"joint_offset_y_{tag}": (joint[8] * 1000, ".1f"),
        f"injected_offset_component_{tag}": (component_mas * 1000, ".1f"),
    }


def jitter_case(data, jitter_uas, sigma_base):
    """Bootstrap precision with extra per-epoch reference jitter in quadrature."""
    idx_od = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")
    inflated = {name: data[name].copy() for name in data.dtype.names}
    sig_new = np.sqrt(data["sigma_ra_mas"] ** 2 + (jitter_uas / 1000.0) ** 2)
    inflated["sigma_ra_mas"] = sig_new
    inflated["sigma_dec_mas"] = sig_new
    best = _fit7(inflated, fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION)
    sig = fit_orbit.bootstrap_samples(inflated, best, n_resamples=N_BOOT).std(axis=0)[idx_od]
    print(f"  jitter {jitter_uas:>5.0f} uas -> sigma(omega_dot) = {sig:.5f} deg/yr ({sig / sigma_base:.2f}x)")
    tag = f"jitter{int(jitter_uas)}"
    return {f"{tag}_sigma": (sig, ".4f"), f"{tag}_ratio": (sig / sigma_base, ".2f")}


def main():
    """Naive bias, joint-fit recovery and precision, and jitter scan."""
    data = fit_orbit.load_observations()
    sigma_base = results_io.read_results("fit_orbit")["sigma_omega_dot_deg_yr"]
    baseline_fit = _fit7(data, fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION)
    results = {"offset_prior_uas": (OFFSET_PRIOR_UAS, ".0f"),
               "offset_prior_uas_conservative": (OFFSET_PRIOR_UAS_CONSERVATIVE, ".0f"),
               "n_boot": N_BOOT}
    for tag, prior_uas in (("fid", OFFSET_PRIOR_UAS), ("cons", OFFSET_PRIOR_UAS_CONSERVATIVE)):
        results.update(offset_case(data, baseline_fit, prior_uas, sigma_base, tag))
    print(f"\n=== Per-epoch reference jitter added in quadrature to the "
          f"{k.GRAVITY_PLUS_ASTROMETRY_MAS * 1000:.0f} uas floor ===")
    for jitter_uas in JITTER_SCAN_UAS:
        results.update(jitter_case(data, jitter_uas, sigma_base))
    results_io.write_results("reference_frame_error", results)


if __name__ == "__main__":
    main()
