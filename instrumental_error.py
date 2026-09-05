"""
Numerical estimate of the leading INSTRUMENTAL systematic error term for
S301's measured apsidal precession -- the second piece of the
comprehensive error budget (random + astrophysical + instrumental) this
project was always meant to include, alongside fit_orbit.py's bootstrap
(random) and extended_mass_error.py's Newtonian-confusion astrophysical
term.

The term quantified here is the astrometric reference-frame tie: the
NIR (GRAVITY/adaptive-optics) frame in which S301's position is measured
is registered to Sgr A*'s true dynamical center via a separate radio
(VLBI) astrometric solution, which itself carries a residual zero-point
offset and a slow linear drift. This is not a hypothetical concern --
it is the SAME real systematic the actual S2 precession-detection paper
identifies as its limiting factor:

  GRAVITY Collaboration (2020, A&A 636, L5, arXiv:2004.07187): "These
  reference frame parameters [x0, y0, vx0, vy0] are now the limiting
  factor in the precision of the detection of the SP [Schwarzschild
  precession] of S2," adopting priors from Plewa et al. (2015):
  x0 = -0.2 +/- 0.2 mas, y0 = 0.1 +/- 0.2 mas,
  vx0 = 0.05 +/- 0.1 mas/yr, vy0 = 0.06 +/- 0.1 mas/yr.

  GRAVITY Collaboration (2022, A&A 657, L12, arXiv:2112.07478), on the
  same systematic: "the NACO-frame zero point and drift on the one
  hand, and the pro- or retrograde precession on the other hand, are
  degenerate" -- i.e., this is not merely additional noise, it can be
  directly reinterpreted by a fit as apsidal precession, biasing
  omega_dot rather than just widening its uncertainty. That is exactly
  the failure mode this script tests for.

Method: this project's own campaign.py/fit_orbit.py pipeline has no
frame-tie nuisance parameters at all (x0, y0, vx0, vy0 are not part of
its 7-parameter model) -- realistic of a first-pass analysis that
assumes a perfect frame tie. We inject a linear frame drift of the same
real magnitude GRAVITY quotes for this exact field into the existing
synthetic dataset (output/synthetic_observations.csv, unchanged
otherwise), refit with the unmodified pipeline, and measure the
resulting BIAS in the recovered omega_dot relative to the drift-free
fit of the identical dataset -- a systematic, not a random-noise,
comparison, following the same differencing methodology as
extended_mass_error.py.

Confusion/crowding, the other named instrumental term, is not
quantified numerically here: GRAVITY Collaboration (2022) states
interferometric astrometry (what this campaign uses) is "less affected
by a factor of several hundred" than single-aperture (NACO-style) AO
imaging by confusion noise, and no paper reviewed gives a citable
quantitative confusion floor for interferometric astrometry specifically
-- see the module docstring discussion in article.tex Section 4 for the
qualitative treatment.
"""

import numpy as np
from scipy.optimize import least_squares

import fit_orbit

# Plewa et al. (2015) reference-frame priors as adopted by GRAVITY
# Collaboration (2020, A&A 636, L5), Sect. on systematic errors --
# real, quoted central values and 1-sigma uncertainties for this exact
# NACO/GRAVITY-to-radio frame tie, not a guessed or illustrative number.
X0_CENTRAL_MAS, X0_SIGMA_MAS = -0.2, 0.2
Y0_CENTRAL_MAS, Y0_SIGMA_MAS = 0.1, 0.2
VX0_CENTRAL_MASYR, VX0_SIGMA_MASYR = 0.05, 0.1
VY0_CENTRAL_MASYR, VY0_SIGMA_MASYR = 0.06, 0.1


def inject_frame_drift(data, x0_mas, y0_mas, vx0_masyr, vy0_masyr, t_ref_yr):
    """Return a copy of the data with a linear reference-frame drift
    added to the observed (not true) RA/Dec: the real astrometric effect
    of an imperfect radio-to-infrared frame tie, which shifts where
    Sgr A* (and hence all offsets measured from it) appears to sit,
    growing linearly over the baseline from t_ref_yr."""
    drifted = {name: data[name].copy() for name in data.dtype.names}
    dt = data["epoch_yr"] - t_ref_yr
    drifted["ra_offset_mas"] = data["ra_offset_mas"] + x0_mas + vx0_masyr * dt
    drifted["dec_offset_mas"] = data["dec_offset_mas"] + y0_mas + vy0_masyr * dt
    return drifted


def refit(data, x0):
    """Thin wrapper around fit_orbit.fit_once for the naive (7-parameter,
    frame-tie-unaware) model."""
    return fit_orbit.fit_once(
        data["epoch_yr"], data["ra_offset_mas"], data["dec_offset_mas"],
        data["sigma_ra_mas"], data["sigma_dec_mas"], x0=x0,
    )


def report_case(label, baseline_fit, drifted_fit, bootstrap_sigma_omega_dot):
    """Print the naive fit's per-parameter bias from the injected frame
    drift, and omega_dot's bias in both sigma-equivalent and
    fraction-of-signal terms."""
    idx_omega_dot = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")
    bias = drifted_fit - baseline_fit
    print(f"\n=== {label} ===")
    for name, b in zip(fit_orbit.PARAM_NAMES, bias):
        print(f"  {name:<18} bias = {b:+.5f}")
    omega_dot_bias = bias[idx_omega_dot]
    print(f"  omega_dot bias vs. bootstrap random-error precision "
          f"({bootstrap_sigma_omega_dot:.4f} deg/yr): "
          f"{abs(omega_dot_bias) / bootstrap_sigma_omega_dot:.2f} sigma-equivalent")
    print(f"  omega_dot bias as fraction of the 1PN signal ({fit_orbit.TRUE_OMEGA_DOT_DEG_YR:.4f} deg/yr): "
          f"{abs(omega_dot_bias) / fit_orbit.TRUE_OMEGA_DOT_DEG_YR * 100:.3f}%")


FRAME_PARAM_NAMES = ["x0_mas", "y0_mas", "vx0_masyr", "vy0_masyr"]
# Same real Plewa et al. (2015) values used to inject the drift, now
# reused as the joint fit's Gaussian priors -- exactly GRAVITY
# Collaboration (2020)'s own stated approach ("adopted priors from
# Plewa et al. 2015") for handling this degeneracy, rather than ignoring
# the frame tie as fit_orbit.py's baseline 7-parameter model does.
FRAME_PRIOR_MEAN = np.array([X0_CENTRAL_MAS, Y0_CENTRAL_MAS, VX0_CENTRAL_MASYR, VY0_CENTRAL_MASYR])
FRAME_PRIOR_SIGMA = np.array([X0_SIGMA_MAS, Y0_SIGMA_MAS, VX0_SIGMA_MASYR, VY0_SIGMA_MASYR])
FRAME_BOUNDS_LO = np.array([-2.0, -2.0, -1.0, -1.0])
FRAME_BOUNDS_HI = np.array([2.0, 2.0, 1.0, 1.0])


def model_observables_with_frame(params11, epochs_yr, t_ref_yr):
    """Same orbital model as fit_orbit.model_observables, plus the linear
    frame-tie drift (x0, y0, vx0, vy0) added on top -- the joint model a
    real analysis (matching GRAVITY Collaboration 2020's own approach)
    would fit, instead of assuming a perfect frame tie."""
    orbit_params, frame_params = params11[:7], params11[7:]
    ra_orbit, dec_orbit, _ = fit_orbit.model_observables(orbit_params, epochs_yr)
    x0_mas, y0_mas, vx0_masyr, vy0_masyr = frame_params
    dt = epochs_yr - t_ref_yr
    return ra_orbit + x0_mas + vx0_masyr * dt, dec_orbit + y0_mas + vy0_masyr * dt


def residuals_with_frame(params11, epochs, ra_obs, dec_obs, sigma_ra, sigma_dec, t_ref_yr):
    """Sigma-normalized astrometric residuals plus Gaussian-prior
    residuals on the 4 frame-tie parameters against their real quoted
    Plewa et al. (2015) uncertainties -- the prior terms are what let the
    fit separate the frame tie from genuine precession at all; without
    them the 11-parameter model is even more degenerate than the
    7-parameter one, not less."""
    ra_model, dec_model = model_observables_with_frame(params11, epochs, t_ref_yr)
    frame_params = params11[7:]
    prior_resid = (frame_params - FRAME_PRIOR_MEAN) / FRAME_PRIOR_SIGMA
    return np.concatenate([
        (ra_model - ra_obs) / sigma_ra,
        (dec_model - dec_obs) / sigma_dec,
        prior_resid,
    ])


def fit_once_with_frame(epochs, ra_obs, dec_obs, sigma_ra, sigma_dec, x0_11, t_ref_yr):
    """Levenberg-Marquardt fit of the 11-parameter joint (orbit +
    frame-tie) model."""
    bounds_lo = np.concatenate([fit_orbit.BOUNDS_LO, FRAME_BOUNDS_LO])
    bounds_hi = np.concatenate([fit_orbit.BOUNDS_HI, FRAME_BOUNDS_HI])
    result = least_squares(
        residuals_with_frame, x0=x0_11, bounds=(bounds_lo, bounds_hi),
        args=(epochs, ra_obs, dec_obs, sigma_ra, sigma_dec, t_ref_yr),
    )
    return result.x


def demonstrate_joint_fit_fix(drifted_data, t_ref_yr, x0_orbit_guess, true_frame_params, label):
    """Fit the drift-injected data with the joint 11-parameter model
    (orbit + frame-tie nuisance parameters with Gaussian priors) instead
    of fit_orbit.py's naive 7-parameter model, and report whether this
    recovers omega_dot close to truth -- the real-world fix GRAVITY
    Collaboration (2020) uses for the same problem with S2."""
    x0_11 = np.concatenate([x0_orbit_guess, FRAME_PRIOR_MEAN])
    fitted_11 = fit_once_with_frame(
        drifted_data["epoch_yr"], drifted_data["ra_offset_mas"], drifted_data["dec_offset_mas"],
        drifted_data["sigma_ra_mas"], drifted_data["sigma_dec_mas"], x0_11, t_ref_yr,
    )
    orbit_fit, frame_fit = fitted_11[:7], fitted_11[7:]
    idx_omega_dot = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")
    omega_dot_true = fit_orbit.TRUE_OMEGA_DOT_DEG_YR
    omega_dot_fit = orbit_fit[idx_omega_dot]

    print(f"\n=== Joint fit with frame-tie nuisance parameters: {label} ===")
    print(f"  true injected frame params:  {dict(zip(FRAME_PARAM_NAMES, true_frame_params))}")
    print(f"  recovered frame params:      {dict(zip(FRAME_PARAM_NAMES, np.round(frame_fit, 4)))}")
    print(f"  omega_dot: truth={omega_dot_true:.4f}, recovered={omega_dot_fit:.4f}, "
          f"residual={omega_dot_fit - omega_dot_true:+.4f} deg/yr "
          f"({abs(omega_dot_fit - omega_dot_true) / omega_dot_true * 100:.3f}% of signal)")
    return fitted_11


def bootstrap_samples_with_frame(data, best_fit_11, t_ref_yr, n_resamples=200):
    """Same bootstrap convention as fit_orbit.bootstrap_samples, but for
    the 11-parameter joint (orbit + frame-tie) model -- gives the
    properly propagated omega_dot uncertainty once the frame-tie
    systematic is marginalized over via its real Gaussian prior, rather
    than either ignoring it (fit_orbit.py's baseline) or leaving it as an
    unmitigated bias (the naive-fit case above)."""
    rng = np.random.default_rng(fit_orbit.RNG_SEED)
    n = len(data["epoch_yr"])
    samples = np.zeros((n_resamples, len(best_fit_11)))
    for k_iter in range(n_resamples):
        idx = rng.integers(0, n, size=n)
        samples[k_iter] = fit_once_with_frame(
            data["epoch_yr"][idx], data["ra_offset_mas"][idx], data["dec_offset_mas"][idx],
            data["sigma_ra_mas"][idx], data["sigma_dec_mas"][idx],
            x0_11=best_fit_11, t_ref_yr=t_ref_yr,
        )
    return samples


def main():
    """Run the naive-vs-joint-fit comparison for both the central and
    conservative Plewa et al. frame-tie cases, reporting the naive bias,
    the joint fit's residual, and its properly propagated uncertainty."""
    data = fit_orbit.load_observations()
    t_ref_yr = data["epoch_yr"].min()  # drift defined as zero at campaign start

    x0_guess = fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION
    baseline_fit = refit(data, x0_guess)

    # Bootstrap random-error precision on omega_dot, for scale comparison
    # -- reuses fit_orbit.py's own established bootstrap convention
    # rather than re-deriving a different uncertainty estimate.
    samples = fit_orbit.bootstrap_samples(data, baseline_fit, n_resamples=200)
    bootstrap_sigma_omega_dot = samples[:, fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")].std()

    print("Baseline (drift-free) fit vs. bootstrap random-error precision established.")
    print(f"  bootstrap sigma(omega_dot) = {bootstrap_sigma_omega_dot:.4f} deg/yr "
          f"(n=200 resamples, for scale)")

    for label, (x0_mas, y0_mas, vx0_masyr, vy0_masyr) in (
        ("Central Plewa et al. (2015) frame-tie values",
         (X0_CENTRAL_MAS, Y0_CENTRAL_MAS, VX0_CENTRAL_MASYR, VY0_CENTRAL_MASYR)),
        ("Conservative (1-sigma-worse, all terms pushed further from zero)",
         (X0_CENTRAL_MAS - X0_SIGMA_MAS,   # -0.2 -> -0.4 mas
          Y0_CENTRAL_MAS + Y0_SIGMA_MAS,   #  0.1 ->  0.3 mas
          VX0_CENTRAL_MASYR + VX0_SIGMA_MASYR,  # 0.05 -> 0.15 mas/yr
          VY0_CENTRAL_MASYR + VY0_SIGMA_MASYR)),  # 0.06 -> 0.16 mas/yr
    ):
        drifted_data = inject_frame_drift(data, x0_mas, y0_mas, vx0_masyr, vy0_masyr, t_ref_yr)
        drifted_fit = refit(drifted_data, x0_guess)
        report_case(label, baseline_fit, drifted_fit, bootstrap_sigma_omega_dot)
        fitted_11 = demonstrate_joint_fit_fix(
            drifted_data, t_ref_yr, x0_guess,
            true_frame_params=(x0_mas, y0_mas, vx0_masyr, vy0_masyr), label=label,
        )
        joint_samples = bootstrap_samples_with_frame(drifted_data, fitted_11, t_ref_yr)
        joint_sigma_omega_dot = joint_samples[:, fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")].std()
        print(f"  properly-propagated omega_dot uncertainty (frame-tie marginalized, "
              f"n=200 resamples): {joint_sigma_omega_dot:.4f} deg/yr "
              f"(vs. {bootstrap_sigma_omega_dot:.4f} deg/yr if the frame tie were perfectly known)")


if __name__ == "__main__":
    main()
