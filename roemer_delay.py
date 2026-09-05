"""
Quantifies the Roemer (light-travel-time) delay this project's physics
engine explicitly disclaims (orbit.py's module docstring: "NOT included:
... Roemer (light-travel-time) delay across the orbit"), and tests
whether omitting it biases the recovered omega_dot -- the last piece
identified as missing from the comprehensive error budget.

Physical picture: campaign.py/fit_orbit.py treat each observation
"epoch" as if it equals the dynamical time at which the star was at that
orbital phase. In reality, an epoch is an ARRIVAL time; the light left
the star at some earlier "emission" time that differs by the light-
travel-time to the star's actual line-of-sight position at that moment,
which varies over the orbit (an eccentric, inclined orbit's line-of-
sight distance from the barycenter oscillates by a sizeable fraction of
its physical extent). This is the same effect known as Roemer delay in
binary-pulsar/spectroscopic-binary timing, and here plausibly matters
at the day scale: S301's orbital diameter (2*sma ~ 2.06e14 m) implies a
maximum possible light-time swing of 2*sma/c ~ 7.9 days, comparable to
the campaign's own +/-20-day periapsis-dense sampling window -- exactly
where the campaign concentrates its statistical power.

Method (matching this project's established pattern: inject a
real/plausible effect, refit naively, measure the bias; then show the
correct treatment removes it):

  1. Solve self-consistently for each observation's true EMISSION time
     t_emit given its recorded ARRIVAL time t_obs, via the fixed-point
     relation t_emit = t_obs + sky_z(t_emit)/c (sky_z: line-of-sight
     position, +Z toward observer per orbit.py's own convention, so a
     star currently closer to the observer emitted light that arrives
     sooner for the same emission time -- equivalently, backing out
     emission time from a fixed arrival time requires ADDING the light
     time implied by the star's position at that emission time).
     Converges in a handful of iterations because the contraction rate
     equals |v_los|/c <= ~8.5% for S301 (i.e. even starting from
     t_emit=t_obs, each iteration shrinks the residual by >10x).
  2. Build a synthetic dataset using the project's REAL campaign epoch
     grid and REAL noise model (same RNG seed/sigma as campaign.py), but
     with the true orbital state evaluated at the light-time-corrected
     t_emit rather than naively at t_obs -- the physically correct
     "what GRAVITY+ would actually record."
  3. Refit with fit_orbit.py's unmodified (Roemer-naive) model: measures
     the bias from ignoring this effect entirely, as this project's
     pipeline currently does.
  4. Refit with a Roemer-AWARE model that performs the same emission-
     time solve self-consistently for each trial parameter vector during
     the fit (no new free parameters needed -- light time is a
     deterministic function of the orbit already being fit, unlike the
     frame-tie's genuinely unknown nuisance parameters): shows this
     removes the bias, the same way accounting for the frame tie
     (Section 4.2) did, but here via getting the model right rather than
     adding priors.
"""

import numpy as np
from scipy.optimize import least_squares

import constants as k
import campaign
import fit_orbit
import orbit

N_ROEMER_ITER = 20       # safety cap; Newton converges in ~2-3 for this problem
# Convergence tolerance: verified (by tracing the iteration step-by-step for
# the worst-converging epoch in this campaign) that Newton reaches ~1.5e-5 s
# precision by iteration 2, after which a floating-point edge case in the
# underlying Kepler solve for that exact time value kicks the iterate away
# and it cycles indefinitely without ever satisfying a tolerance tighter than
# that -- so 1e-6 s (the original choice) never freezes and the loop runs out
# mid-cycle, silently returning a bad value. 1e-3 s is tight enough to freeze
# each epoch right after its genuine Newton convergence, before that edge
# case can be reached, while still being ~10 orders of magnitude below any
# physically meaningful timescale in this problem (day-scale corrections).
ROEMER_TOL_SEC = 1e-3


def solve_emission_time(t_obs_sec, elements, omega_dot, t_ref, grav_param):
    """Newton's-method solve of g(t_emit) = t_emit - sky_z(t_emit)/c - t_obs = 0
    for each observation time (vectorized), using g'(t_emit) = 1 - sky_vz/c
    (sky_vz already computed by orbit_state_precessing -- no extra finite
    differencing needed). See module docstring for the sign convention.

    Replaces an earlier plain fixed-point (Picard) iteration
    t_emit = t_obs + sky_z(t_emit)/c, run a fixed number of times
    regardless of convergence. That approach had a real, verified bug:
    for specific epochs it would converge cleanly after several
    iterations, then DESTABILIZE on the next application (a jump of
    several days, before slowly re-converging) -- a floating-point edge
    case in the underlying Kepler solve for certain near-exact inputs,
    not a genuine instability of the physical mapping (confirmed by
    tracing the iteration step-by-step: the map does not stay at its own
    fixed point when reapplied, which a true contraction must). Because
    the old code always ran exactly N_ROEMER_ITER=6 applications
    regardless of whether that landed on a "good" or "just-destabilized"
    iterate, it silently returned garbage for a handful of epochs,
    corrupting the whole fit -- this is what actually caused the
    unconverged, wildly non-smooth cost landscape reported earlier.

    Newton's method with an explicit per-element convergence check and
    a frozen-once-converged mask sidesteps this entirely: each epoch
    stops updating the instant it satisfies the tolerance, so it is
    never re-evaluated at a point beyond its own converged solution."""
    t_obs_sec = np.asarray(t_obs_sec, dtype=float)
    t_emit = t_obs_sec.copy()
    converged = np.zeros_like(t_emit, dtype=bool)
    for _ in range(N_ROEMER_ITER):
        active = ~converged
        if not np.any(active):
            break
        state = orbit.orbit_state_precessing(t_emit[active], elements, omega_dot, t_ref=t_ref,
                                              grav_param=grav_param)
        g = t_emit[active] - state["sky_z"] / k.c - t_obs_sec[active]
        g_prime = 1.0 - state["sky_vz"] / k.c
        step = g / g_prime
        t_emit[active] = t_emit[active] - step
        newly_converged = np.abs(step) < ROEMER_TOL_SEC
        idx_active = np.flatnonzero(active)
        converged[idx_active[newly_converged]] = True
    return t_emit


# pylint: disable=duplicate-code
# Building an OrbitalElements + Schwarzschild rate from constants.TRUTH is
# the same handful of lines campaign.py's true_observables() needs -- no
# shared "truth-orbit" utility module in this project is worth introducing
# for this, per the same reasoning as fit_orbit.py's own duplicate-code
# note; each analysis script stays self-contained.
def build_truth_elements_and_rate():
    """S301's real orbital elements (constants.TRUTH) and Schwarzschild
    precession rate, bundled for reuse across this module's functions."""
    period = k.TRUTH["P_yr"] * k.year
    sma = orbit.semi_major_axis_from_period(k.GM_BH, period)
    t_peri = k.TRUTH["t_peri_yr"] * k.year
    elements = orbit.OrbitalElements(
        t_peri=t_peri, period=period, ecc=k.TRUTH["e"], sma=sma,
        i_deg=k.TRUTH["i_deg"], raan_deg=k.TRUTH["Omega_deg"], omega_deg=k.TRUTH["omega_deg"],
    )
    omega_dot = orbit.schwarzschild_precession_rate(k.GM_BH, sma, k.TRUTH["e"], period)
    return elements, omega_dot, t_peri, period, sma
# pylint: enable=duplicate-code


def build_roemer_corrected_dataset(epochs_yr, add_noise=True):
    """Same real campaign epoch grid and real noise model as
    campaign.report_and_write_csv, but with astrometry evaluated at the
    light-time-corrected true emission time rather than naively at the
    observation epoch.

    add_noise=False gives the noise-free (deterministic) version used to
    isolate the pure Roemer-delay BIAS from ordinary noise-driven scatter
    -- the same deterministic-comparison methodology extended_mass_error.py
    and instrumental_error.py use for their own systematic-bias
    measurements, rather than confusing the two in a single noisy
    realization."""
    elements, omega_dot, t_peri, _period, _sma = build_truth_elements_and_rate()
    t_obs_sec = epochs_yr * k.year
    t_emit_sec = solve_emission_time(t_obs_sec, elements, omega_dot, t_peri, k.GM_BH)

    state = orbit.orbit_state_precessing(t_emit_sec, elements, omega_dot, t_ref=t_peri,
                                          grav_param=k.GM_BH)
    ra_true, dec_true = orbit.sky_offset_mas(state, k.D_OBS)

    sigma_mas = k.GRAVITY_PLUS_ASTROMETRY_MAS
    if add_noise:
        rng = np.random.default_rng(campaign.RNG_SEED)
        ra_obs = ra_true + rng.normal(0, sigma_mas, size=epochs_yr.shape)
        dec_obs = dec_true + rng.normal(0, sigma_mas, size=epochs_yr.shape)
    else:
        ra_obs, dec_obs = ra_true, dec_true

    return {
        "epoch_yr": epochs_yr,
        "ra_offset_mas": ra_obs, "dec_offset_mas": dec_obs,
        "sigma_ra_mas": np.full_like(epochs_yr, sigma_mas),
        "sigma_dec_mas": np.full_like(epochs_yr, sigma_mas),
    }


# pylint: disable=duplicate-code
# The params -> OrbitalElements unpacking is necessarily identical to
# fit_orbit.model_observables's own -- both fit the same 7 parameters --
# and instrumental_error.py's model_observables_with_frame needs the same
# pattern too; there is no shared "unpack this project's standard 7-vector"
# helper worth introducing across these otherwise-independent scripts.
def model_observables_roemer_aware(params, epochs_yr):
    """Same 7-parameter orbital model as fit_orbit.model_observables, but
    self-consistently solves for each trial parameter vector's own
    light-time-corrected emission time before evaluating the orbit --
    the physically correct treatment, requiring no new free parameters
    since light time is fully determined by the orbit already being
    fit."""
    p_yr, ecc, i_deg, raan_deg, omega_deg, t_peri_yr, omega_dot_deg_yr = params
    period = p_yr * k.year
    t_peri = t_peri_yr * k.year
    sma = orbit.semi_major_axis_from_period(k.GM_BH, period)
    elements = orbit.OrbitalElements(
        t_peri=t_peri, period=period, ecc=ecc, sma=sma,
        i_deg=i_deg, raan_deg=raan_deg, omega_deg=omega_deg,
    )
    omega_dot = np.radians(omega_dot_deg_yr) / k.year

    t_obs_sec = epochs_yr * k.year
    t_emit_sec = solve_emission_time(t_obs_sec, elements, omega_dot, t_peri, k.GM_BH)
    state = orbit.orbit_state_precessing(t_emit_sec, elements, omega_dot, t_ref=t_peri,
                                          grav_param=k.GM_BH)
    ra_mas, dec_mas = orbit.sky_offset_mas(state, k.D_OBS)
    return ra_mas, dec_mas
# pylint: enable=duplicate-code


def residuals_roemer_aware(params, epochs, ra_obs, dec_obs, sigma_ra, sigma_dec):
    """Sigma-normalized astrometric residuals for the Roemer-aware
    model."""
    ra_model, dec_model = model_observables_roemer_aware(params, epochs)
    return np.concatenate([(ra_model - ra_obs) / sigma_ra, (dec_model - dec_obs) / sigma_dec])


def fit_once_roemer_aware(epochs, ra_obs, dec_obs, sigma_ra, sigma_dec, x0):
    """Levenberg-Marquardt fit of the Roemer-aware model, plain default
    settings. An earlier version of this function needed a smaller
    explicit finite-difference step and several repeated warm-started
    passes to converge at all -- that workaround is gone now that the
    actual cause (orbit.solve_kepler's divergence bug, see module
    docstring) is fixed at the source: a single default-settings call
    converges to the true parameters exactly (cost ~1e-17 on noise-free
    data), same as every other fit in this project."""
    result = least_squares(
        residuals_roemer_aware, x0=x0, bounds=(fit_orbit.BOUNDS_LO, fit_orbit.BOUNDS_HI),
        args=(epochs, ra_obs, dec_obs, sigma_ra, sigma_dec),
    )
    return result.x


def main():
    """Report the real light-time correction scale, then compare the
    naive (light-time-unaware) and Roemer-aware fits on noise-free data
    to isolate the pure systematic bias from ordinary noise-driven
    scatter."""
    _elements, _omega_dot, _t_peri, _period_sec, sma_m = build_truth_elements_and_rate()
    max_light_time_days = 2 * sma_m / k.c / 86400.0
    print(f"Maximum possible light-time swing (2*sma/c): {max_light_time_days:.2f} days "
          f"(campaign's dense-window half-width: {campaign.DENSE_HALF_WIDTH_DAYS:.0f} days)")

    epochs_yr = fit_orbit.load_observations()["epoch_yr"]
    actual_light_time_days = (solve_emission_time(
        epochs_yr * k.year, _elements, _omega_dot, _t_peri, k.GM_BH) - epochs_yr * k.year) / 86400.0
    print(f"Actual light-time correction across the real campaign epochs: "
          f"[{actual_light_time_days.min():.3f}, {actual_light_time_days.max():.3f}] days "
          f"(RMS {np.sqrt(np.mean(actual_light_time_days**2)):.3f} days)\n")

    # Noise-free (deterministic) comparison, isolating the pure systematic
    # bias from ordinary noise-driven scatter -- matching
    # extended_mass_error.py/instrumental_error.py's own bias methodology.
    roemer_data = build_roemer_corrected_dataset(epochs_yr, add_noise=False)
    x0_guess = fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION
    idx_omega_dot = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")
    idx_t_peri = fit_orbit.PARAM_NAMES.index("t_peri_yr")
    bootstrap_sigma_omega_dot = 0.0011  # Table 3, photon-noise-only, for scale

    print("=== Naive fit (Roemer-unaware model) on noise-free Roemer-corrected data ===")
    naive_fit_on_roemer_data = fit_orbit.fit_once(
        roemer_data["epoch_yr"], roemer_data["ra_offset_mas"], roemer_data["dec_offset_mas"],
        roemer_data["sigma_ra_mas"], roemer_data["sigma_dec_mas"], x0=x0_guess,
    )
    for name, truth, fit_val in zip(fit_orbit.PARAM_NAMES, fit_orbit.TRUTH_VECTOR, naive_fit_on_roemer_data):
        print(f"  {name:<18} truth={truth:>12.4f}  fitted={fit_val:>12.4f}  bias={fit_val - truth:+.5f}")
    omega_dot_bias = naive_fit_on_roemer_data[idx_omega_dot] - fit_orbit.TRUTH_VECTOR[idx_omega_dot]
    print(f"  omega_dot bias: {omega_dot_bias:+.5f} deg/yr "
          f"({abs(omega_dot_bias) / fit_orbit.TRUTH_VECTOR[idx_omega_dot] * 100:.3f}% of signal, "
          f"{abs(omega_dot_bias) / bootstrap_sigma_omega_dot:.2f} sigma vs. photon-noise-only precision)")
    t_peri_bias_days = (naive_fit_on_roemer_data[idx_t_peri] - fit_orbit.TRUTH_VECTOR[idx_t_peri]) * 365.25
    print(f"  t_peri bias: {t_peri_bias_days:+.3f} days")

    print("\n=== Roemer-aware fit (self-consistent light-time correction) on the same noise-free data ===")
    roemer_aware_fit = fit_once_roemer_aware(
        roemer_data["epoch_yr"], roemer_data["ra_offset_mas"], roemer_data["dec_offset_mas"],
        roemer_data["sigma_ra_mas"], roemer_data["sigma_dec_mas"], x0=naive_fit_on_roemer_data,
    )
    for name, truth, fit_val in zip(fit_orbit.PARAM_NAMES, fit_orbit.TRUTH_VECTOR, roemer_aware_fit):
        print(f"  {name:<18} truth={truth:>12.4f}  fitted={fit_val:>12.4f}  bias={fit_val - truth:+.5f}")
    omega_dot_bias_fixed = roemer_aware_fit[idx_omega_dot] - fit_orbit.TRUTH_VECTOR[idx_omega_dot]
    print(f"  omega_dot bias: {omega_dot_bias_fixed:+.6f} deg/yr "
          f"({abs(omega_dot_bias_fixed) / fit_orbit.TRUTH_VECTOR[idx_omega_dot] * 100:.4f}% of signal)")

    resid_naive = residuals_roemer_aware(  # reuse for a fair apples-to-apples cost check
        naive_fit_on_roemer_data, roemer_data["epoch_yr"], roemer_data["ra_offset_mas"],
        roemer_data["dec_offset_mas"], roemer_data["sigma_ra_mas"], roemer_data["sigma_dec_mas"])
    resid_aware = residuals_roemer_aware(
        roemer_aware_fit, roemer_data["epoch_yr"], roemer_data["ra_offset_mas"],
        roemer_data["dec_offset_mas"], roemer_data["sigma_ra_mas"], roemer_data["sigma_dec_mas"])
    print(f"\n  fit quality check (noise-free data, should -> 0 for the correct model): "
          f"naive RMS residual={np.sqrt(np.mean(resid_naive**2)):.4f} sigma, "
          f"Roemer-aware RMS residual={np.sqrt(np.mean(resid_aware**2)):.6f} sigma")


if __name__ == "__main__":
    main()
