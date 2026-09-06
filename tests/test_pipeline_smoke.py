"""End-to-end smoke test of the campaign -> fit pipeline, kept fast.

Uses campaign.build_epoch_grid_heuristic() (instant) rather than
build_epoch_grid(), whose Fisher-optimal search is the slow part of
campaign.py; a single fit_orbit.fit_once, not the 1000-resample
bootstrap. The point is that the full chain (truth orbit with Roemer
delay -> noise -> 7-parameter Levenberg-Marquardt fit) recovers the
injected precession and period, not to reproduce Table 3's precision.
"""

import numpy as np

import campaign
import constants as k
import fit_orbit
import orbit

SEED = 301


def _noisy_astrometry(epochs):
    """Truth + fixed-seed Gaussian noise at the campaign's per-epoch
    GRAVITY+ precision."""
    ra_true, dec_true, _, _ = campaign.true_observables(epochs)
    sigma = k.GRAVITY_PLUS_ASTROMETRY_MAS
    rng = np.random.default_rng(SEED)
    ra_obs = ra_true + rng.normal(0, sigma, size=epochs.shape)
    dec_obs = dec_true + rng.normal(0, sigma, size=epochs.shape)
    return ra_obs, dec_obs, np.full_like(epochs, sigma)


def test_heuristic_grid_is_the_documented_campaign():
    """131 visibility-filtered epochs spanning 2028.5 to ~2041.7."""
    epochs = campaign.build_epoch_grid_heuristic()
    assert len(epochs) == 131
    assert np.all(np.diff(epochs) > 0)
    assert epochs.min() >= campaign.CAMPAIGN_START_YR
    assert epochs.max() <= campaign.CAMPAIGN_END_YR
    # Everything survives the Mar-Sep Paranal season filter.
    frac = epochs - np.floor(epochs)
    assert np.all((frac >= campaign.VISIBILITY_START_FRAC) & (frac < campaign.VISIBILITY_END_FRAC))


def test_single_fit_recovers_precession_and_period():
    """From the standard perturbed start, one LM fit lands within
    0.01 deg/yr of the injected 1PN rate and 0.01 yr of the period."""
    epochs = campaign.build_epoch_grid_heuristic()
    ra_obs, dec_obs, sigma = _noisy_astrometry(epochs)
    x0 = fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION
    best_fit = fit_orbit.fit_once(epochs, ra_obs, dec_obs, sigma, sigma, x0=x0)

    omega_dot_idx = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")
    period_idx = fit_orbit.PARAM_NAMES.index("P_yr")
    assert abs(best_fit[omega_dot_idx] - fit_orbit.TRUE_OMEGA_DOT_DEG_YR) < 0.01
    assert abs(best_fit[period_idx] - k.TRUTH["P_yr"]) < 0.01
    # The fit actually moved off the perturbed start and toward the truth
    # on every parameter.
    assert np.all(np.abs(best_fit - fit_orbit.TRUTH_VECTOR)
                  < np.abs(fit_orbit.INITIAL_PERTURBATION))
    # Reduced chi^2 of order unity for correctly-modelled Gaussian noise.
    resid = fit_orbit.residuals(best_fit, epochs, ra_obs, dec_obs, sigma, sigma)
    chi2_red = np.sum(resid ** 2) / (resid.size - len(best_fit))
    assert 0.6 < chi2_red < 1.6


def test_roemer_delay_bounded_by_orbit_light_crossing():
    """orbit_state_observed(light_time=True) evaluates the orbit at an
    emission time that differs from the arrival time by at most the
    light-crossing time of the orbit's diameter, 2a/c (~7.9 d for S301);
    and it must actually differ from the light_time=False evaluation."""
    elements, omega_dot = campaign.truth_elements()
    epochs_yr = np.linspace(campaign.CAMPAIGN_START_YR, campaign.CAMPAIGN_END_YR, 50)
    t_obs = epochs_yr * k.year

    t_emit = orbit.solve_emission_time(t_obs, elements, omega_dot, elements.t_peri, k.GM_BH)
    max_shift_days = np.max(np.abs(t_emit - t_obs)) / 86400
    two_a_over_c_days = 2 * elements.sma / k.c / 86400
    assert two_a_over_c_days < 7.95
    assert max_shift_days <= two_a_over_c_days
    assert max_shift_days > 1.0   # the effect is real, days not seconds

    # Self-consistency: t_emit + (distance offset)/c == t_obs.
    state_emit = orbit.orbit_state_precessing(t_emit, elements, omega_dot,
                                              t_ref=elements.t_peri, grav_param=k.GM_BH)
    assert np.allclose(t_emit - state_emit["sky_z"] / k.c, t_obs, rtol=0, atol=1.0)

    with_lt = orbit.orbit_state_observed(t_obs, elements, omega_dot, elements.t_peri,
                                         k.GM_BH, light_time=True)
    without_lt = orbit.orbit_state_observed(t_obs, elements, omega_dot, elements.t_peri,
                                            k.GM_BH, light_time=False)
    ra_lt, dec_lt = orbit.sky_offset_mas(with_lt, k.D_OBS)
    ra_no, dec_no = orbit.sky_offset_mas(without_lt, k.D_OBS)
    assert np.max(np.hypot(ra_lt - ra_no, dec_lt - dec_no)) > k.GRAVITY_PLUS_ASTROMETRY_MAS
    assert np.allclose(with_lt["sky_x"], state_emit["sky_x"])
