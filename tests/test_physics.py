"""Physics regression values the manuscript relies on.

Each target below is a number the paper quotes (or that its argument
depends on), pinned here with a tolerance so a change to the orbital
engine or the photometric model that silently moved it would fail CI.
"""

import numpy as np

import campaign
import constants as k
import orbit
import photometry
import s301_lightcurve


def test_schwarzschild_precession_rate_for_s301():
    """1PN apsidal precession of the published Solution 1 orbit:
    0.2307 deg/yr (about 2.00 deg per 8.68-yr orbit)."""
    elements, omega_dot = campaign.truth_elements()
    rate_deg_yr = np.degrees(omega_dot) * k.year
    assert abs(rate_deg_yr - 0.2307) < 0.001
    # truth_elements() must agree with calling the formula directly.
    direct = orbit.schwarzschild_precession_rate(k.GM_BH, elements.sma, elements.ecc,
                                                 elements.period)
    assert np.isclose(direct, omega_dot, rtol=0, atol=1e-18)


def test_lense_thirring_aligned_spin_per_orbit():
    """Aligned maximal spin (chi=1, xi=0): 0.114 deg/orbit, the discovery
    paper's own quoted '0.11 deg * chi * cos(xi)'."""
    elements, _ = campaign.truth_elements()
    rate = orbit.lense_thirring_apsidal_rate(k.GM_BH, elements.sma, elements.ecc,
                                             elements.period, spin_chi=1.0, xi_deg=0.0)
    per_orbit_deg = np.degrees(rate) * elements.period
    assert abs(per_orbit_deg - 0.114) < 0.002
    # Misaligned by 90 deg the in-plane term vanishes (to floating-point
    # cos(pi/2) ~ 6e-17 of the aligned value); chi scales linearly.
    perpendicular = orbit.lense_thirring_apsidal_rate(k.GM_BH, elements.sma, elements.ecc,
                                                      elements.period, 1.0, 90.0)
    assert abs(perpendicular) < 1e-12 * rate
    half = orbit.lense_thirring_apsidal_rate(k.GM_BH, elements.sma, elements.ecc,
                                             elements.period, 0.5, 0.0)
    assert np.isclose(half, 0.5 * rate)


def test_light_curve_periapsis_amplitude():
    """Peak-to-trough K-band amplitude within +/-20 d of periapsis:
    0.115 mag from the band-integrated model."""
    elements, _ = s301_lightcurve.build_truth_elements()
    t_sec, _, one_plus_z, dmag = s301_lightcurve.compute_light_curve(
        elements, k.T_EFF_K, n_points=4000)
    mask = np.abs(t_sec / 86400) < 20
    amplitude = dmag[mask].max() - dmag[mask].min()
    assert abs(amplitude - 0.115) < 0.003
    # The redshift range the paper quotes for the same passage.
    assert 0.950 < one_plus_z.min() < 0.960
    assert 1.020 < one_plus_z.max() < 1.030


def test_band_integrated_model_is_not_the_bolometric_proxy():
    """photometry.py's amplitude over S301's (1+z) range is ~0.115 mag;
    the old -2.5 log10((1+z)^-3) proxy would give 0.227 mag. Guard
    against anyone reintroducing the proxy."""
    one_plus_z = np.array([0.9552, 1.0, 1.0242])
    dmag = photometry.dmag_from_one_plus_z(one_plus_z)
    span = dmag.max() - dmag.min()
    proxy_span = 7.5 * np.log10(one_plus_z.max() / one_plus_z.min())
    assert abs(span - 0.115) < 0.003
    assert abs(proxy_span - 0.227) < 0.002
    assert span < 0.6 * proxy_span
    # Blueshift (1+z < 1) brightens: more negative delta-mag.
    assert dmag[0] < dmag[1] < dmag[2]


def test_campaign_photometry_uses_shared_model():
    """campaign.true_observables' dmag channel must be the same
    band-integrated model, not a separate implementation."""
    epochs = campaign.build_epoch_grid_heuristic()
    state = campaign.true_state(epochs)
    one_plus_z = orbit.redshift_factor(state, k.Rs)
    _, _, _, dmag = campaign.true_observables(epochs)
    assert np.allclose(dmag, photometry.dmag_from_one_plus_z(one_plus_z))
