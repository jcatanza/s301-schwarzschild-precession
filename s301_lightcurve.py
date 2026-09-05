"""
Simulated near-infrared (Ks-band) light curve of S301 orbiting Sgr A*, as
observed from Earth (ESO Paranal Observatory, Chile).

S301 is a real star: GRAVITY/VLTI collaboration, announced August 2026,
published in Nature as "Discovery of a star sensitive to the spin of
Sgr A*" (arXiv:2607.12664). All orbital elements and stellar parameters
used below are that paper's Solution 1 (see constants.py for the full
values and citations) -- nothing about the orbit geometry here is an
illustrative choice; it's the real published solution.

Physics included:
  1. Full Keplerian orbit (Newton-Raphson eccentric-anomaly solve),
     elements = the real published (P, e, i, Omega, omega, t_peri).
     Semi-major axis a is DERIVED from (GM, P) via Kepler's third law,
     not an independent input -- it should land near the paper's
     687 AU / 83.0 mas as a consistency check.
  2. Exact special-relativistic Doppler factor from the star's true 3D
     velocity (not a small-v approximation) -- captures both classical
     Doppler and transverse (time-dilation) redshift.
  3. Gravitational (Schwarzschild) redshift factor sqrt(1 - Rs/r).
  4. Relativistic intensity boosting via the Lorentz-invariant I_nu/nu^3,
     applied to a blackbody spectrum, integrated over GRAVITY/VLTI's own
     real K-band operating range (1.98-2.40 um) -- the actual instrument
     bandpass S301's published m_K=19.3 was measured in, not a generic
     photometric-system filter.
  5. NOT included: strong-field ray tracing / gravitational lensing
     magnification or occultation by the BH shadow (needs a specific
     line-of-sight alignment not covered by the elements alone), Kerr
     (spin) corrections, or Roemer light-travel-time delay.

ASSUMPTIONS flagged explicitly (not published facts about S301):
  - T_eff = 7300 K: NOT published. Standard-table estimate for the
    paper's published spectral type (F1.5V). R_star = 1.5 R_sun is the
    midpoint of the paper's published 1.4-1.6 R_sun range.
  - Observer distance: Earth's position on its own orbit shifts by
    <<1 AU against the 8277 pc distance to Sgr A* -- utterly negligible,
    so d_obs = R0 to many significant figures. This only sets the
    reported *absolute* apparent magnitude; it has no effect on the
    relativistic delta-mag modulation shape/amplitude computed below
    (that's a flux *ratio*, so a fixed observer distance cancels out).
  - Black hole spin = 0 (Schwarzschild). Measuring the real spin is
    literally what the S301 discovery paper is trying to do with future
    observations -- a Kerr correction is a higher-order effect on top of
    what's modeled here, deliberately out of scope for this script.
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

# pylint: disable=wrong-import-position
# matplotlib's backend must be set (above) before pyplot is imported.
import numpy as np
import matplotlib.pyplot as plt

import constants as k
import orbit
# pylint: enable=wrong-import-position


def planck_bnu(freq, temp):
    """Blackbody spectral radiance B_nu(freq, temp)."""
    x = k.h * freq / (k.kB * temp)
    return (2 * k.h * freq ** 3 / k.c ** 2) / np.expm1(x)


def build_truth_elements():
    """The real published orbit (constants.TRUTH), with the semi-major
    axis derived from (GM_BH, period) rather than taken as an input."""
    period = k.TRUTH["P_yr"] * k.year
    sma = orbit.semi_major_axis_from_period(k.GM_BH, period)
    return orbit.OrbitalElements(
        t_peri=0.0, period=period, ecc=k.TRUTH["e"], sma=sma,
        i_deg=k.TRUTH["i_deg"], raan_deg=k.TRUTH["Omega_deg"], omega_deg=k.TRUTH["omega_deg"],
    ), sma


def compute_light_curve(elements, temp, r_star, d_obs, n_points=40000):
    """Delta-magnitude and redshift time series over one full orbit,
    densely sampled near periapsis (see eccentric_anomaly_grid)."""
    ecc_anom = orbit.eccentric_anomaly_grid(n_points)
    t = orbit.time_since_periapsis(ecc_anom, elements.ecc, elements.period)
    state = orbit.orbit_state(t, elements, k.GM_BH)
    one_plus_z = orbit.redshift_factor(state, k.Rs)

    lam0 = k.KS_LAMBDA0
    sigma_lam = k.KS_FWHM / 2.355
    lam_grid = np.linspace(lam0 - 3 * sigma_lam, lam0 + 3 * sigma_lam, 400)
    freq_grid = k.c / lam_grid
    transmission = np.exp(-0.5 * ((lam_grid - lam0) / sigma_lam) ** 2)

    flux_band = np.zeros_like(t)
    for idx, boost in enumerate(one_plus_z):
        freq_emit = freq_grid * boost
        intensity_emit = planck_bnu(freq_emit, temp)
        intensity_obs = intensity_emit / boost ** 3
        flux_nu_obs = np.pi * (r_star / d_obs) ** 2 * intensity_obs
        flux_band[idx] = np.trapezoid(flux_nu_obs * transmission, freq_grid)
    flux_band = np.abs(flux_band)

    dmag = -2.5 * np.log10(flux_band / np.median(flux_band))
    return t, state, one_plus_z, dmag


def make_plots(t, state, one_plus_z, dmag, outpath):
    """Three-panel figure: full-orbit light curve, periapsis zoom, and the
    relativistic-vs-classical spectral shift comparison."""
    t_days = t / 86400
    fig, axes = plt.subplots(3, 1, figsize=(9, 11), sharex=False)
    fig.suptitle("Simulated as observed by GRAVITY+/VLTI (K-band, 1.98-2.40 um bandpass)",
                 fontsize=10, color="dimgray")

    ax = axes[0]
    ax.plot(t / k.year, dmag, color="#1f77b4", lw=1.2)
    ax.invert_yaxis()
    ax.set_xlabel("Time from periapsis (years)")
    ax.set_ylabel("Delta magnitude (fainter down)")
    ax.set_title("S301 simulated Ks-band light curve -- full 8.68-year orbit")
    ax.grid(alpha=0.3)

    mask = np.abs(t_days) < 20
    ax = axes[1]
    ax.plot(t_days[mask], dmag[mask], color="#d62728", lw=1.4)
    ax.invert_yaxis()
    ax.set_xlabel("Time from periapsis (days)")
    ax.set_ylabel("Delta magnitude (fainter down)")
    ax.set_title("Zoom: periapsis passage (+/-20 days)")
    ax.grid(alpha=0.3)

    ax = axes[2]
    ax.plot(t_days[mask], (one_plus_z[mask] - 1) * 1e2, color="#2ca02c", lw=1.4,
            label="total (SR+GR) redshift z x100")
    ax.plot(t_days[mask], state["beta_r"][mask] * 1e2, color="#9467bd", lw=1.0, ls="--",
            label="classical v_r/c x100")
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xlabel("Time from periapsis (days)")
    ax.set_ylabel("(%)")
    ax.set_title("Spectral shift near periapsis: relativity beats classical Doppler")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    os.makedirs("output", exist_ok=True)
    plt.savefig(outpath, dpi=150)
    return mask


def main():
    """Run the full light-curve simulation and save the figure."""
    elements, sma = build_truth_elements()
    r_p_check = sma * (1 - elements.ecc)

    print(f"Derived semi-major axis a = {sma / k.AU:.1f} AU  "
          f"(published: {k.TRUTH['a_mas'] / 1000 * k.R0_PC:.1f} AU)")
    print(f"Periapsis check         r_p = {r_p_check / k.Rs:.1f} R_S  "
          f"(published: {k.TRUTH['r_p_Rs']} R_S)")
    print(f"Schwarzschild radius     Rs = {k.Rs / 1e3:.3e} km")

    v_p = np.sqrt(k.GM_BH * (2 / r_p_check - 1 / sma))
    print(f"Periapsis speed check   v_p = {v_p / k.c * 100:.2f}% c  "
          f"(published: >8% c, {k.MAX_SPEED_KMS:.0f} km/s)")

    temp = k.T_EFF_K
    r_star = k.R_STAR
    d_obs = k.D_OBS

    t, state, one_plus_z, dmag = compute_light_curve(elements, temp, r_star, d_obs)

    print(f"Max blueshift factor observed (1+z)_min = {one_plus_z.min():.4f}")
    print(f"Max redshift factor observed  (1+z)_max = {one_plus_z.max():.4f}")

    l_star = 4 * np.pi * r_star ** 2 * k.sigma_SB * temp ** 4
    l_sun = 3.828e26
    mbol_sun = 4.74
    mbol = mbol_sun - 2.5 * np.log10(l_star / l_sun)
    m_bol = mbol + 5 * np.log10(d_obs / (10 * k.pc))
    bc_k = -0.5  # rough bolometric correction for T~7300K in K-band, order-of-magnitude only
    m_k_mean = m_bol - bc_k
    print(f"Rough mean apparent K-band mag (order-of-magnitude estimate): "
          f"{m_k_mean:.2f}  (published, GRAVITY: m_K = {k.M_K_APPARENT} +/- {k.SIGMA_M_K})")

    outpath = "output/s301_lightcurve.png"
    mask = make_plots(t, state, one_plus_z, dmag, outpath)
    print(f"Saved figure to {outpath}")
    print(f"Peak-to-trough delta-mag amplitude near periapsis: "
          f"{dmag[mask].max() - dmag[mask].min():.4f} mag")


if __name__ == "__main__":
    main()
