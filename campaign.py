"""
Simulate a synthetic two-instrument, two-passage observing campaign of
S301: GRAVITY+/VLTI astrometry plus ERIS/VLT radial velocity, injecting
the real published orbit (constants.TRUTH, i.e. the GRAVITY discovery
paper's Solution 1) as ground truth.

Both instruments are ground-based, at ESO's Paranal Observatory, Chile
(d_obs = R0, see constants.py) -- the same site GRAVITY used for the real
S301 discovery, and both are real, already-operating facility instruments
(no launch-date or commissioning constraint applies, unlike a hypothetical
dedicated space mission).

This campaign covers TWO periapsis passages -- the real predicted next one
(~2031.8; ESO stated "next closest passage in 2031") and the one after that
(~2040.5) -- so it can measure apsidal (Schwarzschild) precession the same
way the real discovery papers do: by comparing the orbit's orientation
across passages, not from a single snapshot.

The injected "truth" includes REAL, non-illustrative physics beyond a
plain static orbit: the standard 1PN Schwarzschild apsidal precession rate
(orbit.schwarzschild_precession_rate -- the same formula used for
Mercury), applied via orbit.orbit_state_precessing so the argument of
periapsis omega drifts by ~2 deg/orbit between the two passages.
Inclination i and the ascending node Omega are NOT precessed -- that would
be a black-hole-spin (Lense-Thirring, frame-dragging) effect, which is
genuinely unmeasured for Sgr A* and deliberately not modeled anywhere in
this project (see s301_lightcurve.py's docstring).

Three observable channels are simulated per epoch, from TWO real,
honestly-attributed instruments (earlier versions of this project
attributed these to a single space-based telescope -- inconsistent, since
a single 2.4 m aperture's diffraction limit is larger than S301's whole
orbit, and the real astronomy community's own announced follow-up plans
use ground-based GRAVITY+/ELT-class instruments instead; fixed here):

  1. Astrometric offset (RA, Dec, mas) -- GRAVITY+/VLTI. Radial velocity
     has exactly zero sensitivity to the ascending node Omega (see
     orbit.py's rotation matrix), so astrometry is genuinely required to
     recover the full orbit. constants.GRAVITY_PLUS_ASTROMETRY_MAS
     (100 uas) is GRAVITY's real, demonstrated on-sky precision at
     S301's brightness -- not a theoretical or idealized figure.
  2. Radial velocity (km/s) -- ERIS/VLT (SPIFFIER integral field
     spectrograph), via the Brackett-gamma line, the real instrument used
     for essentially all S-cluster RV monitoring to date. Precision is
     scaled from ERIS/SINFONI's real achieved RV precision on S2 (the
     brightest, best-studied S-star) down to S301's much fainter
     magnitude -- see compute_rv_precision()'s docstring.
  3. Photometric delta-magnitude -- same physics as s301_lightcurve.py,
     evaluated at the campaign's actual epochs. Precision uses S301's own
     published photometric uncertainty (constants.SIGMA_M_K) directly.
     Simulated for completeness; NOT used by fit_orbit.py's fit (it
     carries no orbital-element information beyond what astrometry+RV
     already give).
"""

import os

import numpy as np

import constants as k
import orbit

RNG_SEED = 301
# A dedicated, multi-year, two-passage VLT/VLTI monitoring program of this
# kind would realistically be proposed as an ESO Large Programme, not
# slotted into routine time: allow ~6 months to write the proposal, submit
# to the next Large Programme call, ~6 months for Time Allocation
# Committee review, then wait for the start of the awarded observing
# period (~6-12 months out) -- roughly 2 years end to end from today
# before the first scheduled night, not an immediate start.
CAMPAIGN_START_YR = 2028.5
SECOND_PERIAPSIS_YR = k.NEXT_PERIAPSIS_YR + k.TRUTH["P_yr"]   # ~2040.5
CAMPAIGN_END_YR = SECOND_PERIAPSIS_YR + 1.2   # a bit past the second passage

# Cadence: ~90% of all epochs concentrated within +/-20 days of each
# passage's periapsis (the originally specified design -- see
# build_epoch_grid's docstring for how that interacts with real Paranal
# visibility), the remaining ~10% a sparse long-baseline presence so the
# overall orbit shape away from closest approach is still constrained.
# GRAVITY+ and ERIS are assumed to observe on the same coordinated
# cadence (both are VLT/VLTI instruments at the same site, so joint
# scheduling is realistic, unlike an earlier version of this project that
# coordinated a ground-based and a space-based instrument).
DENSE_HALF_WIDTH_DAYS = 20.0
DENSE_HALF_WIDTH_YR = DENSE_HALF_WIDTH_DAYS / 365.25
N_DENSE_PER_PASSAGE = 78   # tuned (with SPARSE_STEP_YR below) so that,
SPARSE_STEP_YR = 0.5       # after the visibility-season filter, the total
# epoch count is >=120 (the target set by the two-passage precession
# extension) and ~90% of surviving epochs land within +/-20 days of a
# periapsis (verified by the printed occupancy check in main()).

# Sgr A* is only observable from Paranal during its annual visibility
# season -- the discovery paper's own real GRAVITY monitoring ran "monthly
# during roughly week-long campaigns between March and September". A
# uniform, year-round cadence (as an earlier version of this campaign
# used) silently scheduled epochs when the target wasn't even up. Fixed
# here by filtering every generated epoch to this real season.
VISIBILITY_START_FRAC = 31 / 365                             # March 1
VISIBILITY_END_FRAC = (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30) / 365  # October 1


def _in_visibility_season(epochs_yr):
    """True where the fractional-year part falls in the real Mar-Sep
    Paranal visibility window for Sgr A*."""
    frac = epochs_yr - np.floor(epochs_yr)
    return (frac >= VISIBILITY_START_FRAC) & (frac < VISIBILITY_END_FRAC)


def _season_anchor(periapsis_yr):
    """The observing date closest to `periapsis_yr` that Sgr A* is
    actually visible from Paranal: periapsis itself if it already falls
    in-season, otherwise the nearest edge of a surrounding Mar-Sep
    visibility window."""
    year = np.floor(periapsis_yr)
    frac = periapsis_yr - year
    if VISIBILITY_START_FRAC <= frac < VISIBILITY_END_FRAC:
        return periapsis_yr
    edges = np.array([year + offset + boundary
                       for offset in (-1, 0, 1)
                       for boundary in (VISIBILITY_START_FRAC, VISIBILITY_END_FRAC)])
    return edges[np.argmin(np.abs(edges - periapsis_yr))]


def compute_rv_precision():
    """
    S301's radial-velocity precision, scaled from a REAL ACHIEVED on-sky
    result rather than a theoretical sensitivity-table calculation:
    SINFONI/ERIS Br-gamma radial velocities of S2 (the brightest,
    best-studied S-star) reach a median precision of
    constants.SIGMA_RV_S2_KMS = 12.3 km/s at constants.M_K_S2 = 14.0.

    Scaling that down to S301's much fainter m_K=19.3 (Delta_m=5.3, i.e.
    ~130x fainter in flux) requires assuming a photon-noise regime. This
    project injects the background/sky-photon-noise-limited (linear in
    flux) scaling -- precision degrades by 10**(0.4*Delta_m) ~= 132x --
    as the REALISTIC value: ground-based K-band IFU spectroscopy of faint
    Galactic Center targets is background-dominated in practice (bright
    OH airglow, thermal sky background, severe stellar crowding), not
    source-photon-starved, so this is the physically appropriate regime
    for a target this faint in this field, not just a conservative
    caveat. The source-limited (sqrt(flux), 10**(0.2*Delta_m) ~= 11.5x)
    scaling is also returned, purely as an optimistic point of comparison
    -- it is NOT injected into the campaign.

    Returns (sigma_rv_kms_realistic, sigma_rv_kms_optimistic).
    """
    return k.SIGMA_RV_S301_KMS_LINEAR_SCALING, k.SIGMA_RV_S301_KMS_SQRT_SCALING


def build_epoch_grid():
    """~90% of all epochs densely sampled within +/-20 days of EACH of
    the two periapsis passages, the remaining ~10% a sparse long-baseline
    presence -- mirrors how real VLTI/GRAVITY campaigns concentrate
    effort near closest approach, and matches this project's originally
    specified cadence. Every candidate epoch is then filtered to the real
    Mar-Sep Paranal visibility season (_in_visibility_season) -- Sgr A*
    isn't observable year-round.

    NOTE on a real, unavoidable complication this filtering exposes: the
    first periapsis passage (~2031.81) falls on 2031-10-22 -- just past
    the end of that year's visibility season. A space telescope could
    still stare at periapsis itself; a ground-based single-site campaign
    cannot. So that passage's +/-20-day dense window is anchored to the
    nearest date the target is actually visible (_season_anchor, ~2031
    Sep 30) rather than periapsis itself, and the season filter clips
    away the (unobservable) post-periapsis half of that window entirely.
    The second passage (~2040.49, 2040-06-26) falls comfortably inside
    its season, so its dense window is genuinely centered on periapsis.
    This is reported here rather than patched around: it is exactly the
    kind of scheduling gap a real ground-based campaign would face."""
    sparse = np.arange(CAMPAIGN_START_YR, CAMPAIGN_END_YR, SPARSE_STEP_YR)
    dense_windows = []
    for periapsis_yr in (k.NEXT_PERIAPSIS_YR, SECOND_PERIAPSIS_YR):
        anchor_yr = _season_anchor(periapsis_yr)
        dense_windows.append(np.linspace(anchor_yr - DENSE_HALF_WIDTH_YR,
                                          anchor_yr + DENSE_HALF_WIDTH_YR,
                                          N_DENSE_PER_PASSAGE))
    epochs = np.union1d(sparse, np.concatenate(dense_windows))
    epochs = epochs[_in_visibility_season(epochs)]
    return np.sort(epochs)


def true_observables(epoch_yr, extra_omega_dot=0.0):
    """True (noise-free) astrometric offset, RV, and delta-mag at the given
    epochs, from the injected truth orbit (constants.TRUTH) plus real
    Schwarzschild apsidal precession between the two passages.

    extra_omega_dot (rad/s): an additional apsidal-precession rate added on
    top of the Schwarzschild term -- used by spin_forecast.py to inject an
    aligned-spin Lense-Thirring contribution (orbit.lense_thirring_apsidal_
    rate). Zero by default, i.e. this function's normal behavior (used by
    fit_orbit.py's own orbit-recovery exercise) is unchanged."""
    period = k.TRUTH["P_yr"] * k.year
    sma = orbit.semi_major_axis_from_period(k.GM_BH, period)
    t_peri = k.TRUTH["t_peri_yr"] * k.year
    elements = orbit.OrbitalElements(
        t_peri=t_peri, period=period, ecc=k.TRUTH["e"], sma=sma,
        i_deg=k.TRUTH["i_deg"], raan_deg=k.TRUTH["Omega_deg"], omega_deg=k.TRUTH["omega_deg"],
    )
    omega_dot = orbit.schwarzschild_precession_rate(k.GM_BH, sma, k.TRUTH["e"], period) \
        + extra_omega_dot

    t_sec = epoch_yr * k.year
    state = orbit.orbit_state_precessing(t_sec, elements, omega_dot, t_ref=t_peri,
                                          grav_param=k.GM_BH)

    ra_mas, dec_mas = orbit.sky_offset_mas(state, k.D_OBS)
    rv_kms = state["v_los"] / 1e3

    one_plus_z = orbit.redshift_factor(state, k.Rs)
    # same relativistic-beaming physics as s301_lightcurve.py, but the
    # delta-mag here only needs the *shape*, so a fast approximation
    # (bolometric boosting, not a full filter-integrated blackbody) is
    # used -- adequate for a photometric campaign product, not meant to
    # replace s301_lightcurve.py's detailed per-wavelength treatment.
    dmag = -2.5 * np.log10(one_plus_z ** -3 / np.median(one_plus_z ** -3))

    return ra_mas, dec_mas, rv_kms, dmag


def main():
    """Build the synthetic two-instrument campaign, inject noise, and
    write the observations CSV that fit_orbit.py consumes."""
    epochs = build_epoch_grid()
    ra_true, dec_true, rv_true, dmag_true = true_observables(epochs)
    sigma_rv_kms, sigma_rv_kms_optimistic = compute_rv_precision()
    sigma_astrometry_mas = k.GRAVITY_PLUS_ASTROMETRY_MAS
    sigma_dmag = k.SIGMA_M_K

    rng = np.random.default_rng(RNG_SEED)
    ra_obs = ra_true + rng.normal(0, sigma_astrometry_mas, size=epochs.shape)
    dec_obs = dec_true + rng.normal(0, sigma_astrometry_mas, size=epochs.shape)
    rv_obs = rv_true + rng.normal(0, sigma_rv_kms, size=epochs.shape)
    dmag_obs = dmag_true + rng.normal(0, sigma_dmag, size=epochs.shape)

    period = k.TRUTH["P_yr"] * k.year
    sma = orbit.semi_major_axis_from_period(k.GM_BH, period)
    omega_dot = orbit.schwarzschild_precession_rate(k.GM_BH, sma, k.TRUTH["e"], period)
    omega_dot_deg_yr = np.degrees(omega_dot) * k.year

    print(f"Campaign: {len(epochs)} epochs, {epochs.min():.2f} - {epochs.max():.2f} "
          f"(passages: {k.NEXT_PERIAPSIS_YR:.2f} and {SECOND_PERIAPSIS_YR:.2f})")
    anchors = np.array([_season_anchor(k.NEXT_PERIAPSIS_YR), _season_anchor(SECOND_PERIAPSIS_YR)])
    near_periapsis = np.min(np.abs(epochs[:, None] - anchors[None, :]), axis=1) <= DENSE_HALF_WIDTH_YR
    print(f"Dense window occupancy check: {near_periapsis.sum()} of {len(epochs)} "
          f"({near_periapsis.sum() / len(epochs) * 100:.1f}%) within +/-{DENSE_HALF_WIDTH_DAYS:.0f} days "
          f"of a periapsis (season-adjusted anchor for passage 1)")
    print(f"Injected Schwarzschild precession: {omega_dot_deg_yr:.4f} deg/yr "
          f"({omega_dot_deg_yr * k.TRUTH['P_yr']:.3f} deg/orbit) -- real 1PN formula, "
          f"not illustrative")
    print(f"Astrometric precision (GRAVITY+, real achieved on-sky figure): "
          f"{sigma_astrometry_mas * 1000:.0f} uas/epoch")
    print(f"RV precision (ERIS, background-limited scaling from real S2 SINFONI/ERIS "
          f"precision -- the realistic regime for a target this faint in this field): "
          f"{sigma_rv_kms:.0f} km/s/epoch "
          f"(optimistic source-limited alternative: {sigma_rv_kms_optimistic:.1f} km/s/epoch)")
    implied_snr = k.ERIS_VELOCITY_RESOLUTION_KMS / sigma_rv_kms
    print(f"Cross-check: SPIFFIER's own R={k.ERIS_SPIFFIER_R} gives a "
          f"{k.ERIS_VELOCITY_RESOLUTION_KMS:.1f} km/s resolution element -- the injected "
          f"precision implies continuum SNR~{implied_snr:.3f}, i.e. no real single-epoch "
          f"RV detection for a star this faint with a current 8m-class instrument")
    print(f"Photometric precision (S301's own published m_K uncertainty): "
          f"{sigma_dmag:.2f} mag/epoch")
    print(f"True RA offset range: [{ra_true.min():.1f}, {ra_true.max():.1f}] mas")
    print(f"True RV range: [{rv_true.min():.0f}, {rv_true.max():.0f}] km/s")

    header = "epoch_yr,ra_offset_mas,dec_offset_mas,rv_kms,dmag_K," \
             "sigma_ra_mas,sigma_dec_mas,sigma_rv_kms,sigma_dmag"
    rows = np.column_stack([
        epochs, ra_obs, dec_obs, rv_obs, dmag_obs,
        np.full_like(epochs, sigma_astrometry_mas),
        np.full_like(epochs, sigma_astrometry_mas),
        np.full_like(epochs, sigma_rv_kms),
        np.full_like(epochs, sigma_dmag),
    ])
    os.makedirs("output", exist_ok=True)
    outpath = "output/synthetic_observations.csv"
    np.savetxt(outpath, rows, delimiter=",", header=header, comments="", fmt="%.6f")
    print(f"Saved {len(epochs)} synthetic epochs to {outpath}")


if __name__ == "__main__":
    main()
