"""
Physical constants and the real, published S301 / Sgr A* catalog values
used throughout this project (s301_lightcurve.py, campaign.py,
fit_orbit.py). Single source of truth so nothing drifts between scripts.

S301 is a genuinely real star, announced August 2026 (ESO press release
eso2612; MPE presskit) and published in Nature as "Discovery of a star
sensitive to the spin of Sgr A*" (arXiv:2607.12664,
DOI 10.1038/s41586-026-10894-w). Orbital elements below are Extended
Data Table 2, Solution 1 of that paper. A second, near-degenerate mirror
solution (SOLUTION_2 below) is also published -- a known ambiguity in
astrometric-only orbit fits -- and is kept here for reference, not used
as the simulation's injected truth.

Sources:
  - https://www.eso.org/public/news/eso2612/
  - https://www.mpe.mpg.de/8219735/news20260819
  - https://arxiv.org/abs/2607.24931 ("S301 and friends")
"""

# ---------------------------------------------------------------- physical constants (SI)
G = 6.674e-11          # m^3 kg^-1 s^-2
c = 2.998e8            # m/s
Msun = 1.989e30        # kg
Rsun = 6.957e8         # m
AU = 1.496e11          # m
pc = 3.086e16          # m
h = 6.626e-34          # J s
kB = 1.381e-23         # J/K
sigma_SB = 5.670e-8    # W m^-2 K^-4
year = 3.15576e7       # s (Julian year)

# ---------------------------------------------------------------- Sgr A* / Galactic Center
M_BH = 4.297e6 * Msun          # arXiv:2607.12664
GM_BH = G * M_BH
Rs = 2 * GM_BH / c ** 2         # Schwarzschild radius
R0_PC = 8277.0                  # distance to Sgr A*, pc (arXiv:2607.12664; "8.3 kpc")

# ---------------------------------------------------------------- S301 orbital elements
# Published uncertainties (sigma_*) are kept alongside each value for
# comparison against this project's own bootstrap-fit uncertainties later.
SOLUTION_1 = {
    "a_mas": 83.0, "sigma_a_mas": 0.7,
    "e": 0.9832, "sigma_e": 0.0010,
    "i_deg": 124.09, "sigma_i_deg": 1.10,
    "Omega_deg": 73.8, "sigma_Omega_deg": 3.5,
    "omega_deg": 293.4, "sigma_omega_deg": 2.2,
    "t_peri_yr": 2023.126, "sigma_t_peri_yr": 0.010,
    "P_yr": 8.68, "sigma_P_yr": 0.11,
    "r_p_Rs": 136,
}

SOLUTION_2 = {
    "a_mas": 83.0, "sigma_a_mas": 0.7,
    "e": 0.9824, "sigma_e": 0.0011,
    "i_deg": 122.84, "sigma_i_deg": 1.12,
    "Omega_deg": 256.9, "sigma_Omega_deg": 3.3,
    "omega_deg": 115.1, "sigma_omega_deg": 2.0,
    "t_peri_yr": 2023.125, "sigma_t_peri_yr": 0.011,
    "P_yr": 8.68, "sigma_P_yr": 0.11,
    "r_p_Rs": 142,
}

# This project's injected "truth" for the campaign/fitting exercise.
TRUTH = SOLUTION_1

# ~2031.8; ESO stated "next closest passage in 2031"
NEXT_PERIAPSIS_YR = TRUTH["t_peri_yr"] + TRUTH["P_yr"]

# ---------------------------------------------------------------- stellar properties
# Spectral type F1.5V and m_K=19.3 are published (arXiv:2607.12664); T_eff
# is NOT published there, so it's a standard-table estimate for an F1.5
# main-sequence star, flagged explicitly as inferred, not measured.
SPECTRAL_TYPE = "F1.5V"
M_K_APPARENT = 19.3
SIGMA_M_K = 0.3
MASS_RANGE_MSUN = (1.1, 1.5)     # published, age-dependent
RADIUS_RANGE_RSUN = (1.4, 1.6)   # published
T_EFF_K = 7300.0                 # NOT published -- standard F1.5V table estimate
MASS_MSUN = sum(MASS_RANGE_MSUN) / 2
R_STAR = (sum(RADIUS_RANGE_RSUN) / 2) * Rsun

MAX_SPEED_KMS = 25000.0   # published; >8% c

# ---------------------------------------------------------------- observer: ESO Paranal Observatory
# All simulated observations in this project are ground-based, from ESO's
# Paranal Observatory, Chile (home of both the VLT and VLTI) -- the same
# site GRAVITY used for the real S301 discovery. Earth's position on its
# orbit shifts by <<1 AU against the 8277 pc distance to Sgr A* -- utterly
# negligible -- so d_obs = R0 to many significant figures.
D_OBS_PC = R0_PC
D_OBS = D_OBS_PC * pc

# GRAVITY/VLTI's actual K-band operating range (not a generic photometric
# filter): 1.98-2.40 um, covering the full atmospheric K window. This is
# the real instrument's own bandpass, used both for S301's published
# m_K=19.3 and for the s301_lightcurve.py light-curve simulation.
# Source: GRAVITY Collaboration, Abuter, R., et al. 2017, A&A, 602, A94,
# "First light for GRAVITY: Phase referencing optical interferometry for
# the Very Large Telescope Interferometer" (arXiv:1705.02345).
# CORRECTED 2026-09-05: previously cited only vaguely as "ESO VLTI/
# GRAVITY instrument page" with no locatable source -- caught during the
# same citation audit that found the SINFONI mis-citation above.
KS_RANGE = (1.98e-6, 2.40e-6)     # m
KS_LAMBDA0 = sum(KS_RANGE) / 2    # m, ~2.19 um
KS_FWHM = KS_RANGE[1] - KS_RANGE[0]   # m, ~0.42 um

# ---------------------------------------------------------------- astrometry: GRAVITY+/VLTI
# Radial velocity has exactly zero sensitivity to the ascending node
# Omega (see orbit.py's rotation matrix -- the line-of-sight row never
# involves raan_deg), so no spectroscopy/photometry-only campaign can
# ever constrain it: astrometry is genuinely required. This project
# attributes it to GRAVITY+, the VLTI interferometer upgrade explicitly
# designed for exactly this kind of faint, close-in S-star follow-up --
# the same class of instrument (an upgrade of GRAVITY itself) that
# discovered S301 in the first place. This is a REAL, on-sky achieved
# performance figure, not a theoretical one: GRAVITY has already
# demonstrated 30-100 uas astrometry on sources down to m_K~20 (combining
# a full night). S301 (m_K=19.3) sits within that demonstrated regime;
# this project uses the more conservative (larger/worse) end of the
# quoted range. Source: arXiv:2301.08071 ("The GRAVITY+ Project").
GRAVITY_PLUS_ASTROMETRY_MAS = 0.1   # 100 uas

# ---------------------------------------------------------------- photometry: same GRAVITY+ data
# S301's own published photometric uncertainty (already defined above as
# SIGMA_M_K=0.3, arXiv:2607.12664) is used directly as this project's
# per-epoch photometric precision -- a real, already-sourced number,
# rather than a freshly derived one. This is the paper's ABSOLUTE
# calibration uncertainty on a single measurement, used here as a
# (probably conservative) stand-in for relative/differential monitoring
# precision between epochs, since no dedicated repeatability figure for
# GRAVITY+ photometry was found.

# ---------------------------------------------------------------- radial velocity: ERIS/VLT
# Attributed to ERIS (VLT UT4's adaptive-optics-fed SPIFFIER integral
# field spectrograph) via the Brackett-gamma line -- the real instrument
# (and successor to SINFONI) used for essentially all S-cluster star RV
# monitoring to date. R=5000 is SPIFFIER's full-K-band-coverage mode
# (as opposed to its R~8000-11000 half-band modes): the full band is
# needed here because S301's own velocity swings by ~8% c, enough to
# Doppler-shift a line most of the way across a half-band setting.
# Source: Davies, R. I., et al. 2023, A&A, 674, A207, "The Enhanced
# Resolution Imager and Spectrograph for the VLT" (arXiv:2304.02343).
ERIS_SPIFFIER_R = 5000
ERIS_VELOCITY_RESOLUTION_KMS = (c / ERIS_SPIFFIER_R) / 1e3   # ~60 km/s per resolution element

# RV precision is NOT derived from a theoretical sensitivity table (as an
# earlier version of this project did for a different, since-abandoned
# instrument choice) -- it's scaled from a real ACHIEVED on-sky result:
# SINFONI radial velocities of S2 (the best-studied, brightest S-star)
# reach a median precision of 12.3 km/s.
# Source: GRAVITY Collaboration, Abuter, R., et al. 2018, A&A, 615, L15,
# "Detection of the gravitational redshift in the orbit of the star S2
# near the Galactic centre massive black hole" (arXiv:1807.09409).
# CORRECTED 2026-09-05: previously mis-cited as arXiv:1709.01598, which
# is an unrelated paper (Nishiyama et al., Subaru/IRCS RV of S2,
# reporting a different, if similar-order, ~13-17 km/s precision) --
# caught during a final citation audit; the 12.3 km/s VALUE itself was
# already correct, only its citation was wrong.
SIGMA_RV_S2_KMS = 12.3
M_K_S2 = 14.0   # S2's real, well-established apparent K magnitude

# Scaling that reference to S301's much fainter m_K=19.3: photon-noise
# precision degrades as some power of the flux ratio, depending on
# whether the noise is source-photon-limited (sqrt(flux) scaling,
# precision ~ 10**(0.2*delta_m)) or background/sky-photon-limited
# (linear-in-flux scaling, precision ~ 10**(0.4*delta_m)).
#
# This project adopts the LINEAR (background-limited) scaling as its
# PRIMARY, actually-injected estimate: ground-based K-band IFU
# spectroscopy of faint Galactic Center targets is well documented in the
# literature as background-dominated (bright, variable OH airglow lines,
# thermal K-band sky background, and severe stellar crowding in this
# specific field) rather than source-photon-starved -- fainter targets in
# this regime are known to need much longer integration precisely because
# the background noise floor doesn't improve with the target's own flux.
# The source-limited (sqrt) scaling is kept as an optimistic alternative
# for comparison, but is NOT what's injected into the synthetic campaign.
_DELTA_M_S301_S2 = M_K_APPARENT - M_K_S2
SIGMA_RV_S301_KMS_LINEAR_SCALING = SIGMA_RV_S2_KMS * 10 ** (0.4 * _DELTA_M_S301_S2)
SIGMA_RV_S301_KMS_SQRT_SCALING = SIGMA_RV_S2_KMS * 10 ** (0.2 * _DELTA_M_S301_S2)

# HARMONI/ELT RV feasibility, extending the SAME background-limited
# regime above -- NOT a new assumption, a consistency extension of it.
# In that regime (background flux per resolution element set by the
# diffraction-limited PSF core, source flux set by collecting area),
# S/N scales as aperture area (D^2) rather than sqrt(area), because a
# larger, equally diffraction-limited aperture both collects more source
# photons AND concentrates them into a smaller solid angle, collecting
# proportionally less background -- the same physical mechanism already
# invoked above to justify linear-in-flux (not sqrt) scaling with
# magnitude. Velocity precision (~ resolution element / S/N) then scales
# as 1/(R * D^2). D_ELT=39 m, D_ERIS=8.2 m (VLT UT); R_HARMONI=17385 is
# the closest real match to GRAVITY's own K-band range (1.98-2.40 um),
# using HARMONI's high-resolution K-short + K-long modes together
# (Thatte, N., et al. 2024, Proc. SPIE, 13096, 1309614, "HARMONI at ELT:
# project status and instrument overview," doi:10.1117/12.3018520).
#
# UNLIKE every other RV number in this project, this one scales a real
# ACHIEVED result to a hypothetical FUTURE instrument, not to a fainter
# target on already-operating hardware -- HARMONI has not observed
# anything yet, so its actual achieved Strehl ratio, sky background, and
# AO correction quality remain unverified assumptions, not measurements.
# Treat this as indicative of the right order of magnitude, not a
# proposal-ready sensitivity estimate.
D_ERIS_M = 8.2
D_ELT_M = 39.0
HARMONI_K_HIGH_R = 17385
_HARMONI_SCALING_VS_ERIS = (ERIS_SPIFFIER_R / HARMONI_K_HIGH_R) * (D_ERIS_M / D_ELT_M) ** 2
SIGMA_RV_S301_KMS_HARMONI_ELT = SIGMA_RV_S301_KMS_LINEAR_SCALING * _HARMONI_SCALING_VS_ERIS
