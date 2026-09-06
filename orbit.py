"""
Shared Keplerian + relativistic orbital-mechanics engine for the S301
simulation suite (s301_lightcurve.py, campaign.py, fit_orbit.py).

Physics: full eccentric-anomaly Kepler solve (Newton-Raphson) -> position
and velocity in the orbital plane -> standard Thiele-Innes-style 3D
rotation via (inclination, RAAN, argument of periapsis) -> line-of-sight
velocity -> exact special-relativistic Doppler factor (from the star's
true 3D speed, not a small-v approximation) combined with Schwarzschild
gravitational redshift sqrt(1 - Rs/r).

Roemer (light-travel-time) delay: orbit_state_observed() solves for the
emission time behind each arrival time self-consistently (Newton), and
campaign.py / fit_orbit.py use it by default, matching the discovery
paper's own S301 orbit fit, which includes the Roemer effect.

NOT included: Kerr (spinning-BH) corrections and strong-field ray
tracing; weak-lensing astrometric deflection is quantified separately
in lensing_error.py as a systematic, not applied in the model.

Sign/orientation convention: standard orbital elements (inclination,
RAAN, argument of periapsis) as used in visual-binary / exoplanet
astrometry. Observer is along +Z; +Z is defined as *toward* the observer,
so increasing Z means decreasing distance (approaching). v_los is
defined positive = receding (astronomy convention), so v_los = -vz.
"""

from typing import NamedTuple

import numpy as np

C_LIGHT = 2.998e8   # m/s -- orbit.py is a standalone, S301-independent engine
                     # (deliberately doesn't import constants.py), so this is
                     # its own copy; must match constants.c.


class OrbitalElements(NamedTuple):
    """Bundled Keplerian elements, kept together so orbit_state() takes a
    manageable number of arguments."""
    t_peri: float
    period: float
    ecc: float
    sma: float
    i_deg: float
    raan_deg: float
    omega_deg: float


def solve_kepler(mean_anomaly, ecc, tol=1e-12, max_iter=200):
    """Bisection-safeguarded Newton-Raphson solve of Kepler's equation
    M = E - e*sin(E) for E.

    mean_anomaly: array (radians), any real value (wrapped internally).
    ecc: scalar eccentricity, 0 <= ecc < 1.

    Plain, unsafeguarded Newton-Raphson from E0=M -- this function's
    original implementation -- CATASTROPHICALLY DIVERGES for a subset of
    ordinary (M, ecc) pairs at high eccentricity: verified directly for
    S301's own e=0.9832 at M=0.2185 rad, an unremarkable value, where the
    naive iteration blows up to E~1e13 and never recovers even given the
    full 200-iteration budget, rather than converging to the true root
    near E=1.1. This wasn't visible in this project's own campaign
    (every real epoch's mean anomaly happened to avoid the divergent
    region for the true orbital elements -- checked directly, residuals
    ~1e-15), but a nonlinear optimizer exploring nearby trial parameters
    (exactly what fitting does) can and does land on such values,
    silently returning wrong positions with no error raised. This is a
    well-documented failure mode of undamped Newton-Raphson on Kepler's
    equation at high eccentricity, not specific to any one script here.

    The fix: g(E) = E - e*sin(E) - M is strictly monotonic increasing
    (g'(E) = 1 - e*cos(E) > 0 for e<1), so a bracket containing the root
    is known a priori from M = E - e*sin(E): E lies in [M-e, M+e] since
    |e*sin(E)| <= e. Each iteration takes a Newton step when it stays
    within the current bracket (fast, quadratic convergence in the
    common case) and falls back to bisection when it would not (which
    can never fail to make progress) -- a standard, globally convergent
    safeguard that costs nothing in the well-behaved case verified
    elsewhere in this project (identical results to machine precision)
    and fixes the divergent one."""
    mean_anomaly = np.asarray(mean_anomaly, dtype=float)
    wrapped = np.mod(mean_anomaly + np.pi, 2 * np.pi) - np.pi

    lo = wrapped - ecc
    hi = wrapped + ecc
    ecc_anom = wrapped.copy()
    converged = np.zeros_like(ecc_anom, dtype=bool)

    for _ in range(max_iter):
        active = ~converged
        if not np.any(active):
            break
        e_act = ecc_anom[active]
        wrapped_act = wrapped[active]
        lo_act, hi_act = lo[active], hi[active]

        resid = e_act - ecc * np.sin(e_act) - wrapped_act
        slope = 1 - ecc * np.cos(e_act)
        newton_candidate = e_act - resid / slope
        # Reject a Newton step that leaves the bracket OR fails to narrow it
        # at all (the degenerate case where the bracket's own midpoint
        # equals the current point, e.g. on the very first iteration,
        # since lo/hi are symmetric about the initial guess) -- either way,
        # fall back to bisection, which always narrows the bracket.
        bisect_candidate = 0.5 * (lo_act + hi_act)
        in_bracket = (newton_candidate > lo_act) & (newton_candidate < hi_act)
        candidate = np.where(in_bracket, newton_candidate, bisect_candidate)

        g_candidate = candidate - ecc * np.sin(candidate) - wrapped_act
        # g monotonic increasing: g(candidate) <= 0 -> root is >= candidate (new lo);
        # g(candidate) >= 0 -> root is <= candidate (new hi).
        neg = g_candidate <= 0
        lo[active] = np.where(neg, candidate, lo_act)
        hi[active] = np.where(neg, hi_act, candidate)
        ecc_anom[active] = candidate

        # Convergence via bracket width, not step size: a bisection step
        # can legitimately leave the point unmoved-looking in edge cases
        # while still having narrowed (or being about to narrow) the
        # bracket, so bracket width is the robust criterion here.
        newly_converged = (hi[active] - lo[active]) < tol
        idx_active = np.flatnonzero(active)
        converged[idx_active[newly_converged]] = True

    # shift back by the same number of 2*pi windings removed above
    return ecc_anom + (mean_anomaly - wrapped)


def eccentric_anomaly_grid(n_points):
    """Uniform-in-E sampling: dt/dE = (1 - e*cos E) * P/2pi is small near
    periapsis (E=0) for high e, so this naturally concentrates time
    resolution right where the fast relativistic stuff happens. Useful for
    smooth plots; NOT used for fitting arbitrary observation epochs (use
    solve_kepler on their actual times instead)."""
    return np.linspace(-np.pi, np.pi, n_points)


def time_since_periapsis(ecc_anom, ecc, period):
    """Convert eccentric anomaly to time-since-periapsis (same units as
    period)."""
    mean_anom = ecc_anom - ecc * np.sin(ecc_anom)
    return mean_anom * period / (2 * np.pi)


def _plane_position_velocity(ecc_anom, ecc, sma, grav_param):
    """Position/velocity in the orbital plane, before 3D rotation."""
    radius = sma * (1 - ecc * np.cos(ecc_anom))
    true_anom = 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(ecc_anom / 2),
                                np.sqrt(1 - ecc) * np.cos(ecc_anom / 2))
    semi_latus = sma * (1 - ecc ** 2)

    x_plane = radius * np.cos(true_anom)
    y_plane = radius * np.sin(true_anom)
    vx_plane = -np.sqrt(grav_param / semi_latus) * np.sin(true_anom)
    vy_plane = np.sqrt(grav_param / semi_latus) * (ecc + np.cos(true_anom))
    return radius, true_anom, x_plane, y_plane, vx_plane, vy_plane


def _rotation_matrix(i_deg, raan_deg, omega_deg):
    """3x2 Thiele-Innes rotation matrix: rows are sky X, sky Y, line-of-sight
    Z; columns are the orbital-plane x, y axes. omega_deg may be an array
    (one value per time sample, for a precessing argument of periapsis) --
    i_deg/raan_deg are then broadcast against it; the result has shape
    (3, 2) for scalar inputs or (3, 2, N) if omega_deg has shape (N,)."""
    incl = np.radians(i_deg)
    raan = np.radians(raan_deg)
    arg_peri = np.radians(omega_deg)
    cos_raan, sin_raan = np.cos(raan), np.sin(raan)
    cos_peri, sin_peri = np.cos(arg_peri), np.sin(arg_peri)
    cos_i, sin_i = np.cos(incl), np.sin(incl)
    return np.array([
        [cos_raan * cos_peri - sin_raan * sin_peri * cos_i,
         -cos_raan * sin_peri - sin_raan * cos_peri * cos_i],
        [sin_raan * cos_peri + cos_raan * sin_peri * cos_i,
         -sin_raan * sin_peri + cos_raan * cos_peri * cos_i],
        [sin_peri * sin_i, cos_peri * sin_i],
    ])


def _rotate_to_sky(x_plane, y_plane, vx_plane, vy_plane, i_deg, raan_deg, omega_deg):
    """Thiele-Innes-style rotation of orbital-plane position/velocity into
    the sky frame: X,Y = sky plane; Z = line of sight, +Z toward observer.
    i_deg, raan_deg, omega_deg are scalars (static orbit)."""
    rot = _rotation_matrix(i_deg, raan_deg, omega_deg)
    sky_pos = np.tensordot(rot, np.stack([x_plane, y_plane]), axes=([1], [0]))
    sky_vel = np.tensordot(rot, np.stack([vx_plane, vy_plane]), axes=([1], [0]))
    return sky_pos[0], sky_pos[1], sky_pos[2], sky_vel[0], sky_vel[1], sky_vel[2]


def _rotate_to_sky_precessing(x_plane, y_plane, vx_plane, vy_plane, i_deg, raan_deg, omega_deg):
    """Same rotation as _rotate_to_sky, but omega_deg is an array (one
    value per time sample, for a precessing argument of periapsis) --
    uses an elementwise batched matrix-vector product instead of
    _rotate_to_sky's shared-matrix tensordot."""
    rot = _rotation_matrix(i_deg, raan_deg, omega_deg)
    sky_pos = np.einsum("ijn,jn->in", rot, np.stack([x_plane, y_plane]))
    sky_vel = np.einsum("ijn,jn->in", rot, np.stack([vx_plane, vy_plane]))
    return sky_pos[0], sky_pos[1], sky_pos[2], sky_vel[0], sky_vel[1], sky_vel[2]


def _finalize_state(radius, true_anom, sky_pos, sky_vel, plane_speed):
    """Bundle the final orbit_state()/orbit_state_precessing() return dict
    from the already-rotated sky position/velocity."""
    sky_x, sky_y, sky_z = sky_pos
    sky_vx, sky_vy, sky_vz = sky_vel
    v_los = -sky_vz  # +Z toward observer => approaching; flip so + = receding
    light_speed = 2.998e8
    beta = plane_speed / light_speed
    beta_r = v_los / light_speed
    gamma = 1 / np.sqrt(1 - beta ** 2)
    return {
        "r": radius, "true_anom": true_anom,
        "sky_x": sky_x, "sky_y": sky_y, "sky_z": sky_z,
        "sky_vx": sky_vx, "sky_vy": sky_vy, "sky_vz": sky_vz,
        "speed": plane_speed, "v_los": v_los, "beta": beta, "beta_r": beta_r, "gamma": gamma,
    }


def orbit_state(t, elements, grav_param):
    """
    Full orbital state at times t, given bundled Keplerian elements.

    t, elements.t_peri, elements.period: consistent time units (e.g. all
    seconds, or all years).
    elements.sma: semi-major axis, consistent length unit with the
    returned position/velocity.
    grav_param: G*M_central, matching the length/time units of sma and
    period.

    Returns dict: r, true_anom, sky_x, sky_y, sky_z (position; x,y=sky
    plane, z=line of sight toward observer), sky_vx, sky_vy, sky_vz
    (velocity), speed, v_los (+ = receding), beta, beta_r, gamma.
    """
    t = np.asarray(t, dtype=float)
    mean_anom = 2 * np.pi * (t - elements.t_peri) / elements.period
    ecc_anom = solve_kepler(mean_anom, elements.ecc)

    radius, true_anom, x_plane, y_plane, vx_plane, vy_plane = _plane_position_velocity(
        ecc_anom, elements.ecc, elements.sma, grav_param)

    sky_x, sky_y, sky_z, sky_vx, sky_vy, sky_vz = _rotate_to_sky(
        x_plane, y_plane, vx_plane, vy_plane,
        elements.i_deg, elements.raan_deg, elements.omega_deg)

    speed = np.sqrt(vx_plane ** 2 + vy_plane ** 2)
    return _finalize_state(radius, true_anom, (sky_x, sky_y, sky_z),
                            (sky_vx, sky_vy, sky_vz), speed)


def orbit_state_precessing(t, elements, omega_dot, t_ref, grav_param):
    """
    Same physics as orbit_state(), but the argument of periapsis precesses
    linearly with time: omega(t) = elements.omega_deg + omega_dot*(t -
    t_ref), i.e. elements.omega_deg is the value AT t_ref (not at
    elements.t_peri). i and raan (Omega) are NOT precessed -- this models
    real, standard-formula Schwarzschild (1PN) apsidal precession only
    (see schwarzschild_precession_rate), never a black-hole-spin/frame-
    dragging (Lense-Thirring, nodal) effect, which is genuinely unmeasured
    for Sgr A* and deliberately not modeled anywhere in this project.

    omega_dot: rad/s (matching t's units, e.g. seconds).
    Returns the same dict shape as orbit_state().
    """
    t = np.asarray(t, dtype=float)
    mean_anom = 2 * np.pi * (t - elements.t_peri) / elements.period
    ecc_anom = solve_kepler(mean_anom, elements.ecc)

    radius, true_anom, x_plane, y_plane, vx_plane, vy_plane = _plane_position_velocity(
        ecc_anom, elements.ecc, elements.sma, grav_param)

    omega_eff_deg = elements.omega_deg + np.degrees(omega_dot * (t - t_ref))
    sky_x, sky_y, sky_z, sky_vx, sky_vy, sky_vz = _rotate_to_sky_precessing(
        x_plane, y_plane, vx_plane, vy_plane,
        elements.i_deg, elements.raan_deg, omega_eff_deg)

    speed = np.sqrt(vx_plane ** 2 + vy_plane ** 2)
    return _finalize_state(radius, true_anom, (sky_x, sky_y, sky_z),
                            (sky_vx, sky_vy, sky_vz), speed)


N_ROEMER_ITER = 20
# Newton converges in ~2-3 iterations (contraction rate |v_los|/c <= 8.5%).
# 1e-3 s is ~10 orders of magnitude below any physical timescale here and
# freezes each epoch right after genuine convergence, avoiding a verified
# floating-point edge case in the Kepler solve at tighter tolerances.
ROEMER_TOL_SEC = 1e-3


def solve_emission_time(t_obs, elements, omega_dot, t_ref, grav_param):
    """Emission time behind each arrival time t_obs (seconds): solves
    g(t_emit) = t_emit - sky_z(t_emit)/c - t_obs = 0 by Newton's method,
    using g' = 1 - sky_vz/c from the same state evaluation. +Z is toward
    the observer, so light from a nearer position arrives sooner for the
    same emission time. Vectorized, with each epoch frozen once its step
    falls below ROEMER_TOL_SEC."""
    t_obs = np.asarray(t_obs, dtype=float)
    t_emit = t_obs.copy()
    converged = np.zeros_like(t_emit, dtype=bool)
    for _ in range(N_ROEMER_ITER):
        active = ~converged
        if not np.any(active):
            break
        state = orbit_state_precessing(t_emit[active], elements, omega_dot, t_ref=t_ref,
                                       grav_param=grav_param)
        g_val = t_emit[active] - state["sky_z"] / C_LIGHT - t_obs[active]
        g_prime = 1.0 - state["sky_vz"] / C_LIGHT
        step = g_val / g_prime
        t_emit[active] = t_emit[active] - step
        newly_converged = np.abs(step) < ROEMER_TOL_SEC
        idx_active = np.flatnonzero(active)
        converged[idx_active[newly_converged]] = True
    return t_emit


def orbit_state_observed(t_obs, elements, omega_dot, t_ref, grav_param, light_time=True):
    """State of the star as SEEN at arrival times t_obs (seconds): with
    light_time=True (default), the orbit is evaluated at the self-
    consistent emission time; with light_time=False it is evaluated
    naively at t_obs, which is what the light-time-unaware model in
    roemer_delay.py uses to quantify the bias of ignoring the effect."""
    t_eval = solve_emission_time(t_obs, elements, omega_dot, t_ref, grav_param) if light_time \
        else np.asarray(t_obs, dtype=float)
    return orbit_state_precessing(t_eval, elements, omega_dot, t_ref=t_ref, grav_param=grav_param)


def schwarzschild_precession_rate(grav_param, sma, ecc, period):
    """Standard 1PN (Schwarzschild) apsidal precession rate -- the same
    formula used for Mercury's perihelion advance, and the effect the real
    S301 discovery papers say has already been measured, seen as its
    "rosette" orbit shape. Real GR physics, not illustrative: computed
    purely from already-known quantities, no assumed/injected inputs.
    Returns domega/dt in rad/s (matching period's time unit)."""
    light_speed = 2.998e8
    delta_omega_per_orbit = 6 * np.pi * grav_param / (light_speed ** 2 * sma * (1 - ecc ** 2))
    return delta_omega_per_orbit / period


def lense_thirring_apsidal_rate(grav_param, sma, ecc, period, spin_chi, xi_deg=0.0):
    """Lense-Thirring (frame-dragging) contribution to the apsidal
    precession rate, for a black-hole spin at angle xi_deg from the
    orbital angular momentum. This is the REAL formula from the S301
    discovery paper's own spin-forecast methods section (arXiv:2607.12664,
    Eq. 1): d(varpi)/dt = 8*pi*chi*cos(xi)*(r_g/p)^1.5 / period, where
    r_g = GM/c^2 and p = a*(1-e^2).

    For xi_deg=0 (spin aligned with the orbital angular momentum -- the
    discovery paper's own illustrative benchmark case), this produces
    ONLY in-plane (apsidal) precession and no nodal drift -- validated
    numerically against the paper's own quoted value for S301: this
    function gives 0.114 deg/orbit at chi=1, xi=0, matching their quoted
    "0.11 deg * chi * cos(xi)" almost exactly.

    NOT implemented: the general misaligned-spin nodal (out-of-plane)
    precession, or the companion paper's (arXiv:2607.24931) Newtonian-
    confusion subtraction methodology -- both explicitly out of scope
    here; see README.md.

    Returns d(varpi)/dt in rad/s (matching period's time unit)."""
    xi = np.radians(xi_deg)
    r_g = grav_param / C_LIGHT ** 2
    p_orb = sma * (1 - ecc ** 2)
    delta_varpi_per_orbit = 8 * np.pi * spin_chi * np.cos(xi) * (r_g / p_orb) ** 1.5
    return delta_varpi_per_orbit / period


def redshift_factor(state, schwarzschild_radius):
    """Combined special-relativistic Doppler + Schwarzschild gravitational
    redshift factor (1+z), given an orbit_state() dict and Schwarzschild
    radius (same length unit as state['r'])."""
    gamma = state["gamma"]
    beta_r = state["beta_r"]
    delta_sr = 1.0 / (gamma * (1 + beta_r))  # v_r>0 (receding) -> delta<1 -> redshift
    grav_factor = np.sqrt(np.clip(1 - schwarzschild_radius / state["r"], 1e-6, None))
    return 1.0 / (delta_sr * grav_factor)


def semi_major_axis_from_period(grav_param, period):
    """Kepler's third law: a^3 = GM*P^2 / (4*pi^2). DERIVED, not assumed --
    the whole point of specifying an orbit by (period, ecc) rather than
    (sma, ecc) directly is that sma then comes out as a physical
    consequence, not an independent input that could silently disagree
    with the period via the assumed grav_param."""
    return (grav_param * period ** 2 / (4 * np.pi ** 2)) ** (1 / 3)


def sky_offset_mas(state, distance):
    """Angular (RA, Dec) offset in mas from an orbit_state()/
    orbit_state_precessing() dict, at observer distance `distance` (same
    length unit as state['r'], e.g. meters)."""
    ra_mas = (state["sky_x"] / distance) * (180 / np.pi) * 3600 * 1000
    dec_mas = (state["sky_y"] / distance) * (180 / np.pi) * 3600 * 1000
    return ra_mas, dec_mas
