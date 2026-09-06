"""
Numerical estimate of the "Newtonian confusion" (extended-mass) systematic
error term for S301's measured apsidal precession -- part of the
comprehensive error budget (instrumental + astrophysical + random) the
project was always meant to include, alongside the random-noise treatment
in fit_orbit.py's bootstrap.

Methodology is the actual post-Newtonian treatment of Rubilar & Eckart
(2001, A&A 374, 95; hereafter RE01), not a recalled/reconstructed formula
-- read from the primary-source PDF in full (Sects. 3, Appendix A) after
an earlier attempt to recall the result from memory produced a formula
that does not appear anywhere in that paper. Specifically:

  - Eq. (9):  Plummer-like density shape rho_alpha(r) = [1+(r/rc)^2]^(-alpha/2)
  - Eq. (10)-(12): rho(r) and enclosed mass M(r), normalized so that
    M(r) -> M as r -> infinity and M(r) -> lambda_p * M as r -> 0, with
    lambda_p the point-mass (BH) fraction and lambda_e = 1 - lambda_p the
    extended fraction.
  - Eq. (13): the *exact* Newtonian potential of a spherically symmetric
    mass distribution, phi(r) = -G * int_r^inf M(zeta)/zeta^2 dzeta
    (not an approximation -- follows directly from the shell theorem via
    Eq. (A.6) below).
  - Eq. (A.6)-(A.7): the post-Newtonian (1PN) equation of motion for a
    spherically symmetric matter distribution,
        dv/dt = -G*M(r)/r^3 * [ (1 + 4*phi/c^2 + v^2/c^2) r - (4/c^2) v (v.r) ]
    which RE01 state (end of Appendix A) reduces to the standard
    Schwarzschild periastron advance, Eq. (1) = 6*pi*G*M/(c^2*a*(1-e^2)),
    for a point mass -- this is the numerical validation check below.

This script does NOT use, assume, or require any radial-velocity data:
it is a pure forward orbital-mechanics calculation (3D position/velocity
state vectors, integrated with a modified acceleration law, periapsis
detected geometrically). RV was permanently dropped from this project's
actual observing campaign (see campaign.py / CLAUDE.md); this systematic-
error estimate concerns the true astrometric precession itself and is
independent of that decision either way.

Two integrations are run from S301's real orbital elements
(constants.TRUTH), both starting at the same true periapsis passage:
  (a) point-mass-only (BH alone) -- validates the integrator against
      orbit.schwarzschild_precession_rate()'s analytic 1PN prediction.
  (b) point mass + Plummer-like extended mass, with the extended mass
      normalized to the real GRAVITY Collaboration (2022, A&A 657, L12)
      bounds on any enclosed extended/"dark" mass component inside the
      S2 apocentre -- M_ext <~ 3000 Msun (1-sigma) and <~ 7500 Msun
      (their conservative 3-sigma limit) -- using that paper's own
      assumed Plummer scale length rc = 0.3 arcsec (~12 mpc at R0) with
      RE01's alpha=5 shape (the classical Plummer profile), so the mass
      bound and the profile it was derived under stay self-consistent.

The extended-mass-induced precession bias is the difference between the
two integrations' measured Delta-omega per orbit, compared against the
fitted omega_dot precision in results/fit_orbit.json and the injected
1PN rate (fit_orbit.TRUE_OMEGA_DOT_DEG_YR). Results go to
results/extended_mass_error.json.
"""

import numpy as np
from scipy.integrate import quad, solve_ivp

import constants as k
import fit_orbit
import orbit
import results_io

ALPHA = 5.0  # alpha=5 in RE01's parametrization IS the classical Plummer
             # profile shape, matching the profile actually assumed below.
# Core radius and extended-mass bound both taken from the same, most
# current source -- GRAVITY Collaboration (2022, A&A 657, L12, arXiv:
# 2112.07478), "Mass distribution in the Galactic Center based on
# interferometric astrometry of multiple stellar orbits" -- verified
# directly against that paper's own text (not a secondary paraphrase;
# an initial secondary-source pass returned three inconsistent numbers,
# 4000/1200/7500 Msun, from three different real papers before this one
# was read directly). Its abstract: "the extended mass component inside
# the S2 apocenter ... must be <~3000 Msun (1-sigma), or <~0.1% of M_BH",
# assuming "a Plummer density profile with scale length 0.3 arcsec".
# Using this paper's own assumed core radius together with its own mass
# bound (rather than pairing a modern bound with RE01's decades-old
# fiducial rc=5.8 mpc) keeps the two numbers self-consistent.
RC_ARCSEC = 0.3
RC = RC_ARCSEC * 4.84814e-6 * k.D_OBS   # arcsec -> rad -> m, at Sgr A*'s distance
M_EXT_1SIGMA = 3000.0 * k.Msun
M_EXT_3SIGMA = 7500.0 * k.Msun          # paper's own "conservative" 3-sigma limit
M_EXT_TOTAL = M_EXT_1SIGMA

R_QUAD_INF = 5000.0 * RC    # numerical stand-in for infinity in the RE01
                             # integrals; rho_alpha ~ r^-alpha at large r
                             # (alpha=5) so the r^2*rho integrand ~ r^-3
                             # and the tail beyond here is negligible


def _rho_alpha(r):
    return (1.0 + (r / RC) ** 2) ** (-ALPHA / 2.0)


def _shape_integrand(r):
    return r ** 2 * _rho_alpha(r)


# RE01 Eq. (11)-(12) normalization integral: int_0^inf zeta^2 rho_alpha(zeta) dzeta
_NORM, _NORM_ERR = quad(_shape_integrand, 0.0, R_QUAD_INF, limit=200)

# Shared log-spaced radial grid for precomputing M(r)/phi(r) once per mass
# model, then interpolating during integration -- calling scipy.quad from
# inside the ODE right-hand side (evaluated ~1e5+ times by an adaptive
# 8th-order integrator) was verified to be too slow (>300s, no result);
# precompute+interpolate is the standard fix and costs one quadrature pass
# instead of one per force evaluation. R_GRID_MIN is set far below any
# radius S301 reaches (periapsis ~0.056 mpc >> R_GRID_MIN).
R_GRID_MIN = 1e-6 * RC
_R_GRID = np.geomspace(R_GRID_MIN, R_QUAD_INF, 4000)


def build_enclosed_mass_fn(m_total, lambda_p):
    """RE01 Eq. (12): M(r)/M = lambda_p + lambda_e * I(r)/I(inf), precomputed
    on _R_GRID and interpolated (linear in log r) for speed. Returns a
    vectorized callable M(r) in kg."""
    lambda_e = 1.0 - lambda_p
    shape_vals = _shape_integrand(_R_GRID)
    cum_shape = np.concatenate([[0.0], np.cumsum(
        0.5 * (shape_vals[1:] + shape_vals[:-1]) * np.diff(_R_GRID))])
    frac_grid = cum_shape / _NORM
    m_grid = m_total * (lambda_p + lambda_e * frac_grid)
    log_r_grid = np.log(_R_GRID)

    def m_of_r(r):
        r = np.atleast_1d(np.asarray(r, dtype=float))
        return np.interp(np.log(np.clip(r, R_GRID_MIN, R_QUAD_INF)), log_r_grid, m_grid)

    return m_of_r


def build_potential_fn(m_of_r, m_total, r_max=R_QUAD_INF):
    """RE01 Eq. (13): phi(r) = -G * int_r^inf M(zeta)/zeta^2 dzeta, exact
    for any spherically symmetric mass distribution (shell theorem).
    Precomputed on _R_GRID via cumulative trapezoid (same speed rationale
    as build_enclosed_mass_fn) plus an analytic tail for [r_max, inf)
    where M(zeta) has already converged to m_total:
    int_{r_max}^inf m_total/zeta^2 dzeta = m_total/r_max."""
    m_grid = m_of_r(_R_GRID)
    g_grid = m_grid / _R_GRID ** 2
    cum_g = np.concatenate([[0.0], np.cumsum(
        0.5 * (g_grid[1:] + g_grid[:-1]) * np.diff(_R_GRID))])  # int_{Rmin}^r g dzeta
    tail_analytic = m_total / r_max
    total_to_rmax = cum_g[-1]  # int_{Rmin}^{r_max} g dzeta
    # int_r^inf M/zeta^2 dzeta = [int_{Rmin}^{r_max} - int_{Rmin}^r] + tail
    phi_grid = -k.G * (total_to_rmax - cum_g + tail_analytic)
    log_r_grid = np.log(_R_GRID)

    def phi_of_r(r):
        r = np.atleast_1d(np.asarray(r, dtype=float))
        out = np.interp(np.log(np.clip(r, R_GRID_MIN, R_QUAD_INF)), log_r_grid, phi_grid)
        out = np.where(r >= r_max, -k.G * m_total / r, out)
        return out

    return phi_of_r


def make_accel_fn(m_of_r, phi_of_r):
    """RE01 Eq. (A.7), full post-Newtonian acceleration for a spherically
    symmetric mass distribution. Falls back exactly to Eq. (A.8) when
    m_of_r/phi_of_r describe a pure point mass."""

    def accel(r_vec, v_vec):
        r = np.linalg.norm(r_vec)
        v2 = float(np.dot(v_vec, v_vec))
        m_enc = float(m_of_r(np.array([r]))[0])
        phi = float(phi_of_r(np.array([r]))[0])
        rdotv = float(np.dot(r_vec, v_vec))
        bracket = (1.0 + 4.0 * phi / k.c ** 2 + v2 / k.c ** 2) * r_vec \
            - (4.0 / k.c ** 2) * v_vec * rdotv
        return -k.G * m_enc / r ** 3 * bracket

    return accel


def integrate_orbits(accel_fn, state0, period, n_orbits=3):
    """Integrate dv/dt = accel_fn(r,v), dx/dt = v, from state0 = (r_vec,
    v_vec) at t=0 for n_orbits periods, detecting periapsis passages via
    the r.v sign change (negative -> positive)."""
    r0, v0 = state0
    y0 = np.concatenate([r0, v0])

    def rhs(_t, y):
        r_vec, v_vec = y[:3], y[3:]
        a_vec = accel_fn(r_vec, v_vec)
        return np.concatenate([v_vec, a_vec])

    def periapsis_event(_t, y):
        return float(np.dot(y[:3], y[3:]))
    periapsis_event.direction = 1.0

    t_span = (0.0, n_orbits * period)
    sol = solve_ivp(rhs, t_span, y0, method="DOP853", events=periapsis_event,
                     rtol=1e-12, atol=1e-2, max_step=period / 20000.0, dense_output=False)
    return sol


def eccentricity_vector(r_vec, v_vec, grav_param):
    """Standard Laplace-Runge-Lenz / eccentricity vector, valid for any
    (possibly perturbed) two-body-like state -- used here purely as a
    convention-independent way to measure the periapsis-direction
    rotation between passages, without needing to replicate orbit.py's
    i/Omega/omega decomposition."""
    r = np.linalg.norm(r_vec)
    v2 = float(np.dot(v_vec, v_vec))
    rdotv = float(np.dot(r_vec, v_vec))
    return ((v2 - grav_param / r) * r_vec - rdotv * v_vec) / grav_param


def measure_delta_omega(sol, grav_param):
    """Signed apsidal shift (rad) between the first two detected periapsis
    events, about the orbit's angular-momentum axis."""
    t_events = sol.t_events[0]
    y_events = sol.y_events[0]
    if len(t_events) < 2:
        raise RuntimeError(f"Only {len(t_events)} periapsis passage(s) detected; "
                            "need >= 2 to measure a shift.")
    r1, v1 = y_events[0][:3], y_events[0][3:]
    r2, v2 = y_events[1][:3], y_events[1][3:]
    e1 = eccentricity_vector(r1, v1, grav_param)
    e2 = eccentricity_vector(r2, v2, grav_param)
    h1 = np.cross(r1, v1)
    h_hat = h1 / np.linalg.norm(h1)
    cross = np.cross(e1, e2)
    delta_omega = np.arctan2(np.dot(cross, h_hat), np.dot(e1, e2))
    return delta_omega, (t_events[1] - t_events[0])


def s301_initial_state():
    """S301's real orbit (constants.TRUTH) as a 3D Cartesian state vector,
    in the same sky/line-of-sight frame orbit.py rotates into (see
    orbit._rotate_to_sky) -- a fixed rotation of an inertial frame, so
    Newtonian/PN dynamics apply identically in it.

    Initialized at APOAPSIS (ecc_anom=pi), not periapsis: starting a
    Cartesian ODE integrator exactly at periapsis of an e=0.9832 orbit is
    numerically unstable -- the acceleration is at its sharpest there,
    and even a tight-tolerance adaptive integrator (DOP853) can take a
    bad first step and pick up spurious energy before its error control
    engages (verified directly: starting at periapsis produced apoapsis
    distances ~4x too large and total energy drift of order unity within
    the first 10% of an orbit). Starting at apoapsis, the smoothest point
    of the orbit, and integrating through the intervening periapsis
    passages instead keeps energy conserved to ~1e-5 over 2+ orbits.

    Returns (r_vec [m], v_vec [m/s], period [s], sma [m])."""
    period = k.TRUTH["P_yr"] * k.year
    sma = orbit.semi_major_axis_from_period(k.GM_BH, period)
    ecc = k.TRUTH["e"]
    # pylint: disable=protected-access
    # orbit.py has no public API for "plane state at a given eccentric
    # anomaly, pre-rotation" -- orbit_state()/orbit_state_precessing() only
    # return post-rotation sky-frame quantities. These two helpers are the
    # only way to get a genuine 3D Cartesian state vector for seeding the
    # ODE integration below.
    _r, _ta, x_p, y_p, vx_p, vy_p = orbit._plane_position_velocity(
        np.array([np.pi]), ecc, sma, k.GM_BH)
    rot = orbit._rotation_matrix(k.TRUTH["i_deg"], k.TRUTH["Omega_deg"], k.TRUTH["omega_deg"])
    # pylint: enable=protected-access
    pos = rot @ np.array([x_p[0], y_p[0]])
    vel = rot @ np.array([vx_p[0], vy_p[0]])
    return pos, vel, period, sma


# pylint: disable=too-many-locals
# main() orchestrates the whole point-mass-vs-extended-mass comparison and
# prints a running diagnostic narrative; splitting it up would scatter that
# narrative across several functions for no real gain in clarity.
def main():
    """Run the point-mass validation, then the 1-sigma and 3-sigma
    extended-mass cases, printing the full comparison narrative."""
    r0, v0, period, sma = s301_initial_state()
    ecc = k.TRUTH["e"]

    print("=== Point-mass-only validation (numerical integrator vs. analytic 1PN) ===")

    def m_point(r):
        return np.full_like(np.atleast_1d(r), k.M_BH, dtype=float)

    def phi_point(r):
        return -k.GM_BH / np.atleast_1d(np.asarray(r, dtype=float))

    accel_point = make_accel_fn(m_point, phi_point)
    sol_point = integrate_orbits(accel_point, (r0, v0), period, n_orbits=2)
    dphi_point, _dt_point = measure_delta_omega(sol_point, k.GM_BH)
    dphi_point_deg = np.degrees(dphi_point)

    analytic_rate = orbit.schwarzschild_precession_rate(k.GM_BH, sma, ecc, period)
    analytic_deg_per_orbit = np.degrees(analytic_rate) * period
    print(f"  numerical Delta-omega (1 orbit, point mass):  {dphi_point_deg:.6f} deg")
    print(f"  analytic 1PN prediction (orbit.py formula):    {analytic_deg_per_orbit:.6f} deg")
    print(f"  fractional agreement: {(dphi_point_deg / analytic_deg_per_orbit - 1) * 100:+.3f}%")

    sigma_fit = results_io.read_results("fit_orbit")["sigma_omega_dot_deg_yr"]
    pub_omega_dot = fit_orbit.TRUE_OMEGA_DOT_DEG_YR  # deg/yr, the 1PN rate the campaign injects
    results = {
        "pointmass_agreement_pct": ((dphi_point_deg / analytic_deg_per_orbit - 1) * 100, ".2f"),
        "rc_mpc": (RC / k.pc * 1e3, ".2f"),
        "rp_mpc": (sma * (1 - ecc) / k.pc * 1e3, ".3f"),
    }
    for tag, label, m_ext in (("one", "1-sigma", M_EXT_1SIGMA),
                              ("three", "3-sigma (conservative)", M_EXT_3SIGMA)):
        print()
        print(f"=== Extended-mass case, {label} bound: RE01 Plummer-like, alpha={ALPHA}, "
              f"rc={RC / k.pc * 1e3:.2f} mpc, M_ext={m_ext / k.Msun:.0f} Msun ===")
        m_total = k.M_BH + m_ext
        lambda_p = k.M_BH / m_total
        print(f"  lambda_p (point-mass fraction) = {lambda_p:.6f}, "
              f"lambda_e = {1 - lambda_p:.3e}")

        m_ext_fn = build_enclosed_mass_fn(m_total, lambda_p)
        phi_ext_fn = build_potential_fn(m_ext_fn, m_total)
        accel_ext = make_accel_fn(m_ext_fn, phi_ext_fn)
        sol_ext = integrate_orbits(accel_ext, (r0, v0), period, n_orbits=2)
        dphi_ext, _dt_ext = measure_delta_omega(sol_ext, k.GM_BH)
        dphi_ext_deg = np.degrees(dphi_ext)
        print(f"  numerical Delta-omega (1 orbit, point + extended): {dphi_ext_deg:.6f} deg")

        bias_deg_per_orbit = dphi_ext_deg - dphi_point_deg
        bias_rate_deg_yr = bias_deg_per_orbit / (period / k.year)
        print(f"  bias per orbit:      {bias_deg_per_orbit:+.6f} deg  "
              f"({bias_deg_per_orbit / analytic_deg_per_orbit * 100:+.3f}% of the 1PN signal)")
        print(f"  bias rate:           {bias_rate_deg_yr:+.6f} deg/yr")
        print(f"  vs. fitted omega_dot precision ({sigma_fit:.4f} deg/yr, Table 3): "
              f"{abs(bias_rate_deg_yr) / sigma_fit:.2f} sigma-equivalent")
        print(f"  vs. injected 1PN omega_dot ({pub_omega_dot:.4f} deg/yr): "
              f"{abs(bias_rate_deg_yr) / pub_omega_dot * 100:.3f}% fractional bias")
        results.update({
            f"bias_rate_{tag}": (bias_rate_deg_yr, ".5f"),
            f"bias_nsigma_{tag}": (abs(bias_rate_deg_yr) / sigma_fit, ".2f"),
            f"bias_pct_{tag}": (abs(bias_rate_deg_yr) / pub_omega_dot * 100, ".3f"),
            f"m_ext_{tag}_msun": (m_ext / k.Msun, ".0f"),
        })
    results_io.write_results("extended_mass_error", results)


if __name__ == "__main__":
    main()
