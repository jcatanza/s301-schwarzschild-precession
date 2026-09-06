"""Kepler-equation solver: global convergence at high eccentricity.

orbit.solve_kepler is a bisection-safeguarded Newton-Raphson. Its
docstring documents the case that made plain Newton-Raphson diverge
(S301's own e=0.9832 at M=0.2185 rad, blowing up to E~1e13); this test
pins that case and sweeps a fixed-seed random sample of (M, e) pairs up
to e=0.9999, well beyond anything the fit explores.
"""

import numpy as np

import orbit


def kepler_residual(ecc_anom, ecc, mean_anom):
    """|E - e sin E - M|, the quantity solve_kepler drives to zero."""
    return np.abs(ecc_anom - ecc * np.sin(ecc_anom) - mean_anom)


def test_documented_divergent_case_converges():
    """The (M, e) pair that broke the unsafeguarded solver."""
    mean_anom, ecc = 0.2185, 0.9832
    ecc_anom = orbit.solve_kepler(np.array([mean_anom]), ecc)
    assert np.isfinite(ecc_anom).all()
    assert kepler_residual(ecc_anom, ecc, mean_anom).max() < 1e-10
    # The true root is near E = 1.1, not 1e13.
    assert 1.0 < ecc_anom[0] < 1.2


def test_random_sweep_worst_residual():
    """2000 fixed-seed (M, e) pairs, e up to 0.9999, M anywhere on the
    real line (the solver wraps internally)."""
    rng = np.random.default_rng(301)
    n_pairs = 2000
    mean_anoms = rng.uniform(-4 * np.pi, 4 * np.pi, size=n_pairs)
    eccs = rng.uniform(0.0, 0.9999, size=n_pairs)
    worst = 0.0
    for mean_anom, ecc in zip(mean_anoms, eccs):
        ecc_anom = orbit.solve_kepler(np.array([mean_anom]), ecc)
        assert np.isfinite(ecc_anom).all()
        worst = max(worst, kepler_residual(ecc_anom, ecc, mean_anom).max())
    assert worst < 1e-9


def test_vectorized_matches_scalar_calls():
    """Passing an array of M gives the same answer as one call per M."""
    ecc = 0.9832
    mean_anoms = np.linspace(-np.pi, np.pi, 257)
    vector = orbit.solve_kepler(mean_anoms, ecc)
    scalar = np.array([orbit.solve_kepler(np.array([m]), ecc)[0] for m in mean_anoms])
    assert np.allclose(vector, scalar, rtol=0, atol=1e-11)


def test_windings_preserved():
    """M and M + 2*pi*n give E and E + 2*pi*n (the solver unwraps)."""
    ecc = 0.9
    base = orbit.solve_kepler(np.array([0.7]), ecc)[0]
    for n_wind in (-3, 1, 5):
        shifted = orbit.solve_kepler(np.array([0.7 + 2 * np.pi * n_wind]), ecc)[0]
        assert abs(shifted - base - 2 * np.pi * n_wind) < 1e-10
