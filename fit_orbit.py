"""
Recover S301's orbital elements from the synthetic noisy GRAVITY+/VLTI
campaign (output/synthetic_observations.csv, produced by campaign.py:
astrometry + photometry, no radial-velocity channel -- see campaign.py's
module docstring for why RV isn't proposed), and compare the recovered
values against both the injected truth and the real published solution
(they're the same thing here: constants.TRUTH IS the paper's Solution 1)
-- this checks that the fitting pipeline itself is sound.

Method: nonlinear least squares (Levenberg-Marquardt via
scipy.optimize.least_squares) on 7 parameters -- the 6 Keplerian elements
(P, e, i, Omega, omega, t_peri) plus omega_dot, the apsidal precession
rate -- fit to the astrometric (RA, Dec) channel across BOTH periapsis
passages in campaign.py's two-passage dataset. The semi-major axis is
NOT an independent fit parameter -- it's derived from (GM_BH, P) via
Kepler's third law (orbit.semi_major_axis_from_period), consistent with
how it's treated everywhere else in this project. omega_deg is the
argument of periapsis AT the fitted t_peri epoch
(orbit.orbit_state_precessing's t_ref); this keeps the model
self-contained -- it never references constants.TRUTH, only its own 7
free parameters.

omega_dot's recovered value is compared against
orbit.schwarzschild_precession_rate's real analytic prediction (~0.23
deg/yr for S301) -- a genuine physics check, not a fit against an
injected/illustrative number: this is the actual real-GR value the
campaign's truth was generated with.

Photometry (dmag_K) is simulated by campaign.py but deliberately not
used here: it's driven by the same (r, v) the astrometric channel
already constrains, so it adds no independent orbital-element
information for this exercise.

Parameter uncertainties are estimated via bootstrap resampling of the
observation epochs (with replacement), refitting each resample --
matching this project's own established bootstrap-uncertainty convention
used elsewhere in this coursework.

VERIFIED, not assumed: an earlier version of this project included a
simulated ERIS RV channel and confirmed directly (by refitting the same
data with and without it) that dropping RV costs nothing -- every
parameter's bootstrap sigma was statistically unchanged (some even
marginally tighter without it). Astrometry alone constrains the full
orbit here, including omega_dot, which is why RV isn't part of this
campaign at all (a decision made on feasibility grounds, not a
limitation discovered after the fact -- see campaign.py).
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

# pylint: disable=wrong-import-position
# matplotlib's backend must be set (above) before pyplot is imported.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from scipy.optimize import least_squares

import campaign
import constants as k
import orbit
# pylint: enable=wrong-import-position

N_BOOTSTRAP = 1000
RNG_SEED = 301

_PERIOD_SEC = k.TRUTH["P_yr"] * k.year
_SMA_M = orbit.semi_major_axis_from_period(k.GM_BH, _PERIOD_SEC)
TRUE_OMEGA_DOT_DEG_YR = np.degrees(
    orbit.schwarzschild_precession_rate(k.GM_BH, _SMA_M, k.TRUTH["e"], _PERIOD_SEC)) * k.year

PARAM_NAMES = ["P_yr", "e", "i_deg", "Omega_deg", "omega_deg", "t_peri_yr", "omega_dot_deg_yr"]
TRUTH_VECTOR = np.array([
    k.TRUTH["P_yr"], k.TRUTH["e"], k.TRUTH["i_deg"],
    k.TRUTH["Omega_deg"], k.TRUTH["omega_deg"], k.TRUTH["t_peri_yr"],
    TRUE_OMEGA_DOT_DEG_YR,
])
# omega_dot has no published-paper sigma to compare against (not reported
# in the press-level sources this project could access) -- NaN prints as
# "n/a" in the comparison table.
PUBLISHED_SIGMA = np.array([
    k.TRUTH["sigma_P_yr"], k.TRUTH["sigma_e"], k.TRUTH["sigma_i_deg"],
    k.TRUTH["sigma_Omega_deg"], k.TRUTH["sigma_omega_deg"], k.TRUTH["sigma_t_peri_yr"],
    np.nan,
])
BOUNDS_LO = np.array([5.0, 0.5, 0.0, -180.0, -180.0, 2015.0, -5.0])
BOUNDS_HI = np.array([15.0, 0.999, 180.0, 540.0, 540.0, 2035.0, 5.0])
INITIAL_PERTURBATION = np.array(
    [0.03 * TRUTH_VECTOR[0], -0.03, 10.0, -15.0, 12.0, 0.05, 0.1])


def load_observations(path="output/synthetic_observations.csv"):
    """Load campaign.py's synthetic observations CSV as a structured array."""
    return np.genfromtxt(path, delimiter=",", names=True)


def model_observables(params, epochs_yr):
    """RA offset (mas), Dec offset (mas), and RV (km/s) predicted at the
    given epochs for the given 7-element parameter vector (6 Keplerian
    elements plus the apsidal precession rate omega_dot). RV is still
    computed and returned -- it's part of the physical orbit model and
    optimal_design.py's Fisher-information calculation uses it purely as
    an analytic cross-check that its own information contribution is
    negligible -- but it is NOT used in residuals()/fit_once() below,
    since campaign.py's data has no RV column at all (no instrument time
    was ever proposed for it -- see campaign.py's module docstring)."""
    p_yr, ecc, i_deg, raan_deg, omega_deg, t_peri_yr, omega_dot_deg_yr = params
    period = p_yr * k.year
    t_peri = t_peri_yr * k.year
    sma = orbit.semi_major_axis_from_period(k.GM_BH, period)
    elements = orbit.OrbitalElements(
        t_peri=t_peri, period=period, ecc=ecc, sma=sma,
        i_deg=i_deg, raan_deg=raan_deg, omega_deg=omega_deg,
    )
    omega_dot = np.radians(omega_dot_deg_yr) / k.year

    t_sec = epochs_yr * k.year
    state = orbit.orbit_state_precessing(t_sec, elements, omega_dot, t_ref=t_peri,
                                          grav_param=k.GM_BH)

    ra_mas, dec_mas = orbit.sky_offset_mas(state, k.D_OBS)
    rv_kms = state["v_los"] / 1e3
    return ra_mas, dec_mas, rv_kms


def residuals(params, epochs, ra_obs, dec_obs, sigma_ra, sigma_dec):
    """Sigma-normalized residuals (model - data) stacked across the
    astrometric (RA, Dec) channels. No RV term -- see module docstring."""
    ra_model, dec_model, _ = model_observables(params, epochs)
    return np.concatenate([
        (ra_model - ra_obs) / sigma_ra,
        (dec_model - dec_obs) / sigma_dec,
    ])


def fit_once(epochs, ra_obs, dec_obs, sigma_ra, sigma_dec, x0):
    """One Levenberg-Marquardt fit of the 6 orbital elements to the given
    (possibly resampled) data, starting from x0. Returns the best-fit
    parameter vector."""
    result = least_squares(
        residuals, x0=x0, bounds=(BOUNDS_LO, BOUNDS_HI),
        args=(epochs, ra_obs, dec_obs, sigma_ra, sigma_dec),
    )
    return result.x


def bootstrap_uncertainty(data, best_fit, n_resamples=N_BOOTSTRAP):
    """Parameter uncertainty via bootstrap resampling of the observation
    epochs (with replacement), refitting from best_fit each time."""
    rng = np.random.default_rng(RNG_SEED)
    n = len(data["epoch_yr"])
    samples = np.zeros((n_resamples, len(best_fit)))
    for k_iter in range(n_resamples):
        idx = rng.integers(0, n, size=n)
        samples[k_iter] = fit_once(
            data["epoch_yr"][idx], data["ra_offset_mas"][idx], data["dec_offset_mas"][idx],
            data["sigma_ra_mas"][idx], data["sigma_dec_mas"][idx], x0=best_fit,
        )
    return samples.std(axis=0)


def print_comparison(best_fit, fit_sigma):
    """Print truth vs. fitted vs. published-uncertainty for each parameter."""
    print(f"\n{'parameter':<18}{'truth':>12}{'fitted':>12}{'fit +/- sigma':>16}"
          f"{'published sigma':>18}")
    for name, truth, fit_val, fit_sig, pub_sig in zip(
            PARAM_NAMES, TRUTH_VECTOR, best_fit, fit_sigma, PUBLISHED_SIGMA):
        pub_sig_str = "n/a (unpublished)" if np.isnan(pub_sig) else f"{pub_sig:.4f}"
        print(f"{name:<18}{truth:>12.4f}{fit_val:>12.4f}{fit_sig:>16.4f}{pub_sig_str:>18}")


def build_plot_epoch_grid(epoch_min, epoch_max, period_yr, t_peri_yr, n_coarse=2000, n_dense=2000):
    """Time grid for smoothly plotting the fitted curve: a coarse
    background grid (fine for the slow apoapsis-side motion) UNION with
    locally dense windows around every periapsis passage the fitted
    ephemeris implies within [epoch_min, epoch_max]. A single
    evenly-spaced-in-time grid badly under-resolves the sharp turn near
    periapsis (most of an ~8.7-year period's worth of angular motion
    happens in a window of days, not years) -- without this, the plotted
    curve visibly cuts corners right where the orbit curves fastest, even
    though the underlying fit (evaluated only at the real data epochs) is
    correct."""
    coarse = np.linspace(epoch_min, epoch_max, n_coarse)
    first_k = int(np.floor((epoch_min - t_peri_yr) / period_yr)) - 1
    last_k = int(np.ceil((epoch_max - t_peri_yr) / period_yr)) + 1
    dense_windows = [coarse]
    for k_orbit in range(first_k, last_k + 1):
        passage_yr = t_peri_yr + k_orbit * period_yr
        if epoch_min - 0.1 <= passage_yr <= epoch_max + 0.1:
            window = np.linspace(passage_yr - 0.1, passage_yr + 0.1, n_dense)
            dense_windows.append(window[(window >= epoch_min) & (window <= epoch_max)])
    return np.sort(np.concatenate(dense_windows))


def _plot_sky_track(fig, ax, data, dense_epochs, ra_fit, dec_fit):
    """Sky-plane track panel. The orbit precesses and is sampled over
    ~1.5 periods, so the 2D track legitimately crosses itself near
    apoapsis -- each precessing loop's far side crosses the other's, real
    orbital geometry rather than a rendering error. A flat single-color
    line still makes those crossings look like stray artifact lines, so
    the path itself is colored by epoch (matching the data/model
    scatter), not just the sparse epoch-colored markers -- that way a
    crossing reads immediately as "two different times," not a glitch."""
    points = np.column_stack([ra_fit, dec_fit]).reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    path_lc = LineCollection(segments, cmap="viridis", norm=Normalize(dense_epochs.min(), dense_epochs.max()),
                              linewidths=0.9, alpha=0.6, zorder=1)
    path_lc.set_array(dense_epochs[:-1])
    ax.add_collection(path_lc)
    ax.plot([], [], color="#440154", lw=1.5, label="fitted orbit (path, colored by epoch)")
    model_scatter = ax.scatter(ra_fit[::15], dec_fit[::15], c=dense_epochs[::15], cmap="viridis",
                                s=8, marker="x", zorder=2, label="fitted orbit (epoch-colored)")
    ax.errorbar(data["ra_offset_mas"], data["dec_offset_mas"],
                xerr=data["sigma_ra_mas"], yerr=data["sigma_dec_mas"],
                fmt="none", ecolor="black", elinewidth=1.0, capsize=2, zorder=2)
    ax.scatter(data["ra_offset_mas"], data["dec_offset_mas"], c=data["epoch_yr"], cmap="viridis",
               s=20, edgecolors="black", linewidths=0.4, zorder=3,
               label="synthetic GRAVITY+ data")
    ax.plot(0, 0, "k*", ms=12, label="Sgr A*", zorder=4)
    ax.set_xlabel("RA offset (mas)")
    ax.set_ylabel("Dec offset (mas)")
    ax.set_title("S301 sky-plane track (GRAVITY+/VLTI astrometry): synthetic data vs. fitted orbit")
    # GRAVITY+'s real ~100 uas precision is ~1500x smaller than this
    # panel's ~150 mas span -- the error bars above are real and drawn at
    # true scale, but will look like a hairline or vanish entirely next
    # to the markers. Say so explicitly rather than leaving it looking
    # like the error bars are simply missing.
    ax.annotate(f"position error bars: +/-{data['sigma_ra_mas'][0] * 1000:.0f} uas per point\n"
                f"(too small to see at this plot's ~150 mas scale)",
                xy=(0.02, 0.02), xycoords="axes fraction", fontsize=7, color="dimgray")
    ax.invert_xaxis()
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.colorbar(model_scatter, ax=ax, label="epoch (year)")


def _plot_photometry(ax, data, dense_epochs, dmag_true_dense):
    """Photometry panel: synthetic data (with error bars) vs. the
    injected truth. Not fit -- see module docstring."""
    ax.errorbar(data["epoch_yr"], data["dmag_K"], yerr=data["sigma_dmag"],
                fmt="o", ms=3, color="#1f77b4", ecolor="#1f77b466",
                label="synthetic GRAVITY+ data")
    ax.plot(dense_epochs, dmag_true_dense, color="gray", lw=1.0,
            label="injected truth (not fit)")
    ax.invert_yaxis()
    ax.set_xlabel("Epoch (year)")
    ax.set_ylabel("Delta magnitude (fainter down)")
    ax.set_title("K-band photometry (GRAVITY+/VLTI): synthetic data (not used in the fit)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)


# pylint: disable=duplicate-code
# The tight_layout/makedirs/savefig/print tail is standard matplotlib
# boilerplate shared with s301_lightcurve.py's make_plots -- there's no
# shared "plotting utilities" module in this project worth introducing
# for four lines of save-figure calls, and orbit.py (the one genuinely
# shared module) is deliberately physics-only.
def make_plots(data, best_fit, outpath="output/fit_orbit.png"):
    """Sky-plane track and photometry: synthetic data vs. the fitted
    orbit. No RV panels -- campaign.py's data has no RV channel (see its
    module docstring for why)."""
    dense_epochs = build_plot_epoch_grid(
        data["epoch_yr"].min(), data["epoch_yr"].max(), best_fit[0], best_fit[5])
    ra_fit, dec_fit, _ = model_observables(best_fit, dense_epochs)
    _, _, _, dmag_true_dense = campaign.true_observables(dense_epochs)

    fig, axes = plt.subplots(2, 1, figsize=(9, 8))
    _plot_sky_track(fig, axes[0], data, dense_epochs, ra_fit, dec_fit)
    _plot_photometry(axes[1], data, dense_epochs, dmag_true_dense)

    plt.tight_layout()
    os.makedirs("output", exist_ok=True)
    plt.savefig(outpath, dpi=150)
    print(f"Saved figure to {outpath}")
# pylint: enable=duplicate-code


def main():
    """Load the campaign, fit the orbit, bootstrap uncertainties, and
    report/plot the result."""
    data = load_observations()
    x0 = TRUTH_VECTOR + INITIAL_PERTURBATION
    print(f"Initial guess (perturbed from truth): {dict(zip(PARAM_NAMES, x0))}")

    best_fit = fit_once(
        data["epoch_yr"], data["ra_offset_mas"], data["dec_offset_mas"],
        data["sigma_ra_mas"], data["sigma_dec_mas"], x0=x0,
    )
    print(f"Running {N_BOOTSTRAP} bootstrap resamples for parameter uncertainties...")
    fit_sigma = bootstrap_uncertainty(data, best_fit)

    print_comparison(best_fit, fit_sigma)
    make_plots(data, best_fit)


if __name__ == "__main__":
    main()
