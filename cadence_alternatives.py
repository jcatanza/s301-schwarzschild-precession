"""
Cadence-alternative tests for the adopted campaign design (Design D,
campaign.build_epoch_grid): every "what if we had scheduled it
differently" number the manuscript quotes is produced here, on the same
noise model, fit, and bootstrap machinery as the headline result, and
written to results/cadence_alternatives.json.

Design D has two kinds of epoch: the periapsis floor (within +/-20 d of
a season-adjusted periapsis anchor -- optimal_design.near_periapsis_mask)
and the "flexible middle" -- the epochs the Fisher-optimal search placed
between the two passages (2030 < t < 2039.5, outside the dense windows).
The tests below move or replace one of those groups at a time, keeping
the total epoch count fixed unless stated otherwise:

  reallocate2032  Move the flexible-middle epochs into the single visible
                  window right after passage 1's solar conjunction
                  (2032.09-2032.74). Tests whether "observe as soon as the
                  star re-emerges" would beat the optimizer's spread.
  complete_p1     Hypothetical: pretend passage 1 were fully observable
                  and centre its dense cluster on the true periapsis
                  (k.NEXT_PERIAPSIS_YR +/- 20 d), ignoring the visibility
                  season. Quantifies what the conjunction costs.
  doubled         Design D rebuilt at twice the epoch budget.
  placement scan  Move the flexible middle into a 0.6-yr visible window
                  at five trial centres between the passages, to show how
                  flat (or not) sigma(omega_dot) is in that placement.
  monte carlo     30 independent noise realizations of Design D, each fit
                  once; the scatter of the best fits across realizations
                  is an independent check on the bootstrap sigmas that
                  fit_orbit.py reports (results/fit_orbit.json).

Baseline sigmas are read from results/fit_orbit.json (1000 bootstrap
resamples of the actual campaign data), never hardcoded.
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np  # noqa: E402  pylint: disable=wrong-import-position

import campaign  # noqa: E402  pylint: disable=wrong-import-position
import constants as k  # noqa: E402  pylint: disable=wrong-import-position
import fit_orbit  # noqa: E402  pylint: disable=wrong-import-position
import optimal_design  # noqa: E402  pylint: disable=wrong-import-position
import results_io  # noqa: E402  pylint: disable=wrong-import-position

N_BOOT = 400
N_BOOT_SCAN = 300
N_MC_REALIZATIONS = 30
MC_SEED_BASE = 1000
FLEXIBLE_START_YR = 2030.0
FLEXIBLE_END_YR = 2039.5
REALLOCATE_WINDOW_YR = (2032.09, 2032.74)   # first visible window after passage 1's conjunction
SCAN_CENTERS_YR = (2032.3, 2033.0, 2034.5, 2036.0, 2037.5)
SCAN_HALF_WIDTH_YR = 0.3
SCAN_TAGS = "abcde"
OMEGA_DOT_INDEX = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")
SHORT_NAME = {"omega_dot_deg_yr": "omega_dot"}


def simulate(epochs, seed=campaign.RNG_SEED):
    """Synthetic GRAVITY+ astrometry at the given epochs: campaign.py's
    noise-free truth plus N(0, GRAVITY_PLUS_ASTROMETRY_MAS) noise, in the
    dict layout fit_orbit.bootstrap_samples consumes. Same RNG seed and
    draw order as campaign.report_and_write_csv, so on Design D's own
    epochs this reproduces the headline dataset."""
    ra_true, dec_true, _, _ = campaign.true_observables(epochs)
    rng = np.random.default_rng(seed)
    sigma = k.GRAVITY_PLUS_ASTROMETRY_MAS
    return {
        "epoch_yr": epochs,
        "ra_offset_mas": ra_true + rng.normal(0, sigma, size=epochs.shape),
        "dec_offset_mas": dec_true + rng.normal(0, sigma, size=epochs.shape),
        "sigma_ra_mas": np.full_like(epochs, sigma),
        "sigma_dec_mas": np.full_like(epochs, sigma),
    }


def fit_data(data):
    """fit_orbit.fit_once from the standard perturbed-truth warm start."""
    return fit_orbit.fit_once(
        data["epoch_yr"], data["ra_offset_mas"], data["dec_offset_mas"],
        data["sigma_ra_mas"], data["sigma_dec_mas"],
        x0=fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION,
    )


def fit_and_bootstrap(epochs, label, n_boot=N_BOOT):
    """Simulate, fit, and bootstrap one cadence; print and return the
    per-parameter bootstrap sigma vector."""
    data = simulate(epochs)
    best = fit_data(data)
    print(f"[{label}] {len(epochs)} epochs, {epochs.min():.2f}-{epochs.max():.2f}; "
          f"running {n_boot} bootstrap resamples...")
    sigma = fit_orbit.bootstrap_samples(data, best, n_resamples=n_boot).std(axis=0)
    print(f"[{label}] sigma(omega_dot) = {sigma[OMEGA_DOT_INDEX]:.5f} deg/yr; "
          f"omega_dot bias = {best[OMEGA_DOT_INDEX] - fit_orbit.TRUE_OMEGA_DOT_DEG_YR:+.5f}")
    return sigma


def flexible_middle_mask(epochs):
    """The epochs the Fisher-optimal search placed between the two
    passages: outside both +/-20-d dense windows and strictly inside
    (FLEXIBLE_START_YR, FLEXIBLE_END_YR)."""
    return (~optimal_design.near_periapsis_mask(epochs)
            & (epochs > FLEXIBLE_START_YR) & (epochs < FLEXIBLE_END_YR))


def passage_one_mask(epochs):
    """Passage 1's dense cluster: within +/-20 d of the season-adjusted
    anchor AND closer to that anchor than to passage 2's."""
    anchors = optimal_design.periapsis_anchors()
    closer_to_one = np.abs(epochs - anchors[0]) < np.abs(epochs - anchors[1])
    return optimal_design.near_periapsis_mask(epochs) & closer_to_one


def visible_window_epochs(center_yr, half_width_yr, n_epochs):
    """n_epochs epochs spread evenly over the visibility-season part of
    [center - half_width, center + half_width]: a fine grid across the
    window is filtered to the Mar-Sep season, then thinned to n_epochs
    evenly spaced survivors."""
    fine = np.linspace(center_yr - half_width_yr, center_yr + half_width_yr, 4000)
    visible = fine[campaign._in_visibility_season(fine)]  # pylint: disable=protected-access
    picks = np.round(np.linspace(0, len(visible) - 1, n_epochs)).astype(int)
    return visible[picks]


def replace_epochs(epochs, drop_mask, new_epochs):
    """Epoch grid with the masked epochs removed and new_epochs added."""
    return np.sort(np.concatenate([epochs[~drop_mask], new_epochs]))


def test_reallocate_2032(base_epochs, base_sigma_omega_dot):
    """(i) Move the flexible middle into the 2032 post-conjunction window."""
    move = flexible_middle_mask(base_epochs)
    new = np.linspace(REALLOCATE_WINDOW_YR[0], REALLOCATE_WINDOW_YR[1], int(move.sum()))
    sigma = fit_and_bootstrap(replace_epochs(base_epochs, move, new), "reallocate2032")
    ratio = sigma[OMEGA_DOT_INDEX] / base_sigma_omega_dot
    print(f"[reallocate2032] moved {move.sum()} epochs; sigma(omega_dot) ratio to base = {ratio:.2f} "
          f"({(ratio - 1) * 100:+.0f}%)")
    return {
        "reallocate_n_moved": int(move.sum()),
        "reallocate_sigma_omega_dot": (sigma[OMEGA_DOT_INDEX], ".5f"),
        "reallocate_ratio": (ratio, ".2f"),
        "reallocate_pct_change": ((ratio - 1) * 100, ".0f"),
    }


def test_complete_passage_one(base_epochs, base_sigma):
    """(ii) Hypothetical: passage 1's cluster centred on the true
    periapsis, ignoring the visibility season."""
    drop = passage_one_mask(base_epochs)
    new = np.linspace(k.NEXT_PERIAPSIS_YR - campaign.DENSE_HALF_WIDTH_YR,
                      k.NEXT_PERIAPSIS_YR + campaign.DENSE_HALF_WIDTH_YR, int(drop.sum()))
    sigma = fit_and_bootstrap(replace_epochs(base_epochs, drop, new), "complete_p1")
    results = {"complete_n_replaced": int(drop.sum())}
    print(f"[complete_p1] replaced {drop.sum()} epochs; sigma ratios to base:")
    for name, sig, base in zip(fit_orbit.PARAM_NAMES, sigma, base_sigma):
        print(f"    {name:<18} {sig:.5f} / {base:.5f} = {sig / base:.2f}")
        results[f"complete_sigma_{SHORT_NAME.get(name, name)}"] = (sig, ".5f")
        results[f"complete_ratio_{SHORT_NAME.get(name, name)}"] = (sig / base, ".2f")
    return results


def test_doubled(n_base, base_sigma_omega_dot):
    """(iii) Design D rebuilt at twice the epoch budget."""
    print(f"[doubled] building the constrained design at n_total={2 * n_base} (slow)...")
    epochs, _ = optimal_design.build_epoch_grid_constrained(n_total=2 * n_base)
    sigma = fit_and_bootstrap(epochs, "doubled")
    ratio = sigma[OMEGA_DOT_INDEX] / base_sigma_omega_dot
    print(f"[doubled] sigma(omega_dot) ratio to base = {ratio:.2f} ({(1 - ratio) * 100:.0f}% gain; "
          f"1/sqrt(2) = {1 / np.sqrt(2):.2f})")
    return {
        "doubled_n": len(epochs),
        "doubled_sigma_omega_dot": (sigma[OMEGA_DOT_INDEX], ".5f"),
        "doubled_ratio": (ratio, ".2f"),
        "doubled_pct_gain": ((1 - ratio) * 100, ".0f"),
    }


def test_placement_scan(base_epochs):
    """(iv) Flexible middle moved into a 0.6-yr visible window at each
    trial centre."""
    move = flexible_middle_mask(base_epochs)
    results = {}
    for tag, center in zip(SCAN_TAGS, SCAN_CENTERS_YR):
        new = visible_window_epochs(center, SCAN_HALF_WIDTH_YR, int(move.sum()))
        sigma = fit_and_bootstrap(replace_epochs(base_epochs, move, new), f"scan_{tag} @ {center:.1f}",
                                  n_boot=N_BOOT_SCAN)
        results[f"scan_{tag}_center"] = (center, ".1f")
        results[f"scan_{tag}_sigma_omega_dot"] = (sigma[OMEGA_DOT_INDEX], ".5f")
    print("[scan] sigma(omega_dot) vs. window centre: " + ", ".join(
        f"{results[f'scan_{t}_center'][0]:.1f}: {results[f'scan_{t}_sigma_omega_dot'][0]:.5f}"
        for t in SCAN_TAGS))
    return results


def test_monte_carlo(base_epochs, base_sigma):
    """(v) Independent-realization cross-check of the bootstrap: fit
    N_MC_REALIZATIONS fresh noise draws of Design D and take the scatter
    of the best fits."""
    print(f"[mc] fitting {N_MC_REALIZATIONS} independent noise realizations of Design D...")
    fits = np.array([fit_data(simulate(base_epochs, seed=MC_SEED_BASE + i))
                     for i in range(N_MC_REALIZATIONS)])
    mc_sigma = fits.std(axis=0, ddof=1)
    results = {"mc_n_realizations": N_MC_REALIZATIONS}
    print(f"{'parameter':<18}{'MC sigma':>12}{'bootstrap':>12}{'MC/boot':>9}{'mean bias':>12}")
    for name, sig, base, bias in zip(fit_orbit.PARAM_NAMES, mc_sigma, base_sigma,
                                     fits.mean(axis=0) - fit_orbit.TRUTH_VECTOR):
        print(f"{name:<18}{sig:>12.5f}{base:>12.5f}{sig / base:>9.2f}{bias:>+12.5f}")
        results[f"mc_sigma_{SHORT_NAME.get(name, name)}"] = (sig, ".5f")
        results[f"mc_over_boot_{SHORT_NAME.get(name, name)}"] = (sig / base, ".2f")
    return results


def main():
    """Run every cadence-alternative test against Design D and write
    results/cadence_alternatives.json."""
    base_epochs = campaign.build_epoch_grid()
    headline = results_io.read_results("fit_orbit")
    base_sigma = np.array([headline[f"sigma_{name}"] for name in fit_orbit.PARAM_NAMES])
    base_sigma_omega_dot = base_sigma[OMEGA_DOT_INDEX]
    print(f"Base design D: {len(base_epochs)} epochs; "
          f"{flexible_middle_mask(base_epochs).sum()} flexible-middle, "
          f"{passage_one_mask(base_epochs).sum()} in passage 1's cluster; headline "
          f"sigma(omega_dot) = {base_sigma_omega_dot:.5f} deg/yr from results/fit_orbit.json")

    results = {"n_boot": N_BOOT, "base_sigma_omega_dot": (base_sigma_omega_dot, ".5f"),
               "n_epochs_base": len(base_epochs)}
    results.update(test_reallocate_2032(base_epochs, base_sigma_omega_dot))
    results.update(test_complete_passage_one(base_epochs, base_sigma))
    results.update(test_placement_scan(base_epochs))
    results.update(test_monte_carlo(base_epochs, base_sigma))
    results.update(test_doubled(len(base_epochs), base_sigma_omega_dot))
    path = results_io.write_results("cadence_alternatives", results)
    print(f"Wrote {len(results)} results to {path}")


if __name__ == "__main__":
    main()
