"""
Fisher-information (D-optimal) design for S301's observing campaign,
replacing campaign.py's "~90% of epochs within +/-20 days of periapsis"
heuristic with a mathematically principled alternative.

Motivation: the heuristic was a reasonable first cut at adapting the
cadence to S301's extreme eccentricity (e=0.9832), but it treats all 7
fit parameters (P, e, i, Omega, omega, t_peri, omega_dot) as equally well
served by periapsis-proximity. Running the actual fit on that design
showed this is false: P_yr and e got MORE precise than an evenly-spread
cadence, but i, Omega, omega, and omega_dot got LESS precise -- because
pinning down the orbit's 3D orientation needs angular diversity across
the orbit's arc, not just a tight cluster at one point on it.

Method: this is a standard problem in optimal experimental design
theory. For a fine grid of candidate epochs (filtered to the real Mar-Sep
Paranal visibility season already established in campaign.py), compute
the Jacobian d(RA,Dec)/d(theta) for all 7 parameters via numerical
central differences around the TRUTH values, using fit_orbit.py's own
model_observables -- the SAME model the actual fit uses, so the design
is honest about what it's optimizing for. (No RV term: campaign.py's
data has no RV channel at all -- see its module docstring for why -- and
RV's Fisher contribution was confirmed negligible before it was dropped,
so excluding it doesn't change any design decision, only the code's
honesty about what's actually collected.) Weight each channel by its
real instrument noise (GRAVITY+ astrometry) to get each
candidate epoch's Fisher-information contribution (a 7x7 matrix per
candidate epoch). Then greedily select N epochs maximizing the
determinant of the summed Fisher matrix -- a D-optimal design, which
minimizes the volume of the joint parameter uncertainty ellipsoid. This
is the standard default criterion in optimal design theory; A-optimal
(minimize sum of variances) or a single-parameter c-optimal design
(minimize omega_dot's variance alone, since that IS the science target)
are the usual alternatives, noted but not implemented here.

The by-product of this approach is verifiable *before* ever running a
fit: the Fisher matrix's inverse is the Cramer-Rao bound, an analytic
lower bound on achievable parameter variance, so designs can be compared
instantly before paying for the 1000-resample bootstrap that empirically
confirms it.

IMPORTANT, real negative result from the first attempt (naive_fisher_
per_epoch/build_epoch_grid_naive below, kept for the comparison): a
Fisher matrix linearized at a single fiducial truth point catastrophically
failed here. It predicted 800-28000x precision improvements by piling
nearly all epochs onto two calendar-adjacent clusters ~18 years from the
reference epoch (2 orbital periods away) and covering NEITHER periapsis
at all. Running the real nonlinear bootstrap fit on that design showed
it was actually 4-600x WORSE than the heuristic on every parameter. The
mechanism: sensitivity to P and t_peri genuinely grows with time
baseline (real physics -- it's why pulsar timing gets great period
precision from long baselines), but a small period error, extrapolated
across >2 periods, compounds into large nonlinear phase drift (a
"cycle-slip"/aliasing hazard). The linearized Fisher information is only
valid for small deviations from the exact truth; it can't see that this
design's information is fragile -- correct only if you already know P to
absurd precision -- while a real nonlinear fit starting from a realistic,
imperfect initial guess has no way to reach that narrow, precisely-tuned
optimum, and gets no orbit-shape information at all from the two
periapsis-free clusters.

The fix (build_epoch_grid_robust, the actual design used below): a
pseudo-Bayesian D-optimal design (Chaloner & Verdinelli 1995) -- instead
of linearizing at one exact truth point, average the per-epoch Fisher
matrix over many parameter draws from a realistic prior (this project's
own published element uncertainties, constants.TRUTH's sigma_* fields).
A candidate epoch whose apparent information is fragile -- large only
for one exact parameter vector, and wildly different for nearby,
equally-plausible vectors -- gets averaged down. Only epochs that stay
informative across the whole plausible parameter range remain
attractive, which directly penalizes reliance on precise, unverifiable
cycle-counting.
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np  # noqa: E402  pylint: disable=wrong-import-position

import campaign  # noqa: E402  pylint: disable=wrong-import-position
import constants as k  # noqa: E402  pylint: disable=wrong-import-position
import fit_orbit  # noqa: E402  pylint: disable=wrong-import-position

CANDIDATE_STEP_DAYS = 1.0
N_MC_DRAWS = 300
MC_RNG_SEED = 301
# Relative finite-difference step per parameter. t_peri_yr gets its own,
# much smaller step: near periapsis the model changes so fast with
# t_peri (a tiny shift in "when periapsis happens" moves the star through
# a large arc) that the default step would smear the numerical
# derivative across a physically significant fraction of the passage
# itself, rather than probing the true local slope.
FD_STEP_FRAC = 1e-6
FD_STEP_FRAC_T_PERI = 1e-7
RIDGE = 1e-6   # tiny numerical-stability prior, not real information


def observable_jacobian(param_vec, epochs_yr):
    """Numerical Jacobian d(RA,Dec)/d(theta) at every epoch in epochs_yr,
    evaluated at param_vec via central finite differences on
    fit_orbit.model_observables (the same model the real fit uses). No
    RV row: campaign.py's data has no RV channel at all (see its module
    docstring), and an earlier version of this project confirmed RV's
    Fisher-information contribution is ~8 orders of magnitude smaller
    than astrometry's anyway (sigma_RV=1621 km/s vs. sigma_astrom=0.1
    mas), so dropping it doesn't change any design decision -- only the
    code's honesty about what's actually being collected. Returns an
    array of shape (n_epochs, 2, n_params)."""
    n_params = len(param_vec)
    jac = np.zeros((len(epochs_yr), 2, n_params))
    for j, val in enumerate(param_vec):
        frac = FD_STEP_FRAC_T_PERI if fit_orbit.PARAM_NAMES[j] == "t_peri_yr" else FD_STEP_FRAC
        step = max(abs(val) * frac, 1e-8)
        p_hi = param_vec.copy()
        p_hi[j] = val + step
        p_lo = param_vec.copy()
        p_lo[j] = val - step
        ra_hi, dec_hi, _ = fit_orbit.model_observables(p_hi, epochs_yr)
        ra_lo, dec_lo, _ = fit_orbit.model_observables(p_lo, epochs_yr)
        jac[:, 0, j] = (ra_hi - ra_lo) / (2 * step)
        jac[:, 1, j] = (dec_hi - dec_lo) / (2 * step)
    return jac


def fisher_per_epoch(jac, sigma_ra, sigma_dec):
    """Per-epoch Fisher-information contribution, treating (RA, Dec) at a
    given epoch as two independent Gaussian-noise channels. Returns an
    array of shape (n_epochs, n_params, n_params)."""
    inv_sigma = np.array([1.0 / sigma_ra, 1.0 / sigma_dec])
    weighted = jac * inv_sigma[None, :, None]   # (n_epochs, 2, n_params)
    return np.einsum("ncp,ncq->npq", weighted, weighted)


def prior_sigma_vector():
    """1-sigma prior width per parameter for the robust design, reusing
    this project's own published element uncertainties
    (fit_orbit.PUBLISHED_SIGMA) wherever they exist. omega_dot has no
    published sigma (it isn't a measured quantity yet) -- fall back to
    half its own true (real 1PN formula) value: generous enough to
    reflect real prior ignorance without being absurd, and explicitly
    NOT an injected/illustrative number, same spirit as the rest of this
    project's honesty conventions."""
    sigma = fit_orbit.PUBLISHED_SIGMA.copy()
    sigma[-1] = 0.5 * fit_orbit.TRUE_OMEGA_DOT_DEG_YR
    return sigma


def robust_fisher_per_epoch(candidates, n_draws=N_MC_DRAWS, rng_seed=MC_RNG_SEED):
    """Pseudo-Bayesian (Chaloner & Verdinelli 1995) per-epoch Fisher
    information: average the naive per-point Fisher matrix over n_draws
    parameter vectors sampled from a realistic prior (prior_sigma_vector)
    around the truth, instead of linearizing at one exact fiducial point.
    This is the fix for naive_fisher_per_epoch's aliasing failure mode
    (see module docstring): a candidate epoch whose apparent information
    depends on knowing P/t_peri to unrealistic precision gets averaged
    down, because nearby-but-different draws see very different (often
    much smaller, sometimes inconsistent) sensitivity there."""
    rng = np.random.default_rng(rng_seed)
    truth = fit_orbit.TRUTH_VECTOR.copy()
    sigma_prior = prior_sigma_vector()
    draws = truth[None, :] + rng.normal(0, 1, size=(n_draws, len(truth))) * sigma_prior[None, :]
    draws[:, 1] = np.clip(draws[:, 1], 0.5, 0.999)   # keep eccentricity physical

    fisher_sum = np.zeros((len(candidates), len(truth), len(truth)))
    for draw in draws:
        jac = observable_jacobian(draw, candidates)
        fisher_sum += fisher_per_epoch(
            jac, k.GRAVITY_PLUS_ASTROMETRY_MAS, k.GRAVITY_PLUS_ASTROMETRY_MAS)
    return fisher_sum / n_draws


def greedy_d_optimal(fisher_cand, n_select, ridge=RIDGE, initial_fisher=None):
    """Greedily select n_select candidate indices maximizing the
    determinant of the cumulative Fisher information matrix (D-optimal
    design -- minimizes the volume of the joint parameter uncertainty
    ellipsoid). `ridge` is a tiny numerical-stability prior so the
    determinant stays well-defined before enough epochs are chosen to
    constrain all parameters -- not real information. `initial_fisher`,
    if given, seeds the cumulative matrix with a fixed contribution
    (e.g. a domain-knowledge periapsis floor -- see
    build_epoch_grid_constrained) that the greedy search then optimally
    supplements, rather than re-discovering from scratch.

    Exact D-optimal subset selection is combinatorially intractable
    (choosing n_select of n_candidates); the greedy sequential-addition
    algorithm used here is the standard, well-established approximation
    in optimal experimental design theory."""
    n_params = fisher_cand.shape[1]
    cumulative = ridge * np.eye(n_params)
    if initial_fisher is not None:
        cumulative = cumulative + initial_fisher
    remaining = list(range(len(fisher_cand)))
    selected = []
    for _ in range(n_select):
        candidate_matrices = cumulative[None, :, :] + fisher_cand[remaining]
        signs, logdets = np.linalg.slogdet(candidate_matrices)
        logdets = np.where(signs > 0, logdets, -np.inf)
        best_local = int(np.argmax(logdets))
        best_idx = remaining.pop(best_local)
        cumulative = cumulative + fisher_cand[best_idx]
        selected.append(best_idx)
    return np.array(selected), cumulative


def build_candidate_grid():
    """Fine, visibility-season-filtered candidate epoch grid spanning
    the same campaign window as campaign.py's heuristic design."""
    step_yr = CANDIDATE_STEP_DAYS / 365.25
    candidates = np.arange(campaign.CAMPAIGN_START_YR, campaign.CAMPAIGN_END_YR, step_yr)
    return candidates[campaign._in_visibility_season(candidates)]  # pylint: disable=protected-access


def build_epoch_grid_naive(n_total=None):
    """NAIVE D-optimal epoch grid (kept only for comparison -- see module
    docstring): greedily selects n_total epochs from a fine visibility-
    filtered candidate grid, maximizing the Fisher-information
    determinant linearized at a single fiducial TRUTH parameter vector.
    This is the version that catastrophically failed the real bootstrap
    test. Returns (epochs, cumulative_fisher_matrix)."""
    if n_total is None:
        n_total = len(campaign.build_epoch_grid_heuristic())

    candidates = build_candidate_grid()
    jac = observable_jacobian(fit_orbit.TRUTH_VECTOR.copy(), candidates)
    fisher_cand = fisher_per_epoch(jac, k.GRAVITY_PLUS_ASTROMETRY_MAS, k.GRAVITY_PLUS_ASTROMETRY_MAS)

    selected_idx, cumulative_fisher = greedy_d_optimal(fisher_cand, n_total)
    return np.sort(candidates[selected_idx]), cumulative_fisher


def build_epoch_grid_robust(n_total=None, n_draws=N_MC_DRAWS):
    """Pseudo-Bayesian D-optimal epoch grid: greedily selects n_total
    epochs from a fine visibility-filtered candidate grid, maximizing the
    determinant of the Fisher matrix AVERAGED over a realistic prior on
    the parameters (robust_fisher_per_epoch), not linearized at one exact
    truth point. This is the fix for build_epoch_grid_naive's aliasing
    failure mode. Returns (epochs, cumulative_fisher_matrix)."""
    if n_total is None:
        n_total = len(campaign.build_epoch_grid_heuristic())

    candidates = build_candidate_grid()
    fisher_cand = robust_fisher_per_epoch(candidates, n_draws=n_draws)

    selected_idx, cumulative_fisher = greedy_d_optimal(fisher_cand, n_total)
    return np.sort(candidates[selected_idx]), cumulative_fisher


PERIAPSIS_FLOOR_FRACTION = 0.4
# Fraction of the epoch budget reserved, up front, for dense coverage
# right at each periapsis passage -- a hard domain-knowledge constraint,
# not something left to the statistical optimality search. Reason:
# Fisher information (even pseudo-Bayesian-averaged over parameter
# uncertainty, build_epoch_grid_robust) only measures NOISE-driven
# uncertainty around a presumed-correct optimum -- it has no way to
# penalize STRUCTURAL/GEOMETRIC non-identifiability from sampling too
# few genuinely distinct orbital phases. Classical orbit determination
# (Gauss's method and its descendants) has always required observations
# spread across an orbit, not just high precision at a couple of points,
# because a highly eccentric orbit's shape is resolved by its curvature
# -- which periapsis alone provides. Both build_epoch_grid_naive and
# build_epoch_grid_robust tried to skip periapsis entirely and failed
# the real bootstrap test as a result (see module docstring). This floor
# encodes that known requirement directly; the Fisher-optimal search is
# then used only for what it's legitimately good at: allocating the
# REMAINING flexible budget to best constrain the angular-orientation
# elements, exactly where the original uniform-allocation heuristic was
# shown to underperform.


def _periapsis_floor_epochs(n_floor):
    """n_floor epochs split evenly between the two periapsis passages'
    +/-20-day dense windows (season-adjusted anchor for passage 1, same
    construction as campaign.py's heuristic design)."""
    # pylint: disable=protected-access
    anchors = [campaign._season_anchor(k.NEXT_PERIAPSIS_YR),
               campaign._season_anchor(k.NEXT_PERIAPSIS_YR + k.TRUTH["P_yr"])]
    n_per_passage = n_floor // 2
    windows = [np.linspace(a - campaign.DENSE_HALF_WIDTH_YR, a + campaign.DENSE_HALF_WIDTH_YR,
                            n_per_passage) for a in anchors]
    epochs = np.concatenate(windows)
    return epochs[campaign._in_visibility_season(epochs)]
    # pylint: enable=protected-access


def build_epoch_grid_constrained(n_total=None, floor_fraction=PERIAPSIS_FLOOR_FRACTION,
                                  n_draws=N_MC_DRAWS):
    """Constrained pseudo-Bayesian D-optimal design: a fixed periapsis
    floor (domain knowledge) plus Fisher-optimal allocation of the
    remaining flexible budget (robust_fisher_per_epoch) across the full
    visibility-filtered candidate grid. Fixes build_epoch_grid_robust's
    remaining failure mode: pure Fisher optimality, even averaged over
    parameter-prior uncertainty, had no mechanism to recognize that
    skipping periapsis entirely breaks geometric identifiability, not
    just precision. Returns (epochs, cumulative_fisher_matrix)."""
    if n_total is None:
        n_total = len(campaign.build_epoch_grid_heuristic())
    n_floor_requested = int(round(n_total * floor_fraction))

    floor_epochs = _periapsis_floor_epochs(n_floor_requested)
    # Passage 1's window straddles the visibility-season boundary (same
    # complication as campaign.py's heuristic design), so the season
    # filter above can cull roughly half of that passage's requested
    # floor points. Rather than hardcode a compensation factor, just
    # measure what actually survived and give the flexible search
    # whatever budget remains, so the total always matches n_total
    # exactly (same "verify by running it" approach campaign.py uses).
    n_flexible = n_total - len(floor_epochs)
    candidates = build_candidate_grid()
    fisher_cand = robust_fisher_per_epoch(candidates, n_draws=n_draws)
    floor_fisher = robust_fisher_per_epoch(floor_epochs, n_draws=n_draws).sum(axis=0)

    selected_idx, cumulative_fisher = greedy_d_optimal(
        fisher_cand, n_flexible, initial_fisher=floor_fisher)

    epochs = np.union1d(floor_epochs, candidates[selected_idx])
    return np.sort(epochs), cumulative_fisher


def crb_sigma(fisher_matrix):
    """Cramer-Rao bound: the analytic best-case parameter sigma implied
    by a Fisher information matrix, sqrt(diag(F^-1)). An instant,
    analytic prediction of what a full bootstrap would find -- doesn't
    replace running the real fit, but lets two designs be compared
    before paying for one."""
    return np.sqrt(np.diag(np.linalg.inv(fisher_matrix)))


def print_crb_comparison(designs):
    """Instant, analytic (Cramer-Rao bound) comparison of several
    designs' predicted parameter precision, before running any fit.
    `designs` is {label: epochs}. Uses the naive (single-point) Fisher
    matrix purely as a reporting metric here, evaluated at the truth --
    NOT the robust design's own selection criterion -- so all designs are
    scored on the same yardstick."""
    truth = fit_orbit.TRUTH_VECTOR.copy()
    sigma_ra = sigma_dec = k.GRAVITY_PLUS_ASTROMETRY_MAS

    crb = {}
    logdet = {}
    for label, epochs in designs.items():
        fisher = fisher_per_epoch(
            observable_jacobian(truth, epochs), sigma_ra, sigma_dec
        ).sum(axis=0) + RIDGE * np.eye(len(truth))
        crb[label] = crb_sigma(fisher)
        _, logdet[label] = np.linalg.slogdet(fisher)

    labels = list(designs.keys())
    print(f"\nCramer-Rao bound comparison ({', '.join(f'{len(designs[l])} {l}' for l in labels)}"
          f" epochs, analytic -- no fit run yet):")
    header = f"{'parameter':<18}" + "".join(f"{l + ' CRB':>18}" for l in labels)
    print(header)
    for i, name in enumerate(fit_orbit.PARAM_NAMES):
        row = f"{name:<18}" + "".join(f"{crb[l][i]:>18.5f}" for l in labels)
        print(row)
    print("log(det Fisher): " + ", ".join(f"{l}={logdet[l]:.2f}" for l in labels)
          + "  (higher = more informative under the naive, truth-linearized metric)")


def fit_and_report(csv_path, label):
    """Load a campaign CSV, run the same fit + bootstrap fit_orbit.py
    uses, and print the truth/fitted/sigma comparison table. Reuses
    fit_orbit's own functions so the fit itself is identical across
    designs -- only the input epochs differ."""
    data = fit_orbit.load_observations(csv_path)
    x0 = fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION
    best_fit = fit_orbit.fit_once(
        data["epoch_yr"], data["ra_offset_mas"], data["dec_offset_mas"],
        data["sigma_ra_mas"], data["sigma_dec_mas"], x0=x0,
    )
    print(f"\n[{label}] Running {fit_orbit.N_BOOTSTRAP} bootstrap resamples...")
    fit_sigma = fit_orbit.bootstrap_samples(data, best_fit).std(axis=0)
    print(f"[{label}] fit results:")
    fit_orbit.print_comparison(best_fit, fit_sigma)
    return best_fit, fit_sigma


def main():
    """Build all four designs (heuristic, naive D-optimal, pseudo-
    Bayesian robust D-optimal, and periapsis-floor-constrained D-optimal)
    at matched epoch count, compare them analytically (Cramer-Rao bound)
    and then empirically (real fits with bootstrap uncertainties) for the
    two designs that matter (heuristic vs. constrained), and report which
    one actually wins on each parameter. See module docstring for why the
    naive and robust designs are included only as documented cautionary
    results -- both tried to skip periapsis entirely and are expected to
    fail the real bootstrap test the same way."""
    heuristic_epochs = campaign.build_epoch_grid_heuristic()
    n_total = len(heuristic_epochs)
    naive_epochs, _ = build_epoch_grid_naive(n_total=n_total)
    robust_epochs, _ = build_epoch_grid_robust(n_total=n_total)
    constrained_epochs, _ = build_epoch_grid_constrained(n_total=n_total)

    print_crb_comparison({
        "heuristic": heuristic_epochs, "naive": naive_epochs,
        "robust": robust_epochs, "constrained": constrained_epochs,
    })

    campaign.report_and_write_csv(
        heuristic_epochs, "output/synthetic_observations_heuristic.csv", label="heuristic")
    campaign.report_and_write_csv(
        constrained_epochs, "output/synthetic_observations_constrained.csv", label="constrained")

    fit_and_report("output/synthetic_observations_heuristic.csv", "heuristic")
    fit_and_report("output/synthetic_observations_constrained.csv", "constrained")


if __name__ == "__main__":
    main()
