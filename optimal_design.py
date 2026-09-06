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

HISTORY, stated honestly. An earlier state of this module -- before the
candidate grid was filtered to the Mar-Sep visibility season, and before
orbit.py's Kepler solver was bracketed -- produced a naive (single-point
linearized) D-optimal design that failed catastrophically: it piled
nearly all epochs onto two calendar-adjacent clusters ~2 orbital
periods from the reference epoch, covered neither periapsis, promised
800-28000x Cramer-Rao gains, and came out 4-600x WORSE than the
heuristic in a real bootstrap fit. The mechanism we diagnosed then is
real physics: sensitivity to P and t_peri grows with baseline, but a
small period error extrapolated over >2 periods compounds into a
nonlinear phase drift (a "cycle-slip"/aliasing hazard) that a Fisher
matrix linearized at the exact truth cannot see. That diagnosis
motivated build_epoch_grid_robust (pseudo-Bayesian D-optimal, Chaloner &
Verdinelli 1995: average the per-epoch Fisher matrix over parameter
draws from the published-uncertainty prior, so fragile information is
averaged down) and, when the robust design in that same code state
still chose zero periapsis epochs, the hard periapsis floor of
build_epoch_grid_constrained.

WHAT THE CURRENT CODE ACTUALLY FINDS (main(), fit_orbit.N_BOOTSTRAP =
1000 resamples per design, all four designs at 131 epochs; every number
is in results/optimal_design.json): the catastrophic failure does NOT
recur. With the season-filtered candidate grid and the fixed solver,
the naive design (B) places 27 of 131 epochs within +/-20 d of a
periapsis anchor across six calendar years, the robust design (C)
places 28 across seven, and the constrained design (D) places 50
(its 39 floor epochs plus 11 the optimizer added on its own). Their
Cramer-Rao bounds on omega_dot are 0.00093 / 0.00099 / 0.00104 deg/yr
against the heuristic's (A) 0.00360 -- a B-over-A gain of 3.9x, not
thousands -- and the bootstrap confirms them closely: sigma(omega_dot)
= 0.00430 (A), 0.00106 (B), 0.00109 (C), 0.00112 (D) deg/yr, with
CRB/bootstrap = 0.93 for D. B, C and D are statistically
indistinguishable on the science target (D/B = 0.95); A is ~4x worse.
On the individual elements B is marginally the tightest on the angles
and D the tightest on P and t_peri (its guaranteed periapsis coverage
buys timing), but all three are within ~30% of each other on every
parameter. Design D is the adopted campaign (campaign.build_epoch_grid)
because its periapsis coverage is guaranteed by construction rather
than being an emergent property of the optimizer that, as the earlier
code state showed, can evaporate when the candidate grid or the model
changes -- not because B or C fail. The earlier "catastrophic failure"
narrative was a pre-fix artifact and is retained above only as history.
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np  # noqa: E402  pylint: disable=wrong-import-position

import campaign  # noqa: E402  pylint: disable=wrong-import-position
import constants as k  # noqa: E402  pylint: disable=wrong-import-position
import fit_orbit  # noqa: E402  pylint: disable=wrong-import-position
import results_io  # noqa: E402  pylint: disable=wrong-import-position

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
    This guards against the cycle-count aliasing hazard described in the
    module docstring: a candidate epoch whose apparent information
    depends on knowing P/t_peri to unrealistic precision gets averaged
    down, because nearby-but-different draws see very different (often
    much smaller, sometimes inconsistent) sensitivity there. (In the
    current, season-filtered code the single-point design no longer
    falls into that trap, so the averaging is insurance rather than a
    fix -- see module docstring.)"""
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
    """NAIVE D-optimal epoch grid (Design B of the manuscript's four-way
    comparison): greedily selects n_total epochs from a fine visibility-
    filtered candidate grid, maximizing the Fisher-information
    determinant linearized at a single fiducial TRUTH parameter vector.
    An earlier code state of this design (candidate grid not yet
    season-filtered, Kepler solver not yet bracketed) failed the real
    bootstrap test catastrophically; the current version places 27 of
    131 epochs within +/-20 d of a periapsis anchor across six calendar
    years and bootstraps as well as the adopted constrained design --
    see module docstring. Returns (epochs, cumulative_fisher_matrix)."""
    if n_total is None:
        n_total = len(campaign.build_epoch_grid_heuristic())

    candidates = build_candidate_grid()
    jac = observable_jacobian(fit_orbit.TRUTH_VECTOR.copy(), candidates)
    fisher_cand = fisher_per_epoch(jac, k.GRAVITY_PLUS_ASTROMETRY_MAS, k.GRAVITY_PLUS_ASTROMETRY_MAS)

    selected_idx, cumulative_fisher = greedy_d_optimal(fisher_cand, n_total)
    return np.sort(candidates[selected_idx]), cumulative_fisher


def build_epoch_grid_robust(n_total=None, n_draws=N_MC_DRAWS):
    """Pseudo-Bayesian D-optimal epoch grid (Design C): greedily selects
    n_total epochs from a fine visibility-filtered candidate grid,
    maximizing the determinant of the Fisher matrix AVERAGED over a
    realistic prior on the parameters (robust_fisher_per_epoch), not
    linearized at one exact truth point -- insurance against the
    aliasing hazard the naive design showed in an earlier code state
    (see module docstring). In the current code it places 28 of 131
    epochs within +/-20 d of a periapsis anchor without being told to.
    Returns (epochs, cumulative_fisher_matrix)."""
    if n_total is None:
        n_total = len(campaign.build_epoch_grid_heuristic())

    candidates = build_candidate_grid()
    fisher_cand = robust_fisher_per_epoch(candidates, n_draws=n_draws)

    selected_idx, cumulative_fisher = greedy_d_optimal(fisher_cand, n_total)
    return np.sort(candidates[selected_idx]), cumulative_fisher


PERIAPSIS_FLOOR_FRACTION = 0.4
# Fraction of the epoch budget REQUESTED, up front, for dense coverage
# within +/-20 d of each periapsis passage (split evenly between the
# two). This is the requested figure, not the realized one: passage 1's
# window straddles the end of the 2031 visibility season, so the season
# filter culls roughly half of that passage's share, and ~30% of the
# budget actually survives as floor epochs (for n_total=131: 52
# requested, 39 survive -- 13 at passage 1, 26 at passage 2). The
# Fisher-optimal search then places some of its own flexible epochs in
# the dense windows as well, so the design's total periapsis occupancy
# ends up at ~38% (n_dense_d in results/optimal_design.json).
#
# Why a hard floor at all, given that the unconstrained designs below
# now also cover periapsis (see module docstring)? Fisher information,
# even pseudo-Bayesian-averaged over parameter uncertainty
# (build_epoch_grid_robust), only measures NOISE-driven uncertainty
# around a presumed-correct optimum -- it has no term that would
# penalize STRUCTURAL/GEOMETRIC non-identifiability from sampling too
# few genuinely distinct orbital phases. Classical orbit determination
# (Gauss's method and its descendants) has always required observations
# spread across an orbit, not just high precision at a couple of points,
# because a highly eccentric orbit's shape is resolved by its curvature
# -- which periapsis alone provides. Whether an unconstrained optimizer
# covers periapsis is therefore an unguaranteed by-product of the
# candidate grid and the model (an earlier code state of this module
# produced designs with zero periapsis epochs; the current one produces
# ~27-28 of 131). The floor encodes the requirement directly, so the
# adopted campaign's periapsis coverage is guaranteed by construction
# rather than by luck; the Fisher-optimal search is then used for what
# it is legitimately good at: allocating the REMAINING flexible budget
# to best constrain the angular-orientation elements, where the
# uniform-allocation heuristic underperforms.


def periapsis_anchors():
    """The two dense-window anchor dates: each periapsis passage,
    season-adjusted (campaign._season_anchor) so passage 1's window --
    which straddles the end of the 2031 visibility season -- is centred
    on the last visible date rather than on the (unobservable) periapsis
    itself. Same construction as campaign.py's heuristic design and its
    printed occupancy check."""
    # pylint: disable=protected-access
    return np.array([campaign._season_anchor(k.NEXT_PERIAPSIS_YR),
                     campaign._season_anchor(campaign.SECOND_PERIAPSIS_YR)])
    # pylint: enable=protected-access


def near_periapsis_mask(epochs):
    """True where an epoch lies within +/-DENSE_HALF_WIDTH_DAYS of either
    periapsis anchor -- the exact test campaign.report_and_write_csv's
    occupancy check applies, factored out so every design is scored
    with the same yardstick."""
    distance = np.min(np.abs(np.asarray(epochs)[:, None] - periapsis_anchors()[None, :]), axis=1)
    return distance <= campaign.DENSE_HALF_WIDTH_YR


def _periapsis_floor_epochs(n_floor):
    """n_floor epochs split evenly between the two periapsis passages'
    +/-20-day dense windows (season-adjusted anchor for passage 1, same
    construction as campaign.py's heuristic design)."""
    n_per_passage = n_floor // 2
    windows = [np.linspace(a - campaign.DENSE_HALF_WIDTH_YR, a + campaign.DENSE_HALF_WIDTH_YR,
                            n_per_passage) for a in periapsis_anchors()]
    epochs = np.concatenate(windows)
    return epochs[campaign._in_visibility_season(epochs)]  # pylint: disable=protected-access


def build_epoch_grid_constrained(n_total=None, floor_fraction=PERIAPSIS_FLOOR_FRACTION,
                                  n_draws=N_MC_DRAWS):
    """Constrained pseudo-Bayesian D-optimal design (Design D, the
    adopted campaign -- campaign.build_epoch_grid): a fixed periapsis
    floor (domain knowledge) plus Fisher-optimal allocation of the
    remaining flexible budget (robust_fisher_per_epoch) across the full
    visibility-filtered candidate grid. The floor guarantees periapsis
    coverage by construction rather than leaving it to the optimizer,
    which has no mechanism to recognize that skipping periapsis would
    break geometric identifiability, not just precision (see the
    PERIAPSIS_FLOOR_FRACTION comment). Returns (epochs,
    cumulative_fisher_matrix)."""
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
    scored on the same yardstick. Returns ({label: crb_sigma_vector},
    {label: log det Fisher}) so main() can record them alongside the
    bootstrap sigmas."""
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
    return crb, logdet


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


DESIGN_TAGS = {"heuristic": "a", "naive": "b", "robust": "c", "constrained": "d"}
# Manuscript labels: Design A = heuristic, B = naive D-optimal, C =
# pseudo-Bayesian robust D-optimal, D = periapsis-floor-constrained
# (the adopted campaign). results/optimal_design.json keys are suffixed
# with these single-letter tags.
OMEGA_DOT_INDEX = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")
RESULT_FMT = {name: ".5f" if name in ("omega_dot_deg_yr", "t_peri_yr") else ".4f"
              for name in fit_orbit.PARAM_NAMES}


def evaluate_design(label, epochs, crb_vec, logdet_val):
    """Full empirical evaluation of one design: generate its synthetic
    campaign with campaign.py's exact noise model, fit + bootstrap it
    with fit_orbit.py's exact machinery, and score its periapsis coverage
    (near_periapsis_mask) and temporal spread (number of distinct
    calendar years occupied). Returns (results dict keyed with the
    design's manuscript tag, bootstrap sigma vector)."""
    tag = DESIGN_TAGS[label]
    csv_path = f"output/synthetic_observations_{label}.csv"
    campaign.report_and_write_csv(epochs, csv_path, label=label)
    best_fit, fit_sigma = fit_and_report(csv_path, label)

    n_dense = int(near_periapsis_mask(epochs).sum())
    years = np.unique(np.floor(epochs).astype(int))
    print(f"[{label}] {n_dense} of {len(epochs)} epochs within "
          f"+/-{campaign.DENSE_HALF_WIDTH_DAYS:.0f} d of a periapsis anchor; "
          f"{len(years)} distinct calendar years: {years.tolist()}")

    results = {
        f"n_epochs_{tag}": len(epochs),
        f"n_dense_{tag}": n_dense,
        f"n_clusters_{tag}": len(years),
        f"bias_omega_dot_{tag}": (best_fit[OMEGA_DOT_INDEX] - fit_orbit.TRUTH_VECTOR[OMEGA_DOT_INDEX], ".4f"),
        f"logdet_{tag}": (logdet_val, ".1f"),
    }
    for name, sigma, crb_val in zip(fit_orbit.PARAM_NAMES, fit_sigma, crb_vec):
        results[f"sigma_{name}_{tag}"] = (sigma, RESULT_FMT[name])
        results[f"crb_{name}_{tag}"] = (crb_val, RESULT_FMT[name])
    return results, fit_sigma


def print_design_summary(designs, results, boot_sigma, crb):
    """One-table recap of the four designs on the quantities the
    manuscript compares them on."""
    print(f"\n{'design':<14}{'tag':>4}{'N':>5}{'n_dense':>9}{'n_years':>9}"
          f"{'boot sig(omega_dot)':>21}{'CRB sig(omega_dot)':>20}{'bias(omega_dot)':>17}")
    for label, epochs in designs.items():
        tag = DESIGN_TAGS[label]
        print(f"{label:<14}{tag.upper():>4}{len(epochs):>5}{results[f'n_dense_{tag}']:>9}"
              f"{results[f'n_clusters_{tag}']:>9}{boot_sigma[tag][OMEGA_DOT_INDEX]:>21.5f}"
              f"{crb[label][OMEGA_DOT_INDEX]:>20.5f}{results[f'bias_omega_dot_{tag}'][0]:>+17.4f}")


def main():
    """Build all four designs (heuristic, naive D-optimal, pseudo-
    Bayesian robust D-optimal, and periapsis-floor-constrained D-optimal)
    at matched epoch count, compare them analytically (Cramer-Rao bound)
    and then empirically -- a real fit plus fit_orbit.N_BOOTSTRAP
    bootstrap resamples for EACH of the four, on synthetic data generated
    with campaign.py's exact noise model -- and record every number in
    results/optimal_design.json. See the module docstring for what the
    comparison actually shows in the current code."""
    heuristic_epochs = campaign.build_epoch_grid_heuristic()
    n_total = len(heuristic_epochs)
    naive_epochs, _ = build_epoch_grid_naive(n_total=n_total)
    robust_epochs, _ = build_epoch_grid_robust(n_total=n_total)
    constrained_epochs, _ = build_epoch_grid_constrained(n_total=n_total)
    designs = {
        "heuristic": heuristic_epochs, "naive": naive_epochs,
        "robust": robust_epochs, "constrained": constrained_epochs,
    }

    crb, logdet = print_crb_comparison(designs)

    results = {"n_boot": fit_orbit.N_BOOTSTRAP,
               "floor_fraction_requested": (PERIAPSIS_FLOOR_FRACTION, ".2f")}
    boot_sigma = {}
    for label, epochs in designs.items():
        design_results, boot_sigma[DESIGN_TAGS[label]] = evaluate_design(
            label, epochs, crb[label], logdet[label])
        results.update(design_results)

    idx = OMEGA_DOT_INDEX
    results.update({
        "floor_fraction_realized": (results["n_dense_d"] / n_total, ".3f"),
        "crb_ratio_omega_dot_b_over_a": (crb["heuristic"][idx] / crb["naive"][idx], ".2f"),
        "boot_ratio_omega_dot_b_over_a": (boot_sigma["a"][idx] / boot_sigma["b"][idx], ".2f"),
        "boot_ratio_omega_dot_d_over_a": (boot_sigma["a"][idx] / boot_sigma["d"][idx], ".2f"),
        "boot_ratio_omega_dot_d_over_b": (boot_sigma["b"][idx] / boot_sigma["d"][idx], ".2f"),
        "crb_over_boot_omega_dot_d": (crb["constrained"][idx] / boot_sigma["d"][idx], ".2f"),
    })
    print_design_summary(designs, results, boot_sigma, crb)
    path = results_io.write_results("optimal_design", results)
    print(f"Wrote {len(results)} results to {path}")


if __name__ == "__main__":
    main()
