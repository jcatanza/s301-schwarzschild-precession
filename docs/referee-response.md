# Response to the internal referee report

**Report:** adversarial self-review of `article.tex` and the repository at
commit `a0049d0` (2026-09-06), recommendation *major revision*.
**Response covers:** commits `08566f6` (manuscript revision draft),
`eb860ee` (editorial pass, new title), `3a78809` (cautious error budget) and
the agent-branch merges between them (`8af286c`, `997c4c6`, `1aaeb0c`).
**Status:** every finding addressed; nothing deferred.

Every number quoted below is produced by a script in this repository and
enters the manuscript as a macro from `numbers.tex` (see M6), so this table
can be re-verified with `./run_all.sh`.

---

## Major issues

| # | Finding | Response | Where |
|---|---------|----------|-------|
| **M1** | The "catastrophic failure" of the Fisher-optimal cadence (Design B piles epochs at the temporal extremes, selects zero periapsis epochs, 4–600× worse) is not reproduced by the submitted code. | **Accepted.** The failure belonged to an earlier code state (no visibility-season filter on the candidate grid, unsafeguarded Kepler solver). `optimal_design.py` now bootstraps all four designs plus the Cramér–Rao bound; B, C and D agree on $\dot\omega$ within bootstrap noise (B/D = 0.95) and all beat the heuristic A (D/A = 3.85). §3.3 tells this honestly and keeps the periapsis floor as insurance, not as a rescue. The abstract, the "design lesson" and the Conclusions no longer claim a catastrophe. | `optimal_design.py`; §3.3 "Cadence design"; Table 2; §6.4 "Lessons" |
| **M2** | The dominant systematic imports the NACO-to-radio frame tie (Plewa et al. 2015), which does not apply to interferometric astrometry referenced directly to Sgr A*. | **Accepted.** `instrumental_error.py` deleted; replaced by `reference_frame_error.py`, which models what does apply: a constant offset between the Sgr A* NIR photocentre and the mass centroid (50 µas prior, 100 µas conservative, anchored to $R_S = 10$ µas and the EHT ring) and per-epoch reference jitter. The offset biases $\dot\omega$ by 0.2–0.4σ if ignored and costs ×1.04 in precision once fitted. The GRAVITY 2020/2022 and discovery-paper statements are quoted in §5.2. The bottom line moved from ×2.5 (frame tie) to ×1.04 (zero point) — the paper was pessimistic, as the report anticipated. | `reference_frame_error.py`; §5.2 "The Sgr A* reference zero point"; Table 4 |
| **M3** | "Demonstrated 100 µas" is the GRAVITY+ forecast; the achieved S301 precision is ≈207 µas. | **Accepted and taken further.** The fiducial per-epoch precision is now 207 µas everywhere (`constants.GRAVITY_PLUS_ASTROMETRY_MAS`); 100 µas is the optimistic sensitivity case. `precision_sensitivity.py` runs both on the same epochs and seed: 0.0023 vs 0.0012 deg/yr (96σ vs 198σ white-noise detections). | `constants.py`; `precision_sensitivity.py`; §3.2; §4.1; Table 4 |
| **M4** | Synthetic photometry used a bolometric proxy, $-2.5\log_{10}(1+z)^{-3}$, with twice the amplitude (0.227 mag) of the validated band-integrated model (0.116 mag). | **Accepted; a real bug.** New `photometry.py` (Planck × $I_\nu/\nu^3$ integrated over the 1.98–2.40 µm top-hat) is the single model used by `s301_lightcurve.py`, `campaign.true_observables` and `photometry_checks.py`. Figure 4 regenerated; the $\chi^2$ check now compares the astrometric orbit's prediction with photometry generated from the same physics ($\chi^2/N = 1.05$, 0.40σ from expectation). | `photometry.py`; `campaign.py`; `photometry_checks.py`; §3.1; §4.3 |
| **M5** | About a dozen quantitative claims had no script (RV added back, joint photometry, 2032 reallocation, complete-passage hypothetical, placement scan, noisy Rømer bootstrap, 262-epoch test). | **Accepted.** Every claim now maps to a script: `cadence_alternatives.py` (reallocation, complete passage 1, doubled budget, placement scan, 30-realization Monte Carlo), `photometry_checks.py`, `rv_channel_test.py`, `roemer_delay.py` (noisy bootstrap), `precision_sensitivity.py`, `mass_distance_test.py`, `lensing_error.py`, `solar_conjunction.py`, `correlated_noise_test.py`, `error_budget.py`. `run_all.sh` runs them in dependency order; the README maps each to its section. | `run_all.sh`; README "Scripts" table |
| **M6** | Stale numbers survived the pipeline rerun in §4.1, §4.2, §4.5, the Table 2 caption and §3.2; three scripts used a 200-resample baseline while Table 3 used 1000. | **Accepted; fixed structurally.** Every script writes `results/<name>.json` through `results_io.write_results`; `make_numbers.py` renders all of them into `numbers.tex` macros that `article.tex` `\input`s. No number is typed by hand (the only hand-typed values are published constants). All σ-multiples now reference the single Table 3 bootstrap σ read from `results/fit_orbit.json`. The Rømer numbers (−0.0023 deg/yr, 1.0σ at 207 µas) and every other stale entry regenerate on each run. | `results_io.py`; `make_numbers.py`; `numbers.tex`; `CLAUDE.md` ground rules |
| **M7** | Holding $M_\bullet$ and $R_0$ fixed over-constrains the fit. | **Accepted.** `mass_distance_test.py` refits with both free under GRAVITY 2022 priors (statistical ⊕ systematic: 0.97% and 0.38%); $e$, $i$, $\omega$ widen substantially, $\dot\omega$ by ×1.04. `error_budget.py` fits orbit + offset + $M$ + $R_0$ jointly (11 parameters) rather than multiplying separate factors. Priors verified against arXiv:2112.07478. | `mass_distance_test.py`; `error_budget.py`; §5.3; §5.9 |
| **M8** | The discovery team's ELT/MICADO radial-velocity plan (1 km/s assumed) was not engaged; ERIS statement not cited. | **Accepted.** §3.2 quotes "Our current ERIS spectroscopy is not deep enough"; §6.2 adds MICADO, states that the 1 km/s figure carries no sensitivity calculation, contrasts it with our HARMONI scaling (~20 km/s), and notes that either would help $P$, $e$, $t_0$ and the mass–distance degeneracy but not $\dot\omega$ (`rv_channel_test.py`: −4.1% at HARMONI precision). §3.5 notes that the discovery fit already includes the Rømer delay and Schwarzschild precession. | `rv_channel_test.py`; §3.2; §6.2 |

## Moderate issues

| # | Finding | Response | Where |
|---|---------|----------|-------|
| 1 | Stochastic (discrete) perturbers absent from the astrophysical budget. | Paragraph added quoting the discovery paper's heavy-tailed-perturbation statement and its $\gtrsim 10^7$ yr relaxation/collision time-scales against the $1.6\times10^3$ yr precession time-scale; flagged as a rare, unbounded event that no smooth budget captures. | §5.4 "Extended mass", second paragraph |
| 2 | Gravitational lensing flagged as "unverified" when a standard estimate exists. | New `lensing_error.py`: exact point-lens primary-image shift applied epoch by epoch; largest shift 22 µas at the far-side periapsis of passage 2, one epoch above 10 µas, $\dot\omega$ bias 0.00σ. Now a row in Table 4. | `lensing_error.py`; §5.6 |
| 3 | Bootstrap uncertainties should be cross-checked against the CRB and a multi-realization Monte Carlo. | CRB reported for Design D in Table 2 (CRB/bootstrap = 0.93 for $\dot\omega$); 30 independent noise realizations in `cadence_alternatives.py` (MC/bootstrap = 0.78); 40 realizations per case in `correlated_noise_test.py`. | Table 2; §3.3; §4.2; §5.7 |
| 4 | Bandpass in code (Gaussian, 1.65–2.73 µm) is not the bandpass in the text (1.98–2.40 µm top-hat). | Code now uses the top-hat over the stated range (`photometry.py`, `BAND_LO_UM`/`BAND_HI_UM` written to `results/s301_lightcurve.json` and quoted by macro). | `photometry.py`; §3.1 |
| 5 | Length and structure: 420-word abstract; 700-word Summary paragraph; project history in §3.3; decorative figures; error budget buried in Results. | Abstract ≤ 250 words; Conclusions as a six-point list; history cut to one paragraph; error budget promoted to its own section (§5) and reordered into modelled → bounded → spin → combined tiers; MNRAS-style sentence-case heads, Acknowledgements and Data Availability sections. Figures 1 and 2 retained by author decision (orientation for a general astronomer reader). | §5; §7; `eb860ee` |

## Beyond the report

The report's central worry — that the precision claim rests on white noise — was pursued further than the report asked, after the author requested a cautious budget:

- **Correlated systematics** (`correlated_noise_test.py`, §5.7): per-observing-run common-mode offsets (60 µas, GRAVITY 2022 calibration systematic) and a red-noise reference term (50 µas, one-season correlation time, GRAVITY 2018 flare-centre bound) on top of the white noise, evaluated from independent realizations. They double $\sigma(\dot\omega)$ (×2.09), supply 77% of the final variance, and the epoch bootstrap applied to a correlated realization reports only 0.71× the true scatter.
- **Unknown-unknowns allowance** (`error_budget.py`, §5.9): inflation by $\sqrt{\chi^2_r}$ = 1.22–1.47 from the reduced $\chi^2$ of real GRAVITY orbit fits (2020: 1.5; 2022: 2.17).
- **Headline:** 0.0074 deg/yr, 3.2% of the signal, a 31σ detection — where the white-noise number alone would have claimed 96σ (or 198σ at the 100 µas forecast the report flagged in M3).

## Repository review

| Item | Response |
|------|----------|
| Licence missing | `LICENSE` (MIT, code) and `LICENSE-MANUSCRIPT.md` (CC BY 4.0); stated in the Data Availability section. |
| Dependencies unpinned | `requirements.txt` pinned, including astropy for `solar_conjunction.py`. |
| No tests | `tests/`: Kepler-solver divergence regression ($M = 0.2185$, $e = 0.9832$), random sweep, vectorized/scalar agreement, 1PN rate, Lense–Thirring per-orbit, light-curve amplitude, band model vs bolometric proxy, shared-photometry-model check, campaign → fit smoke test, Rømer bound. `pytest -q`: 12 passed. |
| Entry point / run order | `run_all.sh` regenerates every result, figure, `numbers.tex`, `article.pdf` (Docker TeX Live) and `arxiv/` in dependency order; `--skip-slow` mode; README script-to-section table. |
| Three diverged manuscript copies | `article.tex` is canonical; `article.md` deleted; README rewritten as a project README; `CLAUDE.md` rewritten; `arxiv/README.md` de-pathed. |
| Lint claim | `pylint --rcfile=.pylintrc *.py` = 10.00 across all modules (`precession_movie.py` refactored). |
| Committed artefacts | `.gitignore` and README now agree: `output/*.png|jpg|mp4`, `article.pdf`, `results/*.json` and `arxiv/` are deliberately committed. |
| Coursework residue | Removed from `orbit.py`, `fit_orbit.py`, `constants.py`. |
| Stale docstrings | `extended_mass_error.py`, `campaign.py`, `optimal_design.py` floor comment (40% requested, ~30% surviving) updated. |

## Minor points

- "Rømer" throughout.
- §2 quotes the periapsis speed from the code (`\lightcurveVpKms` = 25 600 km/s) and substantiates "fastest known star in the Galaxy" against S2, S4714 (Peißker et al. 2020) and S5-HVS1 (Koposov et al. 2020) from primary sources; "most eccentric" is deliberately *not* claimed.
- HARMONI scaling assumption ($S/N \propto D^2$, background-limited, diffraction-limited Strehl) stated in §6.2.
- The ~130 m baseline is now introduced in §3.2 where GRAVITY+ is described.
- Table 4 rows are ordered as the subsections and the note says which rows subsume which; the "SR/GR" caption ambiguity fixed (GR no longer used for both "gravitational redshift" and "general relativity").
- Keywords left as UAT terms while the class is `aastex701`; to be replaced on the switch to the `mnras` class.

## Two report items not adopted

- **Remove Figures 1 and 2.** Retained: the author wants the EHT shadow and the crowded-field image as scale and context for readers outside the S-star community; captions shortened.
- **"Spin-agnostic vs spin-free" phrasing.** The report accepted the distinction; the paper now uses the single term *spin-agnostic*, defined once in §6.3.
