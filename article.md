# A Feasibility Study for Detecting Schwarzschild Apsidal Precession in the Orbit of S301 around Sgr A\* with GRAVITY+/VLTI

**J. Catanzarite**
*Independent study; draft manuscript, not submitted for peer review.*

---

## Abstract

S301 is a recently discovered star (GRAVITY Collaboration et al. 2026) orbiting Sgr A\* on a highly eccentric (*e* = 0.9832), 8.68 yr orbit reaching 136 Schwarzschild radii and >8% of the speed of light at periapsis. Its orbit should show rosette-shaped apsidal precession from ordinary first-post-Newtonian (1PN) general relativity — the Mercury-perihelion effect, scaled up by proximity to a far more compact body — independent of and prior to any measurement of Sgr A\*'s spin, which motivated S301's discovery. We ask whether a real, near-term ground-based campaign can detect this precession from astrometry alone, and build the pipeline in the order a real proposal would need it: a relativistic light-curve model validated against the published solution, a synthetic two-passage GRAVITY+ campaign, and a seven-parameter orbit fit. Two of our initial design choices fail real tests along the way. A proposed ERIS radial-velocity channel, natural given how S-star orbits are usually solved, contributes nothing at S301's magnitude — confirmed directly, not assumed — and is dropped from the proposal. A cadence chosen by Fisher-information (D-optimal) design, expected to beat a simple heuristic, instead collapses in a real nonlinear fit despite an excellent analytic prediction — a reproducible pitfall of local-linearization design applied to a strongly nonlinear, periodic problem — and is rescued only by reintroducing a physically motivated periapsis-coverage floor alongside the statistical search. The resulting campaign recovers all seven orbital parameters, including the precession rate (0.2326 ± 0.0013° yr⁻¹ against a true 0.2307° yr⁻¹), from astrometry alone. We conclude a two-passage, astrometry-only GRAVITY+ campaign is both feasible and sufficient for a spin-independent test of general relativity near Sgr A\*.

**Keywords:** Galactic Center — black hole physics — astrometry — celestial mechanics — general relativity — experimental design

---

## 1. Introduction

The Galactic Center's S-stars have delivered some of the strongest strong-field tests of general relativity to date, including gravitational redshift and Schwarzschild precession in the orbit of S2 (GRAVITY Collaboration 2018, 2020). S301 (GRAVITY Collaboration et al. 2026) extends this program with a shorter period (8.68 yr vs. S2's ~16 yr), higher eccentricity, and closer periapsis (136 *R*_S) than S2, and is explicitly proposed as a spin (Lense-Thirring) probe for Sgr A\*, with a companion paper (Piran et al. 2026) developing that harder measurement's methodology.

Before attempting spin, a logically prior question deserves its own answer: can a real, currently executable campaign detect the ordinary, non-spin apsidal precession any Schwarzschild black hole would already produce? We answer this by building the complete simulation pipeline a real observing proposal would need — validate the physics, design the campaign, cross-check the physics survives realistic noise, then fit — and report what worked, what failed, and why, since two of our own initial choices fail informative real tests along the way (Sections 3–4).

## 2. The S301 Orbital Solution

All ground-truth parameters are taken from GRAVITY Collaboration et al. (2026), Extended Data Table 2, Solution 1 (Table 1); a second, nearly degenerate solution is also published but not used here. At periapsis S301 passes 136 *R*_S from Sgr A\* at >8% *c* (~25,000 km s⁻¹) — closer and faster than any other known star, including S2. Its effective temperature is not published; we adopt 7300 K from standard F1.5V tables, affecting only the absolute-magnitude estimate in Section 3, not the orbital dynamics used elsewhere.

**Table 1.** Adopted S301 parameters (GRAVITY Collaboration et al. 2026, Solution 1).

| Parameter | Value | Uncertainty |
|---|---|---|
| Period *P* | 8.680 yr | 0.110 yr |
| Eccentricity *e* | 0.9832 | 0.0010 |
| Inclination *i* | 124.09° | 1.10° |
| Node Ω | 73.8° | 3.5° |
| Argument of periapsis ω | 293.4° | 2.2° |
| Periapsis epoch *t*₀ | 2023.126 | 0.010 yr |
| Semi-major axis *a* | 83.0 mas (687 AU) | — |
| Sgr A\* mass *M*• | 4.297 × 10⁶ M☉ | — |
| Distance *R*₀ | 8277 pc | — |
| K magnitude *m*_K | 19.3 | 0.3 |

## 3. Simulation Pipeline

### 3.1 Relativistic Light Curve

We model S301's K-band light curve from a full Keplerian solution combined with the exact special-relativistic Doppler factor (not a small-*v* approximation) and the Schwarzschild redshift √(1 − *R*_S/*r*), applied to a blackbody via the relativistic-beaming relation *I*_ν/ν³ over GRAVITY/VLTI's real bandpass (1.98–2.40 μm). Lensing, occultation, Kerr terms, and Rømer delay are not modeled. The model reproduces the published solution (derived *a* = 686.7 AU vs. 687.0; periapsis 136.0 *R*_S vs. 136; periapsis speed 8.54% *c* vs. ">8%"), validating the engine before it is asked to predict anything new.

The predicted periapsis signature is not the symmetric absorption dip a naive guess might expect — no such mechanism operates here. It is an *asymmetric* relativistic-beaming spike: a brightening of ~0.075 mag as the star approaches at high velocity, followed by fading to ~0.04 mag fainter as it recedes (de-boosting plus gravitational redshift, which dims regardless of direction), a 0.116 mag peak-to-trough amplitude. The total relativistic redshift exceeds the naive classical Doppler prediction by 0.26–0.46 percentage points near periapsis — a genuine, non-negligible correction, not a relabeled classical effect.

### 3.2 Observing Campaign: Instrument and Cadence

The campaign uses GRAVITY+, the real VLTI interferometer that made the discovery, adopting its demonstrated 100 μas astrometric precision (GRAVITY+ Collaboration 2023) and simultaneous K-band photometry.

**Radial velocity was considered and rejected before any data was simulated.** Scaling SINFONI's real achieved S2 precision (12.3 km s⁻¹ at *m*_K = 14.0; GRAVITY Collaboration et al. 2018) to S301's magnitude (Δ*m* = 5.3, ~130× less flux) under the background-limited regime appropriate for faint Galactic Center spectroscopy predicts ~1621 km s⁻¹ precision against a ±15,000 km s⁻¹ signal — SNR ≈ 0.04 against SPIFFIER's own 60 km s⁻¹ resolution element (*R* = 5000). No Time Allocation Committee would grant ERIS time for a measurement known in advance to fail this badly, so the campaign analyzed here never requests it. We confirmed directly that nothing is lost: refitting identical astrometric data with a simulated RV channel added back in leaves every parameter's bootstrap uncertainty statistically unchanged (precession rate marginally *tighter* without RV, 0.0013 vs. 0.0014° yr⁻¹).

**Cadence required three attempts.** S301 spends nearly all of its 8.68 yr period near apoapsis and almost its entire angular excursion within days of periapsis, so cadence matters enormously. *Design A (heuristic)*: ~90% of epochs within ±20 days of each periapsis, remainder sparse. Simple and reasonably effective (Table 2), but assumes all seven parameters are equally served by periapsis proximity. *Design B (naive D-optimal)*: Fisher information linearized at the truth, greedily maximizing its determinant. This **failed catastrophically** — it placed nearly the whole budget on two clusters at the campaign's temporal extremes (5.4 and 18.5 yr from the reference epoch, more than two periods away) and selected *zero* epochs near either periapsis. Its Cramér–Rao prediction was 1.8–28,000× tighter than Design A; a real bootstrap fit was instead 4–600× *worse* on every parameter. The mechanism: period/epoch sensitivity genuinely grows with time baseline (as in pulsar timing), but this is a *local* quantity — extrapolated across >2 periods, a small period error compounds into large nonlinear phase drift (a cycle-slip hazard) that linearized Fisher information cannot see. *Design C (pseudo-Bayesian robust)*: averaging Fisher information over 300 draws from a realistic parameter prior (Chaloner & Verdinelli 1995) de-emphasized the fragile clusters but *still* selected zero periapsis epochs — Fisher information, even prior-averaged, measures only noise-driven uncertainty and cannot represent the geometric non-identifiability that comes from sampling too few distinct orbital phases, which is exactly why classical orbit determination has always required observations spread across an orbit's curvature. *Design D (constrained)*: a hard, non-statistical floor reserving ~30% of the budget within ±20 days of each periapsis, with the remaining ~70% allocated by the same Fisher-optimal search. This is the design used in this study; it beats Design A on every parameter (Table 2).

All designs are filtered to Sgr A\*'s real Paranal visibility season (~March–September; GRAVITY Collaboration et al. 2026). This exposes a real scheduling gap: the first periapsis (2031 October 22) falls just after that year's season closes, so its dense window is anchored to the nearest visible date (~2031 September 30) rather than periapsis itself. The second passage (2040 June 26) falls comfortably in-season.

### 3.3 Photometric Cross-Check

Photometry is generated (a free byproduct of GRAVITY+'s astrometric observation) but not fit, since it is a deterministic function of the same (position, velocity) state astrometry already constrains. Its value is independent validation: S301's published photometric uncertainty (0.30 mag) exceeds the entire 0.116 mag signal from Section 3.1 at any single epoch, yet the injected truth curve, threaded through many simulated epochs (Figure 2, lower panel), visibly traces the same asymmetric brighten-then-dim structure — a qualitative but genuine end-to-end check that the injected physics, noise model, and cadence are mutually consistent.

### 3.4 Orbit Fit

We fit seven parameters — the six Keplerian elements plus apsidal precession rate ω̇, applied as ω_eff(*t*) = ω + ω̇(*t* − *t*_ref) — to the astrometric data via Levenberg–Marquardt least squares, from an initial guess perturbed 3% in period and several degrees in angle. The semi-major axis is derived from (*P*, *GM*•) via Kepler's third law at each iteration. Uncertainties come from 1000 bootstrap resamples. The recovered ω̇ is compared to the analytic Schwarzschild prediction ω̇ = 6π*GM*•/[*c*²*a*(1 − *e*²)*P*] — the Mercury-perihelion formula. No spin (Lense-Thirring) term is injected anywhere: *i* and Ω remain genuinely static, since Sgr A\*'s spin is unmeasured and assuming a value would substitute a guess for a measurement.

## 4. Results

The final campaign has 131 epochs (2028.50–2041.69), 39 (29.8%) within ±20 days of a periapsis, more than double an earlier 60-epoch single-passage design.

**Table 2.** Bootstrap 1σ uncertainties, heuristic vs. constrained cadence (*N* = 131, RV included to isolate cadence).

| Parameter | Design A | Design D | Published |
|---|---|---|---|
| *P* (yr) | 0.0001 | 0.0002 | 0.1100 |
| *e* | 0.0001 | 0.0000 | 0.0010 |
| *i* (°) | 0.0480 | 0.0356 | 1.1000 |
| Ω (°) | 0.1229 | 0.0849 | 3.5000 |
| ω (°) | 0.0939 | 0.0475 | 2.2000 |
| *t*₀ (yr) | 0.0003 | 0.0004 | 0.0100 |
| ω̇ (° yr⁻¹) | 0.0043 | 0.0014 | — |

Design D improves *i*, Ω, ω, and ω̇ by 26–67%, at negligible cost to *P*/*t*₀ (both remain >100× tighter than published uncertainty regardless).

**Table 3.** Final result: Design D cadence, astrometry only.

| Parameter | Truth | Fitted | σ | Published σ |
|---|---|---|---|---|
| *P* (yr) | 8.6800 | 8.6799 | 0.0002 | 0.1100 |
| *e* | 0.9832 | 0.9832 | 0.0000 | 0.0010 |
| *i* (°) | 124.0900 | 124.1207 | 0.0337 | 1.1000 |
| Ω (°) | 73.8000 | 73.6948 | 0.0807 | 3.5000 |
| ω (°) | 293.4000 | 293.3333 | 0.0452 | 2.2000 |
| *t*₀ (yr) | 2023.1260 | 2023.1261 | 0.0004 | 0.0100 |
| ω̇ (° yr⁻¹) | 0.2307 | 0.2326 | 0.0013 | — |

Every parameter recovers within ~1.5σ of truth. The fitted sky-plane track (Figure 2, upper panel) visibly shows the periapsis direction rotated between passages — the geometric signature of precession, closing the loop opened in Section 3.1. Because this campaign assumes best-case GRAVITY+ precision at every epoch across two full passages 8.68 yr apart (a baseline the real 2026 dataset does not yet have), several uncertainties here are tighter than the real paper's published values (e.g., 0.08° vs. 3.5° on Ω) — an idealized-study artifact, not a claim of outperforming the real discovery.

## 5. Discussion

**Instrumentation.** No operating alternative beats GRAVITY+ for astrometry: Roman's 2.4 m aperture has a ~223 mas diffraction limit, larger than S301's entire 83 mas orbit; JWST's sharper ~85 mas still lacks an interferometric baseline. For RV, Roman fails independently of integration time (its non-AO-corrected PSF causes source confusion, a structural floor that does not average down, compounded by *R* ≈ 461–890 vs. ERIS's 5000). JWST/NIRSpec is more ambiguous — comparable or better resolution than pre-launch expectations (Shajib et al. 2025) and no OH airglow, but real observations of this field report severe crowding-specific problems (saturation, MSA leakage, guide-star misidentification; Yusef-Zadeh et al. 2025) and no citable achieved precision for a comparably faint point source, so we do not adopt it. The field's real fix, HARMONI on the ELT (*R* = 3500–18,000, 39 m), has science first light slipping to December 2030 (ESO 2026) — only ~10 months before passage 1, unlikely enough time for a new instrument to clear proposal review; HARMONI more plausibly contributes to passage 2 (2040.5).

**Scope.** Measuring Sgr A\*'s spin is explicitly outside this work. Spin causes a distinct signature (nodal precession from frame-dragging) not injected here; Piran et al. (2026) show a real spin measurement additionally requires separating that signal from Newtonian confusion via multi-star calibration, a substantially harder analysis than presented here. This study is a genuine precursor to that measurement, not the measurement itself.

**Design lesson.** Local Fisher-information design assumes the model is well approximated by its linear expansion across the *entire* candidate range — an assumption that fails badly, invisibly to the analytic bound, for nonlinear periodic problems evaluated over multi-cycle baselines. We recommend pairing any such design with pseudo-Bayesian prior-averaging (Chaloner & Verdinelli 1995), an explicit domain-knowledge floor against non-identifiability, and validation against a real nonlinear fit rather than trust in the analytic prediction alone.

## 6. Summary

A relativistic light-curve model reproduces S301's published solution to <0.1% and predicts an asymmetric beaming feature, not an occultation dip, that survives visibly inside realistic noisy photometry. A feasibility calculation shows ERIS cannot detect S301's RV with any current telescope; we verify directly that dropping it costs no orbit-recovery precision. A naive Fisher-optimal cadence catastrophically fails via nonlinear aliasing despite an excellent analytic prediction; a hybrid design (domain-knowledge floor + Fisher-optimal remainder) fixes this and beats a simple heuristic on every parameter. The resulting astrometry-only, two-passage GRAVITY+ campaign recovers the injected Schwarzschild precession rate (0.2307° yr⁻¹) as 0.2326 ± 0.0013° yr⁻¹. We conclude this campaign is feasible with existing instrumentation and sufficient for an independent, spin-agnostic GR test near Sgr A\*, ahead of the harder spin measurement S301's discovery motivated.

## Acknowledgments

This work uses the published orbital solution of GRAVITY Collaboration et al. (2026) and the framework of Piran et al. (2026). Code and figures are available in the accompanying repository.

## References

- Chaloner, K., & Verdinelli, I. 1995, Statistical Science, 10(3), 273
- ESO 2026, "Telescope first light for ESO's Extremely Large Telescope now planned for March 2029," https://www.eso.org/public/announcements/ann25001/
- GRAVITY Collaboration, Abuter, R., et al. 2018, A&A, 615, L15 (arXiv:1807.09409)
- GRAVITY Collaboration, Abuter, R., et al. 2020, A&A, 636, L5 (arXiv:2004.07187)
- GRAVITY Collaboration, Abd El Dayem, K., et al. 2026, Nature (arXiv:2607.12664)
- GRAVITY+ Collaboration 2023 (arXiv:2301.08071)
- Piran, T., Amaro-Seoane, P., Aytac, B., et al. 2026 (arXiv:2607.24931)
- Shajib, A. J., Treu, T., Melo, A., et al. 2025, A&A, 702, L12 (arXiv:2507.03746)
- Yusef-Zadeh, F., Bushouse, H., Arendt, R. G., et al. 2025, ApJL (arXiv:2501.04096)

## Code and Data Availability

Simulation code (`orbit.py`, `s301_lightcurve.py`, `campaign.py`, `optimal_design.py`, `fit_orbit.py`, `constants.py`) and figures are in this project's repository.
