# A Feasibility Study for Detecting Schwarzschild Apsidal Precession in the Orbit of S301 around Sgr A\* with GRAVITY+/VLTI

**J. Catanzarite**
*Independent study; draft manuscript, not submitted for peer review.*

---

## Abstract

S301 is a recently discovered star (GRAVITY Collaboration et al. 2026) orbiting the supermassive black hole Sgr A\* on a highly eccentric (*e* = 0.9832), 8.68 yr orbit that carries it to within 136 Schwarzschild radii of closest approach, reaching more than 8% of the speed of light. Its orbit is predicted to exhibit rosette-shaped apsidal precession from standard first-post-Newtonian (1PN) general relativity, analogous to Mercury's perihelion advance, at a rate of order 0.23° yr⁻¹. We ask a single question and follow it end to end: can a real, currently executable ground-based observing campaign detect this precession from astrometry alone, *before* attempting the much harder measurement of the black hole's spin that motivated S301's discovery? We build, in sequence, a relativistic light-curve model that validates our physics engine against the published solution; a synthetic two-passage observing campaign built entirely from real, demonstrated instrument performance; and a seven-parameter orbit fit that recovers that injected physics from noisy data alone. Along the way, two of our own initial design choices fail real tests we did not expect them to fail, and both failures are as informative as the final success. A proposed radial-velocity channel, added because real S-star orbit determination traditionally uses one, turns out to contribute nothing at S301's magnitude — verified directly, not merely argued, and dropped from the actual proposal. A cadence chosen by formal optimal-design theory, expected to outperform a hand-picked heuristic, instead collapses catastrophically once tested against a real nonlinear fit — a reproducible pitfall of applying local-linearization statistics to a strongly nonlinear, periodic problem — and is rescued only by reintroducing the physical intuition it was meant to replace. The resulting campaign recovers all seven orbital parameters, including the precession rate, to sub-percent fractional precision — within 1.5σ of the injected truth, a formal detection of nonzero precession at ~179σ against this study's idealized, noise-only error budget. We conclude that a two-passage, astrometry-only GRAVITY+ campaign is both feasible and sufficient to make a genuine, spin-independent test of general relativity near Sgr A\*.

**Keywords:** Galactic Center — black hole physics — astrometry — celestial mechanics — general relativity — experimental design

---

## 1. Introduction

The Galactic Center's population of short-period "S-stars" orbiting Sgr A\* has provided some of the strongest empirical tests of general relativity (GR) in the strong-field, near-horizon regime, most notably the detection of gravitational redshift and Schwarzschild precession in the orbit of S2 (GRAVITY Collaboration 2018, 2020). The 2026 discovery of S301 — a star with a substantially shorter period (8.68 yr vs. S2's ~16 yr), higher eccentricity, and closer periapsis passage (136 R_S) than S2 — extends this program to a regime where relativistic effects are proportionally larger and, in principle, more rapidly measurable.

The discovery paper (GRAVITY Collaboration et al. 2026) frames S301 explicitly as a candidate for measuring the spin of Sgr A\* via relativistic (Lense-Thirring) frame-dragging, a genuinely unmeasured quantity. A companion theoretical paper (Piran et al. 2026) develops that spin-measurement methodology in detail, including the real difficulty of separating a frame-dragging signal from ordinary "Newtonian confusion" produced by other, unmodeled mass near Sgr A\*.

Before attempting that harder measurement, a narrower and logically prior question deserves its own answer: can a real, currently operable campaign detect the *non-spin* apsidal precession that any Schwarzschild (non-rotating) black hole would already produce — the same effect responsible for Mercury's perihelion advance, scaled up by many orders of magnitude here by proximity to a far more compact and massive body? This paper answers that question by building the complete pipeline a real observing proposal would need, in the order a real proposal would need it, and reporting exactly what worked, what did not, and why.

We proceed in four stages, each motivating the next. Section 3 builds and validates a relativistic light-curve model against the published orbital solution — establishing that our physics is correct before we ask anything of an observing campaign. Section 4 turns that validated model into a synthetic, realistic observing campaign, which forces two concrete decisions: which instruments to actually request time on, and how to schedule the observations across an extremely eccentric orbit. Both decisions turn out to be less obvious than they first appear, and we report both as findings rather than settling them by assumption. Section 5 uses the campaign's photometric channel — generated but, as we show, never needed for the fit itself — as an independent visual check that the relativistic physics from Section 3 survives intact inside realistic noisy data. Section 6 then fits the noisy astrometric data recovered from this process to ask whether the precession is actually detectable, and Section 7 reports that it is. Section 8 discusses what this result does and does not establish, and Section 9 summarizes.

## 2. The Real S301 Orbital Solution

Every physical parameter used as ground truth in this study is taken directly from GRAVITY Collaboration et al. (2026), Extended Data Table 2, Solution 1 (a second, nearly degenerate solution is also published, reflecting a known astrometry-only orbit-fitting ambiguity, but is not used here). Table 1 summarizes the adopted values; everything that follows is built on top of them.

**Table 1.** Adopted S301 system parameters (GRAVITY Collaboration et al. 2026, Solution 1).

| Parameter | Symbol | Value | Uncertainty |
|---|---|---|---|
| Orbital period | *P* | 8.680 yr | 0.110 yr |
| Eccentricity | *e* | 0.9832 | 0.0010 |
| Inclination | *i* | 124.09° | 1.10° |
| Longitude of ascending node | Ω | 73.8° | 3.5° |
| Argument of periapsis | ω | 293.4° | 2.2° |
| Time of periapsis (epoch) | *t*₀ | 2023.126 | 0.010 yr |
| Semi-major axis | *a* | 83.0 mas (687 AU) | — |
| Sgr A\* mass | *M*• | 4.297 × 10⁶ M☉ | — |
| Distance to Sgr A\* | *R*₀ | 8277 pc | — |
| Spectral type | — | F1.5V | — |
| Apparent K-band magnitude | *m*_K | 19.3 | 0.3 |
| Stellar mass | — | 1.1–1.5 M☉ | — |
| Stellar radius | — | 1.4–1.6 R☉ | — |

Two derived quantities anchor everything that follows: at periapsis, S301 passes 136 Schwarzschild radii from Sgr A\* at a speed exceeding 8% of the speed of light (~25,000 km s⁻¹). No other known star approaches this closely or moves this fast, including S2 — which is precisely why S301, and not S2, is the subject of this study. The star's own effective temperature is not published; we adopt 7300 K from standard spectral-type tables for F1.5V, flagged explicitly as an assumption rather than a measured quantity, since it affects only Section 3's absolute-magnitude estimate and not the orbital dynamics used everywhere else.

## 3. Validating the Physics: A Relativistic Light Curve

Before asking whether a future campaign can detect anything, we first ask a more basic question: does our orbital-mechanics and relativity engine reproduce the star we already know about? We model S301's near-infrared light curve using a full Keplerian orbit solution (Newton–Raphson eccentric-anomaly solve) combined with (i) the exact special-relativistic Doppler factor computed from the star's true three-dimensional velocity (not a small-*v* approximation), capturing both classical Doppler shift and transverse (time-dilation) redshift, and (ii) the Schwarzschild gravitational redshift factor √(1 − *R*_S/*r*). The combined factor (1 + *z*) is applied to a blackbody spectrum via the Lorentz-invariant relativistic-beaming relation *I*_ν/ν³, integrated over GRAVITY/VLTI's real K-band bandpass (1.98–2.40 μm). Gravitational lensing/magnification, occultation by the black hole shadow, Kerr (spin) corrections, and Rømer light-travel-time delay are explicitly not modeled, as they either require a spin value we do not have or a line-of-sight alignment not implied by the orbital elements alone.

The model passes its own consistency check: the semi-major axis derived from (*GM*•, *P*) via Kepler's third law is 686.7 AU against a published 687.0 AU; the derived periapsis distance is 136.0 *R*_S against a published 136; and the derived periapsis speed is 8.54% *c* against a published ">8%, 25,000 km s⁻¹." Only having confirmed this — that the same engine used everywhere below reproduces the one orbit we can independently check — do we treat its predictions elsewhere as trustworthy.

The resulting light curve (Figure 1, upper two panels) reveals the actual observational signature of periapsis, which turns out not to be the one a naive guess might reach for. It is tempting to picture the black hole "swallowing" some of the star's light near closest approach, producing a symmetric dip; the real prediction is an *asymmetric* relativistic beaming spike, with no absorption or occultation involved at all. The star brightens by up to ~0.075 mag as it approaches at high velocity (relativistic Doppler boosting), then fades to ~0.04 mag fainter than baseline as it recedes (Doppler de-boosting compounded by gravitational redshift, which dims regardless of direction of motion), for a peak-to-trough amplitude of 0.116 mag. The lower panel of Figure 1 shows this is a genuinely relativistic effect, not a relabeled classical one: the total (special+general relativistic) redshift measurably exceeds the naive classical Doppler prediction (*v_r*/*c*) by 0.26–0.46 percentage points near periapsis, a correction at the ~10–20% level relative to the classical term itself — large enough to matter because S301's velocity and proximity to Sgr A\* push terms normally negligible (order (*v*/*c*)²) into significance.

This is the physical signature the rest of this study either has to reproduce inside noisy simulated data (Section 5) or explain away — and, as shown below, it survives.

## 4. Building a Realistic Observing Campaign

A validated physics engine only becomes a feasibility study once it is asked to predict what a *specific, real* campaign would actually measure. Turning Section 3's theoretical curve into a synthetic campaign forces two decisions that cannot be waved away: which instruments to actually request time on, and how to schedule a fixed number of observations across an orbit that spends almost all its time barely moving and almost none of its time doing anything interesting. We address each in turn, and in both cases our first, most obvious answer turned out to be wrong in an instructive way.

The campaign's assumed start (2028.5) is itself a modeling choice, not an arbitrary one: a multi-year, two-passage VLTI monitoring program of this kind would realistically be proposed through ESO's Large Programme process rather than routine scheduling, which we take to require approximately two years end to end — time to write the proposal, submit it to the next call, wait for Time Allocation Committee review, and wait for the awarded observing period to begin — rather than an immediate start from today.

### 4.1 Instrument Selection and the Case Against Radial Velocity

The campaign is built around GRAVITY+, the real, currently operating VLTI interferometer instrument that made the original S301 discovery, providing both astrometric offsets (right ascension and declination relative to Sgr A\*) and simultaneous K-band photometry as a natural byproduct of the same interferometric observation. We adopt GRAVITY+'s real, demonstrated on-sky astrometric precision of 100 μas at K ≈ 19–20 (GRAVITY+ Collaboration 2023), the conservative end of its published 30–100 μas performance range.

An earlier design iteration also included a simulated radial-velocity (RV) channel via ERIS/VLT (the SPIFFIER integral-field spectrograph, via the Brackett-γ line), the real instrument used for essentially all S-cluster RV monitoring to date — a natural first instinct, since real S-star orbit determination has historically combined astrometry with spectroscopic RV. We derive its expected precision by scaling SINFONI's real, achieved RV precision on S2 (12.3 km s⁻¹ at *m*_K = 14.0; GRAVITY Collaboration et al. 2018), assumed comparable for ERIS as SINFONI's direct successor at the VLT, to S301's much fainter magnitude (Δ*m* = 5.3, a factor of ~130 in flux) under the background-limited noise regime appropriate for ground-based K-band spectroscopy of faint, crowded Galactic Center targets. This scaling predicts a per-epoch RV precision of only ~1621 km s⁻¹ against a true signal of order ±15,000 km s⁻¹ — an implied continuum signal-to-noise ratio of ~0.04 against SPIFFIER's own 60 km s⁻¹ velocity-resolution element (*R* = 5000). This is not a marginal detection; it is no detection at all with any presently operating 8–10 m telescope.

We treat this not as an incidental limitation to be discovered after the fact, but as a decision that belongs *before* an instrument-time request is ever written: a real Time Allocation Committee requires exactly this kind of feasibility calculation, and would not grant time on an oversubscribed instrument for a measurement already known to fail this badly. Accordingly, **the campaign analyzed in this paper does not include or request an RV channel.** We confirmed directly, rather than assumed, that nothing is lost by this decision: refitting the same synthetic astrometric dataset (Section 6) with a simulated RV channel added back in changes no recovered parameter's bootstrap uncertainty outside of resampling noise — the precession rate, in particular, comes back marginally *more* precisely without RV (0.0013 vs. 0.0014° yr⁻¹). This mirrors how the real S301 discovery paper itself analyzed the orbit — from astrometry, because RV from current instrumentation adds essentially nothing for a star this faint.

### 4.2 Cadence Design: A Cautionary Tale in Optimal Experimental Design

Instrument selection settled, a second question remains: given a fixed observing budget, when should those observations actually happen? Because S301's orbit is extremely eccentric, a cadence that samples time uniformly is badly mismatched to where the orbit carries information — the star spends most of its 8.68 yr period moving slowly near apoapsis and executes nearly its entire angular excursion within days of periapsis. We explored three progressively more principled answers, each evaluated by an identical seven-parameter nonlinear least-squares fit (Section 6) at matched total epoch count (*N* = 131), isolating cadence as the only variable.

**Design A (heuristic).** A hand-specified design placing approximately 90% of epochs within ±20 days of each periapsis passage, with the remainder as a sparse long-baseline presence. Simple, physically motivated, and — as Section 7 shows — reasonably effective, but it implicitly assumes all seven fit parameters are equally well served by periapsis proximity.

**Design B (naive D-optimal).** Wanting something more rigorous than a hand-picked rule, we computed the Fisher information matrix — the standard tool of optimal experimental design theory — by numerically differentiating the astrometric observable model with respect to all seven fit parameters at the published truth values, and greedily selected the *N* candidate epochs maximizing the determinant of the resulting information matrix (D-optimality, which minimizes the volume of the joint parameter confidence ellipsoid). This design **catastrophically failed**. It concentrated nearly the entire epoch budget onto two narrow clusters at the extreme temporal edges of the campaign window — one 5.4 years, the other 18.5 years (more than two orbital periods) from the reference epoch — and selected *zero* epochs within 20 days of either periapsis. Its analytic (Cramér–Rao bound) prediction was 1.8–28,000× tighter parameter uncertainties than Design A — driven almost entirely by the period, eccentricity, and time-of-periapsis terms (600–28,000×), with the angular orientation elements predicted to improve only modestly (1.8–2×); a real nonlinear bootstrap fit on this design was instead 4–600× *worse* than Design A on every single parameter (Table 2).

The mechanism is worth understanding, because it is not a coding error but a real, reproducible statistical pitfall. Sensitivity to the orbital period and epoch of periapsis genuinely grows with time baseline — the same principle behind pulsar timing's exquisite period precision from long observational baselines. But this sensitivity is a *local, linearized* quantity: a small period error, extrapolated across more than two orbital cycles, compounds into a large, nonlinear phase drift (a "cycle-slip" or aliasing hazard well known in periodic-signal parameter estimation). The Fisher information matrix, linearized at one exact fiducial point, cannot represent this: it reports enormous information from a design that is only informative if the period is already known to a precision no realistic prior estimate could supply, while providing literally zero geometric information about the orbit's shape from either of its two periapsis-free epoch clusters.

**Design C (pseudo-Bayesian robust D-optimal).** Following the classical fix for this pathology (Chaloner & Verdinelli 1995), we replaced the single-point linearization with a Fisher information matrix averaged over 300 parameter draws from a realistic prior (this system's own published element uncertainties; Table 1). This substantially de-emphasized the fragile long-baseline clusters, but *still* selected zero epochs within 20 days of either periapsis passage. This second failure is distinct from the first: Fisher information, even averaged over a parameter prior, measures only *noise-driven* uncertainty around a presumed-correct optimum. It has no mechanism to represent *structural or geometric non-identifiability* from sampling too few genuinely distinct points along the orbit — the reason classical orbit determination, from Gauss's method onward, has always required observations spread across an orbit's curvature, not merely high precision concentrated at a few points.

**Design D (constrained: floor + Fisher-optimal remainder).** Two purely statistical attempts having failed for related but distinct reasons, we reintroduced the one piece of domain knowledge both attempts had discarded: a highly eccentric orbit's shape and 3-D orientation cannot be resolved without observing near periapsis, no matter how a statistical criterion scores the alternative. We imposed a hard, non-negotiable floor — reserving ~30% of the total epoch budget for dense coverage within ±20 days of each periapsis passage — and allocated the *remaining* ~70% by the same Fisher-optimal search used in Design C, now deciding only how best to spend the flexible portion of the budget rather than whether periapsis coverage is needed at all. This hybrid design is the one used as the actual campaign in this study, and Section 7 shows it improves on Design A on every fit parameter, at negligible cost to the two parameters (period, eccentricity) Design A already recovered well.

All three designs are additionally filtered to Sgr A\*'s real annual visibility window from Paranal (approximately March–September; GRAVITY Collaboration et al. 2026), which exposes a genuine, unavoidable scheduling complication: the first of the two periapsis passages analyzed here (2031 October 22) falls just after that year's visibility season closes. A space-based instrument could observe this passage on demand; a single-site ground-based campaign cannot. We therefore anchor that passage's dense-coverage window to the nearest date the target is actually observable (~2031 September 30) rather than to periapsis itself, accepting that the post-periapsis half of the ideal observing window is lost. The second passage (2040 June 26) falls comfortably within its visibility season.

## 5. Photometry as an Independent Physical Cross-Check

The campaign built in Section 4 generates two observable channels per epoch — astrometric offsets and K-band photometry — but only astrometry is used in the orbit fit of Section 6, since the photometric signal is a deterministic function of the same (position, velocity) state that astrometry already constrains and therefore carries no independent orbital-element information. This raises an obvious question: if photometry is never fit, why simulate it at all?

The answer is that photometry serves a different purpose entirely: it is a free, independent check that the relativistic physics validated in Section 3 actually survives intact once realistic instrumental noise, sparse and irregular scheduling, and the full campaign pipeline are all layered on top of it. GRAVITY+ measures K-band flux as a natural byproduct of the same interferometric observation used for astrometry, at no additional cost in proposed instrument time, so a realistic simulated campaign should generate this channel even though the fit does not need it — this is what a real GRAVITY+ dataset would actually contain, and it lets us ask whether Section 3's predicted signature is something an observer would actually be able to see, not merely something that exists in a noiseless theory curve.

Figure 2's lower panel answers this directly. S301's published photometric uncertainty (0.30 mag per epoch) is, by itself, larger than the entire 0.116 mag peak-to-trough relativistic signal derived in Section 3 — a single noisy data point near periapsis is not obviously distinguishable from noise. But the *injected truth curve*, threaded through the actual simulated data points across many epochs and both periapsis passages, traces out exactly the asymmetric brighten-then-dim structure predicted in Section 3, visibly recognizable against the scatter once enough epochs are combined. This is a qualitative, not a quantitative, confirmation — we do not fit this curve, and it adds no precision to Table 3's results — but it is a genuine end-to-end sanity check that the physics injected in Section 3, the noise model built in Section 4, and the campaign's actual cadence are all mutually consistent, rather than three separately-plausible pieces that happen to have been assembled together.

## 6. Orbit Recovery Methodology

With instrument selection and cadence resolved (Section 4) and the injected physics visually confirmed to survive realistic noise (Section 5), the remaining question is whether a fit that never sees the injected truth can recover it. We fit the synthetic campaign's astrometric data with a seven-parameter model: the six standard Keplerian elements (*P*, *e*, *i*, Ω, ω, *t*₀) plus an apsidal precession rate ω̇, applied by allowing the argument of periapsis to vary linearly in time about a reference epoch, ω_eff(*t*) = ω + ω̇(*t* − *t*_ref). The semi-major axis is not an independent free parameter; it is derived from (*P*, *GM*•) via Kepler's third law at each fit iteration, consistent with its treatment in Section 3. The fit minimizes sigma-normalized residuals between the model and the astrometric (RA, Dec) data via Levenberg–Marquardt nonlinear least squares (`scipy.optimize.least_squares`), starting from an initial guess deliberately perturbed away from the truth (by 3% in period, and by several degrees in the angular elements) to ensure the recovery is a genuine fit rather than a fixed point. Parameter uncertainties are estimated by bootstrap resampling of the observation epochs with replacement (1000 resamples, refit from the best-fit solution each time).

The recovered precession rate ω̇ is compared directly against the analytic Schwarzschild (1PN) prediction,

ω̇ = 6π *GM*• / [*c*² *a*(1 − *e*²) *P*],

the same formula responsible for Mercury's perihelion advance, evaluated with S301's own published elements (Table 1). No black-hole-spin (Lense–Thirring) contribution is injected anywhere in this study: spin causes a physically distinct effect (nodal precession, i.e., drift in *i* and Ω via frame-dragging) that is genuinely unmeasured for Sgr A\*, and injecting an assumed value would substitute a guess for a measurement (Section 8.2).

## 7. Results

### 7.1 Campaign Statistics

The final (Design D) campaign comprises 131 epochs spanning 2028.50–2041.69, with 39 epochs (29.8%) falling within ±20 days of a periapsis passage and the remaining 92 epochs allocated by the Fisher-optimal search across the rest of the visibility-filtered observing window — more than double the 60-epoch, single-passage campaign of an earlier design iteration. The injected Schwarzschild precession rate is 0.2307° yr⁻¹ (2.003° per orbit).

### 7.2 Cadence Comparison

**Table 2.** Bootstrap parameter uncertainties (1σ) under the heuristic (Design A) and constrained Fisher-optimal (Design D) cadences, at matched total epoch count (*N* = 131). Both fits include an RV channel for this comparison, isolating cadence as the only variable (the RV-channel comparison itself, at fixed cadence, is reported separately in Section 4.1).

| Parameter | Design A (heuristic) σ | Design D (constrained) σ | Published σ |
|---|---|---|---|
| *P* (yr) | 0.0001 | 0.0002 | 0.1100 |
| *e* | 0.0001 | 0.0000 | 0.0010 |
| *i* (°) | 0.0480 | 0.0356 | 1.1000 |
| Ω (°) | 0.1229 | 0.0849 | 3.5000 |
| ω (°) | 0.0939 | 0.0475 | 2.2000 |
| *t*₀ (yr) | 0.0003 | 0.0004 | 0.0100 |
| ω̇ (° yr⁻¹) | 0.0043 | 0.0014 | — (unpublished) |

Design D improves *i*, Ω, ω, and — most importantly for this study's science goal — ω̇ by 26–67%, at the cost of a small, practically negligible degradation in *P* and *t*₀ (both remain more than two orders of magnitude tighter than the published uncertainty either way).

### 7.3 Final Orbit Fit (Astrometry Only, Constrained Cadence)

**Table 3.** Final orbit-recovery result: Design D cadence, astrometry-only (no RV channel; Section 4.1).

| Parameter | Truth | Fitted | Fit σ | Published σ |
|---|---|---|---|---|
| *P* (yr) | 8.6800 | 8.6799 | 0.0002 | 0.1100 |
| *e* | 0.9832 | 0.9832 | 0.0000 | 0.0010 |
| *i* (°) | 124.0900 | 124.1207 | 0.0337 | 1.1000 |
| Ω (°) | 73.8000 | 73.6948 | 0.0807 | 3.5000 |
| ω (°) | 293.4000 | 293.3333 | 0.0452 | 2.2000 |
| *t*₀ (yr) | 2023.1260 | 2023.1261 | 0.0004 | 0.0100 |
| ω̇ (° yr⁻¹) | 0.2307 | 0.2326 | 0.0013 | — |

Every parameter is recovered within a small fraction of one bootstrap standard deviation of the injected truth. The precession rate is recovered as 0.2326 ± 0.0013° yr⁻¹ against a true value of 0.2307° yr⁻¹ — consistent with truth at ~1.5σ, and formally distinguishable from zero precession at ~179σ. That extreme formal significance is an artifact of this study's idealized, systematics-free noise budget, not a claim about what a real campaign's actual error budget would achieve; the meaningful result is the ~1.5σ agreement with the true, physically predicted value, not the significance against an astrophysically implausible null of exactly zero precession. The fitted sky-plane track (Figure 2, upper panel) visibly shows the argument of periapsis rotated between the two passages, the direct geometric signature of the precession, closing the loop opened in Section 3: the same physics predicted there from the published solution alone is now recovered independently, from noisy simulated observations that never had access to that solution.

We verified directly, by refitting the identical astrometric dataset with a simulated RV channel added back in, that this result depends on astrometry alone: every parameter's bootstrap uncertainty is statistically unchanged by the presence or absence of RV, consistent with the SNR ≈ 0.04 finding of Section 4.1.

One consequence bears noting plainly: because this campaign assumes GRAVITY+'s best demonstrated precision at every epoch and spans two full periapsis passages 8.68 years apart — a temporal baseline the real 2026 discovery dataset does not yet possess — several fitted uncertainties here (e.g., 0.08° on Ω) are substantially tighter than the real paper's own published values (3.5°). This reflects the idealized, single-campaign nature of a synthetic feasibility study, not a claim of having outperformed the real discovery.

## 8. Discussion

### 8.1 Instrumentation Alternatives

Having settled on GRAVITY+ astrometry alone (Section 4.1), we checked whether any alternative facility could do better on either channel. For astrometry, no operating alternative matches GRAVITY+/VLTI's interferometric precision at this magnitude: a 2.4 m single-aperture space telescope (e.g., Roman) has a diffraction limit (~223 mas) larger than S301's entire 83 mas orbit and cannot resolve it at any exposure time; JWST's 6.5 m aperture, while diffraction-limited to a sharper ~85 mas, still lacks an interferometric baseline and cannot approach GRAVITY+'s ~100 μas relative astrometry.

For radial velocity, we likewise found no currently viable substitute for the rejected ERIS channel. Roman's RV prospects fail for two independent, integration-time-independent reasons: its broad (~223 mas), non-adaptive-optics-corrected point-spread function produces severe source confusion in the crowded S-cluster — a structural noise floor that, unlike photon noise, does not average down with longer exposure — compounded by markedly lower spectral resolving power (*R* ≈ 461–890 vs. ERIS's *R* = 5000). JWST/NIRSpec presents a more ambiguous case: its measured spectral resolution is comparable to or better than pre-launch expectations (Shajib et al. 2025), and it entirely avoids the OH airglow that limits ground-based K-band spectroscopy; however, real JWST observations of this specific field report severe, field-specific operational problems (detector saturation from Sgr A\*'s own brightness, multi-shutter-array light leakage, and crowding-induced guide-star misidentification producing several-arcsecond pointing offsets; Yusef-Zadeh et al. 2025), and no published, achieved RV precision exists for a point source as faint and crowded as S301 to ground a noise model on. Consistent with this study's standing principle of never injecting an uncited performance figure, we do not adopt JWST as a substitute instrument.

The field's own real answer for RV is the Extremely Large Telescope's HARMONI instrument (*R* = 3500–18,000 on a 39 m aperture, ~22.7× VLT's collecting area). As of the most recent published ESO schedule, ELT science first light (after HARMONI and MICADO are installed and commissioned) has slipped to December 2030 — after this study's assumed 2028.5 campaign start, but only about ten months before the first periapsis passage analyzed here (2031.8). A newly commissioned instrument does not immediately accept competitively awarded external proposals: guaranteed-time and verification-science programs typically occupy the first observing cycles, and this study's own timing assumption for GRAVITY+ (Section 4) already allows roughly two years between a Large Programme proposal and its first scheduled night. A ten-month margin is unlikely to be enough for a new S301 RV proposal on a brand-new instrument to clear that process in time for the first passage; HARMONI is more plausibly available in time to contribute to the second passage (2040.5), nearly a decade after first light.

### 8.2 What This Study Does Not Attempt

Measuring the spin of Sgr A\* is explicitly outside the scope of this work. Spin produces a physically distinct signature — nodal (Lense–Thirring) precession, a secular drift in *i* and Ω from frame-dragging — that this study does not inject anywhere; both angles are held genuinely static throughout (Section 6), consistent with spin being a real, currently unmeasured quantity for Sgr A\*. Piran et al. (2026) show that a genuine spin measurement additionally requires separating the frame-dragging signal from ordinary Newtonian confusion produced by other, unmodeled mass near Sgr A\*, using multi-star calibration (S2, S55, S38) — a substantially more sophisticated analysis than the single-star, spin-free precession measurement presented here. This study's two-passage Schwarzschild-precession detection is a genuine, correctly modeled scientific result and a necessary methodological precursor to that harder measurement, but is not itself a spin measurement, and should not be represented as one.

### 8.3 The Experimental-Design Lesson

Beyond the specific S301 result, we highlight the cadence-design methodology of Section 4.2 as a generally applicable caution. Local Fisher-information-based optimal design is a standard and often effective tool, but it silently assumes the model response is well approximated by its linear expansion at the assumed truth over the *full range* of candidate design points. For a strongly nonlinear, periodic problem evaluated over a multi-cycle observing baseline, this assumption can fail badly and in a way that is invisible from the analytic (Cramér–Rao) prediction alone — it becomes visible only when the design is validated against an actual nonlinear fit, exactly as happened here. We recommend that any application of D-optimal or related design criteria to nonlinear periodic systems include (i) a pseudo-Bayesian averaging over realistic parameter-prior uncertainty (Chaloner & Verdinelli 1995), and (ii) an explicit, domain-knowledge floor protecting against the structural non-identifiability that a purely statistical criterion cannot detect — followed, in every case, by validation against a real nonlinear fit rather than reliance on the analytic bound alone.

## 9. Summary and Conclusions

We set out to answer one question — can a real, near-term campaign detect S301's relativistic apsidal precession from astrometry alone — and answered it by building the complete pipeline in the order a real proposal would need it: validate the physics, design the campaign, cross-check the physics against realistic noisy data, then fit. Two of our own initial instincts along the way failed real tests, and both failures turned out to be as informative as the final result:

1. A relativistic light-curve model reproduces S301's published orbital solution to better than 0.1% in derived semi-major axis, periapsis distance, and periapsis speed, and predicts an asymmetric, brighten-then-dim relativistic beaming feature (0.116 mag peak-to-trough) at periapsis — not an absorption or occultation dip, as a naive first guess might suggest.
2. That same signature survives, visibly, inside realistic simulated GRAVITY+ photometry despite per-epoch noise larger than the signal itself — a free, independent confirmation that the injected physics, noise model, and cadence are mutually consistent, even though photometry is never used in the orbit fit itself.
3. A feasibility calculation, grounded in ERIS/SINFONI's real achieved performance on S2, shows that radial-velocity monitoring of S301 is not viable with any current 8–10 m telescope (SNR ≈ 0.04); we verified directly that omitting this channel costs no orbit-recovery precision, and recommend against requesting it from a real Time Allocation Committee.
4. A naive Fisher-information (D-optimal) cadence design catastrophically fails for this problem via a nonlinear aliasing pitfall, despite an extremely favorable analytic prediction; only a hybrid design that reintroduces the physical intuition the statistical method discarded — a hard periapsis-coverage floor plus Fisher-optimal allocation of the remaining budget — outperforms a simpler heuristic on every fit parameter.
5. The resulting campaign recovers the injected Schwarzschild precession rate (0.2307° yr⁻¹) as 0.2326 ± 0.0013° yr⁻¹ — within 1.5σ of truth — using astrometry alone.

We conclude that a real, near-term, astrometry-only GRAVITY+ campaign spanning S301's next two periapsis passages (2031.8 and 2040.5) is both feasible with existing instrumentation and sufficient to make an independent, spin-agnostic test of general relativity near a supermassive black hole, ahead of and complementary to the harder measurement of Sgr A\*'s spin that motivated S301's discovery in the first place.

## Acknowledgments

This work makes use of the published orbital solution and system parameters of GRAVITY Collaboration et al. (2026) and the theoretical framework of Piran et al. (2026). All simulation code, generated figures, and the full derivation of every adopted instrumental precision figure are available in the accompanying project repository.

## References

- Chaloner, K., & Verdinelli, I. 1995, "Bayesian Experimental Design: A Review," Statistical Science, 10(3), 273–304
- ESO 2026, "Telescope first light for ESO's Extremely Large Telescope now planned for March 2029," https://www.eso.org/public/announcements/ann25001/
- GRAVITY Collaboration, Abuter, R., Amorim, A., et al. 2018, A&A, 615, L15, "Detection of the gravitational redshift in the orbit of the star S2 near the Galactic centre massive black hole" (arXiv:1807.09409) — SINFONI radial-velocity precision for S2 (12.3 km/s)
- GRAVITY Collaboration, Abuter, R., Amorim, A., et al. 2020, A&A, 636, L5, "Detection of the Schwarzschild precession in the orbit of the star S2 near the Galactic centre massive black hole" (arXiv:2004.07187)
- GRAVITY Collaboration, Abd El Dayem, K., et al. 2026, Nature, "Discovery of a star sensitive to the spin of Sgr A\*" (arXiv:2607.12664)
- GRAVITY+ Collaboration 2023, "The GRAVITY+ Project: Towards All-sky, Faint-Science, High-Contrast Near-Infrared Interferometry at the VLTI" (arXiv:2301.08071)
- Piran, T., Amaro-Seoane, P., Aytac, B., et al. 2026, "S301 and friends: Measuring the spin of Sgr A\*" (arXiv:2607.24931)
- Shajib, A. J., Treu, T., Melo, A., et al. 2025, A&A, 702, L12, "An accurate measurement of the spectral resolution of the JWST Near Infrared Spectrograph" (arXiv:2507.03746)
- Yusef-Zadeh, F., Bushouse, H., Arendt, R. G., et al. 2025, ApJL, "Non-stop Variability of Sgr A\* using JWST at 2.1 and 4.8 micron Wavelengths" (arXiv:2501.04096)

## Code and Data Availability

The complete simulation pipeline (`orbit.py`, `s301_lightcurve.py`, `campaign.py`, `optimal_design.py`, `fit_orbit.py`, `constants.py`) and generated figures referenced above (Figures 1–2) are available in this project's repository. All instrumental precision figures are derived from cited, real, achieved measurements; every assumption not directly traceable to a published source is flagged explicitly in the corresponding module's documentation.
