"""
Gravitational-lensing astrometric systematic: how much Sgr A* bends
S301's light, epoch by epoch, and what ignoring it does to the fit.

Geometry: Sgr A* (mass M) is the lens at the origin; the observer is
along +Z (orbit.py convention: +Z toward the observer). When S301 is on
the FAR side of the black hole (sky_z < 0) its light passes the lens on
the way to us and the primary image is displaced outward from the lens
direction; on the near side (sky_z > 0) the light never passes the lens
and the deflection is negligible. For a point lens the primary image of
a source at angular separation theta_s sits at

    theta_+ = [theta_s + sqrt(theta_s^2 + 4 theta_E^2)] / 2,
    theta_E^2 = (4 G M / c^2) * D_LS / (D_L * D_S),

with D_L = R0, D_LS = -sky_z (line-of-sight depth behind the lens) and
D_S = D_L + D_LS. The shift theta_+ - theta_s ~ theta_E^2 / theta_s for
theta_s >> theta_E. At S301's periapsis (D_LS up to ~136 R_S, theta_s ~
1.4 mas) theta_E ~ 0.17 mas, so the shift is of order 20 uas at the few
far-side epochs nearest periapsis and far smaller elsewhere. The
secondary image is fainter by (theta_-/theta_+)^2 ~ (theta_E/theta_s)^4
~ 1e-4 and is ignored.

Method (same inject-and-refit pattern as the other systematics): apply
the exact point-lens shift to the noise-free astrometry, refit with the
lens-unaware model, and report the per-parameter bias against Table 3's
precision; also report the shift distribution over the campaign epochs.
Results go to results/lensing_error.json.
"""

import numpy as np

import campaign
import constants as k
import fit_orbit
import results_io

RAD_TO_MAS = np.degrees(1.0) * 3600.0 * 1000.0


def lensing_shift_mas(state):
    """Astrometric shift (mas) of the primary image, per epoch, from the
    exact point-lens equation; zero for epochs on the near side."""
    theta_s_rad = np.hypot(state["sky_x"], state["sky_y"]) / k.D_OBS
    d_ls = np.clip(-state["sky_z"], 0.0, None)          # depth behind the lens, m
    d_s = k.D_OBS + d_ls
    theta_e_sq = (4.0 * k.GM_BH / k.c ** 2) * d_ls / (k.D_OBS * d_s)
    theta_plus = 0.5 * (theta_s_rad + np.sqrt(theta_s_rad ** 2 + 4.0 * theta_e_sq))
    shift_rad = np.where(d_ls > 0, theta_plus - theta_s_rad, 0.0)
    return shift_rad * RAD_TO_MAS, np.sqrt(theta_e_sq) * RAD_TO_MAS, theta_s_rad * RAD_TO_MAS


def apply_shift(ra_mas, dec_mas, shift_mas):
    """Displace each position radially away from the lens by shift_mas."""
    sep = np.hypot(ra_mas, dec_mas)
    return ra_mas * (1 + shift_mas / sep), dec_mas * (1 + shift_mas / sep)


def main():
    """Shift distribution over the campaign, then the lens-unaware fit bias."""
    data = fit_orbit.load_observations()
    epochs = data["epoch_yr"]
    state = campaign.true_state(epochs)
    shift, theta_e, theta_s = lensing_shift_mas(state)
    far_side = state["sky_z"] < 0
    peri1 = k.NEXT_PERIAPSIS_YR
    print(f"Epochs on the far side of Sgr A*: {far_side.sum()} of {len(epochs)}")
    j_max = int(np.argmax(shift))
    print(f"Max lensing shift: {shift[j_max] * 1000:.1f} uas at epoch {epochs[j_max]:.3f} "
          f"({(epochs[j_max] - peri1) * 365.25:+.0f} d from periapsis 1); "
          f"theta_E there {theta_e[j_max] * 1000:.0f} uas, separation {theta_s[j_max]:.2f} mas")
    print(f"Epochs with shift > 10 uas: {(shift > 0.010).sum()}; > 1 uas: {(shift > 0.001).sum()}; "
          f"RMS over all epochs: {np.sqrt(np.mean(shift ** 2)) * 1000:.2f} uas")

    base = results_io.read_results("fit_orbit")
    sigma_ref = np.array([base[f"sigma_{name}"] for name in fit_orbit.PARAM_NAMES])
    idx_od = fit_orbit.PARAM_NAMES.index("omega_dot_deg_yr")

    ra_true, dec_true, _, _ = campaign.true_observables(epochs)
    ra_lensed, dec_lensed = apply_shift(ra_true, dec_true, shift)
    sig = np.full_like(epochs, k.GRAVITY_PLUS_ASTROMETRY_MAS)
    x0 = fit_orbit.TRUTH_VECTOR + fit_orbit.INITIAL_PERTURBATION
    unlensed_fit = fit_orbit.fit_once(epochs, ra_true, dec_true, sig, sig, x0=x0)
    lensed_fit = fit_orbit.fit_once(epochs, ra_lensed, dec_lensed, sig, sig, x0=unlensed_fit)
    bias = lensed_fit - unlensed_fit
    print("\n=== Lens-unaware 7-parameter fit on lensed, noise-free data ===")
    for name, b_val, s_val in zip(fit_orbit.PARAM_NAMES, bias, sigma_ref):
        print(f"  {name:<18} bias = {b_val:+.6f}  ({abs(b_val) / s_val:.2f} sigma of Table 3 precision)")

    results_io.write_results("lensing_error", {
        "n_far_side": int(far_side.sum()),
        "max_shift_uas": (shift.max() * 1000, ".1f"),
        "max_shift_epoch": (epochs[np.argmax(shift)], ".3f"),
        "max_shift_days_from_peri": ((epochs[np.argmax(shift)] - peri1) * 365.25, ".0f"),
        "theta_e_at_max_uas": (theta_e[np.argmax(shift)] * 1000, ".0f"),
        "n_shift_gt_10uas": int((shift > 0.010).sum()),
        "n_shift_gt_1uas": int((shift > 0.001).sum()),
        "rms_shift_uas": (np.sqrt(np.mean(shift ** 2)) * 1000, ".2f"),
        "omega_dot_bias": (bias[idx_od], ".5f"),
        "omega_dot_bias_nsigma": (abs(bias[idx_od]) / sigma_ref[idx_od], ".2f"),
        "omega_dot_bias_pct": (abs(bias[idx_od]) / fit_orbit.TRUE_OMEGA_DOT_DEG_YR * 100, ".3f"),
        "max_angle_bias_nsigma": (max(abs(bias[j]) / sigma_ref[j] for j in (2, 3, 4)), ".2f"),
    })


if __name__ == "__main__":
    main()
