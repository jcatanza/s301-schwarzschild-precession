"""
Animated visualization of S301's real Schwarzschild apsidal precession,
extrapolated over many more orbits than this project's actual analyzed
campaign (2028.5-2041.69, two passages).

This is NOT a new scientific result and injects nothing illustrative:
every physical quantity -- the published orbital elements, the real 1PN
precession rate (orbit.schwarzschild_precession_rate) -- is the same
ground truth used throughout the rest of this project (constants.TRUTH).
The ONLY thing that goes beyond the actual analyzed study is the time
baseline: N_ORBITS periods (~174 years) instead of the two real passages
the campaign covers, purely so the rosette pattern the precession traces
out is visible to the eye rather than the subtle ~2 deg/orbit rotation
seen between the two loops in fit_orbit.png. This extrapolation is
flagged explicitly on-screen throughout, and the two passages actually
covered by campaign.py's synthetic campaign are highlighted distinctly
from the purely-illustrative orbits around them.

Output: an MP4, intended as an arXiv ancillary file alongside
article.tex (see that file's Figure 2 caption for the reference).
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

# pylint: disable=wrong-import-position
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.collections import LineCollection

import constants as k
import orbit
# pylint: enable=wrong-import-position

N_ORBITS = 20
FRAMES_PER_ORBIT = 60
FPS = 30

# The two passages campaign.py's actual synthetic campaign covers (orbit
# indices 1 and 2 after the published t_peri, matching NEXT_PERIAPSIS_YR
# and NEXT_PERIAPSIS_YR + P_yr in campaign.py) are highlighted distinctly
# from the purely-illustrative orbits around them.
REAL_CAMPAIGN_ORBIT_INDICES = (1, 2)


# pylint: disable=duplicate-code
# Necessarily near-identical to s301_lightcurve.py's own
# build_truth_elements() -- both build the same constants.TRUTH orbit --
# differing only in that this one also returns period; no shared
# "truth-orbit" utility module in this project is worth introducing for
# this handful of lines, per fit_orbit.py's own duplicate-code note.
def build_truth_elements():
    """The real published orbit (constants.TRUTH), matching every other
    script in this project.

    t_peri must be the real absolute periapsis epoch (k.TRUTH["t_peri_yr"]),
    not 0.0: build_trajectory() below offsets its time array by that same
    epoch, and orbit_state_precessing() computes mean anomaly relative to
    elements.t_peri. Using 0.0 here while adding the real epoch to the time
    array introduces a constant ~29-degree mean-anomaly offset (t_peri_sec
    is not an integer number of periods), which desyncs the dense E=0
    sampling region from the actual periapsis -- the true closest approach
    then falls in a sparsely-sampled part of the eccentric-anomaly grid,
    understating how close the trajectory actually comes to Sgr A* by
    roughly an order of magnitude. Caught by checking why the rendered
    track appeared not to enclose the star; verified by comparing the
    grid's minimum r against the analytic periapsis distance a(1-e)."""
    period = k.TRUTH["P_yr"] * k.year
    sma = orbit.semi_major_axis_from_period(k.GM_BH, period)
    return orbit.OrbitalElements(
        t_peri=k.TRUTH["t_peri_yr"] * k.year, period=period, ecc=k.TRUTH["e"], sma=sma,
        i_deg=k.TRUTH["i_deg"], raan_deg=k.TRUTH["Omega_deg"], omega_deg=k.TRUTH["omega_deg"],
    ), sma, period
# pylint: enable=duplicate-code


def build_trajectory():
    """Sky-plane (RA, Dec) trajectory over N_ORBITS periods, densely
    sampled per orbit via the same eccentric-anomaly grid used elsewhere
    in this project (concentrates points near periapsis, avoiding the
    corner-cutting fit_orbit.py's build_plot_epoch_grid warns about).
    Returns (epoch_yr, ra_mas, dec_mas, orbit_index_per_point)."""
    elements, sma, period = build_truth_elements()
    omega_dot = orbit.schwarzschild_precession_rate(k.GM_BH, sma, k.TRUTH["e"], period)
    t_peri_sec = k.TRUTH["t_peri_yr"] * k.year

    ecc_anom = orbit.eccentric_anomaly_grid(FRAMES_PER_ORBIT)
    t_within_orbit = orbit.time_since_periapsis(ecc_anom, k.TRUTH["e"], period)

    t_sec_all = np.concatenate([t_within_orbit + t_peri_sec + n * period for n in range(N_ORBITS)])
    orbit_index = np.repeat(np.arange(N_ORBITS), FRAMES_PER_ORBIT)

    state = orbit.orbit_state_precessing(t_sec_all, elements, omega_dot, t_ref=t_peri_sec,
                                          grav_param=k.GM_BH)
    ra_mas, dec_mas = orbit.sky_offset_mas(state, k.D_OBS)
    epoch_yr = t_sec_all / k.year
    return epoch_yr, ra_mas, dec_mas, orbit_index, omega_dot


# pylint: disable=too-many-locals
# make_movie() sets up the whole animation (figure, artists, update
# closure) in one place; splitting it up would scatter tightly-coupled
# matplotlib setup across functions for no real gain in clarity.
def make_movie(outpath="output/s301_precession.mp4"):
    """Render and save the animation."""
    epoch_yr, ra_mas, dec_mas, orbit_index, omega_dot = build_trajectory()
    omega_dot_deg_yr = np.degrees(omega_dot) * k.year
    is_real_campaign = np.isin(orbit_index, REAL_CAMPAIGN_ORBIT_INDICES)

    fig, ax = plt.subplots(figsize=(8, 8))
    # Additive margin (a fraction of the actual data RANGE), not a multiplicative
    # scaling of the raw min/max -- this orbit's Dec max is small and near zero
    # (apoapsis-to-periapsis asymmetry means the top of the track sits just a
    # few mas above Sgr A*), so a multiplicative margin left almost no room at
    # the top, letting the legend sit directly on top of the periapsis
    # convergence point and making the real, smooth cusp there look "kinked."
    ra_range = np.max(ra_mas) - np.min(ra_mas)
    dec_range = np.max(dec_mas) - np.min(dec_mas)
    margin_frac = 0.08
    ax.set_xlim(np.max(ra_mas) + margin_frac * ra_range,
                np.min(ra_mas) - margin_frac * ra_range)   # inverted RA, like other figures
    ax.set_ylim(np.min(dec_mas) - margin_frac * dec_range,
                np.max(dec_mas) + margin_frac * dec_range)
    ax.set_aspect("equal")
    ax.set_xlabel("RA offset (mas)")
    ax.set_ylabel("Dec offset (mas)")
    # ms=16 here previously drew Sgr A* at a rendered diameter of order 2 mas at
    # this figure's scale -- LARGER than S301's ~1.4 mas periapsis distance, so
    # the star icon itself visually poked outside the orbit's narrow periapsis
    # neck, even though the underlying curve mathematically encloses the focus
    # exactly (a pure rotation about the origin cannot move the origin: see
    # orbit._rotate_to_sky_precessing). A small marker here is load-bearing,
    # not cosmetic -- the periapsis-zoom inset below is what actually lets a
    # viewer confirm enclosure at the relevant scale.
    ax.plot(0, 0, "k*", ms=6, zorder=5, label="Sgr A*")
    ax.grid(alpha=0.3)

    # Fixed-zoom inset on the focus: at the full rosette's scale, S301's
    # ~1.4 mas periapsis distance is far too small to visibly confirm the
    # orbit encloses Sgr A* (this is exactly why that enclosure looked wrong
    # before). This inset re-plots the same trail/comet data, clipped to a
    # small box around the origin, at a scale where the enclosure is
    # unambiguous.
    inset_half_width_mas = 4.0
    ax_inset = ax.inset_axes([0.03, 0.66, 0.31, 0.31])
    ax_inset.set_xlim(inset_half_width_mas, -inset_half_width_mas)
    ax_inset.set_ylim(-inset_half_width_mas, inset_half_width_mas)
    ax_inset.set_aspect("equal")
    ax_inset.plot(0, 0, "k*", ms=10, zorder=5)
    ax_inset.set_title("Zoom on Sgr A*", fontsize=8)
    ax_inset.tick_params(labelsize=6)
    ax_inset.grid(alpha=0.3)
    inset_trail_lc = LineCollection([], cmap="viridis",
                                     norm=plt.Normalize(epoch_yr.min(), epoch_yr.max()),
                                     linewidths=1.3, alpha=0.85, zorder=2)
    ax_inset.add_collection(inset_trail_lc)
    inset_comet, = ax_inset.plot([], [], "o", color="#d62728", ms=5, zorder=6)

    title = ax.set_title("")
    epoch_text = ax.text(0.02, 0.02, "", transform=ax.transAxes, fontsize=9, color="dimgray")
    banner = fig.text(
        0.5, 0.955,
        f"ILLUSTRATIVE: real Schwarzschild precession rate ({omega_dot_deg_yr:.4f} deg/yr) "
        f"extrapolated over {N_ORBITS} orbits (~{N_ORBITS * k.TRUTH['P_yr']:.0f} yr) -- "
        f"only orbits {REAL_CAMPAIGN_ORBIT_INDICES[0]}-{REAL_CAMPAIGN_ORBIT_INDICES[1]} "
        f"(highlighted) are this project's actual analyzed campaign",
        ha="center", va="top", fontsize=8, color="firebrick", wrap=True,
    )

    trail_lc = LineCollection([], cmap="viridis",
                               norm=plt.Normalize(epoch_yr.min(), epoch_yr.max()),
                               linewidths=1.3, alpha=0.85, zorder=2)
    ax.add_collection(trail_lc)
    real_lc = LineCollection([], colors="none", linewidths=3.2, alpha=0.95, zorder=3)
    ax.add_collection(real_lc)
    comet, = ax.plot([], [], "o", color="#d62728", ms=9, zorder=6, label="S301")
    # "lower right" (in data terms: least-negative RA, most-negative Dec) is
    # empty of trajectory for this orbit's geometry, unlike "upper right"
    # which sat directly on top of the periapsis convergence point near
    # Sgr A* -- see the additive-margin fix above for why that mattered.
    ax.legend(loc="lower right", fontsize=8)

    points = np.column_stack([ra_mas, dec_mas]).reshape(-1, 1, 2)
    segments_all = np.concatenate([points[:-1], points[1:]], axis=1)
    real_edge = is_real_campaign[:-1] & is_real_campaign[1:]

    def update(frame):
        end = frame + 1
        trail_lc.set_segments(segments_all[:end])
        trail_lc.set_array(epoch_yr[:end])
        real_segs = segments_all[:end][real_edge[:end]]
        real_lc.set_segments(real_segs)
        real_lc.set_color("#ff7f0e" if len(real_segs) else "none")
        comet.set_data([ra_mas[end - 1]], [dec_mas[end - 1]])
        title.set_text(f"S301 sky-plane track: precession over {N_ORBITS} orbits (real 1PN rate)")
        epoch_text.set_text(f"Epoch: {epoch_yr[end - 1]:.1f}  |  orbit #{orbit_index[end - 1] + 1} "
                             f"of {N_ORBITS}")
        inset_trail_lc.set_segments(segments_all[:end])
        inset_trail_lc.set_array(epoch_yr[:end])
        inset_comet.set_data([ra_mas[end - 1]], [dec_mas[end - 1]])
        return (trail_lc, real_lc, comet, title, epoch_text, banner,
                inset_trail_lc, inset_comet)

    n_frames = len(epoch_yr)
    anim = animation.FuncAnimation(fig, update, frames=n_frames, blit=False, interval=1000 / FPS)

    os.makedirs("output", exist_ok=True)
    writer = animation.FFMpegWriter(fps=FPS, bitrate=2400)
    anim.save(outpath, writer=writer, dpi=130)
    plt.close(fig)
    print(f"Saved {n_frames} frames ({n_frames / FPS:.1f} s at {FPS} fps) to {outpath}")


if __name__ == "__main__":
    make_movie()
