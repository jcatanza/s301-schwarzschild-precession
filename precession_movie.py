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


def build_truth_elements():
    """The real published orbit (constants.TRUTH), matching every other
    script in this project."""
    period = k.TRUTH["P_yr"] * k.year
    sma = orbit.semi_major_axis_from_period(k.GM_BH, period)
    return orbit.OrbitalElements(
        t_peri=0.0, period=period, ecc=k.TRUTH["e"], sma=sma,
        i_deg=k.TRUTH["i_deg"], raan_deg=k.TRUTH["Omega_deg"], omega_deg=k.TRUTH["omega_deg"],
    ), sma, period


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


def make_movie(outpath="output/s301_precession.mp4"):
    """Render and save the animation."""
    epoch_yr, ra_mas, dec_mas, orbit_index, omega_dot = build_trajectory()
    omega_dot_deg_yr = np.degrees(omega_dot) * k.year
    is_real_campaign = np.isin(orbit_index, REAL_CAMPAIGN_ORBIT_INDICES)

    fig, ax = plt.subplots(figsize=(8, 8))
    margin = 1.08
    ax.set_xlim(np.max(ra_mas) * margin, np.min(ra_mas) * margin)   # inverted RA, like other figures
    ax.set_ylim(np.min(dec_mas) * margin, np.max(dec_mas) * margin)
    ax.set_aspect("equal")
    ax.set_xlabel("RA offset (mas)")
    ax.set_ylabel("Dec offset (mas)")
    ax.plot(0, 0, "k*", ms=16, zorder=5, label="Sgr A*")
    ax.grid(alpha=0.3)

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
    ax.legend(loc="upper right", fontsize=8)

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
        return trail_lc, real_lc, comet, title, epoch_text, banner

    n_frames = len(epoch_yr)
    anim = animation.FuncAnimation(fig, update, frames=n_frames, blit=False, interval=1000 / FPS)

    os.makedirs("output", exist_ok=True)
    writer = animation.FFMpegWriter(fps=FPS, bitrate=2400)
    anim.save(outpath, writer=writer, dpi=130)
    plt.close(fig)
    print(f"Saved {n_frames} frames ({n_frames / FPS:.1f} s at {FPS} fps) to {outpath}")


if __name__ == "__main__":
    make_movie()
