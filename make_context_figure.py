"""
Real observational context figure: a genuine near-infrared image of the
Sgr A* region (not a simulation, not this project's own synthetic data),
with S301's orbital scale pointed out relative to the real, already-
labeled S2/Sgr A* zoom inset.

Image: ESO/MPE/S. Gillessen et al., "Image of the Galactic Centre"
(announcement ann17051, released 2017-08-09), NACO/VLT near-infrared
observation. Licensed CC BY 4.0 by ESO's blanket policy for eso.org
public-site images (press releases/announcements/pictures of the week);
credit line reproduced below and in the figure caption per that license's
attribution requirement. Downloaded directly from
https://cdn.eso.org/images/large/ann17051c.jpg.

This image predates S301's announced discovery (August 2026) and was
never intended to show it, so no pixel-precise position for S301 is
plotted -- doing so would require this image's exact plate scale/epoch,
which are not established here. Instead, S301 is pointed out schematically:
its orbit (semi-major axis ~83 mas, comparable to or somewhat larger than
S2's) occupies the same innermost region ESO's own inset already isolates
around Sgr A* and S2, so an arrow to that same inset conveys the real
angular scale honestly without fabricating a precise overlay.
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

# pylint: disable=wrong-import-position
# matplotlib's backend must be set (above) before pyplot is imported --
# matches the same pattern used throughout this project's other plotting
# scripts (s301_lightcurve.py, campaign.py, fit_orbit.py).
import matplotlib.pyplot as plt
# pylint: enable=wrong-import-position

SRC_IMAGE = "output/sgra_eso_ann17051c_large.jpg"
OUT_IMAGE = "output/s301_field_context.png"
CREDIT = "Image credit: ESO/MPE/S. Gillessen et al. (ann17051), CC BY 4.0"


def main():
    """Load the real ESO field image, annotate S301's schematic context
    relative to the existing S2/Sgr A* inset, and save the figure.

    The source image is a network download (see module docstring) and is
    not tracked in the repository, while the annotated OUT_IMAGE is (it is
    a manuscript figure). If the source is absent, say so and leave the
    committed figure untouched rather than failing the whole run_all.sh
    pipeline over an optional re-render."""
    if not os.path.exists(SRC_IMAGE):
        print(f"{SRC_IMAGE} not found; keeping the committed {OUT_IMAGE}. To re-render, download "
              f"https://cdn.eso.org/images/large/ann17051c.jpg to that path and re-run.")
        return
    im = plt.imread(SRC_IMAGE)
    h, w = im.shape[0], im.shape[1]

    fig, ax = plt.subplots(figsize=(9, 7.6))
    ax.imshow(im)
    ax.axis("off")

    # Real, already-present ESO inset box (S2 / Sgr A* zoom), read off the
    # image's own pixel coordinates -- not redrawn, just referenced. Point
    # to its bottom-left corner (empty of the original S2/Sgr A* labels)
    # rather than overlapping them.
    inset_bottom_left = (845, 545)

    ax.annotate(
        "S301 (discovered 2026,\nnot shown here) orbits\nwithin this same\n"
        "innermost ~100 mas region",
        xy=inset_bottom_left, xytext=(0.03, 0.97), textcoords="axes fraction",
        ha="left", va="top", fontsize=11, color="white",
        bbox={"boxstyle": "round,pad=0.4", "fc": "black", "ec": "orange", "alpha": 0.8},
        arrowprops={"arrowstyle": "->", "color": "orange", "lw": 2.2,
                    "connectionstyle": "arc3,rad=0.25"},
    )

    ax.text(0.01, 0.01, CREDIT, transform=ax.transAxes, fontsize=7.5,
             color="white", va="bottom", ha="left",
             bbox={"boxstyle": "square,pad=0.25", "fc": "black", "alpha": 0.6, "ec": "none"})

    fig.tight_layout(pad=0.3)
    fig.savefig(OUT_IMAGE, dpi=200)
    print(f"Saved {OUT_IMAGE} ({w}x{h} source)")


if __name__ == "__main__":
    main()
