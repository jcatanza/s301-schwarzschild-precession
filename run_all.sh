#!/usr/bin/env bash
# Regenerate everything the S301 manuscript depends on, in dependency order:
# every analysis script (each writes results/<name>.json), numbers.tex from
# those JSONs, article.pdf (Docker TeX Live), and the flattened arxiv/
# submission package.
#
# Usage:
#   ./run_all.sh              full pipeline (~1 h; optimal_design.py and the
#                             movie render dominate)
#   ./run_all.sh --skip-slow  skip optimal_design.py, cadence_alternatives.py
#                             and precession_movie.py (~15 min); the committed
#                             results/*.json and output/ artefacts for those
#                             three steps are left in place.
#
# Requires the Python environment in requirements.txt (a .venv/ next to this
# script is activated automatically if present), ffmpeg for the movie, and
# Docker for the PDF and arxiv/ steps (both are skipped, with a notice, if
# `docker` is not on PATH).
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

usage() {
    sed -n '2,18p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

SKIP_SLOW=0
for arg in "$@"; do
    case "$arg" in
        --skip-slow) SKIP_SLOW=1 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "run_all.sh: unknown option '$arg'" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -f .venv/bin/activate ]]; then
    # shellcheck disable=SC1091
    source .venv/bin/activate
fi
PYTHON=${PYTHON:-python}

# Steps that take more than a few minutes; skipped under --skip-slow.
SLOW_STEPS=(optimal_design.py cadence_alternatives.py precession_movie.py)

# The pipeline, in dependency order. campaign.py writes the synthetic CSV
# that fit_orbit.py and every systematics script read; make_numbers.py must
# run last so numbers.tex reflects every results/*.json.
STEPS=(
    s301_lightcurve.py
    campaign.py
    fit_orbit.py
    optimal_design.py
    cadence_alternatives.py
    photometry_checks.py
    rv_channel_test.py
    precision_sensitivity.py
    reference_frame_error.py
    extended_mass_error.py
    confusion_error.py
    lensing_error.py
    roemer_delay.py
    spin_contamination.py
    mass_distance_test.py
    solar_conjunction.py
    error_budget.py
    precession_movie.py
    make_context_figure.py
    make_numbers.py
)

TEX_IMAGE=texlive/texlive:latest
FIGURES=(output/sgra_eht_shadow.jpg output/s301_field_context.png
         output/s301_lightcurve.png output/fit_orbit.png)
MOVIE=output/s301_precession.mp4

is_slow() {
    local step
    for step in "${SLOW_STEPS[@]}"; do
        [[ "$step" == "$1" ]] && return 0
    done
    return 1
}

elapsed() {
    printf '%dm%02ds' $(( $1 / 60 )) $(( $1 % 60 ))
}

run_step() {
    local script=$1 t0 t1
    if (( SKIP_SLOW )) && is_slow "$script"; then
        echo "-- skipping $script (--skip-slow)"
        return
    fi
    echo "== $script"
    t0=$(date +%s)
    "$PYTHON" "$script"
    t1=$(date +%s)
    echo "-- $script done in $(elapsed $(( t1 - t0 )))"
}

docker_tex() {
    # Run a shell command inside TeX Live with the given directory mounted at /work.
    local dir=$1 cmd=$2
    docker run --rm -v "$(cd "$dir" && pwd)":/work -w /work "$TEX_IMAGE" bash -c "$cmd"
}

build_pdf() {
    echo "== article.pdf (Docker $TEX_IMAGE)"
    local t0 t1
    t0=$(date +%s)
    docker_tex . "pdflatex -interaction=nonstopmode article.tex && bibtex article && \
        pdflatex -interaction=nonstopmode article.tex && pdflatex -interaction=nonstopmode article.tex"
    t1=$(date +%s)
    echo "-- article.pdf done in $(elapsed $(( t1 - t0 )))"
}

regenerate_arxiv() {
    # arXiv wants a flat source directory: article.tex with the output/
    # prefix stripped from the four \includegraphics paths, numbers.tex
    # (article.tex \input's it), a pre-built article.bbl (arXiv does not
    # reliably run bibtex), the class file, the figures, and the movie as
    # an ancillary file under anc/.
    echo "== arxiv/ submission package"
    local t0 t1 fig
    t0=$(date +%s)
    mkdir -p arxiv/anc
    rm -f arxiv/*.tex arxiv/*.bbl arxiv/*.cls arxiv/*.bib arxiv/*.png arxiv/*.jpg arxiv/*.pdf
    cp article.tex numbers.tex references.bib arxiv/
    for fig in "${FIGURES[@]}"; do cp "$fig" arxiv/; done
    cp "$MOVIE" arxiv/anc/
    sed -i -E 's#(\\includegraphics(\[[^]]*\])?\{)output/#\1#' arxiv/article.tex
    # Only figure paths must be flattened; the GitHub URL in the Fig. 4 caption
    # legitimately contains "output/".
    if grep -q -E '\\includegraphics(\[[^]]*\])?\{output/' arxiv/article.tex; then
        echo "run_all.sh: arxiv/article.tex still has an output/ figure path -- check the sed above" >&2
        exit 1
    fi
    docker_tex arxiv "pdflatex -interaction=nonstopmode article.tex && bibtex article && \
        pdflatex -interaction=nonstopmode article.tex && pdflatex -interaction=nonstopmode article.tex && \
        cp \"\$(kpsewhich aastex701.cls)\" ."
    rm -f arxiv/article.aux arxiv/article.blg arxiv/article.log arxiv/article.out \
          arxiv/article.pdf arxiv/references.bib
    t1=$(date +%s)
    echo "-- arxiv/ done in $(elapsed $(( t1 - t0 )))"
}

T_START=$(date +%s)
for step in "${STEPS[@]}"; do
    run_step "$step"
done

if command -v docker >/dev/null 2>&1; then
    build_pdf
    regenerate_arxiv
else
    echo "-- docker not found: skipping article.pdf and arxiv/ regeneration"
fi

echo "== all done in $(elapsed $(( $(date +%s) - T_START )))"
