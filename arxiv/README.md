# arXiv submission package

This directory is a ready-to-upload arXiv source package, kept in sync with
the main `../article.tex`. It differs from the main copy only in that
`\includegraphics` paths are flattened (no `output/` prefix), matching the
flat directory structure arXiv expects.

Contents:
- `article.tex` -- the manuscript, with flattened figure paths.
- `article.bbl` -- pre-built bibliography. arXiv's own compile step does not
  reliably run bibtex, so the compiled `.bbl` is submitted directly instead
  of `references.bib` + relying on arXiv to generate it.
- `aastex701.cls` -- bundled in case arXiv's TeX Live version lags the one
  used to build this package.
- `sgra_eht_shadow.jpg`, `s301_field_context.png`, `s301_lightcurve.png`,
  `fit_orbit.png` -- the four figures.
- `anc/s301_precession.mp4` -- ancillary file. arXiv publishes anything in
  an `anc/` subdirectory of the source package as a downloadable file on the
  abstract page, separate from the paper's own figures.

## Rebuilding after `../article.tex` changes

```bash
cd /home/jcatanz/projects/s301
rm -rf arxiv/*.tex arxiv/*.bbl arxiv/*.cls arxiv/*.png arxiv/*.jpg
cp article.tex arxiv/article.tex
cp output/sgra_eht_shadow.jpg output/s301_field_context.png \
   output/s301_lightcurve.png output/fit_orbit.png arxiv/
cp output/s301_precession.mp4 arxiv/anc/
sed -i 's#output/sgra_eht_shadow\.jpg#sgra_eht_shadow.jpg#; \
        s#output/s301_field_context\.png#s301_field_context.png#; \
        s#output/s301_lightcurve\.png#s301_lightcurve.png#; \
        s#output/fit_orbit\.png#fit_orbit.png#' arxiv/article.tex
cp references.bib arxiv/references.bib

cd arxiv
docker run --rm -v "$(pwd)":/work -w /work texlive/texlive:latest bash -c "
  pdflatex -interaction=nonstopmode article.tex &&
  bibtex article &&
  pdflatex -interaction=nonstopmode article.tex &&
  pdflatex -interaction=nonstopmode article.tex &&
  cp /usr/local/texlive/2026/texmf-dist/tex/latex/aastex/aastex701.cls .
"
rm -f article.aux article.blg article.log article.out article.pdf references.bib
```

Verify it compiles standalone before committing -- this package has no
access to the parent directory's files, so a missing image or a stale
`.bbl` would only show up here, not in the main build.
