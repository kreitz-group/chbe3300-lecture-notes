# Figure sources

TikZ/pgfplots sources for figures that were redrawn from hand sketches. Each
`<name>.tex` produces `../<name>.png`; `hetcat-style.tex` is the shared style for
the ch. 10 plots. The curves are computed from the rate laws in the chapter
(parameters are noted in each file), but the axes carry no numbers: these are
schematics.

Build (from this folder, auxiliary files go to a scratch directory):

```
lualatex -output-directory=/tmp/fig <name>.tex
pdftocairo -png -r 300 -singlefile /tmp/fig/<name>.pdf /tmp/fig/<name>
magick /tmp/fig/<name>.png -background white -alpha remove -alpha off \
  -resize '1600x1600>' -strip -define png:compression-level=9 ../<name>.png
```

The PDFs in `/tmp/fig` can also be dropped into the LaTeX notes.
