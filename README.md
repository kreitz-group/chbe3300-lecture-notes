# chbe3300-lecture-notes

Interactive lecture notes for **CHBE 3300 — Chemical Kinetics and Catalysis** at the Georgia
Institute of Technology, built with [Jupyter Book 2](https://next.jupyterbook.org/) / MyST.

📖 **Read the book:** https://kreitz-group.github.io/chbe3300-lecture-notes/

## Repository layout

```
myst.yml                    # book configuration: metadata, math macros, table of contents
index.md                    # landing page
content/                    # one Markdown file per chapter, mirroring the LaTeX source
notebooks/                  # executable computational examples
figures/                    # figures shared with the LaTeX edition
references.bib              # bibliography, shared with the LaTeX edition
environment.yml             # conda environment for running the notebooks locally
requirements.txt            # pip dependencies used by the CI build
.github/workflows/deploy.yml  # builds and publishes to GitHub Pages on push to main
```

The chapter files mirror the LaTeX lecture notes one-to-one. Each stub carries the chapter's
learning objectives and its section headings; HTML comments such as
`<!-- source: Stoichiometry.tex L118 -->` point at the line in the LaTeX source that the section
is to be ported from.

## Building locally

```bash
conda env create -f environment.yml
conda activate chbe3300

jupyter book start          # live preview at http://localhost:3000
jupyter book build --html   # static build into _build/html
```

## Publishing

Pushing to `main` triggers the GitHub Actions workflow, which builds the book and deploys it to
GitHub Pages. This requires **Settings → Pages → Source → GitHub Actions** to be enabled once for
the repository.

## License

Text and figures: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Code: [MIT](LICENSE).
