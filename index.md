---
title: Chemical Kinetics and Catalysis
subtitle: CHBE 3300 · Georgia Institute of Technology
---

These are the lecture notes for **CHBE 3300, Chemical Kinetics and Catalysis**, taught in the
School of Chemical and Biomolecular Engineering at the Georgia Institute of Technology.

Thermodynamics tells us whether a reaction *can* happen. It says nothing about *how fast*.
This course supplies the missing half: how to measure a rate, how to write a rate law, where
that rate law comes from at the molecular level, and how a catalyst changes it. Nearly every
material around you — fuels, plastics, pharmaceuticals, steel, fertilizer, semiconductors — is
made in an engineered chemical reaction, and designing those reactions is what a chemical
engineer does.

## How to use these notes

Each chapter opens with a set of **learning objectives**. Use them as a checklist: if you can do
everything on the list without looking anything up, you are ready for the exam on that chapter.

Many chapters include executable examples. Every page that contains code can be launched
directly in your browser or downloaded as a Jupyter notebook using the {kbd}`↓` and
{kbd}`🚀` buttons at the top right of the page — nothing needs to be installed to read along.
If you would rather run things locally, see [](#running-locally).

## Software

The computational examples use the standard scientific Python stack:

- [NumPy](https://numpy.org/) for array math
- [SciPy](https://scipy.org/) for numerical integration, optimization, and regression
- [Matplotlib](https://matplotlib.org/) for plotting
- [Pint](https://pint.readthedocs.io/) for unit-aware calculations
- [Cantera](https://cantera.org/) for thermodynamic and kinetic data

(running-locally)=
### Running locally

```bash
git clone https://github.com/kreitz-group/chbe3300-lecture-notes.git
cd chbe3300-lecture-notes
conda env create -f environment.yml
conda activate chbe3300
```

To build and preview the book itself:

```bash
jupyter book start
```

## Acknowledgments

The structure of this book is loosely modeled on Kyle Niemeyer's
[Computational Thermodynamics](https://kyleniemeyer.github.io/computational-thermo/), an
excellent example of what interactive course notes can be.

:::{note} Work in progress
These notes are being ported from the LaTeX source used in previous offerings of the course.
Chapters will fill in over the semester. If you find an error, please
[open an issue](https://github.com/kreitz-group/chbe3300-lecture-notes/issues).
:::

## License

The text and figures are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); the code is licensed under the
[MIT License](https://opensource.org/licenses/MIT).
