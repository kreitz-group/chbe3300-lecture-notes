---
title: Introduction
short_title: Introduction
label: ch-introduction
---

:::{admonition} Learning objectives
:class: tip
After completing this chapter, you should be able to:

- Explain why thermodynamic favorability does not guarantee that a reaction proceeds at a useful
  rate, and why reactor design therefore needs kinetics.
- Describe the role of kinetics and catalysis in the chemical process industry.
- Outline the four parts of the course and how each builds on the one before.
:::

## Why kinetics?

Ask what separates a chemical engineer from a mechanical engineer and the answer, stripped down,
is what each of them designs. The mechanical engineer thinks in engines, bearings, and gears. The
chemical engineer thinks in reactors, pipes, and pumps. Our distinguishing skill is that we
*engineer reactions*.

You have already had thermodynamics, so you can decide whether a reaction is favorable — whether
it *can* happen. But thermodynamics is silent on the question that decides whether a process is
viable: **how fast?** A reaction that equilibrium says will go to completion is worthless if it
takes a thousand years to get there. Designing a reactor means knowing how much time the system
needs, and that information comes from **reaction kinetics**.

Nearly everything around you, aside from wood and stone, is the product of an engineered chemical
reaction: fuels, plastics, pharmaceuticals, steel, paper, Gorilla Glass, transistors, and fermented
foods such as beer, bread, and cheese. This course covers the principles required to engineer
essentially all of it.

### The case of ammonia synthesis

Ammonia synthesis is the example we will return to throughout these notes, and it is worth stating
plainly why. It is among the most consequential chemical processes ever developed, and one that
almost nobody thinks about. Without it, only about half of the people alive today could be fed;
the rest would be competing for a badly limited food supply. The process turns the nitrogen in
ordinary air into fertilizer.

$$
\ce{N2 + 3 H2 <=> 2 NH3}
$$

Thermodynamics says this reaction is favorable. Yet it does not happen at any useful rate on its
own — it needs a **catalyst**. Part of this course is devoted to what a catalyst is and how it
works.

The equation above also hides a great deal. Written this way it suggests that a nitrogen molecule
and three hydrogen molecules collide and rearrange in one step. They do not: a large number of
electrons have to be moved around, and the transformation proceeds through a sequence of elementary
steps — a **reaction mechanism**. Uncovering and using such mechanisms is another thread running
through the course.

### Beyond the petrochemical industry

The field of kinetics grew out of the petrochemical industry, and most examples and textbooks —
including many in these notes — still reflect that origin. The applications, however, are far
broader. The concepts we cover apply equally to transformations in the natural world, including
geochemical and biological ones. Your chemical engineering training carries directly into those
areas.

## Course mechanics

Course logistics — schedule, assessment, office hours, and the syllabus itself — live on
**Canvas**, which is the authoritative source and is updated during the semester. These notes cover
the technical content only.

Two points about the notes themselves are worth stating here:

- **The course is taught in MATLAB.** Worked examples are given as MATLAB code, and MATLAB is what
  you will use for homework and exams. Georgia Tech provides it to all students through the
  [campus license](https://software.oit.gatech.edu/).
- **Each chapter opens with learning objectives.** Treat them as a checklist. If you can do
  everything on a chapter's list without looking anything up, you are ready for the exam on that
  chapter.

## What this course covers

The material is organized in four parts, each building on the one before.

**Reaction stoichiometry and thermodynamics.** The bookkeeping that every reacting system obeys:
how to encode a reaction mathematically, how to describe a mixture, and how to track its
composition as the reaction proceeds ([](#ch-stoichiometry)). Then the thermodynamic limits —
what equilibrium permits, and how it shifts with temperature and pressure
([](#ch-thermodynamics)).

**Reactor balances and rate laws.** Mass balances for the ideal reactor types
([](#ch-mass-balances)), the rate laws that go into them ([](#ch-reaction-kinetics)), what changes
when several reactions run at once ([](#ch-multiple-reactions)), and how rate laws are extracted
from experimental data in the first place ([](#ch-experiments)).

**Molecular basis of reaction rates.** Where a rate constant comes from — collisions, energy
barriers, and transition states ([](#ch-microscopic)) — and how elementary steps assemble into the
mechanisms behind observed global rate laws ([](#ch-mechanisms)).

**Catalysis.** How a solid catalyst changes a reaction and how to write rate laws for catalytic
cycles ([](#ch-hetcat)), and how nature solves the same problem with enzymes
([](#ch-biocatalysis)).
