---
title: Reaction Kinetics
short_title: Reaction kinetics
label: ch-reaction-kinetics
---

<!-- LaTeX source: ReactionKinetics.tex -->

:::{admonition} Learning objectives
:class: tip
After completing this chapter, you should be able to:

- Define the rate of reaction $r$ and the species production rate $r_i$, and relate them through
  the stoichiometric coefficient via $r_i = \nu_i r$.
- State Boudart's five rules for reaction kinetics and apply them to write power-law rate
  expressions.
- Use the Arrhenius equation to relate rate constants at different temperatures and extract
  $E\un{a}$ and $A$ from a $\ln k$ vs. $1/T$ plot.
- Distinguish elementary reactions from global (lumped) reactions and recognize when the law of
  mass action applies.
- Derive and use integrated rate laws for irreversible zero-, first-, and second-order reactions
  in a constant-volume batch reactor.
- Recast batch-reactor solutions in dimensionless form using the Damköhler number $Da_I$ and the
  conversion $X$ (or remaining fraction $f$).
- Apply the Euler–Cauchy method to integrate a batch-reactor ODE numerically and compare against
  the analytical solution.
- Solve a reversible first-order reaction analytically and recover the equilibrium constant as
  $K = k\un{fwd} / k\un{rev}$.
:::

## Rate of reaction

<!-- source: ReactionKinetics.tex L4-85 -->

## Boudart's rules for reaction kinetics

<!-- source: ReactionKinetics.tex L85 -->

## Integrated rate laws

<!-- source: ReactionKinetics.tex L287 -->

### Irreversible, first-order reaction

<!-- source: ReactionKinetics.tex L297 -->

### Irreversible, zero-order reaction

<!-- source: ReactionKinetics.tex L366 -->

### Irreversible, second-order reaction

<!-- source: ReactionKinetics.tex L397 -->

### Dimensionless equations

<!-- source: ReactionKinetics.tex L440 -->

### Second-order irreversible reactions, first order in both reactants

<!-- source: ReactionKinetics.tex L511 -->

## Interlude: solving ordinary differential equations numerically

<!-- source: ReactionKinetics.tex L581 -->

See [](../notebooks/batch-reactor-first-order.ipynb) for a worked numerical example.

## Reversible, first-order reaction

<!-- source: ReactionKinetics.tex L634 -->

## Summary

<!-- source: ReactionKinetics.tex L738 -->
