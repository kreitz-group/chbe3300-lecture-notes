---
title: Stoichiometry and Definitions
short_title: Stoichiometry
label: ch-stoichiometry
---

<!-- LaTeX source: Stoichiometry.tex -->
<!-- NOTE on porting: mhchem here is math-only, so \ce{} in prose must be wrapped in $...$.
     And \un{} expands to _{\textrm{...}}, which puts \ce{} into text mode and breaks KaTeX --
     write species subscripts as _{\ce{N2}}, never \un{\ce{N2}}. -->

:::{admonition} Learning objectives
:class: tip
After completing this chapter, you should be able to:

- Identify reactants, products, and stoichiometric coefficients for any chemical reaction, and
  translate freely between the equation form and the canonical form $0 = \sum_i \nu_i A_i$.
- Distinguish stoichiometric numbers (positive) from stoichiometric coefficients (signed).
- Describe a reaction mixture using the appropriate extensive and intensive quantities
  ($n_i$, $m_i$, $c_i$, $y_i$, $p_i$).
- Use the extent of reaction $\xi$ to track composition changes for both single and multiple
  reactions.
- Identify the limiting reactant of a reaction system and compute conversion and the maximum
  extent.
- Set up and solve for the equilibrium composition of a closed system from a given equilibrium
  constant.
:::

## Reactions and stoichiometric coefficients

<!-- source: Stoichiometry.tex L16 -->

Stoichiometry defines a set of rules that a reaction mixture follows over the course of reaction.
A *reaction mixture* is the summary of all species present in a reaction system — a chemical
reactor, for instance. The reason a reaction mixture obeys stoichiometric rules is the conservation
of the number of atoms of every element over the course of the reaction.

A generic chemical reaction is written

$$
\underbrace{\ce{\alpha A + \beta B}}_{\text{reactants}}\ \ce{->}\
\underbrace{\ce{\gamma C + \delta D}}_{\text{products}} ,
$$ (eq-generic-reaction)

where A, B, C, and D are the components, or *species*. The species on the left-hand side are the
**reactants** and those on the right-hand side are the **products**; $\alpha$, $\beta$, $\gamma$,
and $\delta$ are the **stoichiometric numbers**.

The same reaction can be written in a purely mathematical form, using a set of species $A_i$, in
the *canonical form*

$$
0 = \sum_{i=1}^{N\un{species}} \nu_i A_i ,
$$ (eq-canonical-form)

where the $\nu_i$ are the **stoichiometric coefficients**. This is the form that lets us evaluate
chemical reactions mathematically. The stoichiometric coefficients carry a sign:

$$
\begin{aligned}
\nu_i &< 0 \quad \text{for reactants,} \\
\nu_i &> 0 \quad \text{for products.}
\end{aligned}
$$

The difference between stoichiometric numbers and stoichiometric coefficients is precisely this
sign: stoichiometric *numbers* are always positive, whereas stoichiometric *coefficients* may be
positive or negative. Applying this to [](#eq-generic-reaction) gives $\nu\un{A} = -\alpha$,
$\nu\un{B} = -\beta$, $\nu\un{C} = \gamma$, and $\nu\un{D} = \delta$.

### Example: ammonia synthesis

The reaction equation for $\ce{NH3}$ synthesis is

$$
\ce{N2 + 3 H2 <=> 2 NH3} ,
$$ (eq-nh3-rxn)

or, with $A_1 = \ce{N2}$, $A_2 = \ce{H2}$, and $A_3 = \ce{NH3}$,

$$
A_1 + 3\,A_2 \eq 2\,A_3 ,
$$

with the stoichiometric coefficients $\nu_1 = -1$, $\nu_2 = -3$, $\nu_3 = 2$. The same reaction may
equally well be written as any of

$$
\begin{aligned}
&0 = -\ce{N2} - 3\,\ce{H2} + 2\,\ce{NH3} \\
&\ce{N2} + 3\,\ce{H2} = 2\,\ce{NH3} \\
&0 = -\tfrac{1}{2}\,\ce{N2} - \tfrac{3}{2}\,\ce{H2} + \ce{NH3} .
\end{aligned}
$$

::::{admonition} Discussion
:class: seealso
Does it matter how you write the reaction equation?

:::{dropdown} Answer
Not really. How you write the reaction equation is somewhat arbitrary and largely a matter of
personal preference. What matters is that the *ratios* of the stoichiometric coefficients are
always the same:

$$
\frac{\nu_{\ce{N2}}}{\nu_{\ce{NH3}}} = \frac{-1}{2} = \frac{-1/2}{1} = -\frac{1}{2} .
$$
:::
::::

### Multiple reactions

In many chemical systems several reactions occur at the same time. It is then convenient to write
the chemical reactions in linear-algebra notation, using a stoichiometry vector $\vect{\nu}$ and a
species vector $\vect{A}$:

$$
\begin{aligned}
&\vect{\nu} \cdot \vect{A}^{\T} = 0 , \\
\text{with}\quad & \vect{\nu} = \left(\nu_1\ \nu_2\ \dots\ \nu_n\right) \\
\text{and}\quad & \vect{A} = \left(A_1\ A_2\ \dots\ A_n\right) .
\end{aligned}
$$

These matrix notations are the ones you will actually use when solving kinetic problems with
numerical tools.

## Quantities used to describe a reaction mixture

<!-- source: Stoichiometry.tex L86 -->

Before we can quantify reaction progress, we need a vocabulary for describing the mixture itself.
Quantities used to characterize a reaction mixture fall into two classes: **extensive** ones, which
are proportional to system size, and **intensive** ones, which are independent of it.

### Extensive quantities

- Number of moles of species $i$: $n_i$, in mol.
- Mass of species $i$: $m_i$, in g.

The two are related by the molar mass $M_i$ of the species,

$$
m_i = n_i M_i .
$$

The total moles and total mass of the mixture are

$$
n = \sum_{i=1}^{N\un{species}} n_i ,
\qquad
m = \sum_{i=1}^{N\un{species}} m_i .
$$

A chemical reaction changes the moles of each species and may change the total moles $n$ — for
example, $\ce{N2 + 3 H2 <=> 2 NH3}$ converts four moles of reactants into two moles of product.
The total mass $m$, however, is always conserved.

Extensive quantities scale with the system: doubling the size of the reactor doubles $n$ and $m$.

### Intensive quantities

Intensive quantities are independent of system size and are therefore the natural way to describe
*compositions*, and to express rate laws and equilibrium relations:

$$
\begin{aligned}
\text{molar concentration:}\quad & c_i = \frac{n_i}{V} \quad \text{in}\ \mathrm{mol\,L^{-1}}, \\
\text{mole fraction:}\quad & y_i = \frac{n_i}{n} = \frac{n_i}{\sum_k n_k} , \\
\text{partial pressure (gas phase):}\quad & p_i = y_i\, p ,
\end{aligned}
$$

where $V$ is the reaction volume and $p$ the total pressure. By construction $\sum_i y_i = 1$ and
$\sum_i p_i = p$, the latter being Dalton's law.

## Extent of reaction

<!-- source: Stoichiometry.tex L118 -->

In kinetics we are interested in the *change* of a system brought about by chemical reaction. Stay
with the ammonia synthesis example, $\ce{N2 + 3 H2 <=> 2 NH3}$.

Because the number of atoms of each element is conserved, the ratio of the change in molar amount
to the stoichiometric coefficient is the same for every species as the reaction proceeds:

$$
\frac{\mathrm{d}n_{\ce{N2}}}{\nu_{\ce{N2}}}
= \frac{\mathrm{d}n_{\ce{H2}}}{\nu_{\ce{H2}}}
= \frac{\mathrm{d}n_{\ce{NH3}}}{\nu_{\ce{NH3}}}
= \frac{\mathrm{d}n_i}{\nu_i}
= \mathrm{d}\xi .
$$

The quantity $\xi$ is the **extent of reaction**. Integrating gives

$$
\xi = \frac{n_i - n\un{i,0}}{\nu_i} ,
$$

where $n\un{i,0}$ is the number of moles at the beginning of the reaction. The useful consequence
is that the entire progress of a single reaction is described by one number. Rearranged,

$$
n_i = n\un{i,0} + \nu_i \xi ,
$$ (eq-extent-integrated)

which lets us follow the progress of a reaction in a closed system — for instance, to obtain the
composition of an equilibrium mixture from a given starting composition.

::::{admonition} In-class exercise
:class: note
Suppose we start with 1 mol of $\ce{N2}$ and 2 mol of $\ce{H2}$, and we convert 0.5 mol of
$\ce{N2}$. How many moles of $\ce{H2}$ are left: 1.5 mol, 1 mol, or 0.5 mol?

:::{dropdown} Solution
Apply [](#eq-extent-integrated) to $\ce{N2}$:

$$
\xi = \frac{n_{\ce{N2}} - n_{\ce{N2},0}}{\nu_{\ce{N2}}}
    = \frac{0.5\ \mathrm{mol} - 1\ \mathrm{mol}}{-1}
    = 0.5\ \mathrm{mol} .
$$

Then for $\ce{H2}$,

$$
n_{\ce{H2}} = 2\ \mathrm{mol} + (-3)(0.5\ \mathrm{mol}) = 0.5\ \mathrm{mol} .
$$

Each mole of $\ce{N2}$ that reacts consumes three moles of $\ce{H2}$, so the answer is 0.5 mol.
:::
::::

We will use the extent of reaction later to define the rate of reaction. Note that $\xi$ is an
*extensive* property and therefore depends on system size.

For multiple reactions the definition generalizes directly,

$$
n_i = n\un{i,0} + \sum_{j=1}^{N\un{reactions}} \nu_{i,j}\, \xi_j ,
$$

with a separate extent of reaction for each independent reaction.

Another useful quantity is the **conversion**, an intensive measure of reaction progress:

$$
X = \frac{n\un{i,0} - n_i}{n\un{i,0}} = 1 - \frac{n_i}{n\un{i,0}} ,
$$

which lies between 0 and 1.

When reactants are not present in stoichiometric proportions, one of them runs out first — the
**limiting reactant**. The limiting reactant sets the bounds on the extent of reaction: the maximum
extent corresponds to complete conversion ($X = 1$) of the limiting reactant,

$$
\xi\un{max} = -\frac{n\un{lim,0}}{\nu\un{lim}} .
$$ (eq-xi-max)

::::{admonition} Discussion
:class: seealso
In the example above — 1 mol $\ce{N2}$ and 2 mol $\ce{H2}$ — can we achieve complete conversion of
*both* reactants?

:::{dropdown} Answer
No. Stoichiometry requires three moles of $\ce{H2}$ for every mole of $\ce{N2}$. With 1 mol
$\ce{N2}$ and only 2 mol $\ce{H2}$ on hand, $\ce{H2}$ runs out first: it is the limiting reactant.
The maximum extent is

$$
\xi\un{max} = -\frac{2\ \mathrm{mol}}{-3} = \tfrac{2}{3}\ \mathrm{mol} ,
$$

at which point all $\ce{H2}$ has been consumed and $1 - 2/3 = 1/3$ mol of $\ce{N2}$ remains
unreacted.
:::
::::

## Equilibrium composition

<!-- source: Stoichiometry.tex L209 -->

The extent of reaction lets us derive the composition of a mixture from a starting composition and
a given extent. Most frequently the method is used the other way around: to find the composition at
equilibrium for a given equilibrium constant. Equilibrium and thermodynamics are treated in detail
in [](#ch-thermodynamics).

### Single-reaction example: ammonia synthesis

<!-- source: Stoichiometry.tex L216 -->

Return to ammonia synthesis,

$$
\ce{3 H2 + N2 <=> 2 NH3} ,
$$

under the conditions $T = 450\ \mathrm{K}$, $p = 4.25\ \mathrm{bar}$, and
$K\un{p} = 1.397\ \mathrm{bar^{-2}}$. The equilibrium constant is defined by the **law of mass
action**,

$$
K\un{p} = \frac{p_{\ce{NH3}}^2}{p_{\ce{N2}}\, p_{\ce{H2}}^3}
        = \frac{1}{p^2}\, \frac{y_{\ce{NH3}}^2}{y_{\ce{N2}}\, y_{\ce{H2}}^3} ,
$$ (eq-nh3-kp)

where we have used $p_i = p\, y_i$. Initially the reactor contains

$$
n_{\ce{N2},0} = 1\ \mathrm{mol}, \qquad
n_{\ce{H2},0} = 2\ \mathrm{mol}, \qquad
n_{\ce{NH3},0} = 0\ \mathrm{mol} .
$$

The equilibrium composition is most easily worked out in a table.

:::{table} Determination of the equilibrium composition of $\ce{NH3}$ synthesis using the extent of reaction.
:label: tab-eq-nh3-synthesis

| $i$ | $A_i$ | $n\un{i,0}$ (mol) | $n_i = n\un{i,0} + \nu_i \xi$ | $y_i = n_i / n$ |
|:---:|:-----:|:-----------------:|:-----------------------------:|:---------------:|
| 1 | $\ce{N2}$  | 1 | $1 - \xi$   | $\dfrac{1-\xi}{3-2\xi}$   |
| 2 | $\ce{H2}$  | 2 | $2 - 3\xi$  | $\dfrac{2-3\xi}{3-2\xi}$  |
| 3 | $\ce{NH3}$ | 0 | $2\xi$      | $\dfrac{2\xi}{3-2\xi}$    |
|   | $\sum$     | 3 | $3 - 2\xi$  | 1 |
:::

Substituting the mole fractions from [](#tab-eq-nh3-synthesis) into [](#eq-nh3-kp),

$$
K\un{p} = \frac{1}{p^2}\,
\frac{\left(\dfrac{2\xi}{3-2\xi}\right)^{\!2}}
     {\dfrac{1-\xi}{3-2\xi}\left(\dfrac{2-3\xi}{3-2\xi}\right)^{\!3}} .
$$

Inserting the given numbers and clearing denominators leaves

$$
\underbrace{25}_{K\un{p}\,p^2}\left(1-\xi\right)\left(2-3\xi\right)^3
= \left(3-2\xi\right)^2\left(2\xi\right)^2 .
$$

This is a *quartic* equation — degree four on both sides. It cannot be solved with the quadratic
formula, so we must fall back on a numerical method such as Newton–Raphson.

We can, however, bound the solution before searching. The limiting reactant is $\ce{H2}$, so the
maximum extent of reaction follows from [](#eq-xi-max):

$$
0 = 2 - 3\xi \quad \Rightarrow \quad \xi\un{max} = \tfrac{2}{3}\ \mathrm{mol} .
$$

Solving numerically on $[0, \xi\un{max}]$ gives $\xi = 0.453\ \mathrm{mol}$.

### Multi-reaction example: methanol synthesis

<!-- source: Stoichiometry.tex L276 -->

For $M$ reactions we write one stoichiometric equation per reaction $j$,

$$
0 = \sum_{i=1}^{N\un{species}} \nu_{i,j} A_i , \qquad j = 1, \dots, M ,
$$

and describe the progress of the system with one extent per reaction,

$$
n_i = n\un{i,0} + \sum_{j=1}^{M} \nu_{i,j}\, \xi_j ,
$$

or, in matrix notation,

$$
\vect{n} = \vect{n}_0 + \mtrx{N}\, \vect{\xi} ,
$$

where $\mtrx{N}$ is the matrix of stoichiometric coefficients and $\vect{\xi}$ the vector of
extents.

Take methanol synthesis as the example:

$$
\begin{aligned}
(1)\quad & \ce{CO + 2 H2 <=> CH3OH} \\
(2)\quad & \ce{CO2 + 3 H2 <=> CH3OH + H2O} \\
(3)\quad & \ce{H2O + CO <=> CO2 + H2}
\end{aligned}
$$

The first two are the key reactions of the system. The third is a linear combination of the first
two and therefore need not be considered separately — it introduces no new chemistry and no new
degree of freedom. We consequently need two equilibrium constants,

$$
K\un{p,1} = \frac{p_{\ce{CH3OH}}}{p_{\ce{CO}}\, p_{\ce{H2}}^2} ,
\qquad
K\un{p,2} = \frac{p_{\ce{CH3OH}}\, p_{\ce{H2O}}}{p_{\ce{CO2}}\, p_{\ce{H2}}^3} .
$$

Assume we start with $n_{\ce{CO},0} = 1\ \mathrm{mol}$, $n_{\ce{CO2},0} = 0.1\ \mathrm{mol}$,
and $n_{\ce{H2},0} = 2\ \mathrm{mol}$. In analogy with the single-reaction case, a table gives
the $n_i$ and $y_i$.

:::{table} Determination of the equilibrium composition of $\ce{CH3OH}$ synthesis using one extent of reaction per independent reaction.
:label: tab-eq-ch3oh-synthesis

| $i$ | $A_i$ | $n\un{i,0}$ (mol) | $n_i = n\un{i,0} + \sum_j \nu_{i,j}\xi_j$ | $y_i = n_i / n$ |
|:---:|:-----:|:-----------------:|:----------------------------------------:|:---------------:|
| 1 | $\ce{CO}$    | 1   | $1 - \xi_1$              | $\dfrac{1-\xi_1}{3.1-2\xi_1-2\xi_2}$              |
| 2 | $\ce{CO2}$   | 0.1 | $0.1 - \xi_2$            | $\dfrac{0.1-\xi_2}{3.1-2\xi_1-2\xi_2}$            |
| 3 | $\ce{H2}$    | 2   | $2 - 2\xi_1 - 3\xi_2$    | $\dfrac{2-2\xi_1-3\xi_2}{3.1-2\xi_1-2\xi_2}$      |
| 4 | $\ce{CH3OH}$ | 0   | $\xi_1 + \xi_2$          | $\dfrac{\xi_1+\xi_2}{3.1-2\xi_1-2\xi_2}$          |
| 5 | $\ce{H2O}$   | 0   | $\xi_2$                  | $\dfrac{\xi_2}{3.1-2\xi_1-2\xi_2}$                |
|   | $\sum$       | 3.1 | $3.1 - 2\xi_1 - 2\xi_2$  | 1 |
:::

This nonlinear system of two equations in two unknowns is too complex to solve by hand, so we rely
on a numerical solution — in MATLAB, for example — to determine the equilibrium composition for a
known pair $K\un{p,1}$, $K\un{p,2}$.

## Summary

<!-- source: Stoichiometry.tex L328 -->

- A chemical reaction is encoded compactly as $0 = \sum_i \nu_i A_i$, with $\nu_i < 0$ for
  reactants and $\nu_i > 0$ for products. The unsigned counterparts are the stoichiometric numbers.
- A reaction mixture is described by extensive ($n_i$, $m_i$) and intensive ($c_i$, $y_i$, $p_i$)
  quantities. Total mass is conserved by every reaction; total moles need not be.
- For a single reaction, the extent of reaction $\xi$ tracks composition via
  $n_i = n\un{i,0} + \nu_i \xi$. Conversion is $X_i = 1 - n_i/n\un{i,0}$, and the limiting reactant
  sets $\xi\un{max} = -n\un{lim,0}/\nu\un{lim}$.
- For $M$ reactions, $n_i = n\un{i,0} + \sum_j \nu_{i,j}\xi_j$, equivalently
  $\vect{n} = \vect{n}_0 + \mtrx{N}\vect{\xi}$. Only linearly independent reactions need their own
  extent.
- The equilibrium composition follows from a stoichiometric table together with the law of mass
  action; the resulting equation is generally nonlinear and solved numerically. Always bound the
  solution to the physically allowed range $[0, \xi\un{max}]$ before searching.
