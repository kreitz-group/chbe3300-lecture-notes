---
title: Extracting Rate Laws from Experiments
short_title: Experimental kinetics
label: ch-experiments
---

<!-- LaTeX source: Experiments.tex -->
<!-- Porting notes: mhchem is math-only here, so \ce{} in prose needs $...$; and \un{} expands to
     _{\textrm{...}}, so species subscripts must be written _{\ce{CH4}}, never \un{\ce{CH4}}.
     Nested directives need a longer outer fence (:::: around :::).
     The nicematrix environments in the source have no KaTeX equivalent: the two annotated
     stoichiometric matrices, which carry species row labels and reaction column headers, are
     rendered as MyST tables; the plain numeric matrices stay as pmatrix. -->

:::{admonition} Learning objectives
:class: tip
After completing this chapter, you should be able to:

- Explain why a kinetic experiment measures the extent of reaction, or a related observable, rather
  than the reaction rate itself.
- State the rules a kinetic measurement must satisfy, and distinguish chemical from physical
  (*in situ* / *in operando*) measurement methods.
- Apply the integral method: postulate a rate law, integrate the batch-reactor mass balance,
  linearize, and extract $k$ and $n$ graphically.
- Apply the differential method: estimate $\mathrm{d}c\un{A}/\mathrm{d}t$ from concentration data
  using forward, central, or higher-order finite-difference formulas, then linearize
  $\ln(-\mathrm{d}c\un{A}/\mathrm{d}t)$ vs. $\ln c\un{A}$ to read off $n$ and $k$.
- Linearize the Arrhenius equation and extract $E\un{a}$ and $A$ from a $\ln k$ vs. $1/T$ plot, or
  from rate constants at two temperatures.
- Set up and solve a nonlinear regression problem to fit kinetic parameters from
  concentration-versus-time data.
- Apply the isolation (excess-of-reactants) method to decouple the reaction orders of two
  reactants.
- Apply the method of initial rates to extract $n$ and $k$ when reverse reactions or product
  inhibition would otherwise distort the analysis.
- Identify common pitfalls in kinetic experiments: side reactions, transport limitations,
  temperature gradients, analytical accuracy.
- Use the rank of the stoichiometric-coefficient matrix to identify a set of key species and key
  reactions, and reconstruct the changes in the non-key species via
  $\Delta\vect{n}_2 = \mtrx{N}_{2,1}\,\mtrx{N}_{1,1}^{-1}\,\Delta\vect{n}_1$.
:::

Until now we have assumed that the reaction rate can be written as

$$
r = k\, c\un{A}^{n} c\un{B}^{m} \dots
$$ (eq-power-law)

This is the power-law approximation, and the parameters $k$, $n$, $m$, … are typically determined
from experiments. Their values are set by Nature for a given reaction, but Nature does not tell us
what they are; we have to measure them in the lab.

:::{admonition} Discussion
:class: seealso
Show of hands: who has done experiments in a lab? And collected concentration-versus-time data?
:::

We now turn to the question of how to perform kinetic experiments and how to determine kinetic
parameters by regression.

**The reaction rate cannot be measured directly.** What is measured is the extent of reaction, or a
related observable, and the rate is then inferred from the data.

We focus on batch reactors here. Another common class is flow reactors — continuous stirred-tank or
plug-flow — but those have different material balances and are covered in ChBE 4320. Batch reactors
are widely used in practice for kinetic studies.

:::{admonition} Discussion
:class: seealso
Why are batch reactors so commonly used for kinetic studies?
:::

- They are cheap and easy to control.
- Only a small quantity of material is needed.
- They eliminate transport limitations — we want to see chemical kinetics, not mass-transfer
  artifacts.
- Their material balances are simple, which helps the analysis.

The overall procedure to determine reaction kinetics is as follows:

- Collect $\xi$ vs. $t$ data — recall that the reaction rate is
  $r = \frac{1}{V}\frac{\mathrm{d}\xi}{\mathrm{d}t}$.
- Fit a material balance to the data.
- Extract the rate constant $k$, the reaction orders, and the temperature dependence.

:::{admonition} Discussion
:class: seealso
Which methods can we use to measure the extent of reaction?

- Spectroscopy or mass spectrometry of specific components
- Titration or chromatography
- Following the change in a physical property, such as conductivity or refractive index
- Following pressure in constant-volume systems
- Following volume in constant-pressure systems
:::

Regardless of the method, there are a few rules that any kinetic measurement must follow:

- The measurement must not disturb the system under investigation.
- The measurement must be representative of the system as a whole.
- The method must provide a true measure of the extent of reaction.

Methods are classified as either chemical or physical.

**Chemical methods** draw a sample for analysis, which raises a few issues: drawing a sample can
change the volume of the system, and the reaction continues to proceed in the sample unless it is
rapidly quenched — chilled or chemically deactivated to stop the reaction. Very fast reactions
cannot be studied this way. Chemical methods typically do not give a continuous reading, since each
sample takes time to analyze, although individual measurements can be very fast.

**Physical methods** measure a property of the reacting system without removing material. Such
measurements can be performed *in situ* (Latin for "in place", meaning inside the reactor at
reaction conditions) or *in operando* ("while operating", meaning during the actual operation of
the reactor). The reaction vessel must be accessible to the probe, often optically, but the
measurement can then be performed continuously. The challenge is to find a physical property whose
change is large enough to track accurately.

The second step in extracting kinetic parameters is to fit a material balance to the data. There
are two classical ways to do this: the *integral method* and the *differential method*.

:::{figure} ../figures/ExpMethods.png
:label: fig-exp-methods
:alt: Hand-drawn schematic in three parts. On the left, the differential equation dc_A/dt equals a function of c_A, braced underneath and labelled integral method. An arrow points to a boxed centre panel labelled c_A versus t, braced and labelled experiment. A second arrow points to the same differential equation on the right, braced and labelled differential method. The integral method runs from the rate law to the data; the differential method runs from the data back to the rate law.
:width: 75%

Illustration of the integral and differential methods for determining kinetic parameters from
experiments.
:::

## Integral method

<!-- source: Experiments.tex L102 -->

In the integral method we postulate a rate law, integrate it analytically, and then fit the rate
constant and reaction orders so that the integrated form agrees with the data.

**Procedure.** Postulate a rate law, or take one from a mechanistic analysis as done in the
previous two chapters:

$$
r\un{A} = \nu\un{A} r = -k\, c\un{A}^{n} .
$$ (eq-integral-postulate)

Solve the batch-reactor mass balance for $c\un{A}$ vs. $t$:

$$
\begin{aligned}
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} &= -k\, c\un{A}^{n} \\
\int_{c\un{A,0}}^{c\un{A}} \frac{\mathrm{d}c\un{A}'}{(c\un{A}')^{n}} &= -k\int_0^t \mathrm{d}t' \\
\frac{1}{1-n}\left(c\un{A}^{1-n} - c\un{A,0}^{1-n}\right) &= -kt .
\end{aligned}
$$ (eq-integral-general)

[](#eq-integral-general) is only valid for $n \neq 1$; the case $n = 1$ gives the familiar
exponential derived earlier. Linearize the integrated form and plot the experimental data
accordingly:

$$
\begin{aligned}
c\un{A}^{1-n} &= c\un{A,0}^{1-n} - (1-n)\, k\, t \\
y &= b + m\, x .
\end{aligned}
$$ (eq-integral-linearized)

Try different values of $n$ until the data fall on a straight line. Once the reaction order is
identified, the rate constant follows from the slope.

:::{figure} ../figures/IntegralGraphMethod.png
:label: fig-integral
:alt: Hand-drawn plot with c_A to the power one minus n on the vertical axis against time on the horizontal axis. One set of data points falls on a straight descending line, annotated with y-intercept equal to c_A0 to the power one minus n and slope equal to minus one minus n times k. A second set of points, marked with crosses, curves away from the line and is annotated as indicating an incorrect rate law, so try again.
:width: 75%

Graphical illustration of the integral method.
:::

:::{admonition} Discussion
:class: seealso
What is the main practical drawback of the integral method?
:::

The integration is often difficult, and several reaction orders may need to be tested before the
data linearize, which is time-consuming.

## Differential method

<!-- source: Experiments.tex L145 -->

The second classical approach is the differential analysis of kinetic data. This method is very
useful for complex rate laws, where the integral form is hard to obtain analytically. The procedure
is:

- Determine $\mathrm{d}c\un{A}/\mathrm{d}t$ from the $c\un{A}$ vs. $t$ data.
- Postulate a rate law — or obtain it from a mechanistic analysis using the PSSA, QEA, or RDS
  assumption — linearize it, and plot the data appropriately.

:::{figure} ../figures/DiffMethod.png
:label: fig-differential
:alt: Hand-drawn plot of the concentration of A against time. Experimental points marked with crosses fall along a decaying curve, and short straight tangent lines are drawn through several of the points, annotated slope at each point.
:width: 55%

Reaction rates can be determined from concentration-versus-time data by differentiation.
:::

The derivative can be approximated by a finite-difference quotient. The species production rate is

$$
r_i = \frac{1}{\nu_i}\frac{\Delta c_i}{\Delta t} .
$$ (eq-species-rate-diff)

The simplest finite-difference scheme is the forward difference,

$$
f'(x) = \frac{f(x+h) - f(x)}{h} ,
$$ (eq-forward-difference)

and a more accurate alternative is the symmetric, or central, difference quotient,

$$
f'(x) = \frac{f(x+h) - f(x-h)}{2h} .
$$ (eq-central-difference)

Higher-order schemes that use additional data points — the five-point central-difference formula,
or the backward-difference formula — provide even better estimates.

Once the rate has been extracted from the experimental data, the reaction order and the rate
constant can be determined by linear regression, by nonlinear regression, or graphically. We start
with the graphical approach.

**Example.**

$$
\begin{aligned}
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} &= -k\, c\un{A}^n \\
\underbrace{\ln\left(-\frac{\mathrm{d}c\un{A}}{\mathrm{d}t}\right)}_{y}
&= \underbrace{n \ln(c\un{A})}_{mx} + \underbrace{\ln(k)}_{b} ,
\end{aligned}
$$ (eq-differential-linearized)

the linearized form. Using the differential method, both $n$ and $k$ can be extracted from a single
plot, whereas the integral method generally requires several attempts.

:::{figure} ../figures/DifferentialMethod.png
:label: fig-differential-linearized
:alt: Hand-drawn plot of the natural logarithm of minus dc_A by dt on the vertical axis against the natural logarithm of c_A on the horizontal axis. Scattered data points fall along a straight rising line, annotated with slope equal to n and y-intercept equal to the natural logarithm of k.
:width: 65%

Graphical illustration of the differential method.
:::

## Temperature dependence

<!-- source: Experiments.tex L204 -->

The temperature dependence of the rate constant is captured by the Arrhenius equation,

$$
k = A \exp\left(\frac{-E\un{a}}{RT}\right) .
$$ (eq-arrhenius-exp)

By performing experiments at several temperatures and extracting $k(T)$ from each, the Arrhenius
parameters can be determined by the same linearization strategy. Taking the logarithm,

$$
\ln(k) = -\frac{E\un{a}}{RT} + \ln(A) ,
$$ (eq-arrhenius-linear)

so a plot of $\ln k$ vs. $1/T$ has slope $-E\un{a}/R$ and intercept $\ln A$. If $k$ is known at only
two temperatures, the activation energy follows directly from

$$
\ln\left(\frac{k_2}{k_1}\right) = \frac{E\un{a}}{R}\left[\frac{1}{T_1} - \frac{1}{T_2}\right] .
$$ (eq-arrhenius-two-point)

## Nonlinear regression

<!-- source: Experiments.tex L229 -->

The procedures above are the classical, graphical approach. They predate modern computing and are
still very useful as sanity checks, because the linearization makes the parameter dependence
visually transparent. In a modern lab, however, parameters are typically extracted by nonlinear
regression directly against the ODE system.

For the simple reaction $\ce{A -> products}$ with first-order kinetics in A,

$$
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t}
= -A \exp\left(\frac{-E\un{a}}{RT}\right) c\un{A}^{\,n} .
$$ (eq-nonlinear-ode)

A residual function is defined that compares the simulated and measured concentrations, or rates,
at each experimental data point, and the parameters are found by minimizing the sum of squared
residuals:

$$
\min_{A,\,E\un{a},\,n} \sum_i
\underbrace{\left[r\un{exp,i} - r(A, E\un{a}, n)\right]^2}_{\text{objective function}} .
$$ (eq-nonlinear-objective)

This is the standard approach for complex rate expressions. A wide range of nonlinear optimizers is
available in MATLAB or Python — Levenberg–Marquardt, BFGS, and others. In most cases you will need
to provide parameter bounds and reasonable initial guesses, since the objective function often has
many local minima.

## Excess of reactants methods

<!-- source: Experiments.tex L254 -->

The earlier discussion focused on the simple case $\ce{A -> products}$, but for most reactions of
practical interest the rate depends on more than one reactant,

$$
\ce{A + B -> products} .
$$

For such cases the integral method can be used in principle, but exploring the joint $(n,m)$
parameter space by trial and error is tedious; the differential method, combined with a clever
choice of initial conditions, is more efficient.

Assume a power-law rate expression,

$$
-\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} = k\, c\un{A}^n\, c\un{B}^m .
$$ (eq-isolation-general)

The idea is to run the experiment with one reactant in large excess, $c\un{B,0} \gg c\un{A,0}$. Then
$c\un{B} \approx c\un{B,0}$ throughout the reaction, and the rate expression collapses to a
single-reactant form,

$$
-\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} = \underbrace{k'}_{k\, c\un{B,0}^m}\, c\un{A}^n ,
$$ (eq-isolation-a)

from which $n$ can be determined by the differential method. A second experiment with A in excess,
so that $c\un{A} \approx c\un{A,0}$, gives

$$
-\frac{\mathrm{d}c\un{B}}{\mathrm{d}t} = \underbrace{k''}_{k\, c\un{A,0}^n}\, c\un{B}^m ,
$$ (eq-isolation-b)

from which $m$ is obtained. Two experiments are therefore needed to determine the orders $n$ and
$m$ separately.

:::{figure} ../figures/Excess.png
:label: fig-excess
:alt: Hand-drawn plot of the natural logarithm of minus dc_A by dt on the vertical axis, against the natural logarithm of c_A and of c_B on the horizontal axis. Two straight lines of different slope are drawn through two separate sets of data points, the upper labelled c_A and the lower labelled c_B, annotated that the orders n and m are extracted from different experiments.
:width: 75%

Graphical illustration of the excess-of-reactant (isolation) method for determining the kinetic
parameters.
:::

With $n$ and $m$ in hand, the rate constant $k$ in the original two-reactant expression
$-r\un{A} = k\, c\un{A}^n\, c\un{B}^m$ follows from any single rate measurement at known $c\un{A}$
and $c\un{B}$:

$$
k = \frac{-r\un{A}}{c\un{A}^n\, c\un{B}^m} .
$$ (eq-isolation-k)

The same procedure extends to more complex rate laws. Linearization is not strictly required:
nonlinear regression on the original, un-linearized data typically gives better parameter estimates,
because linearization can place disproportionate weight on certain regions of the data and distort
the fit. Take care when the data span only a narrow range of concentration.

## Method of initial rates

<!-- source: Experiments.tex L305 -->

The differential method is attractive because a single experiment is in principle enough to extract
$k$ and $n$.

:::{admonition} Discussion
:class: seealso
What can go wrong when applying the differential method to a real reaction?
:::

The method becomes unreliable when the rate is significantly affected by the reverse reaction or by
product inhibition: the apparent rate then slows with increasing extent of reaction even when the
underlying forward kinetics is, say, first order. A good rule of thumb is to keep the kinetic
measurement away from equilibrium.

In the **method of initial rates**, the rate is evaluated only at the very start of the reaction, at
low conversion, where the reverse-reaction and product-inhibition contributions are negligible. At
the initial conditions,

$$
-r\un{A,0} = k\, c\un{A,0}^n .
$$ (eq-initial-rate)

Repeating the experiment with different initial concentrations and plotting
$\ln\left(-r\un{A,0}\right)$ vs. $\ln c\un{A,0}$ then yields the reaction order from the slope and
the rate constant from the intercept. Note the minus sign: $r\un{A,0}$ is the *production* rate of a
reactant and is therefore negative, so it is $-r\un{A,0}$ that we take the logarithm of.

:::{admonition} Caution
:class: warning
The methods in this chapter extract empirical rate-law parameters but cannot tell us what happens
at the atomic level. They do not verify a proposed reaction mechanism on their own. Discriminating
between candidate mechanisms generally requires complementary information, such as
electronic-structure calculations of the underlying elementary steps ([](#ch-microscopic)).
:::

## Experimental data

<!-- source: Experiments.tex L330 -->

The quality of the experimental data is critical. A few common pitfalls that recur in the
literature are worth keeping in mind.

- Make sure that the reaction being studied is actually the one of interest — watch for side
  reactions, reactions in series, or the reactor itself acting as a catalyst.
- Experimental conditions must be clearly specified and measured: reactant purity, temperature,
  catalyst mass, and so on.
- Minimize temperature effects. An isothermal reactor is the goal, but endothermic or exothermic
  reactions can produce internal gradients.
- Avoid mass-transfer limitations and concentration inhomogeneities.
- Verify the accuracy of the analytical method.
- Check reproducibility across multiple runs.

## Key reactions and key species

<!-- source: Experiments.tex L343 -->

**Example: steam reforming.** Steam reforming of methane is a major industrial process for the
production of synthesis gas, $\ce{CO}$ and $\ce{H2}$. It involves a complex reaction network:

$$
\begin{aligned}
\ce{CH4 + H2O &<=> CO + 3 H2} \\
\ce{CO + H2O &<=> CO2 + H2} \\
\ce{CH4 &<=> C + 2 H2} \\
\ce{2 CH4 &<=> C2H6 + H2} \\
\ce{C + H2O &<=> CO + H2} \\
\ce{2 CO &<=> C + CO2}
\end{aligned}
$$

A natural question is whether it is necessary to measure all of the species in the reaction mixture
in order to extract the reaction kinetics. The answer is no: in general only a small subset of
species needs to be measured, and the rest can be reconstructed from stoichiometry. The question
becomes which species to measure, and how many.

The matrix of stoichiometric coefficients $\mtrx{N}$ provides the answer. It has one column per
reaction and one row per species, with entries $\nu_{i,j}$. When the system has more reactions than
independent ones — that is, when the reactions are not all linearly independent — the rank of
$\mtrx{N}$ tells us how many species must be measured.

:::{table} The stoichiometric-coefficient matrix $\mtrx{N}$ for the steam-reforming network, with one column per reaction $j$ and one row per species.
:label: tab-n-matrix

| Species | 1 | 2 | 3 | 4 | 5 | 6 |
|:--------|--:|--:|--:|--:|--:|--:|
| $\ce{H2}$   |  3 |  1 |  2 |  1 |  1 |  0 |
| $\ce{CO}$   |  1 | -1 |  0 |  0 |  1 | -2 |
| $\ce{H2O}$  | -1 | -1 |  0 |  0 | -1 |  0 |
| $\ce{CH4}$  | -1 |  0 | -1 | -2 |  0 |  0 |
| $\ce{C}$    |  0 |  0 |  1 |  0 | -1 |  1 |
| $\ce{CO2}$  |  0 |  1 |  0 |  0 |  0 |  1 |
| $\ce{C2H6}$ |  0 |  0 |  0 |  1 |  0 |  0 |
:::

We can determine a set of **key species** and a set of **key reactions**, both fixed by the rank
$R_{\nu}$ of the matrix. To find out which reactions and species matter, use Gaussian elimination
to triangulate the matrix and read off its rank, which gives

:::{table} The stoichiometric matrix after Gaussian elimination. Three rows vanish, so the rank is 4.
:label: tab-n-matrix-triangulated

| Species | 1 | 2 | 3 | 4 | 5 | 6 |
|:--------|--:|--:|--:|--:|--:|--:|
| $\ce{H2}$   | 3 | 1              | 2              | 1              | 1  | 0  |
| $\ce{CO}$   | 0 | $-\frac{4}{3}$ | $-\frac{2}{3}$ | $-\frac{1}{3}$ | $\frac{2}{3}$ | -2 |
| $\ce{H2O}$  | 0 | 0              | 1              | $\frac{1}{2}$  | -1 | 1  |
| $\ce{CH4}$  | 0 | 0              | 0              | $-\frac{3}{2}$ | 0  | 0  |
| $\ce{C}$    | 0 | 0              | 0              | 0              | 0  | 0  |
| $\ce{CO2}$  | 0 | 0              | 0              | 0              | 0  | 0  |
| $\ce{C2H6}$ | 0 | 0              | 0              | 0              | 0  | 0  |
:::

For this system the rank is 4, which means that

- only the first four reactions are linearly independent, and
- only the conversion of $\ce{H2}$, $\ce{CO}$, $\ce{H2O}$, and $\ce{CH4}$ needs to be measured. The
  rest can be determined from the stoichiometry.

Using the rank, the matrix of stoichiometric coefficients can be organized so that the key species
and key reactions form the first $R_{\nu}$ rows and columns. The matrix then divides into four
sub-matrices. Similarly, we can construct a vector containing the change in the number of moles and
partition it into key and non-key species.

:::{figure} ../figures/matrix.png
:label: fig-n-partition
:alt: Diagram of a large rectangle representing the matrix N, divided by one vertical and one horizontal line into four blocks labelled N_1,1 top left, N_1,2 top right, N_2,1 bottom left and N_2,2 bottom right. The columns are indexed 1 to R_v to M and the rows 1 to R_v to N, so the top-left block holds the key species and key reactions. Beside it a tall narrow rectangle representing the vector n is divided horizontally into an upper block n_1 and a lower block n_2.
:width: 60%

Partitioning of the stoichiometric-coefficient matrix into key and non-key species and reactions.
:::

If the rows and columns are sorted appropriately, we can invert $\mtrx{N}_{1,1}$:

$$
\mtrx{N}_{1,1}^{-1} = \begin{pmatrix}
0 & \frac{1}{2} & -\frac{1}{2} & 0 \\
0 & -\frac{1}{2} & -\frac{1}{2} & 0 \\
\frac{2}{3} & -\frac{1}{2} & \frac{7}{6} & \frac{1}{3} \\
-\frac{1}{3} & 0 & -\frac{1}{3} & -\frac{2}{3}
\end{pmatrix} .
$$

We can then calculate the vector of the change in molar amounts of the species,

$$
\Delta\vect{n}_1 = \mtrx{N}_{1,1}\, \vect{\xi} , \qquad
\Delta\vect{n}_2 = \mtrx{N}_{2,1}\, \vect{\xi} ,
$$

and determine the extent-of-reaction vector $\vect{\xi}$ from the known $\Delta\vect{n}_1$,

$$
\vect{\xi} = \mtrx{N}_{1,1}^{-1}\, \Delta\vect{n}_1 .
$$

Thus we can determine the conversion of the non-key species via

$$
\Delta\vect{n}_2 = \mtrx{N}_{2,1}\, \mtrx{N}_{1,1}^{-1}\, \Delta\vect{n}_1 ,
$$ (eq-nonkey)

where

$$
\mtrx{N}_{2,1} = \begin{pmatrix}
0 & 0 & 1 & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}
\quad \rightarrow \quad
\mtrx{N}_{2,1}\mtrx{N}_{1,1}^{-1} = \begin{pmatrix}
\frac{2}{3} & -\frac{1}{2} & \frac{7}{6} & \frac{1}{3} \\
0 & -\frac{1}{2} & -\frac{1}{2} & 0 \\
-\frac{1}{3} & 0 & -\frac{1}{3} & -\frac{2}{3}
\end{pmatrix} .
$$

Finally, we can determine the amounts of the non-key species:

$$
\begin{aligned}
\Delta n_{\ce{C}} &= \tfrac{2}{3}\Delta n_{\ce{H2}} - \tfrac{1}{2}\Delta n_{\ce{CO}}
  + \tfrac{7}{6}\Delta n_{\ce{H2O}} + \tfrac{1}{3}\Delta n_{\ce{CH4}} \\
\Delta n_{\ce{CO2}} &= 0\,\Delta n_{\ce{H2}} - \tfrac{1}{2}\Delta n_{\ce{CO}}
  - \tfrac{1}{2}\Delta n_{\ce{H2O}} + 0\,\Delta n_{\ce{CH4}} \\
\Delta n_{\ce{C2H6}} &= -\tfrac{1}{3}\Delta n_{\ce{H2}} - 0\,\Delta n_{\ce{CO}}
  - \tfrac{1}{3}\Delta n_{\ce{H2O}} - \tfrac{2}{3}\Delta n_{\ce{CH4}}
\end{aligned}
$$

This method is particularly useful when some species are difficult to detect: choose the key species
to be the ones that are easy and accurate to measure, and reconstruct the rest. The method does
require knowledge of all reactions that can proceed; if a relevant reaction is missing from
$\mtrx{N}$, the wrong species can be selected as key.

An alternative method that does not require this knowledge uses the **element-species matrix**,
constructed from the elemental composition of each species. The key components are then chosen to be
those that are easiest to detect, and conservation of elements provides the additional constraints
needed to reconstruct the rest.

## Summary

<!-- source: Experiments.tex L488 -->

- Kinetic experiments measure the extent of reaction, or a related observable, not the rate
  directly. Methods are classified as chemical (sampling) or physical (*in situ* / *in operando*).
- Integral method: postulate $r = -k c\un{A}^n$, integrate to get [](#eq-integral-general),
  linearize via [](#eq-integral-linearized), and try values of $n$ until the data fall on a straight
  line.
- Differential method: estimate $\mathrm{d}c\un{A}/\mathrm{d}t$ using forward,
  [](#eq-forward-difference), or central, [](#eq-central-difference), differences, then linearize
  $\ln(-\mathrm{d}c\un{A}/\mathrm{d}t)$ vs. $\ln c\un{A}$, [](#eq-differential-linearized), to read
  off both $n$ and $k$ from a single plot.
- Arrhenius linearization, [](#eq-arrhenius-linear), gives $E\un{a}$ from the slope of $\ln k$ vs.
  $1/T$; for two temperatures the activation energy follows directly from
  [](#eq-arrhenius-two-point).
- Modern practice is nonlinear regression directly on the ODE system,
  [](#eq-nonlinear-objective), using Levenberg–Marquardt, BFGS, or similar. Provide bounds and good
  initial guesses; expect multiple local minima.
- Isolation method: with $c\un{B,0} \gg c\un{A,0}$ the rate law collapses to single-reactant form,
  [](#eq-isolation-a), so $n$ is recovered by the differential method; a second experiment with
  $c\un{A,0} \gg c\un{B,0}$, [](#eq-isolation-b), recovers $m$.
- Method of initial rates: evaluate $-r\un{A,0}$ at low conversion, [](#eq-initial-rate), to avoid
  distortion from reverse reactions and product inhibition.
- Pitfalls in experimental data: side reactions, transport limitations, temperature gradients,
  analytical accuracy, reproducibility.
- Empirical rate-law fitting cannot verify a reaction mechanism; complementary tools —
  electronic-structure calculations, isotope labeling, spectroscopy — are needed for mechanistic
  discrimination ([](#ch-microscopic), [](#ch-mechanisms)).
- For a reaction network with more reactions than independent ones, the rank $R_\nu$ of the
  stoichiometric matrix sets the number of key species that must be measured. Non-key species follow
  from $\Delta\vect{n}_2 = \mtrx{N}_{2,1}\,\mtrx{N}_{1,1}^{-1}\,\Delta\vect{n}_1$, [](#eq-nonkey).
