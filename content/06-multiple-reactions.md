---
title: Kinetics of Multiple Reactions
short_title: Multiple reactions
label: ch-multiple-reactions
---

<!-- LaTeX source: multiple_reactions.tex -->
<!-- Porting notes: mhchem is math-only here, so \ce{} in prose needs $...$; and \un{} expands to
     _{\textrm{...}}, so species subscripts must be written _{\ce{Br2}}, never \un{\ce{Br2}}.
     Radical dots use radicals use \ce{Br\bullet}, which puts the dot on the math axis;
     \ce{Br^{.}} sets it as a superscript and \ce{Br{.}} drops it to the baseline.
     Nested directives need a longer outer fence (:::: around :::). -->

:::{admonition} Learning objectives
:class: tip
After completing this chapter, you should be able to:

- Solve the batch-reactor mass balances for two parallel first-order reactions and recognize that the product selectivity reduces to the ratio of rate constants.
- Define and compute the conversion $X$, the selectivity $S$, and the yield $Y$
- Apply the pseudo-steady-state approximation (PSSA) and the quasi-equilibrium approximation (QEA), recognize the timescale conditions under which each is valid, and use them to reduce an ODE system to closed-form rate expressions.
- Apply the PSSA to a radical chain mechanism (HBr synthesis) to derive an apparent rate law and identify the apparent rate constant and apparent activation energy.
- Apply the QEA to a pre-equilibrium mechanism ($\ce{NO + Br2}$) to derive a closed-form rate expression in terms of the equilibrium constant and a slow-step rate constant.
- Apply the rate-determining step (RDS) assumption and recognize when it can be invoked safely.
:::

We have looked into various simple reaction systems, considering only a single irreversible or reversible reaction. Unfortunately chemistry is not that simple, and we have to deal with much more complex systems in which multiple reactions occur simultaneously. We can broadly distinguish four types of reaction system: parallel, series, independent, and complex.

In a system with **parallel reactions**, a reactant is converted to two possible products by different pathways:

$$
\ce{A ->[$k_1$] B} \qquad \ce{A ->[$k_2$] C} .
$$

An example is the epoxidation of ethylene to ethylene oxide, a key intermediate in the production of ethylene glycol,

$$
\ce{C2H4 + 0.5 O2 -> C2H4O} ,
$$

where it is necessary to avoid the complete combustion of ethylene to $\ce{CO2}$,

$$
\ce{C2H4 + 3 O2 -> 2 CO2 + 2 H2O} .
$$

In addition to parallel reactions we frequently have **reactions in series**. The second step may be an unwanted reaction of our product to an undesired side product, or a necessary intermediate that is later converted to the target product:

$$
\ce{A -> B -> C} .
$$

If reactions occur at the same time but share no reactants or products, they are **independent**:

$$
\ce{A -> B + C} \qquad \ce{D -> E + F} .
$$

The cracking of crude oil is an example of independent reactions,

$$
\begin{aligned}
\ce{C15H32 &-> C12H26 + C3H6} \\
\ce{C8H18 &-> C6H14 + C2H4} .
\end{aligned}
$$

Solving the material balance for independent reactions is straightforward: we can consider each reaction in isolation and apply the rules derived before. By contrast, when two reactions *do* share a reactant or product they are **coupled**, or *dependent*, and the mass balances must be solved together. A reversible reaction $\ce{A <=> B}$ is the simplest example: the forward and reverse steps share both species and cannot be analyzed in isolation.

The last category is the **complex** reaction system, combining the previous types: multiple parallel reactions, reactions in series, and possibly some independent reactions. An example is the partial oxidation of n-butane ($\ce{C4H10}$) to maleic anhydride ($\ce{C4H2O3}$) over a vanadium phosphate oxide (VPO) catalyst, a major industrial chemical produced at roughly 3 million tonnes per year. Maleic anhydride is used as an intermediate to make polyester resins, copolymers, and lubricants.

:::{figure} ../figures/MA_system.png
:label: fig-ma-network
:alt: Reaction network diagram. n-butane at top left has three outgoing arrows: one straight across the top to MA, maleic anhydride, and two descending to AcA, acetic acid, and AcrA, acrylic acid. AcA and AcrA each have two outgoing arrows, to CO and to CO2, crossing over one another. MA at top right also has two arrows down to CO and CO2. Below the network is the skeletal structure of maleic anhydride, a five-membered ring with one ring oxygen, a carbon-carbon double bond, and two carbonyl oxygens.
:width: 60%

Complex reaction network of maleic anhydride (MA) synthesis from n-butane by partial oxidation, which leads to various side products and undesired parallel products, including acrylic acid (AcrA) and acetic acid (AcA).
:::

## Parallel reactions, irreversible

<!-- source: multiple_reactions.tex L70 -->

For some multi-reaction systems it is possible to derive analytical solutions to the material balances of the batch reactor. Consider a simple network of two parallel reactions with first-order kinetics,

$$
r_1 = k_1 c\un{A} , \qquad r_2 = k_2 c\un{A} ,
$$

which has the mass balances

$$
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} = -k_1 c\un{A} - k_2 c\un{A} , \qquad
\frac{\mathrm{d}c\un{B}}{\mathrm{d}t} = k_1 c\un{A} , \qquad
\frac{\mathrm{d}c\un{C}}{\mathrm{d}t} = k_2 c\un{A} .
$$

The concentration function of A is straightforward,

$$
\frac{c\un{A}}{c\un{A,0}} = \exp\left[-(k_1 + k_2)t\right] .
$$ (eq-parallel-ca)

Inserting this into the differential equation for B, separating the variables, and integrating gives

$$
\frac{c\un{B}}{c\un{A,0}} = \frac{c\un{B,0}}{c\un{A,0}}
+ \frac{k_1}{k_1 + k_2}\left\{1 - \exp\left[-(k_1 + k_2)t\right]\right\} ,
$$ (eq-parallel-cb)

and the same procedure for C yields

$$
\frac{c\un{C}}{c\un{A,0}} = \frac{c\un{C,0}}{c\un{A,0}}
+ \frac{k_2}{k_1 + k_2}\left\{1 - \exp\left[-(k_1 + k_2)t\right]\right\} .
$$ (eq-parallel-cc)

:::{figure} ../figures/CprofileParallel.png
:label: fig-parallel-reactions
:alt: Hand-drawn plot of concentration against time for two parallel reactions with k1 equal to twice k2. The curve for A starts at its initial value and decays to zero. The curves for B and C both rise from zero and level off, with B levelling off at roughly twice the plateau value of C.
:width: 70%

Concentration profiles for two parallel first-order reactions $\ce{A -> B}$ and $\ce{A -> C}$ in a constant-volume batch reactor.
:::

Taking the ratio of [](#eq-parallel-cb) and [](#eq-parallel-cc) gives a simple measure of how much B is produced relative to C:

$$
\frac{c\un{B} - c\un{B,0}}{c\un{C} - c\un{C,0}} = \frac{k_1}{k_2} .
$$ (eq-parallel-ratio)

The product ratio depends only on the ratio of rate constants, not on $c\un{A,0}$ or on time. We revisit this result in the next section, where selectivity is defined.

:::{admonition} Live example
:class: seealso
Drag the sliders to change the two rate constants:
the plateaus of B and C move with $k_1/(k_1 + k_2)$ and $k_2/(k_1 + k_2)$, but their ratio is always $k_1/k_2$. The MATLAB code that produces the same result follows underneath.
:::

```{marimo} python
import marimo as mo
import numpy as np
import matplotlib.pyplot as plt
```

```{marimo} python
par_k1 = mo.ui.slider(
    start=0.2,
    stop=5,
    step=0.1,
    value=2.2,
    label="Rate constant k1 (1/s)",
    show_value=True,
)
par_k2 = mo.ui.slider(
    start=0.2,
    stop=5,
    step=0.1,
    value=1.0,
    label="Rate constant k2 (1/s)",
    show_value=True,
)
mo.vstack([par_k1, par_k2])
```

```{marimo} python
par_cA0 = 1.0  # mol/m^3
par_cB0 = 0.0  # mol/m^3
par_cC0 = 0.0  # mol/m^3

par_ksum = par_k1.value + par_k2.value
par_time = np.linspace(0, 5, 200)

# Integrated material balances for A, B, and C
par_cA = par_cA0 * np.exp(-par_ksum * par_time)
par_cB = par_cB0 + par_k1.value / par_ksum * (1 - np.exp(-par_ksum * par_time)) * par_cA0
par_cC = par_cC0 + par_k2.value / par_ksum * (1 - np.exp(-par_ksum * par_time)) * par_cA0

mo.md(
    f"$k_1/k_2 = {par_k1.value / par_k2.value:.2f}$. "
    f"At long times $c_\\mathrm{{B}} \\to {par_k1.value / par_ksum * par_cA0:.3f}$ and "
    f"$c_\\mathrm{{C}} \\to {par_k2.value / par_ksum * par_cA0:.3f}\\ \\mathrm{{mol\\,m^{{-3}}}}$; "
    f"99 % of A is consumed after $t = {np.log(100) / par_ksum:.2f}\\ \\mathrm{{s}}$."
)
```

```{marimo} python
fig_par, ax_par = plt.subplots(figsize=(5.5, 3.9))

ax_par.plot(par_time, par_cA, lw=2.0, label="A")
ax_par.plot(par_time, par_cB, lw=2.0, label="B")
ax_par.plot(par_time, par_cC, lw=2.0, label="C")

ax_par.set_xlim(0, 5)
ax_par.set_ylim(0, 1)
ax_par.set_xlabel("time (s)", fontsize=13)
ax_par.set_ylabel(r"concentration (mol m$^{-3}$)", fontsize=13)
ax_par.tick_params(labelsize=12, width=1.2, length=5)
for sp_par in ax_par.spines.values():
    sp_par.set_linewidth(1.2)
ax_par.legend(loc="center right", fontsize=13, frameon=False)

fig_par.tight_layout()
fig_par
```

In MATLAB, set `k1` and `k2` to the values you want and run:

```matlab
cA0 = 1;  % mol/m3
cB0 = 0;  % mol/m3
cC0 = 0;  % mol/m3

k1 = 2.2;  % 1/s
k2 = 1;    % 1/s

% Integrated material balances for the species A, B, and C
cA = @(time) ...
    cA0*exp(-(k1+k2)*time);

cB = @(time) ...
    cB0 + k1/(k1+k2)*(1-exp(-(k1+k2)*time))*cA0;

cC = @(time) ...
    cC0 + k2/(k1+k2)*(1-exp(-(k1+k2)*time))*cA0;

time = linspace(0, 5, 100);

figure('Units','centimeters','Position',[5 5 14 10])
hold on

plot(time, cA(time), 'LineWidth', 2.0, 'LineStyle','-','DisplayName','A')
plot(time, cB(time), 'LineWidth', 2.0, 'LineStyle','-','DisplayName','B')
plot(time, cC(time), 'LineWidth', 2.0, 'LineStyle','-','DisplayName','C')

xlim([0 5])
ylim([0 1])

xlabel('$\mathrm{time\ (s)}$','Interpreter','latex','FontSize',16)
ylabel('$\mathrm{concentration\ \left(mol\,m^{-3}\right)}$','Interpreter','latex','FontSize',16)

set(gca, ...
    'FontName','lmodern', ...
    'FontSize',16, ...
    'LineWidth',1.5, ...
    'TickLength',[0.015 0.015], ...
    'Box','on')

legend('Interpreter','latex','Location','east','FontSize',16)

grid off
hold off
```

:::{tip} Check yourself
With $k_1 = 2.2\ \mathrm{s^{-1}}$ and $k_2 = 1\ \mathrm{s^{-1}}$ the plateaus are
$c\un{B} = 0.688$ and $c\un{C} = 0.313\ \mathrm{mol\,m^{-3}}$, and their ratio is 2.2 at every instant. Check both against [](#eq-parallel-cb), [](#eq-parallel-cc), and [](#eq-parallel-ratio).
:::

## Selectivity and yield

<!-- source: multiple_reactions.tex L123 -->

In a system with multiple reactions, achieving a high conversion of the reactant is no longer
enough. Recall that the conversion of a reactant $i$ and the unconverted fraction $f$ were defined
as

$$
X = \frac{n\un{i,0} - n_i}{n\un{i,0}} , \qquad f = 1 - X = \frac{n_i}{n\un{i,0}} .
$$ (eq-conversion-def)

A high conversion is necessary, but what we really want is for the converted reactant to end up in
the desired product. If side products are formed, downstream separations are needed, which directly
affects the economics of the process. We therefore try to boost the production of the target
product while suppressing side products.

The quantity that captures how a converted reactant is partitioned among the possible products is
the **selectivity**. Consider a reactant A that can be converted into either B or C. The
selectivity of producing B from A is defined as

$$
S\un{B,A} = \frac{n\un{B} - n\un{B,0}}{n\un{A} - n\un{A,0}}\, \frac{\nu\un{A}}{\nu\un{B}} .
$$ (eq-selectivity-def)

The ratio of moles is already dimensionless; the stoichiometric factor $\nu\un{A}/\nu\un{B}$ is what
normalizes $S\un{B,A}$ so that it is positive and bounded between 0 and 1. The same definition
extends to any number of products formed from a common reactant.

**Example.** Consider methanol synthesis from synthesis gas,

$$
\ce{CO + 2 H2 -> CH3OH} ,
$$

with the side reactions

$$
\begin{aligned}
\ce{2 CO + 4 H2 &-> CH3OCH3 + H2O} \\
\ce{3 CO + 4 H2 &-> CH3COOCH3 + H2O}
\end{aligned}
$$

forming dimethyl ether ($\ce{CH3OCH3}$) and methyl acetate ($\ce{CH3COOCH3}$). The following molar
amounts were measured after a reaction time $t_1$.

:::{table} Molar amounts (in mol) measured for the methanol-synthesis example, used to calculate selectivity.
:label: tab-selectivity-example

| $t$   | $\ce{CO}$ | $\ce{H2}$ | $\ce{CH3OH}$ | $\ce{CH3OCH3}$ | $\ce{CH3COOCH3}$ | $\ce{H2O}$ |
|:-----:|:---------:|:---------:|:------------:|:--------------:|:----------------:|:----------:|
| 0     | 1         | 2.2       | 0            | 0              | 0                | 0          |
| $t_1$ | 0.78      | 1.8       | 0.1          | 0.03           | 0.02             | 0.05       |
:::

From this we can calculate the selectivity of forming each product from the reactant — for
instance, the selectivity of forming methanol from CO:

$$
\begin{aligned}
S\un{CH3OH,CO} &= \frac{0.1\ \mathrm{mol} - 0}{0.78\ \mathrm{mol} - 1\ \mathrm{mol}}\,\frac{-1}{1}
= 0.455 \\
S\un{CH3OCH3,CO} &= \frac{0.03\ \mathrm{mol} - 0}{0.78\ \mathrm{mol} - 1\ \mathrm{mol}}\,\frac{-2}{1}
= 0.273 \\
S\un{CH3COOCH3,CO} &= \frac{0.02\ \mathrm{mol} - 0}{0.78\ \mathrm{mol} - 1\ \mathrm{mol}}\,\frac{-3}{1}
= 0.273
\end{aligned}
$$

The sum over all selectivities equals 1 and indicates how the converted reactant is partitioned
among the products.

For systems where the rates rather than the integrated amounts are known, it is more convenient to
work with the **point selectivity**, or instantaneous selectivity, defined as the ratio of the
production rates of two species $i$ and $j$:

$$
S\un{ij} = \frac{r_i}{r_j} .
$$ (eq-point-selectivity)

The integral selectivity in [](#eq-selectivity-def) characterizes the cumulative behavior over a
finite reaction time and is the natural quantity when batch yields are reported; the point
selectivity in [](#eq-point-selectivity) characterizes the local behavior at a given concentration
and is the natural quantity when designing or comparing reactor configurations.

Another important quantity is the **yield** $Y$, defined as

$$
Y\un{B} = \frac{n\un{B,0} - n\un{B}}{n\un{A,0}}\, \frac{\nu\un{A}}{\nu\un{B}} ,
$$ (eq-yield-def)

$$
Y\un{B} = X\un{A} \cdot S\un{B,A} .
$$ (eq-yield-xs)

The yield is the fraction of the initial reactant A that has been converted into product B. For the
methanol example,

$$
Y\un{CH3OH} = \frac{0 - 0.1\ \mathrm{mol}}{1\ \mathrm{mol}}\,\frac{-1}{1} = 0.1 ,
$$

which equals the product of selectivity and conversion, as [](#eq-yield-xs) states. Analogously, an
**instantaneous yield** can be defined as

$$
Y\un{B} = \frac{r\un{B}}{r\un{A}}\, \frac{\nu\un{A}}{\nu\un{B}} .
$$ (eq-inst-yield)

## Reactions in series, irreversible

<!-- source: multiple_reactions.tex L206 -->

Consider the reaction

$$
\ce{A -> B -> C} .
$$

:::{admonition} Discussion
:class: seealso
What is our intuition for the shape of the concentration profiles? What do they depend on?
:::

:::{figure} ../figures/SketchSeriesReactions.png
:label: fig-reactions-in-series
:alt: Hand-drawn plot of concentration against time for a series reaction. The curve for A starts at its initial value and decays towards zero. The curve for B rises from zero, passes through a clear maximum at intermediate time, and then decays back towards zero. The curve for C starts at zero, stays flat during an induction period, then rises with an inflection and levels off near the initial concentration of A.
:width: 75%

Sketch of the concentration profiles for reactions in series.
:::

There are many examples of this type of reaction sequence. A series reaction poses a design
challenge, because the intermediate B is itself reacted away to C, so we have to optimize the time
at which we stop the reaction in order to maximize B. We therefore need to understand reactions in
series quantitatively.

Let the reaction rates be given by simple first-order expressions in the corresponding reactants,

$$
r_1 = k_1 c\un{A} ,
$$ (eq-series-r1)

$$
r_2 = k_2 c\un{B} .
$$ (eq-series-r2)

The material balances are

$$
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} = \nu\un{A} r_1 = -k_1 c\un{A} ,
$$ (eq-series-dca)

$$
\frac{\mathrm{d}c\un{B}}{\mathrm{d}t} = \nu\un{B,1} r_1 + \nu\un{B,2} r_2 = k_1 c\un{A} - k_2 c\un{B} ,
$$ (eq-series-dcb)

$$
\frac{\mathrm{d}c\un{C}}{\mathrm{d}t} = \nu\un{C} r_2 = k_2 c\un{B} .
$$ (eq-series-dcc)

The solution for A is the familiar first-order decay,

$$
c\un{A} = c\un{A,0}\exp\left(-k_1 t\right) .
$$ (eq-series-ca)

Inserting this into the material balance for B gives

$$
\frac{\mathrm{d}c\un{B}}{\mathrm{d}t} + k_2 c\un{B} = k_1 c\un{A,0}\exp\left(-k_1 t\right) .
$$ (eq-series-dcb-nonsep)

The variables can no longer be separated, so we need a more general technique. The standard tool
for first-order linear ODEs is the **integrating-factor method**, which applies to any equation of
the form

$$
y' + p(t)\, y = f(t) ,
$$

where here $y = c\un{B}(t)$. Comparing with [](#eq-series-dcb-nonsep),

$$
p(t) = k_2 , \qquad f(t) = k_1 c\un{A,0}\exp(-k_1 t) .
$$

We multiply both sides by an integrating factor $I(t)$,

$$
I\, y' + I\, p(t)\, y = I\, f(t) ,
$$

chosen so that the left-hand side is the derivative of $I y$. This is achieved by

$$
I(t) = \exp\int p(t)\, \mathrm{d}t = \exp(k_2 t) .
$$ (eq-integrating-factor)

Inserting this,

$$
\exp(k_2 t)\, y' + k_2 \exp(k_2 t)\, y = \exp(k_2 t)\, k_1 c\un{A,0}\exp(-k_1 t) .
$$

The left-hand side is now $\mathrm{d}\!\left[\exp(k_2 t)\, y\right]/\mathrm{d}t$ — use the product
rule to confirm — so we can integrate directly:

$$
\begin{aligned}
\int \frac{\mathrm{d}\!\left[\exp(k_2 t)\, y\right]}{\mathrm{d}t}\, \mathrm{d}t
&= \int k_1 c\un{A,0}\exp\left[(k_2 - k_1)t\right] \mathrm{d}t \\
\exp(k_2 t)\, y &= \frac{k_1 c\un{A,0}}{k_2 - k_1}\exp\left[(k_2 - k_1)t\right] + C .
\end{aligned}
$$

The constant of integration $C$ is pinned by the initial condition. Dividing through by
$\exp(k_2 t)$ and rearranging,

$$
c\un{B} = \exp(-k_2 t)\left[\frac{k_1 c\un{A,0}}{k_2 - k_1}\exp\left[(k_2 - k_1)t\right] + C\right] .
$$

Applying the initial condition $c\un{B}(0) = c\un{B,0}$,

$$
\begin{aligned}
c\un{B,0} &= \frac{k_1 c\un{A,0}}{k_2 - k_1} + C \\
C &= c\un{B,0} - \frac{k_1 c\un{A,0}}{k_2 - k_1} ,
\end{aligned}
$$

and inserting and rearranging leads to

$$
c\un{B} = \frac{k_1 c\un{A,0}}{k_2 - k_1}
\left[\exp(-k_1 t) - \exp(-k_2 t)\right] + c\un{B,0}\exp(-k_2 t) .
$$ (eq-series-cb)

:::{admonition} Discussion
:class: seealso
Does [](#eq-series-cb) make sense? Check the limits at $t = 0$ and $t \to \infty$.
:::

For species C we could integrate $\mathrm{d}c\un{C}/\mathrm{d}t = k_2 c\un{B}$ directly, but that is
tedious. A cleaner shortcut uses the overall species balance: summing [](#eq-series-dca),
[](#eq-series-dcb), and [](#eq-series-dcc) gives

$$
\begin{aligned}
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} + \frac{\mathrm{d}c\un{B}}{\mathrm{d}t}
+ \frac{\mathrm{d}c\un{C}}{\mathrm{d}t} &= 0 \\
\frac{\mathrm{d}\left(c\un{A} + c\un{B} + c\un{C}\right)}{\mathrm{d}t} &= 0 \\
c\un{A} + c\un{B} + c\un{C} &= c\un{A,0} + c\un{B,0} + c\un{C,0} .
\end{aligned}
$$

Therefore

$$
\begin{aligned}
c\un{C} &= c\un{A,0} + c\un{B,0} + c\un{C,0} - c\un{A} - c\un{B} \\
c\un{C} &= c\un{A,0} + c\un{B,0} + c\un{C,0} - c\un{A,0}\exp\left(-k_1 t\right)
- \frac{k_1 c\un{A,0}}{k_2 - k_1}\left[\exp(-k_1 t) - \exp(-k_2 t)\right]
- c\un{B,0}\exp(-k_2 t) .
\end{aligned}
$$

For $c\un{B,0} = c\un{C,0} = 0$ this reduces to

$$
c\un{C} = c\un{A,0}\left(1 - \frac{k_2}{k_2 - k_1}\exp\left(-k_1 t\right)
+ \frac{k_1}{k_2 - k_1}\exp\left(-k_2 t\right)\right) .
$$ (eq-series-cc)

:::{admonition} Live example
:class: seealso
Drag the sliders to change $k_1$ and $k_2$ and watch the maximum of B move. Setting $k_2$ equal to
$k_1$ is allowed — the expressions above are then indeterminate, and the code switches to the
limiting form $c\un{B} = c\un{A,0}\, k_1 t\, \exp(-k_1 t)$. The MATLAB code that produces the same
result follows underneath.
:::

```{marimo} python
ser_k1 = mo.ui.slider(
    start=0.2,
    stop=5,
    step=0.1,
    value=1.0,
    label="Rate constant k1 (1/s)",
    show_value=True,
)
ser_k2 = mo.ui.slider(
    start=0.2,
    stop=10,
    step=0.1,
    value=2.0,
    label="Rate constant k2 (1/s)",
    show_value=True,
)
mo.vstack([ser_k1, ser_k2])
```

```{marimo} python
ser_cA0 = 1.0  # mol/m^3
ser_cB0 = 0.0  # mol/m^3
ser_k1v = ser_k1.value
ser_k2v = ser_k2.value
ser_time = np.linspace(0, 5, 300)

# Integrated material balances for A, B, and C
ser_cA = ser_cA0 * np.exp(-ser_k1v * ser_time)

if np.isclose(ser_k1v, ser_k2v):
    ser_cB = ser_cA0 * ser_k1v * ser_time * np.exp(-ser_k2v * ser_time)
    ser_t_max = 1 / ser_k1v
else:
    ser_cB = ser_k1v * ser_cA0 / (ser_k2v - ser_k1v) * (
        np.exp(-ser_k1v * ser_time) - np.exp(-ser_k2v * ser_time)
    ) + ser_cB0 * np.exp(-ser_k2v * ser_time)
    ser_t_max = np.log(ser_k2v / ser_k1v) / (ser_k2v - ser_k1v)

ser_cC = ser_cA0 + ser_cB0 - ser_cA - ser_cB

if np.isclose(ser_k1v, ser_k2v):
    ser_cB_max = ser_cA0 * ser_k1v * ser_t_max * np.exp(-ser_k2v * ser_t_max)
else:
    ser_cB_max = ser_k1v * ser_cA0 / (ser_k2v - ser_k1v) * (
        np.exp(-ser_k1v * ser_t_max) - np.exp(-ser_k2v * ser_t_max)
    )

mo.md(
    f"$k_2/k_1 = {ser_k2v / ser_k1v:.2f}$. "
    f"B peaks at $t_\\mathrm{{max}} = {ser_t_max:.3f}\\ \\mathrm{{s}}$ with "
    f"$c_\\mathrm{{B,max}} = {ser_cB_max:.3f}\\ \\mathrm{{mol\\,m^{{-3}}}}$, "
    f"i.e. {100 * ser_cB_max / ser_cA0:.1f} % of the initial A."
)
```

```{marimo} python
fig_ser, ax_ser = plt.subplots(figsize=(5.5, 3.9))

ax_ser.plot(ser_time, ser_cA, lw=2.0, label="A")
ax_ser.plot(ser_time, ser_cB, lw=2.0, label="B")
ax_ser.plot(ser_time, ser_cC, lw=2.0, label="C")
ax_ser.axvline(ser_t_max, color="0.6", lw=1.0, ls=":")

ax_ser.set_xlim(0, 5)
ax_ser.set_ylim(0, 1)
ax_ser.set_xlabel("time (s)", fontsize=13)
ax_ser.set_ylabel(r"concentration (mol m$^{-3}$)", fontsize=13)
ax_ser.tick_params(labelsize=12, width=1.2, length=5)
for sp_ser in ax_ser.spines.values():
    sp_ser.set_linewidth(1.2)
ax_ser.legend(loc="center right", fontsize=13, frameon=False)

fig_ser.tight_layout()
fig_ser
```

In MATLAB, set `k1` and `k2` to the values you want and run:

```matlab
cA0 = 1;  % mol/m3
cB0 = 0;  % mol/m3
cC0 = 0;  % mol/m3

k1 = 1;  % 1/s
k2 = 2;  % 1/s

% Integrated material balances for the species A, B, and C
cA = @(time) ...
    cA0.*exp(-k1.*time);

if k1 == k2
    cB = @(time) ...
        cA0.*k1.*time.*exp(-k2.*time);
else
    cB = @(time) ...
        k1.*cA0/(k2-k1).*(exp(-k1.*time)-exp(-k2.*time)) + cB0.*exp(-k2.*time);
end

cC = @(time) ...
    cA0 + cB0 + cC0 - cA(time) - cB(time);

time = linspace(0, 5, 100);

figure('Units','centimeters','Position',[5 5 14 10])
hold on

plot(time, cA(time), 'LineWidth', 2.0, 'LineStyle','-','DisplayName','A')
plot(time, cB(time), 'LineWidth', 2.0, 'LineStyle','-','DisplayName','B')
plot(time, cC(time), 'LineWidth', 2.0, 'LineStyle','-','DisplayName','C')

xlim([0 5])
ylim([0 1])

xlabel('$\mathrm{time\ (s)}$','Interpreter','latex','FontSize',16)
ylabel('$\mathrm{concentration\ \left(mol\,m^{-3}\right)}$','Interpreter','latex','FontSize',16)

set(gca, ...
    'FontName','lmodern', ...
    'FontSize',16, ...
    'LineWidth',1.5, ...
    'TickLength',[0.015 0.015], ...
    'Box','on')

legend('Interpreter','latex','Location','east','FontSize',16)

grid off
hold off
```

:::{tip} Check yourself
Set $\mathrm{d}c\un{B}/\mathrm{d}t = 0$ in [](#eq-series-cb) and show that B peaks at
$t\un{max} = \ln(k_2/k_1)/(k_2 - k_1)$. For $k_1 = 1\ \mathrm{s^{-1}}$ and $k_2 = 2\ \mathrm{s^{-1}}$
this gives $t\un{max} = \ln 2 = 0.693\ \mathrm{s}$ and $c\un{B,max} = 0.25\ \mathrm{mol\,m^{-3}}$,
which the script should reproduce. What happens to $t\un{max}$ and $c\un{B,max}$ as $k_2$ grows?
:::

The derivation was tedious, but the result is worth keeping. For different ratios of rate constants
the concentration profiles look very different, and we can leverage the rate-constant ratio to
simplify our kinetic analysis — the topic of the next section.

## Timescales

<!-- source: multiple_reactions.tex L348 -->

We looked at the simple example $\ce{A -> B -> C}$, but in reality reactions are typically more
complicated. Some are very fast while others are slow; some reach equilibrium almost immediately
while others remain far from it for a long time. In such systems the rate constants typically
differ by several orders of magnitude, which leads to scenarios where some species have very high
concentrations while others are barely present.

The disparate timescales make it numerically challenging to solve complex reaction networks:
timescales differing by orders of magnitude lead to **stiff ODE** systems, and the time step
required for stability is set by the fastest mode.

The same separation of timescales, however, can also make our lives easier. It allows us to invoke
approximations that simplify the rate laws or the system of ODEs and, in many cases, give physical
insight into the reaction mechanism. We introduce two such approximations: the
*pseudo-steady-state assumption* and the *quasi-equilibrium assumption*.

### The pseudo-steady-state approximation

<!-- source: multiple_reactions.tex L363 -->

Take the series example from the previous section and write the mass balances,

$$
\begin{aligned}
&\ce{A ->[$k_1$] B ->[$k_2$] C} \\
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} &= -k_1 c\un{A} \\
\frac{\mathrm{d}c\un{B}}{\mathrm{d}t} &= k_1 c\un{A} - k_2 c\un{B} \\
\frac{\mathrm{d}c\un{C}}{\mathrm{d}t} &= k_2 c\un{B} .
\end{aligned}
$$

For concreteness, take the parameter set

$$
k_1 = 1\ \mathrm{s^{-1}} , \qquad
k_2 = 1\ \mathrm{s^{-1}} , \qquad
c\un{A,0} = 1\ \mathrm{mol\,L^{-1}} ,
$$

which gives the concentration profile in [](#fig-cprofile-series).

:::{figure} ../figures/reactions_in_series.png
:label: fig-cprofile-series
:alt: Plot of concentration in mol per cubic metre against time in seconds, from 0 to 10, for the series reaction with k1 equal to k2 equal to 1 per second. A decays exponentially from 1 to zero. B rises to a maximum of about 0.37 near t equals 1 second and then decays back to zero. C rises smoothly from zero and approaches 1.
:width: 60%

Concentration profile for reactions in series.
:::

If we increase $k_2$ while keeping everything else fixed, we approach the limiting case

$$
\ce{A ->[slow] B ->[fast] C} .
$$

The mass balances are unchanged,

$$
\begin{aligned}
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} &= -k_1 c\un{A}
\quad \Rightarrow \quad c\un{A} = c\un{A,0}\exp\left(-k_1 t\right) \\
\frac{\mathrm{d}c\un{B}}{\mathrm{d}t} &= k_1 c\un{A} - k_2 c\un{B} \\
\frac{\mathrm{d}c\un{C}}{\mathrm{d}t} &= k_2 c\un{B} ,
\end{aligned}
$$

and the analytical solution for $c\un{B}$ derived in the previous section, with $c\un{B,0} = 0$, is

$$
c\un{B} = \frac{k_1 c\un{A,0}}{k_2 - k_1}\left[\exp(-k_1 t) - \exp(-k_2 t)\right] .
$$

In the limit $k_2 \gg k_1$ the second exponential decays very quickly and the prefactor becomes
$k_1 c\un{A,0}/k_2$, so

$$
\begin{aligned}
c\un{B} &= \frac{k_1 c\un{A,0}}{k_2}\exp(-k_1 t) \\
c\un{B} &= \frac{k_1}{k_2} c\un{A} ,
\end{aligned}
$$ (eq-pssa-cb)

a much simpler closed form.

:::{figure} ../figures/reactions_in_series_varied_k2.png
:label: fig-pssa-var-k2
:alt: Three panels sharing a time axis from 0 to 10 seconds, for k2 values of 1, 2, 5, 10, 50 and 100 per second. Left, the concentration of A: all six curves coincide, since A does not depend on k2. Middle, the concentration of B: the peak shrinks and shifts to earlier times as k2 increases, from about 0.37 at k2 equals 1 down to almost nothing at k2 equals 100. Right, the concentration of C: all curves rise to 1, faster as k2 increases.
:width: 100%

Concentration profiles with varied $k_2$.
:::

If $k_2$ is very large, B is produced slowly by the first reaction and consumed almost immediately
by the second. This is typical of short-lived, high-energy intermediates. Consequently the
concentration of B stays very small and is approximately constant over time, after a brief
induction period. We arrive at the same conclusion by invoking

$$
r\un{B} = \frac{\mathrm{d}c\un{B}}{\mathrm{d}t} \approx 0 ,
$$ (eq-pssa-def)

the **pseudo-steady-state approximation** (PSSA), also called the quasi-steady-state approximation
(QSSA) or the Bodenstein approximation.

The PSSA does *not* mean that $c\un{B}$ itself is zero. If $c\un{B} = 0$ then no C would ever be
formed, which is clearly wrong. We also do not have a true steady state, since that would require
*all* concentrations to be constant in time. The correct interpretation is that the **net** rate of
production of B is approximately zero, because as quickly as B is produced it is consumed by the
next reaction.

Applying the PSSA to the series example,

$$
\begin{aligned}
\frac{\mathrm{d}c\un{B}}{\mathrm{d}t} &= k_1 c\un{A} - k_2 c\un{B} = 0 \\
0 &= k_1 c\un{A} - k_2 c\un{B} \\
c\un{B} &= \frac{k_1}{k_2} c\un{A} ,
\end{aligned}
$$

which matches the limiting form in [](#eq-pssa-cb). Inserting this into the mass balance for C,

$$
\begin{aligned}
\frac{\mathrm{d}c\un{C}}{\mathrm{d}t} &= k_2 \frac{k_1}{k_2} c\un{A} \\
\frac{\mathrm{d}c\un{C}}{\mathrm{d}t} &= k_1 c\un{A} = c\un{A,0}\, k_1 \exp\left(-k_1 t\right) \\
c\un{C} &= c\un{A,0}\left[1 - \exp\left(-k_1 t\right)\right] .
\end{aligned}
$$ (eq-pssa-cc)

The two derivations give the same result.

:::{figure} ../figures/reactions_in_series_PSSA.png
:label: fig-comparison-ode-pssa
:alt: Four panels comparing the full ODE solution, drawn as solid lines, with the PSSA, drawn as dashed lines, for k2 equal to 1, 5, 10 and 100 per second. At k2 equals 1 the PSSA curve for C rises far too early and clearly disagrees with the ODE solution. By k2 equals 5 the two are close, and at k2 equals 10 and 100 they are indistinguishable on this scale.
:width: 100%

Comparison of the ODE system with the PSSA for various values of $k_2$.
:::

The PSSA agrees with the full ODE solution increasingly well as $k_2$ becomes large, as expected.
The approximation still fails at very short times, however, during the initial fast transient.

:::{admonition} Live example
:class: seealso
Solid lines are the numerical solution of the full ODE system, dashed lines the PSSA result from
[](#eq-pssa-cb) and [](#eq-pssa-cc). Increase $k_2$ to see the two converge, and shorten the time
window to zoom in on the transient at the start, where the PSSA is always off. The MATLAB code that
produces the same result follows underneath.
:::

```{marimo} python
pssa_k2 = mo.ui.slider(
    steps=[1, 2, 5, 10, 20, 50, 100],
    value=10,
    label="Rate constant k2 (1/s)",
    show_value=True,
)
pssa_t_end = mo.ui.slider(
    steps=[0.2, 0.5, 1, 2, 5, 10],
    value=10,
    label="Time window (s)",
    show_value=True,
)
mo.vstack([pssa_k2, pssa_t_end])
```

```{marimo} python
pssa_cA0 = 1.0  # mol/m^3
pssa_k1 = 1.0   # 1/s
pssa_k2v = float(pssa_k2.value)

# Full ODE system, integrated with a fixed-step Runge-Kutta scheme (ode45 in MATLAB)
pssa_t = np.linspace(0, pssa_t_end.value, 4000)
pssa_h = pssa_t[1] - pssa_t[0]


def pssa_model(z):
    """Right-hand side of the material balances for [A, B, C]."""
    return np.array([-pssa_k1 * z[0], pssa_k1 * z[0] - pssa_k2v * z[1], pssa_k2v * z[1]])


pssa_z = np.empty((len(pssa_t), 3))
pssa_z[0] = [pssa_cA0, 0.0, 0.0]
for pssa_i in range(len(pssa_t) - 1):
    pssa_s1 = pssa_model(pssa_z[pssa_i])
    pssa_s2 = pssa_model(pssa_z[pssa_i] + 0.5 * pssa_h * pssa_s1)
    pssa_s3 = pssa_model(pssa_z[pssa_i] + 0.5 * pssa_h * pssa_s2)
    pssa_s4 = pssa_model(pssa_z[pssa_i] + pssa_h * pssa_s3)
    pssa_z[pssa_i + 1] = pssa_z[pssa_i] + pssa_h / 6 * (pssa_s1 + 2 * pssa_s2 + 2 * pssa_s3 + pssa_s4)

# PSSA: cB = k1/k2 cA, cC = cA0 (1 - exp(-k1 t))
pssa_cA = pssa_cA0 * np.exp(-pssa_k1 * pssa_t)
pssa_cB = pssa_k1 / pssa_k2v * pssa_cA
pssa_cC = pssa_cA0 * (1 - np.exp(-pssa_k1 * pssa_t))

pssa_dev = np.abs(pssa_cC - pssa_z[:, 2])
pssa_i_dev = int(np.argmax(pssa_dev))

mo.md(
    f"$k_2/k_1 = {pssa_k2v / pssa_k1:.0f}$. The PSSA predicts "
    f"$c_\\mathrm{{B}} = {pssa_k1 / pssa_k2v:.3f}\\, c_\\mathrm{{A}}$; the largest gap in "
    f"$c_\\mathrm{{C}}$ over this window is {pssa_dev[pssa_i_dev]:.3f} mol m$^{{-3}}$ at "
    f"$t = {pssa_t[pssa_i_dev]:.2f}\\ \\mathrm{{s}}$."
)
```

```{marimo} python
fig_pssa, ax_pssa = plt.subplots(figsize=(5.5, 3.9))

# Zoomed-in windows drop A, which is identical in both models, so B and C fill the axes.
pssa_show_A = pssa_t_end.value >= 1
if pssa_show_A:
    ax_pssa.plot(pssa_t, pssa_z[:, 0], lw=2.0, color="C0", label="A")
    pssa_y_max = 1.0
else:
    pssa_y_max = 1.1 * max(pssa_z[:, 1:].max(), pssa_cC.max())
ax_pssa.plot(pssa_t, pssa_z[:, 1], lw=2.0, color="C1", label="B")
ax_pssa.plot(pssa_t, pssa_z[:, 2], lw=2.0, color="C2", label="C")
ax_pssa.plot(pssa_t, pssa_cB, lw=2.0, ls="--", color="C1", label="B (PSSA)")
ax_pssa.plot(pssa_t, pssa_cC, lw=2.0, ls="--", color="C2", label="C (PSSA)")

ax_pssa.set_xlim(0, pssa_t_end.value)
ax_pssa.set_ylim(0, pssa_y_max)
ax_pssa.set_xlabel("time (s)", fontsize=13)
ax_pssa.set_ylabel(r"concentration (mol m$^{-3}$)", fontsize=13)
ax_pssa.tick_params(labelsize=12, width=1.2, length=5)
for sp_pssa in ax_pssa.spines.values():
    sp_pssa.set_linewidth(1.2)
ax_pssa.legend(loc="best", fontsize=12, frameon=False)

fig_pssa.tight_layout()
fig_pssa
```

In MATLAB, solve the ODE system with `ode45` and overlay the PSSA. Save this as a script file,
since it ends with two local functions, and set `k2` to the value you want:

```matlab
z0 = [1; 0; 0];            % [A; B; C], mol/m3
t  = linspace(0, 10, 4000); % s

A0 = 1;   % mol/m3
k1 = 1;   % 1/s
k2 = 10;  % 1/s

figure('Units','centimeters','Position',[5 5 14 10])
hold on

% Full ODE system
[t_sol, z] = ode45(@(tt,zz) model(tt, zz, k1, k2), t, z0);

% Pseudo-steady-state approximation
[A_pssa, B_pssa, C_pssa] = analytical(A0, t_sol, k1, k2);

plot(t_sol, z(:,1), 'LineWidth', 2.0, 'DisplayName','A')
plot(t_sol, z(:,2), 'LineWidth', 2.0, 'DisplayName','B')
plot(t_sol, z(:,3), 'LineWidth', 2.0, 'DisplayName','C')
plot(t_sol, B_pssa, 'LineWidth', 2.0, 'LineStyle','--', 'DisplayName','B (PSSA)')
plot(t_sol, C_pssa, 'LineWidth', 2.0, 'LineStyle','--', 'DisplayName','C (PSSA)')

xlim([0 10])
ylim([0 1])

xlabel('$\mathrm{time\ (s)}$','Interpreter','latex','FontSize',16)
ylabel('$\mathrm{concentration\ \left(mol\,m^{-3}\right)}$','Interpreter','latex','FontSize',16)

set(gca, ...
    'FontName','lmodern', ...
    'FontSize',16, ...
    'LineWidth',1.5, ...
    'TickLength',[0.015 0.015], ...
    'Box','on')

legend('Interpreter','latex','Location','east','FontSize',16)

grid off
hold off

function dzdt = model(~, z, k1, k2)
A = z(1);
B = z(2);
C = z(3);

dAdt = -k1*A;
dBdt =  k1*A - k2*B;
dCdt =  k2*B;

dzdt = [dAdt; dBdt; dCdt];
end

function [A, B, C] = analytical(A0, t, k1, k2)
A = A0 .* exp(-k1.*t);
B = (k1./k2) .* A;
C = A0 .* (1 - exp(-k1.*t));
end
```

:::{tip} Check yourself
For $k_1 = 1\ \mathrm{s^{-1}}$ and $k_2 = 10\ \mathrm{s^{-1}}$ the PSSA overestimates $c\un{C}$ by
at most $0.077\ \mathrm{mol\,m^{-3}}$, at $t \approx 0.26\ \mathrm{s}$. At
$k_2 = 100\ \mathrm{s^{-1}}$ the gap shrinks to about $0.01\ \mathrm{mol\,m^{-3}}$. Extract these
numbers from your own `ode45` solution. Why does the PSSA always run *ahead* of the true solution
for C?
:::

:::{figure} ../figures/reactions_in_series_initial.png
:label: fig-dev-pssa-ode
:alt: Plot zoomed into the first 0.2 seconds, with concentration from 0 to 0.1. The PSSA curve for B, dashed, jumps immediately to its plateau value of about 0.01, while the ODE curve rises to it over roughly 0.02 seconds. Correspondingly the PSSA curve for C runs above the ODE curve, the gap opening during the transient and then staying roughly constant.
:width: 60%

Deviation between the PSSA and the ODE system at the start of the reaction.
:::

The PSSA reduces the order of the ODE system: an ODE is replaced by an algebraic equation.

:::{admonition} Discussion
:class: seealso
Why does this matter? Modern ODE solvers are fast and accurate, so what is the practical benefit of
the PSSA?
:::

There are two reasons. First, when reaction networks become very large — full mechanisms with
hundreds or thousands of species — reducing the number of ODEs cuts simulation cost dramatically.
Second, when timescales differ by many orders of magnitude the resulting stiff ODE system requires
extremely small time steps for stability. A system whose timescales span $10^{10}$ requires at
least $10^{10}$ steps to integrate the slow process accurately, whereas the PSSA bypasses the fast
mode entirely. We always prefer closed-form analytical rate expressions when we can get them,
because they also make it possible to determine rate equations and reaction orders directly from
experiments.

### The quasi-equilibrium approximation

<!-- source: multiple_reactions.tex L480 -->

Strictly speaking, every elementary reaction in a mechanism is reversible. Treating an elementary
step as irreversible is only justified when the equilibrium lies far toward the products,
$K \gg 1$. We now relax that assumption and consider the fully reversible series

$$
\ce{A} \eqa{k\un{+1}}{k\un{-1}} \ce{B} \eqa{k\un{+2}}{k\un{-2}} \ce{C} .
$$

Consider a scenario with the following values:

- $k\un{+1} = 1\ \mathrm{s^{-1}}$
- $K\un{1} = 2$, so $k\un{-1} = k\un{+1}/K\un{1} = 0.5\ \mathrm{s^{-1}}$
- $k\un{+2} = 1\ \mathrm{s^{-1}}$
- $K\un{2} = 5$, so $k\un{-2} = 0.2\ \mathrm{s^{-1}}$
- $c\un{A,0} = 1\ \mathrm{mol\,L^{-1}}$, $c\un{B,0} = 0.2\ \mathrm{mol\,L^{-1}}$

All rate laws are first order, leading to the concentration profiles in
[](#fig-reversible-series).

:::{figure} ../figures/reactions_in_series_w_eq.png
:label: fig-reversible-series
:alt: Plot of concentration against time from 0 to 5 seconds for the reversible series. A decays from 1 towards about 0.12 without reaching zero. B starts at 0.2, rises to a maximum of about 0.4 near t equals 0.8 seconds, then falls back to about 0.21. C rises from zero and levels off near 0.87. All three approach non-zero equilibrium values.
:width: 60%

Reversible reactions in series.
:::

As with the PSSA we can reduce the full mechanism based on a separation of timescales, but the
underlying assumption is different. Suppose
$k\un{+1}, k\un{-1} \gg k\un{+2}, k\un{-2}$. This means the first reaction equilibrates very
quickly, while the second remains far from equilibrium for a long time.

**Example: water in tanks.**

:::{figure} ../figures/WaterTanks.png
:label: fig-water-tanks
:alt: Hand-drawn analogy showing three tall water tanks side by side, connected in sequence. Between the first and second tank a high opening labelled k1 lets water flow forward and a lower opening labelled k minus 1 lets it flow back. The second and third tanks are connected the same way by openings labelled k2 and k minus 2. The water levels represent the species concentrations.
:width: 60%

Water-tank analogy for reactions in series.
:::

The full system of ODEs is

$$
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} = -k\un{+1}c\un{A} + k\un{-1}c\un{B} = -r_1 ,
$$ (eq-qea-dca)

$$
\frac{\mathrm{d}c\un{B}}{\mathrm{d}t}
= k\un{+1}c\un{A} - k\un{-1}c\un{B} - k\un{+2}c\un{B} + k\un{-2}c\un{C} = r_1 - r_2 ,
$$ (eq-qea-dcb)

$$
\frac{\mathrm{d}c\un{C}}{\mathrm{d}t} = k\un{+2}c\un{B} - k\un{-2}c\un{C} = r_2 .
$$ (eq-qea-dcc)

:::{figure} ../figures/reactions_in_series_w_eq_various_k1.png
:label: fig-var-k1-qea
:alt: Three panels sharing a time axis from 0 to 5 seconds, sweeping k1 over 1, 2, 5, 10, 50 and 100 per second. Left, A: the larger k1, the faster the initial drop, after which all curves relax to the same value near 0.1. Middle, B: for small k1 the curve rises gradually to a broad maximum, while for large k1 it jumps almost vertically to about 0.78 at t equals 0 and then decays; all curves converge to about 0.2. Right, C: all curves rise to about 0.9, faster for larger k1.
:width: 100%

Variation in the rate constant $k\un{+1}$.
:::

The first reaction is rapidly equilibrated from the start. At equilibrium,

$$
r_1 = r\un{+1} - r\un{-1} = 0 ,
$$

so we introduce the constraint

$$
\begin{aligned}
r_1 = 0 &= k\un{+1}c\un{A} - k\un{-1}c\un{B} \\
0 &= K\un{1}c\un{A} - c\un{B}
\qquad \text{with} \qquad K\un{1} = \frac{k\un{+1}}{k\un{-1}} .
\end{aligned}
$$ (eq-qea-constraint)

Adding the mass balances for A and B in [](#eq-qea-dca) and [](#eq-qea-dcb) eliminates $r_1$:

$$
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} + \frac{\mathrm{d}c\un{B}}{\mathrm{d}t}
= -r_1 + r_1 - r_2 = -r_2 .
$$

So far we have made no approximations. Differentiating the equilibrium constraint
[](#eq-qea-constraint) with respect to time gives

$$
0 = K\un{1}\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} - \frac{\mathrm{d}c\un{B}}{\mathrm{d}t} .
$$

We now have two equations in two derivatives. Eliminating $\mathrm{d}c\un{B}/\mathrm{d}t$ gives

$$
\begin{aligned}
\frac{\mathrm{d}c\un{B}}{\mathrm{d}t} &= K\un{1}\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} \\
K\un{1}\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} + \frac{\mathrm{d}c\un{A}}{\mathrm{d}t} &= -r_2 ,
\end{aligned}
$$

so the reduced system is

$$
\begin{aligned}
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} &= \frac{-r_2}{1 + K\un{1}} \\
&= \frac{-1}{1 + K\un{1}}\left(k\un{+2}c\un{B} - k\un{-2}c\un{C}\right) \\
&= \frac{-1}{1 + K\un{1}}\left(k\un{+2}K\un{1}c\un{A} - k\un{-2}c\un{C}\right) \\
c\un{B} &= K\un{1}c\un{A} \\
\frac{\mathrm{d}c\un{C}}{\mathrm{d}t} &= r_2 = k\un{+2}K\un{1}c\un{A} - k\un{-2}c\un{C} .
\end{aligned}
$$ (eq-qea-reduced)

We have again reduced the order of the ODE system by one. The concentration of B has been replaced
by the equilibrium relation $c\un{B} = K\un{1}c\un{A}$. [](#fig-comparison-qea-ode) shows that the
simplified model agrees with the full model whenever the timescale separation holds. The initial
condition for the reduced system must of course be re-evaluated; in the QEA we impose
$c\un{B,0} = K\un{1}c\un{A,0}$ rather than the original $c\un{B,0}$.

:::{figure} ../figures/reactions_in_series_w_eq_QEA.png
:label: fig-comparison-qea-ode
:alt: Three panels comparing the full ODE solution at k1 equal to 1, 10 and 100 per second against the QEA, drawn as black squares. For A, B and C alike the k1 equals 1 curve sits well away from the QEA markers, while the k1 equals 10 and 100 curves lie on top of them over almost the whole time range.
:width: 95%

Comparison of the full ODE system with the QEA for the first reaction.
:::

It is tempting to ask: if $r_1 = 0$, why go to all this trouble? Why not simply plug $r_1 = 0$
directly into the original ODE system to obtain

$$
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} = -r_1 = 0 , \qquad
\frac{\mathrm{d}c\un{B}}{\mathrm{d}t} = r_1 - r_2 = -r_2 , \qquad
\frac{\mathrm{d}c\un{C}}{\mathrm{d}t} = r_2 \, ?
$$

This is clearly wrong. With these equations $c\un{A}$ would never change at all, while $c\un{B}$
would deplete and $c\un{C}$ would grow without any reactant ever being consumed — mass is not
conserved. The resolution is that $r_1$ does *not* actually vanish. Looking at the mass balance for
A,

$$
\begin{aligned}
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} &= -r_1 \\
&= -k\un{+1}c\un{A} + k\un{-1}c\un{B} \\
&= \underbrace{k\un{-1}}_{\rightarrow \infty}
   \underbrace{\left(c\un{B} - K\un{1}c\un{A}\right)}_{\rightarrow 0} \neq 0 .
\end{aligned}
$$

The bracketed term is the equilibrium constraint, which goes to zero, but the prefactor goes to
infinity in the fast-equilibrium limit. The product is finite and nonzero. This is exactly the rate
that drives the slow approach to equilibrium and feeds the second reaction.

The same reduction works in the opposite limit,
$k\un{+1}, k\un{-1} \ll k\un{+2}, k\un{-2}$, where the second reaction equilibrates rapidly.

:::{figure} ../figures/reactions_in_series_w_eq_various_k2.png
:label: fig-qea-var-k2
:alt: Three panels sharing a time axis from 0 to 5 seconds, sweeping k2 over 1, 2, 5, 10, 50 and 100 per second. Left, A: the curves are close together, decaying to about 0.1, slightly faster for larger k2. Middle, B: for k2 equals 1 the curve peaks near 0.4, and as k2 increases the peak flattens away entirely, the largest values starting low and rising monotonically to the common plateau near 0.19. Right, C: all curves rise to about 0.92, faster for larger k2.
:width: 95%

Reversible reactions in series with varied $k\un{+2}$.
:::

:::{figure} ../figures/reactions_in_series_w_eq_QEA_k2.png
:label: fig-comparison-qea-ode-var-k2
:alt: Three panels comparing the full ODE solution at k2 equal to 1, 10 and 100 per second against the QEA, drawn as black squares. As before the k2 equals 1 curve deviates clearly while the k2 equals 10 and 100 curves track the QEA markers closely for A, B and C.
:width: 95%

Comparison of the full ODE system with the QEA for the second reaction.
:::

At small times the QEA and the ODE system still deviate from each other, because equilibrium has
not yet been reached.

:::{figure} ../figures/reactions_in_series_w_eq_initial.png
:label: fig-dev-qea-ode
:alt: Plot of the concentration of B against time, zoomed into the first 0.5 seconds. The k2 equals 10 curve falls steeply from its initial value of 0.2 to a minimum near 0.13 and then turns upward. The k2 equals 100 curve falls much faster and further, to about 0.048, before rising. The dotted QEA line starts lower, at about 0.034, and rises smoothly, converging with the k2 equals 100 curve only after the transient.
:width: 60%

Comparison of the concentration profile of B from the ODE system and the QEA for large
$k\un{+2}$.
:::

## Application of PSSA and QEA to derive global rate expressions

<!-- source: multiple_reactions.tex L635 -->

### PSSA: HBr synthesis

<!-- source: multiple_reactions.tex L638 -->

**Example.** As a first application, consider the photochemical production of hydrogen bromide from
the elements in the gas phase,

$$
\ce{H2 + Br2 -> 2 HBr} .
$$ (eq-hbr-overall)

:::{admonition} Discussion
:class: seealso
Is [](#eq-hbr-overall) an elementary reaction? Why or why not?
:::

The reactions considered up to this point have been simplified cases, mostly stoichiometric or
global reactions, each summarizing an underlying mechanism. We return to mechanisms in detail in
[](#ch-mechanisms). For [](#eq-hbr-overall), if it were elementary the law of mass action would
give

$$
r = k\, c_{\ce{H2}}\, c_{\ce{Br2}} ,
$$

first order in each reactant. Kinetic experiments, however, show that
$r \propto c_{\ce{Br2}}^{1/2}$ at low conversion, a clear sign that the overall reaction does not
happen as written. Instead it proceeds via a radical chain mechanism that starts with the
dissociation of bromine into bromine radicals:

$$
\text{Initiation/termination:}\quad \ce{Br2} \eqa{k\un{+1}}{k\un{-1}} \ce{2 Br\bullet} ,
\qquad r_1 = k\un{+1}c_{\ce{Br2}} - k\un{-1}c_{\ce{Br\bullet}}^2
$$ (eq-hbr-step1)

$$
\text{Propagation:}\quad \ce{Br\bullet + H2} \eqa{k\un{+2}}{k\un{-2}} \ce{HBr + H\bullet} ,
\qquad r_2 = k\un{+2}c_{\ce{Br\bullet}}c_{\ce{H2}} - k\un{-2}c_{\ce{HBr}}c_{\ce{H\bullet}}
$$ (eq-hbr-step2)

$$
\text{Propagation:}\quad \ce{H\bullet + Br2 ->[$k_3$] HBr + Br\bullet} ,
\qquad r_3 = k_3\, c_{\ce{H\bullet}}c_{\ce{Br2}}
$$ (eq-hbr-step3)

These are elementary reactions; together they form a reaction mechanism, and summing the steps
gives back the overall stoichiometric expression [](#eq-hbr-overall).

Suppose we are interested in the production rate of HBr because we have experimental data for it.
$\ce{HBr}$ is produced in steps [](#eq-hbr-step2) and [](#eq-hbr-step3), so

$$
r_{\ce{HBr}} = \frac{\mathrm{d}c_{\ce{HBr}}}{\mathrm{d}t}
= \underbrace{k\un{+2}c_{\ce{Br\bullet}}c_{\ce{H2}}
- k\un{-2}c_{\ce{HBr}}c_{\ce{H\bullet}}}_{r_2}
+ \underbrace{k_3\, c_{\ce{H\bullet}}c_{\ce{Br2}}}_{r_3} .
$$ (eq-rhbr-full)

The radical species $\ce{Br\bullet}$ and $\ce{H\bullet}$ are short-lived intermediates, so we apply the
PSSA,

$$
\frac{\mathrm{d}c_{\ce{Br\bullet}}}{\mathrm{d}t}
= \frac{\mathrm{d}c_{\ce{H\bullet}}}{\mathrm{d}t} = 0 .
$$

Writing the mass balance for each radical species,

$$
\begin{aligned}
\frac{\mathrm{d}c_{\ce{Br\bullet}}}{\mathrm{d}t} &= 0
= \underbrace{2k\un{+1}c_{\ce{Br2}} - 2k\un{-1}c_{\ce{Br\bullet}}^2}_{2r_1}
- \underbrace{\left(k\un{+2}c_{\ce{Br\bullet}}c_{\ce{H2}}
  - k\un{-2}c_{\ce{HBr}}c_{\ce{H\bullet}}\right)}_{r_2}
+ \underbrace{k_3\, c_{\ce{H\bullet}}c_{\ce{Br2}}}_{r_3} \\
\frac{\mathrm{d}c_{\ce{H\bullet}}}{\mathrm{d}t} &= 0
= \underbrace{k\un{+2}c_{\ce{Br\bullet}}c_{\ce{H2}}
  - k\un{-2}c_{\ce{HBr}}c_{\ce{H\bullet}}}_{r_2}
- \underbrace{k_3\, c_{\ce{H\bullet}}c_{\ce{Br2}}}_{r_3} .
\end{aligned}
$$

One way to solve for $c_{\ce{Br\bullet}}$ is simply to add the two equations and work through the
algebra. A cleaner argument follows from the structure of the mechanism: under the PSSA the
propagation steps [](#eq-hbr-step2) and [](#eq-hbr-step3) produce and consume $\ce{Br\bullet}$ in
equal measure, so they cannot be the steps that set $c_{\ce{Br\bullet}}$. Only the
initiation/termination step [](#eq-hbr-step1) alters the radical population, and at steady state
the termination rate must equal the initiation rate. Applying this to $\ce{Br\bullet}$,

$$
\begin{aligned}
2k\un{+1}c_{\ce{Br2}} &= 2k\un{-1}c_{\ce{Br\bullet}}^2 \\
c_{\ce{Br\bullet}} &= \sqrt{\frac{k\un{+1}c_{\ce{Br2}}}{k\un{-1}}} .
\end{aligned}
$$ (eq-hbr-brrad)

Plugging this into the steady-state equation for $\ce{H\bullet}$ gives

$$
c_{\ce{H\bullet}} = \frac{k\un{+2}\sqrt{k\un{+1}c_{\ce{Br2}}/k\un{-1}}\, c_{\ce{H2}}}
{k\un{-2}c_{\ce{HBr}} + k_3\, c_{\ce{Br2}}} .
$$ (eq-hbr-hrad)

Substituting [](#eq-hbr-brrad) and [](#eq-hbr-hrad) into [](#eq-rhbr-full) and simplifying, the
algebra collapses into the closed-form rate law

$$
r_{\ce{HBr}} = \frac{2k\un{+2}\sqrt{k\un{+1}/k\un{-1}}\, c_{\ce{H2}}\sqrt{c_{\ce{Br2}}}}
{1 + \dfrac{k\un{-2}c_{\ce{HBr}}}{k_3\, c_{\ce{Br2}}}} .
$$ (eq-hbr-rate-full)

Experiments show that the rate is one-half order in $\ce{Br2}$ at low conversions. At low
conversion $c_{\ce{HBr}} \to 0$, and [](#eq-hbr-rate-full) simplifies to

$$
r_{\ce{HBr}} = 2k\un{+2}\sqrt{\frac{k\un{+1}}{k\un{-1}}}\, c_{\ce{H2}}\sqrt{c_{\ce{Br2}}} .
$$ (eq-hbr-rate-low)

The individual rate constants cannot be measured independently, so we lump them into an apparent
rate constant

$$
k\un{app} = k\un{+2}\sqrt{\frac{k\un{+1}}{k\un{-1}}} .
$$ (eq-hbr-kapp)

A temperature-dependence measurement of $k\un{app}$ would yield an apparent activation energy

$$
E\un{a,app} = E\un{a,+2} + \frac{E\un{a,+1} - E\un{a,-1}}{2} ,
$$ (eq-hbr-eapp)

so that

$$
r_{\ce{HBr}} = 2k\un{app}\, c_{\ce{H2}}\sqrt{c_{\ce{Br2}}} ,
$$

which is exactly the half-order dependence observed experimentally. This does not prove the assumed
mechanism is correct, only that it is consistent with the experimental rate law.

The full empirical rate expression determined by Bodenstein in 1906 is

$$
r_{\ce{HBr}} = \frac{k_a\, c_{\ce{H2}}\sqrt{c_{\ce{Br2}}}}
{1 + k_b\, c_{\ce{HBr}}/c_{\ce{Br2}}} ,
$$ (eq-hbr-bodenstein)

where $k_a$ and $k_b$ are empirical constants. Comparing with [](#eq-hbr-rate-full) gives
$k_a = 2k\un{+2}\sqrt{k\un{+1}/k\un{-1}}$ and $k_b = k\un{-2}/k_3$.

### Why are no other reactions relevant in the mechanism?

<!-- source: multiple_reactions.tex L754 -->

This is an excellent question, and one we revisit rigorously in [](#ch-mechanisms). Identifying the
relevant reaction pathways within an enormous chemical reaction space is a research-level problem
in its own right. To describe the flame of a candle, for instance, roughly $10^4$ reactions are
typically needed, which means starting from a candidate space of perhaps $10^6$ possible reactions.

For the HBr system specifically, a useful rule-of-thumb estimate explains why we restrict the
mechanism to $\ce{Br2}$ dissociation rather than including $\ce{H2}$ or $\ce{HBr}$ dissociation,

$$
\ce{H2 -> 2 H\bullet} \qquad \ce{HBr -> H\bullet + Br\bullet} .
$$

The bond-dissociation enthalpies — the reaction enthalpies of these elementary steps — are

- $\ce{H2 -> 2 H\bullet}$: $\Delta\un{r}H = 436\ \mathrm{kJ\,mol^{-1}}$
- $\ce{HBr -> H\bullet + Br\bullet}$: $\Delta\un{r}H = 399\ \mathrm{kJ\,mol^{-1}}$
- $\ce{Br2 -> 2 Br\bullet}$: $\Delta\un{r}H = 194\ \mathrm{kJ\,mol^{-1}}$

:::{figure} ../figures/EaBDE.png
:label: fig-bde
:alt: Hand-drawn energy diagram against a reaction coordinate showing two endothermic pathways from a common reactant level up to two different product levels. Each pathway has a barrier drawn just above its own product level, annotated as the barrierless limit, with a vertical arrow marking the reaction enthalpy. The pathway ending at the higher product level necessarily has the higher barrier, annotated higher activation energy and therefore slower.
:width: 75%

Bond-dissociation enthalpy as a lower bound on the activation barrier of an endothermic
dissociation step.
:::

For an endothermic reaction the reaction enthalpy is a lower bound on the activation barrier, so the
$\ce{Br2}$ dissociation has by far the lowest barrier and is much faster than $\ce{H2}$ or
$\ce{HBr}$ dissociation. When the thermodynamics are not available — some intermediates are not in
any database — electronic-structure calculations are required to estimate the barriers and rate
constants directly.

### QEA: $\ce{NO + Br2}$ synthesis of nitrosyl bromide

<!-- source: multiple_reactions.tex L783 -->

**Example.** Consider the reaction

$$
\ce{2 NO + Br2 <=> 2 NOBr} ,
$$

where $\ce{NOBr}$ is nitrosyl bromide. The reaction proceeds via the mechanism

$$
\ce{NO + Br2} \eqa{k\un{+1}}{k\un{-1}} \ce{NOBr2} \qquad \text{very fast}
$$ (eq-nobr-step1)

$$
\ce{NOBr2 + NO ->[$k_2$] 2 NOBr} \qquad \text{slow}
$$ (eq-nobr-step2)

The production rate of $\ce{NOBr}$ is

$$
\frac{\mathrm{d}c_{\ce{NOBr}}}{\mathrm{d}t} = r_{\ce{NOBr}}
= 2k_2\, c_{\ce{NOBr2}}\, c_{\ce{NO}} .
$$ (eq-nobr-rate-raw)

The intermediate concentration $c_{\ce{NOBr2}}$ is hard to measure, but the quasi-equilibrium
approximation for the fast step [](#eq-nobr-step1) relates it directly to the reactant
concentrations:

$$
\begin{aligned}
r_1 &= r\un{+1} - r\un{-1} = 0 \\
k\un{+1}\, c_{\ce{NO}}\, c_{\ce{Br2}} &= k\un{-1}\, c_{\ce{NOBr2}} \\
c_{\ce{NOBr2}} &= K\un{1}\, c_{\ce{NO}}\, c_{\ce{Br2}} .
\end{aligned}
$$ (eq-nobr-intermediate)

Substituting [](#eq-nobr-intermediate) into [](#eq-nobr-rate-raw) gives the closed-form rate

$$
r_{\ce{NOBr}} = 2k_2 K\un{1}\, c_{\ce{NO}}^2\, c_{\ce{Br2}} .
$$ (eq-nobr-rate)

:::{admonition} Discussion
:class: seealso
What is the rate of reaction $r$ associated with this rate law?
:::

Using $r = r_{\ce{NOBr}}/\nu_{\ce{NOBr}}$ with $\nu_{\ce{NOBr}} = 2$,

$$
r = \frac{r_{\ce{NOBr}}}{2} = k_2 K\un{1}\, c_{\ce{NO}}^2\, c_{\ce{Br2}} .
$$

### The rate-determining step assumption

<!-- source: multiple_reactions.tex L827 -->

Consider a more complex reaction sequence in which one step is much slower than the others. That
slowest step is the bottleneck, or **rate-determining step** (RDS), and the overall rate is set
entirely by its rate. This is the limiting form of the QEA: all the other steps are assumed to be
in quasi-equilibrium relative to the slow step.

Consider the reaction

$$
\ce{NO2 + CO -> NO + CO2} .
$$

Experiments show that the rate is second order in $\ce{NO2}$ and independent of $\ce{CO}$. The
mechanism is

$$
\ce{2 NO2 ->[$k_1$] NO + NO3\bullet} \qquad \text{very slow}
$$ (eq-no2-step1)

$$
\ce{NO3\bullet + CO} \eqa{k\un{2+}}{k\un{2-}} \ce{NO2 + CO2} \qquad \text{fast}
$$ (eq-no2-step2)

where $\ce{NO3\bullet}$ is the nitrate radical. If the first step is rate-determining, the overall rate
is approximated directly from [](#eq-no2-step1),

$$
r = r_1 = k_1\, c_{\ce{NO2}}^2 ,
$$ (eq-no2-rate)

in agreement with the experimentally observed rate law. The remaining reactions are then assumed to
be in quasi-equilibrium.

The RDS assumption must be applied with care, particularly when the slow step is itself reversible.
If the reverse rate of the slow step is comparable to its forward rate near equilibrium, the simple
approximation $r \approx r\un{forward}$ breaks down. The net rate
$r\un{forward} - r\un{reverse}$ can be very different from $r\un{forward}$ alone, even when the
forward step is the slowest in the mechanism. In such cases the full quasi-equilibrium analysis is
required.

## Summary

<!-- source: multiple_reactions.tex L861 -->

- Multi-reaction systems fall into four broad classes: parallel ($\ce{A -> B}$, $\ce{A -> C}$),
  series ($\ce{A -> B -> C}$), independent (no shared species), and complex (combinations of the
  above).
- For two parallel first-order reactions in a constant-volume batch reactor,
  $c\un{A} = c\un{A,0}\exp[-(k_1+k_2)t]$, [](#eq-parallel-ca), and the product ratio is the
  rate-constant ratio $k_1/k_2$, [](#eq-parallel-ratio).
- Conversion, [](#eq-conversion-def), selectivity, [](#eq-selectivity-def), and yield,
  [](#eq-yield-def), are linked by $Y\un{B} = X\un{A}\cdot S\un{B,A}$, [](#eq-yield-xs). Integral
  and point versions answer different questions: integral for cumulative batch yields, point for
  the local rate behavior used in reactor design.
- For a series first-order reaction $\ce{A -> B -> C}$, the integrating-factor method gives
  $c\un{B}$, [](#eq-series-cb), and the species balance gives $c\un{C}$, [](#eq-series-cc). The
  intermediate B goes through a maximum, so the optimal stop-time is finite.
- PSSA: when an intermediate is short-lived ($k\un{out} \gg k\un{in}$), set
  $\mathrm{d}c\un{B}/\mathrm{d}t \approx 0$, [](#eq-pssa-def). It is the *net* rate of production
  that is approximately zero, not $c\un{B}$ itself.
- QEA: when one step equilibrates much faster than the others, impose the equilibrium constraint
  $c\un{B} = K\un{1}c\un{A}$, [](#eq-qea-constraint); the rate of the fast step does *not* vanish —
  it is the indeterminate form $k\un{-1} \to \infty$ times
  $(K\un{1}c\un{A} - c\un{B}) \to 0$.
- Both PSSA and QEA reduce the order of the ODE system by replacing one ODE with an algebraic
  constraint — valuable both for stiff numerical problems and for deriving closed-form rate laws to
  compare with experiment.
- Applying the PSSA to the radical chain mechanism of HBr synthesis recovers the experimentally
  observed half-order dependence on $\ce{Br2}$, [](#eq-hbr-rate-low), and predicts an apparent
  activation energy $E\un{a,app} = E\un{a,+2} + (E\un{a,+1} - E\un{a,-1})/2$, [](#eq-hbr-eapp).
- Applying the QEA to a fast pre-equilibrium step (the $\ce{NO + Br2}$ mechanism) gives the
  closed-form rate $r_{\ce{NOBr}} = 2k_2 K\un{1} c_{\ce{NO}}^2 c_{\ce{Br2}}$, [](#eq-nobr-rate).
- The rate-determining-step assumption is the limiting form of the QEA: the slowest step sets the
  overall rate, with all other steps in quasi-equilibrium. Use with caution when the slow step is
  reversible.
