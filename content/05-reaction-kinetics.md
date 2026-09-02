---
title: Reaction Kinetics
short_title: Reaction kinetics
label: ch-reaction-kinetics
---

<!-- LaTeX source: ReactionKinetics.tex -->
<!-- Porting notes: mhchem is math-only here, so \ce{} in prose needs $...$; and \un{} expands to
     _{\textrm{...}}, so species subscripts must be written _{\ce{CH4}}, never \un{\ce{CH4}}.
     Nested directives need a longer outer fence (:::: around :::). -->

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

<!-- source: ReactionKinetics.tex L18-84 -->

For the start of our journey into reaction kinetics we will focus exclusively on **homogeneous
reactions**, where the reaction occurs in just one phase — in the liquid or the gas phase, for
instance. Consider such a homogeneous system in an isothermal batch reactor at constant pressure,
in which a single reaction occurs, say the isomerization $\ce{A -> B}$.

The progress of this reaction is described by the extent of reaction $\xi$, defined in
[](#ch-stoichiometry) as

$$
\mathrm{d}\xi = \frac{\mathrm{d}n_i}{\nu_i} .
$$

Since this quantity describes the progress of the reaction towards equilibrium, the rate at which
the reaction proceeds is simply the change of the extent of reaction with time,

$$
\text{rate of reaction} = \frac{\mathrm{d}\xi}{\mathrm{d}t} .
$$

The units would be mol s$^{-1}$, making this an extensive property that depends on system size.
Almost always we prefer an intensive quantity on a volume basis, so the **rate of reaction** $r$ is
defined as

$$
r = \frac{1}{V}\frac{\mathrm{d}\xi}{\mathrm{d}t} ,
\qquad r\ [=]\ \mathrm{mol\,m^{-3}\,s^{-1}} .
$$

Note that $r$ is always positive, that the reaction rate is independent of any particular species
participating in the reaction, and that the time derivative is essential, since the rate changes as
a function of time.

:::{admonition} Discussion
:class: seealso
How do we determine the change in the extent of reaction with time?
:::

Replacing $\mathrm{d}\xi$ with its definition gives the reaction rate as

$$
r = \frac{1}{\nu_i}\frac{1}{V}\frac{\mathrm{d}n_i}{\mathrm{d}t} .
$$

The production or destruction rate of a specific species, $r_i$, is then

$$
r_i = \frac{1}{V}\frac{\mathrm{d}n_i}{\mathrm{d}t} ,
$$

and it is related to the reaction rate through the stoichiometric coefficient:

$$
r_i = \nu_i r .
$$ (eq-ri-nui-r)

**Example: steam reforming of methane**, a large industrial process for producing $\ce{H2}$:

$$
\begin{aligned}
&\ce{CH4 + H2O -> CO + 3 H2} , \\
&r = \frac{r_i}{\nu_i} , \\
&\frac{r_{\ce{CH4}}}{-1} = \frac{r_{\ce{H2O}}}{-1}
 = \frac{r_{\ce{CO}}}{1} = \frac{r_{\ce{H2}}}{3} = r .
\end{aligned}
$$

:::{admonition} In-class exercise
:class: note
If $r_{\ce{H2}} = 18\ \mathrm{mol\,m^{-3}\,s^{-1}}$, what is $r_{\ce{CH4}}$?
:::

We can also formulate the reaction rate in terms of molar concentrations. Assuming $V$ is constant,
this recovers the batch-reactor material balance introduced in [](#ch-mass-balances):

$$
r_i = \frac{\mathrm{d}c_i}{\mathrm{d}t} = \nu_i r .
$$

We still need to know what $r$ *is*, and to express it in a form that shows how it depends on
composition. We could also normalize by the surface area or the mass of a catalyst, but we will
return to that when we discuss heterogeneous catalysis.

## Boudart's rules for reaction kinetics

<!-- source: ReactionKinetics.tex L85 -->

So far nothing has been said about what $r$ actually is — we have only shuffled things around. To
address the reaction rate, we introduce five rules, sometimes referred to as **Boudart's rules**.
Michel Boudart was a famous chemical engineering professor at Stanford University. He neither
invented nor discovered these rules, but he formalized them in his textbooks.

Consider a general chemical reaction,

$$
\begin{aligned}
&\ce{\nu_1 A_1 + \nu_2 A_2 <=> \nu_3 A_3 + \nu_4 A_4} , \\
\text{or mathematically}\quad & 0 = \sum_{i=1}^{N\un{species}} \nu_i A_i .
\end{aligned}
$$

To begin, observe that the reaction rate $r$ is a very complex and unknown function depending on
temperature $T$, composition $c_i$, pressure $p$, and many other things:

$$
r = f\left(c_i, T, p\right) .
$$

Typically $r$ has to be determined experimentally. Note that these are *rules*, not laws like those
of thermodynamics — which is to say it is acceptable to break one.

1. The rate of reaction decreases monotonically with increasing extent of reaction.

2. The rate of reaction in a single direction can be written as

   $$
   r = \underbrace{k(T)}_{\text{rate constant}} f(c_i) ,
   $$

   the **rate constant assumption**: the dependence is split into a rate constant that depends on
   temperature and an unknown function that depends on composition.

3. The rate constant has an **Arrhenius** form,

   $$
   k(T) = A \exp\left(\frac{-E\un{a}}{RT}\right) ,
   $$ (eq-arrhenius)

   where $A$ is the **pre-exponential coefficient** and $E\un{a}$ the **activation energy**. Named
   after Svante Arrhenius, who received the Nobel Prize in 1903. It was known that increasing the
   temperature increases the rate of reaction; he observed that the dependence is exponential.

4. The functional form of the concentration dependence can be approximated by a **power law**,
   which works quite frequently:

   $$
   f(c_i) = \prod_i c_i^{\alpha_i}
   \qquad \text{or} \qquad
   f(p_i) = \prod_i p_i^{\alpha_i} .
   $$

5. The net rate of reaction is the difference between the forward and the reverse reaction,

   $$
   r = r\un{fwd} - r\un{rev} .
   $$

Let us look at these rules in more detail.

### 1. The rate decreases monotonically with extent of reaction

When we start a reaction we are typically far from equilibrium, so the driving force is large and
the reaction is fast. As reactants are consumed and we approach equilibrium, the reverse reaction
kicks in and rates slow down. The result is an exponential decay in the rate.

:::{figure} ../figures/Monotonic.png
:label: fig-monotonic
:alt: Two hand-drawn plots. Left, reaction rate on the vertical axis against extent of reaction on the horizontal axis: the curve starts high at the initial composition and decays exponentially, flattening onto a horizontal line marked equilibrium. A note below reads that as the reaction progresses from initial to final state it slows down, because reactants become depleted and the reverse reaction contributes. Right, labelled autocatalytic, rate against extent of reaction: the curve starts near zero and accelerates upward.
:width: 95%

Left: normal kinetics, with the rate decreasing monotonically as the extent of reaction increases.
Right: an autocatalytic reaction $\ce{A + B -> B + B}$, where the rate *accelerates* with extent of
reaction.
:::

There is one exception to this rule: **autocatalytic reactions**. A reaction such as
$\ce{A + B -> B + B}$ accelerates with extent of reaction.

:::{admonition} Discussion
:class: seealso
Can anyone name an example of an autocatalytic reaction?
:::

### 2. The rate constant assumption

The reaction rate equation is a complex function, so we try to factor out some — if not all — of
the temperature dependence into a prefactor, our rate constant. We will say much more about rate
constants when we look into transition state theory in [](#ch-microscopic).

### 3. The rate constant has an Arrhenius form

Most of you have seen the Arrhenius expression, [](#eq-arrhenius), before. It can be linearized as

$$
\ln k = \ln A - \frac{E\un{a}}{R}\frac{1}{T} .
$$ (eq-arrhenius-linearized)

:::{figure} ../figures/Arrhenius.png
:label: fig-arrhenius
:alt: Two hand-drawn plots. Left, ln k on the vertical axis against 1 over T on the horizontal axis: a straight line of negative slope, annotated slope equals minus E_a over R, with the low end of the horizontal axis marked hotter and the high end marked colder. Right, k on the vertical axis against ln p on the horizontal axis: a curve that rises steeply at low pressure and then flattens to a plateau.
:width: 95%

Left: the linearized Arrhenius form, $\ln k$ vs. $1/T$, with slope $-E\un{a}/R$ and intercept
$\ln A$. Right: the weak pressure dependence of the rate constant, $k$ vs. $\ln p$.
:::

Activation energies are typically between 40 and 200 kJ mol$^{-1}$. You can also use the Arrhenius
equation to find the ratio between two rate constants $k_1$ and $k_2$ at two temperatures $T_1$ and
$T_2$:

$$
\frac{k_2}{k_1} = \exp\left[-\frac{E\un{a}}{R}\left(\frac{1}{T_2} - \frac{1}{T_1}\right)\right] .
$$ (eq-arrhenius-ratio)

A rule of thumb is that a change of 10 K doubles the reaction rate, for
$E\un{a} \approx 60\ \mathrm{kJ\,mol^{-1}}$. Assuming the rate constant depends only on temperature
is not quite correct — it has a weak pressure dependence as well, shown in the right panel of
[](#fig-arrhenius) — but we will not cover that here.

:::{admonition} Discussion
:class: seealso
If $k$ is a function of temperature, why do we call it a *constant*? What is it constant with
respect to?
:::

A second question: why do we call it an activation *energy* if $E\un{a}/R$ has units of
temperature?

:::{figure} ../figures/ActivationEnergy.png
:label: fig-activation-energy
:alt: Hand-drawn energy diagram against a reaction coordinate. A flat line on the left labelled energy of reactants rises to a single peak and falls to a lower flat line on the right labelled energy of products. A double-headed vertical arrow spans from the reactant level to the peak, labelled activation energy.
:width: 75%

$E\un{a}/R$ has units of temperature.
:::

It is related to some barrier of activation. We leave this vague intentionally, because we will go
into it when we discuss transition state theory in [](#ch-microscopic). Spoiler: it has to do with
free energy.

### 4. The power law approach

We assume the unknown concentration dependence of the rate equation takes the power-law form

$$
f(c_i) = \prod_i c_i^{\alpha_i} ,
$$

where the $\alpha_i$ can be any real numbers — positive, negative, or fractional. For an
irreversible reaction $\nu_1 A_1 + \nu_2 A_2 \to \dots$ we can write

$$
r = k\, c_1^{\alpha_1} c_2^{\alpha_2} ,
$$

where $\alpha_i$ is the **reaction order** with respect to species $i$ and
$\alpha = \sum_i \alpha_i$ is the **total reaction order**. The $\alpha_i$ are *not* the
stoichiometric coefficients; they are empirical properties determined from experiments.

Consider the units of $k$. The rate $r$ is in mol m$^{-3}$ s$^{-1}$ almost always, and the
concentrations are in mol m$^{-3}$ raised to the total order $\alpha$. Thus

$$
k\ [=]\ \left(\frac{\text{volume}}{\text{mole}}\right)^{\alpha-1}\frac{1}{\text{time}} .
$$

**Examples:**

1. Decomposition of acetaldehyde, $\ce{CH3CHO -> CH4 + CO}$, with
   $r = k\, c_{\ce{CH3CHO}}^{1.5}$.
2. Ammonia synthesis, $\ce{N2 + 3 H2 -> 2 NH3}$, with
   $r = k\, c_{\ce{N2}}\, c_{\ce{H2}}^{2.25}\, c_{\ce{NH3}}^{-1.5}$.
3. Decomposition of diethyl ether, $\ce{(C2H5)2O -> C2H6 + CH3CHO}$, with
   $r = k\, c_{\ce{(C2H5)2O}}$.

Sometimes species occur in a mixture but do not appear in the rate expression; these are **inert**
species, though they still affect the concentrations of the mixture. Sometimes species appear in
the power law without being part of the reaction — a catalyst or an inhibitor, for example. If
$\alpha_i > 0$ for a product we call the reaction **autocatalytic**, since the product accelerates
the rate as it is formed. If $\alpha_i < 0$ for a product, that product is an **inhibitor** and
slows the reaction down as it progresses. In example 2 above, $\ce{NH3}$ inhibits the reaction.

There is one exception, where the reaction orders *are* the stoichiometric numbers: **elementary
reactions**. Elementary reactions are the individual reaction steps that make up a reaction
mechanism, and we will study them in detail in [](#ch-mechanisms). We will not dive into the fifth
Boudart rule yet, leaving it open until we reach reversible reactions in
[](#sec-reversible-first-order).

### The law of mass action

The stoichiometric coefficients used so far are *net* stoichiometric coefficients. Consider a
general reaction,

$$
\ce{A + B <=> C + D} ,
$$

with net stoichiometric coefficients

$$
\nu\un{A} = -1 , \qquad \nu\un{B} = -1 , \qquad \nu\un{C} = 1 , \qquad \nu\un{D} = 1 .
$$

We can define separate stoichiometric numbers for the forward and reverse directions,

$$
\nu_i' = \text{forward} , \qquad\qquad \nu_i'' = \text{reverse} ,
$$

which are essentially the stoichiometric numbers. The net coefficient is the difference between the
two,

$$
\nu\un{A} = \nu\un{A}'' - \nu\un{A}' .
$$

The **law of mass action** applies when the reaction orders are the stoichiometric coefficients,

$$
\alpha_i^{\text{fwd}} \to \nu_i' , \qquad\qquad \alpha_i^{\text{rev}} \to \nu_i'' ,
$$

so that the rate expression for the equilibrium reaction becomes

$$
r = r\un{fwd} - r\un{rev}
  = k\un{fwd} \prod_i c_i^{\nu_i'} - k\un{rev} \prod_i c_i^{\nu_i''} .
$$ (eq-mass-action)

This approach is called **mass action kinetics**. As noted, it typically does not work for global
or lumped reactions, but it does work for elementary reactions.

<!-- 
## Integrated rate laws

source: ReactionKinetics.tex L287
Reaction rate laws can be rather complex functions and frequently have to be determined from
experiments. For some simple kinetics, however, we can find an analytical solution to the material
balance of the batch reactor. Such a solution gives the temporal concentration profile, and we can
use it to design batch reactors. We will look at a few rate laws, starting with the simplest case.

### Irreversible, first-order reaction

source: ReactionKinetics.tex L297
We assume an irreversible reaction in a batch reactor,

$$
\ce{A -> B} ,
$$

an isomerization, for example, and take it to be described by a first-order power law,

$$
r = k\, c\un{A} .
$$

Implementing this in the batch-reactor mole balance for species A gives

$$
\begin{aligned}
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} &= \nu\un{A}\, r \\
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} &= -k\, c\un{A} .
\end{aligned}
$$

This is an ordinary differential equation, and it can be solved by separating the variables and
integrating. The initial condition is a starting concentration $c\un{A,0}$ at $t = 0$:

$$
\begin{aligned}
\int_{c\un{A,0}}^{c\un{A}} \frac{1}{c\un{A}}\, \mathrm{d}c\un{A} &= -k \int_0^t \mathrm{d}t \\
\ln\left(\frac{c\un{A}}{c\un{A,0}}\right) &= -kt
\quad \Rightarrow \quad
c\un{A} = c\un{A,0} \exp\left(-kt\right) .
\end{aligned}
$$ (eq-first-order-solution)

The units of $k$ for a first-order reaction are s$^{-1}$. Similarly, we can solve for $c\un{B}$
from its material balance,

$$
\frac{\mathrm{d}c\un{B}}{\mathrm{d}t} = \nu\un{B}\, r = k\, c\un{A} .
$$

To solve this we need $c\un{B}$ in terms of $c\un{A}$, which the overall mass balance supplies:

$$
c\un{A} + c\un{B} = c\un{A,0} + c\un{B,0}
\quad \Rightarrow \quad
c\un{A} = c\un{A,0} + c\un{B,0} - c\un{B} .
$$

Substituting and solving the integral leads to

$$
c\un{B} = c\un{A,0}\left(1 - \exp\left(-kt\right)\right) .
$$

Alternatively we could have used the mass balance directly, since we already have an expression for
the concentration of A:

$$
c\un{B} = c\un{A,0} + c\un{B,0} - c\un{A}
\quad \Rightarrow \quad
c\un{B} = c\un{A,0} + c\un{B,0} - c\un{A,0}\exp\left(-kt\right) .
$$

So always be smart and think first, before launching into a more complex integration.

With the rate constant set by the Arrhenius equation, $k(T) = A\,\mathrm{e}^{-E\un{a}/RT}$, the
concentration profiles depend on temperature only through $k$ — but that dependence is exponential,
and it is worth developing intuition for how violent it is. Drag the temperature slider below and
watch the profiles: over a span of 80 K the reaction goes from barely proceeding to essentially
complete within the same 200 s window.

:::{admonition} Live example
:class: seealso
The figure below is interactive and runs entirely in your browser — nothing to install, though it
takes a few seconds to wake up the first time. The MATLAB code that produces the same result
follows underneath.
:::

```{marimo} python
import marimo as mo
import numpy as np
import matplotlib.pyplot as plt

cA0 = 1.0    # mol/m^3, initial concentration of A
cB0 = 0.0    # mol/m^3, initial concentration of B
Ea = 60.0    # kJ/mol, activation energy
A_pre = 1e6  # 1/s, pre-exponential factor
R = 8.314    # J/(mol K), gas constant


def rate_constant(temp):
    """Arrhenius rate constant k(T) in 1/s, with temp in K."""
    return A_pre * np.exp(-Ea * 1e3 / (R * temp))
```

```{marimo} python
T = mo.ui.slider(
    start=340,
    stop=420,
    step=1,
    value=380,
    label="Temperature / K",
    show_value=True,
)
T
```

```{marimo} python
k = rate_constant(T.value)
time = np.linspace(0, 200, 200)
cA = cA0 * np.exp(-k * time)
cB = cB0 + cA0 * (1 - np.exp(-k * time))

mo.md(
    f"At $T = {T.value}\\ \\mathrm{{K}}$: "
    f"$k = {k:.3e}\\ \\mathrm{{s^{{-1}}}}$, "
    f"half-life $t_{{1/2}} = {np.log(2) / k:.1f}\\ \\mathrm{{s}}$."
)
```

```{marimo} python
fig, ax = plt.subplots(figsize=(5.5, 3.9))

ax.plot(time, cA, lw=2.0, label="A")
ax.plot(time, cB, lw=2.0, label="B")

ax.set_xlim(0, 200)
ax.set_ylim(0, 1)
ax.set_xlabel("time (s)", fontsize=13)
ax.set_ylabel(r"concentration (mol m$^{-3}$)", fontsize=13)
ax.tick_params(labelsize=12, width=1.2, length=5)
for spine in ax.spines.values():
    spine.set_linewidth(1.2)
ax.legend(loc="center right", fontsize=13, frameon=False)

fig.tight_layout()
fig
```

The slider is a convenience for exploring the temperature dependence in the browser; the
calculation behind it is the one you will write yourself. In MATLAB, set `T` to the temperature you
want and run:

```matlab
cA0 = 1;      % mol/m^3
cB0 = 0;      % mol/m^3
Ea  = 60;     % kJ/mol
R   = 8.314;  % J/mol/K

T = 380;      % K -- change this to move along the slider

k = @(temp) 1e6*exp(-Ea*1e3/R/temp);   % 1/s

cA = @(time) cA0*exp(-k(T)*time);
cB = @(time) cB0 + cA0*(1 - exp(-k(T)*time));

time = linspace(0, 200, 100);

figure('Units','centimeters','Position',[5 5 14 10])
hold on

plot(time, cA(time), 'LineWidth', 2.0, 'LineStyle','-','DisplayName','A')
plot(time, cB(time), 'LineWidth', 2.0, 'LineStyle','-','DisplayName','B')

xlim([0 200])
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
At $T = 380$ K the script should give $k = 5.65\times10^{-3}\ \mathrm{s^{-1}}$ and a half-life of
about 123 s. Reproduce those numbers before trusting anything else the script tells you.
:::

A useful quantity is the **half-life** of the reactant, the time at which half of the initial
reactant concentration has been converted to products:

$$
\begin{aligned}
\frac{c\un{A}}{c\un{A,0}} &= \exp\left(-kt\right) \\
\frac{1}{2} &= \exp\left(-k\, t_{1/2}\right) \\
t_{1/2} &= \frac{\ln 2}{k} .
\end{aligned}
$$ (eq-half-life-first-order)

The main takeaway is that **the half-life does not depend on the initial concentration** for a
first-order reaction. Starting from any concentration of A, it always takes $\ln(2)/k$ to convert
half of it. You will work with this quantity in the homework assignment.

### Irreversible, zero-order reaction

source: ReactionKinetics.tex L366
There are of course many other reaction orders, and we give two more examples. The first is a
zero-order rate law.

:::{admonition} Discussion
:class: seealso
What is the rate expression $r$ if the reaction is zero-order?
:::

$$
r = k , \qquad k\ [=]\ \mathrm{mol\,m^{-3}\,s^{-1}} .
$$ (eq-zero-order-rate)

For a zero-order reaction the rate does not depend on the concentration of the reactant. Following
the same approach as before,

$$
\begin{aligned}
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} &= \nu\un{A} r \\
                                      &= -k ,
\end{aligned}
$$

and integration leads to

$$
c\un{A} = c\un{A,0} - kt .
$$ (eq-zero-order-solution)

An example is the catalytic decomposition of dinitrogen monoxide ($\ce{N2O}$, laughing gas) over a
hot platinum wire,

$$
\ce{2 N2O ->[Pt] 2 N2 + O2} .
$$

Strictly speaking this is a heterogeneous reaction — it occurs on the Pt surface, not in the bulk
gas phase. We use it here as a clean illustration of a zero-order rate, since the rate depends only
on the surface area of the hot Pt wire and not on the gas-phase concentration.

An interesting feature of a zero-order reaction is that we can determine exactly when all of A has
been converted:

$$
t\un{end} = \frac{c\un{A,0}}{k} .
$$ (eq-zero-order-tend)

### Irreversible, second-order reaction

source: ReactionKinetics.tex L397
The last example we discuss is second-order reactions, of which there are two types. Consider

$$
\ce{A -> products} ,
\qquad\text{for example}\qquad
\ce{2 NO2 ->[$\Delta$] 2 NO + O2} ,
$$

which could have the power law

$$
r = k\, c\un{A}^2 \qquad \text{(Case I)} ,
$$

and would be second order with respect to the concentration of A. Keep an eye on the stoichiometric
coefficient when you use a real example: for the $\ce{NO2}$ decomposition above
$\nu_{\ce{NO2}} = -2$, so the mole balance picks up an extra factor of two relative to the generic
$\nu\un{A} = -1$ case worked out below.

Alternatively, a second-order rate law can arise for

$$
\ce{A + B -> C} .
$$

Most reactions are not simply A converting to B; more typically two reactants combine to form the
product, which frequently follows

$$
r = k\, c\un{A} c\un{B} \qquad \text{(Case II)} .
$$

This reaction is first order with respect to A and first order with respect to B, but the total
order is two.

For Case I the procedure is exactly as before, and fairly straightforward:

$$
\begin{aligned}
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} &= \nu\un{A} r \\
                                      &= -k\, c\un{A}^2 ,
\qquad\text{with}\qquad k\ [=]\ \mathrm{m^3\,mol^{-1}\,s^{-1}} .
\end{aligned}
$$

Integration is a little more difficult, but still easy enough:

$$
\begin{aligned}
\int_{c\un{A,0}}^{c\un{A}} \frac{1}{c\un{A}^2}\, \mathrm{d}c\un{A} &= -k\int_0^t \mathrm{d}t \\
c\un{A} &= \frac{c\un{A,0}}{1 + k t\, c\un{A,0}} .
\end{aligned}
$$ (eq-second-order-casei)

The solution generalizes to an $n$-th order reaction, $r = k c\un{A}^n$; you can find the
derivation and result in the Hill & Root textbook. Before moving to Case II, let us take some time
to plug numbers into the material balances we have derived.

### Dimensionless equations

source: ReactionKinetics.tex L440
You can already see that there are many rate laws out there, each with very different rate
constants. That makes it difficult to compare different reaction orders, starting concentrations,
and rate constants, since they carry different units. This is why the chemical engineering
community developed dimensionless numbers.

:::{admonition} Discussion
:class: seealso
Have some of you already heard about dimensionless numbers? Can you name one?
:::

To define them, look at the batch-reactor material balance for the first-order rate law,

$$
c\un{A} = c\un{A,0} \exp\left(-kt\right) .
$$

Dividing by $c\un{A,0}$ makes everything dimensionless:

$$
f = \frac{c\un{A}}{c\un{A,0}} = \exp\left(-kt\right) .
$$

The quantity $f$ is the fraction of reactant not yet converted,

$$
f = 1 - X = \frac{n_i}{n\un{i,0}} ,
\qquad\text{and for } V = \text{const.} \qquad
f = \frac{c_i}{c\un{i,0}} ,
$$

related to the conversion $X$ through

$$
X = \frac{n\un{i,0} - n_i}{n\un{i,0}} .
$$

The term in the exponential is already dimensionless, since $k\ [=]\ \mathrm{s^{-1}}$. This
dimensionless quantity is the **Damköhler number** $Da_I$ of the first kind,

$$
Da_I = \frac{\text{reaction time}}{\text{time constant of the reaction}} .
$$

The Damköhler number can be determined from the reaction rate at initial conditions $r_0$ and the
initial concentration of the limiting component $c\un{1,0}$, giving the mathematical definition

$$
Da_I = \frac{r_0 t}{c\un{1,0}} = c\un{1,0}^{\,n-1} k t ,
\qquad\text{with}\qquad r_0 = k\, c\un{1,0}^{\,n} .
$$ (eq-damkohler)

The reaction time is $t$, and the time constant of the reaction is $c\un{1,0}/r_0$ — the time it
would take to consume all of the limiting reactant if the rate stayed at its initial value. The
useful thing is that $Da_I$ describes the extent of reaction independently of the details of the
kinetics. Looking at the expression above, the right-hand side contains $Da_I = kt$ for a
first-order reaction, so the dimensionless material balance solves to

$$
c\un{A} = c\un{A,0}\exp\left(-kt\right)
\quad \Rightarrow \quad
f\un{A} = \exp\left(-Da_I\right) .
$$ (eq-f-da-first-order)

:::{admonition} Discussion
:class: seealso
How long does it take for a reaction to reach $f = 0$, or a conversion $X = 1$, when $n \geq 1$?
:::

:::{figure} ../figures/rate_law_f_Da.png
:label: fig-n-da-orders
:alt: Plot of the remaining fraction f, equal to c_A over c_A0, on the vertical axis from 0 to 1, against the first Damkohler number on the horizontal axis from 0 to 2. Five curves for different reaction orders all start together at f equals 1. In order of decreasing steepness they are n equals minus one half, which reaches zero at about Da 0.67; n equals 0, a straight line reaching zero at Da 1; n equals one half, reaching zero near Da 2; n equals 1, decaying exponentially to about 0.13 at Da 2; and n equals 2, the shallowest, still at about 0.33 at Da 2. Only the orders below one reach f equals zero at finite Da.
:width: 80%

Fractional conversion of species A as a function of the $Da_I$ number for various reaction orders
of irreversible reactions.
:::

### Second-order irreversible reactions, first order in both reactants

source: ReactionKinetics.tex L511
For Case II the solution of the material balance is more complicated. Formulate the mass balances
for the two reactants:

$$
\begin{aligned}
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} &= -k\, c\un{A} c\un{B} \\
\frac{\mathrm{d}c\un{B}}{\mathrm{d}t} &= -k\, c\un{A} c\un{B} .
\end{aligned}
$$

This is a little more challenging to integrate, since both $c\un{A}$ and $c\un{B}$ depend on time.
There are different ways to approach it.

**1.** Use the stoichiometry of the reaction. If the concentration of A is reduced by $x$, then B is
reduced by the same amount:

$$
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t}
= -k\left(c\un{A,0} - x\right)\left(c\un{B,0} - x\right)
= -\frac{\mathrm{d}x}{\mathrm{d}t} .
$$

At any time $c\un{A} = c\un{A,0} - x$, so $\mathrm{d}c\un{A}/\mathrm{d}t = -\mathrm{d}x/\mathrm{d}t$.

**2.** Subtract the material balance of B from that of A, as we did above:

$$
\frac{\mathrm{d}\left(c\un{A} - c\un{B}\right)}{\mathrm{d}t} = 0 .
$$

Both concentrations always change in the same way, so their difference is constant. We can
therefore solve for the concentration of B in terms of A,

$$
c\un{B} = c\un{A} - c\un{A,0} + c\un{B,0} ,
$$

which leads to

$$
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t}
= -k\, c\un{A}\left(c\un{A} - c\un{A,0} + c\un{B,0}\right) .
$$

Separating the variables and integrating,

$$
\int_{c\un{A,0}}^{c\un{A}}
\frac{\mathrm{d}c\un{A}}{c\un{A}\left(c\un{A} - c\un{A,0} + c\un{B,0}\right)} = -kt .
$$

This integral is rather challenging, but the method of partial fractions gives

$$
\frac{1}{c\un{B,0} - c\un{A,0}}\left[
\int_{c\un{A,0}}^{c\un{A}} \frac{\mathrm{d}c\un{A}}{c\un{A}}
- \int_{c\un{A,0}}^{c\un{A}} \frac{\mathrm{d}c\un{A}}{c\un{A} + c\un{B,0} - c\un{A,0}}
\right] = -kt ,
$$

resulting finally in

$$
\frac{1}{c\un{B,0} - c\un{A,0}}\left[
\ln\left(\frac{c\un{A}}{c\un{A,0}}\right)
- \ln\left(\frac{c\un{A} + c\un{B,0} - c\un{A,0}}{c\un{B,0}}\right)
\right] = -kt .
$$ (eq-second-order-caseii)

Phew — that is a lot of work for a rather simple rate law. The result is more complicated than for
Case I, since it depends on *both* initial concentrations. We can simplify by defining the ratio of
starting concentrations,

$$
\kappa\un{B} = \frac{c\un{B,0}}{c\un{A,0}} ,
$$ (eq-kappa-b)

and recasting the result in dimensionless form using $Da_I$. The fraction of reactant A remaining
is then

$$
f = \frac{c\un{A}}{c\un{A,0}}
  = \frac{\kappa\un{B} - 1}
         {\kappa\un{B}\exp\left[\left(\kappa\un{B} - 1\right)Da_I\right] - 1} .
$$ (eq-f-da-caseii)

:::{figure} ../figures/second_order_rate_law_f_Da.png
:label: fig-n-da-caseii
:alt: Plot of the remaining fraction f, equal to c_A over c_A0, on the vertical axis from 0 to 1, against the first Damkohler number on the horizontal axis from 0 to 2. Four curves for different initial reactant ratios all start at f equals 1 and decay. The larger the excess of B, the faster the decay: kappa_B equals 0.2 falls only to about 0.61 at Da 2; kappa_B equals 0.8 falls to about 0.43; kappa_B equals 2 falls to about 0.07; and kappa_B equals 5 has essentially reached zero by Da 1.5.
:width: 80%

Fractional conversion of species A as a function of the $Da_I$ number for a second-order rate law
$r = k c\un{A} c\un{B}$, for various initial ratios of the reactants $\kappa$.
:::

In dimensional form, the same result can be written as

$$
c\un{A} = \left(c\un{A,0} - c\un{B,0}\right)
\left[1 - \frac{c\un{B,0}}{c\un{A,0}}
\exp\left(\left(c\un{B,0} - c\un{A,0}\right)kt\right)\right]^{-1} .
$$ (eq-caseii-dimensional)

The full step-by-step derivation is uploaded to Canvas; you can also try integrating it yourself.

## Interlude: solving ordinary differential equations numerically

source: ReactionKinetics.tex L581
There are other rate laws for which we can find analytical solutions. In most cases, however, we
have to rely on numerical methods to solve our material balances — especially for more complex
systems with multiple reactions. Luckily there is a wide range of numerical tools available, and
you will start practicing with some of them in the next homework.

You will be using MATLAB in the homework assignment. When working with MATLAB you will frequently
use the `ode45` method to solve your system of ordinary differential equations. You can use that
solver without knowing the fundamentals behind it, but it is of course better if you do. Here we
look at another method, somewhat easier, used long before computers existed: the **Euler–Cauchy
method**.

$$
y(x_{i+1}) = y(x_i) + \int_{x_i}^{x_{i+1}} f(x, y(x))\, \mathrm{d}x .
$$

Our equation is

$$
\underbrace{\frac{\mathrm{d}c\un{A}}{\mathrm{d}t}}_{y'}
= \underbrace{-k\, c\un{A}}_{f(x,y)} .
$$

The integral can be approximated using the rectangle method,

$$
\int_{x_i}^{x_{i+1}} f(x, y(x))\, \mathrm{d}x \approx h \cdot f(x_i, y_i) ,
\qquad\text{with}\quad h = x_{i+1} - x_i .
$$

:::{figure} ../figures/Eulercauchy.png
:label: fig-euler-cauchy
:alt: Two hand-drawn panels. Left, y against x: between grid lines at x_i and x_i+1 two nearly straight lines rise from the same point at y_i, the upper one labelled true value and the lower one the Euler step, ending slightly below it at y_i+1. Right, the integrand f against x over the same interval: the curve rises from f_i to f_i+1, and the area beneath it from x_i to x_i+1 is hatched as a rectangle of height f_i, annotated that the geometrical integral is the slope.
:width: 80%

Graphical illustration of the Euler–Cauchy method. Over one step the integrand is held at its value
at $x_i$, so the area under the curve is approximated by a rectangle of width $h$.
:::

The approximate solution is therefore carried out via

$$
y_{i+1} = y_i + h \cdot f(x_i, y_i) ,
$$ (eq-euler-cauchy)

where the step size $h$ has to be chosen. Take the scenario
$c\un{A,0} = 1\ \mathrm{mol\,m^{-3}}$, $k = 0.02\ \mathrm{s^{-1}}$, and a step size of
$h = 20\ \mathrm{s}$. We can do the first step manually:

$$
\begin{aligned}
y_{i+1} &= y_i + h \cdot f(x_i, y_i)
&&\text{with}\quad f(x_i, y_i) = -k\, c\un{A} \\
\text{At } t = 0: \quad & c\un{A,i} = c\un{A,0} \\
\text{At time step 1:} \quad & c\un{A,i+1} = c\un{A,0} + h\left(-k\, c\un{A,0}\right) \\
& c\un{A,i+1} = 1\ \mathrm{mol\,m^{-3}}
  + 20\ \mathrm{s}\left(-0.02\ \mathrm{s^{-1}} \cdot 1\ \mathrm{mol\,m^{-3}}\right)
  = 0.600\ \mathrm{mol\,m^{-3}} .
\end{aligned}
$$

The exact solution at $t = 20\ \mathrm{s}$ would be

$$
\begin{aligned}
c\un{A} &= c\un{A,0}\exp\left(-kt\right) \\
c\un{A} &= 1\ \mathrm{mol\,m^{-3}}
  \exp\left(-0.02\ \mathrm{s^{-1}} \cdot 20\ \mathrm{s}\right)
  = 0.670\ \mathrm{mol\,m^{-3}} .
\end{aligned}
$$

(sec-reversible-first-order)=
## Reversible, first-order reaction

source: ReactionKinetics.tex L634
So far we have only looked at irreversible reactions. Basically all reactions, however, are
equilibrium reactions; in some cases the equilibrium simply lies so far to one side that we can
treat the reaction as irreversible. For reversible reactions it is also possible to derive an
analytical expression for simple systems.

Assume a reversible reaction in a batch reactor,

$$
\ce{A <=> B} ,
$$

with the following rate laws, assuming mass action kinetics:

$$
r\un{fwd} = k\un{fwd}\, c\un{A} , \qquad\qquad r\un{rev} = k\un{rev}\, c\un{B} .
$$

The net rate of reaction is the difference between the forward and reverse rates — the fifth
Boudart rule:

$$
r = r\un{fwd} - r\un{rev} = k\un{fwd}\, c\un{A} - k\un{rev}\, c\un{B} .
$$

Implementing this in the batch-reactor material balance,

$$
\begin{aligned}
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t}
  &= \nu\un{A} \underbrace{\left(r\un{fwd} - r\un{rev}\right)}_{r} \\
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t}
  &= \nu\un{A}\left(k\un{fwd}\, c\un{A} - k\un{rev}\, c\un{B}\right) \\
\frac{\mathrm{d}c\un{B}}{\mathrm{d}t}
  &= \nu\un{B}\left(r\un{fwd} - r\un{rev}\right) \\
\frac{\mathrm{d}c\un{B}}{\mathrm{d}t}
  &= \nu\un{B}\left(k\un{fwd}\, c\un{A} - k\un{rev}\, c\un{B}\right) .
\end{aligned}
$$

These ODEs are again challenging to solve, since each depends on the concentration of both A and B.
The mass balance simplifies matters. Mass is constant at all times, and assuming constant volume,

$$
c\un{A} + c\un{B} = c\un{A,0} + c\un{B,0}
\quad \Rightarrow \quad
c\un{B} = c\un{A,0} - c\un{A}
\quad \text{if } c\un{B,0} = 0 .
$$

The mass balance for A is therefore

$$
\begin{aligned}
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t}
  &= -k\un{fwd}\, c\un{A} + k\un{rev}\left(c\un{A,0} - c\un{A}\right) \\
  &= -\left(k\un{fwd} + k\un{rev}\right)c\un{A} + k\un{rev}\, c\un{A,0} ,
\end{aligned}
$$

and the analytical solution to this differential equation is

$$
c\un{A} = \frac{k\un{rev}
  + k\un{fwd}\exp\left(-\left(k\un{fwd} + k\un{rev}\right)t\right)}
  {k\un{fwd} + k\un{rev}}\, c\un{A,0} .
$$ (eq-reversible-first-order-solution)

The concentration profile for B follows from the mass balance:

$$
\begin{aligned}
c\un{B} &= c\un{A,0} - c\un{A}
 = c\un{A,0} - \frac{k\un{rev}
   + k\un{fwd}\exp\left(-\left(k\un{fwd} + k\un{rev}\right)t\right)}
   {k\un{fwd} + k\un{rev}}\, c\un{A,0} \\
c\un{B} &= c\un{A,0}\left[1 - \frac{k\un{rev}
   + k\un{fwd}\exp\left(-\left(k\un{fwd} + k\un{rev}\right)t\right)}
   {k\un{fwd} + k\un{rev}}\right] .
\end{aligned}
$$

:::{figure} ../figures/EquilProfile.png
:label: fig-reversible-reaction
:alt: Hand-drawn plot of concentration against time. Starting from c_A0 at the top left, a solid curve for A falls and levels off; a solid curve for B rises from zero and levels off at the same value, the two meeting at the equilibrium concentration marked on the vertical axis. This symmetric pair is the case where the forward and reverse rate constants are equal. A second, dashed pair of curves shows the case where the forward rate constant exceeds the reverse one: A falls further and B rises higher, so equilibrium lies toward the product side.
:width: 75%

Concentration profiles for a reversible first-order reaction $\ce{A <=> B}$ in a constant-volume
batch reactor, approaching the equilibrium values set by $K = k\un{fwd}/k\un{rev}$.
:::

:::{admonition} Discussion
:class: seealso
What happens to the concentration profiles as $t \to \infty$?
:::

For $t \to \infty$, $\exp(-\infty) \to 0$, and the equation reduces to

$$
c\un{A,eq} = \frac{k\un{rev}\, c\un{A,0}}{k\un{fwd} + k\un{rev}} ,
$$

and

$$
c\un{B,eq} = \frac{k\un{fwd}\, c\un{A,0}}{k\un{fwd} + k\un{rev}} .
$$

Dividing these by each other leads to

$$
K = \frac{c\un{B,eq}}{c\un{A,eq}} = \frac{k\un{fwd}}{k\un{rev}} .
$$ (eq-k-from-kratios)

The concentration profiles at equilibrium depend only on the ratio $k\un{fwd}/k\un{rev}$.

:::{admonition} Discussion
:class: seealso
What *is* this ratio — have we seen it before in this course?
:::

Consider a general reversible chemical reaction that follows mass action kinetics, where the
reaction orders are the stoichiometric coefficients. At equilibrium the extent of reaction no
longer changes,

$$
\frac{\mathrm{d}\xi}{\mathrm{d}t} = r = 0 ,
$$

so the net reaction rate is zero, which means the forward and reverse rates are balanced:

$$
\begin{aligned}
r\un{fwd} &= r\un{rev} \\
k\un{fwd} \prod_i c_i^{\nu_i'} &= k\un{rev} \prod_i c_i^{\nu_i''} \\
\frac{k\un{fwd}}{k\un{rev}} &= \frac{\prod_i c_i^{\nu_i''}}{\prod_i c_i^{\nu_i'}} \\
\frac{k\un{fwd}}{k\un{rev}} &= \prod_i c_i^{\nu_i} \\
\frac{k\un{fwd}}{k\un{rev}} &= K\un{c} .
\end{aligned}
$$ (eq-kc-from-kratios)

The ratio of the rate constants is in fact the concentration equilibrium constant $K\un{c}$ from
[](#ch-thermodynamics) — kinetics and thermodynamics are linked through this relation. On this
note, the fifth Boudart rule is actually a law: **all reactions are equilibrium reactions and
reversible.**

## Summary

source: ReactionKinetics.tex L738
- The rate of reaction $r = (1/V)\,\mathrm{d}\xi/\mathrm{d}t$ is intensive and species-independent;
  the species production rate is $r_i = \nu_i r$, [](#eq-ri-nui-r).
- Boudart's five rules: (1) the rate decreases monotonically with extent of reaction, autocatalytic
  reactions being the exception; (2) the rate-constant assumption $r = k(T) f(c_i)$; (3) the
  Arrhenius form $k = A\exp(-E\un{a}/RT)$, [](#eq-arrhenius); (4) power-law concentration
  dependence $f(c_i) = \prod_i c_i^{\alpha_i}$; (5) net rate $= r\un{fwd} - r\un{rev}$.
- The reaction orders $\alpha_i$ are empirical and generally differ from the stoichiometric
  coefficients — except for elementary reactions, where the law of mass action gives
  $\alpha_i = \nu_i'$ forward and $\alpha_i = \nu_i''$ reverse, [](#eq-mass-action).
- Integrated rate laws for a constant-volume batch reactor:
  $c\un{A} = c\un{A,0}\exp(-kt)$ for first order, [](#eq-first-order-solution), with
  $t_{1/2} = \ln 2/k$, [](#eq-half-life-first-order); $c\un{A} = c\un{A,0} - kt$ for zero order,
  [](#eq-zero-order-solution); and $c\un{A} = c\un{A,0}/(1 + k c\un{A,0} t)$ for second order
  Case I, [](#eq-second-order-casei).
- The Damköhler number $Da_I = c\un{1,0}^{\,n-1} k t$, [](#eq-damkohler), collapses kinetic data of
  like reaction order onto a single dimensionless curve of $f$ vs. $Da_I$,
  [](#eq-f-da-first-order).
- For systems without an analytical solution, the Euler–Cauchy method approximates the ODE by
  $y_{i+1} = y_i + h \cdot f(x_i, y_i)$, [](#eq-euler-cauchy).
- For a reversible first-order reaction $\ce{A <=> B}$, the equilibrium ratio
  $c\un{B,eq}/c\un{A,eq} = k\un{fwd}/k\un{rev} = K$, [](#eq-k-from-kratios); for general mass action
  kinetics, $k\un{fwd}/k\un{rev} = K\un{c}$, [](#eq-kc-from-kratios). Kinetics and thermodynamics
  are linked through the rate-constant ratio.
 -->

