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

For the irreversible reaction $\ce{A -> B}$ in a constant-volume batch reactor, the mass balance
and its solution are

$$
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} = -k\,c\un{A},
\qquad
c\un{A}(t) = c\un{A,0}\,\mathrm{e}^{-kt},
\qquad
c\un{B}(t) = c\un{B,0} + c\un{A,0}\left(1 - \mathrm{e}^{-kt}\right),
$$

with the rate constant set by the Arrhenius equation, $k(T) = A\,\mathrm{e}^{-E\un{a}/RT}$.

The concentration profiles depend on temperature only through $k$ — but that dependence is
exponential, and it is worth developing intuition for how violent it is. Drag the temperature
slider below and watch the profiles: over a span of 80 K the reaction goes from barely proceeding
to essentially complete within the same 200 s window.

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

## Reversible, first-order reaction

<!-- source: ReactionKinetics.tex L634 -->

## Summary

<!-- source: ReactionKinetics.tex L738 -->
