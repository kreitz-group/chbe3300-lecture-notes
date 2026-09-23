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

The examples in this chapter all use the same measurements: the concentration of A in an isothermal
batch reactor, sampled every $0.1\ \mathrm{h}$ for two hours, from runs at $50$, $65$, and
$80\ \mathrm{^\circ C}$. The first two examples use only the run at $65\ \mathrm{^\circ C}$.

:::{admonition} Live example
:class: seealso
Choose a postulated reaction order $n$. The plot shows the data in the linearized form
[](#eq-integral-linearized) — for $n = 1$ that is $\ln(c\un{A,0}/c\un{A})$ vs. $t$ — together with a
straight-line fit. Only the right order puts the points on the line. For a wrong order the points
curve systematically around the fit, even when $R^2$ still looks respectable. The MATLAB code that
produces the same result follows underneath.
:::

```{marimo} python
import marimo as mo
import numpy as np
import matplotlib.pyplot as plt
```

```{marimo} python
exp_time = np.round(np.arange(21) * 0.1, 1)  # h
exp_conc = np.array([  # mol/m^3, one row per temperature
    [100, 97.14, 94.36, 91.66, 89.04, 86.50, 84.02, 81.62, 79.29, 77.02, 74.82,
     72.68, 70.60, 68.58, 66.62, 64.71, 62.86, 61.07, 59.32, 57.62, 55.98],
    [100, 90.48, 81.86, 74.06, 67.01, 60.63, 54.85, 49.63, 44.90, 40.63, 36.76,
     33.26, 30.09, 27.22, 24.63, 22.29, 20.16, 18.24, 16.51, 14.93, 13.51],
    [100, 73.29, 53.71, 39.37, 28.85, 21.14, 15.50, 11.38, 8.32, 6.10, 4.47,
     3.28, 2.40, 1.76, 1.29, 0.95, 0.69, 0.51, 0.37, 0.27, 0.20],
])
exp_temps = np.array([50.0, 65.0, 80.0]) + 273  # K
exp_R = 8.314  # J/(mol K)
```

```{marimo} python
int_n = mo.ui.slider(
    steps=[0, 0.5, 1, 1.5, 2, 3],
    value=0,
    label="Postulated reaction order n",
    show_value=True,
)
int_n
```

```{marimo} python
int_order = float(int_n.value)
int_conc = exp_conc[1]  # run at 65 °C

# Linearized integrated rate law, y = m*t + b
if int_order == 1:
    int_y = np.log(int_conc[0] / int_conc)
    int_ylabel = r"$\ln(c_\mathrm{A0}/c_\mathrm{A})$"
    int_kunit = r"\mathrm{h^{-1}}"
else:
    int_y = int_conc ** (1 - int_order)
    int_ylabel = rf"$c_\mathrm{{A}}^{{{1 - int_order:g}}}$"
    int_kunit = rf"(\mathrm{{mol\,m^{{-3}}}})^{{{1 - int_order:g}}}\,\mathrm{{h^{{-1}}}}"

int_m, int_b = np.polyfit(exp_time, int_y, 1)
int_fit = int_m * exp_time + int_b
int_r2 = 1 - np.sum((int_y - int_fit) ** 2) / np.sum((int_y - int_y.mean()) ** 2)
int_k = int_m if int_order == 1 else -int_m / (1 - int_order)

mo.md(f"$n = {int_order:g}$: $k = {int_k:.4g}\\ {int_kunit}$, $R^2 = {int_r2:.4f}$.")
```

```{marimo} python
fig_int, ax_int = plt.subplots(figsize=(5.5, 3.9))

ax_int.plot(exp_time, int_y, "o", ms=6, label="Exp. data")
ax_int.plot(exp_time, int_fit, lw=2.0, label="Fit")

ax_int.set_xlim(0, 2.1)
ax_int.set_xlabel("time (h)", fontsize=13)
ax_int.set_ylabel(int_ylabel, fontsize=13)
ax_int.tick_params(labelsize=12, width=1.2, length=5)
for sp_int in ax_int.spines.values():
    sp_int.set_linewidth(1.2)
ax_int.legend(loc="best", fontsize=13, frameon=False)

fig_int.tight_layout()
fig_int
```

In MATLAB, set `n` to the order you want to test and run:

```matlab
% Experimental data
time = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, ...
        1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2];     % h
conc = [100, 90.48, 81.86, 74.06, 67.01, 60.63, 54.85, 49.63, ...
        44.90, 40.63, 36.76, 33.26, 30.09, 27.22, 24.63, 22.29, ...
        20.16, 18.24, 16.51, 14.93, 13.51];                       % mol/m3

n = 0;  % postulated reaction order

% Plot the data in the linearized form that belongs to this order
if n == 1
    y = log(conc(1)./conc);   % ln(cA0/cA) = k*t
    y_label = '$\mathrm{ln(c_{A0}/c_A)\ \left(1\right)}$';
else
    y = conc.^(1-n);          % cA^(1-n) = cA0^(1-n) - (1-n)*k*t
    y_label = sprintf('$\\mathrm{c_A^{%g}}$', 1-n);
end

% Linear regression y = m*x + b
mdl = fitlm(time, y);
b = mdl.Coefficients.Estimate(1);  % the built-in function puts the intercept first
m = mdl.Coefficients.Estimate(2);

if n == 1
    k = m;
else
    k = -m/(1-n);
end
fprintf('n = %g: k = %.4g, R^2 = %.4f\n', n, k, mdl.Rsquared.Ordinary)

y_fitted = m*time + b;

figure('Units','centimeters','Position',[5 5 14 10])
hold on

plot(time, y, ...
    'LineWidth', 2.0, 'LineStyle','None','DisplayName','Exp. Data',...
    'Marker','o','MarkerFaceColor','auto')

plot(time, y_fitted, ...
    'LineWidth', 2.0, 'LineStyle','-','DisplayName','Fit')

xlim([0 2.1])

xlabel('$\mathrm{time\ (h)}$','Interpreter','latex','FontSize',16)
ylabel(y_label,'Interpreter','latex','FontSize',16)

set(gca, ...
    'FontName','lmodern', ...
    'FontSize',16, ...
    'LineWidth',1.5, ...
    'TickLength',[0.015 0.015], ...
    'Box','on')

legend('Interpreter','latex','Location','best','FontSize',16)

grid off
hold off
```

:::{tip} Check yourself
Only $n = 1$ gives $R^2 = 1.0000$, with $k = 1.00\ \mathrm{h^{-1}}$. Notice that $n = 0.5$ and
$n = 1.5$ still reach $R^2 = 0.982$, so a high $R^2$ alone does not confirm a rate law. The
systematic curvature of the points around the line is the better test.
:::

## Differential method

<!-- source: Experiments.tex L145 -->

The second classical approach is the differential analysis of kinetic data. This method is very
useful for complex rate laws, where the integral form is hard to obtain analytically. The procedure
is:

- Determine $\mathrm{d}c\un{A}/\mathrm{d}t$ from the $c\un{A}$ vs. $t$ data.
- Postulate a rate law — or obtain it from a mechanistic analysis using the PSSA or the QEA —
  linearize it, and plot the data appropriately.

:::{figure} ../figures/DiffMethod.png
:label: fig-differential
:alt: Hand-drawn plot of the concentration of A against time. Experimental points marked with crosses fall along a decaying curve, and short straight tangent lines are drawn through several of the points, annotated slope at each point.
:width: 55%

Reaction rates can be determined from concentration-versus-time data by differentiation.
:::

The derivative can be approximated by a finite-difference quotient, which gives the reaction rate as

$$
r = \frac{1}{\nu_i}\frac{\Delta c_i}{\Delta t} .
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

For an $n$th-order rate law, taking the logarithm of the mass balance gives

$$
\begin{aligned}
\frac{\mathrm{d}c\un{A}}{\mathrm{d}t} &= -k\, c\un{A}^n \\
\underbrace{\ln\left(-\frac{\mathrm{d}c\un{A}}{\mathrm{d}t}\right)}_{y}
&= \underbrace{n \ln(c\un{A})}_{mx} + \underbrace{\ln(k)}_{b} ,
\end{aligned}
$$ (eq-differential-linearized)

a straight line in $\ln(c\un{A})$. Using the differential method, both $n$ and $k$ can be extracted from a single
plot, whereas the integral method generally requires several attempts.

:::{figure} ../figures/DifferentialMethod.png
:label: fig-differential-linearized
:alt: Hand-drawn plot of the natural logarithm of minus dc_A by dt on the vertical axis against the natural logarithm of c_A on the horizontal axis. Scattered data points fall along a straight rising line, annotated with slope equal to n and y-intercept equal to the natural logarithm of k.
:width: 65%

Graphical illustration of the differential method.
:::

:::{admonition} Live example
:class: seealso
The same run at $65\ \mathrm{^\circ C}$, analysed with the differential method. Rates are computed
with the forward difference, [](#eq-forward-difference), called the Newton method in the code, and
with the central difference, [](#eq-central-difference). Each set of rates is then fitted with
[](#eq-differential-linearized). Use the slider to space the samples further apart: both schemes
keep returning $n = 1$, but the forward difference underestimates $k$ more and more. It assigns the
slope of the secant over $[t, t + \Delta t]$ to the left end of the interval, where the rate is
highest. The MATLAB code that produces the same result follows underneath.
:::

```{marimo} python
dif_dt = mo.ui.slider(
    steps=[0.1, 0.2, 0.5],
    value=0.1,
    label="Sampling interval Δt (h)",
    show_value=True,
)
dif_dt
```

```{marimo} python
dif_every = int(round(dif_dt.value / 0.1))
dif_t = exp_time[::dif_every]
dif_c = exp_conc[1][::dif_every]  # run at 65 °C

# Forward (Newton) difference: rate at t[i] from points i and i+1
dif_r_fwd = -(dif_c[1:] - dif_c[:-1]) / (dif_t[1:] - dif_t[:-1])
# Symmetric difference: rate at t[i] from points i-1 and i+1
dif_r_sym = -(dif_c[2:] - dif_c[:-2]) / (dif_t[2:] - dif_t[:-2])

# Linearized rate law ln(r) = n ln(c) + ln(k)
dif_n_fwd, dif_lnk_fwd = np.polyfit(np.log(dif_c[:-1]), np.log(dif_r_fwd), 1)
dif_n_sym, dif_lnk_sym = np.polyfit(np.log(dif_c[1:-1]), np.log(dif_r_sym), 1)

mo.md(
    f"Newton method: $n = {dif_n_fwd:.3f}$, "
    f"$k = {np.exp(dif_lnk_fwd):.3f}\\ \\mathrm{{h^{{-1}}}}$. "
    f"Symmetric difference: $n = {dif_n_sym:.3f}$, "
    f"$k = {np.exp(dif_lnk_sym):.3f}\\ \\mathrm{{h^{{-1}}}}$."
)
```

```{marimo} python
fig_dif, (ax_dif1, ax_dif2) = plt.subplots(1, 2, figsize=(9, 3.6))

ax_dif1.plot(dif_t[:-1], dif_r_fwd, "o", ms=6, color="C0", label="Newton method")
ax_dif1.plot(dif_t[1:-1], dif_r_sym, "s", ms=6, color="C1", label="Sym. diff. method")
ax_dif1.set_xlim(0, 2.1)
ax_dif1.set_ylim(0, 100)
ax_dif1.set_xlabel("time (h)", fontsize=13)
ax_dif1.set_ylabel(r"$r$ (mol m$^{-3}$ h$^{-1}$)", fontsize=13)
ax_dif1.legend(loc="upper right", fontsize=11, frameon=False)

dif_x = np.log(dif_c)
ax_dif2.plot(np.log(dif_c[:-1]), np.log(dif_r_fwd), "o", ms=6, color="C0")
ax_dif2.plot(np.log(dif_c[1:-1]), np.log(dif_r_sym), "s", ms=6, color="C1")
ax_dif2.plot(dif_x, dif_n_fwd * dif_x + dif_lnk_fwd, "-", lw=2.0, color="C0")
ax_dif2.plot(dif_x, dif_n_sym * dif_x + dif_lnk_sym, "--", lw=2.0, color="C1")
ax_dif2.set_xlabel(r"$\ln(c_\mathrm{A})$", fontsize=13)
ax_dif2.set_ylabel(r"$\ln(r_\mathrm{A})$", fontsize=13)

for ax_d in (ax_dif1, ax_dif2):
    ax_d.tick_params(labelsize=12, width=1.2, length=5)
    for sp_dif in ax_d.spines.values():
        sp_dif.set_linewidth(1.2)

fig_dif.tight_layout()
fig_dif
```

In MATLAB, set `every` to the sampling interval you want and run:

```matlab
% Experimental data
time = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, ...
        1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2];     % h
conc = [100, 90.48, 81.86, 74.06, 67.01, 60.63, 54.85, 49.63, ...
        44.90, 40.63, 36.76, 33.26, 30.09, 27.22, 24.63, 22.29, ...
        20.16, 18.24, 16.51, 14.93, 13.51];                       % mol/m3

every = 1;  % use every n-th sample: 1 -> 0.1 h, 2 -> 0.2 h, 5 -> 0.5 h
t = time(1:every:end);
c = conc(1:every:end);

% Calculate the rates using numerical differentiation
r_Newton = zeros(1, numel(t)-1);  % forward (Newton) difference
r_SDM    = zeros(1, numel(t)-2);  % symmetric difference

for index = 1:numel(t)-1
    r_Newton(index) = -(c(index+1)-c(index))/(t(index+1)-t(index));
end

for index = 2:numel(t)-1
    r_SDM(index-1) = -(c(index+1)-c(index-1))/(t(index+1)-t(index-1));
end

% Linearized rate law: ln(r) = n*ln(cA) + ln(k)
mdl_Newton = fitlm(log(c(1:end-1)), log(r_Newton));
mdl_SDM    = fitlm(log(c(2:end-1)), log(r_SDM));

n_Newton = mdl_Newton.Coefficients.Estimate(2);
k_Newton = exp(mdl_Newton.Coefficients.Estimate(1));
n_SDM    = mdl_SDM.Coefficients.Estimate(2);
k_SDM    = exp(mdl_SDM.Coefficients.Estimate(1));

fprintf('Newton method:     n = %.3f, k = %.3f 1/h\n', n_Newton, k_Newton)
fprintf('Sym. diff. method: n = %.3f, k = %.3f 1/h\n', n_SDM, k_SDM)

% Rates vs. time
figure('Units','centimeters','Position',[5 5 14 10])
hold on

plot(t(1:end-1), r_Newton, ...
    'LineWidth', 2.0, 'LineStyle','None','DisplayName','Newton Method',...
    'Marker','o')
plot(t(2:end-1), r_SDM, ...
    'LineWidth', 2.0, 'LineStyle','None','DisplayName','Sym. Diff. Method',...
    'Marker','s')

xlim([0 2.1])
ylim([0 100])

xlabel('$\mathrm{time\ (h)}$','Interpreter','latex','FontSize',16)
ylabel('$\mathrm{r\ \left(mol\,m^{-3}\,h^{-1}\right)}$','Interpreter','latex','FontSize',16)

set(gca, ...
    'FontName','lmodern', ...
    'FontSize',16, ...
    'LineWidth',1.5, ...
    'TickLength',[0.015 0.015], ...
    'Box','on')

legend('Interpreter','latex','Location','northeast','FontSize',16)

grid off
hold off

% ln(r) vs. ln(cA) with the linear fits
figure('Units','centimeters','Position',[5 5 14 10])
hold on

plot(log(c(1:end-1)), log(r_Newton), ...
    'LineWidth', 2.0, 'LineStyle','None','DisplayName','Newton Method',...
    'Marker','o')
plot(log(c(2:end-1)), log(r_SDM), ...
    'LineWidth', 2.0, 'LineStyle','None','DisplayName','Sym. Diff. Method',...
    'Marker','s')
plot(log(c), n_Newton*log(c) + log(k_Newton), ...
    'LineWidth', 2.0, 'LineStyle','-','HandleVisibility','off')
plot(log(c), n_SDM*log(c) + log(k_SDM), ...
    'LineWidth', 2.0, 'LineStyle','--','HandleVisibility','off')

xlabel('$\mathrm{ln(c_A)}$','Interpreter','latex','FontSize',16)
ylabel('$\mathrm{ln(r_A)}$','Interpreter','latex','FontSize',16)

set(gca, ...
    'FontName','lmodern', ...
    'FontSize',16, ...
    'LineWidth',1.5, ...
    'TickLength',[0.015 0.015], ...
    'Box','on')

legend('Interpreter','latex','Location','northwest','FontSize',16)

grid off
hold off
```

:::{tip} Check yourself
At $\Delta t = 0.1\ \mathrm{h}$ the Newton method gives $k = 0.953$ and the symmetric difference
gives $1.003\ \mathrm{h^{-1}}$; at $\Delta t = 0.5\ \mathrm{h}$ they give $0.788$ and $1.043$. For
first-order data, $c\un{A} = c\un{A,0}e^{-kt}$, show that the forward difference returns
$(1 - e^{-k\Delta t})/\Delta t$ instead of $k$ and the central difference returns
$\sinh(k\Delta t)/\Delta t$. Evaluate both with $k = 1\ \mathrm{h^{-1}}$ from the integral method.
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

:::{admonition} Live example
:class: seealso
All three runs, at $50$, $65$, and $80\ \mathrm{^\circ C}$, analysed in two linear steps. First, the
differential method at each temperature: rates from central differences (MATLAB's `gradient`, which
falls back to one-sided differences at the two end points), then a fit of
[](#eq-differential-linearized) that gives $n$ and $\ln k$. Second, a fit of $\ln k$ vs. $1000/T$,
[](#eq-arrhenius-linear), that gives $E\un{a}$ from the slope and $A$ from the intercept. This
example has no controls; it is computed from the data in your browser, and the MATLAB code that
produces the same result follows underneath.
:::

```{marimo} python
arr_lnc = np.log(exp_conc)
arr_lnr = np.log(-np.gradient(exp_conc, exp_time, axis=1))

# Differential method at each temperature: ln(r) = n ln(c) + ln(k)
arr_fits = np.array([np.polyfit(arr_lnc[i], arr_lnr[i], 1) for i in range(3)])
arr_order, arr_lnk = arr_fits[:, 0], arr_fits[:, 1]

# Arrhenius fit: ln(k) = -(Ea/R) (1/T) + ln(A), with x = 1000/T
arr_x = 1000 / exp_temps
arr_slope, arr_int = np.polyfit(arr_x, arr_lnk, 1)
arr_Ea = -arr_slope * exp_R * 1000  # J/mol
arr_A = np.exp(arr_int)  # 1/h

mo.md(
    "Reaction orders: "
    + ", ".join(f"${n:.3f}$" for n in arr_order)
    + "; rate constants: "
    + ", ".join(f"${np.exp(lk):.3f}$" for lk in arr_lnk)
    + r"$\ \mathrm{h^{-1}}$. "
    + f"$E_\\mathrm{{a}} = {arr_Ea / 1000:.2f}\\ \\mathrm{{kJ\\,mol^{{-1}}}}$, "
    + f"$A = {arr_A:.2e}\\ \\mathrm{{h^{{-1}}}}$."
)
```

```{marimo} python
fig_arr, (ax_arr1, ax_arr2) = plt.subplots(1, 2, figsize=(9, 3.6))

for arr_i in range(3):
    ax_arr1.plot(arr_lnc[arr_i], arr_lnr[arr_i], "o", ms=5, color=f"C{arr_i}",
                 label=f"{exp_temps[arr_i] - 273:.0f} °C")
    ax_arr1.plot(arr_lnc[arr_i], arr_order[arr_i] * arr_lnc[arr_i] + arr_lnk[arr_i],
                 "-", lw=1.5, color=f"C{arr_i}")
ax_arr1.set_xlim(-2, 5)
ax_arr1.set_ylim(-2, 6)
ax_arr1.set_xlabel(r"$\ln(c)$", fontsize=13)
ax_arr1.set_ylabel(r"$\ln(r)$", fontsize=13)
ax_arr1.legend(loc="lower right", fontsize=11, frameon=False)

ax_arr2.plot(arr_x, arr_lnk, "ko", ms=7, label=r"$\ln(k)$")
ax_arr2.plot(arr_x, arr_slope * arr_x + arr_int, "k-", lw=1.5, label="Arrhenius fit")
ax_arr2.set_xlim(2.8, 3.15)
ax_arr2.set_ylim(-1.5, 1.5)
ax_arr2.set_xlabel(r"$1000/T$ (K$^{-1}$)", fontsize=13)
ax_arr2.set_ylabel(r"$\ln(k)$", fontsize=13)
ax_arr2.legend(loc="upper right", fontsize=11, frameon=False)

for ax_a in (ax_arr1, ax_arr2):
    ax_a.tick_params(labelsize=12, width=1.2, length=5)
    for sp_arr in ax_a.spines.values():
        sp_arr.set_linewidth(1.2)

fig_arr.tight_layout()
fig_arr
```

In MATLAB:

```matlab
% Time in hours
time = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, ...
        1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2];

% Concentration c in mol/m^3, one row per temperature
conc = [100, 97.14, 94.36, 91.66, 89.04, 86.50, 84.02, 81.62, 79.29, 77.02, 74.82, 72.68, 70.60, 68.58, 66.62, 64.71, 62.86, 61.07, 59.32, 57.62, 55.98;
        100, 90.48, 81.86, 74.06, 67.01, 60.63, 54.85, 49.63, 44.90, 40.63, 36.76, 33.26, 30.09, 27.22, 24.63, 22.29, 20.16, 18.24, 16.51, 14.93, 13.51;
        100, 73.29, 53.71, 39.37, 28.85, 21.14, 15.50, 11.38, 8.32, 6.10, 4.47, 3.28, 2.40, 1.76, 1.29, 0.95, 0.69, 0.51, 0.37, 0.27, 0.20];

% Temperature in Kelvin
temps = [50, 65, 80] + 273;

%% Rates by numerical differentiation
rate   = zeros(size(conc));
lnrate = zeros(size(conc));
lnconc = zeros(size(conc));

for i = 1:size(conc,1)
    rate(i,:)   = -gradient(conc(i,:), time);  % central differences inside, one-sided at the ends
    lnrate(i,:) = log(rate(i,:));
    lnconc(i,:) = log(conc(i,:));
end

%% Fit the linearized rate law ln(r) = n*ln(c) + ln(k) at each temperature
lnk   = zeros(size(conc,1),1);
order = zeros(size(conc,1),1);

for i = 1:size(conc,1)
    p = polyfit(lnconc(i,:), lnrate(i,:), 1);
    order(i) = p(1);   % slope = reaction order
    lnk(i)   = p(2);   % intercept = ln(k)
end

disp(table((temps-273)', order, exp(lnk), 'VariableNames', {'T_degC','n','k_per_h'}))

%% Arrhenius fit: ln(k) = -(Ea/R)*(1/T) + ln(A), with x = 1000/T
x = 1000 ./ temps;
y = lnk';

p = polyfit(x, y, 1);   % p(1) = slope, p(2) = intercept

R  = 8.314;             % J/(mol K)
Ea = -p(1) * R * 1000;  % J/mol (the factor 1000 undoes x = 1000/T)
A  = exp(p(2));         % same units as k, 1/h

fprintf('Ea = %.2f kJ/mol\n', Ea/1000);
fprintf('A  = %.2e 1/h\n', A);

%% Plotting
figure('Position',[100 100 1000 400])
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

% ln(r) vs ln(c) with the fitted lines
nexttile
hold on; grid on
for i = 1:size(conc,1)
    h = plot(lnconc(i,:), lnrate(i,:), 'o', 'LineWidth', 2, 'MarkerSize', 8, ...
        'DisplayName', sprintf('%d °C', temps(i)-273));
    plot(lnconc(i,:), order(i)*lnconc(i,:) + lnk(i), '-', ...
        'Color', h.Color, 'LineWidth', 1.5, 'HandleVisibility', 'off')
end
legend('Location','best')
xlim([-2 5])
ylim([-2 6])
xlabel('$\ln(c)$','Interpreter','latex')
ylabel('$\ln(r)$','Interpreter','latex')
set(gca,'FontSize',14)

% Arrhenius plot
nexttile
hold on; grid on
plot(1000./temps, lnk, 'ko', 'LineWidth', 2, 'MarkerSize', 9, 'DisplayName', 'ln(k)')
plot(x, polyval(p, x), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Arrhenius fit')
legend('Location','best')
xlim([2.8 3.15])
ylim([-1.5 1.5])
xlabel('$1000/T\ \mathrm{(K^{-1})}$','Interpreter','latex')
ylabel('$\ln(k)$','Interpreter','latex')
set(gca,'FontSize',14)
```

:::{tip} Check yourself
The orders come out at $n = 0.99$ rather than exactly 1, a small bias from the one-sided
differences at the ends of each run. The Arrhenius fit gives
$E\un{a} = 74.5\ \mathrm{kJ\,mol^{-1}}$ and $A = 3.4 \times 10^{11}\ \mathrm{h^{-1}}$. Check the
activation energy with the two-point formula, [](#eq-arrhenius-two-point), using the rate constants
at $50$ and $80\ \mathrm{^\circ C}$.
:::

## Nonlinear regression

<!-- source: Experiments.tex L229 -->

The procedures above are the classical, graphical approach. They predate modern computing and are
still very useful as sanity checks, because the linearization makes the parameter dependence
visually transparent. In a modern lab, however, parameters are typically extracted by nonlinear
regression directly against the ODE system.

For the simple reaction $\ce{A -> products}$ with an $n$th-order rate law in A, the batch-reactor
mass balance is

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
available — Nelder–Mead, Levenberg–Marquardt, BFGS, and others; in MATLAB, `fminsearch` minimizes a
scalar objective and `lsqnonlin` minimizes a sum of squared residuals. In most cases you will need
to provide parameter bounds and reasonable initial guesses, since the objective function often has
many local minima.

In the two examples below the residuals are concentrations: each iteration solves the batch-reactor
mass balance, [](#eq-nonlinear-ode) with $n = 1$, and compares the simulated $c\un{A}(t)$ with the
measured one.

:::{admonition} Live example
:class: seealso
A fit to a single run. With only one temperature, the data contain no information about the
temperature dependence, so $A$ and $E\un{a}$ cannot both be determined: $E\un{a}$ is fixed and only
$A$ is fitted, with `fminsearch`. Choose the run, then change the fixed activation energy. $A$ moves
by orders of magnitude, but the fit and its sum of squared residuals do not change at all. Every
pair $(A, E\un{a})$ that gives the same $k(T)$ describes this run equally well. The MATLAB code that
produces the same result follows underneath.
:::

```{marimo} python
from scipy.optimize import minimize_scalar

nls_run = mo.ui.dropdown(
    options={"50 °C": 0, "65 °C": 1, "80 °C": 2},
    value="50 °C",
    label="Run used for the fit",
)
nls_Ea = mo.ui.slider(
    steps=[60, 65, 70, 75, 80, 85, 90],
    value=75,
    label="Fixed activation energy Ea (kJ/mol)",
    show_value=True,
)
mo.vstack([nls_run, nls_Ea])
```

```{marimo} python
nls_i = nls_run.value
nls_c = exp_conc[nls_i]
nls_T = exp_temps[nls_i]
nls_Ea_J = nls_Ea.value * 1e3  # J/mol


def nls_sim(k0):
    """Solution of dc/dt = -k0 exp(-Ea/RT) c at the sample times (what ode45 computes)."""
    return nls_c[0] * np.exp(-k0 * np.exp(-nls_Ea_J / (exp_R * nls_T)) * exp_time)


def nls_cost(log_k0):
    return np.sum((nls_sim(10**log_k0) - nls_c) ** 2)


# Minimize over log10(k0) so the search is well scaled for any Ea
nls_opt = minimize_scalar(nls_cost, bounds=(0, 20), method="bounded",
                          options={"xatol": 1e-10})
nls_k0 = 10**nls_opt.x
nls_fit = nls_sim(nls_k0)

mo.md(
    f"$E_\\mathrm{{a}} = {nls_Ea.value:.2f}\\ \\mathrm{{kJ\\,mol^{{-1}}}}$ (fixed), "
    f"$A = {nls_k0:.2e}\\ \\mathrm{{h^{{-1}}}}$, "
    f"sum of squared residuals $= {nls_opt.fun:.4f}\\ (\\mathrm{{mol\\,m^{{-3}}}})^2$."
)
```

```{marimo} python
fig_nls, (ax_nls1, ax_nls2) = plt.subplots(1, 2, figsize=(9, 3.6))
nls_col = f"C{nls_i}"
nls_label = f"{nls_T - 273:.0f} °C"

ax_nls1.plot(exp_time, nls_fit, "-", lw=2.0, color=nls_col, label=f"{nls_label}, fit")
ax_nls1.plot(exp_time, nls_c, "o", ms=6, mfc="none", color=nls_col, label=f"{nls_label}, exp")
ax_nls1.set_xlim(0, 2.1)
ax_nls1.set_ylim(0, 105)
ax_nls1.set_xlabel("time (h)", fontsize=13)
ax_nls1.set_ylabel(r"$c$ (mol m$^{-3}$)", fontsize=13)
ax_nls1.legend(loc="upper right", fontsize=11, frameon=False)

ax_nls2.plot([0, 100], [0, 100], ":", color="0.5", lw=1.0)
ax_nls2.plot(nls_c, nls_fit, "o", ms=6, mfc="none", color=nls_col, label=nls_label)
ax_nls2.set_xlim(0, 105)
ax_nls2.set_ylim(0, 105)
ax_nls2.set_xlabel(r"$c_\mathrm{exp}$ (mol m$^{-3}$)", fontsize=13)
ax_nls2.set_ylabel(r"$c_\mathrm{fit}$ (mol m$^{-3}$)", fontsize=13)
ax_nls2.legend(loc="lower right", fontsize=11, frameon=False)

for ax_n in (ax_nls1, ax_nls2):
    ax_n.tick_params(labelsize=12, width=1.2, length=5)
    for sp_nls in ax_n.spines.values():
        sp_nls.set_linewidth(1.2)

fig_nls.tight_layout()
fig_nls
```

In MATLAB, save this as a script file, since it ends with local functions, and set
`data_set_eval` and `Ea_init` to the values you want:

```matlab
% Time in hours
time = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, ...
        1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2];

% Concentration c in mol/m^3, one row per temperature
conc = [100, 97.14, 94.36, 91.66, 89.04, 86.50, 84.02, 81.62, 79.29, 77.02, 74.82, 72.68, 70.60, 68.58, 66.62, 64.71, 62.86, 61.07, 59.32, 57.62, 55.98;
        100, 90.48, 81.86, 74.06, 67.01, 60.63, 54.85, 49.63, 44.90, 40.63, 36.76, 33.26, 30.09, 27.22, 24.63, 22.29, 20.16, 18.24, 16.51, 14.93, 13.51;
        100, 73.29, 53.71, 39.37, 28.85, 21.14, 15.50, 11.38, 8.32, 6.10, 4.47, 3.28, 2.40, 1.76, 1.29, 0.95, 0.69, 0.51, 0.37, 0.27, 0.20];

% Temperature in Kelvin
temps = [50, 65, 80] + 273;

%% Parameters
% Select which dataset to use for fitting
% 1: T = 50 °C, 2: T = 65 °C, 3: T = 80 °C
data_set_eval = 1;

% Initial guesses for the parameters before optimization
% k0: pre-exponential factor (Arrhenius), units: 1/h
% Ea: activation energy, units: J/mol (fixed during optimization)
k0_init = 3.5e11;
Ea_init = 75e3;

%% Optimization using fminsearch (Nelder-Mead)
t_data = time;
c_data = conc(data_set_eval, :);
T_data = temps(data_set_eval);

% Since Ea is fixed, only k0 is optimized. The anonymous function @(p)
% fixes Ea and the data while letting p = k0 vary; fminsearch minimizes
% the scalar cost returned by cost_Temp.
p0 = k0_init;
p_opt = fminsearch(@(p) cost_Temp(p, Ea_init, t_data, c_data, T_data), p0);

k0_fit = p_opt(1);
Ea_fit = Ea_init;   % Ea was not optimized, so we keep its initial value

fprintf('Ea = %.2f kJ/mol (fixed)\n', Ea_fit/1000);
fprintf('A  = %.2e 1/h\n', k0_fit);
fprintf('SSE = %.4f (mol/m3)^2\n', cost_Temp(k0_fit, Ea_fit, t_data, c_data, T_data));

%% Simulate with the fitted parameters
c_fit = sim_exp_Temp(t_data, c_data(1), k0_fit, Ea_fit, T_data);

label_T = [num2str(T_data - 273), ' °C'];

%% Plotting
figure('Position', [100, 100, 1000, 500]);

% Left: concentration vs time
subplot(1, 2, 1);
hold on; grid on; box on
plot(t_data, c_fit, '-', 'LineWidth', 2, 'DisplayName', [label_T, ', fit']);
plot(t_data, c_data, 'o', 'MarkerSize', 9, 'DisplayName', [label_T, ', exp']);
xlabel('$t \;/\; \mathrm{h}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$c \;/\; \mathrm{mol \; m^{-3}}$', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'best');
hold off;

% Right: parity plot (simulated vs experimental concentration)
subplot(1, 2, 2);
hold on; grid on; box on
plot(c_data, c_fit, 'o', 'LineWidth', 2, 'MarkerSize', 9, 'DisplayName', label_T);
plot([0 100], [0 100], 'k:', 'HandleVisibility', 'off');
xlabel('$c_\mathrm{exp} \;/\;\mathrm{mol\; m^{-3}}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$c_\mathrm{fit} \;/\;\mathrm{mol\; m^{-3}}$', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'best');
hold off;

function dcdt = balance_Temp(~, c, k0, Ea, T)
    % Batch-reactor mass balance for a first-order Arrhenius reaction:
    % dc/dt = -k0 * exp(-Ea / RT) * c
    R = 8.314;  % universal gas constant, J/(mol K)
    dcdt = -k0 * exp(-Ea / (R * T)) * c;
end

function sim = sim_exp_Temp(time, c_init, k0, Ea, T)
    % Solves the ODE and returns the simulated concentrations at the
    % experimental time points as a row vector, matching the data layout.
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [~, sol] = ode45(@(t, c) balance_Temp(t, c, k0, Ea, T), time, c_init, opts);
    sim = sol';
end

function cost = cost_Temp(p, Ea, time, data, T)
    % Sum of squared residuals between simulation and experiment
    k0 = p(1);
    sim = sim_exp_Temp(time, data(1), k0, Ea, T);
    residuals = sim - data;
    cost = sum(residuals.^2);
end
```

:::{tip} Check yourself
For the run at $50\ \mathrm{^\circ C}$, $E\un{a} = 75\ \mathrm{kJ\,mol^{-1}}$ gives
$A = 3.91 \times 10^{11}\ \mathrm{h^{-1}}$ and $E\un{a} = 60\ \mathrm{kJ\,mol^{-1}}$ gives
$A = 1.47 \times 10^{9}\ \mathrm{h^{-1}}$. Show that both pairs give the same rate constant,
$k(323\ \mathrm{K}) = 0.290\ \mathrm{h^{-1}}$.
:::

To separate $A$ from $E\un{a}$, the fit must see several temperatures at once. The residuals of all
runs are concatenated into one vector and minimized together with `lsqnonlin`.

:::{admonition} Live example
:class: seealso
A global fit of $A$ and $E\un{a}$ to all three runs simultaneously. Untick a temperature to leave
its run out of the fit: its curve is then a *prediction* from the other two, drawn dashed. How well a
fitted model predicts a run it has never seen is a much stronger test than how well it fits the
runs it was fitted to. The MATLAB code that produces the same result follows underneath.
:::

```{marimo} python
from scipy.optimize import least_squares

glb_use = mo.ui.array(
    [mo.ui.checkbox(value=True, label=f"{t - 273:.0f} °C") for t in exp_temps]
)
mo.hstack([mo.md("Runs used in the fit:"), glb_use.hstack(gap=1.5)], justify="start")
```

```{marimo} python
glb_sets = [i for i in range(3) if glb_use.value[i]]


def glb_sim(k0, Ea, i):
    """Solution of dc/dt = -k0 exp(-Ea/RT) c for run i (what ode45 computes)."""
    return exp_conc[i, 0] * np.exp(-k0 * np.exp(-Ea / (exp_R * exp_temps[i])) * exp_time)


def glb_resid(p):
    # p = [log10(k0), Ea in kJ/mol]; both are of order 10-100, so the fit is well scaled
    return np.concatenate([glb_sim(10 ** p[0], p[1] * 1e3, i) - exp_conc[i] for i in glb_sets])


if len(glb_sets) >= 2:
    glb_opt = least_squares(glb_resid, x0=[np.log10(3.1e11), 75.0],
                            bounds=([0, 5], [13, 150]), xtol=1e-12, ftol=1e-12)
    glb_k0, glb_Ea = 10 ** glb_opt.x[0], glb_opt.x[1] * 1e3
    glb_msg = mo.md(
        f"$E_\\mathrm{{a}} = {glb_Ea / 1000:.2f}\\ \\mathrm{{kJ\\,mol^{{-1}}}}$, "
        f"$A = {glb_k0:.2e}\\ \\mathrm{{h^{{-1}}}}$, "
        f"sum of squared residuals $= {np.sum(glb_opt.fun ** 2):.4f}\\ (\\mathrm{{mol\\,m^{{-3}}}})^2$."
    )
else:
    glb_k0 = glb_Ea = None
    glb_msg = mo.md(
        "Select at least two runs: a single temperature cannot separate $A$ from $E_\\mathrm{a}$."
    )
glb_msg
```

```{marimo} python
fig_glb, (ax_glb1, ax_glb2) = plt.subplots(1, 2, figsize=(9, 3.6))
ax_glb2.plot([0, 100], [0, 100], ":", color="0.5", lw=1.0)

for glb_i in range(3):
    glb_col = f"C{glb_i}"
    glb_label = f"{exp_temps[glb_i] - 273:.0f} °C"
    ax_glb1.plot(exp_time, exp_conc[glb_i], "o", ms=5, mfc="none", color=glb_col)
    if glb_k0 is not None:
        glb_fit = glb_sim(glb_k0, glb_Ea, glb_i)
        glb_used = glb_i in glb_sets
        ax_glb1.plot(exp_time, glb_fit, "-" if glb_used else "--", lw=2.0, color=glb_col,
                     label=f"{glb_label}, {'fit' if glb_used else 'predicted'}")
        ax_glb2.plot(exp_conc[glb_i], glb_fit, "o", ms=5, mfc="none", color=glb_col,
                     label=glb_label)

ax_glb1.set_xlim(0, 2.1)
ax_glb1.set_ylim(0, 105)
ax_glb1.set_xlabel("time (h)", fontsize=13)
ax_glb1.set_ylabel(r"$c$ (mol m$^{-3}$)", fontsize=13)
ax_glb2.set_xlim(0, 105)
ax_glb2.set_ylim(0, 105)
ax_glb2.set_xlabel(r"$c_\mathrm{exp}$ (mol m$^{-3}$)", fontsize=13)
ax_glb2.set_ylabel(r"$c_\mathrm{fit}$ (mol m$^{-3}$)", fontsize=13)
if glb_k0 is not None:
    ax_glb1.legend(loc="upper right", fontsize=10, frameon=False)
    ax_glb2.legend(loc="lower right", fontsize=10, frameon=False)

for ax_g in (ax_glb1, ax_glb2):
    ax_g.tick_params(labelsize=12, width=1.2, length=5)
    for sp_glb in ax_g.spines.values():
        sp_glb.set_linewidth(1.2)

fig_glb.tight_layout()
fig_glb
```

In MATLAB, save this as a script file and list the runs to fit in `fit_sets`. Note that the listing
fits $\log_{10} A$ rather than $A$: $A \approx 10^{11}$ and $E\un{a} \approx 10^{5}$ differ by six
orders of magnitude, and on that scale `lsqnonlin` stalls before it has moved $A$ from its initial
guess. Rescaling parameters to similar magnitudes is a routine part of setting up a fit.

```matlab
% Time in hours
time = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, ...
        1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2];

% Concentration c in mol/m^3, one row per temperature
conc = [100, 97.14, 94.36, 91.66, 89.04, 86.50, 84.02, 81.62, 79.29, 77.02, 74.82, 72.68, 70.60, 68.58, 66.62, 64.71, 62.86, 61.07, 59.32, 57.62, 55.98;
        100, 90.48, 81.86, 74.06, 67.01, 60.63, 54.85, 49.63, 44.90, 40.63, 36.76, 33.26, 30.09, 27.22, 24.63, 22.29, 20.16, 18.24, 16.51, 14.93, 13.51;
        100, 73.29, 53.71, 39.37, 28.85, 21.14, 15.50, 11.38, 8.32, 6.10, 4.47, 3.28, 2.40, 1.76, 1.29, 0.95, 0.69, 0.51, 0.37, 0.27, 0.20];

% Temperature in Kelvin
temps = [50, 65, 80] + 273;

fit_sets = [1 2 3];  % experiments used in the fit; the others are predicted

%% Setup
c_inits = conc(:,1)';  % initial concentrations of all experiments

% Initial guesses for both free parameters
k0_init = 3.1e11;  % 1/h
Ea_init = 75e3;    % J/mol

% k0 (~1e11) and Ea (~1e5) differ by six orders of magnitude, which stalls
% the optimizer. Fitting log10(k0) instead puts both on a similar scale.
p0 = [log10(k0_init), Ea_init];
lb = [log10(1),     5e3];    % k0 >= 1 1/h,    Ea >= 5 kJ/mol
ub = [log10(10e12), 150e3];  % k0 <= 1e13 1/h, Ea <= 150 kJ/mol

%% Optimization using lsqnonlin (least squares)
% All experiments are fitted at once by concatenating their residuals into
% a single vector; lsqnonlin minimizes the sum of its squares.
opts = optimoptions('lsqnonlin', 'Display', 'iter');
p_opt = lsqnonlin(@(p) resid_Temp_multi(p, time, c_inits(fit_sets), ...
                  conc(fit_sets,:), temps(fit_sets)), p0, lb, ub, opts);

k0_fit = 10^p_opt(1);  % back-transform from log10 scale
Ea_fit = p_opt(2);

fprintf('Ea = %.2f kJ/mol\n', Ea_fit/1000);
fprintf('A  = %.2e 1/h\n', k0_fit);

%% Simulate every experiment with the fitted parameters
sim_res = zeros(size(conc));
for i = 1:size(conc, 1)
    sim_res(i,:) = sim_exp_Temp(time, c_inits(i), k0_fit, Ea_fit, temps(i));
end

%% Plotting
figure('Position', [100, 100, 1400, 600]);
colors = lines(3);

% Left: concentration vs time. Solid lines are fits, dashed are predictions.
subplot(1, 2, 1);
hold on; grid on;
for i = 1:size(conc, 1)
    label_T = [num2str(temps(i)-273), ' °C'];
    if ismember(i, fit_sets)
        style = '-';  tag = ', fit';
    else
        style = '--'; tag = ', predicted';
    end
    plot(time, sim_res(i,:), style, 'Color', colors(i,:), 'LineWidth', 2, ...
        'DisplayName', [label_T, tag]);
    plot(time, conc(i,:), 'o', 'Color', colors(i,:), 'MarkerSize', 9, ...
        'DisplayName', [label_T, ', exp']);
end
xlabel('$t \quad \mathrm{(h)}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$c \quad \mathrm{(mol\,m^{-3})}$', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'best');
hold off;

% Right: parity plot (simulated vs experimental concentration)
subplot(1, 2, 2);
hold on; grid on;
for i = 1:size(conc, 1)
    plot(conc(i,:), sim_res(i,:), 'o', 'Color', colors(i,:), 'LineWidth', 2, ...
        'MarkerSize', 9, 'DisplayName', [num2str(temps(i)-273), ' °C']);
end
plot([0 100], [0 100], 'k:', 'HandleVisibility', 'off');
xlabel('$c_\mathrm{exp} \quad \mathrm{(mol\,m^{-3})}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$c_\mathrm{fit} \quad \mathrm{(mol\,m^{-3})}$', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'best');
hold off;

function dcdt = balance_Temp(~, c, k0, Ea, T)
    % Batch-reactor mass balance for a first-order Arrhenius reaction:
    % dc/dt = -k0 * exp(-Ea / RT) * c
    R = 8.314;  % universal gas constant, J/(mol K)
    dcdt = -k0 * exp(-Ea / (R * T)) * c;
end

function sim = sim_exp_Temp(time, c_init, k0, Ea, T)
    % Solves the ODE for one experiment and returns the simulated
    % concentrations at the experimental time points as a row vector.
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [~, sol] = ode45(@(t, c) balance_Temp(t, c, k0, Ea, T), time, c_init, opts);
    sim = sol';
end

function res = resid_Temp_multi(p, time, c_inits, data, T)
    % Residuals of all experiments, concatenated into one row vector.
    %   p(1) = log10(k0), p(2) = Ea
    %   data : experimental concentrations, one row per experiment
    k0 = 10^p(1);
    Ea = p(2);
    res = [];
    for i = 1:length(T)
        sim = sim_exp_Temp(time, c_inits(i), k0, Ea, T(i));
        res = [res, sim - data(i,:)];
    end
end
```

:::{tip} Check yourself
The global fit to all three runs gives $E\un{a} = 74.93\ \mathrm{kJ\,mol^{-1}}$ and
$A = 3.80 \times 10^{11}\ \mathrm{h^{-1}}$, close to the two-step linearized result in the
previous section ($74.5\ \mathrm{kJ\,mol^{-1}}$). Fitted to the $50$ and
$65\ \mathrm{^\circ C}$ runs only, it predicts the $80\ \mathrm{^\circ C}$ run almost exactly.
:::

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
electronic-structure calculations of the underlying elementary steps.
<!-- Restore the link ([](#ch-microscopic)) when chapter 8 is released. -->

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
  discrimination.
  <!-- Restore the links ([](#ch-microscopic), [](#ch-mechanisms)) when chapters 8 and 9 are released. -->

- For a reaction network with more reactions than independent ones, the rank $R_\nu$ of the
  stoichiometric matrix sets the number of key species that must be measured. Non-key species follow
  from $\Delta\vect{n}_2 = \mtrx{N}_{2,1}\,\mtrx{N}_{1,1}^{-1}\,\Delta\vect{n}_1$, [](#eq-nonkey).
