---
title: Chemical Thermodynamics
short_title: Thermodynamics
label: ch-thermodynamics
---

:::{admonition} Learning objectives
:class: tip
After completing this chapter, you should be able to:

- Derive the equilibrium constant $K = \exp(-\Delta\un{r} G^0/RT)$ from the chemical-potential
  criterion $\sum_i \nu_i \mu_i = 0$.
- Distinguish the dimensionless equilibrium constant $K$ from the dimensional special forms
  $K\un{p}$ and $K\un{c}$, and convert between them.
- Compute $\Delta\un{r} G^0$, $\Delta\un{r} H^0$, and $\Delta\un{r} S^0$ from tabulated formation
  properties.
- State Kirchhoff's law and apply the Ulich approximations to obtain $\Delta\un{r} H^0(T)$ when
  only data at the reference temperature are available.
- Use the van't Hoff equation to shift an equilibrium constant from one temperature to another,
  and read $\Delta\un{r} H^0$ from the slope of a $\ln K$ vs. $1/T$ plot.
:::

:::{admonition} Activity
:class: warning
Who knows about the iron-making process? Sketch out what you think happens inside a blast furnace,
and what the overall reaction is.
:::

The overall reaction inside the furnace is

$$
\ce{Fe2O3 + 3 CO -> 2 Fe + 3 CO2} ,
$$

where $\ce{Fe2O3}$ is hematite. At the end of the 19th century a large number of these blast
furnaces were constructed, and engineers had by then developed the measurement tools to record
concentrations at various points along the length of the reactor. They measured a great deal of
unreacted CO at the exhaust and concluded that the reactor needed to be improved so that more of
the ore could be converted.

The engineers' approach was to make the reactor longer, improving the contacting time between gas
and ore. They built furnaces taller and taller, some reaching more than 30 m. To their
disappointment, conversion did not improve.

Le Châtelier was the first to recognize that this reaction is limited by thermodynamic equilibrium,
so increasing the residence time cannot push the system past its equilibrium composition. A few
years later Boudouard characterized the side reaction that is also active under blast-furnace
conditions and now bears his name,

$$
\ce{2 CO <=> CO2 + C} .
$$

As noted earlier, thermodynamics tells us nothing about *how fast* a system approaches
equilibrium — that is the job of kinetics. But we still need to know *where* the equilibrium lies
before we can sensibly evaluate kinetics, because every reaction rate slows as the system
approaches its equilibrium composition. Thermodynamics therefore sits at the core of everything we
do in this course, kinetics included.

The takeaway is simple: **always look at the equilibrium first.** If a reactor already operates
close to its thermodynamic limit, no improvement in catalyst or reactor design will significantly
raise conversion. Equilibrium is the ultimate ceiling on reactor performance.

## Equilibrium constant

Let us continue with the ammonia synthesis,

$$
\ce{N2 + 3 H2 <=> 2 NH3} .
$$

We can use the law of mass action to write down a partial-pressure equilibrium constant

$$
K\un{p} = \frac{p_{\ce{NH3}}^2}{p_{\ce{N2}}\, p_{\ce{H2}}^3} ,
\qquad
[K\un{p}] = \mathrm{bar^{-2}} ,
$$

which carries units, since the exponents on the partial pressures do not balance. We will show
below that the *true* equilibrium constant $K$ is dimensionless, and that $K\un{p}$ is one of
several special forms of $K$.

:::{admonition} Discussion
:class: seealso
At equilibrium there is no net reaction, so by definition the left-hand side of the law of mass
action equals a constant. Where does that constant come from, and why does it depend on temperature
but not on the initial composition?
:::

Most of you have already learned **Le Châtelier's principle**: if a system at equilibrium is
disturbed by changing the conditions, the position of equilibrium shifts to counteract the change
and reestablish equilibrium. For the ammonia synthesis, raising the total pressure shifts the
reaction toward the product side, because the right-hand side has fewer moles of gas.

Some of you may also have seen the equilibrium constant written as

$$
K = \exp\left(\frac{-\Delta\un{r} G^0}{RT}\right) ,
$$ (eq-k-def)

where $\Delta\un{r} G^0$ is the Gibbs free energy change of the reaction at the standard pressure
$p^0 = 1\ \mathrm{bar}$. There is no IUPAC-defined *standard temperature*, but most thermodynamic
data are tabulated at the agreed-upon reference temperature
$T\un{R} = 298.15\ \mathrm{K}$.

At first glance it is not obvious how the partial-pressure form and the exponential form are
connected. The remainder of this section shows the link.

Consider a closed system holding a mixture — here $\ce{N2}$, $\ce{H2}$, and $\ce{NH3}$ — undergoing
a reversible change.

:::{figure} ../figures/Nh3-system.png
:label: fig-nh3-system
:alt: Sketch of two irregular blobs representing a closed system. The left blob is labelled "system" and contains N2, H2 and NH3; an arrow labelled "dG with dn_i at constant T, p" points to the right blob, which contains the same three species after a change in composition. A second label identifies the system's Gibbs free energy G.
:width: 70%

Example: $\ce{NH3}$ system undergoing a change in composition (at fixed $T$, $p$) due to a chemical
reaction.
:::

This system can be described by the Gibbs free energy, which is the maximum non-expansion work the
system can do. We are interested in how the system changes over the course of reaction, and we
typically assume this happens at constant temperature and pressure. In that case

$$
\mathrm{d}G = \sum_i \mu_i\, \mathrm{d}n_i ,
$$

where $\mu_i$ is the chemical potential of component $i$, a combination of the first and second
laws of thermodynamics. The chemical potential is the partial molar Gibbs free energy,

$$
\mu_i \equiv \left(\frac{\partial G}{\partial n_i}\right)_{T,\,p,\,n_{j \neq i}} .
$$

For a single reaction at constant temperature and pressure, the change in moles of each species is
fixed by stoichiometry through the extent of reaction $\xi$ introduced in
[](#ch-stoichiometry): $\mathrm{d}n_i = \nu_i\, \mathrm{d}\xi$. Hence

$$
\mathrm{d}G = \sum_i \left(\nu_i \mu_i\right) \mathrm{d}\xi ,
$$

and therefore

$$
\left(\frac{\partial G}{\partial \xi}\right)_{T,\,p} = \sum_i \left(\nu_i \mu_i\right)
= \Delta\un{r} G .
$$

At equilibrium the Gibbs energy is at a minimum, so its derivative with respect to the extent of
reaction must vanish:

$$
\sum_i \left(\nu_i \mu_i\right) = 0 = \Delta\un{r} G .
$$ (eq-equilibrium-criterion)

:::{figure} ../figures/Gibbs.png
:label: fig-gibbs-energy-extent
:alt: Hand-drawn plot of Gibbs energy on the vertical axis against extent of reaction on the horizontal axis for the example reaction A in equilibrium with B. The curve is convex with a single minimum. On the descending left branch the slope is negative, annotated Delta_r G < 0, forward reaction is spontaneous. On the ascending right branch the slope is positive, annotated Delta_r G > 0, reverse reaction is spontaneous. At the minimum the tangent is horizontal, annotated Delta_r G = 0, equilibrium.
:width: 80%

Gibbs free energy of reaction as a function of the extent of the reaction.
:::

At equilibrium, then, the sum of the chemical potentials weighted by their stoichiometric
coefficients must be zero. We do not use this equation as it stands, because the chemical potential
is expressed in terms of different quantities:

$$
\mu_i = \mu_i^0 + RT \ln a_i = \Delta\un{f} G_i^0 + RT \ln a_i ,
$$

where $\Delta\un{f} G_i^0$ is the Gibbs free energy of formation of species $i$ and $a_i$ is the
**activity**, essentially an effective concentration or pressure. The term $RT \ln a_i$ describes
the pressure dependence of the chemical potential. Gibbs free energies of formation are tabulated
values — see, for example, the NIST Chemistry WebBook or the Active Thermochemical Tables.

For gases, which are the majority in this course, the activity is defined as

$$
a_i = \frac{f_i}{p^0} = \varphi_i\, y_i\, \frac{p}{p^0} ,
$$

based on the **fugacity** $f_i$ of a species. The name derives from the Latin *fugare*, "to flee",
and reflects the tendency of a species to escape its own phase. Here $p^0 = 1\ \mathrm{bar}$ is the
standard pressure and $\varphi_i$ is the fugacity coefficient, which equals one for an ideal gas.
In most cases we will assume ideal-gas behavior, so that the activity is described by the partial
pressure.

:::{figure} ../figures/Fugacity.png
:label: fig-fugacity
:alt: Hand-drawn plot of chemical potential mu_i on the vertical axis against pressure p on the horizontal axis, with the vertical axis origin marked as the Gibbs free energy of formation. A solid curve labelled "perfect" rises and flattens with increasing pressure. A dashed curve lying below it in the low-pressure region is labelled "attractive interactions"; a second dashed curve lying above it at high pressure is labelled "repulsive interactions". The standard pressure p-nought is marked on the pressure axis.
:width: 80%

Chemical potential as a function of pressure for an ideal and a non-ideal gas.
:::

Multiplying the chemical-potential equation by $\nu_i$ and summing gives

$$
\begin{aligned}
\sum_i \left(\nu_i \mu_i\right) &= \sum_i \nu_i \Delta\un{f} G_i^0 + RT \sum_i \nu_i \ln a_i , \\
0 &= \Delta\un{r} G^0 + RT \ln \prod_i a_i^{\nu_i} ,
\end{aligned}
$$

with

$$
K = \prod_i a_i^{\nu_i} .
$$

If $K > 1$ the equilibrium lies on the product side; if $K < 1$ it lies on the side of the
reactants. Rearranging leads back to [](#eq-k-def),

$$
K = \exp\left(\frac{-\Delta\un{r} G^0}{RT}\right) .
$$

For an ideal gas the fugacity is simply the partial pressure, so

$$
K = \left(\frac{1}{p^0}\right)^{\sum_i \nu_i} \prod_i p_i^{\nu_i} .
$$

$K$ is the true, dimensionless equilibrium constant. In practice we often do not work with $K$
directly, but with defined *special* equilibrium constants such as

$$
K\un{p} = \prod_i p_i^{\nu_i} , \qquad K\un{c} = \prod_i c_i^{\nu_i} .
$$

These special equilibrium constants can carry units. They all describe the same physical
equilibrium, just expressed in different variables. The two forms are connected through the
ideal gas law $p_i = c_i RT$:

$$
K\un{c} = K\un{p} \left(RT\right)^{-\sum_i \nu_i} .
$$

::::{admonition} Example 3.1 — Three constants, one equilibrium
:class: note
:label: ex-three-constants
In [](#ch-stoichiometry) we used $K\un{p} = 1.397\ \mathrm{bar^{-2}}$ for ammonia synthesis at
450 K without saying where the units came from. Report the same equilibrium as the dimensionless
$K$ and as $K\un{c}$.

:::{dropdown} Solution
Everything follows from the sum of the stoichiometric coefficients. For
$\ce{N2 + 3 H2 <=> 2 NH3}$,

$$
\sum_i \nu_i = 2 - 1 - 3 = -2 ,
$$

which is exactly why $K\un{p}$ carries $\mathrm{bar^{-2}}$: the pressure exponents fail to cancel
by two.

**The dimensionless constant.** Rearranging $K = (1/p^0)^{\sum_i \nu_i} K\un{p}$,

$$
K = \left(p^0\right)^{2} K\un{p}
  = (1\ \mathrm{bar})^2 \times 1.397\ \mathrm{bar^{-2}}
  = 1.397 .
$$

The number is unchanged and the units are gone. This is a coincidence of choosing
$p^0 = 1\ \mathrm{bar}$ — numerically invisible, but conceptually the step that makes $K$ a proper
dimensionless argument for $\ln$ in $\Delta\un{r}G^0 = -RT\ln K$.

**The concentration constant.** With $\sum_i \nu_i = -2$ the conversion factor is $(RT)^{+2}$:

```matlab
Kp = 1.397;          % bar^-2
R  = 0.08314;        % bar L /(mol K)   <- units chosen to match Kp
T  = 450;            % K

K  = Kp * 1^2                  % dimensionless, p0 = 1 bar
Kc = Kp * (R*T)^2              % Kc = 1955 L^2/mol^2
```

so $K\un{c} = 1.96 \times 10^{3}\ \mathrm{L^2\,mol^{-2}}$.

The three numbers — 1.397, 1.397, and 1955 — describe one and the same physical equilibrium. Only
$K$ is dimensionless, and only $K$ belongs in $\exp(-\Delta\un{r}G^0/RT)$. Note also that $K\un{c}$
is the one that changes if you quote it in $\mathrm{m^6\,mol^{-2}}$ instead
($1.96\times10^{-3}$), which is a common source of factor-of-$10^6$ errors.
:::
::::

The Gibbs free energy of reaction is derived from the enthalpy of reaction and the entropy of
reaction. Alternatively it can be obtained directly from the Gibbs free energies of formation of
the individual species:

$$
\begin{aligned}
\Delta\un{r} G^0 &= \Delta\un{r} H^0 - T \Delta\un{r} S^0 , \\
\Delta\un{r} G^0 &= \sum_i \nu_i\, \Delta\un{f} G_i^0 , \\
\Delta\un{r} H^0 &= \sum_i \nu_i\, \Delta\un{f} H_i^0 , \\
\Delta\un{r} S^0 &= \sum_i \nu_i\, S_i^0 .
\end{aligned}
$$

In most data sources, Gibbs free energies of formation are not reported directly. Instead, tables
list enthalpies of formation at standard pressure and the reference temperature alongside standard
entropies, and $\Delta\un{r} G^0$ is reconstructed from the equations above. Determining accurate
enthalpies of formation is harder than it sounds.

::::{admonition} Example 3.2 — From a data table to an equilibrium constant
:class: note
:label: ex-drg-from-tables
Compute $\Delta\un{r} H^0$, $\Delta\un{r} S^0$, $\Delta\un{r} G^0$, and $K$ at
$T\un{R} = 298.15\ \mathrm{K}$ for $\ce{N2 + 3 H2 <=> 2 NH3}$, using only tabulated formation
properties.

| species | $\Delta\un{f} H^0$ (kJ mol⁻¹) | $S^0$ (J mol⁻¹ K⁻¹) |
|:---|---:|---:|
| $\ce{N2}$  | 0      | 191.61 |
| $\ce{H2}$  | 0      | 130.68 |
| $\ce{NH3}$ | −45.90 | 192.77 |

:::{dropdown} Solution
Both reaction properties are stoichiometry-weighted sums, with $\nu_i$ *signed*. The elements are
in their reference states, so their enthalpies of formation are zero by definition — but their
entropies are emphatically not.

$$
\begin{aligned}
\Delta\un{r} H^0 &= 2(-45.90) - (0) - 3(0) = -91.80\ \mathrm{kJ\,mol^{-1}} , \\
\Delta\un{r} S^0 &= 2(192.77) - 191.61 - 3(130.68) = -198.11\ \mathrm{J\,mol^{-1}\,K^{-1}} .
\end{aligned}
$$

Both are negative, and that is the whole story of ammonia synthesis in two numbers: the reaction
releases heat, but it converts four moles of gas into two, so it is strongly entropy-disfavored.
Which term wins depends on temperature.

At the reference temperature the enthalpy term still dominates:

$$
\Delta\un{r} G^0 = \Delta\un{r} H^0 - T\Delta\un{r} S^0
= -91\,800 - 298.15\,(-198.11) = -32.73\ \mathrm{kJ\,mol^{-1}} .
$$

```matlab
R  = 8.314;    TR = 298.15;                 % J/(mol K), K

nu   = [-1, -3,  2];                        % N2, H2, NH3
dfH  = [ 0,  0, -45.90e3];                  % J/mol
S0   = [191.61, 130.68, 192.77];            % J/(mol K)

drH = nu*dfH.';                             % -91.80 kJ/mol
drS = nu*S0.';                              % -198.11 J/(mol K)
drG = drH - TR*drS                          % -32.73 kJ/mol
K   = exp(-drG/(R*TR))                      % K = 5.43e5
```

giving $K(298.15\ \mathrm{K}) = 5.43 \times 10^{5}$.

Two checks worth making a habit. Tabulations that *do* list $\Delta\un{f} G^0$ give
$-16.4\ \mathrm{kJ\,mol^{-1}}$ for $\ce{NH3}$, so $\Delta\un{r} G^0 = 2(-16.4) =
-32.8\ \mathrm{kJ\,mol^{-1}}$ — agreement to better than a percent. And note how brutally the
exponential amplifies: an error of $5\ \mathrm{kJ\,mol^{-1}}$ in $\Delta\un{r} G^0$, well within
the spread between data sources, moves $K$ by a factor of $\exp(5000/RT) \approx 7.5$. This is why
the remark above — that determining accurate enthalpies of formation is harder than it sounds —
matters in practice.
:::
::::

## Temperature dependence of the equilibrium constant

Since enthalpies of formation, entropies, and Gibbs free energies of formation are tabulated only
at the reference temperature of 298.15 K, we need to be able to calculate the enthalpy of formation
at other temperatures. For enthalpies this conversion is straightforward using the heat capacity at
constant pressure, $c\un{p}$. As with the properties above, the heat capacity change of the
reaction is the stoichiometry-weighted sum over species,

$$
\Delta\un{r} c\un{p}(T) = \sum_{i=1}^{N\un{species}} \nu_i\, c\un{p,i} .
$$

Combined with the standard enthalpy of reaction $\Delta\un{r} H^0(T\un{R})$, the enthalpy of
reaction at any temperature follows from **Kirchhoff's law**:

$$
\Delta\un{r} H^0(T) = \Delta\un{r} H^0(T\un{R})
+ \int_{T\un{R}}^{T} \Delta\un{r} c\un{p}(T)\, \mathrm{d}T .
$$ (eq-kirchhoff)

Heat capacities are themselves temperature dependent. We will not go into their measurement or
derivation, which would take us too deep into physical chemistry. As engineers we simply rely on
data carefully tabulated by others, often supplied as polynomials in temperature.

The heat-capacity integral is usually much smaller in magnitude than
$\Delta\un{r} H^0(T\un{R})$ itself. Chemical engineers therefore introduce approximations to
$\Delta\un{r} c\un{p}(T)$ that make temperature-dependent enthalpies of reaction — and hence
equilibrium constants — easier to determine. The first commonly used approximation is

$$
\Delta\un{r} c\un{p}(T) = 0 \qquad \text{(first Ulich approximation)} ,
$$

in which case the enthalpy of reaction does not change with temperature:
$\Delta\un{r} H^0(T) = \Delta\un{r} H^0(T\un{R})$.

A slightly more accurate but still tractable approximation assumes a constant,
temperature-independent $\Delta\un{r} c\un{p}$,

$$
\Delta\un{r} c\un{p}(T) = \text{const} \qquad \text{(second Ulich approximation)} ,
$$

in which case Kirchhoff's law integrates trivially to

$$
\Delta\un{r} H^0(T) = \Delta\un{r} H^0(T\un{R})
+ \Delta\un{r} c\un{p}\left(T - T\un{R}\right) .
$$

The most accurate but most laborious option is to use the full temperature-dependent
$\Delta\un{r} c\un{p}(T)$, typically supplied as a polynomial fit. **In this course the first Ulich
approximation is used by default.** For the gas-phase reactions we encounter, the integral over
$\Delta\un{r} c\un{p}(T)$ is usually small compared with $\Delta\un{r} H^0(T\un{R})$ itself, so the
simpler form gives an answer that is good enough for reactor-design purposes.

Whenever we are interested in the temperature dependence of the enthalpy or the Gibbs free energy
of reaction, what we ultimately want is the position of the equilibrium and the equilibrium
constant. We have already defined the true equilibrium constant at the reference temperature,

$$
K(T\un{R}) = \exp\left(\frac{-\Delta\un{r} G^0(T\un{R})}{R\,T\un{R}}\right) ,
$$

but we typically need it at some other temperature $T$. The most accurate route is to compute
$\Delta\un{r} G^0(T)$ from $\Delta\un{r} H^0(T)$ and $\Delta\un{r} S^0(T)$ with their full
temperature dependencies, then substitute. That is exact, but tedious.

The easiest simplification is to invoke the first Ulich approximation,
$\Delta\un{r} c\un{p}(T) = 0$, so that $\Delta\un{r} H^0$ and $\Delta\un{r} S^0$ are both treated as
temperature-independent. Combining the Gibbs–Helmholtz relation with
$K = \exp(-\Delta\un{r} G^0/RT)$,

$$
\begin{aligned}
K(T) &= \exp\left(\frac{-\Delta\un{r} G^0(T)}{RT}\right) , \\
\ln K(T) &= \frac{-\Delta\un{r} G^0(T)}{RT}
          = \frac{-\left(\Delta\un{r} H^0 - T\Delta\un{r} S^0\right)}{RT} , \\
\ln K(T) &= -\frac{\Delta\un{r} H^0}{RT} + \frac{\Delta\un{r} S^0}{R} .
\end{aligned}
$$

Although the enthalpy of reaction has been frozen at its reference value, the Gibbs free energy of
reaction — and therefore $K$ — still depends on $T$ through the explicit $1/T$ factor.
Differentiating with respect to $T$, holding $\Delta\un{r} H^0$ and $\Delta\un{r} S^0$
temperature-independent, yields

$$
\frac{\mathrm{d}\ln K}{\mathrm{d}T} = \frac{\Delta\un{r} H^0}{RT^2} ,
$$ (eq-vant-hoff)

the **van't Hoff equation**, which describes the temperature dependence of the equilibrium
constant. Under our assumptions $K(T)$ depends only on the standard enthalpy of reaction, which is
widely tabulated, so the equation is extremely useful for converting an equilibrium constant from
one temperature to another with minimal thermodynamic input.

Integrating [](#eq-vant-hoff),

$$
\int_{K(T_1)}^{K(T_2)} \mathrm{d}\ln K
= \int_{T_1}^{T_2} \frac{\Delta\un{r} H^0}{RT^2}\, \mathrm{d}T ,
$$

gives

$$
\begin{aligned}
\ln\left(\frac{K(T_2)}{K(T_1)}\right)
  &= -\frac{\Delta\un{r} H^0}{R}\left(\frac{1}{T_2} - \frac{1}{T_1}\right) , \\
K(T_2) &= K(T_1) \exp\left[-\frac{\Delta\un{r} H^0}{R}
          \left(\frac{1}{T_2} - \frac{1}{T_1}\right)\right] .
\end{aligned}
$$ (eq-vant-hoff-integrated)

The van't Hoff equation is also linear when plotted as $\ln K$ against $1/T$, with slope
$-\Delta\un{r} H^0/R$.

:::{figure} ../figures/Vanthoff.png
:label: fig-vant-hoff
:alt: Two hand-drawn plots side by side, each with ln K on the vertical axis and 1/T on the horizontal axis. The left plot, labelled endothermic, shows a straight line with negative slope. The right plot, labelled exothermic, shows a straight line with positive slope. On both, a small slope triangle marks the gradient, annotated minus Delta H over R.
:width: 100%

Van't Hoff diagram for an endothermic and an exothermic reaction.
:::

::::{admonition} Example 3.3 — Shifting $K$ to reactor temperature, and what it costs
:class: note
Ammonia synthesis runs nowhere near 298 K. Starting from $K(298.15\ \mathrm{K}) = 5.43\times10^5$
found in [Example 3.2](#ex-drg-from-tables), use the van't Hoff equation to obtain $K$ at 450 K,
and compare with the $K\un{p} = 1.397\ \mathrm{bar^{-2}}$ that [](#ch-stoichiometry) asked you to
take on trust.

:::{dropdown} Solution
Apply [](#eq-vant-hoff-integrated) with $\Delta\un{r} H^0 = -91.80\ \mathrm{kJ\,mol^{-1}}$ held
constant — the first Ulich approximation, which is this course's default:

```matlab
R = 8.314;    T1 = 298.15;    T2 = 450;      % J/(mol K), K, K

K1  = 5.432e5;                                % from Example 3.2
drH = -91.80e3;                               % J/mol, first Ulich: constant

K2 = K1*exp(-drH/R*(1/T2 - 1/T1))             % K2 = 2.03
```

$K$ collapses from $5.4\times10^5$ to about **2.03** over 152 K. That is the exponential in
[](#eq-vant-hoff) doing its work, and it is the central engineering tension of the Haber–Bosch
process: the reaction is exothermic, so heating it to get a workable *rate* destroys the
equilibrium *yield*. Everything about the industrial process — the pressure, the recycle loop, the
catalyst — is a response to that conflict.

**Now the comparison.** Chapter 2 quoted $K\un{p} = 1.397\ \mathrm{bar^{-2}}$, equivalently
$K = 1.397$ from [Example 3.1](#ex-three-constants). We obtained 2.03 — high by 45 %. The
approximation is what failed, not the arithmetic. Recall $\Delta\un{r} c\un{p} \approx
-45.5\ \mathrm{J\,mol^{-1}\,K^{-1}}$ for this reaction, which is not negligible over 152 K. Redoing
it with the second Ulich approximation, $\Delta\un{r}c\un{p} = \text{const}$, so that

$$
\Delta\un{r} H^0(T_2) = \Delta\un{r} H^0(T\un{R}) + \Delta\un{r} c\un{p}(T_2 - T\un{R}) ,
\qquad
\Delta\un{r} S^0(T_2) = \Delta\un{r} S^0(T\un{R}) + \Delta\un{r} c\un{p} \ln\frac{T_2}{T\un{R}} ,
$$

gives $\Delta\un{r}G^0(450) = -1.13\ \mathrm{kJ\,mol^{-1}}$ and $K = 1.35$ — now within 4 % of the
quoted value. The residual gap is the temperature dependence of $\Delta\un{r} c\un{p}$ itself,
which the full polynomial treatment captures and which yields the 1.397 used in chapter 2.

So the three levels of approximation introduced above land at

| treatment | $K(450\ \mathrm{K})$ | error |
|:---|---:|---:|
| first Ulich, $\Delta\un{r}c\un{p} = 0$ | 2.03 | +45 % |
| second Ulich, $\Delta\un{r}c\un{p} = \text{const}$ | 1.35 | −3 % |
| full $\Delta\un{r}c\un{p}(T)$ polynomial | 1.397 | — |

The lesson is not that the first Ulich approximation is bad — it got the collapse from $10^5$ to
order unity, which is the engineering conclusion — but that you should know its accuracy before
quoting a conversion to three digits.
:::
::::

## Summary

- At constant $T$ and $p$, equilibrium is fixed by $\sum_i \nu_i \mu_i = 0$. Substituting
  $\mu_i = \mu_i^0 + RT\ln a_i$ yields the dimensionless equilibrium constant
  $K = \prod_i a_i^{\nu_i} = \exp(-\Delta\un{r} G^0/RT)$, [](#eq-k-def).
- For ideal gases $a_i = p_i/p^0$, and $K$ relates to the special forms
  $K\un{p} = \prod_i p_i^{\nu_i}$ and $K\un{c} = \prod_i c_i^{\nu_i}$ via
  $K = (1/p^0)^{\sum_i \nu_i} K\un{p}$ and $K\un{c} = K\un{p}(RT)^{-\sum_i \nu_i}$. The special
  forms can carry units; $K$ itself is always dimensionless.
- $\Delta\un{r} G^0$, $\Delta\un{r} H^0$, and $\Delta\un{r} S^0$ are reconstructed from tabulated
  formation properties as stoichiometry-weighted sums. Only $\Delta\un{f} H^0$ and $S^0$ are
  typically reported, and $\Delta\un{r} G^0 = \Delta\un{r} H^0 - T\Delta\un{r} S^0$ closes the loop.
- Kirchhoff's law, [](#eq-kirchhoff), shifts $\Delta\un{r} H^0$ from the reference temperature to
  any $T$ using $\Delta\un{r} c\un{p}(T)$. The first Ulich approximation
  ($\Delta\un{r} c\un{p} = 0$) is the default in this course; the second
  ($\Delta\un{r} c\un{p} = \text{const}$) is a refinement; the full polynomial form is used only
  when high accuracy is required.
- Under the first Ulich approximation, the van't Hoff equation
  $\mathrm{d}\ln K/\mathrm{d}T = \Delta\un{r} H^0/(RT^2)$, [](#eq-vant-hoff), integrates to
  [](#eq-vant-hoff-integrated), allowing $K$ to be transferred between two temperatures using only
  $\Delta\un{r} H^0$. A plot of $\ln K$ vs. $1/T$ has slope $-\Delta\un{r} H^0/R$.
