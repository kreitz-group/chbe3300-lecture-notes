---
title: Heterogeneous Catalysis
short_title: Heterogeneous catalysis
label: ch-hetcat
---

<!-- LaTeX source: HeterogeneousCatalysis.tex -->
<!-- Porting notes: mhchem is math-only here, so \ce{} in prose needs $...$; and \un{} expands to
     _{\textrm{...}}, so species subscripts must be written _{\ce{CO}}, never \un{\ce{CO}}.
     Surface sites are \ce{\ast}, adsorbates \ce{A\ast}, adjacent pairs \ce{\ast\ast}. Distinct
     site types are written \ce{\ast_1} / \ce{\ast_2} (subscripted) rather than the source's
     "\ast 1", which mhchem would read as a stoichiometric coefficient.
     Nested directives need a longer outer fence (:::: around :::). -->

:::{admonition} Learning objectives
:class: tip
After completing this chapter, you should be able to:

- List the seven steps of heterogeneous catalysis and identify which are transport steps and which
  are chemical steps.
- Distinguish physisorption from chemisorption by binding strength and by the effect on the
  adsorbate, and recognize dissociative adsorption.
- State the Langmuir assumptions and derive the Langmuir isotherm for molecular, competitive, and
  dissociative adsorption.
- Predict how the equilibrium coverage responds to temperature and to the adsorption enthalpy.
- Use chemisorption to determine the active-site count, metallic surface area, mean particle size,
  and dispersion of a supported catalyst.
- Distinguish Langmuir–Hinshelwood from Eley–Rideal steps and write the corresponding rate
  expressions in surface concentrations or coverages.
- Derive an LHHW rate law from a catalytic mechanism by combining a rate-determining step,
  quasi-equilibrated steps, and a site balance.
- Extract apparent reaction orders and apparent activation energies from a derived rate law, and
  explain the rate maximum in $\ce{CO}$ oxidation.
- Use the most abundant reaction intermediate (MARI) assumption to simplify a site balance.
- Compute the Thiele modulus and effectiveness factor, and explain how pore diffusion lowers the
  apparent activation energy.
:::

So far we have looked only at homogeneous reactions in the gas or liquid phase, which often proceed
through radical intermediates. We now shift our attention to catalytic reactions.

:::{admonition} Discussion
:class: seealso
What is a catalyst? What properties make a material a catalyst?
:::

A catalyst accelerates a reaction by providing an alternative reaction pathway with a lower
activation barrier, and is regenerated at the end of each catalytic cycle. Compared with
uncatalyzed reactions, catalysts can boost rates by many orders of magnitude and can also steer
selectivity toward a desired product. Compared with homogeneous catalysis, heterogeneous catalysis
has the practical advantage that the catalyst is in a different phase from the reactants and
products, which makes separating it from the reaction mixture far easier.

Our main focus is on understanding the kinetics of catalytic reactions: we aim to develop catalytic
reaction mechanisms and then derive simplified rate laws from them. To do so, we first have to
understand the seven steps of catalysis.

Let us assume that our catalyst is a spherical, porous particle surrounded by a gas phase. A real
catalyst is typically a few weight percent of small metal or metal-oxide nanoparticles dispersed
over a high-surface-area porous support, rather than a solid block of metal — this is what makes
the surface-area-per-gram and the active-site density so high.

The following steps occur:

1. The reactants $\ce{A1}$ and $\ce{A2}$ diffuse from the bulk of the fluid to the outer surface of
   the catalyst. This is called **film diffusion**, because there is a stagnant film around the
   catalyst particles and we use the film model to describe the process.
2. The reactants continue and diffuse from the outer surface to the inner surface of the catalyst,
   where the metal nanoparticles are finely dispersed inside the pores. This is **pore diffusion**.
   Pores are typically 1–100 nm across. If the pores are really small, below about 5 nm, we also
   have to consider Knudsen diffusion, where the molecules collide more often with the wall than
   with each other.
3. The reactants **adsorb** on the catalyst surface, on the so-called active site.
4. The adsorbed intermediates **react** to form the adsorbed products.
5. The adsorbed products **desorb** from the catalyst surface.
6. The products diffuse through the pores back to the outer catalyst surface.
7. The products diffuse through the film back to the bulk of the fluid phase.

We will mostly look at steps 3–5, and only return to the diffusion steps at the end of the chapter.
The magic of catalytic reactions happens on the catalyst surface, and to explain that magic we have
to understand adsorption and desorption.

:::{figure} ../figures/catalysis.png
:label: fig-catalysis-7steps
:alt: Cross-section diagram. On the left is the gas phase, separated by a dashed line marking the film diffusion layer from a hatched region on the right representing the catalyst, into which a narrow pore extends. Arrows numbered 1 and 2 run from the gas through the film and into the pore; the label 3,4,5 sits at the pore where adsorption, surface reaction and desorption occur; arrows numbered 6 and 7 run back out through the pore and the film to the gas.
:width: 55%

Schematic illustration of a porous catalyst particle, highlighting the seven steps of heterogeneous
catalysis. Adapted from Hagen (2015).
:::

Reaction rates for heterogeneously catalyzed reactions are usually expressed with respect to the
surface area or the mass of the catalyst,

$$
r = \frac{1}{A}\frac{\mathrm{d}\xi}{\mathrm{d}t} .
$$

## Adsorption

<!-- source: HeterogeneousCatalysis.tex L71 -->

The first step in the catalytic cycle is the adsorption of the reactants on the catalyst surface.
Gas molecules interact with solid surfaces through collisions of the gas with the surface, and we
can quantify the collision frequency using kinetic gas theory, much as in collision theory. The
collision flux $Z\un{W}$ for a specific molecule is

$$
Z\un{W} = \frac{p}{\left(2\pi m k\un{B}T\right)^{1/2}} ,
$$ (eq-collision-flux)

where $m$ is the mass of the molecule. For $\ce{N2}$ at room temperature and 1 bar, $Z\un{W}$ is
around $3\times10^{27}\ \mathrm{m^{-2}\,s^{-1}}$. On a metal surface area of 1 m² we have around
$10^{19}$ surface atoms, give or take, so each metal atom is struck by roughly $10^{8}$ gas-phase
molecules per second.

Adsorption describes the process by which an atom or molecule binds to the catalyst surface. The
adsorbing molecule or atom is the **adsorbate**, and the surface is the **adsorbent**. Absorption,
by contrast, refers to a species being taken up into the bulk of the adsorbent — a different
phenomenon.

There are two types of adsorption: physisorption and chemisorption.

:::{admonition} Discussion
:class: seealso
Who knows what the difference is?
:::

Before explaining the two, let us define our adsorbent, the catalyst surface. Catalysts are
typically metal or metal-oxide surfaces, and a catalyst surface has many different active sites. We
usually assume that each active site can adsorb one adsorbate. We denote a surface site with an
asterisk, $\ast$, and the catalyst has a fixed total number of active sites.

**Physisorption** refers to a weak interaction of an adsorbate with a surface, through mostly
electrostatic or van der Waals forces. These interactions do not substantially change the electronic
configuration of the adsorbate, so physisorption is typically not strong enough to activate the
adsorbate toward a chemical reaction:

$$
\begin{aligned}
\ce{Ar + \ast &<=> Ar\ast} \\
\ce{H2 + \ast &<=> H2\ast}
\end{aligned}
$$

**Chemisorption**, in contrast, leads to very strong interactions between adsorbate and catalyst,
and to the formation of strong bonds with the surface, on the order of covalent bond strengths. CO,
for example, adsorbs very strongly on metal surfaces such as Pt or Ni. Chemisorption can also lead
to dissociation of the adsorbate, which we call **dissociative adsorption**:

$$
\begin{aligned}
\ce{CO + \ast &<=> CO\ast} \\
\ce{H2 + 2 \ast &<=> 2 H\ast}
\end{aligned}
$$

Physisorption and chemisorption can be distinguished by their binding strength. As a rule of thumb,
physisorbed molecules have adsorption enthalpies in the range $-5$ to
$-40\ \mathrm{kJ\,mol^{-1}}$, while chemisorbed species typically have $-40$ to
$-500\ \mathrm{kJ\,mol^{-1}}$; the boundary between the two regimes is not sharp. Dissociative
adsorption steps typically have an activation barrier.

:::{figure} ../figures/adsorption.png
:label: fig-pes-h2-diss
:alt: Potential-energy diagram against distance from a hatched Ni surface on the left. Curve 1, for molecular H2, comes in from the right and dips into a shallow well of depth delta-H-P just below zero. Curve 2, for two H atoms, starts at the H plus H level a distance E_D above the H2 level and falls steeply into a much deeper well of depth delta-H-C at point B, close to the surface. The two curves cross at point A slightly above zero; the height of that crossing above the H2 level is marked E_A, the activation barrier for dissociation.
:width: 60%

Lennard–Jones-style potential-energy diagram for the dissociative adsorption of $\ce{H2}$ on a metal
surface. The shallow well at long distance corresponds to molecular, physisorbed $\ce{H2}$; the
deeper well at short distance corresponds to two atomically chemisorbed H atoms. The crossing of the
two curves sets the activation barrier for dissociation.
:::

### Langmuir isotherm

<!-- source: HeterogeneousCatalysis.tex L131 -->

We are interested in the kinetics of this adsorption process, starting with molecular adsorption.
First, some assumptions about our catalyst and adsorbate. Assume the catalyst surface looks like a
checkerboard, where every square is an active site, with a finite total number of available surface
sites $N\un{available}$. During adsorption the surface fills with adsorbates, occupying
$N\un{occupied}$ sites.

We often express the concentration of occupied or vacant sites as a **fractional coverage**,

$$
\theta = \frac{\text{number of occupied sites}}{\text{number of available surface sites}}
= \frac{N\un{occupied}}{N\un{available}} .
$$ (eq-coverage-def)

Here $\theta = 1$ means the surface is completely occupied, $\theta = 0$ that it is empty, and
$\theta = 0.5$ that it is 50 % occupied and 50 % vacant. In the literature you will also see the
coverage defined by volume, since we work with gases,

$$
\theta = \frac{V}{V_\infty} ,
$$

where $V$ is the volume of adsorbed gas and $V_\infty$ the volume corresponding to a full
**monolayer** — all adsorption sites occupied, nothing else fits.

For the derivation of the adsorption equilibrium we need some assumptions about the nature of the
active sites. These were first introduced by Irving Langmuir:

- Adsorption cannot proceed beyond a monolayer coverage.
- All adsorption sites are equivalent.
- A molecule can only adsorb at a vacant site.
- There are no interactions between neighboring sites.

## Molecular adsorption

<!-- source: HeterogeneousCatalysis.tex L162 -->

Consider the reaction

$$
\ce{A + \ast ->[$k\un{ads}$][$k\un{des}$] A\ast} ,
$$

which could be the adsorption of $\ce{CO}$ on a Pt surface. We can write the rate of adsorption as

$$
r\un{ads} = \frac{\mathrm{d}N_{\ce{A\ast}}}{\mathrm{d}t} = k\un{ads}\, p\un{A}\, N_{\ce{\ast}} ,
$$

using the total number of moles, so the units of $k\un{ads}$ differ from those for non-catalyzed
reactions. The catalyst has a finite number of available surface sites $N\un{\ast,0}$, which we
typically do not know and which is hard to calculate.

Dividing both sides by the total number of available sites gives

$$
r\un{ads} = \frac{\mathrm{d}\theta_{\ce{A\ast}}}{\mathrm{d}t}
= k\un{ads}\, p\un{A}\, \theta_{\ce{\ast}} ,
$$

where $\theta_{\ce{A\ast}}$ is the coverage of species A and $\theta_{\ce{\ast}}$ the fraction of
vacant sites. With $p\un{A}$ in pressure units, $k\un{ads}$ now has units of bar$^{-1}$ or
Pa$^{-1}$. In catalysis it is very common to write the mass balance in terms of coverages.

For the desorption reaction

$$
\ce{A\ast ->[$k\un{des}$] A + \ast}
$$

we can likewise write

$$
r\un{des} = k\un{des}\, \theta_{\ce{A\ast}} .
$$

For the overall adsorption/desorption equilibrium,

$$
\ce{A + \ast} \eqa{k\un{ads}}{k\un{des}} \ce{A\ast} ,
$$

the overall change in the coverage of $\ce{A\ast}$ is

$$
\frac{\mathrm{d}\theta_{\ce{A\ast}}}{\mathrm{d}t} = r\un{ads} - r\un{des}
= k\un{ads}\, p\un{A}\, \theta_{\ce{\ast}} - k\un{des}\, \theta_{\ce{A\ast}} ,
$$

and at equilibrium

$$
\begin{aligned}
\frac{\mathrm{d}\theta_{\ce{A\ast}}}{\mathrm{d}t} &= 0 \\
k\un{ads}\, p\un{A}\, \theta_{\ce{\ast}} &= k\un{des}\, \theta_{\ce{A\ast}} \\
K\un{ads} = \frac{k\un{ads}}{k\un{des}}
&= \frac{\theta_{\ce{A\ast}}}{p\un{A}\, \theta_{\ce{\ast}}} .
\end{aligned}
$$

Most adsorption/desorption reactions are equilibrated, since they proceed on much faster timescales
than the surface reactions. Our goal is an expression for the coverage of the adsorbate in terms of
pressure and equilibrium constant, and for that we make a **site balance**: the fractional coverages
of occupied and vacant sites must sum to one,

$$
\theta_{\ce{\ast}} + \theta_{\ce{A\ast}} = 1
\quad \Rightarrow \quad
\theta_{\ce{\ast}} = 1 - \theta_{\ce{A\ast}} .
$$

Substituting,

$$
\begin{aligned}
K\un{ads} &= \frac{\theta_{\ce{A\ast}}}{p\un{A} - p\un{A}\theta_{\ce{A\ast}}} \\
\theta_{\ce{A\ast}} &= \frac{K\un{ads}\, p\un{A}}{1 + K\un{ads}\, p\un{A}} ,
\end{aligned}
$$ (eq-langmuir-isotherm)

which is the **Langmuir adsorption isotherm**.

:::{figure} ../figures/monolayer.png
:label: fig-langmuir-iso
:alt: Hand-drawn plot of the coverage of A against the partial pressure of A. The curve rises steeply and linearly from the origin, then bends over and flattens, approaching a horizontal dashed line labelled monolayer coverage without crossing it.
:width: 55%

Langmuir adsorption isotherm: equilibrium coverage $\theta_{\ce{A\ast}}$ as a function of the partial
pressure $p\un{A}$ for molecular adsorption. The coverage is linear at low pressure and approaches
monolayer saturation, $\theta = 1$, at high pressure.
:::

At low partial pressures of A, $K\un{ads}p\un{A} \ll 1$ and thus
$\theta_{\ce{A\ast}} \to K\un{ads}p\un{A}$. At high pressures, $K\un{ads}p\un{A} \gg 1$ and the
coverage approaches 1.

:::{admonition} Discussion
:class: seealso
What is an isotherm?
:::

It is called an isotherm because we determine the equilibrium coverage of $\ce{A\ast}$ at fixed
temperature — which sets $K\un{ads}$ — as a function of pressure.

The balance can equally be written in concentrations instead of coverages, which you will also find
frequently in the literature:

$$
\frac{[\ce{A\ast}]}{[\ast]_0} = \frac{K\un{ads}\, p\un{A}}{1 + K\un{ads}\, p\un{A}} ,
\qquad
[\ce{A\ast}] = \frac{K\un{ads}\, p\un{A}\,[\ast]_0}{1 + K\un{ads}\, p\un{A}} ,
$$

where $[\ce{A\ast}]$ is the concentration of adsorbed A and $[\ast]$ the concentration of vacant
sites. These concentrations are based on the surface, with units of mol m$^{-2}$. The catalyst has a
total concentration of available surface sites $[\ast]_0$, also often called the surface site
density $\Gamma$. Similarly, the rates are

$$
r\un{ads} = k\un{ads}\, p\un{A}\,[\ce{\ast}] , \qquad
r\un{des} = k\un{des}\,[\ce{A\ast}] .
$$

:::{admonition} Discussion
:class: seealso
What happens if I increase the temperature at fixed pressure — does the coverage increase or
decrease?
:::

We can answer this by replacing the equilibrium constant with its thermodynamic definition:

$$
\begin{aligned}
\theta_{\ce{A\ast}}
&= \frac{\exp\left(\frac{-\Delta\un{ads}G}{RT}\right) p\un{A}}
{1 + \exp\left(\frac{-\Delta\un{ads}G}{RT}\right) p\un{A}} \\
\theta_{\ce{A\ast}}
&= \frac{\exp\left(\frac{-\Delta\un{ads}H + T\Delta\un{ads}S}{RT}\right) p\un{A}}
{1 + \exp\left(\frac{-\Delta\un{ads}H + T\Delta\un{ads}S}{RT}\right) p\un{A}} .
\end{aligned}
$$

The adsorption enthalpy is usually exothermic, $\Delta\un{ads}H < 0$, so $-\Delta\un{ads}H$ in the
exponent is positive and the equilibrium constant decreases with increasing $T$. The coverage
therefore **decreases** with temperature, which also makes intuitive sense: higher temperature means
more thermal energy and more desorption. Conversely, making $\Delta\un{ads}H$ more exothermic — by
switching to a more strongly binding species, say — increases the coverage at fixed conditions.

## Competitive adsorption

<!-- source: HeterogeneousCatalysis.tex L266 -->

In most cases we cannot do much chemistry with just one adsorbate; we typically have two reactants.
In ammonia synthesis, for example, both $\ce{H2}$ and $\ce{N2}$ have to adsorb, and the two
adsorbates compete for adsorption sites. Consider a toy system of adsorbates A and B, assuming each
surface site can be occupied by either $\ce{A\ast}$ or $\ce{B\ast}$ but not both:

$$
\begin{aligned}
\ce{A + \ast &<=> A\ast} \\
\ce{B + \ast &<=> B\ast}
\end{aligned}
$$

Assuming equilibrium,

$$
K\un{1} = \frac{\theta_{\ce{A\ast}}}{p\un{A}\theta_{\ce{\ast}}} , \qquad
K\un{2} = \frac{\theta_{\ce{B\ast}}}{p\un{B}\theta_{\ce{\ast}}} .
$$

The **site balance** is now

$$
\begin{aligned}
1 &= \theta_{\ce{A\ast}} + \theta_{\ce{B\ast}} + \theta_{\ce{\ast}} \\
1 &= K\un{1}p\un{A}\theta_{\ce{\ast}} + K\un{2}p\un{B}\theta_{\ce{\ast}} + \theta_{\ce{\ast}} \\
1 &= \theta_{\ce{\ast}}\left(1 + K\un{1}p\un{A} + K\un{2}p\un{B}\right) \\
\theta_{\ce{\ast}} &= \frac{1}{1 + K\un{1}p\un{A} + K\un{2}p\un{B}} .
\end{aligned}
$$

To get the coverages of $\ce{A\ast}$ and $\ce{B\ast}$,

$$
\theta_{\ce{A\ast}} = K\un{1}p\un{A}\theta_{\ce{\ast}}
= \frac{K\un{1}p\un{A}}{1 + K\un{1}p\un{A} + K\un{2}p\un{B}} ,
\qquad
\theta_{\ce{B\ast}} = \frac{K\un{2}p\un{B}}{1 + K\un{1}p\un{A} + K\un{2}p\un{B}} .
$$

This generalizes naturally to $N$ adsorbates,

$$
\theta_{X_i} = \frac{K_i\, p_{X_i}}{1 + \sum_{j=1}^{N} K_j\, p_{X_j}} .
$$ (eq-competitive-isotherm)

We made several idealizing assumptions initially. On real catalysts the active sites are
heterogeneous, with different binding strengths; adsorbates can interact with each other, often
destabilizing one another; and it is possible to put more adsorbates on a surface than a single
monolayer. Many adsorption-isotherm models have been developed to account for these non-idealities.
The most prominent is the **Brunauer–Emmett–Teller (BET)** isotherm, which is *the* standard method
for determining the surface area of catalyst particles, via physisorption of $\ce{N2}$ at the
boiling point of liquid $\ce{N2}$. Other common types are the Temkin and Freundlich isotherms, but
we will not focus on those — the Langmuir isotherm is the standard approach and typically works well
enough.

## Dissociative adsorption

<!-- source: HeterogeneousCatalysis.tex L318 -->

Most molecules dissociate upon adsorption — $\ce{H2}$, for instance, because it leads to much more
strongly bound adsorbates. It is worth pausing on how remarkable that is: $\ce{H2}$ is very stable,
with a gas-phase dissociation energy of $436\ \mathrm{kJ\,mol^{-1}}$, yet it decomposes readily at
50 °C over most transition metals.

Dissociative adsorption combines adsorption and reaction into a single elementary step,

$$
\ce{A2 + 2 \ast <=> 2 A\ast} .
$$

For dissociative adsorption to occur, however, we need two *adjacent* surface sites, not just any
two sites. We can rewrite the reaction to reflect that pair,

$$
\ce{A2 + \ast\ast <=> A\ast\ast A} .
$$

We need an expression for the density of adjacent vacant pairs and adjacent occupied pairs. This
time we use concentrations. Lattice statistics gives

$$
[\ast\ast] = \frac{z}{2}\frac{[\ast]^2}{[\ast]_0} , \qquad
[\ce{A}\ast\ast\ce{A}] = \frac{z}{2}\frac{[\ce{A}\ast]^2}{[\ast]_0} ,
$$

where $z$ is the coordination number of the site and the factor of $1/2$ avoids overcounting
indistinguishable neighboring sites.

Consider $[\ast\ast]$, the number density of adjacent pairs of vacant sites. The probability of
randomly selecting a vacant site on the lattice is $[\ast]/[\ast]_0$, so the expected number of such
sites across the surface is $[\ast]_0 \cdot ([\ast]/[\ast]_0)$. After choosing one vacant site, we
consider the probability that an adjacent site is also vacant. That probability is again
$[\ast]/[\ast]_0$, but now the relevant pool is limited to the $z$ adjacent sites, giving
$z \cdot ([\ast]/[\ast]_0)$. Assuming these two events are independent, we multiply the
probabilities to obtain $[\ast\ast] = z \cdot ([\ast]^2/[\ast]_0)$, and a final factor of $1/2$
avoids double-counting indistinguishable pairs. The derivation for $[\ce{A}\ast\ast\ce{A}]$ is
similar.

The rates of adsorption and desorption are

$$
\begin{aligned}
r\un{ads} &= k\un{ads}\, p_{\ce{A2}}[\ast\ast]
= k\un{ads}\, p_{\ce{A2}}\frac{z}{2}\frac{[\ast]^2}{[\ast]_0} \\
r\un{des} &= k\un{des}\,[\ce{A}\ast\ast\ce{A}]
= k\un{des}\frac{z}{2}\frac{[\ce{A}\ast]^2}{[\ast]_0} ,
\end{aligned}
$$

and again at equilibrium,

$$
\begin{aligned}
r\un{ads} &= r\un{des} \\
k\un{ads}\, p_{\ce{A2}}\frac{z}{2}\frac{[\ast]^2}{[\ast]_0}
&= k\un{des}\frac{z}{2}\frac{[\ce{A}\ast]^2}{[\ast]_0} \\
K\un{ads} &= \frac{[\ce{A}\ast]^2}{p_{\ce{A2}}[\ast]^2} .
\end{aligned}
$$

The statistical factors cancel in the equilibrium constant, so it would not have mattered whether we
considered the distribution of sites. It would, however, have led to adsorption and desorption
*rates* that are overestimated.

To derive the coverage of A as a function of pressure, we solve for the unknown concentration of
vacant sites:

$$
\begin{aligned}
[\ast]_0 &= [\ast] + [\mathrm{A}\ast] \quad \rightarrow \quad
[\ast] = [\ast]_0 - [\mathrm{A}\ast] \qquad \text{(site balance)} \\
k\un{ads}\, p_{\ce{A2}}\left([\ast]_0 - [\mathrm{A}\ast]\right)^2
&= k\un{des}[\mathrm{A}\ast]^2 \\
\sqrt{K\un{ads}\, p_{\ce{A2}}}\left([\ast]_0 - [\mathrm{A}\ast]\right) &= [\mathrm{A}\ast] \\
[\ast]_0 &= [\mathrm{A}\ast]\left(\frac{1}{\sqrt{K\un{ads}\, p_{\ce{A2}}}} + 1\right) \\
\theta_{\ce{A\ast}} &= \frac{\sqrt{K\un{ads}\, p_{\ce{A2}}}}{1 + \sqrt{K\un{ads}\, p_{\ce{A2}}}} .
\end{aligned}
$$ (eq-dissociative-isotherm)

At low pressures the dissociative-adsorption coverage rises more rapidly with $p$ than the
molecular-adsorption coverage — square root versus linear — while at higher pressures molecular
adsorption has the stronger pressure dependence.

Typically the factor $z/2$ is simply absorbed into the rate constant:

$$
r\un{ads} = k'\un{ads}\, p_{\ce{A2}}\frac{[\ast]^2}{[\ast]_0} , \qquad
r\un{des} = k'\un{des}\frac{[\ce{A}\ast]^2}{[\ast]_0} .
$$

## Chemisorption as a catalyst characterization technique

<!-- source: HeterogeneousCatalysis.tex L386 -->

We defined the coverage as $\theta = N\un{occupied}/N\un{available}$, but typically we do not know
the number of available active sites.

:::{admonition} Discussion
:class: seealso
How could we determine that?
:::

We can use chemisorption itself for the measurement. You keep track of how much gas you have exposed
to your catalyst, either by a static volumetric method or by a pulse-chemisorption technique at
fixed pressure. If you know the total amount of, say, $\ce{H2}$ adsorbed, you can convert it to an
active-site count by assuming one H atom occupies one active site — usually true, though other
molecules can have different occupation stoichiometries. The amount of available sites is then

$$
n\un{available} = 2q_{\ce{H2}} , \qquad q_{\ce{H2}}\ [=]\ \mathrm{mol\,g^{-1}} ,
$$ (eq-site-count)

where the factor of 2 comes from the dissociative adsorption of $\ce{H2}$. We work on a per-gram
basis to make the quantities intensive. Another common probe molecule is $\ce{CO}$, but its
adsorption stoichiometry is less well defined and depends on the metal.

We can also use this site count to calculate the **metallic surface area**. The area occupied per
metal atom is known — for Ni, $\sigma\un{Ni} = 6.51\ \text{Å}^2$ — so

$$
a\un{metal} = 2q_{\ce{H2}}\, N\un{A}\, \sigma\un{metal} .
$$ (eq-metal-area)

With assumptions about particle shape, typically spherical, we can also estimate an average particle
diameter. For transition metals,

$$
d\un{metal}\ [\mathrm{nm}]
= \frac{6000}{a\un{metal}\ [\mathrm{m^2\,g^{-1}}]\ \times\ \rho\un{metal}\ [\mathrm{g\,cm^{-3}}]} .
$$ (eq-particle-size)

The numerical constant 6000 packages the geometric factor 6 — the surface-to-volume ratio of a
sphere — together with the unit conversion to nanometers. Watch the units: the constant is only 6000
when $a\un{metal}$ is in m² g$^{-1}$ and the density in g cm$^{-3}$. For Ni, with
$\rho = 8.9\ \mathrm{g\,cm^{-3}}$ and $a\un{metal} = 10\ \mathrm{m^2\,g^{-1}}$, this gives
$d\un{metal} = 67\ \mathrm{nm}$.

The **dispersion** characterizes the fraction of metal atoms that are surface atoms. A dispersion of
100 % means every metal atom is accessible to the adsorbate, and therefore active. In practice this
is rarely the case: for supported transition-metal catalysts the dispersion is typically 10–30 %,
meaning a large fraction of the metal is not used.

## Reactions on the surface

<!-- source: HeterogeneousCatalysis.tex L423 -->

Before discussing reactions on surfaces, recall our definition of the reaction rate. For a
homogeneous reaction in a gas or liquid phase we defined it based on the extent of reaction,

$$
r = \frac{1}{V}\frac{\mathrm{d}\xi}{\mathrm{d}t} , \qquad
r\ [=]\ \mathrm{mol\,m^{-3}\,s^{-1}} .
$$

For surface reactions it is more convenient to work with a rate based on the surface area of the
catalyst, since that is where the reaction happens:

$$
r' = \frac{1}{A\un{cat}}\frac{\mathrm{d}\xi}{\mathrm{d}t} , \qquad
r'\ [=]\ \mathrm{mol\,m^{-2}\,s^{-1}} .
$$

You can convert between the two using the catalyst surface area per volume,
$a\un{cat} = A\un{cat}/V\un{reactor}$, so that $r = r'\, a\un{cat}$.

There are two major surface-reaction types: **Langmuir–Hinshelwood** (LH) and **Eley–Rideal** (ER).
The key distinction is the state of the reacting species. In an LH step, *both* reactants are
adsorbed on the surface and react across adjacent sites. In an ER step, an adsorbed species reacts
directly with a gas-phase molecule that has not adsorbed.

### Langmuir–Hinshelwood reactions

<!-- source: HeterogeneousCatalysis.tex L447 -->

For Langmuir–Hinshelwood reactions we assume the surface is well mixed, and that the rate depends on
average concentrations or coverages rather than on the local spatial arrangement of adsorbates.

By far the most common surface reaction type is a dissociation step; the reverse of a dissociation
is an association reaction:

$$
\ce{AB\ast + \ast <=> A\ast + B\ast} .
$$

An adsorbate dissociates over an empty surface site to produce two adsorbates. Using surface
concentrations,

$$
r' = k\un{+}[\ce{AB\ast}][\ce{\ast}] - k\un{-}[\ce{A\ast}][\ce{B\ast}] .
$$

Consider the units. Surface concentrations are in mol m$^{-2}$ and the rate is in
mol m$^{-2}$ s$^{-1}$, so

$$
\mathrm{mol\,m^{-2}\,s^{-1}} = [k]\cdot \mathrm{mol\,m^{-2}} \cdot \mathrm{mol\,m^{-2}} ,
$$

and therefore $k$ must be in m² mol$^{-1}$ s$^{-1}$. Again, it is often more common to work with
coverages,

$$
r'' = k''\un{+}\theta_{\ce{AB\ast}}\theta_{\ce{\ast}}
- k''\un{-}\theta_{\ce{A\ast}}\theta_{\ce{B\ast}} .
$$

Another LH reaction type is an abstraction reaction,

$$
\ce{AB\ast + C\ast <=> A\ast + CB\ast} ,
$$

with rate

$$
r'' = k''\un{+}\theta_{\ce{AB\ast}}\theta_{\ce{C\ast}}
- k''\un{-}\theta_{\ce{A\ast}}\theta_{\ce{CB\ast}} .
$$

The reverse of an abstraction is also an abstraction. Abstraction reactions are rather rare in the
catalysis literature — dissociation steps make up roughly 80–90 % of the elementary steps in
catalytic mechanisms. They are, however, extremely relevant in the combustion community; think of
our branching reactions. Their role in catalysis is an active research direction.

### Eley–Rideal reactions

<!-- source: HeterogeneousCatalysis.tex L486 -->

In an Eley–Rideal reaction, a gas-phase molecule reacts with an adsorbed species to form one or two
adsorbed species,

$$
\begin{aligned}
\ce{A\ast + B &<=> AB\ast} \\
\ce{A\ast + BC &<=> AB\ast + C\ast}
\end{aligned}
$$

Invoking mass-action kinetics,

$$
r' = k\un{+}[\ce{A\ast}]\, p\un{BC} - k\un{-}[\ce{AB\ast}][\ce{C\ast}] .
$$

As with LH reactions this can equivalently be written in coverages; always pay attention to the
units of the rate constant.

:::{admonition} Discussion
:class: seealso
Can you think of a case where Eley–Rideal reactions might be important?
:::

When the surface is highly covered, the ER step might be the dominant pathway. In practice, however,
ER reactions are quite rare, because catalysts are rarely fully covered, and it is still debated
whether genuine Eley–Rideal steps occur at all. In this course we therefore focus on
Langmuir–Hinshelwood reactions; when in doubt, assume an LH mechanism.

## Reaction mechanisms of catalytic reactions

<!-- source: HeterogeneousCatalysis.tex L507 -->

We now have everything we need to construct reaction mechanisms for heterogeneously catalyzed
reactions. We know how to describe adsorption and desorption of reactants and products — more
precisely, we assume throughout that adsorption and desorption are quasi-equilibrated — and the
surface reactions proceed according to a Langmuir–Hinshelwood mechanism. These rate expressions were
pioneered by Hougen and Watson, so the rate laws we obtain are called
**Langmuir–Hinshelwood–Hougen–Watson (LHHW)** rate expressions.

### Simple mechanism

<!-- source: HeterogeneousCatalysis.tex L515 -->

Start with a simple example: a gas-phase molecule A adsorbs and then isomerizes to form a product P
which desorbs into the gas phase,

$$
\begin{aligned}
\ce{A + \ast} &\eqa{k\un{1+}}{k\un{1-}} \ce{A\ast} \\
\ce{A\ast} &\ce{->[$k\un{2}$] P + \ast} .
\end{aligned}
$$

We assume the second step is irreversible and forms a gas-phase product directly. Our goal is a rate
expression for product P; the mass balance is

$$
r = k\un{2}[\ce{A\ast}] .
$$

Our recurring goal is to replace the concentration of something we cannot directly measure — an
adsorbate coverage — with something we can measure, a partial pressure or gas-phase concentration.
Assuming the first step is quasi-equilibrated, this is a simple molecular adsorption, so

$$
[\ce{A\ast}] = \frac{K\un{1}p\un{A}\Gamma}{1 + K\un{1}p\un{A}} ,
$$

and substituting into the rate expression,

$$
r = \frac{k\un{2}K\un{1}p\un{A}\Gamma}{1 + K\un{1}p\un{A}} .
$$ (eq-simple-mech-rate)

Consider two limiting cases:

$$
\begin{aligned}
\text{Weak adsorption, low } p\un{A}: \quad
&K\un{1}p\un{A} \ll 1 \;\Rightarrow\; r\un{A} = k\un{2}K\un{1}p\un{A}\Gamma \\
\text{Strong adsorption, high } p\un{A}: \quad
&K\un{1}p\un{A} \gg 1 \;\Rightarrow\; r\un{A} = k\un{2}\Gamma
\end{aligned}
$$

:::{figure} ../figures/MonomolecularReaction.png
:label: fig-simple-mech
:alt: Hand-drawn plot of the rate of A against the partial pressure of A. The curve rises linearly from the origin, annotated as proportional to k_s times K_ads times p_A, then bends over and saturates at a horizontal dashed line annotated as proportional to k_s. A dashed straight line continues the initial slope to show the departure from linearity.
:width: 50%

Reaction rate for the simple A → P mechanism as a function of $p\un{A}$. At low pressure the rate is
first order in $p\un{A}$, the weak-adsorption limit; at high pressure it saturates at
$k\un{2}\Gamma$, the strong-adsorption limit.
:::

What is the reaction order with respect to A in each limit? The formal definition of the apparent
order is

$$
\alpha\un{A,app} = p\un{A}\frac{\partial \ln r}{\partial p\un{A}} ,
$$ (eq-apparent-order)

and the corresponding apparent activation energy is

$$
E\un{a,app} = RT^2 \frac{\partial \ln k\un{app}}{\partial T} .
$$ (eq-apparent-ea)

The temperature dependence of the elementary rate constant follows Arrhenius and that of the
equilibrium constant follows van 't Hoff,

$$
k\un{2} = A\un{2}\exp\left(\frac{-E\un{a,2}}{RT}\right) , \qquad
K\un{1} = \exp\left(\frac{-\Delta G\un{ads,1}}{RT}\right)
= \exp\left(\frac{\Delta S\un{ads,1}}{R}\right)\exp\left(\frac{-\Delta H\un{ads,1}}{RT}\right) .
$$

**Case 1: strong adsorption ($K\un{1}p\un{A} \gg 1$).** The rate simplifies to $r = k\un{2}\Gamma$,
so $k\un{app} = k\un{2}$. The apparent order is

$$
\alpha\un{A,app} = p\un{A}\frac{\partial \ln(k\un{2}\Gamma)}{\partial p\un{A}} = 0 ,
$$

and the apparent activation energy is simply the elementary barrier,

$$
E\un{a,app} = RT^2 \frac{\partial}{\partial T}\!\left(\ln A\un{2} - \frac{E\un{a,2}}{RT}\right)
= RT^2\left(\frac{E\un{a,2}}{RT^2}\right) = E\un{a,2} .
$$

**Case 2: weak adsorption ($K\un{1}p\un{A} \ll 1$).** The rate simplifies to
$r = k\un{2}K\un{1}p\un{A}\Gamma$, so $k\un{app} = k\un{2}K\un{1}$ and the apparent order is 1.
Applying [](#eq-apparent-ea),

$$
E\un{a,app} = E\un{a,2} + \Delta H\un{ads,1} , \qquad
A\un{app} = A\un{2}\exp\left(\frac{\Delta S\un{ads,1}}{R}\right) .
$$ (eq-weak-ads-ea)

Since the adsorption enthalpy is usually exothermic, $\Delta H\un{ads,1} < 0$, the apparent
activation energy in this limit is *lower* than the elementary activation energy of the
rate-determining surface step.

### LH reaction kinetics: CO oxidation

<!-- source: HeterogeneousCatalysis.tex L585 -->

Now a more complex system: $\ce{CO}$ oxidation on Pd. This is probably the best-studied catalytic
reaction; Gerhard Ertl received his Nobel Prize in part for studying it. The mechanism is

$$
\begin{aligned}
\ce{CO + \ast} &\eqa{k\un{1+}}{k\un{1-}} \ce{CO\ast} \\
\ce{O2 + 2 \ast} &\eqa{k\un{2+}}{k\un{2-}} \ce{2 O\ast} \\
\ce{CO\ast + O\ast} &\ce{->[$k\un{3}$] CO2 + 2 \ast} .
\end{aligned}
$$

We assume the third step, the surface reaction, is the rate-determining step and proceeds
irreversibly, and that the two adsorption steps are quasi-equilibrated. Our rate of interest is the
$\ce{CO2}$ production rate,

$$
r_{\ce{CO2}} = k\un{3}[\ce{CO\ast}][\ce{O\ast}] .
$$

The task is to relate the coverages to measurable quantities. For the quasi-equilibrated reactions
this is straightforward:

$$
\begin{aligned}
K\un{1} = \frac{[\ce{CO\ast}]}{p_{\ce{CO}}[\ce{\ast}]}
\qquad &\text{and} \qquad
K\un{2} = \frac{[\ce{O\ast}]^2}{p_{\ce{O2}}[\ce{\ast}]^2} \\
[\ce{CO\ast}] = K\un{1}p_{\ce{CO}}[\ce{\ast}]
\qquad &\text{and} \qquad
[\ce{O\ast}] = \sqrt{K\un{2}p_{\ce{O2}}}\,[\ce{\ast}]
\end{aligned}
$$

The next step is the site balance,

$$
\begin{aligned}
\Gamma &= [\ce{CO\ast}] + [\ce{O\ast}] + [\ce{\ast}] \\
\Gamma &= K\un{1}p_{\ce{CO}}[\ce{\ast}] + \sqrt{K\un{2}p_{\ce{O2}}}[\ce{\ast}] + [\ce{\ast}] \\
[\ce{\ast}] &= \frac{\Gamma}{1 + K\un{1}p_{\ce{CO}} + \sqrt{K\un{2}p_{\ce{O2}}}} .
\end{aligned}
$$

Inserting this into the expressions for $[\ce{CO\ast}]$ and $[\ce{O\ast}]$,

$$
[\ce{CO\ast}] = \frac{K\un{1}p_{\ce{CO}}\,\Gamma}
{1 + K\un{1}p_{\ce{CO}} + \sqrt{K\un{2}p_{\ce{O2}}}} , \qquad
[\ce{O\ast}] = \frac{\sqrt{K\un{2}p_{\ce{O2}}}\,\Gamma}
{1 + K\un{1}p_{\ce{CO}} + \sqrt{K\un{2}p_{\ce{O2}}}} ,
$$

and finally into the rate expression,

$$
r_{\ce{CO2}}
= \frac{k\un{3}K\un{1}K\un{2}^{1/2}\, p_{\ce{CO}}\, p_{\ce{O2}}^{1/2}\, \Gamma^2}
{\left(1 + K\un{1}p_{\ce{CO}} + \sqrt{K\un{2}p_{\ce{O2}}}\right)^2} .
$$ (eq-co-oxidation-rate)

:::{admonition} Discussion
:class: seealso
What happens to the rate as $p_{\ce{CO}}$ or $p_{\ce{O2}} \to \infty$?
:::

The rate approaches zero in either limit, because the surface saturates with one species and the
other can no longer find a vacant site.

What if $K\un{1}$ is very large, because $\ce{CO\ast}$ binds very strongly? In that limit the rate
simplifies to

$$
r_{\ce{CO2}} = \frac{k\un{3}\sqrt{K\un{2}p_{\ce{O2}}}\,\Gamma^2}{K\un{1}p_{\ce{CO}}} .
$$

The apparent reaction order with respect to $\ce{CO}$ is then $-1$, and with respect to $\ce{O2}$ it
is $+1/2$.

:::{admonition} Discussion
:class: seealso
Does that make sense?
:::

Yes: if $K\un{1} \gg 1$ then almost the entire surface is covered by $\ce{CO}$, and adding more
$\ce{CO}$ only displaces $\ce{O\ast}$ and slows the rate. The analogous case for strong adsorption
of $\ce{O2}$ works the same way. The point is that these orders are exactly what you would measure
experimentally under different conditions.

:::{figure} ../figures/LHReaction.png
:label: fig-lh-co2
:alt: Hand-drawn plot of the rate of A against the partial pressure of A. The curve rises linearly from the origin, annotated rate proportional to p_A, passes through a sharp maximum annotated with the condition that K_ads,A times p_A equals K_ads,B times p_B, then falls away, annotated rate proportional to one over p_A.
:width: 50%

CO-oxidation rate as a function of $p_{\ce{CO}}$ at fixed $p_{\ce{O2}}$ in a Langmuir–Hinshelwood
mechanism. The rate has a maximum because the two reactants compete for the same sites: at low
$p_{\ce{CO}}$ the rate is limited by $\ce{CO\ast}$ coverage, while at high $p_{\ce{CO}}$ the surface
is poisoned by $\ce{CO}$ and starves $\ce{O\ast}$.
:::

## Eley–Rideal reaction mechanism

<!-- source: HeterogeneousCatalysis.tex L652 -->

In an Eley–Rideal mechanism only one of the reactants is adsorbed on the catalyst surface; the other
reacts directly with the adsorbed species without adsorbing first. Such a sequence is

$$
\begin{aligned}
\ce{A + \ast} &\eqa{k\un{1+}}{k\un{1-}} \ce{A\ast} \\
\ce{A\ast + B} &\ce{->[$k\un{2}$] C\ast}
\end{aligned}
$$

The derivation follows the same steps as before. Assume the adsorption reaction is equilibrated and
the second reaction is rate-limiting, so

$$
r\un{A} = k\un{2}[\ce{A\ast}]\, p\un{B} .
$$

From the equilibrium assumption,
$K\un{ads,A} = [\ce{A\ast}]/(p\un{A}[\ce{\ast}])$, and the site balance for this system is
$\Gamma = [\ce{A\ast}] + [\ce{\ast}]$. Therefore the overall rate is

$$
r\un{A} = \frac{k\un{2}K\un{ads,A}\, p\un{A}\, p\un{B}\, \Gamma}{1 + K\un{ads,A}\, p\un{A}} .
$$ (eq-er-rate)

:::{figure} ../figures/ERReaction.png
:label: fig-er-rate-pa
:alt: Hand-drawn plot of the rate of A against the partial pressure of A. The curve rises from the origin, annotated as proportional to k_s times K_ads,A times p_A times p_B, then bends over and saturates at a plateau annotated as proportional to k_s times p_B. Unlike the Langmuir-Hinshelwood case there is no maximum; the curve is monotonic.
:width: 50%

Reaction rate of an Eley–Rideal reaction as a function of $p\un{A}$ at fixed $p\un{B}$. Unlike the LH
case, the rate keeps increasing with $p\un{A}$ and saturates at $k\un{2}p\un{B}\Gamma$ rather than
passing through a maximum, because B is not competing for surface sites.
:::

## Most abundant reaction intermediate

<!-- source: HeterogeneousCatalysis.tex L689 -->

Typically multiple adsorbates compete for the various adsorption sites on the catalyst surface, and
we have seen that small changes in the equilibrium constant can drastically change whether A or B
covers the surface. Quite often one adsorbate binds more tightly than all the others and covers most
of the surface. We call it the **most abundant reaction intermediate (MARI)**; you will also see
MASI, for most abundant surface intermediate, which means the same thing.

The MARI concept lets us greatly simplify the kinetic expression, because we no longer have to
consider all quasi-equilibrated steps. Let us demonstrate this for ammonia synthesis in the
Haber–Bosch process, with the mechanism

$$
\begin{aligned}
\ce{N2 + 2 \ast &<=> 2 N\ast} \\
\ce{H2 + 2 \ast &<=> 2 H\ast} \\
\ce{N\ast + H\ast &<=> NH\ast + \ast} \\
\ce{NH\ast + H\ast &<=> NH2\ast + \ast} \\
\ce{NH2\ast + H\ast &<=> NH3 + 2 \ast}
\end{aligned}
$$

We assume the dissociative adsorption of $\ce{N2}$ is the rate-determining step, that all other
elementary reactions are quasi-equilibrated, and that $\ce{H\ast}$ is the MARI.

:::{admonition} Discussion
:class: seealso
How could we arrive at this assumption?
:::

Experimentally — via infrared spectroscopy of the surface under reaction conditions, for instance —
or computationally, from first-principles calculations of binding energies.

With these two assumptions we can lump multiple elementary steps together to simplify the algebra.
First, our reaction rate is limited by the rate-determining step,

$$
r = r\un{1} = k\un{1}p_{\ce{N2}}[\ast]^2 .
$$

Our only unknown is the concentration of vacant sites, which we obtain from the site balance,

$$
\Gamma = [\ce{H\ast}] + [\ast] .
$$

We do not need the concentrations of the other adsorbates, because they are very small compared with
$[\ce{H\ast}]$. For the dissociative adsorption of $\ce{H2}$ under quasi-equilibrium,

$$
K\un{2} = \frac{[\ce{H\ast}]^2}{p_{\ce{H2}}[\ast]^2}
\quad \Rightarrow \quad
[\ce{H\ast}] = [\ast]\sqrt{K\un{2}p_{\ce{H2}}} ,
$$

so the site balance gives

$$
[\ast] = \frac{\Gamma}{1 + \sqrt{K\un{2}p_{\ce{H2}}}} ,
$$

and finally

$$
r = \frac{k\un{1}p_{\ce{N2}}\Gamma^2}{\left(1 + \sqrt{K\un{2}p_{\ce{H2}}}\right)^2} .
$$ (eq-mari-h)

This shows that we need not account for any reaction that is neither the rate-determining step nor
the step that forms the MARI. Sometimes we do have to consider multiple steps, because the MARI is
not formed through an adsorption step. We could equally have lumped all the remaining steps together
and written the mechanism as

$$
\begin{aligned}
\ce{N2 + 2 \ast &<=> 2 N\ast} \\
\ce{H2 + 2 \ast &<=> 2 H\ast} \\
\ce{N\ast + 3 H\ast &<=> NH3 + 4 \ast}
\end{aligned}
$$

### Second example: $\ce{N\ast}$ as the MARI

<!-- source: HeterogeneousCatalysis.tex L753 -->

Consider another scenario for ammonia synthesis,

$$
\begin{aligned}
\ce{N2 + 2 \ast &<=> 2 N\ast} \\
\ce{N\ast + 3/2 H2 &<=> NH3 + \ast}
\end{aligned}
$$

where the last step is again not elementary but lumped. Assume the first reaction is the RDS and
that $\ce{N\ast}$ is the MARI — which could happen if you change the catalyst. The rate of the RDS
is

$$
r = r\un{1} = k\un{1}p_{\ce{N2}}[\ast]^2 ,
$$

and the site balance is $\Gamma = [\ast] + [\ce{N\ast}]$. To eliminate the concentration of vacant
sites we invoke quasi-equilibrium of the second reaction,

$$
K\un{2} = \frac{p_{\ce{NH3}}[\ast]}{[\ce{N\ast}]\, p_{\ce{H2}}^{3/2}}
\quad \Rightarrow \quad
[\ce{N\ast}] = \frac{p_{\ce{NH3}}[\ast]}{p_{\ce{H2}}^{3/2}K\un{2}} ,
$$

so that

$$
[\ast] = \frac{\Gamma}{1 + \dfrac{p_{\ce{NH3}}}{p_{\ce{H2}}^{3/2}K\un{2}}} ,
$$

and accordingly

$$
r = \frac{k\un{1}p_{\ce{N2}}\Gamma^2}
{\left(1 + \dfrac{p_{\ce{NH3}}}{p_{\ce{H2}}^{3/2}K\un{2}}\right)^2} .
$$ (eq-mari-n)

This example emphasizes that it does not matter what happens after the RDS — you can simply merge
everything together.

## Multiple active sites

<!-- source: HeterogeneousCatalysis.tex L790 -->

There are frequently multiple types of active site on a catalyst surface. Some adsorbates only
adsorb on one site type, others on a different type. Consider the mechanism

$$
\begin{aligned}
\ce{A + \ast_1 &<=> A\ast_1} \\
\ce{B + \ast_2 &<=> B\ast_2} \\
\ce{A\ast_1 + B\ast_2 &<=> C\ast_1 + \ast_2} \\
\ce{C\ast_1 &<=> C + \ast_1} ,
\end{aligned}
$$

where $\ast_1$ is the active site for species A and C, and $\ast_2$ the active site for species B.

Assume the surface reaction is the RDS and irreversible, so

$$
r = k\un{3}[\ce{A\ast_1}][\ce{B\ast_2}] .
$$

Assuming quasi-equilibrium for all other steps,

$$
\begin{aligned}
K\un{1} = \frac{[\ce{A\ast_1}]}{p\un{A}[\ast_1]}
\quad &\Rightarrow \quad [\ce{A\ast_1}] = K\un{1}p\un{A}[\ast_1] \\
K\un{2} = \frac{[\ce{B\ast_2}]}{p\un{B}[\ast_2]}
\quad &\Rightarrow \quad [\ce{B\ast_2}] = K\un{2}p\un{B}[\ast_2] \\
K\un{4} = \frac{p\un{C}[\ast_1]}{[\ce{C\ast_1}]}
\quad &\Rightarrow \quad [\ce{C\ast_1}] = \frac{p\un{C}}{K\un{4}}[\ast_1] .
\end{aligned}
$$

We now have to perform *two* site balances,

$$
\Gamma\un{1} = [\ast_1] + [\ce{A\ast_1}] + [\ce{C\ast_1}] , \qquad
\Gamma\un{2} = [\ast_2] + [\ce{B\ast_2}] ,
$$

which give

$$
[\ast_1] = \frac{\Gamma\un{1}}{1 + K\un{1}p\un{A} + \dfrac{p\un{C}}{K\un{4}}} , \qquad
[\ast_2] = \frac{\Gamma\un{2}}{1 + K\un{2}p\un{B}} .
$$

Inserting into the rate expression, we arrive at

$$
r = \frac{k\un{3}K\un{1}p\un{A}K\un{2}p\un{B}\Gamma\un{1}\Gamma\un{2}}
{\left(1 + K\un{1}p\un{A} + \dfrac{p\un{C}}{K\un{4}}\right)\left(1 + K\un{2}p\un{B}\right)} .
$$ (eq-two-site-rate)

## Langmuir–Hinshelwood–Hougen–Watson rate equations

<!-- source: HeterogeneousCatalysis.tex L828 -->

We have now derived a few LH-type rate expressions. Do they respect thermodynamic equilibrium? They
do not, because we have so far written the surface reaction as irreversible. Real rates also
decrease monotonically as the system approaches the equilibrium composition.

To handle both issues at once, Hougen and Watson generalized many of these LH derivations and
arrived at a generic rate equation of the form

$$
r = \frac{(\text{kinetic term})\,(\text{driving force / displacement from equilibrium})}
{(\text{adsorption term})^n} ,
$$

for example, for an $\ce{A + B <=> C}$ reaction,

$$
r\un{A} = \frac{\left(k\un{s}K\un{ads,A}K\un{ads,B}\right)
\left(p\un{A}p\un{B} - \dfrac{p\un{C}p\un{0}}{K}\right)}
{\left(1 + K\un{ads,A}p\un{A} + K\un{ads,B}p\un{B} + K\un{ads,C}p\un{C}\right)^2} .
$$ (eq-lhhw)

The numerator factor is the displacement from equilibrium,

$$
p\un{A}p\un{B} - p\un{A}^{\ast}p\un{B}^{\ast}
= p\un{A}p\un{B} - \frac{p\un{C}p\un{0}}{K} ,
\qquad \text{with} \qquad
K = \frac{p\un{C}p\un{0}}{p\un{A}p\un{B}} ,
$$

where the starred pressures denote equilibrium values. Note that $K$ here is the *dimensionless*
equilibrium constant of [](#ch-thermodynamics), not the dimensional
$K\un{p} = \prod_i p_i^{\nu_i}$; the factor $p\un{0}$, the reference pressure, is what makes it so
for this reaction, where $\sum_i \nu_i = -1$.

LHHW expressions are widely used in the catalysis community and in industrial applications.

## Rate constants for catalyzed reactions

<!-- source: HeterogeneousCatalysis.tex L850 -->

We now have an expression for the rate of $\ce{CO2}$ formation, but we still do not know the
numerical values of $k\un{3}$, $K\un{1}$, $K\un{2}$, and so on.

:::{admonition} Discussion
:class: seealso
How do we derive these? And how could we confirm the apparent reaction order?
:::

Either by fitting to experimental data or by computing them directly from quantum-chemical
calculations. To confirm an apparent order, perform experiments in which only the partial pressure
of the species of interest is varied.

In general you will not be able to distinguish all the individual parameters from a fit, so we group
them together,

$$
r_{\ce{CO2}} = \frac{k\, p_{\ce{CO}}\, p_{\ce{O2}}^{1/2}}
{\left(1 + K\un{1}p_{\ce{CO}} + \sqrt{K\un{2}p_{\ce{O2}}}\right)^2} ,
$$

and then regress the lumped parameters against experimental data.

:::{admonition} Discussion
:class: seealso
But how do I know my reaction-rate expression is the right one? I could have assumed a different
MARI or a different RDS.
:::

You usually derive 20–30 candidate LHHW expressions based on different assumptions about MARI and
RDS, informed by experiments and the literature, and then see which one fits best. LHHW expressions
work very well in practice and are routinely used in catalysis research — but you cannot obtain
mechanistic insight from a fitted LHHW expression alone. If it fits the data and lets us design a
reactor, we have achieved what we set out to do.

:::{admonition} Discussion
:class: seealso
How do we achieve a full mechanistic understanding?
:::

Not from fitting a rate expression to experiments alone. Instead you need **microkinetic modeling**:
determining the rate constants and equilibrium constants from first principles, typically via
transition state theory and quantum-chemical calculations.

## Mass transfer

<!-- source: HeterogeneousCatalysis.tex L882 -->

Heterogeneous catalysis is really a combination of microkinetics and mass transfer; together these
define the macrokinetics. Reactants must reach the active sites through external (film) and internal
(pore) diffusion before they can react, and products must diffuse out by the same route. In the
isothermal case the highest local rate occurs at the outside of the catalyst particle, where the
reactant concentration is highest. As a result only part of the catalyst volume is fully utilized:
most supported catalysts operate at less than 100 % effectiveness, and the effective utilization
drops further as the intrinsic rate becomes faster.

:::{figure} ../figures/CProfiles.png
:label: fig-c-profile
:alt: Hand-drawn concentration profile across a catalyst particle, with the particle centre in the middle and its outer surface at plus and minus d_p over 2, surrounded by a film layer of thickness delta. A black dashed curve dips only slightly from the bulk value toward the centre. A blue curve, for a faster intrinsic rate, drops noticeably across the film and again inside the particle, reaching its minimum at the centre; the blue point where it meets the outer surface is annotated as the maximum rate.
:width: 55%

Concentration profile of a reactant inside a porous catalyst particle when diffusion is fast
compared with reaction. The concentration is nearly uniform from the outer surface to the center,
and the entire particle is utilized.
:::

:::{figure} ../figures/CProfilesLimited.png
:label: fig-c-profile-limited
:alt: Hand-drawn concentration profile across a catalyst particle. A green horizontal line at the top marks the ideal case with no gradient. The black curve falls slightly across the film on each side, labelled film diffusion, then plunges steeply inside the particle, labelled pore diffusion, reaching zero well before the centre. The central region where the concentration is zero is hatched in red and labelled dead zone.
:width: 55%

Concentration profile when pore diffusion is slow compared with reaction. The reactant is consumed
before it can reach the center, so the inner core of the particle sees almost no reactant and
contributes little to the overall rate.
:::

To derive the governing equation, we balance diffusion against reaction inside the particle. For a
first-order reaction $r = k\,c$,

$$
\underbrace{D\un{eff}\frac{\mathrm{d}^2 c}{\mathrm{d}x^2}}_{\text{diffusion}}
= \underbrace{k\, c}_{\text{reaction}} .
$$ (eq-diffusion-reaction)

Measuring $x$ from the centre of the particle, the boundary conditions are

$$
c(L) = c\un{s} \quad \text{(concentration at outer surface)} , \qquad
\left.\frac{\mathrm{d}c}{\mathrm{d}x}\right|_{x=0} = 0 \quad \text{(symmetry at center)} .
$$

Integration yields

$$
\frac{c(x)}{c\un{s}} = \frac{\cosh\left(\Phi\, x/L\right)}{\cosh \Phi} ,
$$

which is flat at the centre and rises to $c\un{s}$ at the outer surface, as it must. Here $\Phi$ is
the **Thiele modulus**,

$$
\Phi = L\sqrt{\frac{k}{D\un{eff}}}
= \sqrt{\frac{\text{time constant of diffusion}}{\text{time constant of reaction}}} ,
$$ (eq-thiele)

since the diffusion time constant is $L^2/D\un{eff}$ and the reaction time constant is $1/k$, so
that $\Phi^2 = kL^2/D\un{eff}$. The characteristic length $L$ is $V\un{cat}/A\un{cat}$ for a generic
particle, $d\un{p}/6$ for a sphere, and $d/2$ for a plate.

We define the catalyst **effectiveness factor** as

$$
\eta = \frac{\bar{r}}{r\un{max}} = \frac{k\,\bar{c}}{k\, c\un{s}} ,
$$

and for a first-order reaction in a slab the integration gives

$$
\eta = \frac{\tanh \Phi}{\Phi} ,
$$ (eq-effectiveness)

which is approximately $1/\Phi$ for $\Phi \gtrsim 2$.

:::{figure} ../figures/ApparenEaMassTransfer.png
:label: fig-ea-masstransfer
:alt: Hand-drawn Arrhenius plot of the logarithm of the effective rate constant against one over temperature. The line has three segments of increasing steepness from left to right. The shallowest, at the left, is labelled film diffusion limited with an apparent barrier of 5 to 20 kilojoules per mole. The middle segment is labelled pore diffusion limited with slope corresponding to E_a over 2. The steepest, at the right, is the kinetic regime with the full E_a.
:width: 55%

Arrhenius plot, $\ln r$ vs. $1/T$, showing the transition from kinetic to mass-transfer-limited
operation. At low temperature the apparent activation energy equals the intrinsic value; at high
temperature the rate becomes diffusion-controlled and the apparent activation energy drops to half
the intrinsic value, and eventually to that of the diffusion process.
:::

At high temperatures the intrinsic rate constant grows exponentially with $T$ while the diffusion
coefficient grows only weakly, so the reaction becomes mass-transfer limited and the apparent
activation energy drops below the intrinsic value.

To approach $\eta \to 1$ we need to keep $\Phi$ small, which means small particles, small intrinsic
rate constants, or large effective diffusion coefficients. Practically this is achieved with small
catalyst particles, with core-shell particles where only the outer shell is catalytically active, or
with monolith structures. A common rule of thumb is to aim for $\Phi < 0.3$.

:::{figure} ../figures/CoreShell.png
:label: fig-coreshell
:alt: Two hand-drawn sketches. On the left, a core-shell particle: an irregular black blob labelled solid, surrounded by a hatched blue rim labelled catalyst, so only the thin outer layer is active. On the right, a monolith: a square channel drawn in cross-section with a hatched wall, the wall labelled catalyst, so the reacting gas flows down the open centre past a thin active layer.
:width: 55%

Core-shell catalyst design, left: a non-active core surrounded by a thin catalytically active shell.
The short diffusion path through the shell keeps $\Phi$ small and the effectiveness factor near
unity. Right: a monolith channel, which achieves the same thing with a thin active layer on the
channel wall.
:::

## Summary

<!-- source: HeterogeneousCatalysis.tex L950 -->

- A catalytic cycle has seven steps: film diffusion in, pore diffusion in, adsorption, surface
  reaction, desorption, pore diffusion out, film diffusion out. Steps 3–5 are the chemistry; the
  rest is transport.
- Physisorption is weak, roughly $-5$ to $-40\ \mathrm{kJ\,mol^{-1}}$, and does not activate the
  adsorbate; chemisorption is strong, roughly $-40$ to $-500\ \mathrm{kJ\,mol^{-1}}$, and can be
  dissociative.
- Under the Langmuir assumptions — monolayer, equivalent sites, adsorption only at vacant sites, no
  lateral interactions — the equilibrium coverage is
  $\theta_{\ce{A\ast}} = K\un{ads}p\un{A}/(1 + K\un{ads}p\un{A})$ for molecular adsorption,
  [](#eq-langmuir-isotherm); $K_i p_i/(1 + \sum_j K_j p_j)$ for competitive adsorption,
  [](#eq-competitive-isotherm); and
  $\sqrt{K\un{ads}p}/(1 + \sqrt{K\un{ads}p})$ for dissociative adsorption,
  [](#eq-dissociative-isotherm).
- Because adsorption is exothermic, coverage *decreases* with increasing temperature.
- Chemisorption of a probe molecule such as $\ce{H2}$ counts active sites, [](#eq-site-count), from
  which the metallic surface area, [](#eq-metal-area), the mean particle size,
  [](#eq-particle-size), and the dispersion follow. Supported transition-metal catalysts typically
  achieve only 10–30 % dispersion.
- Langmuir–Hinshelwood steps have both reactants adsorbed; Eley–Rideal steps have one reactant in
  the gas phase. LH is the default assumption — genuine ER steps are rare.
- The recipe for an LHHW rate law is always the same: write the rate of the rate-determining step,
  use quasi-equilibrium to replace unmeasurable coverages with partial pressures, and close the
  system with a site balance.
- Competition for sites is what makes catalytic rate laws non-monotonic. In $\ce{CO}$ oxidation the
  rate passes through a maximum, [](#eq-co-oxidation-rate): strong $\ce{CO}$ binding gives an
  apparent order of $-1$ in $\ce{CO}$ and $+1/2$ in $\ce{O2}$.
- Apparent orders and apparent activation energies are limit-dependent. In the weak-adsorption limit
  $E\un{a,app} = E\un{a,2} + \Delta H\un{ads,1}$, [](#eq-weak-ads-ea), which is *lower* than the
  intrinsic barrier because the adsorption enthalpy is negative.
- A fitted LHHW expression that reproduces the data is enough to design a reactor, but it does not
  establish the mechanism; that requires microkinetic modelling from first principles.
- Pore diffusion limits utilization of the particle. The Thiele modulus
  $\Phi = L\sqrt{k/D\un{eff}}$, [](#eq-thiele), sets the effectiveness factor
  $\eta = \tanh\Phi/\Phi$, [](#eq-effectiveness), which is $\approx 1/\Phi$ for $\Phi \gtrsim 2$. Aim
  for $\Phi < 0.3$ via small particles, core-shell designs, or monoliths. Under strong
  pore-diffusion limitation the apparent activation energy falls to about half its intrinsic value.
