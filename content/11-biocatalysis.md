---
title: Biocatalysis
short_title: Biocatalysis
label: ch-biocatalysis
---

<!-- LaTeX source: Biocatalysis.tex -->
<!-- Porting notes: mhchem is math-only here, so \ce{} in prose needs $...$; and \un{} expands to
     _{\textrm{...}}, so species subscripts must be written _{\ce{S}}, never \un{\ce{S}}.
     mhchem arrow labels contain their own $, so ->[$k$] must live in a $$ block, never in inline
     math -- inline it terminates the math early and dumps raw macro source onto the page.
     Nested directives need a longer outer fence (:::: around :::). -->

:::{admonition} Learning objectives
:class: tip
After completing this chapter, you should be able to:

- Explain what distinguishes an enzyme from the heterogeneous catalysts of [](#ch-hetcat), and list
  the advantages and the main drawback of enzymatic catalysis.
- Classify enzymes by the reaction type they catalyze.
- Write the governing ODEs for the enzyme–substrate mechanism and reduce them with the PSSA on the
  enzyme–substrate complex.
- Derive the Michaelis–Menten rate law
  $r = r\un{max}[\ce{S}]/(K\un{M} + [\ce{S}])$ using an enzyme balance, in direct analogy with the
  catalyst site balance.
- Identify the first-order and zeroth-order limits of the Michaelis–Menten expression and interpret
  $r\un{max}$ and $K\un{M}$ physically.
- Read $K\un{M}$ off a rate-versus-substrate plot as the concentration at half the maximum rate.
:::

We have just wrapped up our discussion of heterogeneous catalysts, and talked about how to use
microkinetic modeling to identify the best catalyst for $\ce{NH3}$ synthesis. As you hopefully
remember, the industrial process operates between 400 and 500 °C and between 150 and 300 bar, and
only works with a catalyst — without one the reaction would be far too slow. We need a great deal of
energy to break the $\ce{N#N}$ triple bond to make $\ce{NH3}$, which then ends up as fertilizer and
eventually on our plate.

And yet nature has figured out how to fix nitrogen at room temperature and ambient pressure. How?
With an enzyme: nitrogenase does what no metal catalyst is capable of.

## Enzymes

<!-- source: Biocatalysis.tex L25 -->

Enzymes are catalysts that occur in nature and are essential for life. In their natural state they
are *homogeneous* catalysts: dissolved in the same, usually aqueous, phase as their substrate.
Industrially they are sometimes immobilized on a support to become heterogeneous, but that is a
deliberate engineering choice.

:::{admonition} Discussion
:class: seealso
Who has taken biochemistry or already knows about enzymes?
:::

We will focus on enzymes as catalysts. In this regard they are truly amazing structures, and
superior to most other catalysts:

- very specific and selective,
- give stereochemically pure products,
- work at ambient conditions.

The main downside is that they are fragile.

Enzymes are proteins with a very specific function: to catalyze a reaction. Like all proteins they
are made up of amino acids, nature's building blocks. They are classified by the reaction type they
catalyze:

:::{table} Classification of enzymes by the reaction type they catalyze.
:label: tab-enzyme-classes

| Class | What it does | Generic reaction |
|:------|:-------------|:-----------------|
| **oxidoreductases** | transfer an electron   | $\ce{A^- + B -> A + B^-}$ |
| **transferases**    | transfer functional groups | $\ce{AX + B -> A + BX}$ |
| **hydrolases**      | use water to cleave a bond | $\ce{AB + H2O -> AOH + BH}$ |
| **lyases**          | break a bond           | $\ce{AB -> A + B}$ |
| **isomerases**      | isomerize              | $\ce{A -> B}$ |
| **ligases**         | join two molecules     | $\ce{A + B -> AB}$ |
:::

## Reaction kinetics of enzyme-catalyzed reactions

<!-- source: Biocatalysis.tex L53 -->

We do not have time to dive deep into the chemistry of enzymes. In this course we are mostly
interested in the kinetics of enzyme-catalyzed reactions and what influences them.

Some terminology:

- **enzyme** — the catalyst, E
- **substrate** — the reactant we are trying to convert, S
- **enzyme–substrate complex** — substrate bound to the enzyme, ES
- **product** — the result of the reaction, P

The reaction sequence for an enzymatic reaction is

$$
\begin{aligned}
\ce{E + S} &\eqa{k\un{1+}}{k\un{1-}} \ce{ES} \\
\ce{ES} &\ce{->[$k\un{2}$] E + P} .
\end{aligned}
$$ (eq-enzyme-mechanism)

:::{admonition} Discussion
:class: seealso
What are the governing equations?
:::

$$
\begin{aligned}
\frac{\mathrm{d}[\ce{E}]}{\mathrm{d}t} &= -r\un{1} + r\un{2}
= -k\un{1+}[\ce{E}][\ce{S}] + k\un{1-}[\ce{ES}] + k\un{2}[\ce{ES}] \\
\frac{\mathrm{d}[\ce{S}]}{\mathrm{d}t} &= -r\un{1}
= -k\un{1+}[\ce{E}][\ce{S}] + k\un{1-}[\ce{ES}] \\
\frac{\mathrm{d}[\ce{ES}]}{\mathrm{d}t} &= r\un{1} - r\un{2}
= k\un{1+}[\ce{E}][\ce{S}] - k\un{1-}[\ce{ES}] - k\un{2}[\ce{ES}] \\
\frac{\mathrm{d}[\ce{P}]}{\mathrm{d}t} &= r\un{2} = k\un{2}[\ce{ES}]
\end{aligned}
$$ (eq-enzyme-odes)

:::{admonition} Discussion
:class: seealso
Which assumption could we invoke to simplify this system?
:::

We assume that the enzyme–substrate complex is present in small, slowly changing amounts compared
with the free substrate and the total enzyme, so we can invoke the PSSA on $\ce{ES}$. Note that we
are not claiming $k_2$ is large; the PSSA holds because $[\ce{ES}]$ is small and its time derivative
is negligible compared with the individual production and consumption terms. We also assume the
second reaction is irreversible.

Invoking the PSSA for $\ce{ES}$ leads to

$$
\begin{aligned}
\frac{\mathrm{d}[\ce{ES}]}{\mathrm{d}t} &= 0
= k\un{1+}[\ce{E}][\ce{S}] - k\un{1-}[\ce{ES}] - k\un{2}[\ce{ES}] \\
k\un{1+}[\ce{E}][\ce{S}] &= \left(k\un{1-} + k\un{2}\right)[\ce{ES}] .
\end{aligned}
$$

Imagine we have a fixed number of enzymes. We do not know how many of them are bound and how many
are free. This is very similar to the catalyst surface, where we have a fixed number of surface
sites that can be either occupied or vacant. Making a balance over all enzymes leads to

$$
[\ce{E}]_0 = [\ce{E}] + [\ce{ES}]
\quad \Rightarrow \quad
[\ce{E}] = [\ce{E}]_0 - [\ce{ES}] ,
$$ (eq-enzyme-balance)

which we can plug back into the PSSA result:

$$
\begin{aligned}
k\un{1+}[\ce{E}][\ce{S}] &= \left(k\un{1-} + k\un{2}\right)[\ce{ES}] \\
k\un{1+}\left([\ce{E}]_0 - [\ce{ES}]\right)[\ce{S}]
&= \left(k\un{1-} + k\un{2}\right)[\ce{ES}] \\
k\un{1+}[\ce{E}]_0[\ce{S}] - k\un{1+}[\ce{ES}][\ce{S}]
&= \left(k\un{1-} + k\un{2}\right)[\ce{ES}] \\
k\un{1+}[\ce{E}]_0[\ce{S}]
&= \left(k\un{1-} + k\un{2} + k\un{1+}[\ce{S}]\right)[\ce{ES}] \\
[\ce{ES}] &= \frac{k\un{1+}[\ce{E}]_0[\ce{S}]}{k\un{1-} + k\un{2} + k\un{1+}[\ce{S}]} \\
[\ce{ES}] &= \frac{[\ce{E}]_0[\ce{S}]}
{\left(\dfrac{k\un{1-} + k\un{2}}{k\un{1+}}\right) + [\ce{S}]} .
\end{aligned}
$$

Finally we plug this into the expression for the rate of product formation,

$$
r = \frac{\mathrm{d}[\ce{P}]}{\mathrm{d}t}
= \frac{k\un{2}[\ce{E}]_0[\ce{S}]}
{\left(\dfrac{k\un{1-} + k\un{2}}{k\un{1+}}\right) + [\ce{S}]} .
$$

A couple of substitutions make this easier to read. First the denominator,

$$
K\un{M} = \frac{k\un{1-} + k\un{2}}{k\un{1+}} ,
$$ (eq-km)

where $K\un{M}$ is the **Michaelis–Menten constant**, named after Leonor Michaelis and Maud Menten
for their work on unraveling the kinetics of enzyme-catalyzed reactions. Both rate constants in the
numerator, $k\un{1-}$ and $k\un{2}$, are first order with units s$^{-1}$, while $k\un{1+}$ in the
denominator is second order with units m³ mol$^{-1}$ s$^{-1}$. As a result $K\un{M}$ has units of
concentration, mol m$^{-3}$.

Next, a maximum rate. Taking our rate expression for the product,
$\mathrm{d}[\ce{P}]/\mathrm{d}t = k\un{2}[\ce{ES}]$:

:::{admonition} Discussion
:class: seealso
What is the highest reaction rate we can achieve?
:::

The maximum occurs when 100 % of the enzyme is bound as the ES complex, so that
$[\ce{ES}] = [\ce{E}]_0$, giving

$$
r\un{max} = k\un{2}[\ce{E}]_0 .
$$ (eq-rmax)

Substituting into the general expression,

$$
r = \frac{k\un{2}[\ce{E}]_0[\ce{S}]}
{\left(\dfrac{k\un{1-} + k\un{2}}{k\un{1+}}\right) + [\ce{S}]}
= \frac{r\un{max}[\ce{S}]}{K\un{M} + [\ce{S}]} .
$$ (eq-michaelis-menten)

The rate now depends only on the substrate concentration and on two lumped parameters, $r\un{max}$
and $K\un{M}$.

## Michaelis–Menten cases

<!-- source: Biocatalysis.tex L148 -->

Let us look at the limiting cases.

### Case 1: low substrate, $[\ce{S}] \ll K\un{M}$

$$
r = \frac{r\un{max}[\ce{S}]}{K\un{M} + [\ce{S}]}
\;\xrightarrow{[\ce{S}] \ll K\un{M}}\;
\frac{r\un{max}}{K\un{M}}[\ce{S}] .
$$

At low substrate concentration the reaction is first order in S: vacant enzyme is plentiful, and the
rate is limited by how quickly substrate finds it.

:::{figure} ../figures/10_Michaelis_Menten_1.png
:label: fig-michaelis-menten-1
:alt: Slide headed "Case I, low substrate concentration". It states the Michaelis-Menten rate law and its limit as the substrate concentration becomes much smaller than K_M, namely r_max over K_M times the substrate concentration. Below, a plot of rate against substrate concentration shows a short blue straight line rising from the origin, annotated in red as first-order in substrate.
:width: 60%

Michaelis–Menten rate $r$ as a function of substrate concentration $[\ce{S}]$ in the low-$[\ce{S}]$
limit. The rate rises linearly with $[\ce{S}]$ with slope $r\un{max}/K\un{M}$ — first-order
kinetics.
:::

### Case 2: high substrate, $[\ce{S}] \gg K\un{M}$

$$
r \;\xrightarrow{[\ce{S}] \gg K\un{M}}\; r\un{max} .
$$

At high substrate concentration the rate becomes independent of $[\ce{S}]$ — zeroth order. This
makes sense: with abundant substrate, essentially all enzyme is tied up as ES, and the turnover step
$k\un{2}[\ce{ES}]$ sets the rate.

:::{figure} ../figures/10_Michaelis_Menten_2.png
:label: fig-michaelis-menten-2
:alt: Slide headed "Case II, high substrate concentration". It states the Michaelis-Menten rate law and its limit as the substrate concentration becomes much larger than K_M, namely r_max. Below, a plot of rate against substrate concentration shows the short rising blue line at low concentration and, separately at high concentration, a flat blue line sitting on the dashed r_max level, annotated in red as zeroth-order in substrate.
:width: 60%

Michaelis–Menten rate $r$ in the high-$[\ce{S}]$ limit. The rate saturates at
$r\un{max} = k\un{2}[\ce{E}]_0$ — zeroth-order kinetics. Adding more substrate does not increase the
rate, because all enzyme is already bound.
:::

### Case 3: intermediate substrate

For $[\ce{S}] \sim K\un{M}$ the rate is a smooth saturating function connecting the two limits, with
the exact shape set by the values of $r\un{max}$ and $K\un{M}$.

One more special point on the curve is worth highlighting. Let $[\ce{S}]_{1/2}$ denote the substrate
concentration at which $r = r\un{max}/2$. Then

$$
\begin{aligned}
\frac{1}{2}r\un{max} &= \frac{r\un{max}\,[\ce{S}]_{1/2}}{K\un{M} + [\ce{S}]_{1/2}} \\
\frac{1}{2}\bigl(K\un{M} + [\ce{S}]_{1/2}\bigr) &= [\ce{S}]_{1/2} \\
K\un{M} &= [\ce{S}]_{1/2} .
\end{aligned}
$$ (eq-km-half-max)

The Michaelis–Menten constant $K\un{M}$ is therefore the substrate concentration at which the enzyme
operates at half its maximum rate. This gives $K\un{M}$ a direct experimental meaning, and a
convenient way to read it off a rate-versus-$[\ce{S}]$ plot.

:::{figure} ../figures/10_Michaelis_Menten_3.png
:label: fig-michaelis-menten-3
:alt: Slide headed "Typical Michaelis-Menten curve", stating the rate law. Below, a plot of rate against substrate concentration shows a single blue curve rising from the origin and saturating at a dashed line marked r_max. A second dashed line at r_max over 2 meets the curve at the substrate concentration marked K_M on the horizontal axis. The initial rising portion is annotated in red as first-order in substrate and the saturated portion as zeroth-order in substrate.
:width: 60%

Full Michaelis–Menten curve. The rate approaches $r\un{max}$ at high substrate concentration; the
half-maximum rate $r\un{max}/2$ is reached at $[\ce{S}] = K\un{M}$, giving a graphical
interpretation of the Michaelis–Menten constant.
:::

## Summary

<!-- source: Biocatalysis.tex L197 -->

- Enzymes are protein catalysts. They are highly specific and selective, give stereochemically pure
  products, and work at ambient conditions — nitrogenase fixes $\ce{N2}$ where the Haber–Bosch
  process needs hundreds of degrees and hundreds of bar. Their main drawback is fragility.
- In their natural state enzymes are homogeneous catalysts; immobilizing them on a support to make
  them heterogeneous is an engineering choice, not a property of the enzyme.
- The enzyme mechanism $\ce{E + S <=> ES}$ followed by $\ce{ES -> E + P}$,
  [](#eq-enzyme-mechanism), is reduced by applying the PSSA to $\ce{ES}$ together with an enzyme
  balance $[\ce{E}]_0 = [\ce{E}] + [\ce{ES}]$, [](#eq-enzyme-balance). This is the exact analogue of
  the site balance used for a catalyst surface.
- The result is the Michaelis–Menten rate law
  $r = r\un{max}[\ce{S}]/(K\un{M} + [\ce{S}])$, [](#eq-michaelis-menten), with
  $r\un{max} = k\un{2}[\ce{E}]_0$, [](#eq-rmax), and
  $K\un{M} = (k\un{1-} + k\un{2})/k\un{1+}$, [](#eq-km). The whole rate collapses onto two lumped
  parameters.
- At low substrate the reaction is first order in $\ce{S}$ with slope $r\un{max}/K\un{M}$; at high
  substrate it saturates at $r\un{max}$ and is zeroth order, because essentially all enzyme is bound
  as $\ce{ES}$.
- $K\un{M}$ has units of concentration and equals the substrate concentration at which the enzyme
  runs at half its maximum rate, [](#eq-km-half-max), which is how it is read off experimentally.
