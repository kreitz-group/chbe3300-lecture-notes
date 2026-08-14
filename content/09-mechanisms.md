---
title: Reaction Mechanisms
short_title: Mechanisms
label: ch-mechanisms
---

<!-- LaTeX source: Mechanisms.tex -->
<!-- Porting notes: mhchem is math-only here, so \ce{} in prose needs $...$; and \un{} expands to
     _{\textrm{...}}, so species subscripts must be written _{\ce{H2}}, never \un{\ce{H2}}.
     Radicals use \ce{Cl\bullet} and energized molecules \ce{A\ast}, both on the math axis.
     The source colour-codes the H2/O2 steps as branching, propagating and terminating, then
     repeats the same reactions three more times regrouped by colour. Colour alone is not
     accessible, so the 19 steps appear once, in a table with an explicit Role column. -->

:::{admonition} Learning objectives
:class: tip
After completing this chapter, you should be able to:

- Distinguish an overall (stoichiometric) reaction from the mechanism of elementary steps that
  underlies it, and apply the law of mass action to each elementary step.
- Explain why a unimolecular reaction requires bimolecular activation, and derive the
  Lindemann–Hinshelwood rate law and its high- and low-pressure limits.
- Define a radical and classify the steps of a radical chain mechanism as initiation, propagation,
  or termination.
- Identify the common propagation reaction types — abstraction, addition to a double bond,
  $\beta$-scission — and write their reverse reactions.
- Describe the catalytic chlorine cycle that destroys stratospheric ozone and explain why a single
  radical destroys many ozone molecules.
- Classify the steps of the $\ce{H2}$/$\ce{O2}$ mechanism as branching, propagating, or
  terminating, and use the competition between them to explain the shape of the explosion-limit
  (Z) curve.
- Explain the difference between high- and low-temperature combustion chemistry and the role of the
  alkylperoxy radical in auto-ignition.
:::

We have discussed reaction mechanisms here and there, especially when we talked about the PSSA and
the QEA. Let us revisit a few of the key concepts to make sure everyone is well versed in the
terminology going forward.

As a running example, take the production of phosgene, $\ce{COCl2}$. The overall reaction equation
is

$$
\ce{Cl2 + CO <=> COCl2} .
$$

This overall stoichiometric reaction is not how the reaction actually happens. Instead we have a set
of elementary reactions that make up the reaction mechanism:

$$
\begin{aligned}
\ce{Cl2} &\eqa{k\un{+1}}{k\un{-1}} \ce{2 Cl\bullet} \\
\ce{Cl\bullet + CO} &\eqa{k\un{+2}}{k\un{-2}} \ce{COCl\bullet} \\
\ce{COCl\bullet + Cl2} &\ce{->[$k_3$] COCl2 + Cl\bullet}
\end{aligned}
$$

We use the law of mass action to express the reaction rates.

:::{admonition} Discussion
:class: seealso
What was that again?
:::

The law of mass action states that, for an elementary step, the reaction orders equal the
stoichiometric coefficients,

$$
\alpha_i^{\text{fwd}} = \nu_i' , \qquad \alpha_i^{\text{rev}} = \nu_i'' ,
$$

so the rate expression becomes

$$
r = r\un{fwd} - r\un{rev}
  = k\un{fwd}\prod_i c_i^{\nu_i'} - k\un{rev}\prod_i c_i^{\nu_i''} .
$$

This approach is called **mass action kinetics**. The net rate of reaction is the difference between
forward and reverse reaction, which is also known as microscopic reversibility.

Mass action kinetics is always connected to equilibrium. At equilibrium the extent of reaction no
longer changes, $\mathrm{d}\xi/\mathrm{d}t = r = 0$, so the overall reaction rate is zero and the
forward and reverse rates are balanced:

$$
\begin{aligned}
r\un{fwd} &= r\un{rev} \\
k\un{fwd}\prod_i c_i^{\nu_i'} &= k\un{rev}\prod_i c_i^{\nu_i''} \\
\frac{k\un{fwd}}{k\un{rev}} &= \frac{\prod_i c_i^{\nu_i''}}{\prod_i c_i^{\nu_i'}} .
\end{aligned}
$$

With the net stoichiometric coefficient $\nu_i = \nu_i'' - \nu_i'$, this collapses to

$$
\frac{k\un{fwd}}{k\un{rev}} = \prod_i c_i^{\nu_i} = K\un{c} .
$$ (eq-mech-kc)

Writing a reaction as irreversible carries implications about the parameters of the mechanism, and
is mostly done for instructional purposes.

Those are the key assumptions worth recalling. What follows looks first at unimolecular reactions,
then at radical chain reactions in the context of atmospheric chemistry, explosions and combustion,
and polymerization. Towards the end we look at how mechanisms are constructed. We then move on to
heterogeneous catalysis — my favorite topic of this course.

## Unimolecular reactions

<!-- source: Mechanisms.tex L77 -->

Many gas-phase reactions follow first-order kinetics, such as the isomerization of cyclopropane to
propene,

$$
\ce{cyclo-C3H6 -> CH3CHCH2} , \qquad r = k\,[\ce{cyclo-C3H6}] .
$$

:::{admonition} Discussion
:class: seealso
Does anyone see what the issue could be with first-order kinetics? Think about our grain of sand in
the box.
:::

The molecule has to somehow gather enough energy to overcome the activation barrier and cross the
potential energy surface, and that energy can only be supplied through collision with other
molecules. Collision is a bimolecular step, so these unimolecular reactions have bimolecular
elementary steps. How does that fit together?

The mechanism for this activation was proposed and explained by Frederick Lindemann and Cyril
Hinshelwood, and is called the **Lindemann–Hinshelwood mechanism**:

$$
\ce{A + A} \eqa{k\un{+a}}{k\un{-a}} \ce{A\ast + A} ,
\qquad r\un{a} = k\un{+a}[\mathrm{A}]^2 - k\un{-a}[\ce{A\ast}][\ce{A}] .
$$

Here $\ce{A\ast}$ is an *energized* molecule, not a transition-state activated complex. The reverse
reaction is simply the energized molecule losing energy upon collision with another molecule. If the
energized molecule has enough energy, it decomposes or rearranges in a unimolecular step,

$$
\ce{A\ast ->[$k\un{d}$] P} , \qquad r\un{d} = k\un{d}[\ce{A\ast}] .
$$

We can formulate the net rate of formation of the energized molecule,

$$
\begin{aligned}
\frac{\mathrm{d}[\ce{A\ast}]}{\mathrm{d}t} &= r\un{a} - r\un{d} \\
\frac{\mathrm{d}[\ce{A\ast}]}{\mathrm{d}t}
&= k\un{+a}[\mathrm{A}]^2 - k\un{-a}[\ce{A\ast}][\ce{A}] - k\un{d}[\ce{A\ast}] .
\end{aligned}
$$

What can we invoke about our energized molecule? It is probably short-lived, so the PSSA applies:

$$
\begin{aligned}
\frac{\mathrm{d}[\ce{A\ast}]}{\mathrm{d}t} &\approx 0
= k\un{+a}[\mathrm{A}]^2 - k\un{-a}[\ce{A\ast}][\ce{A}] - k\un{d}[\ce{A\ast}] \\
[\ce{A\ast}] &= \frac{k\un{+a}[\mathrm{A}]^2}{k\un{d} + k\un{-a}[\mathrm{A}]} .
\end{aligned}
$$

We are interested in the rate law for formation of the product P,

$$
r\un{P} = \frac{\mathrm{d}[\mathrm{P}]}{\mathrm{d}t} = k\un{d}[\ce{A\ast}]
= \frac{k\un{d}k\un{+a}[\mathrm{A}]^2}{k\un{d} + k\un{-a}[\mathrm{A}]} .
$$ (eq-lindemann)

This is clearly not first order. Did we do something wrong? No — the math is correct. Let us consider
two limiting scenarios.

At **high pressure** of A, collisions between A molecules are very frequent — remember the collision
count for $\ce{N2}$ at ambient conditions — so we can assume $k\un{-a}[\ce{A}] \gg k\un{d}$, and the
$k\un{d}$ in the denominator is negligible:

$$
\frac{\mathrm{d}[\mathrm{P}]}{\mathrm{d}t}
= \frac{k\un{d}k\un{+a}}{k\un{-a}}[\ce{A}] = k\un{app}[\ce{A}] .
$$ (eq-lindemann-high-p)

At **very low pressure** of A, the unimolecular decay is faster than deactivation,
$k\un{-a}[\ce{A}] \ll k\un{d}$, and the kinetics switch to second order:

$$
\frac{\mathrm{d}[\mathrm{P}]}{\mathrm{d}t} = k\un{+a}[\ce{A}]^2 .
$$ (eq-lindemann-low-p)

## Radical reactions

<!-- source: Mechanisms.tex L133 -->

We have already talked about radicals and used them in several reaction mechanisms, but have not
actually explained what they are. A **radical** is a species that possesses an unpaired electron. We
use $\bullet$ to denote them — $\ce{Cl\bullet}$, $\ce{H\bullet}$, $\ce{OH\bullet}$, and so on.

:::{admonition} Discussion
:class: seealso
What about $\ce{O2}$?
:::

Oxygen has a triplet and a singlet electronic structure. The most stable form is the triplet, in
which there is an unpaired electron on each oxygen atom; the first excited state is the singlet.

Radicals are highly reactive and very short lived, yet they are crucial to many reaction mechanisms,
including

- atmospheric chemistry,
- combustion,
- cracking,
- polymer synthesis, and
- biological processes.

Although all radicals are reactive, they are not equally so. Comparatively, $\ce{O2}$ is quite
stable; $\ce{NO}$ and $\ce{NO2}$ persist for a long time; $\ce{HO2\bullet}$ is less reactive than
most other radicals; and $\ce{OH\bullet}$ is highly reactive.

Most reaction mechanisms involving radicals are so-called **radical chain reactions**. In these the
radical acts rather like a catalyst: it is consumed but regenerated, so it can undergo another
cycle. There are always three steps in a chain reaction — initiation, propagation, and termination.
We have already looked at the chain mechanism for the formation of HBr in
[](#ch-multiple-reactions). Let us look more closely at each step.

### Initiation

The first step is always the creation of radicals, mostly by homolytic cleavage of a sigma bond.

:::{admonition} Discussion
:class: seealso
What does *homolytic* mean?
:::

Homolytic cleavage means the bond breaks evenly, so each fragment keeps one electron. Breaking a
bond does not just happen on its own: we have to inject energy into the molecule, either through
radiation or through collision.

$$
\begin{aligned}
\ce{H2 + M &-> 2 H\bullet + M} \\
\ce{H2O2 + M &-> 2 OH\bullet + M} \\
\ce{Cl2 &->[$h\nu$] 2 Cl\bullet}
\end{aligned}
$$

These are dissociation reactions. In some special cases other reaction types can lead to initiation,

$$
\begin{aligned}
\ce{CH4 + O2 &-> CH3\bullet + HO2\bullet} \\
\ce{NH3 + O &-> NH2\bullet + OH\bullet} .
\end{aligned}
$$

### Propagation

Many reaction types exist for propagation, but two are especially important.

**Abstraction reactions** transfer an atom, usually H, from a closed-shell molecule to a radical,
regenerating a different radical:

$$
\begin{aligned}
\ce{Cl\bullet + CH4 &-> HCl + CH3\bullet} \\
\ce{OH\bullet + CH4 &-> H2O + CH3\bullet} \\
\ce{H\bullet + Cl2 &-> HCl + Cl\bullet}
\end{aligned}
$$

The reverse of an abstraction reaction is also an abstraction.

A second class is **addition to a double bond**,

$$
\ce{CH3\bullet + CH2CH2 -> CH3CH2CH2\bullet} ,
$$

whose reverse is a $\beta$-scission. Various positions in the molecule have names: $\alpha$,
$\beta$, $\gamma$, and so on. The atom that carries the radical is the $\alpha$ position, and
$\beta$ is the next neighbor. A $\beta$-scission breaks a bond at the $\beta$ position,

$$
\ce{CH3O\bullet -> H\bullet + CH2O} .
$$

### Termination

Most termination reactions are simply the reverse of an initiation reaction,

$$
\begin{aligned}
\ce{2 H\bullet + M &-> H2 + M} \\
\ce{2 Cl\bullet + M &-> Cl2 + M} .
\end{aligned}
$$

Other examples include wall reactions, or conversion to a more stable, less reactive radical.
Converting to a stable radical effectively terminates the chain, because the resulting radical no
longer propagates efficiently:

$$
\ce{R\bullet + wall -> product} , \qquad
\ce{H2O2 + H\bullet -> HO2\bullet + H2} .
$$

These steps are the ones we typically invoke the PSSA on in radical-chain mechanisms.

### Atmospheric chemistry

<!-- source: Mechanisms.tex L236 -->

Radicals are hugely important in atmospheric chemistry. Unfortunately, they can also play a very
harmful role.

Most of the Earth's ozone is in the stratosphere, which begins about 10–15 km above the surface. The
ozone is concentrated in a layer roughly 15–35 km high, which we call the ozone layer. How is this
ozone formed? Through natural processes in the stratosphere: $\ce{O2}$ decomposes under UV-C
irradiation to oxygen radicals,

$$
\ce{O2 ->[$h\nu$] 2 O\bullet} ,
$$

where UV-C light covers wavelengths in the range 100–280 nm. The oxygen radicals are highly reactive
and combine with diatomic oxygen, with a third body M carrying away the excess energy, to form
ozone:

$$
\ce{O\bullet + O2 + M -> O3 + M} .
$$

The production of ozone is balanced by its destruction. Ozone absorbs ultraviolet light — UV-B, in
the range 280–315 nm — and decomposes, which is what shields us from most of the harmful UV
radiation. This is a fine balance between production and destruction.

Closer to the ground, in the troposphere, ozone has a low abundance, roughly 100 ppb. Tropospheric
ozone can be formed through pollutants such as $\ce{NO_x}$, $\ce{CO}$, and various hydrocarbons, as
well as during thunderstorms, and it is also emitted by trees.

We have interfered with the stratospheric balance through chemicals released into the environment,
in particular the chlorofluorocarbons (CFCs).

:::{admonition} Discussion
:class: seealso
Who has heard of CFCs? Why did we start using these materials?
:::

They have been used in many processes and were once praised as an excellent chemical solution to
many problems, because they are mostly non-toxic, non-flammable, and stable. They can decompose
under UV radiation in the stratosphere, however:

$$
\ce{CCl3F ->[$h\nu$] CCl2F\bullet + Cl\bullet} .
$$

CFCs also contribute to the greenhouse effect. The chlorine or bromine radicals released this way
are bad for the ozone layer: they destroy ozone through a radical chain reaction,

$$
\begin{aligned}
\ce{Cl\bullet + O3 &-> ClO\bullet + O2} \\
\ce{ClO\bullet + O\bullet &-> Cl\bullet + O2} \\ \hline
\ce{O3 + O\bullet &-> 2 O2}
\end{aligned}
$$

The radicals act almost like a catalyst. One chlorine or bromine radical can destroy 100,000 ozone
molecules before it recombines with another radical. We will not derive the rate law using the PSSA
here; you can try this at home if you want.

CFCs have led to the destruction of the ozone layer and the formation of an ozone hole, increasing
the risk of skin cancer. This discovery led to the Montreal Protocol in 1987. Chemists identified
CFCs as the main cause of ozone-layer destruction by unraveling the detailed kinetics of the
chlorine cycle, which led to the phase-out of many CFCs. The ozone layer is now recovering.

### Oxidation of $\ce{H2}$

<!-- source: Mechanisms.tex L285 -->

Most of you have done this reaction at some point in high school or college. Definitely do not
attempt it at home. If we mix $\ce{H2}$ and $\ce{O2}$ and provide a spark or flame, the mixture can
explode. Explosions are extremely fast radical reactions. Whether the mixture explodes depends on
the temperature — it is stable at room temperature — and also on the pressure.

:::{admonition} Discussion
:class: seealso
How does the explosion limit depend on pressure?
:::

A simple expectation would be that the explosion limit, in temperature, decreases with increasing
pressure. Reality is more complex, however, and we observe a so-called **Z-curve**.

Let us see how a detailed understanding of the reaction mechanism explains that Z-curve. First we
need the mechanism. Fortunately this is one of the best-studied systems, and we know that the
relevant set of stable species and radicals is $\ce{H2}$, $\ce{O2}$, $\ce{H2O}$, $\ce{H2O2}$,
$\ce{H\bullet}$, $\ce{O\bullet}$, $\ce{OH\bullet}$, and $\ce{HO2\bullet}$.

The system has 19 reactions, which group into branching, propagating, and terminating steps.
Initiation reactions count as branching, since they increase the radical pool.

:::{table} The 19-step $\ce{H2}$/$\ce{O2}$ mechanism. Steps 1–9 are the core $\ce{H2}$, $\ce{O2}$, $\ce{H2O}$, $\ce{H\bullet}$, $\ce{O\bullet}$, $\ce{OH\bullet}$ mechanism; steps 10–19 are the $\ce{H2O2}$ and $\ce{HO2\bullet}$ submechanism. The role column replaces the colour coding used in the printed notes.
:label: tab-h2-o2-mechanism

| # | Reaction | Role |
|--:|:---------|:-----|
| 1  | $\ce{H2 + M <=> 2 H\bullet + M}$ | branching |
| 2  | $\ce{O2 + M <=> 2 O\bullet + M}$ | branching |
| 3  | $\ce{H2 + O2 <=> H\bullet + HO2\bullet}$ | branching |
| 4  | $\ce{H\bullet + O2 <=> OH\bullet + O\bullet}$ | branching |
| 5  | $\ce{O\bullet + H2 <=> OH\bullet + H\bullet}$ | branching |
| 6  | $\ce{H2 + OH\bullet <=> H\bullet + H2O}$ | propagating |
| 7  | $\ce{OH\bullet + OH\bullet <=> O\bullet + H2O}$ | terminating |
| 8  | $\ce{H\bullet + O\bullet + M <=> OH\bullet + M}$ | terminating |
| 9  | $\ce{H\bullet + OH\bullet + M <=> H2O + M}$ | terminating |
| 10 | $\ce{H\bullet + O2 + M <=> HO2\bullet + M}$ | the odd one out |
| 11 | $\ce{HO2\bullet + H\bullet <=> OH\bullet + OH\bullet}$ | propagating |
| 12 | $\ce{HO2\bullet + O\bullet <=> O2 + OH\bullet}$ | terminating |
| 13 | $\ce{HO2\bullet + OH\bullet <=> O2 + H2O}$ | terminating |
| 14 | $\ce{HO2\bullet + HO2\bullet <=> H2O2 + O2}$ | terminating |
| 15 | $\ce{H2O2 + M <=> 2 OH\bullet + M}$ | branching |
| 16 | $\ce{H2O2 + H\bullet <=> H2O + OH\bullet}$ | propagating |
| 17 | $\ce{H2O2 + H\bullet <=> HO2\bullet + H2}$ | propagating |
| 18 | $\ce{H2O2 + O\bullet <=> HO2\bullet + OH\bullet}$ | branching |
| 19 | $\ce{H2O2 + OH\bullet <=> HO2\bullet + H2O}$ | propagating |
:::

That is seven branching steps (1, 2, 3, 4, 5, 15, 18), six terminating steps (7, 8, 9, 12, 13, 14),
five propagating steps (6, 11, 16, 17, 19), and one that does not fit the pattern — step 10.

The most important reaction in the whole sequence is the branching step
$\ce{H\bullet + O2 -> OH\bullet + O\bullet}$, reaction 4. This is probably *the* best-studied
gas-phase reaction. How is the initial $\ce{H\bullet}$ produced? At low temperatures through
$\ce{H2 + O2}$, and at high temperatures through $\ce{H2 + M}$. Once started, the radical
concentration explodes quickly.

Now back to the Z-curve. Fix a temperature and look at what happens as we increase the pressure.

At **low pressures**, collisions between molecules are rare. Bimolecular branching has a
second-order pressure dependence, while wall collisions are comparatively faster. The radicals are
lost to the walls before they can branch, so the mixture does not explode.

As we **increase the pressure**, the bimolecular branching and propagation steps overtake wall loss,
thanks to their second-order pressure dependence. In this regime the sequence
$4 + 5 + 2\times 6$ converts one $\ce{H\bullet}$ into three $\ce{H\bullet}$ radicals — this is the
first explosion limit.

Increasing $p$ **further** activates the termolecular step
$\ce{H\bullet + O2 + M -> HO2\bullet + M}$, step 10. $\ce{HO2\bullet}$ is comparatively stable and
effectively terminates the chain, so the mixture stops exploding. This is the second limit.

Increasing $p$ **further still** leads to chain branching once again, this time through a peroxide
pathway, steps 15 and 18. Combining these regimes reproduces the full Z-curve.

### Combustion

<!-- source: Mechanisms.tex L376 -->

The $\ce{H2}$/$\ce{O2}$ mechanism is one of the simplest combustion mechanisms — really. Everything
else is much more complicated, and so is the next example. Let us talk about combustion, and more
specifically about the combustion process in the engine of a car.

:::{admonition} Discussion
:class: seealso
Who has a car? What kind of engine does it have? What is the research octane number (RON) of the
fuel you put in it? How is the octane number determined?
:::

In addition to diesel and spark-ignition gasoline engines, there is another engine type: the
low-temperature auto-ignition engine. It has very high efficiency and produces only low amounts of
$\ce{NO_x}$ and soot. We really need to understand the kinetics in order to design these engines.

:::{admonition} Discussion
:class: seealso
Which molecules are in your gasoline?
:::

Ignition is highly sensitive to the fuel composition. Developing a reaction mechanism for these
systems is hugely complicated: there are simply too many reactions to consider, and the number grows
combinatorially with the number of species.

Let us look in more detail at the sequence of steps in an auto-ignition engine. Spontaneous ignition
is undesirable because it leads to rough operation and can damage the engine. Why does the mixture
ignite when we compress it? Thermal and pressure effects. After the end of the compression stroke it
typically takes a short time before the mixture explodes; this delay is the **ignition delay**.

:::{admonition} Discussion
:class: seealso
How does the ignition delay change with temperature — does it increase or decrease?
:::

We might expect the delay to shrink with hotter mixtures. The observed trend is more nuanced.

We can split the ignition-delay diagram into a hot zone and a cold zone, corresponding to two
different reaction mechanisms. The hot zone is where diesel and gasoline engines operate.
High-temperature combustion chemistry is comparatively simple: the fuel breaks down quickly through
$\beta$-scission, and our $\ce{H2}$/$\ce{O2}$ chemistry determines the overall rate. Most fuels
share the same high-temperature characteristics, which is why one can, for example, run a diesel
engine on sunflower oil without much difference.

At low temperatures, fuels behave more individually, and the oxidation depends heavily on the
structure of the alkane. The key reaction intermediate here is the **alkylperoxy radical**. This
structural sensitivity is why RON depends on the alkane structure. Following one specific example,
propane, as a representative alkane:

1. The alkane reacts with $\ce{O2}$ to form an alkylperoxy radical
   ($\ce{R-O-O\bullet}$, written $\ce{ROO\bullet}$).
2. The alkylperoxy radical isomerizes to a hydroperoxyalkyl radical ($\ce{\bullet QOOH}$).
3. A second $\ce{O2}$ addition gives a hydroperoxyalkylperoxy radical ($\ce{\bullet OOQOOH}$).
4. This species rearranges into a keto-hydroperoxide ($\mathrm{OQ'OOH}$).
5. Decomposition of the keto-hydroperoxide gives a formyl-alkoxy radical
   ($\mathrm{OQ'O}\bullet$) and $\ce{OH\bullet}$, providing the radicals that eventually ignite the
   mixture.

## Polymerization

<!-- source: Mechanisms.tex L412 -->

:::{note} Covered in lecture
This section has no written notes in the LaTeX source — the material is delivered from slides. As a
side note, polymer synthesis also proceeds through radical reactions.
:::

## Constructing reaction mechanisms

<!-- source: Mechanisms.tex L417 -->

This is our last topic on homogeneous reaction mechanisms. We have looked at several different
systems: atmospheric chemistry, including the chlorine cycle and the formation of HBr; explosions,
in the $\ce{H2}$/$\ce{O2}$ system; and the combustion of fuels. The Rice–Herzfeld mechanism for
pyrolysis reactions came up in the homework.

In each of these cases we were handed a list of elementary reactions and told that this is how the
chemistry happens. But how do we actually know that the reaction proceeds this way and not some
other way? The answer is that we have advanced considerably in building reaction mechanisms from the
ground up, and in determining the kinetic and thermophysical parameters from quantum-chemical
computations. This has led to the area of **predictive kinetics**. Detailed mechanistic knowledge is
very powerful.

## Summary

<!-- source: Mechanisms.tex L429 -->

- An overall stoichiometric equation is a bookkeeping statement, not a description of what happens.
  The mechanism is the underlying set of elementary steps, and only for those does the law of mass
  action apply, with $\alpha_i^{\text{fwd}} = \nu_i'$ and $\alpha_i^{\text{rev}} = \nu_i''$. At
  equilibrium the same kinetics gives back $k\un{fwd}/k\un{rev} = K\un{c}$, [](#eq-mech-kc).
- Unimolecular reactions are first order only in appearance: the molecule must be energized by
  collision first. The Lindemann–Hinshelwood mechanism plus the PSSA on the energized molecule
  $\ce{A\ast}$ gives
  $\mathrm{d}[\mathrm{P}]/\mathrm{d}t = k\un{d}k\un{+a}[\mathrm{A}]^2/(k\un{d} + k\un{-a}[\mathrm{A}])$,
  [](#eq-lindemann), which is first order at high pressure, [](#eq-lindemann-high-p), and second
  order at low pressure, [](#eq-lindemann-low-p).
- A radical carries an unpaired electron. Radical chain mechanisms always consist of initiation —
  usually homolytic bond cleavage, requiring radiation or collision — propagation (abstraction,
  addition to a double bond, $\beta$-scission), and termination (typically the reverse of
  initiation, wall reactions, or conversion to an unreactive radical).
- In the stratosphere, $\ce{O2}$ photolysis and $\ce{O\bullet + O2 + M}$ form the ozone layer.
  Chlorine radicals released from CFCs destroy ozone catalytically — one radical accounts for
  roughly 100,000 ozone molecules — which is what motivated the Montreal Protocol.
- The $\ce{H2}$/$\ce{O2}$ system separates into branching, propagating, and terminating steps,
  [](#tab-h2-o2-mechanism). The competition between them explains the Z-shaped explosion limit: wall
  termination dominates at low $p$; bimolecular branching (steps $4 + 5 + 2\times 6$, which turns
  one $\ce{H\bullet}$ into three) wins at intermediate $p$; the termolecular step
  $\ce{H\bullet + O2 + M -> HO2\bullet + M}$ terminates the chain at higher $p$; and the peroxide
  pathway branches again at higher $p$ still.
- High-temperature combustion is dominated by $\beta$-scission and $\ce{H2}$/$\ce{O2}$ chemistry, so
  most fuels behave alike. Low-temperature auto-ignition proceeds through the alkylperoxy radical
  and is highly structure-sensitive, which is why the octane number depends on the shape of the
  alkane.
