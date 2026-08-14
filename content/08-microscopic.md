---
title: Microscopic Details of Chemical Reactions
short_title: Microscopic details
label: ch-microscopic
---

<!-- LaTeX source: Microscopic.tex -->
<!-- Porting notes: mhchem is math-only here, so \ce{} in prose needs $...$; and \un{} expands to
     _{\textrm{...}}, so species subscripts must be written _{\ce{N2}}, never \un{\ce{N2}}.
     Radicals use \ce{F\bullet}; \Tilde is not a KaTeX command, so number densities use \tilde.
     Nested directives need a longer outer fence (:::: around :::). -->

:::{admonition} Learning objectives
:class: tip
After completing this chapter, you should be able to:

- Explain why the rate constant of an elementary step is generally not directly measurable, and
  what is gained by knowing it.
- Derive the collision frequency of a bimolecular gas-phase reaction from kinetic gas theory, and
  compute mean and relative molecular speeds from the Maxwell–Boltzmann distribution.
- State the collision-theory rate constant, identify the energetic and steric factors, and explain
  why hard-sphere collision theory fails for reactions with strict orientation requirements.
- Count the degrees of freedom of a molecule and identify the dimensionality of its potential
  energy surface.
- Characterize minima and first-order saddle points by the sign of the energy derivatives, and
  identify the reaction coordinate.
- Explain why chemical reactions are rare events on the timescale of molecular vibrations.
- Derive the Eyring equation from the quasi-equilibrium assumption between reactants and the
  transition state, and evaluate $K^\ddagger$ from molecular partition functions.
- Recast transition state theory in thermodynamic form and relate $E\un{a}$ and $A$ to
  $\Delta\un{act}H^\ddagger$ and $\Delta\un{act}S^\ddagger$.
:::

In this chapter we shift gears and look at the fundamental level on which chemical reactions
proceed. We have already talked a little about reaction mechanisms, but not very rigorously. Most
of the reactions and examples we have dealt with were stoichiometric reaction equations, such as
the methanation reaction $\ce{CO + 3 H2 -> CH4 + H2O}$. These simply summarize the stoichiometry
that we observe experimentally, and each has an underlying sequence of elementary reactions — the
reaction mechanism. We have seen that we can simplify such a mechanism using the PSSA or the QEA to
derive global rate expressions that we can fit to experimental data. It is not possible, however,
to measure the rate constants of individual elementary reactions directly in the general case; a
few reactions are accessible with very specialized experimental setups.

:::{admonition} Discussion
:class: seealso
Do we actually need to know all the individual rate constants?
:::

It depends what you want to do. For designing a reactor based on a reaction that works in the
laboratory, typically not — you can often fit a simple power-law equation. A reasonable
understanding of the mechanism helps to develop suitable rate expressions that can be fitted to
experiments, and these work very well for catalysis.

Do not think of these elementary steps as happening in sequence. Rather, they all happen at the
same time.

So what is the appeal of knowing all the individual rate constants? You have achieved a full
mechanistic understanding of the reaction system. Without any fits, you can describe how the
reaction will behave at any given condition, and you can design the reactor or catalyst to work at
the limit of what is possible. To figure this out, physical chemists have turned to developing
theoretical methods to elucidate the parameters. We will thus leave the realm of chemical
engineering and follow our chemistry colleagues on this quest to understand rate constants. This
goes into the weeds of physical chemistry.

There are two main theories for determining rate constants: collision theory and transition state
theory. We start with collision theory.

## Collision theory

<!-- source: Microscopic.tex L49 -->

Collision theory is the simplest theory for deriving the rate constant of an elementary reaction.
For a bimolecular reaction to occur, two reactant molecules must collide:

$$
\ce{A + B -> C} , \qquad r = k\, c\un{A} c\un{B} ,
$$

with the Arrhenius rate constant

$$
k = A \exp\left(\frac{-E\un{a}}{RT}\right) .
$$

Collision theory relies on kinetic gas theory.

:::{admonition} Discussion
:class: seealso
Who has heard of kinetic gas theory before?
:::

We will not have time to go into depth on the underlying methods, but we will provide enough
explanation to follow the derivation. Kinetic gas theory makes a few key assumptions that
essentially treat molecules as an ideal gas:

- Hard spheres with elastic collisions, conserving momentum and kinetic energy.
- Particles are very small, with volume $\approx 0$.
- No attractive or repulsive interactions.

A chemical reaction only occurs if the molecules collide, so we need to figure out how often
molecules collide — the collision frequency — and we want a theory that predicts it from the
properties of the colliding molecules. The following thought experiment provides a route. Imagine
that all molecules are stationary except one, which moves in a straight line and sweeps out a
cylinder, colliding with every molecule whose center lies within its collision diameter. A collision
happens whenever two centers come closer than $d\un{AB}$, defined as

$$
d\un{AB} = \frac{1}{2}\left(d\un{A} + d\un{B}\right) ,
$$

where $d\un{A}$ is the collision diameter of A and $d\un{B}$ that of B.

:::{figure} ../figures/Cylinder.png
:label: fig-cylinder
:alt: Hand-drawn sketch of the swept collision cylinder. Two horizontal lines mark the top and bottom of the cylinder, with a small blob labelled d_A moving along inside it and a larger blob labelled d_B straddling the upper boundary. A vertical arrow at the right marks the cylinder radius, annotated as d_A plus d_B over 2 equals r. A brace under the cylinder gives its length as delta t times the mean relative speed.
:width: 55%

Schematic of the swept-out collision cylinder. The radius of the cylinder is $d\un{AB}$, the sum of
the molecular radii.
:::

From this we define a collision cross section,

$$
\sigma\un{AB} = \pi d\un{AB}^2 .
$$

The molecule travels with an average relative speed $\langle v\un{rel}\rangle$, so the volume of the
cylinder swept during a time interval $\Delta t$ is

$$
V\un{c} = \sigma\un{AB}\langle v\un{rel}\rangle \Delta t .
$$

During this time, the number of collisions $N\un{AB}$ that one A molecule experiences with B depends
on the number density of species B, $\tilde{c}\un{B}$:

$$
\begin{aligned}
N\un{AB} &= V\un{c}\,\tilde{c}\un{B}
 = \sigma\un{AB}\langle v\un{rel}\rangle \Delta t\, \tilde{c}\un{B} \\
z\un{AB} &= \sigma\un{AB}\langle v\un{rel}\rangle \tilde{c}\un{B} \qquad \text{(per unit time)} .
\end{aligned}
$$

Here $z\un{AB}$ is in s$^{-1}$ and $\tilde{c}\un{B}$ is in molecules per cubic meter. We must also
count collisions of B molecules with A, which add to the total, so the overall collision frequency
per unit volume is

$$
Z\un{AB} = \sigma\un{AB}\langle v\un{rel}\rangle\, \tilde{c}\un{A}\, \tilde{c}\un{B} .
$$

The number density relates to the molar concentration by $\tilde{c}\un{A} = c\un{A} N\un{A}$, where
$N\un{A}$ is the Avogadro constant.

:::{admonition} Discussion
:class: seealso
What do we mean by the relative speed of the molecule? How does it relate to the average speeds of
A and B individually?
:::

We have a very large ensemble of molecules. At a concentration of 1 mol m$^{-3}$ we have
$6.023\times10^{23}$ molecules per cubic meter. Not all of them move at the same speed, however,
just as people do not all walk at the same pace. There is a distribution of molecular speeds,
described by the **Maxwell–Boltzmann distribution**,

$$
f(v) = 4\pi \left(\frac{m}{2\pi k\un{B}T}\right)^{\frac{3}{2}} v^2
\exp\left(\frac{-mv^2}{2k\un{B}T}\right) ,
$$ (eq-maxwell-boltzmann)

where $k\un{B}$ is the Boltzmann constant and $m$ the molecular mass. The Boltzmann constant is
$k\un{B} = 1.381\times10^{-23}\ \mathrm{J\,K^{-1}}$, and $f(v)$ is the probability density for
finding a molecule with speed $v$ at temperature $T$.

:::{admonition} Discussion
:class: seealso
How does the Boltzmann constant relate to the ideal gas constant?
:::

The connection is $R = k\un{B} N\un{A}$.

:::{figure} ../figures/Boltzmann.png
:label: fig-boltzmann
:alt: Hand-drawn plot of the speed distribution f of v against speed v, showing three skewed curves. The leftmost and tallest, with the narrowest spread, is labelled low temperature and high mass. A middle curve is labelled intermediate. The rightmost, lowest and broadest curve is labelled high temperature and low mass. All three start at the origin and decay to zero at large v.
:width: 55%

Maxwell–Boltzmann speed distribution $f(v)$ for different temperatures and molecular masses.
Increasing $T$ broadens and shifts the distribution to higher $v$; increasing $m$ does the opposite.
:::

To determine the average speed of species A, we compute the first moment of the distribution,

$$
\langle v\rangle = \int_0^{\infty} v\, f(v)\, \mathrm{d}v
= \left(\frac{8k\un{B}T}{\pi m}\right)^{1/2} .
$$

Since two species are reacting, we need the average relative speed of A and B. Kinetic gas theory
together with the Maxwell–Boltzmann distribution yields

$$
\langle v\un{rel}\rangle = \left(\frac{8k\un{B}T}{\pi m\un{AB}}\right)^{1/2} ,
$$ (eq-v-rel)

where $m\un{AB}$ is the reduced mass,

$$
m\un{AB} = \frac{m\un{A} m\un{B}}{m\un{A} + m\un{B}} ,
$$

which reduces a two-body problem to a one-body problem.

**Example: average speed of $\ce{N2}$ molecules.** Let us calculate the average speed of $\ce{N2}$
molecules at 25 °C, or 298 K. A nitrogen molecule weighs 28.02 u, converted to kg by multiplying by
$1.661\times10^{-27}\ \mathrm{kg\,u^{-1}}$. For two identical molecules the reduced mass is
$m\un{AA} = \tfrac{1}{2}m\un{A}$, so

$$
\langle v\un{rel}\rangle = \left(
\frac{8 \cdot 1.381\times10^{-23}\ \mathrm{J\,K^{-1}} \cdot 298\ \mathrm{K}}
{\pi \cdot 14.01\ \mathrm{u} \cdot 1.661\times10^{-27}\ \mathrm{kg\,u^{-1}}}
\right)^{1/2} = 671\ \mathrm{m\,s^{-1}} .
$$

:::{admonition} Discussion
:class: seealso
What is the speed of sound in air at sea level? Compare it to $\langle v\un{rel}\rangle$ above.
:::

For reference, the speed of sound in air at sea level is approximately 343 m s$^{-1}$.

Returning to the collision frequency, we substitute $\langle v\un{rel}\rangle$ and obtain the total
collision frequency per unit volume,

$$
Z\un{AB} = \sigma\un{AB} \left(\frac{8k\un{B}T}{\pi m\un{AB}}\right)^{1/2}
c\un{A} c\un{B} N\un{A}^2 .
$$ (eq-collision-frequency)

**Example: collision frequency of $\ce{N2}$ molecules.** Consider the same temperature, 298 K, and a
pressure of 1 bar, at which the molar concentration of $\ce{N2}$ is 40 mol m$^{-3}$. The collision
cross section is $\sigma = 0.43\ \mathrm{nm^2} = 0.43\times10^{-18}\ \mathrm{m^2}$. For collisions
of identical molecules the formula carries an extra factor of $\tfrac{1}{2}$ to avoid
double-counting indistinguishable pairs, which combines with the reduced mass $m\un{AA} = m/2$ to
give

$$
Z_{\ce{N2},\ce{N2}} = \sigma \left(\frac{4k\un{B}T}{\pi m_{\ce{N2}}}\right)^{1/2}
c_{\ce{N2}}^2 N\un{A}^2 = 8.4\times10^{34}\ \mathrm{m^{-3}\,s^{-1}} ,
$$

which is enormous. Even in 1 cm³ that is $8.4\times10^{28}$ collisions per second.

Comparing with the rate expression from the start of the chapter, the molar second-order rate
constant predicted by collision theory is

$$
k = N\un{A}\, \sigma \left(\frac{8k\un{B}T}{\pi m\un{AB}}\right)^{1/2} .
$$

The Avogadro factor converts the per-molecule formula into a molar rate constant, with units of
m³ mol$^{-1}$ s$^{-1}$, or conventionally L mol$^{-1}$ s$^{-1}$. We see that $k \sim \sqrt{T}$,
which does *not* have the Arrhenius temperature dependence we expect from rate constants in
practice.

:::{admonition} Discussion
:class: seealso
What did we neglect that would give us the missing exponential temperature dependence?
:::

:::{figure} ../figures/Collision.png
:label: fig-collision-reactive-nonreactive
:alt: Two hand-drawn cartoons. On the left, labelled non-reactive, arrows for A and BC converge on a jagged POW burst and emerge again unchanged as A and BC. On the right, labelled reactive, arrows for A and BC converge on a similar burst but emerge as the rearranged products AB and C.
:width: 55%

Non-reactive (left) and reactive (right) collisions. A reaction requires the kinetic energy along
the line of impact to exceed a threshold $E\un{c}$.
:::

We have neglected the energy distribution of the molecules. Not every collision carries enough
energy to lead to a chemical reaction; only those above a threshold $E\un{c}$ are reactive. The
energy distribution is also Boltzmann-distributed. We jump straight to the result, since the full
derivation is involved; the underlying assumption is that the collision cross section depends on
the collision energy.

$$
k = N\un{A}\,
\underbrace{\sigma \left(\frac{8k\un{B}T}{\pi m\un{AB}}\right)^{1/2}}_{\text{collision frequency factor}}
\underbrace{\exp\left(\frac{-E\un{c}}{k\un{B}T}\right)}_{\text{energetic requirement}}
\cdot \underbrace{P}_{\text{steric factor}} .
$$ (eq-collision-theory-k)

We see that $k \sim \sqrt{T}\exp\left(-E\un{c}/k\un{B}T\right)$, which is similar in form to the
Arrhenius expression.

:::{table} Comparison of rate constants predicted by collision theory with experimental values. Collision theory captures the right order of magnitude for radical reactions but is far off for the $\ce{ClO}$ self-reaction, where steric requirements matter.
:label: tab-examples-collision-theory

| Reaction | Theory (L mol$^{-1}$ s$^{-1}$) | Experiment (L mol$^{-1}$ s$^{-1}$) |
|:---------|-------------------------------:|-----------------------------------:|
| $\ce{H + C2H6 -> C2H5 + H2}$ | $37.3\times10^{10}$ | $9.6\times10^{10}$ |
| $\ce{2 ClO -> Cl2 + O2}$     | $2.5\times10^{10}$  | $6.3\times10^{7}$  |
:::

The steric factor $P$ accounts for the fact that only certain collision orientations lead to
reaction. Head-on and glancing collisions of atoms or polyatomic molecules differ in their
reactivity, for instance. The theory works reasonably well for radical reactions — low energy
threshold, modest steric requirements — and is theoretically exact at very low pressures. In other
cases the rate constant can be off by many orders of magnitude, because the theory has substantial
limitations. First, we must know the energy threshold and the steric factor, both of which require
experimental input or more sophisticated theories. A more fundamental omission is the neglect of
attractive and repulsive interactions between the molecules.

:::{figure} ../figures/LJP.png
:label: fig-lennard-jones
:alt: Two hand-drawn plots of energy against the A-to-B separation. On the left, labelled potential energy diagram, the curve falls steeply from high energy at short separation, passes through a well below zero at intermediate separation, and rises back to zero at large separation. On the right, labelled hard sphere, the energy is a vertical wall at the separation d_AB and exactly zero everywhere beyond it, with no well.
:width: 55%

Lennard–Jones (or Morse) interaction potential between two molecules. The well at intermediate
distance reflects attractive interactions; the steep rise at short distance reflects repulsive
interactions. Hard-sphere collision theory ignores both features.
:::

[](#fig-lennard-jones) shows the molecule–molecule interaction potential, often approximated as a
Lennard–Jones or Morse potential. There are more advanced theories that incorporate these
interactions, but they go well beyond undergraduate kinetics. We will instead turn to transition
state theory. First, we need a brief detour into potential energy surfaces, which we have so far
only mentioned in passing.

## Potential energy surface

<!-- source: Microscopic.tex L265 -->

We now look more closely at the potential energy diagram and explain what it actually represents.
To do this, we first need the concept of a harmonic oscillator.

:::{admonition} Discussion
:class: seealso
What is a harmonic oscillator? You should have encountered this in introductory physics.
:::

The physical picture is a mass connected to a spring.

:::{figure} ../figures/Hookes_law.png
:label: fig-hookes-law
:alt: Diagram of a mass on a horizontal spring attached to a wall, drawn in three stacked states. In the top state the spring is stretched and a red arrow labelled F points right. In the middle state the mass sits at the equilibrium position marked 0 on the x axis. In the bottom state the spring is compressed and a red arrow points left, indicating the restoring force always acts back toward equilibrium.
:width: 45%

Mass on a spring obeying Hooke's law. Perturbing the mass from its equilibrium position $x_0$
generates a restoring force $F = -k(x - x_0)$, leading to oscillatory motion at a single
characteristic frequency.
:::

What happens if we perturb this system and push the mass to the right? When released, the spring
pulls it back. Kinetic and potential energy, the latter stored in the spring, are continually
exchanged. This motion has a characteristic frequency set by the spring constant and the mass. We
avoid the full mathematical treatment here and trust that you have seen it before; it is also not
strictly necessary for what follows.

The key point is that molecules behave very much like this mass-spring system: each bond can be
approximated as a harmonic oscillator.

:::{figure} ../figures/harmonic_oscillator_h2.png
:label: fig-h2-ho
:alt: Slide showing an H-H molecule drawn as two blue spheres joined by a bond, with an arrow to a mass on a spring at its equilibrium position. Below, a plot of potential energy against H-H distance shows a blue anharmonic curve that rises steeply at short distance and flattens at long distance, overlaid with a red parabola. An inset magnifies the region near the minimum, where the two curves are almost indistinguishable, annotated that for small displacements the harmonic approximation is quite good.
:width: 70%

Harmonic-oscillator approximation for the $\ce{H2}$ bond. The energy is parabolic in the bond
displacement around the equilibrium bond length.
:::

There is an equilibrium bond length; if we perturb the bond, the restoring force pulls it back
toward equilibrium.

A brief comment on the status of this model. Undergraduate courses introduce many models that are
merely stepping stones to more complete theories. The harmonic oscillator is not one of those.
Entire advanced textbooks are written on how the model can be used in kinetics, and it is used
routinely in research. You have to go quite deep into edge cases before a more complex anharmonic
model is needed.

With the harmonic oscillator we can now build a potential energy surface. Just as energy is stored
in the spring, energy is stored in each bond, and the surface tells us how that energy varies with
molecular geometry.

Let us look at the $\ce{H2}$ molecule.

:::{admonition} Discussion
:class: seealso
What is a degree of freedom?
:::

Each atom can move in three orthogonal directions in space, giving three degrees of freedom (DOF)
per atom. For an $N$-atom molecule there are thus $3N$ DOF in total; for $\ce{H2}$, six.

Does the potential energy surface also have $3N$ dimensions? Not quite — we have to subtract the
collective motions that do not change the molecular shape. Translation of the whole molecule and
rotation about its center of mass do not change the potential energy. Only changes in the H–H
distance do, and the H–H distance is an internal degree of freedom.

:::{admonition} Discussion
:class: seealso
What kind of motion is described by changes in the H–H distance?
:::

The answer is vibration, and this internal degree of freedom gives one vibrational mode. For
$\ce{H2}$, then, we have a one-dimensional potential energy surface.

Now consider a slightly more complex molecule, $\ce{H2O}$, in vacuum.

:::{admonition} Discussion
:class: seealso
$\ce{H2O}$ has three atoms, so how many DOF? Of these, how many are translational, rotational, and
internal?
:::

$\ce{H2O}$ has $3N = 9$ DOF: three translational and three rotational, since it is a nonlinear
molecule, leaving $3N - 6 = 3$ internal degrees of freedom. We can choose these three internal
coordinates as

- the O–H$_1$ distance,
- the O–H$_2$ distance,
- the H–O–H angle.

We can only visualize them one at a time, by holding the other two fixed. These internal coordinates
correspond to three vibrational modes: symmetric O–H stretch, asymmetric O–H stretch, and bending.

:::{figure} ../figures/triatomic_pes.png
:label: fig-h2o-pes
:alt: Slide with three plots of potential energy, one against the A-B distance, one against the B-C distance, and one against the A-B-C angle. Each shows a blue curve with a minimum, overlaid near the bottom of the well with a red dashed parabola representing the harmonic approximation; the angle plot additionally shows a second maximum and a further descent. A note in red observes that a three-dimensional surface is already very hard to visualize.
:width: 70%

Schematic slices of the $\ce{H2O}$ potential energy surface, one internal coordinate at a time. The
full surface is three-dimensional.
:::

The full potential energy surface can become very complicated, with one dimension for each internal
degree of freedom. Methane, for example, has five atoms, giving $3N - 6 = 9$ internal dimensions.
What we care about on the surface are small perturbations around equilibrium.

Let us look at a general one-dimensional potential energy surface. The key features are the
stationary points: the minima and the maxima.

:::{figure} ../figures/stationary_points_1.png
:label: fig-1d-pes
:alt: Slide showing a wavy one-dimensional curve with two minima separated by a maximum. Red starburst markers sit at the two minima and the maximum, labelled stationary points, where the annotation notes the force, minus the gradient of U, is zero. Two orange markers on the rising and falling flanks are labelled inflection points.
:width: 70%

Generic one-dimensional potential energy surface highlighting the stationary points: minima, which
are stable structures, and the local maximum, which is the transition state.
:::

The minima are the bottoms of wells. Any small perturbation away from a minimum is restored back to
the equilibrium structure. At a minimum,

$$
\begin{aligned}
\frac{\partial E}{\partial r_i} &= 0 \qquad \text{for all } i \\
\frac{\partial^2 E}{\partial r_i^2} &> 0 \qquad \text{for all } i .
\end{aligned}
$$ (eq-minimum-conditions)

:::{admonition} Discussion
:class: seealso
What is the physical meaning of the derivative of the energy with respect to position?
:::

The answer is the force: at a stationary point the force on every atom vanishes. The Hessian, the
matrix of second derivatives, is positive definite.

The $r_i$ are our internal degrees of freedom, and minima correspond to stable molecules. In two
dimensions we can construct a contour plot of the surface. A useful analogy is a topographic map of
the kind used for hiking, where contour lines connect points of equal elevation; here they connect
points of equal potential energy.

:::{figure} ../figures/minima_contour_1.png
:label: fig-2d-contour-1
:alt: Slide showing a two-dimensional contour plot with axes r_i and r_j. Concentric blue closed contours shrink toward a single red starburst marker at the centre, marking one minimum. A note in red suggests thinking of topological maps used for hiking.
:width: 70%

Two-dimensional contour plot of a potential energy surface with a single minimum. Contour lines join
points of equal energy; closely spaced contours indicate a steep wall.
:::

The other stationary points we care about are first-order saddle points: maxima along one coordinate
but minima along all others.

:::{figure} ../figures/pes_minima.png
:label: fig-1d-minima
:alt: Slide headed "we associate minima with molecules", showing a single blue parabola-like well with a red starburst at its base, annotated a stable molecule. Beside it the conditions are written out: the first derivative of U with respect to every coordinate is zero, and the second derivative is positive for every coordinate, so any small perturbation returns to the equilibrium structure.
:width: 70%

A minimum on a one-dimensional surface: the first derivative vanishes and every second derivative is
positive, so the structure is stable.
:::

At a first-order saddle point,

$$
\begin{aligned}
\frac{\partial E}{\partial r_i} &= 0 \qquad \text{for all } i \\
\frac{\partial^2 E}{\partial r_j^2} &< 0 \qquad \text{for exactly one } j \\
\frac{\partial^2 E}{\partial r_i^2} &> 0 \qquad \text{for all } i \neq j .
\end{aligned}
$$ (eq-saddle-conditions)

The specific degree of freedom $j$ along which the surface is curved downward is the **reaction
coordinate**.

:::{figure} ../figures/pes_saddlepoint.png
:label: fig-1d-saddle
:alt: Slide headed that the minimum energy path also contains a stationary point. A blue one-dimensional curve runs from a deep well up over a rounded maximum, marked with a red starburst, and down into a second shallower well. Beside it the saddle-point conditions are written out: the first derivative vanishes for every coordinate, the second derivative is negative for exactly one coordinate j, and positive for all the others.
:width: 70%

First-order saddle point: a maximum along the reaction coordinate but a minimum along all orthogonal
directions.
:::

We can also examine the relative depth and width of the wells, whose features map onto thermodynamic
quantities:

- **Well depth:** a deeper well is lower in enthalpy.
- **Well width:** a wider well is higher in entropy, since more vibrational states are accessible.

At equilibrium every system minimizes its Gibbs free energy, $\Delta G = \Delta H - T\Delta S$. At
low temperature enthalpy dominates; at high temperature the entropic term $-T\Delta S$ becomes
important. A few illustrative cases:

- **Same depth, same width.** For an ensemble of molecules we expect a 50/50 split at any
  temperature.
- **One deeper, same width.** At low temperature equilibrium heavily favors the deeper well. At
  higher temperatures some population is in the shallower well, but the deeper well always
  dominates.
- **Same depth, one wider.** At low temperature, a 50/50 split. As temperature increases the
  equilibrium shifts toward the wider, higher-entropy well.
- **One deeper, one wider.** At low temperature the population is concentrated in the deep well; at
  high temperature much of it shifts to the wide well, the split depending on relative depth and
  width.
- **Different barriers.** Two wells of similar depth, separated by barriers of different heights.
  Starting from all population in the higher well, both systems eventually equilibrate toward the
  lower well, but the system with the lower barrier gets there first.

### Rare events

<!-- source: Microscopic.tex L435 -->

We now want to understand how fast molecules go from one well into the other, which is what we mean
by a chemical reaction. Before we begin, two big-picture statements that may seem confusing at
first.

**Rate as events per time.**

$$
\text{rate} = \frac{\text{events}}{\text{time}} .
$$

This should be familiar by now. Equivalently, a rate is a number of molecules multiplied by the
frequency with which the event occurs in each one.

**Reactions are rare events.** Recall the mental picture of an $\ce{H2O}$ molecule as a small set of
mass-and-spring systems joined together. The molecule is constantly vibrating and jiggling.

- Most of the time, displacements from equilibrium are small.
- Infrequently, energy concentrates spontaneously into one mode, leading to a larger displacement.
- Very rarely, enough energy concentrates in one mode that the bond effectively breaks.

That last bullet is, in essence, a chemical reaction. Most molecular vibrations have wavenumbers in
the range 100–3000 cm$^{-1}$. Converted via the speed of light, these correspond to frequencies
between $3\times10^{12}$ and $9\times10^{13}\ \mathrm{s^{-1}}$, which is extremely fast. Reactions,
by contrast, are much slower. For a reaction occurring on a timescale of 1 s$^{-1}$ — already fast —
there are roughly 10 to 100 trillion vibrations per reactive event. For comparison, the odds of
winning the Powerball are roughly one in 300 million; reactive events are one in 100 trillion, six
orders of magnitude rarer. So yes, reactions truly are rare events at the molecular level. The
reason we still observe them in bulk is simply that we have $6.023\times10^{23}$ molecules per mole.

## Transition state theory

<!-- source: Microscopic.tex L466 -->

We now have a working understanding of the potential energy surface and how it reflects the
structure of the reactants. We will leverage this to derive a more comprehensive theory of the rate
constant. Consider the elementary reaction

$$
\ce{A + B ->[$k$] P} .
$$

Picture a 2D potential energy surface with a dividing surface separating reactants and products.
Returning to our hiking analogy: most hikers crossing a mountain ridge will pass over the lowest
point of the ridge, the saddle point. Only a few will take a longer, higher path. Molecules behave
similarly — most reactive trajectories pass through, or very near, the saddle point.

:::{figure} ../figures/06_2well_contour.png
:label: fig-2d-contour
:alt: Slide showing a two-dimensional contour plot with two sets of concentric blue contours, the left basin labelled R for reactant and the right labelled P for product, joined through a narrow neck where the contours pinch together at the saddle point. To the right, a small one-dimensional profile shows the same thing as two wells separated by a barrier, with a red arrow arcing from the reactant well over the barrier into the product well.
:width: 80%

Two-well contour plot of a 2D potential energy surface. Reactant and product basins are separated by
a dividing surface that passes through the first-order saddle point along the minimum-energy path.
:::

We can extract the 1D potential energy profile along the minimum-energy path connecting the two
basins, as in [](#fig-1d-pes). The reaction coordinate corresponds to one specific vibrational-like
mode of the system.

:::{figure} ../figures/Tst.png
:label: fig-tst
:alt: Hand-drawn energy profile along a reaction coordinate. A curve rises from a well labelled A plus BC over a single rounded barrier and down into a well labelled AB plus C. The top of the barrier is labelled activated complex, written as A, dots, B, dots, C, indicating the partially formed and partially broken bonds.
:width: 55%

Schematic of transition state theory: reactants in quasi-equilibrium with the transition state
$\ce{ABC^{\ddagger}}$, which decomposes irreversibly into products with frequency $\nu^\ddagger$
along the reaction coordinate.
:::

In transition state theory we assume there is a well-defined transition state connecting reactants
and products along the minimum-energy path on the reaction coordinate, and a dividing surface at the
saddle point that, once crossed, leads irreversibly to products — no recrossing. We can rewrite the
reaction as

$$
\ce{A + BC} \eqa{K^\ddagger}{} \ce{ABC^{\ddagger} ->[$k^\ddagger$] AB + C} ,
$$

where $K^\ddagger$ is the equilibrium constant for transition-state formation,

$$
K^\ddagger = \frac{[\mathrm{ABC}^{\ddagger}]\, c\un{0}}{[\mathrm{A}][\mathrm{BC}]} .
$$

The net rate of product formation is $r = k^\ddagger [\mathrm{ABC}^{\ddagger}]$.

The concentration of the transition state is not directly observable. We assume quasi-equilibrium
between reactants and the transition state, however, which lets us express
$[\mathrm{ABC}^\ddagger]$ in terms of reactant concentrations. For an ideal gas the standard
concentration is $c\un{0} = p\un{0}/RT$, so

$$
r = k^\ddagger K^\ddagger [\mathrm{A}][\mathrm{BC}]\frac{1}{c\un{0}}
  = k^\ddagger K^\ddagger [\mathrm{A}][\mathrm{BC}]\frac{RT}{p\un{0}} .
$$

Recalling that the reaction coordinate is one of the vibrational-like modes of the internal degrees
of freedom, this mode has a characteristic vibrational frequency $\nu^\ddagger$ along which the
transition state decomposes. We can therefore identify $k^\ddagger = \nu^\ddagger$, giving

$$
r = \underbrace{\nu^\ddagger K^\ddagger}_{k}[\mathrm{A}][\mathrm{BC}]\frac{1}{c\un{0}} .
$$

Invoking equilibrium for a system that is by definition far from equilibrium may seem inconsistent.
The justification is a separation of timescales: molecular relaxations within the reactant basin are
much faster than reactive crossings of the barrier. By the time a reactive trajectory reaches the
saddle point, the system has forgotten how it got there, and the local distribution at the
transition state is well approximated by an equilibrium distribution with the reactants. This is
what allows us to use an equilibrium constant to describe a kinetic process.

### Equilibrium constant from molecular partition functions

<!-- source: Microscopic.tex L539 -->

The next task is to evaluate $K^\ddagger$. We previously defined the general equilibrium constant as

$$
K^\ddagger = \frac{a_{\mathrm{ABC}^\ddagger}}{a\un{A}\, a\un{BC}} ,
$$

and it can also be written using molecular partition functions per unit volume,

$$
K^\ddagger = \frac{N\un{A}\, q_{\mathrm{ABC}^\ddagger}}{q\un{A} q\un{BC}} ,
$$

where the $q$ are molecular partition functions per unit volume and $N\un{A}$ is the Avogadro
constant. We will not work through the formal derivation here.

:::{admonition} Discussion
:class: seealso
What is a partition function, intuitively?
:::

Partition functions come from statistical thermodynamics. For an ensemble of molecules that can
occupy a set of discrete energy states, the partition function is essentially a count of accessible
states, weighted by their Boltzmann factors. The fraction of molecules in a given energy state is
given by the Boltzmann distribution, which we covered earlier in the semester.

Partition functions are generally defined with respect to the lowest-energy state of each species,
so we have to add an explicit term to the equilibrium constant accounting for the energy difference
between transition state and reactants at 0 K:

$$
K^\ddagger = \frac{N\un{A}\, q_{\mathrm{ABC}^\ddagger}}{q\un{A} q\un{BC}}
\exp\left(-\frac{\Delta E\un{0}^\ddagger}{RT}\right) ,
$$

where $\Delta E\un{0}^\ddagger$ is that 0 K energy difference. The total partition function
factorizes into contributions from translation, rotation, vibration, and electronic states,

$$
q = q\un{trans}\, q\un{rot}\, q\un{vib}\, q\un{el} ,
$$

and $q\un{el}$ can typically be set to 1. It is useful to express each contribution as a product of
factors per degree of freedom. For a nonlinear $N$-atom molecule,

$$
q\un{trans} = f\un{trans}^3 , \qquad
q\un{rot} = f\un{rot}^3 , \qquad
q\un{vib} = f\un{vib}^{3N-6} .
$$

For a linear molecule, replace $f\un{rot}^3$ by $f\un{rot}^2$ and $f\un{vib}^{3N-6}$ by
$f\un{vib}^{3N-5}$. Standard expressions for the partition functions are

$$
\begin{aligned}
q\un{trans} &= \left(\frac{2\pi m k\un{B}T}{h^2}\right)^{\frac{3}{2}}
&&\approx 10^{26}\ \text{to}\ 10^{28}\ \mathrm{L^{-1}} \\
q\un{rot} &= \frac{8\pi^2 I k\un{B}T}{\sigma\un{sym}h^2} &&\text{(linear)} \\
q\un{rot} &= \frac{8\pi^2(8\pi^3 I\un{1}I\un{2}I\un{3})^{1/2}(k\un{B}T)^{3/2}}{\sigma\un{sym}h^3}
&&\text{(nonlinear)}, \approx 10\ \text{to}\ 100 \\
q\un{vib} &= \prod_i \left(1 - \exp\left(-\frac{h\nu_i}{k\un{B}T}\right)\right)^{-1}
&&\approx 1\ \text{to}\ 10
\end{aligned}
$$

The Planck constant is $h = 6.626\times10^{-34}\ \mathrm{J\,s}$ and the Boltzmann constant is
$k\un{B} = 1.381\times10^{-23}\ \mathrm{J\,K^{-1}}$. The approximate ranges given are for
temperatures between 300 K and 500 K.

We now identify the vibrational mode $\nu^\ddagger$ of the reaction coordinate. This is a
large-amplitude motion, so its vibrational frequency is much lower than the other modes, and we
factor it out separately. For a single low-frequency mode,

$$
q\un{vib} = \frac{1}{1 - \exp\left(-h\nu^\ddagger/k\un{B}T\right)} ,
$$

and in the low-frequency limit $h\nu^\ddagger/k\un{B}T \ll 1$ we can use $\exp(-x) \approx 1 - x$,
giving

$$
q\un{vib} = \frac{k\un{B}T}{h\nu^\ddagger} .
$$

The partition function of the transition state is factored to separate this projected-out degree of
freedom,

$$
q_{\mathrm{AB}^\ddagger} = \left(q\un{vib}\right)\un{TST}\, q_{\mathrm{AB}^\ddagger}'
\quad \text{with} \quad
q_{\mathrm{AB}^\ddagger}' = q\un{rot}\, q\un{trans}\, q\un{vib}' ,
$$

where $q\un{vib}'$ contains only the remaining $3N-7$ vibrational modes, or $3N-6$ for a linear
transition state. Inserting these into the rate expression, the unknown frequency $\nu^\ddagger$
cancels:

$$
k = \nu^\ddagger \frac{k\un{B}T}{h\nu^\ddagger}\left(K^\ddagger\right)'\frac{RT}{p\un{0}}
  = \frac{k\un{B}T}{h}\left(K^\ddagger\right)'\frac{RT}{p\un{0}} .
$$ (eq-eyring)

This is the **Eyring equation**, the central result of transition state theory.

#### Example: TST rate constant for $\ce{F\bullet + H2 -> HF + H\bullet}$

We apply transition state theory to the gas-phase reaction
$\ce{F\bullet + H2 -> HF + H\bullet}$ at 300 K, starting from

$$
k = \frac{k\un{B}T}{h}\left(K^\ddagger\right)'\frac{RT}{p\un{0}} .
$$

**Translational partition functions.**

$$
\begin{aligned}
\frac{(q\un{trans})_{\mathrm{AB}^\ddagger}}{(q\un{trans})\un{A}\,(q\un{trans})\un{B}}
&= \frac{\left(\dfrac{2\pi m\un{TS} k\un{B}T}{h^2}\right)^{\frac{3}{2}}}
{\left(\dfrac{2\pi m\un{F} k\un{B}T}{h^2}\right)^{\frac{3}{2}}
 \left(\dfrac{2\pi m_{\ce{H2}} k\un{B}T}{h^2}\right)^{\frac{3}{2}}} \\
&= \frac{9.407\times10^{28}\ \mathrm{L^{-1}}}
{8.086\times10^{28}\ \mathrm{L^{-1}} \cdot 2.795\times10^{27}\ \mathrm{L^{-1}}}
= 4.162\times10^{-28}\ \mathrm{L} .
\end{aligned}
$$

The only inputs needed are the masses, available from the periodic table. Translational partition
functions are large.

**Rotational partition functions.** For a monatomic species $q\un{rot} = 1$. For a linear species,

$$
q\un{rot} = \frac{8\pi^2 I k\un{B}T}{\sigma\un{sym}h^2} .
$$

The moment of inertia is computed from the reduced mass and bond distance,

$$
I_{\ce{H2}} = \frac{M\un{H}M\un{H}}{M\un{H} + M\un{H}}\, d\un{H-H}^2
= 4.60\times10^{-48}\ \mathrm{kg\,m^2} ,
$$

using $d\un{H-H} = 0.74\ \text{Å}$ and the mass of the hydrogen atom. The symmetry number of
$\ce{H2}$ is $\sigma\un{sym} = 2$, since the molecule has two indistinguishable orientations around
the center of mass. The transition-state moment of inertia is taken from a quantum-chemical
calculation, $I\un{TS} = 1.234\times10^{-46}\ \mathrm{kg\,m^2}$, with $\sigma\un{sym} = 1$ (no
symmetry). The ratio of rotational partition functions is then

$$
\frac{(q\un{rot})_{\mathrm{AB}^\ddagger}}{(q\un{rot})\un{A}\,(q\un{rot})\un{B}}
= \frac{\dfrac{8\pi^2 I\un{TS}\, k\un{B}T}{\sigma\un{sym,TS}\, h^2}}
{1 \cdot \dfrac{8\pi^2 I_{\ce{H2}}\, k\un{B}T}{\sigma_{\mathrm{sym},\ce{H2}}\, h^2}}
= \frac{91.917}{1 \cdot 1.707} = 53.851 .
$$

**Vibrational partition functions.** The monatomic F has no vibrational modes. $\ce{H2}$ has a
stretching frequency of 4395 cm$^{-1}$, which gives a partition function that is essentially 1. The
linear three-atom transition state would normally have $3N - 5 = 4$ vibrational modes, but we
projected out one as the reaction coordinate, leaving three modes from quantum-chemical
calculations — 4007 cm$^{-1}$, 398 cm$^{-1}$, and 398 cm$^{-1}$ — giving
$(q\un{vib})_{\mathrm{AB}^\ddagger} = 1.378$. Hence

$$
\frac{(q\un{vib})_{\mathrm{AB}^\ddagger}}{(q\un{vib})\un{A}\,(q\un{vib})\un{B}}
= \frac{1.378}{1 \cdot 1} = 1.378 .
$$

**Assembly.** Putting it together,

$$
k\un{TST} = \frac{k\un{B}T}{h}\, N\un{A} \cdot 4.162\times10^{-28}\ \mathrm{L}
\cdot 53.851 \cdot 1.378
\exp\left(-\frac{\Delta E\un{0}^\ddagger}{RT}\right) .
$$

The barrier $\Delta E\un{0}^\ddagger$ comes from quantum-mechanical calculations; using
$\Delta E\un{0}^\ddagger = 6\ \mathrm{kJ\,mol^{-1}}$,

$$
k\un{TST} = 1.05\times10^{10}\ \mathrm{L\,mol^{-1}\,s^{-1}} , \qquad
k\un{exp} = 1.4\times10^{10}\ \mathrm{L\,mol^{-1}\,s^{-1}} .
$$

The experimental value from the NIST database is
$k = 1.4\times10^{10}\ \mathrm{L\,mol^{-1}\,s^{-1}}$, so the theory lands within about 25 % of
experiment — excellent agreement, given that the only empirical input was the barrier height.

Because it assumes that every trajectory crossing the dividing surface goes on to form products,
transition state theory is a strict upper bound on the rate constant *for a given potential energy
surface*: recrossing can only reduce the true rate. In practice a calculated rate constant can still
fall below experiment, as it does here, because the computed barrier $\Delta E\un{0}^\ddagger$
carries its own uncertainty, and because quantum tunnelling through the barrier — neglected entirely
in the treatment above — speeds up reactions involving light atoms such as H.

## Thermodynamic considerations

<!-- source: Microscopic.tex L672 -->

Transition state theory can also be cast in a different notation, one that exposes the thermodynamic
meaning of the activation parameters. The equilibrium constant for transition-state formation can be
written in terms of the Gibbs free energy of activation. In the derivation above we worked with a
potential energy surface; what we actually want is a free energy surface.

:::{admonition} Discussion
:class: seealso
What is the difference between a potential energy surface and a free energy surface? Hint: entropy.
:::

The rate constant depends directly on the change in free energy on going from reactants to the
transition state,

$$
k\un{TST} = \frac{k\un{B}T}{h}K^\ddagger
= \frac{k\un{B}T}{h}\exp\left(\frac{-\Delta\un{act}G^\ddagger}{RT}\right) .
$$

Splitting the free energy into entropic and enthalpic contributions via $G = H - TS$,

$$
\begin{aligned}
k\un{TST} &= \frac{k\un{B}T}{h}
\exp\left(\frac{-\left[\Delta\un{act}H^\ddagger - T\Delta\un{act}S^\ddagger\right]}{RT}\right) \\
k\un{TST} &= \frac{k\un{B}T}{h}
\exp\left(\frac{\Delta\un{act}S^\ddagger}{R}\right)
\exp\left(\frac{-\Delta\un{act}H^\ddagger}{RT}\right) .
\end{aligned}
$$ (eq-tst-thermo)

Suppose for a moment that we know the free energy of activation. One useful payoff of transition
state theory is that it gives a microscopic interpretation of the empirical Arrhenius expression.
Equating the two,

$$
\begin{aligned}
k\un{TST} &= k\un{Arrhenius} \\
\frac{k\un{B}T}{h}\exp\left(\frac{-\Delta\un{act}G^\ddagger}{RT}\right)
&= A\exp\left(\frac{-E\un{a}}{RT}\right) \\
\frac{k\un{B}T}{h}\exp\left(\frac{\Delta\un{act}S^\ddagger}{R}\right)
\exp\left(\frac{-\Delta\un{act}H^\ddagger}{RT}\right)
&= A\exp\left(\frac{-E\un{a}}{RT}\right) ,
\end{aligned}
$$

which implies, up to small temperature-dependent corrections,

$$
E\un{a} \approx \Delta\un{act}H^\ddagger ,
$$ (eq-tst-arrhenius-ea)

$$
A \approx \frac{k\un{B}T}{h}\exp\left(\frac{\Delta\un{act}S^\ddagger}{R}\right) .
$$ (eq-tst-arrhenius-a)

In other words, the activation energy reflects the change in enthalpy between reactants and
transition state, while the pre-exponential factor reflects the change in entropy. The entropy in
turn is derived from the partition functions, closing the loop with the Eyring derivation above.

## Summary

<!-- source: Microscopic.tex L713 -->

- Collision theory treats molecules as hard spheres and builds the rate constant from the collision
  frequency $Z\un{AB} = \sigma\un{AB}\langle v\un{rel}\rangle \tilde{c}\un{A}\tilde{c}\un{B}$,
  [](#eq-collision-frequency), with
  $\langle v\un{rel}\rangle = (8k\un{B}T/\pi m\un{AB})^{1/2}$, [](#eq-v-rel), evaluated from the
  Maxwell–Boltzmann distribution, [](#eq-maxwell-boltzmann), and the reduced mass $m\un{AB}$.
- Collisions alone give $k \sim \sqrt{T}$, not an Arrhenius form. Adding the Boltzmann fraction of
  collisions above a threshold energy $E\un{c}$ and a steric factor $P$ recovers the familiar shape,
  [](#eq-collision-theory-k). The theory works for radical reactions but can be orders of magnitude
  off when orientation matters, because it ignores the intermolecular potential entirely.
- An $N$-atom molecule has $3N$ degrees of freedom; subtracting translation and rotation leaves
  $3N-6$ internal coordinates, or $3N-5$ if linear, and the potential energy surface has one
  dimension per internal coordinate.
- Minima, where all second derivatives are positive, [](#eq-minimum-conditions), are stable
  structures; first-order saddle points, with exactly one negative second derivative,
  [](#eq-saddle-conditions), are transition states, and the direction of negative curvature is the
  reaction coordinate. Well depth maps to enthalpy and well width to entropy.
- Reactions are rare events: molecular vibrations occur at $10^{12}$–$10^{14}\ \mathrm{s^{-1}}$, so
  a reaction on a one-second timescale takes many trillions of vibrations per reactive event. We
  observe them in bulk only because a mole contains $6.022\times10^{23}$ molecules.
- Transition state theory assumes a quasi-equilibrium between reactants and the transition state and
  no recrossing of the dividing surface. Identifying $k^\ddagger$ with the frequency $\nu^\ddagger$
  along the reaction coordinate makes that frequency cancel, giving the Eyring equation
  $k = (k\un{B}T/h)(K^\ddagger)'(RT/p^0)$, [](#eq-eyring).
- $\left(K^\ddagger\right)'$ is evaluated from molecular partition functions, factorized as
  $q = q\un{trans}q\un{rot}q\un{vib}q\un{el}$, together with the 0 K energy difference
  $\Delta E\un{0}^\ddagger$. The transition-state vibrational partition function excludes the
  reaction coordinate, leaving $3N-7$ modes, or $3N-6$ if linear.
- In thermodynamic form,
  $k\un{TST} = (k\un{B}T/h)\exp(\Delta\un{act}S^\ddagger/R)\exp(-\Delta\un{act}H^\ddagger/RT)$,
  [](#eq-tst-thermo), which identifies $E\un{a} \approx \Delta\un{act}H^\ddagger$,
  [](#eq-tst-arrhenius-ea), and
  $A \approx (k\un{B}T/h)\exp(\Delta\un{act}S^\ddagger/R)$, [](#eq-tst-arrhenius-a). This gives the
  empirical Arrhenius parameters a microscopic meaning.
