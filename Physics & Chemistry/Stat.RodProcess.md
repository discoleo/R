

# A. Linear Rod Process #

Inspired by the discrete hard rod process (see Ref. 1). However, that process introduces discontinuities during the swap of rods. Although it is a "discrete" model, I am rather unhappy with such discontinuities.

1. CRM: Chiara Franceschini: Some results for hard rods of inhomogeneous size.
  > https://www.youtube.com/watch?v=wmUlfSVO1Fo

## Description ##

Inhomogeneous hard rods interacting inertially during collisions:

**Setting:**
- 1D tube of infinite length;
- hard rods are placed inside this tube: it is assumed that there is NO friction;

**Rods:**
- N rods of random length: e.g. uniform distribution between (L1, L2);
- the rods are placed randomly at positions between [0, Xf]:
  - the centre of the rod indicates the position of the rod;
  - rods do NOT intersect;
  - TODO: work out overlaps between rods in a way that preserves the properties of the underlying distribution;
  - TODO: clarify usefulness of centre vs start of rod;
- each rod has an initial velocity, randomly selected between [-vm, vm]:
  - rods with negative velocity move in opposite direction;
  - velocities are randomly distributed, e.g. uniformly, or based on a truncated Gaussian distribution;

**Events:**
- the rods can collide: momentum and kinetic energy are preserved during each collision;
- the "mass" is considered proportional to the rod length;
- momentum ~ l * v; energy ~ l * v^2, where l = length of the specific rod;

**Variants:**
1. No damping during each collision;
2. Constant damping during each collision;
3. Damping proportional to abs(diff(v)) during each collision;

For simple damping:\
E{after collision} = (1 - damping) * E{before collision};


## Outcomes ##

1. Distribution of Escape velocities: e.g. the rods can only escape (when there is NO damping);
  - Escape = NO further collisions possible;

2. Proportion of rods with arbitrarily small velocity (in the variants with damping);

3. Phase Transition:
- a phase transition is likely in the variants with damping;
- the number of rods with velocity close to 0 may increase abruptly at a certain critical "density" (number) of rods;
- number N of rods at which a phase transition is expected, dependent on the value of the damping factor and the parameters of the initial distributions;

4. Time Measures:
- Median Escape Time: at which 50% of rods have escaped (or, for the variants with damping: 50% of rods in the escaping compartment);
- Median Slowdown Time: at which the average velocity of the rods has decreased below some critical value (in the model with damping), e.g. below 50% of the average initial velocity;
> Note: the median velocity may behave badly due to the different sizes of the compartments that will arise (Escape compartment vs quasi-Stationary compartment);


====

# B. Circular Rod Process #

Circular rods are placed in a large circular tube of radius R.

For simplicity, the N rods can be placed equidistantly. The remaining features remain the same as in the Linear rod model: rod lengths and initial velocities follow some given distributions.

## Outcomes ##

1. Time Measures:
- Mean/Median Time until the velocities have almost-converged to a uniform velocity: e.g. variance of velocities less than a critical value or, alternatively, until the variance decreased to 50% of the initial variance;
- Median Slowdown Time: for the models with damping (similar to the Linear rod model);

2. Terminal velocity

## Open Questions ##

1. Do the rods group/aggregate into 2 separate compartments?
> The group-velocity of the 2 compartments has to be equal for the 2 separate compartments to be stable.

2. Is it possible to separate the rods into more than 2 large compartments? Which values of the initial parameters would induce such behaviour?

3. Median Time of survival of the compartments


====

# C. Mathematical Alternatives #

The formulas for the energy and for the momentum can be varied, e.g. E = l * h^3, or p = l^2 * v, where l = length of rod. Note: constants have been omitted from any of the formulas.

Although these formulas may not represent any real physical phenomena, the resulting mathematical models may be still interesting.


====

# D. Spin Models #

Circular disks can be used instead of rods. The disks can spin (positive or negative spin). Conservation of spin has to be added to the formulas.

Note: These models are more complicated and will be moved to a separate file (the disks are not rods anyway).
