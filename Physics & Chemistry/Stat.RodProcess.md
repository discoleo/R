

# The Rod Process #

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


### Outcomes ###

1. Distribution of Escape velocities: e.g. the rods can only escape (when there is NO damping);

2. Proportion of rods with arbitrarily small velocity (in the variants with damping);

3. Phase Transition:
- a phase transition is likely in the variants with damping;
- the number of rods with velocity close to 0 may increase abruptly at a certain critical "density" (number) of rods;
- number N of rods at which a phase transition is expected, dependent on the value of the damping factor and the parameters of the initial distributions;

4. Time Measures:
- Median Escape Time: at which 50% of rods have escaped;
- Median Slowdown Time: at which the average velocity of the rods has decreased below some critical value (in the model with damping), e.g. below 50% of the average initial velocity;
