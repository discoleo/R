
# Prediction of Protein Structure

## Energy Minimization

Can be done in 3 steps:
- Step 1: Backbone
  - rotations around the backbone bonds subject to some constraints;
- Step 2: minimize surface exposure of non-polar AA;
- Step 3: maximize H-bonds, van der Waals bonds and other ionic interactions between AAs;

**Notes:**
- most proteins can be viewed as linear chains;
- disulfide bonds may pose some challenges, but the chains are still mostly linear;


### Backbone

Backbone: ... - N(H) - C(HR) - C(=O) - ...;

- optimize rotation around each bond;
- intercalated bonds: "V__|", with slight perturbations around the minimum energy configuration;
- reduce surface exposure of non-polar AA;
- maximize H-bonds;

### Constraints:
- non-overlapping atoms, considering the Lennard-Jones potential;
  - can be done initially with a discrete set of precomputed distances for each atom-pair;
- non-polar AA: forbidden on the surface;
  - e.g. Ala, Val, Leu, Ile, Phe, Tyr (with the exception of the - OH), Trp (parts of it);
- polar AA: maintain on the surface (including around pockets);
  - Asp, Glu, Arg, Lys, His;
  - OH: Ser, Tre, -OH of Tyr;

### Precomputed Structures:
- His, Trp, Phy, Tyr;
- Pro: partly; can pre-compute a number of variants;

### Challenges:
- challenging optimization problem: but not due to any quantum-mechanics!
- Cys, Met: S is a larger atom and more variability possible;
- disulfide bonds: additional constraints;
- Pro: can use a number of precomputed structural variants;

