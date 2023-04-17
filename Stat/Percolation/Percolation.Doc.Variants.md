
# Percolation Variants

Intended as subsequent projects:
- not covered by the current project;


## A. Radial Channels

The channels run concentrically around each other (similar to Dante's Inferno):
- Center Channel: source of the liquid;
- each channel has a number of blocks: either various homogeneous distributions, or variable distributions (e.g. as functions of the radial distance, or of sqrt(radial distance));
- each wall has a number of pores;


### Measures

For many of the measured parameters it is meaningful to compute:
- Specialists in Combinatorics: the median and various percentiles (e.g. 0.75 and 0.90);
- Frequentists: expectation;

1. The number of levels the system percolates;
2. Proportion of filled surface (radial length) for each level;
3. Proportion of Back-fill: sub-channels filled from a level further away;


## B. Linear Channels around a Cylinder

Similar to **A** and to the **Priginal**, but encircling a cylinder:
- Inflow: can be the whole Top channel or a specified column;

Variant: Continuous Case
- let the number of Blocks per channel = pb * n and number Pores = pp * n, with pp < pb;
- Limit Case as n tends to Infinity;
- Other examples: Poisson process with lambda = px * n;
