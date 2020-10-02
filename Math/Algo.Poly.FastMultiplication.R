
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Multiplication:
### Fast Multiplication Algorythm
###
### draft v.0.1


####################
### Introduction ###

### Terminology:

# - let P(x) be a polynomial of order n;
# - let r[i] be the roots of this polynomial;
# - we do NOT need to compute explicitly these roots;

### Task:

# Compute (P(x))^k, where k >> n;
# - computation should be as efficient as possible;


### Solution:

# P(x) = Prod(x - r[i])
# =>
# (P(x))^k = Prod((x - r[i])^k)

# (x - r[i])^k are binomial expansions;
# - we can compute one binomial expansion and then replicate it n times;
# - compute the parametric product of these n polynomials (n << k);
# - decompose the resulting coefficents into elementary polynomials based on r[i];
# - compute the values using the values of the elementary polynomials;

