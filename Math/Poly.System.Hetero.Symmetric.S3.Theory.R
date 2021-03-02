
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Theory
###
### draft v.0.1a


### Hetero-Symmetric S3 System

# - let P(x, y, z) = polynomial;
### System:
# P(x, y, z) = 0
# P(y, z, x) = 0
# P(z, x, y) = 0

### Trivial solution: x = y = z;
### Non-Trivial system:
# - where at least 2 of (x, y, z) are different;

### Theorem
# The non-trivial system is always decomposable into a simpler system: P[3] o P[simpler];

# - let (x0, y0, z0) be a root tuple;
# - then (y0, z0, x0) & (z0, x0, y0) are also roots;
# Therefore:
# - for each pair of roots, there are other 2 pairs,
#   which have the same value for: S = x0 + y0 + z0.
# - the order of the polynomial in S = Order(non-trivial system) / 3;
