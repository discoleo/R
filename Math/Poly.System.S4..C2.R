########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S4: C2-Hetero-Symmetric
###
### draft v.0.1a


##############
### Theory ###
##############

### System:
# x1^n1 + x2^n1 + y1^n1 + y2^n1 = R1
# x1^n2*y1^n3 + x2^n2*y2^n3 = R2
# x1*x2*y1 + x1*x2*y2 + x1*y1*y2 + x2*y1*y2 = R3
# x1*x2*y1*y2 = R4

### Symmetries:
# - if (x1, x2, y1, y2) is a solution,
#   so is also the C2 permutation: (x2, x1, y2, y1);

##################
### Basic Type ###
##################

### Symmetric Eq 2
# n2 = n3 = n
# (x1*y1)^n + (x2*y2)^n = R2
# Note:
# - has additional symmetry;
# - if (x1, x2, y1, y2) is a solution,
#   so is also the permutation: (y1, y2, x1, x2) & its C2 permutation;

### Solution:
# (Eq 2) * (x1*y1)^n =>
(x1*y1)^(2*n) - R2*(x1*y1)^n + R2^n # = 0
# - solve for (x1*y1) & (x2*y2);

### Eq 1 & Eq 3:

### Case 1: n1 = 1;
(x1 + y1) + (x2 + y2) - R1 # = 0
(x2*y2)*(x1 + y1) + (x1*y1)*(x2 + y2) - R3 # = 0
# - solve for: (x1 + y1) & (x2 + y2);

### Solve 2x2:
# x1 + y1 = s1;
# x1 * y1 = p1;
# and:
# x2 + y2 = s2;
# x2 * y2 = p2;

# TODO: implement & check;


### Case 2: n1 > 1;
# TODO

