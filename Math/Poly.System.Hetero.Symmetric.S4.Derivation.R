########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
###  == Derivation ==
###
### draft v.0.1b


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R
# - e.g. round0(), round0.p(),
#   solve.EnAll(), solveEn();


########################

###############
### Order 2 ###
###############

### V3:
### x1^2 + b*x2*x3*x4 = R
x1^2 + b*x2*x3*x4 # - R
x2^2 + b*x1*x3*x4 # - R
x3^2 + b*x1*x2*x4 # - R
x4^2 + b*x1*x2*x3 # - R

### Solution:

### Diff Eq[i] - Eq[i+1] =>
(x1 - x2)*(x1 + x2 - b*x3*x4) # = 0
(x2 - x3)*(x2 + x3 - b*x1*x4) # = 0
(x3 - x4)*(x3 + x4 - b*x1*x2) # = 0
(x4 - x1)*(x1 + x4 - b*x2*x3) # = 0
### Diff Eq[i] - Eq[i+2], etc. =>
(x1 - x3)*(x1 + x3 - b*x2*x4) # = 0
(x2 - x4)*(x2 + x4 - b*x1*x3) # = 0
# Case: x[i] != x[j]: Sum =>
3*S - b*E2 # = 0

### Sum =>
S^2 - 2*E2 + b*E3 - 4*R # = 0

### Sum(x1*...) =>
(x1^3 + x2^3 + x3^3 + x4^3) + 4*b*E4 - R*S # = 0
S^3 - 3*E2*S + 3*E3 + 4*b*E4 - R*S # = 0

### Diff(x[i]*Eq[i] - x[i+1]*Eq[i+1]) =>
(x1 - x2)*(x1^2 + x2^2 + x1*x2 - R) # = 0
(x2 - x3)*(x2^2 + x3^2 + x2*x3 - R) # = 0
# ...

### Sum(Diff 1)^2 vs Sum(Diff 2):
x1^2 + x2^2 + 2*x1*x2 - b^2*(x3*x4)^2 # = 0
3*S^2 - 6*E2 + 2*E2 - b^2*(E2^2 - 2*S*E3 + 2*E4) # = 0
### Sum(Diff 2):
3*S^2 - 5*E2 - 6*R # = 0
# =>
E2 - b^2*(E2^2 - 2*S*E3 + 2*E4) + 6*R = 0
# TODO:
# - verify!
# - analyze existence of solutions;
# - possibility to simplify solution using this Diff;


### =>
#   b*E2 = 3*S
# - b^2*E3 = b*S^2 - 6*S - 4*b*R
# - 4*b^3*E4 = b^2*S^3 - 12*b*S^2 - b^2*R*S + 18*S + 12*b*R

### Sum(x1^2*...) =>
(x1^4 + x2^4 + x3^4 + x4^4) + b*E4*S - R*(S^2 - 2*E2) # = 0
S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 - 4*E4 + b*E4*S - R*(S^2 - 2*E2) # = 0
S^4 - R*S^2 + 2*E2^2 - 4*E2*S^2 + 2*R*E2 + 4*E3*S + b*E4*S - 4*E4 # = 0
#
3*b^3*S^4 - 48*b^2*S^3 - 3*b^3*R*S^2 + 102*b*S^2 + 72*b^2*R*S + 72*S + 48*b*R # = 0


#########
### [old] [probably cyclic redundancy]

### Sum(x2*x3*x4*...) =>
# - possibly/probably redundant!
E4*S + b*Sum((x2*x3*x4)^2 ) - R*E3 # = 0
E4*S + b*(E3^2 - 2*E4*E2) - R*E3 # = 0
E4*S - 2*b*E2*E4 + b*E3^2 - R*E3 # = 0

### SEq1:
(10*R*S^2 - 14*R*S^3*b + 40*R^2*S*b - S^4 + S^5*b) +
(- 22*R*S*b^2 + 32*R*b + 6*S - 9*S^2*b + 4*S^3*b^2)*E3^1 +
(3*S*b^3 - 14*b^2)*E3^2

### SEq2:
(- 8*R*S^3 + 12*R*S^4*b - 8*R^2*S - 24*R^2*S^2*b - 32*R^3*b + S^5 - S^6*b) +
(8*R - 28*R*S*b + 10*R*S^2*b^2 + 32*R^2*b^2 - 8*S^2 + 10*S^3*b - 3*S^4*b^2)*E3^1 +
(- 10*R*b^3 + 9*S*b^2 - S^2*b^3 - 8*b)*E3^2 +
(b^4)*E3^3


### E3 =>
Subst = 560*R*S^2 - 780*R*S^3*b - 224*R*S^4*b^2 + 213*R*S^5*b^3 - 9*R*S^6*b^4 + 3024*R^2*S*b - 564*R^2*S^2*b^2 - 768*R^2*S^3*b^3 + 24*R^2*S^4*b^4 + 816*R^3*S*b^3 - 16*R^3*S^2*b^4 + 3136*R^3*b^2 - 56*S^4 + 36*S^5*b + 37*S^6*b^2 - 18*S^7*b^3 + S^8*b^4;
Subst = - Subst;

E3Div = - 324*R*S*b^2 - 256*R*S^2*b^3 + 249*R*S^3*b^4 - 5*R*S^4*b^5 + 1008*R*b - 252*R^2*S*b^4 + 4*R^2*S^2*b^5 - 1408*R^2*b^3 + 336*S - 188*S^2*b - 240*S^3*b^2 + 199*S^4*b^3 - 48*S^5*b^4 + S^6*b^5;

### =>
(10*R*S^2 - 14*R*S^3*b + 40*R^2*S*b - S^4 + S^5*b) +
(- 22*R*S*b^2 + 32*R*b + 6*S - 9*S^2*b + 4*S^3*b^2) * (Subst/E3Div) +
(3*S*b^3 - 14*b^2) * (Subst/E3Div)^2

### Eq:
b^5*S^10 +
	+ 2*b^4*S^9 - (2*b^3 + 6*R*b^5)*S^8  + (19*b^2 - 110*R*b^4)*S^7 +
	+ b*(-34 + 16*R*b^2 + 9*R^2*b^4)*S^6 + (-56 - 432*R*b^2 + 523*R^2*b^4)*S^5 +
	+ b*(568*R + 2026*R^2*b^2 - 4*R^3*b^4)*S^4 + (1120*R + 2472*R^2*b^2 - 596*R^3*b^4)*S^3 +
	- (2624*b*R^2 + 6608*b^3*R^3)*S^2 + (-3584*R^2 - 13184*b^2*R^3 + 256*b^4*R^4)*S +
	+ -7168*b*R^3 + 256*b^3*R^4
# (b*S^3 + 4*S^2 - 64*R) * (b*S + 1) * P[6]
### b*S + 1: Solution to "distinct" system (but x3 == x4);
### P[6]: Solution to degenerate System
# P[6] = (b*S^2 - 2*S - 4*b*R) * P[4]
(- 28*R + R^2*b^2) +
(- 24*R*b)*S^1 +
(7 - 2*R*b^2)*S^2 +
(- b)*S^3 +
(b^2)*S^4

##################

### Debug:
R = 2;
b = 3;
#
x1 =  1.6128981492115229;
x2 = -0.5852711593612232;
x3 = x4 = x2;


### Test
x1^2 + b*x2*x3*x4 # - R
x2^2 + b*x1*x3*x4 # - R
x3^2 + b*x1*x2*x4 # - R
x4^2 + b*x1*x2*x3 # - R


#######################
#######################



###############
### Order 3 ###
###############

### V3:
### x1^3 + b*x2*x3*x4 = R
x1^3 + b*x2*x3*x4 # - R
x2^3 + b*x1*x3*x4 # - R
x3^3 + b*x1*x2*x4 # - R
x4^3 + b*x1*x2*x3 # - R

### Solution:

### Diff Eq[i] - Eq[i+1] =>
(x1 - x2)*(x1^2 + x2^2 + x1*x2 - b*x3*x4) # = 0
# ...
# Case: x[i] != x[j]: Sum =>
3*S^2 - 5*E2 - b*E2 # = 0

### Sum =>
S^3 - 3*E2*S + 3*E3 + b*E3 - 4*R # = 0

### Sum(x1*...) =>
(x1^4 + x2^4 + x3^4 + x4^4) + 4*b*E4 - R*S # = 0
S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 - 4*E4 + 4*b*E4 - R*S # = 0

### Diff(x[i]*Eq[i] - x[i+1]*Eq[i+1]) =>
(x1 - x2)*(x1^3 + x2^3 + x1^2*x2 + x1*x2^2 - R) # = 0
(x2 - x3)*(x2^3 + x3^3 + x2^2*x3 + x2*x3^2 - R) # = 0
# ...
# TODO: sum(Diff);

