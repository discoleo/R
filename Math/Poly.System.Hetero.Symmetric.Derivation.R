########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S2
### Heterogenous Symmetric
###
### Derivation of Formulas
###
### draft v.0.3a


### Systems:
# - are described in file:
#   Poly.System.Hetero.Symmetric.R;

# this file contains the detailed derivations;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R;
# e.g. round0(), round0.p;


##########################

##########################
### Polynomial Systems ###
##########################

###############
### Order 3 ###
###############

### x^3 + b*y

# x^3 + b1*y = R
# y^3 + b1*x = R

### Extensions:
# E1: x^3 + b2*x*y + b1*y = R; [P3 => P6]
# E2: x^3 + b3*(x*y)^2 + b2*x*y + b1*y = R; [P4 => P8]
# - are discussed in a separate section;

### Solution:

### Diff =>
# x^3 - y^3 - b1*(x-y) = 0
# (x - y)*(x^2 + y^2 + x*y - b1) = 0
# => x = y *OR* x^2 + y^2 + x*y - b1 = 0;
# =>
# (x+y)^2 - x*y - b1 = 0
# x*y = (x+y)^2 - b1;
# x*y = S^2 - b1;

### Sum =>
# (x+y)^3 - 3*x*y*(x+y) + b1*(x+y) = 2*R
# S^3 - 3*(S^2 - b1)*S + b1*S - 2*R = 0
# S^3 - 2*b1*S + R = 0

### Classic Polynomial:
### Derivation:
# b1*y = R - x^3
# =>
# (R - x^3)^3 / b1^3 + b1*x - R = 0
# (R - x^3)^3 + b1^4*x - R*b1^3
# (x^3 - R)^3 - b1^4*x + R*b1^3
# (x^3 + b1*x - R - b1*x)^3 - b1^4*x + R*b1^3
# (x^3 + b1*x - R)^3 - 3*b1*x*(x^3 + b1*x - R)^2 + 3*b1^2*x^2*(x^3 + b1*x - R) - b1^3*x^3 - b1^4*x + R*b1^3
# let: p = (x^3 + b1*x - R)
# p^3 - 3*b1*x*p^2 + 3*b1^2*x^2*p - b1^3*p
# p*(p^2 - 3*b1*x*p + 3*b1^2*x^2 - b1^3)
# (x^3 + b1*x - R)*(x^6 - b1*x^4 - 2*R*x^3 + b1^2*x^2 + b1*R*x + R^2 - b1^3)
# =>
(x^3 - b[1]/2 * x - R)^2 + 3/4 * b[1]^2*x^2 - b[1]^3

err = x^6 - b[1]*x^4 - 2*R*x^3 + b[1]^2*x^2 + b[1]*R*x + R^2 - b[1]^3
round0(err)


