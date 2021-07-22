########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### DE Systems: Polynomial
###
### draft v.0.1a

#############
### Types ###
#############

### Simple:
# Level n: TODO;
### Hetero-Symmetric:
# Simple Order n: TODO;
### Others:
# TODO


###############
### History ###
###############


### draft v.0.1a:
# - system:
#   n*(R - c0*g)*df + c0*f*dg = f*dR;


########################
########################

########################
### Hetero-Symmetric ###
########################

### Derived from:
# f^n + c0*g = R
# g^n + c0*f = R
# where (f, g) = functions in x;
# R = R(x) is a given parameter function;
# c0 = numeric constant;

### D & Mult(f*...) =>
n*(R - c0*g)*df + c0*f*dg - f*dR # = 0
n*(R - c0*f)*dg + c0*g*df - g*dR # = 0

### Examples:
# R = x + b0;
n*(x + b0 - c0*g)*df + c0*f*dg - f # = 0
n*(x + b0 - c0*f)*dg + c0*g*df - g # = 0
# =>
c0*f*dg - c0*g*df + n*(x + b0)*df - f # = 0
c0*g*df - c0*f*dg + n*(x + b0)*dg - g # = 0

### TODO:
# - check;

