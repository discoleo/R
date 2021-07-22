########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### DE Systems: Polynomial
###
### draft v.0.1b

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


### draft v.0.1a - v.0.1b:
# - system:
#   n*(R - c0*g)*df + c0*f*dg = f*dR;
# - system:
#   n*(R1 - c0*g)*(g*df + f*dg) + c0*f*g*dg = f*g*dR1;


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
c0*f*dg - n*c0*g*df + n*(x + b0)*df - f # = 0
c0*g*df - n*c0*f*dg + n*(x + b0)*dg - g # = 0

### TODO:
# - check;


#########################
#########################

### Derived from:
# (f*g)^n + c0*g = R1
# (f*g)^n + c0*f = R2
# where (f, g) = functions in x;
# R1, R2 = is a given parameter function;
# c0 = numeric constant;

### D & Mult(f*g*...) =>
n*(R1 - c0*g)*(g*df + f*dg) + c0*f*g*dg - f*g*dR1 # = 0
n*(R2 - c0*f)*(g*df + f*dg) + c0*f*g*df - f*g*dR2 # = 0

### Examples:
# R1 = R2 = x + b0;
n*(x + b0 - c0*g)*(g*df + f*dg) + c0*f*g*dg - f*g # = 0
n*(x + b0 - c0*f)*(g*df + f*dg) + c0*f*g*df - f*g # = 0
# =>
n*c0*g^2*df + (n-1)*c0*f*g*dg - n*(x + b0)*(g*df + f*dg) + f*g # = 0
n*c0*f^2*dg + (n-1)*c0*f*g*df - n*(x + b0)*(g*df + f*dg) + f*g # = 0

### TODO:
# - check;

