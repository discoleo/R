########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Linear ODEs - Mixed / Products
##
## draft v.0.1a



####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### y = B(x) * Log(P1(x)) * Log(P2(x))

### Simple: LOG()^2

### y = x * log(x)^2

# Check:
# for Homogenous: p0 = 0;
x = sqrt(3); p0 = -1/2; params = list(x=x, p0=p0);
e = expression(x * log(x)^2 + p0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*d2y - x*dy + y - 2*x - p0 # = 0


# D =>
x*dy - y - 2*x*log(x) + p0 # = 0

# D2 =>
x*d2y - 2*log(x) - 2 # = 0
x^2*d2y - x*dy + y - 2*x - p0 # = 0

