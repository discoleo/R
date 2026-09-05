#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Power of Log
##
## draft v.0.1a



####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

##############
###  SQRT  ###

### y = x^p ^ SQRT(x + k*log(x))

# Check:
p = 1/3; # p = -1; # p = -1/2;
k = sqrt(2);
x = sqrt(3);
params = list(x=x, k=k, p=p);
e = expression(x^p * sqrt(x + k*log(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*x*dy - 2*p*y - x^(2*p) * (x+k) / y # = 0
2*x*y*dy - 2*p*y^2 - x^(2*p) * (x+k) # = 0

# D2 =>
2*x*y*d2y + 2*x*dy^2 - (4*p-2) * y*dy +
	- 2*p*x^(2*p-1) * (x+k) - x^(2*p) # = 0
# Entanglement:
2*x^2*y*d2y + 2*x^2*dy^2 - 2*(2*p-1)*x * y*dy +
	- 2*(2*p+1) * (x*y*dy - p*y^2) + k*x^(2*p) # = 0

### ODE:
x^2 * y*d2y + x^2 * dy^2 - 4*p*x * y*dy +
	+ p*(2*p+1) * y^2 + k/2 * x^(2*p) # = 0

### Special Cases:

### p = -1/2;
x^3 * y*d2y + x^3 * dy^2 + 2*x^2 * y*dy + k/2 # = 0


