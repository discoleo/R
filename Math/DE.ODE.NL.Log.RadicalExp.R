########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log Compositions
## Log( Radical(Exp) )
##
## draft v.0.1a


### Examples:

# x^5*d2y + (x*dy - y)^3 - (x + 2)*x^3 * (x*dy - y) = 0;


####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### y = x * log(sqrt(exp(2*x) + b0) - exp(x))
# Note: Very simple variant;

# Check:
k = -1/5; b0 = 2/3;
# Note: b0 has no impact on ODE;
x = sqrt(3); params = list(x=x, b0=b0);
e = expression(x * log(sqrt(exp(2*x) + b0) - exp(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*dy - y + x^2 * exp(x) / sqrt(exp(2*x) + b0) # = 0

# D2 =>
x*d2y + (x^2 + 2*x) * exp(x) / sqrt(exp(2*x) + b0) +
	- x^2 * exp(3*x) / sqrt(exp(2*x) + b0)^3 # = 0

### ODE:
x^5*d2y + (x*dy - y)^3 - (x + 2)*x^3 * (x*dy - y) # = 0

