########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log of Log
##
## draft v.0.1a


####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### y = x * log(x + k*log(x))

# Check:
k = sqrt(2);
x = sqrt(3);
params = list(x=x, k=k);
e = expression(x * log(x + k*log(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - y - x * (x+k) / (x + k*log(x)) # = 0

# D2 =>
x*d2y + (x+k)^2 / (x + k*log(x))^2 +
	- (2*x+k) / (x + k*log(x)) # = 0

### ODE:
x^3*(x+k) * d2y + (x+k)*(x*dy - y)^2 - x*(2*x+k) * (x*dy - y) # = 0

