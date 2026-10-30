########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log of Atan
##
## draft v.0.1a


### Examples:



####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### y = x * log(atan(x/k))
# Note: Very simple variant;

# Check:
k = 1/5;
x = sqrt(3); params = list(x=x, k=k);
e = expression(x * log(atan(x/k)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*(x^2+k^2)*dy - (x^2+k^2)*y - k*x^2 / atan(x/k) # = 0

# D2 =>
x*(x^2+k^2)*d2y + 2*x^2*dy - 2*x*y +
	- 2*k*x / atan(x/k) + k^2*x^2 / (x^2+k^2) / atan(x/k)^2 # = 0
x^2*(x^2+k^2)*d2y - 2*k^2*x*dy + 2*k^2*y +
	+ k^2*x^3 / (x^2+k^2) / atan(x/k)^2 # = 0
x^3*(x^2+k^2) * d2y - 2*k^2*x^2 * dy + 2*k^2*x * y +
	+ (x^2+k^2)*(x*dy - y)^2 # = 0

### ODE:
x^3*(x^2+k^2) * d2y + x^2*(x^2+k^2) * dy^2 - 2*x*(x^2+k^2) * y*dy +
	+ (x^2+k^2) * y^2 - 2*k^2*x^2 * dy + 2*k^2*x * y # = 0

