#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Atan of SQRT(Exp)
##
## draft v.0.1a


### Examples:

# k*x^5 * d2y + (x*dy - y)^3 + k*x^3*(k*x - 2) * (x*dy - y) = 0;


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### y = x * Atan(Sqrt(Exp(k*x) - 1))
# Note: Exp(x/k) works better;

# Check:
k = 1/5;
x = sqrt(3); params = list(x=x, k=k);
e = expression(x * atan(sqrt(exp(2*k*x) - 1)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*dy - y - k*x^2 / sqrt(exp(2*k*x) - 1) # = 0
# =>
k^2*x^4 / (x*dy - y)^2 - exp(2*k*x) + 1 # = 0

# D2 =>
x*d2y - 2*k*x / sqrt(exp(2*k*x) - 1) +
	+ k^2*x^2 * exp(2*k*x) / sqrt(exp(2*k*x) - 1)^3 # = 0
k*x^5*d2y - 2*k*x^4*dy + 2*k*x^3*y + exp(2*k*x) * (x*dy - y)^3 # = 0
k*x^5*d2y - 2*k*x^4*dy + 2*k*x^3*y +
	+ (k^2*x^4 / (x*dy - y)^2 + 1) * (x*dy - y)^3 # = 0

### ODE:
k*x^5 * d2y + (x*dy - y)^3 + k*x^3*(k*x - 2) * (x*dy - y) # = 0

