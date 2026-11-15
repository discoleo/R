#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Fraction of Log * Exp
##
## v.0.1a


### Examples:

# x*y*d2y - 2*x*dy^2 - (2*k*x-1)*y*dy + k*(k*x-1)*y^3 - k*(k*x - 1)*y^2 = 0;


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#####################

### y = 1 / (Exp(k*x) * Log(x) + 1)

# Check:
k = 1/5;
x = sqrt(3); params = list(x=x, k=k);
e = expression(1 / (exp(k*x) * log(x) + 1))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
dy + (k/y - k + exp(k*x)/x) * y^2 # = 0
x*dy + (exp(k*x)*y - k*x*y + k*x) * y # = 0

# D2 =>
x*d2y + dy + (exp(k*x)*y - k*x*y + k*x) * dy +
	+ (exp(k*x)*dy + k*exp(k*x)*y - k*x*dy - k*y + k) * y # = 0
x*y*d2y - x * dy^2 - k*x * y^2*dy - (k*x-1)*y*dy +
	+ k*((k*x-1)*y - k*x + 1) * y^2 + exp(k*x) * y^2 * dy # = 0

### ODE:
x * y*d2y - 2*x * dy^2 - (2*k*x-1) * y*dy +
	+ k*(k*x-1) * y^3 - k*(k*x - 1) * y^2 # = 0

