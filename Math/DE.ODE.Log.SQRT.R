#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODE: Coupled SQRT * Log(SQRT)
##
## draft v.0.1a


### Coupled SQRT * Log(SQRT)
# y = SQRT(P(x)^2 + b) * Log(SQRT(P(x)^2 + b) - P(x))


### Examples:

# (x^2-1)^2 * d2y - (x^2-1)*(2*x^2 + x - 2) * dy + (x^4 + x^3 - x^2 - x + 2) * y = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### Mixed w. Exp

### y = sqrt(x^2 + b) * log(sqrt(x^2 + b) - x) * exp(k*x)

# Check:
k = 1/sqrt(5); b = exp(1/3); # b = -1 + 0i; k = 1;
x = sqrt(3); params = list(x=x, k=k, b=b);
e = expression(sqrt(x^2 + b) * log(sqrt(x^2 + b) - x) * exp(k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2 + b) * dy - k*(x^2 + b) * y - x*y + (x^2 + b)*exp(k*x) # = 0

# D2 =>
(x^2 + b) * d2y - (k*x^2 - x + b*k) * dy - (2*k*x + 1) * y +
	+ k*(x^2 + b)*exp(k*x) + 2*x*exp(k*x) # = 0
(x^2 + b) * d2y - (2*k*x^2 - x + 2*b*k) * dy +
	+ (k^2*x^2 - k*x + b*k^2 - 1) * y + 2*x*exp(k*x) # = 0

### ODE:
(x^2+b)^2 * d2y - (x^2+b)*(2*k*x^2 + x + 2*b*k) * dy +
	+ ((x^2+b)*(k^2*x^2 + k*x + b*k^2 + 1) - 2*b) * y # = 0

### Special Cases:

### Case: b = -1; k = 1;
(x^2-1)^2 * d2y - (x^2-1)*(2*x^2 + x - 2) * dy +
	+ (x^4 + x^3 - x^2 - x + 2) * y # = 0

