########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Linear ODEs - Trig w. Radicals
##
## draft v.0.1a

### Linear ODE: Trig w. Radicals
# y = B1(x) * sin(SQRT(R(x))) + B2(x) * SQRT(R(x)) * cos(R(x));

### Theory
# D(Sin(SQRT(...))) => SQRT * Cos(SQRT(...));
# D(SQRT * Cos(SQRT(...))) => Sin(SQRT(...)) + SQRT * Cos(SQRT(...));


####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################


#####################

### y = x * sin(k*SQRT(x)) + c1 * SQRT(x) * cos(k*SQRT(x));

# Check:
k = sqrt(2); c1 = -1/3;
# k = sqrt(2); c1 = 1/k; # k = 1i * sqrt(3); c1 = 1/k;
x = sqrt(3);
params = list(x=x, k=k, c1=c1);
e = expression(x * sin(k*sqrt(x)) + c1*sqrt(x)*cos(k*sqrt(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*x*dy + (c1*k-2) * x * sin(k*sqrt(x)) +
	- (k*x + c1) * sqrt(x) * cos(k*sqrt(x)) # = 0

# System:
(k*x + c1*(c1*k-1)) * x * sin(k*sqrt(x)) # ==
- 2*c1*x * dy + (k*x + c1) * y;
#
(k*x + c1*(c1*k-1)) * sqrt(x) * cos(k*sqrt(x)) # ==
2*x * dy + (c1*k-2) * y;

# D2 =>
4*x^2*d2y + 4*x*dy +
	+ (k*(k*x + c1) + 2*(c1*k-2)) * x * sin(k*sqrt(x)) +
	+ ((c1*k-5)*k*x - c1) * sqrt(x) * cos(k*sqrt(x)) # = 0
# Substitution =>

### ODE:
4*x^2*(k*x + c1*(c1*k-1)) * d2y +
	- 2*x * (3*k*x + c1*(c1*k-1)) * dy +
	+ (k^3*x^2 + k*(c1^2*k^2 - 3*c1*k + 6)*x + 2*c1*(c1*k - 1)) * y # = 0

### Special Cases:

# Case: c1*k = 1;
4*x^3 * d2y - 6*x^2 * dy + x*(k^2*x + 4) * y # = 0

