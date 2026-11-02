########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Trig( Exp * SQRT )
##
## draft v.0.1a


### NL ODEs: Trigonometric-Type
#   Subtype: Trig( EXP * SQRT )
# - Inhomogeneous NL ODE;


### Examples:

# 4*(x+b)*(2*k*(x+b) + 1) * (y*d2y - dy^2) +
#  - 2*((2*k*(x+b) + 1)^2 - 2) * (y*dy + k*y^2) + (2*k*(x+b) + 1)^3 = 0;


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


########################

### Inversely-Coupled Exp

### y = sin(exp(k*x) * sqrt(x^n + b)) * exp(-k*x)

# Check:
e = expression(sin(exp(k*x) * sqrt(x^n + b)) * exp(-k*x));
k = sqrt(5);
n = 1/sqrt(2); b = 2/3; # n = 1;
x = sqrt(3);
params = list(x=x, k=k);
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
2*(x^n + b)*dy + 2*k*(x^n + b)*y +
	- (2*k*(x^n + b) + n*x^(n-1)) * sqrt(x^n + b) * cos(exp(k*x) * sqrt(x^n + b)) # = 0

# D2 =>
4*(x^n + b) * d2y + 2*(2*k*x^n + n*x^(n-1) + 2*b*k) * dy + 2*k*n*x^(n-1) * y +
	- 2*(2*k*n*x^(n-1) + n*(n-1)*x^(n-2)) * sqrt(x^n + b) * cos(exp(k*x) * sqrt(x^n + b)) +
	+ (2*k*(x^n + b) + n*x^(n-1))^2 * exp(k*x) * sin(exp(k*x) * sqrt(x^n + b)) # = 0
4*(x^n + b) * y*d2y + 2*(2*k*x^n + n*x^(n-1) + 2*b*k) * y*dy + 2*k*n*x^(n-1) * y^2 +
	- 4*(x^n + b)*(2*k*n*x^(n-1) + n*(n-1)*x^(n-2)) * y*(dy + k*y) / (2*k*(x^n + b) + n*x^(n-1)) +
	+ (2*k*(x^n + b) + n*x^(n-1))^2 * sin(exp(k*x) * sqrt(x^n + b))^2 # = 0
4*(x^n + b) * y*d2y + 2*(2*k*x^n + n*x^(n-1) + 2*b*k) * y*dy + 2*k*n*x^(n-1) * y^2 +
	- 4*(x^n + b)*(2*k*n*x^(n-1) + n*(n-1)*x^(n-2)) * y*(dy + k*y) / (2*k*(x^n + b) + n*x^(n-1)) +
	- 4*(x^n + b)*(dy + k*y)^2 + (2*k*(x^n + b) + n*x^(n-1))^2 # = 0

# TODO: simplify;


### Special Cases:

### Case: n = 1;
4*(x+b)*(2*k*(x+b) + 1) * (y*d2y - dy^2) +
	- 2*((2*k*(x+b) + 1)^2 - 2) * (y*dy + k*y^2) +
	+ (2*k*(x+b) + 1)^3 # = 0

