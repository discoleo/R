########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Trig w. Radicals
##
## draft v.0.1a



####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


######################
######################

### Example:
# y = sin(k*(x^n + 1)^(1/n))

# Check:
n = 3^(2/3);
k = sqrt(5);
x = 2^(2/3);
params = list(x=x, n=n, k=k);
e = expression(sin(k*(x^n + 1)^(1/n)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^n + 1)*dy - k*x^(n-1) * cos(k*(x^n + 1)^(1/n)) * (x^n + 1)^(1/n) # = 0

# D2 =>
(x^n + 1)^2 * d2y + n*x^(n-1)*(x^n + 1) * dy +
	- k*(n-1)*x^(n-2)*(x^n + 1) * cos(k*(x^n + 1)^(1/n)) * (x^n + 1)^(1/n) +
	+ k^2*x^(2*n-2) * sin(k*(x^n + 1)^(1/n)) * (x^n + 1)^(2/n) +
	- k*x^(2*n-2) * cos(k*(x^n + 1)^(1/n)) * (x^n + 1)^(1/n) # = 0
x*(x^n + 1)^2 * d2y + n*x^n*(x^n + 1) * dy +
	- (x^n + 1) * (n*x^n + n-1) * dy +
	+ k^2*x^(2*n-1) * (x^n + 1)^(2/n) * y # = 0
x*(x^n+1) * d2y - (n-1) * dy +
	+ x*(x^n+1) * y / (1-y^2) * dy^2 # = 0

### ODE:
x*(x^n+1) * (1-y^2)*d2y - (n-1)*(1-y^2)*dy +
	+ x*(x^n+1) * y * dy^2 # = 0

