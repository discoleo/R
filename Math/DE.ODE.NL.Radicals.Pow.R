#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Radicals: Power
##
## draft v.0.1b

### Log to Power:
# y = P(x)^n + P(x)^(j*n)


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### ODE: Order 1

### y = (x^2 + k)^n + (x^2 + k)^(2*n)

# Check:
k = 1/sqrt(5);
n = 1/sqrt(7);
x = sqrt(3); params = list(x=x, k=k);
e = expression((x^2 + k)^n + (x^2 + k)^(2*n))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2 + k)*dy - 2*n*x*((x^2 + k)^n + 2*(x^2 + k)^(2*n)) # = 0

# System:
R = (x^2 + k)^n; # =>
R^2 + R - y # = 0
2*n*x*(2*R^2 + R) - (x^2 + k)*dy # = 0
# =>
2*n*x*R + (x^2 + k)*dy - 4*n*x*y  # = 0
# =>
4*n^2*x^2*R^2 + 4*n^2*x^2*R - 4*n^2*x^2*y # = 0
((x^2 + k)*dy - 4*n*x*y)^2 +
	- 2*n*x*((x^2 + k)*dy - 4*n*x*y) - 4*n^2*x^2*y # = 0

### ODE:
(x^2+k)^2 * dy^2 - 8*n*x*(x^2+k) * y*dy + 16*n^2*x^2 * y^2 +
	- 2*n*x*((x^2+k)*dy - 2*n*x*y) # = 0


#########################
#########################

### ODE: Order 2

### y = (x^2 + k)^n + (x^2 + k)^(2*n) + (x^2 + k)^(3*n)
# Note: Order 1 will be presented separately;

# Check:
k = 1/sqrt(5);
n = 1/sqrt(7);
x = sqrt(3); params = list(x=x, k=k);
e = expression((x^2 + k)^n + (x^2 + k)^(2*n) + (x^2 + k)^(3*n))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2 + k)*dy - 2*n*x*((x^2 + k)^n + 2*(x^2 + k)^(2*n) + 3*(x^2 + k)^(3*n)) # = 0

# D2 =>
(x^2 + k)*d2y + 2*x*dy +
	- 4*n^2*x^2*((x^2 + k)^n + 4*(x^2 + k)^(2*n) + 9*(x^2 + k)^(3*n)) / (x^2 + k) +
	- 2*n*((x^2 + k)^n + 2*(x^2 + k)^(2*n) + 3*(x^2 + k)^(3*n)) # = 0

### System:
# - there are distinct options regarding dependent radicals;
R = (x^2 + k)^n; R3 = (x^2 + k)^(3*n); # =>
R^2 + R - y + R3 # = 0
2*n*x*(2*R^2 + R + 3*R3) - (x^2+k)*dy # = 0
(x^2+k)^2*d2y + 2*x*(x^2+k)*dy +
	- 4*n^2*x^2*(4*R^2 + R + 9*R3) +
	- 2*n*(x^2+k)*(2*R^2 + R + 3*R3) # = 0

# TODO: solve system;

