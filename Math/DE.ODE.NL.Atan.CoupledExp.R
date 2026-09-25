########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs - ATAN w. Coupled Exp
##
## draft v.0.1b


### Examples:

# y*d2y - k*dy - 4*k^2 * y^2 + k^2 = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y = cosh(k*x) * atan(exp(k*x))
# - Simple example;

# Check:
k = sqrt(3); # k = 1;
x = 3^(3/5);
params = list(x=x, k=k);
e = expression(cosh(k*x) * atan(exp(k*x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*dy - 2*k*sinh(k*x)* atan(exp(k*x)) - k # = 0

# Linear System:
# [actually NOT needed]
exp(k*x) * atan(exp(k*x)) # ==
(2*dy - k + 2*k*y) / (2*k);
#
exp(-k*x) * atan(exp(k*x)) # ==
- (2*dy - k - 2*k*y) / (2*k);

# D2 =>
2*d2y - 2*k^2 * y - k^2*sinh(k*x) / cosh(k*x) # = 0
4*d2y - 4*k^2 * y - k*(2*dy - k) / y # = 0

### ODE:
4*y*d2y - 2*k*dy - 4*k^2 * y^2 + k^2 # = 0

# Variant:
kd = k / 2;
y*d2y - kd*dy - 4*kd^2 * y^2 + kd^2 # = 0


#########################

### Extension: y = cosh(k*x^n) * atan(exp(k*x^n))

# Check:
k = sqrt(3); # k = 1;
n = 2/5;
x = 3^(3/5);
params = list(x=x, k=k);
e = expression(cosh(k*x^n) * atan(exp(k*x^n)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*dy - 2*k*n*x^(n-1)*sinh(k*x^n) * atan(exp(k*x^n)) - k*n*x^(n-1) # = 0

# D2 =>
2*d2y - 2*k^2*n^2*x^(2*n-2)*cosh(k*x^n) * atan(exp(k*x^n)) +
	- k^2*n^2*x^(2*n-2)*sinh(k*x^n) / cosh(k*x^n) +
	- 2*k*n*(n-1)*x^(n-2)*sinh(k*x^n) * atan(exp(k*x^n)) +
	- k*n*(n-1)*x^(n-2) # = 0
2*x*d2y - 2*k^2*n^2*x^(2*n-1) * y +
	- k*n*x^n / 2 * 2*k*n*x^(n-1)*sinh(k*x^n) * atan(exp(k*x^n)) / y +
	- (n-1) * (2*dy - k*n*x^(n-1)) +
	- k*n*(n-1)*x^(n-1) # = 0

### ODE:
4*x*y*d2y - 4*(n-1) * y*dy - 2*k*n*x^n * dy +
	- 4*k^2*n^2*x^(2*n-1) * y^2 + k^2*n^2*x^(2*n-1) # = 0

