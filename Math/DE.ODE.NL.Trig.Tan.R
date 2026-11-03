########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Trig: Tan
##
## draft v.0.1a


### NL ODEs: Trigonometric
#   Trivial Tan


### Examples

# x*d2y - 2*K*x^n * y*dy - (n-1)*dy = 0;
# d2y - K*x^(n-2) * (2*K*x^n * y^3 + (n-1)*y^2 + 2*K*x^n * y + n-1) = 0;


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


########################

### Simple Components

### y = tan(k*x^n)

# Check:
e = expression(tan(k*x^n));
k = sqrt(5); n = 1/sqrt(2);
x = sqrt(3);
params = list(x=x, k=k);
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
dy - k*n*x^(n-1) * (1 + y^2) # = 0

# D2 =>
d2y - k*n*(n-1)*x^(n-2) * (1 + y^2) +
	- 2*k*n*x^(n-1) * y*dy # = 0

### Variant 1:
d2y - k*n*(n-1)*x^(n-2) * (1 + y^2) +
	- 2*k^2*n^2*x^(2*n-2) * y*(1 + y^2) # = 0

### ODE:
d2y - k*n*x^(n-2) * (2*k*n*x^n * y^3 + (n-1)*y^2 + 2*k*n*x^n * y + n-1) # = 0


### Variant 2:
d2y - k*n*(n-1)*x^(n-2) * (1 + y^2) +
	- 2*k^2*n^2*x^(2*n-2) * y*(1 + y^2) # = 0

### ODE:
x*d2y - 2*k*n*x^n * y*dy - (n-1)*dy # = 0

