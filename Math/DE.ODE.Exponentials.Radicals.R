########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs - Exponentials w. Radicals
##
## draft v.0.1a



####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### Simple: SQRT

### y = exp(x + k*sqrt(x + b0))

x = sqrt(3); k = -1/5; b0 = 2/3; params = list(x=x, k=k, b0=b0);
e = expression(exp(x + k*sqrt(x + b0)))[[1]];
#
y   = exp(x + k*sqrt(x + b0))
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
2*dy - (2 + k/sqrt(x + b0)) * y # = 0

# D2 =>
4*d2y - 2*(2 + k/sqrt(x + b0)) * dy +
	+ k/(x + b0) / sqrt(x + b0) * y # = 0
4*d2y - 8*dy + (4 - k^2/(x+b0)) * y +
	+ k/(x + b0) / sqrt(x + b0) * y # = 0

### ODE:
4*(x + b0)*d2y - (8*x + 8*b0 - 2)*dy +
	+ (4*x - k^2 + 4*b0 - 2) * y # = 0

### Special Cases:

### k = sqrt(4*b0-2)
x = sqrt(3); b0 = 2/3; k = sqrt(4*b0-2);
params = list(x=x, k=k, b0=b0);
# Re-run assignments from above;
4*(x + b0)*d2y - (8*x + 8*b0 - 2)*dy + 4*x*y # = 0


### Variant:
### y = x * exp(x + k*sqrt(x + b0))

x = sqrt(3); k = -1/5; b0 = 2/3; params = list(x=x, k=k, b0=b0);
e = expression(x * exp(x + k*sqrt(x + b0)))[[1]];
#
y   = x * exp(x + k*sqrt(x + b0));
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
2*x*dy - (2 + 2*x + k*x/sqrt(x + b0)) * y # = 0

# D2 =>
4*x*d2y - (2 + k/sqrt(x + b0)) * (2 + 2*x + k*x/sqrt(x + b0)) * y +
	- k*(1 + b0/(x + b0))/sqrt(x + b0) * y - 4*y # = 0
4*x*d2y - k*(4*x + 3 + b0/(x + b0)) / sqrt(x + b0) * y +
	- (k^2 + 8 + 4*x - k^2*b0/(x + b0)) * y # = 0
4*x*(x + b0)*d2y - k*((4*x + 3)*(x + b0) + b0) / sqrt(x + b0) * y +
	- ((k^2 + 8 + 4*x)*(x + b0) - k^2*b0) * y # = 0

### ODE:
4*x^2*(x + b0)*d2y - 2*x*(4*x^2 + (4*b0+3)*x + 4*b0) * dy +
	+ (4*x^3 + (4*b0 + 6 - k^2)*x^2 + 8*b0*x + 6*x + 8*b0) * y # = 0


### Special Cases:

### k = sqrt(4*b0 + 6)
x = sqrt(3); b0 = 2/3; k = sqrt(4*b0 + 6);
params = list(x=x, k=k, b0=b0);
# Re-run assignments from above;
4*x^2*(x + b0)*d2y - 2*x*(4*x^2 + (4*b0+3)*x + 4*b0) * dy +
	+ (4*x^3 + 8*b0*x + 6*x + 8*b0) * y # = 0

### & b0 = -3/4
x = sqrt(5); b0 = -3/4; k = sqrt(4*b0 + 6); # k = sqrt(3);
params = list(x=x, k=k, b0=b0);
# Re-run assignments from above;
x^2*(4*x - 3)*d2y - 2*x*(4*x^2 - 3) * dy + (4*x^3 - 6) * y # = 0


##############################

### y = exp(x + sqrt(x^2 + k))

x = - sqrt(3); k = -1/5; params = list(x=x, k=k);
e = expression(exp(x + sqrt(x^2 + k)))[[1]];
#
y   = exp(x + sqrt(x^2 + k))
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
dy - (1 + x/sqrt(x^2 + k)) * y # = 0

# D2 =>
d2y - (1 + x/sqrt(x^2 + k)) * dy +
	- (1 - x^2/(x^2 + k)) / sqrt(x^2 + k) * y # = 0
d2y - 2*dy + k/(x^2 + k) * y +
	- (1 - x^2/(x^2 + k)) / sqrt(x^2 + k) * y # = 0
d2y - 2*dy + k/(x^2 + k) * y +
	- k/(x^2 + k) / sqrt(x^2 + k) * y # = 0
x*d2y - 2*x*dy + k*x/(x^2 + k) * y +
	- k/(x^2 + k) * (dy - y) # = 0

### ODE:
x*(x^2 + k)*d2y - (2*x*(x^2 + k) + k)*dy + k*(x+1)*y # = 0


###################

### y = exp(x^2 + sqrt(x^2 + k))

x = - sqrt(3); k = -2/5; params = list(x=x, k=k);
e = expression(exp(x^2 + sqrt(x^2 + k)))[[1]];
#
y   = exp(x^2 + sqrt(x^2 + k))
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
dy - x*(2 + 1/sqrt(x^2 + k)) * y # = 0

# D2 =>
d2y - x*(2 + 1/sqrt(x^2 + k)) * dy +
	- 2*y - (1 - x^2/(x^2 + k)) / sqrt(x^2 + k) * y # = 0
d2y - 4*x*dy + (4*x^2 - 1)*y + k/(x^2 + k) * y +
	- 2*y - (1 - x^2/(x^2 + k)) / sqrt(x^2 + k) * y # = 0
x*d2y - 4*x^2*dy - k/(x^2 + k) * dy +
	+ x*(4*x^2 - 3)*y + 3*k*x/(x^2 + k) * y # = 0
x*(x^2 + k)*d2y - (4*x^2*(x^2 + k) + k)*dy +
	+ x*((4*x^2 - 3)*(x^2 + k) + 3*k)*y # = 0

### ODE:
x*(x^2 + k)*d2y - (4*x^2*(x^2 + k) + k)*dy +
	+ x^3*(4*x^2 + 4*k - 3)*y # = 0

### Special Cases:

### k = 3/4
x = - sqrt(3); k = 3/4; params = list(x=x, k=k);
# Re-run assignments above;
x*(4*x^2 + 3)*d2y - (16*x^4 + 12*x^2 + 3)*dy + 16*x^5*y # = 0
	
