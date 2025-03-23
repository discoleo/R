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
	
