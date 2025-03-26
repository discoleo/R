########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log w. Radicals
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

### y = log(sqrt(x + b0) + d)

# Check:
x = sqrt(3); b0 = 2/5; d = 1/3;
params = list(x=x, b0=b0, d=d);
e = expression(log(sqrt(x + b0) + d))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);
d3y = eval(D(D(D(e, "x"), "x"), "x"), params);


# D =>
2*(sqrt(x + b0) + d)*dy - 1/sqrt(x + b0) # = 0
#
2*(x + b0 + d*sqrt(x + b0))*dy - 1 # = 0

# D2 =>
4*(x + b0 + d*sqrt(x + b0))*d2y +
	+ 2 * dy + 1/sqrt(x + b0)^2 # = 0

### ODE:
2*(x + b0)*d2y + 2*(x + b0)*dy^2 + dy # = 0


### Variant:

# D3 =>
2*(x + b0)*d3y + 4*(x+b0)*dy*d2y + 3*d2y + 2*dy^2 # = 0
2*(x + b0)^2*d3y + 4*(x+b0)^2*dy*d2y + (x + b0)*d2y - dy # = 0

### ODE:
z = dy; dz = d2y; d2z = d3y;
2*(x + b0)^2*d2z + 4*(x+b0)^2*z*dz + (x + b0)*dz - z # = 0


#####################

### y = x * log(sqrt(x + b0) + d)

# Check:
x = sqrt(3); b0 = 2/5; d = 1/3;
params = list(x=x, b0=b0, d=d);
e = expression(x * log(sqrt(x + b0) + d))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*x*(sqrt(x + b0) + d)*dy - 2*(sqrt(x + b0) + d)*y - x^2/sqrt(x + b0) # = 0
#
2*x*(x + b0 + d*sqrt(x + b0))*dy - 2*(x + b0 + d*sqrt(x + b0))*y - x^2 # = 0
d*sqrt(x + b0)*(2*x*dy - 2*y) # =
- 2*(x + b0)*(x*dy - y) + x^2;

# D2 =>
2*x*(x + b0 + d*sqrt(x + b0))*d2y +
	+ x*(2 + d / sqrt(x + b0))*dy +
	- (2 + d / sqrt(x + b0))*y - 2*x # = 0
2*x*(x + b0)*(sqrt(x + b0) + d)*d2y +
	+ x*(2*sqrt(x + b0) + d)*dy +
	- (2*sqrt(x + b0) + d)*y - 2*sqrt(x + b0)*x # = 0
2*x*(x + b0)*(- 2*(x + b0)*(x*dy - y) + x^2 + d^2*(2*x*dy - 2*y))*d2y +
	- x*((4*x + 4*b0 - 2*d^2)*(x*dy - y) - 2*x^2)*dy +
	+ ((4*x + 4*b0 - 2*d^2)*(x*dy - y) - 2*x^2)*y +
	+ 4*x*(x + b0)*(x*dy - y) - 2*x^3 # = 0

### ODE:
4*x^2*(x + b0)*(x + b0 - d^2)*dy*d2y +
	- 4*x*(x + b0)*(x + b0 - d^2)*y*d2y - 2*x^3*(x + b0)*d2y +
	+ 2*x^2*(2*x + 2*b0 - d^2)*dy^2 +
	- 4*x*(2*x + 2*b0 - d^2)*y*dy - x^2*(6*x + 4*b0)*dy +
	+ 2*(2*x + 2*b0 - d^2)*y^2 + x*(6*x + 4*b0)*y + 2*x^3 # = 0

