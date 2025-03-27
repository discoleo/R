########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Linear ODEs - Log w. Radicals
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

x = sqrt(3); b0 = -sqrt(2); d = 1/3; params = list(x=x, b0=b0, d=d);
e = expression(log(sqrt(x + b0) + d))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
2*dy - 1/(sqrt(x + b0) + d) /sqrt(x + b0) # = 0
2*dy - (x + b0 - d*sqrt(x + b0))/(x + b0 - d^2) / (x + b0) # = 0
2*(x + b0)*(x + b0 - d^2)*dy - (x + b0 - d*sqrt(x + b0)) # = 0

# D2 =>
4*(x + b0)*(x + b0 - d^2)*d2y + 4*(2*x + 2*b0 - d^2)*dy +
	- (2 - d/sqrt(x + b0)) # = 0
4*(x + b0)^2*(x + b0 - d^2)*d2y + 4*(x + b0)*(2*x + 2*b0 - d^2)*dy +
	- (2*(x + b0) - d*sqrt(x + b0)) # = 0

### ODE:
4*(x + b0)*(x + b0 - d^2)*d2y + 2*(3*x + 3*b0 - d^2)*dy - 1 # = 0


