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

