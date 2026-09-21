########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs - Coupled SQRT
##
## draft v.0.1a



####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### Type: 1 Component



#########################

### Type: 2 Components

### y = B1(x) * SQRT(P1(x)) * ATAN(SQRT(P1(x))) + B2(x) * SQRT(P2(x)) * ATAN(SQRT(P2(x)))

# Check:
c1 = sqrt(2);
x = 3^(3/5);
params = list(x=x, c1=c1);
e = expression(sqrt(x^2-1) * atan(sqrt(x^2-1)) + c1*sqrt(x^4-1) * atan(sqrt(x^4-1)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x^2-1)*(x^4-1)*dy +
	- x^2*(x^4-1) * sqrt(x^2-1) * atan(sqrt(x^2-1)) +
	- 2*c1*x^4*(x^2-1) * sqrt(x^4-1) * atan(sqrt(x^4-1)) +
	- (2*c1+1)*(x^2-1)*(x^4-1) # = 0

# System:
div = 2*x^4*(x^2-1) - x^2*(x^4-1);
sqrt(x^2-1) * atan(sqrt(x^2-1)) # ==
- (x*(x^2-1)*(x^4-1)*dy - 2*x^4*(x^2-1)*y - (2*c1+1)*(x^2-1)*(x^4-1)) / div;
#
c1*sqrt(x^4-1) * atan(sqrt(x^4-1)) # ==
(x*(x^2-1)*(x^4-1)*dy - x^2*(x^4-1)*y - (2*c1+1)*(x^2-1)*(x^4-1)) / div;

# D2 =>
x*(x^2-1)*(x^4-1) * d2y + (7*x^6 - 5*x^4 - 3*x^2 + 1) * dy +
	- x^3*(x^4-1) * sqrt(x^2-1) * atan(sqrt(x^2-1)) / (x^2-1) +
	- (6*x^5-2*x) * sqrt(x^2-1) * atan(sqrt(x^2-1)) +
	- 4*c1*x^7*(x^2-1) * sqrt(x^4-1) * atan(sqrt(x^4-1)) / (x^4-1) +
	- 2*c1*(6*x^5-4*x^3) * sqrt(x^4-1) * atan(sqrt(x^4-1)) +
	- x*(x^4-1) - 4*c1*x^3*(x^2-1) - (2*c1+1)*(6*x^5 - 4*x^3 - 2*x) # = 0

# TODO


####################

### Mixed Components

