#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Exp of Radicals
##
## draft v.0.1b

### Theory

# NL ODE derived from:
# y = C1(x) * Exp(B1(x) * P(x)^p + B0(x));
# y = C1(x) * P(x)^p2 * Exp(B1(x) * P(x)^p1 + B0(x));


### Examples

# 3*(x+1) * y*d2y - 3*(x+1) * dy^2 + 2*y*dy - 2*k*y^2 = 0;
# n*(x+1) * y*d2y - n*(x+1) * dy^2 + (n-1)*y*dy - (n-1)*k*y^2 = 0;


####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

######################
###  EXP(Radical)  ###

### y = x^p * exp((x+1)^(1/3) + k*x)

# Check:
p = 1/3; # p = 0;
k = sqrt(2);
x = sqrt(3);
params = list(x=x, k=k, p=p);
e = expression(x^p * exp((x+1)^(1/3) + k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
3*x*(x+1)*dy - 3*p*(x+1)*y - x*((x+1)^(1/3) + 3*k*(x+1)) * y # = 0
3*x*(x+1)*dy - 3*(x+1)*(k*x + p)*y - x*(x+1)^(1/3) * y # = 0

# D2 =>
9*x*(x+1)^2 * d2y +
	- 9*(x+1)*(k*x*(x+1) + (p-2)*x + p-1) * dy +
	- 9*(x+1)*(2*k*x + p+k) * y +
	- 3*x*(x+1) * (x+1)^(1/3) * dy +
	- (4*x+3) * (x+1)^(1/3) * y # = 0

### ODE:
3*x^2*(x+1) * y*d2y - 3*x^2*(x+1) * dy^2 + 2*x^2 * y*dy +
	- (2*k*x^2 - p*x - 3*p) * y^2 # = 0


### Special Cases:

### p = 0;
3*(x+1) * y*d2y - 3*(x+1) * dy^2 + 2*y*dy - 2*k*y^2 # = 0


### Variant:
### y = x^p * exp((x+1)^(1/4) + k*x)

# Check:
p = 1/3; # p = 0;
k = sqrt(2); n = 5; # n = -1/5;
x = sqrt(3);
params = list(x=x, k=k, p=p);
e = expression(x^p * exp((x+1)^(1/n) + k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
n*x^2*(x+1) * y*d2y - n*x^2*(x+1) * dy^2 + (n-1)*x^2 * y*dy +
	- ((n-1)*k*x^2 - p*x - n*p) * y^2 # = 0


### Special Cases:

### p = 0;
n*(x+1) * y*d2y - n*(x+1) * dy^2 + (n-1)*y*dy - (n-1)*k*y^2 # = 0

