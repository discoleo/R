########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log in Fractions
##
## draft v.0.1a


####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### y = 1 / (x + k*log(x))

# Check:
k = sqrt(2);
x = sqrt(3);
params = list(x=x, k=k);
e = expression(1 / (x + k*log(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy + (x+k) / (x + k*log(x))^2 # = 0
x*dy + (x+k)*y^2 # = 0

# D2 =>
x*d2y + 2*(x+k)*y*dy + y^2 + dy # = 0
# Substitution y^2 =>

### ODE:
x*(x+k) * d2y + 2*(x+k)^2 * y*dy + k*dy # = 0

#######################

### Partial Generalization
# y = x^p / (x + k*log(x))

# Check:
p = 1/3; # p = -1;
k = sqrt(2);
x = sqrt(3);
params = list(x=x, k=k, p=p);
e = expression(x^p / (x + k*log(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x^(p+1) * dy - p*x^p * y + (x+k)*y^2 # = 0

# D2 & Substitution =>

### ODE:
x^(p+1)*(x+k) * d2y + 2*(x+k)^2 * y*dy + k*x^p * dy +
	- (p^2-p)*x^p * y - k*p^2*x^(p-1) * y # = 0

### Special Cases:

# Case: p = -1;
x^2*(x+k) * d2y + 2*x^2(x+k)^2 * y*dy + k*x * dy - (2*x + k) * y # = 0

