########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log of Log
##
## draft v.0.1a


####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### y = x * log(x + k*log(x))

# Check:
k = sqrt(2);
x = sqrt(3);
params = list(x=x, k=k);
e = expression(x * log(x + k*log(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - y - x * (x+k) / (x + k*log(x)) # = 0

# D2 =>
x*d2y + (x+k)^2 / (x + k*log(x))^2 +
	- (2*x+k) / (x + k*log(x)) # = 0

### ODE:
x^3*(x+k) * d2y + (x+k)*(x*dy - y)^2 - x*(2*x+k) * (x*dy - y) # = 0

#######################

# Simple Variant:
# y = log(x + k*log(x))

# Check:
k = sqrt(2);
x = sqrt(3);
params = list(x=x, k=k);
e = expression(log(x + k*log(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - (x+k) / (x + k*log(x)) # = 0

# D2 =>
x^2*d2y + x*dy - x / (x + k*log(x)) +
	+ (x+k)^2 / (x + k*log(x))^2 # = 0

### ODE:
x^2*(x+k) * d2y + x^2*(x+k) * dy^2 + k*x * dy # = 0


##########
# Variant:
# y = log(x + k*log(x)) + c1*log(x)

# Check:
k = sqrt(2);
c1 = - sqrt(5); # c1 = -1; # c1 = 1/2;
x = sqrt(3);
params = list(x=x, k=k, c1=c1);
e = expression(log(x + k*log(x)) + c1*log(x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - (x+k) / (x + k*log(x)) - c1 # = 0

# D2 =>
x^2*d2y + x*dy - x * (x*dy-c1) / (x+k) +
	+ (x*dy - c1)^2 # = 0
x^2*(x+k) * d2y + (x+k)*(x*dy - c1)^2 +
	+ k*x * dy + c1*x # = 0

### ODE:
x^2*(x+k) * d2y + x^2*(x+k) * dy^2 +
	- x*(2*c1*x + 2*k*c1 - k) * dy + c1*(c1+1)*x + k*c1^2 # = 0

### Special Cases:

### c1 = -1;
x^2*(x+k) * d2y +
	+ x^2*(x+k) * dy^2 + x*(2*x + 3*k) * dy + k # = 0

### c1 = 1/2;
4*x^2*(x+k) * d2y +
	+ 4*x^2*(x+k) * dy^2 - 4*x^2 * dy + 3*x + k # = 0

