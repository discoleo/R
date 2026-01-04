########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Linear ODEs - Mixed / Products
##
## draft v.0.1a


### Examples:

# x^2*d2y - (2*p-1)*x*dy + p^2*y - 2*x^p - p^2*c0 # = 0


####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### Type: LOG()^2
# y = B(x) * Log(P(x))^2 + P(x)

### y = x * log(x)^2

# Check:
# for Homogenous: p0 = 0;
x = sqrt(3); p0 = -1/2; params = list(x=x, p0=p0);
e = expression(x * log(x)^2 + p0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*d2y - x*dy + y - 2*x - p0 # = 0


# D =>
x*dy - y - 2*x*log(x) + p0 # = 0

# D2 =>
x*d2y - 2*log(x) - 2 # = 0
x^2*d2y - x*dy + y - 2*x - p0 # = 0


############################

### y = (x+b0) * log(x+b0)^2
# - just a shift;

# Check:
# for Homogenous: p0 = 0;
x = sqrt(3); b0 = 1; p0 = -1/2; params = list(x=x, b0=b0, p0=p0);
e = expression((x+b0) * log(x+b0)^2 + p0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
(x+b0)^2*d2y - (x+b0)*dy + y - 2*(x+b0) - p0 # = 0


# D =>
dy - log(x+b0)^2 - 2*log(x+b0) # = 0
(x+b0)*dy - y - 2*(x+b0)*log(x+b0) + p0 # = 0

# D2 =>
(x+b0)*d2y - 2*log(x+b0) - 2 # = 0
(x+b0)^2*d2y - ((x+b0)*dy - y + p0) - 2*(x+b0) # = 0

############

### y = x^p * log(x)^2
# - Simple with exponent;

# Check:
# for Homogenous: p0 = 0;
x = sqrt(3); p = 2/3; c0 = -1/2; params = list(x=x, p=p, c0=c0);
e = expression(x^p * log(x)^2 + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*d2y - (2*p-1)*x*dy + p^2*y - 2*x^p - p^2*c0 # = 0


# D =>
x*dy - p*y - 2*x^p*log(x) + p*c0 # = 0

# D2 =>
x*d2y - (p-1)*dy - 2*p*x^(p-1)*log(x) - 2*x^(p-1) # = 0
x^2*d2y - (p-1)*x*dy - p*(x*dy - p*y + p*c0) - 2*x^p # = 0


###############

### y = x * log(x+b0)^2

# Check:
# for Homogenous: p0 = 0;
x = sqrt(3); b0 = 1; p0 = -1/2; params = list(x=x, b0=b0, p0=p0);
e = expression(x * log(x+b0)^2 + p0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*(x+b0)^2*d2y - x*(x+2*b0)*(x+b0)*dy + (x+2*b0)*(x+b0)*y +
	- 2*x^3 - p0*(x+2*b0)*(x+b0) # = 0


# D =>
x*(x+b0)*dy - (x+b0)*y - 2*x^2*log(x+b0) + p0*(x+b0) # = 0

# D2 =>
x*(x+b0)*d2y + x*dy - y - 4*x*log(x+b0) - 2*x^2/(x+b0) + p0 # = 0
x^2*(x+b0)*d2y + x^2*dy - x*y - 2*(x*(x+b0)*dy - (x+b0)*y + p0*(x+b0)) +
	- 2*x^3/(x+b0) + p0*x # = 0
x^2*(x+b0)*d2y - x*(x+2*b0)*dy + (x+2*b0)*y +
	- 2*x^3/(x+b0) - p0*(x+2*b0) # = 0

