########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Linear ODEs - Squares
##
## draft v.0.1j


### Base-Solution:
# Simple: 1 Component
# y = B(x) * Log(P(x))^2 + F0(x);
# y = B(x) * Atan(P(x))^2 + F0(x);
# w. Radicals:
# y = B(x) * Log(sqrt(P(x)^2 + b0) - P(x))^2 + F0(x);
# y = B(x) * Atan(sqrt(P(x)))^2 + F0(x);

# Full: 2 Components
# y = B1(x) * Log(P(x))^2 + B2(x) * Log(P(x)) + F0(x);
# y = B1(x) * Atan(P(x))^2 + B2(x) * Atan(P(x)) + F0(x);


### Examples:

# x^2*d2y - (2*p-1)*x*dy + p^2*y - 2*x^p - p^2*c0 # = 0
# 2*x^2*(x-1)*d2y - x*(x-2)*dy + (x-2)*y - (c0+1)*x + 2*c0 # = 0
# 2*x^2*(x-1)*d2y + x^2*dy - x*(y - c0) - (x + 2) # = 0


####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### Type: LOG()^2

### Components: 1
# y = B(x) * Log(P(x))^2 + F0(x)

### y = x * log(x)^2

# Check:
# for Quasi-Homogenous: c0 = 0;
x = sqrt(3); c0 = -1/2; params = list(x=x, c0=c0);
e = expression(x * log(x)^2 + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*d2y - x*dy + y - 2*x - c0 # = 0


# D =>
x*dy - y - 2*x*log(x) + c0 # = 0

# D2 =>
x*d2y - 2*log(x) - 2 # = 0
x^2*d2y - x*dy + y - 2*x - c0 # = 0


####################

### y = x * log(x)^2 + c2*x^2 + c1*x + c0;
# - more complicated F0(x);

# Check:
# for Quasi-Homogenous: c0 = 0;
# Note: c1 has NO impact;
x = sqrt(3); c2 = -1/3; c1 = -5/2; c0 = -1/2;
params = list(x=x, c2=c2, c1=c1, c0=c0);
e = expression(x * log(x)^2 + c2*x^2 + c1*x + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*d2y - x*dy + y - c2*x^2 - 2*x - c0 # = 0


# D =>
x*dy - y - 2*x*log(x) - c2*x^2 + c0 # = 0

# D2 =>
x*d2y - 2*log(x) - 2*c2*x - 2 # = 0
x^2*d2y - x*dy + y - c2*x^2 - 2*x - c0 # = 0


### Variant: with Coef c1;
x = sqrt(3); c2 = -1/3; c1 = -2/5; c0 = -1/2;
params = list(x=x, c2=c2, c1=c1, c0=c0);
e = expression(c1/2 * x * log(x)^2 + c2*x^2 + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### IDE:
x^2*d2y - x*dy + y - (c2*x^2 + c1*x + c0) # = 0


############################

### y = (x+b0) * log(x+b0)^2
# - just a shift;

# Check:
# for Homogenous: c0 = 0;
x = sqrt(3); b0 = 1; c0 = -1/2;
params = list(x=x, b0=b0, c0=c0);
e = expression((x+b0) * log(x+b0)^2 + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
(x+b0)^2*d2y - (x+b0)*dy + y - 2*(x+b0) - c0 # = 0


# D =>
dy - log(x+b0)^2 - 2*log(x+b0) # = 0
(x+b0)*dy - y - 2*(x+b0)*log(x+b0) + c0 # = 0

# D2 =>
(x+b0)*d2y - 2*log(x+b0) - 2 # = 0
(x+b0)^2*d2y - ((x+b0)*dy - y + c0) - 2*(x+b0) # = 0


### Variant:
e = expression((x+b0)/x * log(x+b0)^2 + c0)[[1]];
y   = eval(e, params); dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x*(x+b0)^2*d2y + (x+b0)*(x+2*b0)*dy - b0*(y - c0) - 2*(x+b0) # = 0


############

### y = x^p * log(x)^2
# - Simple with exponent;

# Check:
# for Quasi-Homogenous: c0 = 0;
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


#################

#################
### Components: 2

# y = b1(x) * log(P(x)) + B2(x) * Log(P(x))^2 + F0(x)

### y = x^2 * log(x)^2 + x*log(x)

# Check:
# for Quasi-Homogenous: c0 = 0;
x = sqrt(3); c0 = -1/2; params = list(x=x, c0=c0);
e = expression(x^2 * log(x)^2 + x*log(x) + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*(2*x-1)*d2y - 2*x*(3*x-1)*dy + 2*(4*x-1)*(y - c0) - 4*x^3 + 6*x^2 - x # = 0


# D =>
dy - 2*x*log(x)^2 - (2*x+1)*log(x) - 1 # = 0
x*dy - 2*(y - c0) - x*(2*x-1)*log(x) - x # = 0

# D2 =>
x*d2y - dy - (4*x-1)*log(x) - 2*x # = 0
x^2*(2*x-1)*d2y - x*(2*x-1)*dy - (4*x-1)*(x*dy - 2*(y - c0) - x) - 2*x^2*(2*x-1) # = 0
x^2*(2*x-1)*d2y - 2*x*(3*x-1)*dy + 2*(4*x-1)*(y - c0) - 4*x^3 + 6*x^2 - x # = 0


### Variant: y = c2*x^2 * log(x)^2 + x*log(x)

# Check:
# for Quasi-Homogenous: c0 = 0;
x = sqrt(3); c0 = 1/2; c2 = 1/sqrt(2);
# c2 = 1/2; # c2 = -1/2;
params = list(x=x, c0=c0, c2=c2);
e = expression(c2 * x^2 * log(x)^2 + x*log(x) + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*(2*c2*x-1)*d2y - 2*x*(3*c2*x-1)*dy + 2*(4*c2*x-1)*(y - c0) +
	- 4*c2^2*x^3 + 6*c2*x^2 - x # = 0


# D =>
dy - 2*c2*x*log(x)^2 - (2*c2*x+1)*log(x) - 1 # = 0
x*dy - 2*(y - c0) - x*(2*c2*x-1)*log(x) - x # = 0

# D2 =>
x*d2y - dy - (4*c2*x-1)*log(x) - 2*c2*x # = 0
x^2*(2*c2*x-1)*d2y - 2*x*(3*c2*x-1)*dy + 2*(4*c2*x-1)*(y - c0) +
	- 4*c2^2*x^3 + 6*c2*x^2 - x # = 0

### Special Cases:

### Case: c2 = 1/2
x^2*(x-1)*d2y - x*(3*x-2)*dy + 2*(2*x-1)*(y - c0) +
	- x^3 + 3*x^2 - x # = 0

### Case: c2 = -1/2
x^2*(x+1)*d2y - x*(3*x+2)*dy + 2*(2*x+1)*(y - c0) +
	+ x^3 + 3*x^2 + x # = 0


############
### Generic:

### y = B2(x) * log(P(x))^2 + B1(x)*log(P(x))

# Check:
# for Quasi-Homogenous: c0 = 0;
x = sqrt(3); c0 = -1/2; params = list(x=x, c0=c0);
eb2 = expression(x^2 + x/3 + 1/5)[[1]];
eb1 = expression(x^2 + x/5 + 2/3)[[1]];
# eb1 = eb2; # Case: b1 == b2;
epx = expression(x^2 + x/2 + 3/5)[[1]];
# eb1 = eb2; epx = expression(x)[[1]];
# eb1 = eb2 = expression(sqrt(x))[[1]]; epx = expression(x)[[1]];
# eb1 = eb2 = expression(x^(1/3))[[1]]; epx = expression(x)[[1]]; p = 1/3;
e = expression(b2 * log(px)^2 + b1 * log(px) + c0)[[1]];
e[[2]][[2]][[3]][[2]][[2]] = epx;
e[[2]][[3]][[3]][[2]] = epx;
e[[2]][[2]][[2]] = eb2;
e[[2]][[3]][[2]] = eb1;

#
y  = eval(e, params); dy = eval(D(e, "x"), params);
px = eval(epx, params); dp  = eval(D(epx, "x"), params);
b1 = eval(eb1, params); db1 = eval(D(eb1, "x"), params);
b2 = eval(eb2, params); db2 = eval(D(eb2, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);
d2p = eval(D(D(epx, "x"), "x"), params);
d2b1 = eval(D(D(eb1, "x"), "x"), params);
d2b2 = eval(D(D(eb2, "x"), "x"), params);

### ODE:
# - see formula below;


# D =>
px*dy - db2*px * log(px)^2 - (2*b2*dp + db1*px) * log(px) - b1*dp # = 0
b2*px * dy - db2*px * (y - c0) +
	- (b2*(2*b2*dp + db1*px) - b1*db2*px) * log(px) - b1*b2*dp # = 0

# D2 =>
b2*px * d2y + b2*dp * dy - (db2*dp + d2b2*px) * (y - c0) +
	+ (d2b2*b1*px - d2b1*b2*px + db2*(b1-4*b2)*dp - db1*b2*dp - 2*b2^2*d2p) * log(px) +
	- b1*b2*d2p - 2*db1*b2*dp - 2*b2^2*dp^2 / px # = 0
b2*px^2*((b2*db1 - b1*db2)*px + 2*b2^2*dp) * d2y +
	+ b2*px * ((b1*d2b2 - b2*d2b1)*px^2 - 4*b2*db2*px*dp + 2*b2^2*dp^2 - 2*b2^2*px*d2p) * dy +
	- px*(b2*(db1*d2b2 - db2*d2b1)*px^2 - 4*b2*db2^2*px*dp + 2*b2^2*d2b2*px*dp +
		+ 2*b2^2*db2*dp^2 - 2*b2^2*db2*px*d2p) * (y - c0) +
	- (2*db1*b2*(b2*db1 - b1*db2) + b1*b2*(b1 - b2)*d2b2) * px^2*dp +
	- b2*(b1 - 6*b2)*(b1*db2 - b2*db1) * px*dp^2 - 4*b2^4*dp^3 +
	- b1*b2*(b2*db1 - b1*db2) * px^2*d2p # = 0

# TODO:
# - gain insight from formula;

### Special Cases:

### Case: b1 = b2
b2^2*px^2*dp * d2y +
	- b2*px * (2*db2*px*dp - b2*dp^2 + b2*px*d2p) * dy +
	+ px * ((2*db2^2 - b2*d2b2)*px*dp - b2*db2*dp^2 + b2*db2*px*d2p) * (y - c0) +
	- 2*b2^3*dp^3 # = 0

### Case: b1 = b2; px = x;
b2^2*x^2 * d2y - b2*x * (2*db2*x - b2) * dy +
	+ x*((2*db2^2 - b2*d2b2)*x - b2*db2) * (y - c0) - 2*b2^3 # = 0

### Case: b1 = b2 = sqrt(x); px = x;
4*x^2*d2y + (y - c0) - 8*sqrt(x) # = 0
### Case: b1 = b2 = x^p; px = x;
x^2*d2y - (2*p-1)*x*dy + p^2*(y - c0) - 2*x^p # = 0


#######################
#######################

### Type: ATAN()^2
# y = B(x) * Atan(P(x))^2 + P0(x)

### y = x * atan(x)^2

# Check:
# for Quasi-Homogenous: c0 = 0;
x = sqrt(3); c0 = -1/2; params = list(x=x, c0=c0);
e = expression(x * atan(x)^2 + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*(x^2+1)^2*d2y - 2*x*(x^2+1)*dy + 2*(x^2+1)*y - 2*x^3 - 2*c0*(x^2+1) # = 0


# D =>
x*dy - y - 2*x^2*atan(x) / (x^2+1) + c0 # = 0
x*(x^2+1)*dy - (x^2+1)*y - 2*x^2*atan(x) + c0*(x^2+1) # = 0

# D2 =>
x*(x^2+1)*d2y + 2*x^2*dy - 2*x*y +
	- 4*x*atan(x) - 2*x^2/(x^2+1) + 2*c0*x # = 0
x^2*(x^2+1)*d2y - 2*x*dy + 2*y - 2*x^3/(x^2+1) - 2*c0 # = 0


### y = (x^2+1) * atan(x)^2

# Check:
# for Quasi-Homogenous: c0 = 0;
x = sqrt(3); c0 = -1/2; params = list(x=x, c0=c0);
e = expression((x^2+1) * atan(x)^2 + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
(x^2+1)^2*d2y - 2*x*(x^2+1)*dy + 2*(x^2-1)*y - 2*(c0 + 1)*x^2 + 2*(c0 - 1) # = 0


# D =>
(x^2+1)*dy - 2*x*y - 2*(x^2+1)*atan(x) + 2*x*c0 # = 0

# D2 =>
(x^2+1)*d2y - 2*y - 4*x*atan(x) + 2*c0 - 2 # = 0
(x^2+1)^2*d2y - 2*x*(x^2+1)*dy + 2*(x^2-1)*y - 2*(c0 + 1)*x^2 + 2*(c0 - 1) # = 0


#####################
#####################

### Simple Radical & Function

# y = B1(x) * R(x)^p * log(P(x))^2
# Linear ODE: Order 3 OR w. Radical;

### Simple Example
# y = sqrt(x^2 + b0) * log(x^2 + b0)^2

# Check:
# for Quasi-Homogenous: c0 = 0;
x = sqrt(3); b0 = 1/3; c0 = -1/2; params = list(x=x, c0=c0);
e = expression(sqrt(x^2+b0) * log(x^2+b0)^2 + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
# see below;


# D =>
(x^2+b0)*dy - x*(y - c0) - 4*x * sqrt(x^2+b0) * log(x^2+b0) # = 0

# D2 =>
x*(x^2+b0)^2*d2y - (x^2+b0)^2*dy + x^3 * (y - c0) +
	- 8*x^3 * sqrt(x^2+b0) # = 0


#####################
#####################

### w. Radicals

### Log & Radicals

### y = x * log(sqrt(x^2 + b0) - x)^2

# Check:
# for Homogenous: c0 = 0;
x = sqrt(3); b0 = 1/3; c0 = -1/2;
params = list(x=x, b0=b0, c0=c0);
e = expression(x * log(sqrt(x^2 + b0) - x)^2 + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*(x^2+b0)*d2y - x*(x^2+2*b0)*dy + (x^2+2*b0)*y +
	- 2*x^3 - c0*(x^2+2*b0) # = 0

# D =>
x*dy - (y - c0) - 2*x^2 * log(sqrt(x^2+b0) - x) * (x/sqrt(x^2+b0) - 1) / (sqrt(x^2+b0) - x) # = 0
b0*x*dy - b0*(y - c0) - 2*x^2 * log(sqrt(x^2+b0) - x) * (x/sqrt(x^2+b0) - 1) * (sqrt(x^2+b0) + x) # = 0
x*(x^2+b0)*dy - (x^2+b0)*(y - c0) + 2*x^2 * sqrt(x^2+b0) * log(sqrt(x^2+b0) - x) # = 0

# D2 =>
x*(x^2+b0)*d2y + 2*x^2*dy - 2*x*(y - c0) +
	+ (4*x + 2*x^3/(x^2+b0)) * sqrt(x^2+b0) * log(sqrt(x^2+b0) - x) - 2*x^2 # = 0
x^2*(x^2+b0)*d2y + 2*x^3*dy - 2*x^2*y +
	- (2 + x^2/(x^2+b0)) * (x*(x^2+b0)*dy - (x^2+b0)*(y - c0)) - 2*x^3 + 2*c0*x^2 # = 0
x^2*(x^2+b0) * d2y - x*(x^2+2*b0) * dy + (x^2+2*b0) * y +
	- 2*x^3 - c0*(x^2+2*b0) # = 0

####################

### y = x * log(sqrt(x^4 + b0) - x^2)^2
# - Effect of Power inside SQRT;

# Check:
# for Homogenous: c0 = 0;
x = sqrt(3); b0 = 1/3; c0 = -1/2;
params = list(x=x, b0=b0, c0=c0);
e = expression(x * log(sqrt(x^4 + b0) - x^2)^2 + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*(x^4+b0)*d2y - x*(x^4+3*b0)*dy + (x^4+3*b0)*y +
	- 8*x^5 - c0*(x^4+3*b0) # = 0

# D =>
x*dy - (y - c0) - 4*x^3 * log(sqrt(x^4+b0) - x^2) * (x^2/sqrt(x^4+b0) - 1) / (sqrt(x^4+b0) - x^2) # = 0
x*(x^4+b0)*dy - (x^4+b0)*(y - c0) + 4*x^3 * sqrt(x^4+b0) * log(sqrt(x^4+b0) - x^2) # = 0

# D2 =>
x*(x^4+b0)*d2y + 4*x^4*dy - 4*x^3*(y - c0) +
	+ (3 + 2*x^4/(x^4+b0)) * 4*x^2 * sqrt(x^4+b0) * log(sqrt(x^4+b0) - x^2) - 8*x^4 # = 0
x*(x^4+b0)*d2y + 4*x^4*dy - 4*x^3*(y - c0) +
	- (3*(x^4+b0) + 2*x^4) * (x*dy - (y - c0))/x - 8*x^4 # = 0


### Gen: y = x * log(sqrt(x^(2*n) + b0) - x^n)^2
# - Effect of Power inside SQRT;

# Check:
# for Homogenous: c0 = 0;
x = sqrt(3); n = 3/5;  b0 = 1/3; c0 = -1/2;
params = list(x=x, b0=b0, c0=c0);
e = expression(x * log(sqrt(x^(2*n) + b0) - x^n)^2 + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*(x^(2*n)+b0)*d2y - x*(x^(2*n)+(n+1)*b0)*dy + (x^(2*n)+(n+1)*b0)*y +
	- 2*n^2 * x^(2*n+1) - c0*(x^(2*n)+(n+1)*b0) # = 0


### y = x * log(sqrt((x + b1)^2 + b0) - x - b1)^2
# - Slightly more complicated SQRT;

# Check:
# for Quasi-Homogenous: c0 = 0;
x = sqrt(3); b0 = 1/3; b1 = 3/4; c0 = -1/2;
params = list(x=x, b0=b0, c0=c0);
e = expression(x * log(sqrt((x + b1)^2 + b0) - x - b1)^2 + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*((x + b1)^2 + b0)*d2y - x*((x + b1)*(x + 2*b1) + 2*b0)*dy +
	+ ((x + b1)*(x + 2*b1) + 2*b0)*y +
	- 2*x^3 - c0*((x + b1)*(x + 2*b1) + 2*b0) # = 0


# D =>
x*dy - (y - c0) - 2*x^2 * log(sqrt((x + b1)^2 + b0) - x - b1) *
	((x+b1)/sqrt((x + b1)^2 + b0) - 1) / (sqrt((x + b1)^2 + b0) - x - b1) # = 0
x*dy - (y - c0) + 2*x^2 * log(sqrt((x + b1)^2 + b0) - x - b1) / sqrt((x + b1)^2 + b0) # = 0

# D2 =>
x*d2y + 4*x * log(sqrt((x + b1)^2 + b0) - x - b1) / sqrt((x + b1)^2 + b0) +
	- 2*x^2 * (x + b1) * log(sqrt((x + b1)^2 + b0) - x - b1) / sqrt((x + b1)^2 + b0)^3 +
	- 2*x^2 / ((x + b1)^2 + b0) # = 0
x^2*((x + b1)^2 + b0)*d2y - x*((x + b1)*(x + 2*b1) + 2*b0)*dy +
	+ ((x + b1)*(x + 2*b1) + 2*b0)*y +
	- 2*x^3 - c0*((x + b1)*(x + 2*b1) + 2*b0) # = 0


### Variant: b0 < 0;
# - same ODE: see above;
# for Quasi-Homogenous: c0 = 0;
x = sqrt(3); b0 = -1/5; b1 = 3/4; c0 = -1/2;
params = list(x=x, b0=b0, c0=c0);
e = expression(x * log(x + b1 - sqrt((x + b1)^2 + b0))^2 + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);



######################

### ATAN( Radicals )

### y = x * atan(sqrt(x+b0))^2

# Check:
# for Quasi-Homogenous: c0 = 0;
x = sqrt(3); c0 = -1/2;
b0 = 1/3; # b0 = -1;
params = list(x=x, b0=b0, c0=c0);
e = expression(x * atan(sqrt(x+b0))^2 + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
# see below;

# D =>
x*dy - y - x^2 * atan(sqrt(x+b0)) / sqrt(x+b0) / (x+b0+1) + c0 # = 0
x*(x+b0+1)*dy - (x+b0+1)*y - x^2 * atan(sqrt(x+b0)) / sqrt(x+b0) + c0*(x+b0+1) # = 0

# D2 =>
x*(x+b0+1)*d2y + x*dy - y +
	- 2*x * atan(sqrt(x+b0)) / sqrt(x+b0) +
	- 1/2 * x^2 / (x+b0) / (x+b0+1) +
	+ 1/2 * x^2 * atan(sqrt(x+b0)) / sqrt(x+b0)^3 + c0 # = 0
2*x^2*(x+b0+1)*d2y + 2*x^2*dy - 2*x*y +
	- 4*(x*(x+b0+1)*dy - (x+b0+1)*y + c0*(x+b0+1)) +
	- x^3 / (x+b0) / (x+b0+1) +
	+ x / (x+b0) * (x*(x+b0+1)*dy - (x+b0+1)*y + c0*(x+b0+1)) + 2*c0*x # = 0
2*x^2*(x+b0)*(x+b0+1)*d2y +
	- x * (x^2 + 5*b0*x + 3*x + 4*b0 + 4*b0^2) * dy +
	+ (x^2 + 3*x + 5*b0*x + 4*b0 + 4*b0^2) * y +
	- x^3 / (x+b0+1) - c0*(x^2 + 3*x + 5*b0*x + 4*b0 + 4*b0^2) # = 0
2*x^2*(x+b0)*(x+b0+1)^2*d2y +
	- x*(x+b0+1) * (x^2 + 5*b0*x + 3*x + 4*b0 + 4*b0^2) * dy +
	+ (x+b0+1) * (x^2 + 3*x + 5*b0*x + 4*b0 + 4*b0^2) * y +
	- c0*(x+b0+1) * (x^2 + 3*x + 5*b0*x + 4*b0 + 4*b0^2) - x^3 # = 0


# Special Cases:

# Case: b0 = -1
2*x^2*(x-1)*d2y - x*(x-2)*dy + (x-2)*y - (c0+1)*x + 2*c0 # = 0

#################

#################
### Components: 2

### y = x * atan(sqrt(x+b0))^2 + c1 * sqrt(x+b0) * atan(sqrt(x+b0))

# Check:
# for Quasi-Homogenous: c0 = 0;
x = sqrt(3); c1 = 2/5; c0 = -1/2; # c1 = 2;
b0 = 1/3; # b0 = -1; c1 = 2;
params = list(x=x, b0=b0, c0=c0);
e = expression(x * atan(sqrt(x+b0))^2 + c1 * sqrt(x+b0) * atan(sqrt(x+b0)) + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
# see below;

# D =>
(x+b0)*(x+b0+1)*dy - (x+b0)*(x+b0+1)*atan(sqrt(x+b0))^2 - x*sqrt(x+b0)*atan(sqrt(x+b0)) +
	- c1/2 * (x+b0+1) * sqrt(x+b0) * atan(sqrt(x+b0)) - c1/2 * (x+b0) # = 0
2*x*(x+b0)*(x+b0+1)*dy - 2*(x+b0)*(x+b0+1)*(y - c0) +
	- (2*x^2 - c1*(x+2*b0)*(x+b0+1)) * sqrt(x+b0) * atan(sqrt(x+b0)) - c1*x*(x+b0) # = 0

# D2 =>
4*x*(x+b0)^2*(x+b0+1) * d2y + 4*x*(x+b0)*(2*x+2*b0+1) * dy +
	- 4*(x+b0)*(2*x+2*b0+1) * (y - c0) +
	- (2*x*(5*x+4*b0) - c1*x*(5*x+7*b0+3) - c1*b0*(6*x+8*b0+4)) * sqrt(x+b0) * atan(sqrt(x+b0)) +
	- 2*x^2 * (x+b0) / (x+b0 + 1) - 3*c1*x*(x+b0) # = 0
# Substitution:
2*x*(x+b0)^2*(x+b0+1)*(2*x^2 - c1*(x+2*b0)*(x+b0+1)) * d2y +
	+ x*(x+b0) * ((c1-2)*x^3 + 2*b0*c1*x^2 - 10*b0*x^2 + 2*c1*x^2 - 6*x^2 +
		- 8*b0*x - 8*b0^2*x + c1*x + 2*b0*c1*x + b0^2*c1*x) * dy +
	- (x+b0) * ((c1-2)*x^3 + (2*b0*c1 - 10*b0 + 2*c1 - 6)*x^2 +
		+ b0^2*c1*x + 2*b0*c1*x - 8*b0^2*x - 8*b0*x + c1*x) * (y - c0) +
	- 2*x^4 * (x+b0) / (x+b0 + 1) +
	- c1*x * (x+b0) * ((c1 - 3)*x^2 + 2*b0*(c1-3)*x + b0*c1*(b0-1)) # = 0


### Special Cases:

### Case: c1 = 2
2*(x+b0)^2*(x+b0+1)^2 * ((3*b0+1)*x + 2*b0^2 + 2*b0) * d2y +
	+ x*(x+b0)*(x+b0+1) * ((3*b0 + 1)*x + 3*b0^2 + 2*b0 - 1) * dy +
	- (x+b0)*(x+b0+1) * ((3*b0 + 1)*x + 3*b0^2 + 2*b0 - 1) * (y - c0) +
	- (x+b0) * ((3*b0 + 1)*x^2 + 4*b0*x - 2*b0^3 + 2*b0) # = 0

### Case: c1 = 2; b1 = -1;
2*x^2*(x-1)*d2y + x^2*dy - x*(y - c0) - (x + 2) # = 0

