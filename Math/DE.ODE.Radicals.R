########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Linear ODEs - Radicals
##
## draft v.0.1d


### Theory

# y = B1(x) * R(x)^(1/n) * FUN(P(x)) + B2(x) * R(x)^(1/n);
# - where FUN = Log, Atan;
# - the Radical in the first term entangles with the free radical;
# y = B1(x) * R(x)^(1/n) * SIN(P(x)) + B2(x) * R(x)^(1/n) * COS(P(x));
# - the coefficients of SIN & COS entangle with each-other;


### Examples:

# x*(x^2+b0)^2 * d2y - b0*(x^2+b0) * dy + x^3 * (y - c0) = 0;


####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### Radicals w. Log
# From: y = B1(x) * R(x)^(1/n) * log(P(x)) + B2(x) * R(x)^(1/n);


### y = x^2 * (x^2+b0)^(1/3) * log(x) + x * (x^2+b0)^(1/3);

# Check:
x = sqrt(3); b0 = sqrt(2); c0 = -1/3; params = list(x=x, b0=b0, c0=c0);
e = expression(x^2 * (x^2+b0)^(1/3) * log(x) + x * (x^2+b0)^(1/3) + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
9*x^2*(x-1)*(x^2+b0)^2 * d2y +
	- 3*x*(x^2+b0)*(13*x^3 - 10*x^2 + 9*b0*x - 6*b0) * dy +
	+ (64*x^5 - 40*x^4 + 84*b0*x^3 - 42*b0*x^2 + 36*b0^2*x - 18*b0^2) * (y - c0)  # = 0


# D =>
3*(x^2+b0)*dy - (8*x^2 + 6*b0) * x * (x^2+b0)^(1/3) * log(x) +
	- (3*x^3 + 5*x^2 + 3*b0*x + 3*b0) * (x^2+b0)^(1/3) # = 0
3*x*(x^2+b0)*dy - (8*x^2 + 6*b0) * (y - c0) +
	- 3*x*(x^3 - x^2 + b0*x - b0) * (x^2+b0)^(1/3) # = 0

# D2 =>
3*x*(x^2+b0)*d2y + (x^2-3*b0)*dy - 16*x * (y - c0) +
	- (14*x^3 - 11*x^2 + 6*b0*x - 3*b0) * (x^2+b0)^(1/3) # = 0
9*x^2*(x-1)*(x^2+b0)^2*d2y + 3*x*(x-1)*(x^2+b0)*(x^2-3*b0)*dy +
	- 48*x^2*(x-1)*(x^2+b0) * (y - c0) +
	- (14*x^3 - 11*x^2 + 6*b0*x - 3*b0) * (3*x*(x^2+b0)*dy - (8*x^2 + 6*b0) * (y - c0)) # = 0
9*x^2*(x-1)*(x^2+b0)^2 * d2y +
	- 3*x*(x^2+b0)*(13*x^3 - 10*x^2 + 9*b0*x - 6*b0) * dy +
	+ (64*x^5 - 40*x^4 + 84*b0*x^3 - 42*b0*x^2 + 36*b0^2*x - 18*b0^2) * (y - c0)  # = 0


### Special Cases:

### Case: b0 = -1
x = sqrt(3); b0 = -1; c0 = -1/3; params = list(x=x, b0=b0, c0=c0);
e = expression(x^2 * (x^2+b0)^(1/3) * log(x) + x * (x^2+b0)^(1/3) + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
9*x^2*(x-1)^2*(x+1)^2 * d2y +
	- 3*x*(x-1)*(x+1)*(13*x^2 + 3*x - 6) * dy +
	+ (64*x^4 + 24*x^3 - 60*x^2 - 18*x + 18) * (y - c0)  # = 0

###################
###################

### w. Trig

### y = x * (x^2+b0)^(1/3) * sin(k*x^2);

# Check:
x = sqrt(3); k = 2/5; b0 = sqrt(2); c0 = -1/3;
params = list(x=x, k=k, b0=b0, c0=c0);
e = expression(x * (x^2+b0)^(1/3) * sin(k*x^2) + c0)[[1]];
#
y   = eval(e, params); dy = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
9*x^2*(x^2+b0)^2*d2y - 3*x*(x^2+b0)*(13*x^2 + 9*b0) * dy +
	+ (36*k^2*x^8 + 72*k^2*b0*x^6 + (36*k^2*b0^2+55)*x^4 +
		+ 66*b0*x^2 + 27*b0^2) * (y - c0) # = 0


# D =>
3*x*(x^2+b0)*dy - (5*x^2+3*b0) * (y - c0) - 6*k*x^3*(x^2+b0) * (x^2+b0)^(1/3) * cos(k*x^2) # = 0

# D2 =>
9*x^2*(x^2+b0)^2*d2y - 3*x*(x^2+b0)*(13*x^2 + 9*b0) * dy +
	+ (36*k^2*x^8 + 72*k^2*b0*x^6 + (36*k^2*b0^2+55)*x^4 + 66*b0*x^2 + 27*b0^2) * (y - c0) # = 0


#################

### Type: Sin(Log(P(x)))
# y = (x^2+b0)^(1/3) * sin(k*log(x));

# Check:
x = sqrt(3); k = 2/5; b0 = sqrt(2); c0 = -1/3;
params = list(x=x, k=k, b0=b0, c0=c0);
e = expression((x^2+b0)^(1/3) * sin(k*log(x)) + c0)[[1]];
#
y   = eval(e, params); dy = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
9*x^2*(x^2+b0)^2 * d2y - 3*x*(x^2+b0)*(x^2-3*b0) * dy +
	+ (9*k^2*(x^2+b0)^2 + 4*x^4 - 12*b0*x^2) * (y - c0) # = 0


# D =>
3*x*(x^2+b0)*dy - 2*x^2 * (y - c0) +
	- 3*k*(x^2+b0) * (x^2+b0)^(1/3) * cos(k*log(x)) # = 0

# D2 =>
9*x^2*(x^2+b0)^2*d2y - 3*x*(x^2+b0)*(x^2-3*b0)*dy +
	+ (9*k^2*(x^2+b0)^2 + 4*x^4 - 12*b0*x^2) * (y - c0) # = 0


#################

### y = (x^2+b0)^p * sin(k*log(x^2+b0));

# Check:
x = sqrt(3); p = 2/3; k = 2/5; b0 = sqrt(2); c0 = -1/3;
# p = 1/4; k = sqrt(3)/4;
params = list(x=x, k=k, b0=b0, c0=c0);
e = expression((x^2+b0)^p * sin(k*log(x^2+b0)) + c0)[[1]];
#
y   = eval(e, params); dy = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x*(x^2+b0)^2 * d2y - ((4*p-1)*x^2 + b0)*(x^2+b0) * dy +
	+ 4*(k^2 + p^2)*x^3 * (y - c0) # = 0

# D =>
(x^2+b0)*dy - 2*p*x * (y - c0) - 2*k*x * (x^2+b0)^p * cos(k*log(x^2+b0)) # = 0

# D2 =>
(x^2+b0)^2 * d2y - 2*(p-1)*x*(x^2+b0) * dy +
	+ (4*k^2*x^2 - 2*p*(x^2+b0)) * (y - c0) +
	- 2*k*(2*p*x^2 + x^2+b0) * (x^2+b0)^p * cos(k*log(x^2+b0)) # = 0
x*(x^2+b0)^2 * d2y - ((4*p-1)*x^2 + b0)*(x^2+b0) * dy +
	+ 4*(k^2 + p^2)*x^3 * (y - c0) # = 0


### Special Cases:

### Case: p = 1/4; k = sqrt(3)/4;
x*(x^2+b0)^2 * d2y - b0*(x^2+b0) * dy + x^3 * (y - c0) # = 0

