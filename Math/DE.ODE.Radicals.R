########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Linear ODEs - Radicals
##
## draft v.0.1k


### Theory

### Simple Radicals:
# y = B1(x) * R1(x)^p1 + B2(x) * R2(x)^p2;

### Non-Trig:
# y = B1(x) * R(x)^p * FUN(P(x)) + B2(x) * R(x)^p;
# - where FUN = Log, Atan;
# - the Radical in the first term entangles with the free radical;
# - Power p = (usually) non-integer value;

### w. Trig:
# y = B1(x) * R(x)^p * SIN(P(x)) + B2(x) * R(x)^p * COS(P(x));
# y = B1(x) * R(x)^p * SIN(k*log(P(x))) + B2(x) * R(x)^p * COS(k*log(P(x)));
# - the coefficients of SIN & COS entangle with each-other;
# - Power p = non-integer value;

### Nested-Radicals
# y = B(x) * sqrt(sqrt(P(x)^2 + b0) - P(x));

# Note:
# - Radicals persist across all levels of differentiation:
# => Homogenous ODE;


### Examples:

# 16*x*(x+b0) * d2y + 8*(2*x + b0) * dy - (y - c0) = 0;
# 16*x^2*(x+b0) * d2y + 8*x*(3*x + 2*b0) * dy - b0 * (y - c0) = 0;

# Gen: (x^2+b0) * d2y + x * dy - k^2 * (y - c0) = 0
# Special Cases & Variants:
# 4*(x^2+b0) * d2y + 4*x * dy - (y - c0) = 0;
# x*(x^2+b0)^2 * d2y - b0*(x^2+b0) * dy + x^3 * (y - c0) = 0;
# 4*(x^2+b0)^2 * d2y - 4*x*(x^2+b0) * dy + (3*x^2 - 5*b0) * (y - c0) = 0;


####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### Simple Radicals

# y = B1(x) * R1(x)^k1 + B2(x) * R2(x)^k2

### y = (x2+b0)^k1 + c1 * (x^2+b0)^k2

# Check:
x = sqrt(3); k1 = 1/3; k2 = 3/5; b0 = sqrt(2); c1 = 2/5; c0 = -1/3;
# k2 = 1 - k1; # k1 = (1-1i*sqrt(3))/2; k2 = 1 - k1;
params = list(x=x, k1=k1, k2=k2, b0=b0, c1=c1, c0=c0);
e = expression((x^2+b0)^k1 + c1 * (x^2+b0)^k2 + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x*(x^2+b0)^2*d2y +
	- (x^2+b0)*(2*(k1+k2)*x^2 - x^2 + b0) * dy +
	+ 4*k1*k2*x^3 * (y - c0) # = 0


# D =>
(x^2+b0)*dy - 2*k1*x*(x^2+b0)^k1 - 2*c1*k2*x*(x^2+b0)^k2 # = 0

# Isolation of Radicals:
2*(k1-k2)*x*(x^2+b0)^k1 # ==
(x^2+b0)*dy - 2*k2*x*(y - c0);
#
2*c1*(k1-k2)*x*(x^2+b0)^k2 # ==
- ((x^2+b0)*dy - 2*k1*x*(y - c0));

# D2 =>
(x^2+b0)^2*d2y + 2*x*(x^2+b0)*dy +
	- 2*k1*(2*k1*x^2 + x^2 + b0) * (x^2+b0)^k1 +
	- 2*c1*k2*(2*k2*x^2 + x^2 + b0) * (x^2+b0)^k2 # = 0
(k1-k2)*x*(x^2+b0)^2*d2y +
	- (k1-k2)*(x^2+b0)*(2*(k1+k2)*x^2 - x^2 + b0) * dy +
	+ 4*k1*k2*(k1-k2)*x^3 * (y - c0) # = 0

### Special Cases:

### Case: k1 + k2 = 1
x*(x^2+b0)^2 * d2y - (x^2+b0)^2 * dy - 4*k1*(k1-1)*x^3 * (y - c0) # = 0
### Case: k1*(k1-1) = -1;
x*(x^2+b0)^2 * d2y - (x^2+b0)^2 * dy + 4*x^3 * (y - c0) # = 0


#############

### Radicals:
### Polynomial Fractions

### P(x): Order 1
# y = (1+d0/(x+b0))^k1 + c1 * (1+d0/(x+b0))^k2

# Check:
x = sqrt(5); k1 = 1/3; k2 = 3/5; b0 = sqrt(2); d0 = sqrt(3); c1 = 2/5; c0 = -1/3;
# k2 = 1 - k1;
params = list(x=x, k1=k1, k2=k2, b0=b0, c1=c1, c0=c0, d0=d0);
e = expression((1+d0/(x+b0))^k1 + c1 * (1+d0/(x+b0))^k2 + c0)[[1]];
#
y   = eval(e, params); dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
(x+b0)^2*(x+b0+d0)^2 * d2y +
	+ (x+b0)*(x+b0+d0)*(2*(x+b0+d0) + d0*(k1+k2-1)) * dy +
	+ k1*k2*d0^2 * (y - c0) # = 0

# D =>
(x+b0)*(x+b0+d0)*dy + k1*d0*(1+d0/(x+b0))^k1 +
	+ c1*k2*d0*(1+d0/(x+b0))^k2 # = 0

# Isolation of Radicals:
(k1-k2)*d0*(1+d0/(x+b0))^k1 # ==
- ((x+b0)*(x+b0+d0)*dy + k2*d0*(y - c0));
#
(k1-k2)*c1*d0*(1+d0/(x+b0))^k2 # ==
((x+b0)*(x+b0+d0)*dy + k1*d0*(y - c0));

# D2 =>
(x+b0)^2*(x+b0+d0)^2*d2y + (x+b0)*(x+b0+d0)*(2*x+2*b0+d0)*dy +
	- k1^2*d0^2 * (1+d0/(x+b0))^k1 +
	- c1*k2^2*d0^2 * (1+d0/(x+b0))^k2 # = 0
(k1-k2)*(x+b0)^2*(x+b0+d0)^2*d2y +
	+ (k1-k2)*(x+b0)*(x+b0+d0)*(2*x+2*b0+d0)*dy +
	+ (k1^2-k2^2)*d0*(x+b0)*(x+b0+d0)*dy +
	+ k1*k2*(k1-k2)*d0^2 * (y - c0) # = 0

### Special Cases:

### Case: k1 + k2 = 1
(x+b0)^2*(x+b0+d0)^2 * d2y + 2*(x+b0)*(x+b0+d0)^2 * dy +
	+ k1*k2*d0^2 * (y - c0) # = 0
# Reformulated:
# d0 = d - b; # Solve: k1 + k2 = 1; k1*k2 = k / (d-b)^2;
b = b0; d = b0+d0; k = k1*k2*d0^2;
(x+b)^2*(x+d)^2 * d2y + 2*(x+b)*(x+d)^2 * dy + k * (y - c0) # = 0


#################
### P(x): Order 2
# y = ((x^2+b0)/(x+d0))^k1 + c1 * ((x^2+b0)/(x+d0))^k2

# Check:
x = sqrt(5); k1 = 1/3; k2 = 3/5; b0 = sqrt(2); d0 = sqrt(3); c1 = 2/5; c0 = -1/3;
# b0 = - d0^2; # Trivial
# k = 3; k1 = (1 - sqrt(1 - 4*k + 0i))/2; k2 = 1 - k1;
params = list(x=x, b0=b0, c1=c1, c0=c0, d0=d0);
e = expression(((x^2+b0)/(x+d0))^k1 + c1 * ((x^2+b0)/(x+d0))^k2 + c0)[[1]];
#
y   = eval(e, params); dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
(x+d0)^2*(x^2+b0)^2*(x^2 + 2*d0*x - b0) * d2y +
	- (x+d0)*(x^2+b0) *
		((k1+k2-1)*(x^2 + 2*d0*x - b0)^2 + 2*(d0^2 + b0)*(x^2+b0)) * dy +
	+ k1*k2*(x^2 + 2*d0*x - b0)^3 * (y - c0) # = 0

# D =>
(x+d0)*(x^2+b0)*dy +
	- k1 * (x^2 + 2*d0*x - b0) * ((x^2+b0)/(x+d0))^k1 +
	- c1*k2 * (x^2 + 2*d0*x - b0) * ((x^2+b0)/(x+d0))^k2 # = 0

# Isolation of Radicals:
(k1-k2) * (x^2 + 2*d0*x - b0) * ((x^2+b0)/(x+d0))^k1 # ==
((x+d0)*(x^2+b0)*dy - k2*(x^2 + 2*d0*x - b0) * (y - c0))
#
c1*(k1-k2) * (x^2 + 2*d0*x - b0) * ((x^2+b0)/(x+d0))^k2 # ==
- ((x+d0)*(x^2+b0)*dy - k1*(x^2 + 2*d0*x - b0) * (y - c0))

# D2 =>
(x+d0)^2*(x^2+b0)^2 * d2y +
	+ (x+d0)*(x^2+b0)*(3*x^2 + 2*d0*x + b0)*dy +
	- k1^2 * (x^2 + 2*d0*x - b0)^2 * ((x^2+b0)/(x+d0))^k1 +
	- 2*k1 * (x + d0)^2*(x^2+b0) * ((x^2+b0)/(x+d0))^k1 +
	- c1*k2^2 * (x^2 + 2*d0*x - b0)^2 * ((x^2+b0)/(x+d0))^k2 +
	- 2*c1*k2 * (x + d0)^2*(x^2+b0) * ((x^2+b0)/(x+d0))^k2 # = 0
(k1-k2)*(x+d0)^2*(x^2+b0)^2*(x^2 + 2*d0*x - b0) * d2y +
	- (k1-k2)*(x+d0)*(x^2+b0) *
		((k1+k2-1)*(x^2 + 2*d0*x - b0)^2 + 2*(d0^2 + b0)*(x^2+b0)) * dy +
	+ k1*k2*(k1-k2) * (x^2 + 2*d0*x - b0)^3 * (y - c0) # = 0

### Special Cases:

### Trivial Case: b0 = -d0^2
(x-d0)^2 * d2y - (k1+k2-1)*(x-d0) * dy + k1*k2 * (y - c0) # = 0

### Case: k1 + k2 = 1
k = 3; k1 = (1 - sqrt(1 - 4*k + 0i))/2; k2 = 1 - k1;
(x+d0)^2*(x^2+b0)^2*(x^2 + 2*d0*x - b0) * d2y +
	- 2*(d0^2 + b0)*(x+d0)*(x^2+b0)^2 * dy +
	+ k*(x^2 + 2*d0*x - b0)^3 * (y - c0) # = 0


########################
########################

### Radicals & Functions

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


############
### Radicals w. Simple Log
# From: y = B2(x) * R(x)^k + B1(x) * log(P(x));


### y = (x^2+b0)^k + c1*x * log(x);

# Check:
x = sqrt(3); k = 1/3; b0 = sqrt(2); c1 = 2/5; c0 = -1/3;
# k = 1/2; # k = 1/4;
params = list(x=x, k=k, b0=b0, c1=c1, c0=c0);
e = expression((x^2+b0)^k + c1*x*log(x) + c0)[[1]];
#
y   = eval(e, params); dy = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x*(x^2+b0)*((2*k-1)*x^2-b0) * d2y +
	- 2*k*x^2 * ((2*k-1)*x^2+b0) * dy +
	+ 2*k*x * ((2*k-1)*x^2 + b0) * (y - c0) +
	+ c1 * ((2*k-1)^2*x^4 + 2*b0*x^2 + b0^2) # = 0

# D =>
(x^2+b0)*dy - 2*k*x*(y - c0) + c1*((2*k-1)*x^2-b0)*log(x) - c1*(x^2+b0) # = 0

# D2 =>
x*(x^2+b0) * d2y - 2*(k-1)*x^2 * dy - 2*k*x * (y - c0) +
	+ 2*(2*k-1)*x^2 * c1*log(x) +
	+ c1*((2*k-3)*x^2-b0) # = 0


### Special Cases:

### Case: k = 1/2
x*(x^2+b0) * d2y + 2*k*x^2 * dy - 2*k*x * (y - c0) +
	- c1 * (2*x^2 + b0) # = 0

### Case: k = 1/4
2*x*(x^2+b0)*(x^2 + 2*b0) * d2y +
	- x^2*(x^2 - 2*b0) * dy + x*(x^2 - 2*b0) * (y - c0) +
	- c1 * (x^4 + 8*b0*x^2 + 4*b0^2) # = 0


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
params = list(x=x, k=k, p=p, b0=b0, c0=c0);
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

### Case: p = 1/4;
k0 = -1;
k = (k0/4 - 1/4^2); k = ifelse(k >= 0, sqrt(k), sqrt(k + 0i));
params = list(x=x, k=k, p=1/4, b0=b0, c0=c0);
# Note: must re-initialize above;
x*(x^2+b0)^2 * d2y - b0*(x^2+b0) * dy + k0*x^3 * (y - c0) # = 0


#########################
#########################

### Type: Nested Radicals

# Simple: 1 Component
# y = B(x) * (sqrt(P(x)^2+b0) - P(x))^k or
# y = B(x) * sqrt(P(x)^2+b0) * (sqrt(P(x)^2+b0) - P(x))^k
# where k = Rational / Real / Complex;

### Example:
# y = sqrt(sqrt(x^2+b0) - x)

# Check:
x = sqrt(3); b0 = sqrt(2); c0 = -1/3;
# p = 1/4; k = sqrt(3)/4;
params = list(x=x, b0=b0, c0=c0);
e = expression(sqrt(sqrt(x^2+b0) - x) + c0)[[1]];
#
y   = eval(e, params); dy = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
4*(x^2+b0)*d2y + 4*x*dy - (y - c0) # = 0


# D =>
2*(x^2+b0)*dy + sqrt(x^2+b0) * sqrt(sqrt(x^2+b0) - x) # = 0

# D2 =>
4*(x^2+b0)*d2y + 8*x*dy +
	+ 2*x / (x^2+b0) * sqrt(x^2+b0) * sqrt(sqrt(x^2+b0) - x) +
	- sqrt(sqrt(x^2+b0) - x) # = 0
4*(x^2+b0)*d2y + 4*x*dy - (y - c0) # = 0


### Variant:

### y = sqrt(x^2+b0) * sqrt(sqrt(x^2+b0) - x)

# Check:
x = sqrt(3); b0 = sqrt(2); c0 = -1/3;
# p = 1/4; k = sqrt(3)/4;
params = list(x=x, k=k, p=p, b0=b0, c0=c0);
e = expression(sqrt(x^2+b0) * sqrt(sqrt(x^2+b0) - x) + c0)[[1]];
#
y   = eval(e, params); dy = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
4*(x^2+b0)^2*d2y - 4*x*(x^2+b0)*dy + (3*x^2 - 5*b0)*(y - c0) # = 0


# D =>
2*(x^2+b0)*dy - 2*x * (y - c0) + (x^2+b0) * sqrt(sqrt(x^2+b0) - x) # = 0

# D2 =>
4*(x^2+b0)*d2y + 4*x*dy - 5*(y - c0) +
	+ 4*x * sqrt(sqrt(x^2+b0) - x) # = 0
4*(x^2+b0)^2*d2y - 4*x*(x^2+b0)*dy + (3*x^2 - 5*b0)*(y - c0) # = 0


### Generalisation of Power

### y = sqrt(sqrt(x^(2*n)+b0) - x^n)

# Check:
x = sqrt(3); n = 2/5; b0 = sqrt(2); c0 = -1/3;
# n = 1/2; # n = -1/2;
params = list(x=x, k=k, n=n, b0=b0, c0=c0);
e = expression(sqrt(sqrt(x^(2*n)+b0) - x^n) + c0)[[1]];
#
y   = eval(e, params); dy = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
4*x*(x^(2*n)+b0) * d2y + 4*(x^(2*n) - (n-1)*b0) * dy +
	- n^2*x^(2*n-1) * (y - c0) # = 0

# D =>
2*(x^(2*n)+b0)*dy + n*x^(n-1) * sqrt(x^(2*n)+b0) * sqrt(sqrt(x^(2*n)+b0) - x^n) # = 0

# D2 =>
4*x*(x^(2*n)+b0) * d2y + 4*(x^(2*n) - (n-1)*b0) * dy +
	- n^2*x^(2*n-1) * (y - c0) # = 0

### Special Cases:

### Case: n = 1/2
16*x*(x+b0) * d2y + 8*(2*x + b0) * dy - (y - c0) # = 0

### Case: n = -1/2
bi = 1/b0;
16*x^2*(x+bi) * d2y + 8*x*(3*x + 2*bi) * dy - bi * (y - c0) # = 0


##################

### Outer Radical: Higher Power
# y = (sqrt(x^2+b0) - x)^k

# Check:
x = sqrt(3); k = 1/3; b0 = sqrt(2); c0 = -1/3;
# k = 1i/3;
params = list(x=x, k=k, b0=b0, c0=c0);
e = expression((sqrt(x^2+b0) - x)^k + c0)[[1]];
#
y   = eval(e, params); dy = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
(x^2+b0)*d2y + x*dy - k^2 * (y - c0) # = 0


# D =>
dy - k * (x/sqrt(x^2+b0) - 1) * (sqrt(x^2+b0) - x)^k / (sqrt(x^2+b0) - x) # = 0
(x^2+b0)*dy + k * sqrt(x^2+b0) * (sqrt(x^2+b0) - x)^k # = 0

# D2 =>
(x^2+b0)*d2y + 2*x*dy +
	+ k*x * sqrt(x^2+b0) / (x^2+b0) * (sqrt(x^2+b0) - x)^k +
	- k^2 * (sqrt(x^2+b0) - x)^k # = 0
(x^2+b0)*d2y + x*dy - k^2 * (y - c0) # = 0


##############
### Gen: Power
# y = (sqrt(x^(2*n)+b0) - x^n)^k

# Check:
x = sqrt(3); n = 2/5; k = 1/3; b0 = sqrt(2); c0 = -1/3;
# n = 2; k = 1i/2;
params = list(x=x, n=n, k=k, b0=b0, c0=c0);
e = expression((sqrt(x^(2*n)+b0) - x^n)^k + c0)[[1]];
#
y   = eval(e, params); dy = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x*(x^(2*n) + b0) * d2y + (x^(2*n) - (n-1)*b0) * dy +
	- k^2*n^2*x^(2*n-1) * (y - c0) # = 0


# D =>
(x^(2*n)+b0)*dy + k*n*x^(n-1) * sqrt(x^(2*n)+b0) * (sqrt(x^(2*n)+b0) - x^n)^k # = 0

# D2 =>
x*(x^(2*n) + b0) * d2y + (x^(2*n) - (n-1)*b0) * dy +
	- k^2*n^2*x^(2*n-1) * (y - c0) # = 0

### Special Cases:

### Case: n = 2
x*(x^4 + b0) * d2y + (x^4 - b0) * dy - 4*k^2*x^3 * (y - c0) # = 0

