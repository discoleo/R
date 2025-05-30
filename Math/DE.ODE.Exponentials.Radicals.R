########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs - Exponentials w. Radicals
##
## draft v.0.1f


### Types:

# 1. Simple SQRT:
#    y = G(x) * exp(B(x) * sqrt(P(x)) + B0(x));
# 2. Double SQRT:
#    y = G(x) * sqrt(P(x)) * exp(B(x) * sqrt(P(x) + B0(x));
# Note:
# - sqrt(P(x)) is the same with the one in the exponential;
# - P, B, G: polynomials (or polynomial fractions);


####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### Simple: SQRT

### y = exp(x + k*sqrt(x + b0))

x = sqrt(3); k = -1/5; b0 = 2/3; params = list(x=x, k=k, b0=b0);
e = expression(exp(x + k*sqrt(x + b0)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
2*dy - (2 + k/sqrt(x + b0)) * y # = 0

# D2 =>
4*d2y - 2*(2 + k/sqrt(x + b0)) * dy +
	+ k/(x + b0) / sqrt(x + b0) * y # = 0
4*d2y - 8*dy + (4 - k^2/(x+b0)) * y +
	+ k/(x + b0) / sqrt(x + b0) * y # = 0

### ODE:
4*(x + b0)*d2y - (8*x + 8*b0 - 2)*dy +
	+ (4*x - k^2 + 4*b0 - 2) * y # = 0

### Special Cases:

### k = sqrt(4*b0-2)
x = sqrt(3); b0 = 2/3; k = sqrt(4*b0-2);
params = list(x=x, k=k, b0=b0);
# Re-run assignments from above;
4*(x + b0)*d2y - (8*x + 8*b0 - 2)*dy + 4*x*y # = 0


### Variant:
### y = x * exp(x + k*sqrt(x + b0))

x = sqrt(3); k = -1/5; b0 = 2/3; params = list(x=x, k=k, b0=b0);
e = expression(x * exp(x + k*sqrt(x + b0)))[[1]];
#
y   = x * exp(x + k*sqrt(x + b0));
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
2*x*dy - (2 + 2*x + k*x/sqrt(x + b0)) * y # = 0

# D2 =>
4*x*d2y - (2 + k/sqrt(x + b0)) * (2 + 2*x + k*x/sqrt(x + b0)) * y +
	- k*(1 + b0/(x + b0))/sqrt(x + b0) * y - 4*y # = 0
4*x*d2y - k*(4*x + 3 + b0/(x + b0)) / sqrt(x + b0) * y +
	- (k^2 + 8 + 4*x - k^2*b0/(x + b0)) * y # = 0
4*x*(x + b0)*d2y - k*((4*x + 3)*(x + b0) + b0) / sqrt(x + b0) * y +
	- ((k^2 + 8 + 4*x)*(x + b0) - k^2*b0) * y # = 0

### ODE:
4*x^2*(x + b0)*d2y - 2*x*(4*x^2 + (4*b0+3)*x + 4*b0) * dy +
	+ (4*x^3 + (4*b0 + 6 - k^2)*x^2 + 8*b0*x + 6*x + 8*b0) * y # = 0


### Special Cases:

### k = sqrt(4*b0 + 6)
x = sqrt(3); b0 = 2/3; k = sqrt(4*b0 + 6);
params = list(x=x, k=k, b0=b0);
# Re-run assignments from above;
4*x^2*(x + b0)*d2y - 2*x*(4*x^2 + (4*b0+3)*x + 4*b0) * dy +
	+ (4*x^3 + 8*b0*x + 6*x + 8*b0) * y # = 0

### & b0 = -3/4
x = sqrt(5); b0 = -3/4; k = sqrt(4*b0 + 6); # k = sqrt(3);
params = list(x=x, k=k, b0=b0);
# Re-run assignments from above;
x^2*(4*x - 3)*d2y - 2*x*(4*x^2 - 3) * dy + (4*x^3 - 6) * y # = 0


#################

### y = exp(1/x + k*sqrt(x + b0))

x = sqrt(3); k = -1/5; b0 = 2/3; params = list(x=x, k=k, b0=b0);
e = expression(exp(1/x + k*sqrt(x + b0)))[[1]];
#
y   = eval(e, params)
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
2*x^2*dy + (2 - k*x^2/sqrt(x + b0)) * y # = 0

# D2 =>
4*x^2*d2y + 8*x*dy + (4 - 2*k*x^2/sqrt(x + b0)) * dy +
	- (3*k*x + k*b0 - k*b0^2 / (x + b0)) / sqrt(x + b0) * y # = 0
4*x^4*(x + b0)*d2y + 8*x^3*(x + b0)*dy +
	- x^2 * (3*k*x^2 + 4*k*b0*x - 4*k*x - 4*k*b0) / sqrt(x + b0) * y +
	- (k^2*x^4 + 4*x + 4*b0) * y # = 0

### ODE:
4*x^4*(x + b0)*d2y + 2*x^2*(x^2 + 4*x + 4*b0)*dy +
	- (k^2*x^4 + 6*x^2 + (8*b0-4)*x - 4*b0) * y # = 0

### Special Case: b0 = 1;
x = sqrt(3); k = -3/5; params = list(x=x, k = 2*k, b0 = 1);
# Re-run assignments from above;
2*x^4*(x + 1)*d2y + x^2*(x + 2)^2*dy - (2*k^2*x^4 + 3*x^2 + 2*x - 2) * y # = 0


##############################
##############################

### y = exp(x + sqrt(x^2 + b0))

x = - sqrt(3); b0 = -3/5; params = list(x=x, b0=b0);
e = expression(exp(x + sqrt(x^2 + b0)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
dy - (1 + x/sqrt(x^2 + b0)) * y # = 0

# D2 =>
d2y - (1 + x/sqrt(x^2 + b0)) * dy +
	- (1 - x^2/(x^2 + b0)) / sqrt(x^2 + b0) * y # = 0
d2y - 2*dy + b0/(x^2 + b0) * y +
	- (1 - x^2/(x^2 + b0)) / sqrt(x^2 + b0) * y # = 0
d2y - 2*dy + b0/(x^2 + b0) * y +
	- b0/(x^2 + b0) / sqrt(x^2 + b0) * y # = 0
x*d2y - 2*x*dy + b0*x/(x^2 + b0) * y +
	- b0/(x^2 + b0) * (dy - y) # = 0

### ODE:
x*(x^2 + b0)*d2y - (2*x*(x^2 + b0) + b0)*dy + b0*(x+1)*y # = 0


###################

### y = exp(x^2 + sqrt(x^2 + k))

x = - sqrt(3); k = -2/5; params = list(x=x, k=k);
e = expression(exp(x^2 + sqrt(x^2 + k)))[[1]];
#
y   = exp(x^2 + sqrt(x^2 + k))
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
dy - x*(2 + 1/sqrt(x^2 + k)) * y # = 0

# D2 =>
d2y - x*(2 + 1/sqrt(x^2 + k)) * dy +
	- 2*y - (1 - x^2/(x^2 + k)) / sqrt(x^2 + k) * y # = 0
d2y - 4*x*dy + (4*x^2 - 1)*y + k/(x^2 + k) * y +
	- 2*y - (1 - x^2/(x^2 + k)) / sqrt(x^2 + k) * y # = 0
x*d2y - 4*x^2*dy - k/(x^2 + k) * dy +
	+ x*(4*x^2 - 3)*y + 3*k*x/(x^2 + k) * y # = 0
x*(x^2 + k)*d2y - (4*x^2*(x^2 + k) + k)*dy +
	+ x*((4*x^2 - 3)*(x^2 + k) + 3*k)*y # = 0

### ODE:
x*(x^2 + k)*d2y - (4*x^2*(x^2 + k) + k)*dy +
	+ x^3*(4*x^2 + 4*k - 3)*y # = 0

### Special Cases:

### k = 3/4
x = - sqrt(3); k = 3/4; params = list(x=x, k=k);
# Re-run assignments above;
x*(4*x^2 + 3)*d2y - (16*x^4 + 12*x^2 + 3)*dy + 16*x^5*y # = 0


########
### Gen:

### y = exp(x + k*sqrt(x^2 + b1*x + b0))

x = sqrt(3); k = - sqrt(2); b0 = -3/5; b1 = 2/3;
params = list(x=x, k=k, b0=b0, b1=b1);
e = expression(exp(x + k*sqrt(x^2 + b1*x + b0)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
2*dy - (2 + k*(2*x + b1)/sqrt(x^2 + b1*x + b0)) * y # = 0

# D2 =>
4*d2y - (2 + k*(2*x + b1)/sqrt(x^2 + b1*x + b0))^2 * y +
	- 4*k / sqrt(x^2 + b1*x + b0) * y +
	+ k*(2*x + b1)^2 / sqrt(x^2 + b1*x + b0)^3 * y # = 0
4*(x^2 + b1*x + b0)*d2y +
	- k*(8*x^3 + 12*b1*x^2 + (4*b1^2 + 8*b0)*x +
		- b1^2 + 4*b0*b1 + 4*b0) / sqrt(x^2 + b1*x + b0) * y +
	- (4*(k^2 + 1)*x^2 + 4*b1*(k^2+1)*x + k^2*b1^2 + 4*b0) * y # = 0

### ODE:
4*(2*x + b1)*(x^2 + b1*x + b0)*d2y +
	- 2*(8*x^3 + 12*b1*x^2 + (4*b1^2 + 8*b0)*x +
		- b1^2 + 4*b0*b1 + 4*b0) * dy +
	- (8*(k^2-1)*x^3 + 12*(k^2-1)*b1*x^2 + (6*k^2*b1^2 - 4*b1^2 - 8*b0)*x +
		+ k^2*b1^3 + 2*b1^2 - 4*b0*b1 - 8*b0) * y # = 0


### Special Cases:

### Case: k = 1; b0 = 0;
# b1 => 2*b1;
x = sqrt(3); b1 = 2/5; params = list(x=x, b1 = 2*b1, k = 1, b0 = 0);
# Re-run assignments above;
x*(x + b1)*(x + 2*b1)*d2y +
	- (2*x^3 + 6*b1*x^2 + 4*b1^2*x - b1^2) * dy +
	- b1^2*(x + b1 + 1) * y # = 0

### Case: b1 = -1;
x*(x - 1)*(x - 2)*d2y +
	- (2*x^3 - 6*x^2 + 4*x - 1) * dy - x*y # = 0


############
### Variant:

### y = exp(k*sqrt(x^2 + b1*x + b0))

x = sqrt(3); k = -1/5; b0 = 2/3; b1 = -1/sqrt(5);
params = list(x=x, k=k, b0=b0, b1=b1);
e = expression(exp(k*sqrt(x^2 + b1*x + b0)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


### D =>
2*dy - k*(2*x + b1) / sqrt(x^2 + b1*x + b0) * y # = 0
# =>
2*sqrt(x^2 + b1*x + b0)*dy - k*(2*x + b1)*y # = 0

### D2 =>
4*(x^2 + b1*x + b0)*d2y - 2*k*(2*x + b1) * sqrt(x^2 + b1*x + b0) * dy +
	- k * (4*b0 - b1^2) / sqrt(x^2 + b1*x + b0) * y # = 0
4*(x^2 + b1*x + b0)*d2y - k^2*(2*x + b1)^2*y +
	- k * (4*b0 - b1^2) / sqrt(x^2 + b1*x + b0) * y # = 0


### ODE:
4*(2*x + b1)*(x^2 + b1*x + b0)*d2y - 2*(4*b0 - b1^2)*dy +
	- k^2*(2*x + b1)^3*y # = 0

# b1 => 2*b1;
params = list(x=x, k=k, b0=b0, b1 = 2*b1);
# Re-run assignments above;
(x + b1)*(x^2 + 2*b1*x + b0)*d2y - (b0 - b1^2)*dy +
	- k^2*(x + b1)^3*y # = 0

# Note: Case b0 = b1^2 is trivial;


###################

### y = exp(k*sqrt(x^2 + b1*x + b0) + x)

x = sqrt(3); k = -1/5; b0 = 2/3; b1 = -1/sqrt(5);
params = list(x=x, k=k, b0=b0, b1=b1);
e = expression(exp(k*sqrt(x^2 + b1*x + b0) + x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


### D =>
2*dy - k*(2*x + b1) / sqrt(x^2 + b1*x + b0) * y - 2*y # = 0
# =>
2*sqrt(x^2 + b1*x + b0)*(dy - y) - k*(2*x + b1)*y # = 0

### D2 =>
4*(x^2 + b1*x + b0)*(d2y - dy) - k^2*(2*x + b1)^2*y  +
	- 2*k*(2*x + b1) * sqrt(x^2 + b1*x + b0) * y +
	- k*(4*b0 - b1^2) / sqrt(x^2 + b1*x + b0) * y # = 0
4*(x^2 + b1*x + b0)*(d2y - dy) - k^2*(2*x + b1)^2 * y  +
	- 2*((4*b0 - b1^2) + 2*(2*x + b1)*(x^2 + b1*x + b0)) *
		(dy - y) / (2*x + b1) # = 0

### ODE:
4*(2*x + b1)*(x^2 + b1*x + b0)*d2y +
	- 2*(8*x^3 + 12*b1*x^2 + 4*(b1^2 + 2*b0)*x +
		- b1^2 + 4*b0*b1 + 4*b0) * dy +
	- (8*(k^2 - 1)*x^3 + 12*b1*(k^2 - 1)*x^2 +
		+ 6*b1^2*(k^2 - 1)*x + 2*(b1^2 - 4*b0)*x +
		+ b1^3*k^2 + 2*b1^2 - 4*b0*(b1 + 2)) * y # = 0


### Special Cases:

### k = +/- 1; b1 = -2;
params = list(x=x, b0=b0, b1 = -2, k = 1);
# Re-run assignments above;
(x - 1)*(x^2 - 2*x + b0)*d2y +
	- (2*x^3 - 6*x^2 + 2*(b0 + 2)*x - (b0 + 1)) * dy +
	+ (b0 - 1)*x*y # = 0

# Shift: z = x - 1;
z = x - 1;
z*(z^2 + b0 - 1)*d2y +
	- (2*z^3 + 2*(b0 - 1)*z + b0 - 1) * dy +
	+ (b0 - 1)*(z + 1)*y # = 0
# b0 = 7/8 & Shift:
params = list(x = x+1, b0 = 7/8, b1 = -2, k = 1);
# Re-run assignments above;
x*(8*x^2 - 1)*d2y - (2*x - 1)*(8*x^2 + 4*x + 1)*dy - (x+1)*y # = 0


### b0 = 0
params = list(x=x, k=k, b1=b1, b0 = 0);
# Re-run assignments above;
4*x*(2*x + b1)*(x + b1)*d2y +
	- 2*(8*x^3 + 12*b1*x^2 + 4*b1^2*x - b1^2) * dy +
	- ((k^2 - 1)*(8*x^3 + 12*b1*x^2 + 6*b1^2*x) +
		+ 2*b1^2 * x + b1^3*k^2 + 2*b1^2) * y # = 0


###################

###################
### Double SQRT ###

### y = sqrt(x + b0) * exp(k*sqrt(x + b0))

x = sqrt(3); k = -1/5; b0 = 2/3; params = list(x=x, k=k, b0=b0);
e = expression(sqrt(x + b0) * exp(k*sqrt(x + b0)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


### D =>
2*dy - y / (x+b0) - k * exp(k*sqrt(x + b0)) # = 0

### D2 =>
4*(x+b0)^2*d2y - 2*(x+b0)*dy - (k^2*(x+b0) - 2) * y # = 0


### ODE:
4*(x+b0)^2*d2y - 2*(x+b0)*dy - (k^2*(x+b0) - 2) * y # = 0


######################

### y = 1/sqrt(x + b0) * exp(k*sqrt(x + b0))

x = sqrt(3); k = -1/5; b0 = 2/3; params = list(x=x, k=k, b0=b0);
e = expression(1/sqrt(x + b0) * exp(k*sqrt(x + b0)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


### D =>
2*(x+b0)*dy + y - k * exp(k*sqrt(x + b0)) # = 0

### D2 =>
4*(x+b0)*d2y + 6*dy - k^2 * y # = 0

### ODE:
4*(x+b0)*d2y + 6*dy - k^2*y # = 0


####################

### y = sqrt(x + b0) * exp(k / sqrt(x + b0))

x = sqrt(3); k = -1/5; b0 = 2/3; params = list(x=x, k=k, b0=b0);
e = expression(sqrt(x + b0) * exp(k / sqrt(x + b0)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


### D =>
2*(x+b0)*dy - y + k * exp(k / sqrt(x + b0)) # = 0

### D2 =>
4*(x+b0)*d2y + 2*dy - k^2 / (x+b0)^2 * y # = 0


### ODE:
4*(x+b0)^3 * d2y + 2*(x+b0)^2 * dy - k^2 * y # = 0

