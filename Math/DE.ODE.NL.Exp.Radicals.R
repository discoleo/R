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

# library(deSolve)

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
### y = x^p * exp((x+1)^r + k*x)

# Check:
p = 1/3; # p = 0;
k = sqrt(2); n = 5; # n = -1/5;
x = sqrt(3);
params = list(x=x, k=k, p=p, r = 1/n);
e = expression(x^p * exp((x+1)^r + k*x))[[1]];
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


######################

### Type: Coupled

### y = (x+b)^(-r) * exp(k*(x+b)^r)

# Check:
r = 1/5;
k = sqrt(2); b = sqrt(5);
x = sqrt(3);
params = list(x=x, k=k, b=b, r=r);
e = expression((x+b)^(-r) * exp(k*(x+b)^r))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x+b)*dy + r*y - k*r*exp(k*(x+b)^r) # = 0

# D2 =>
(x+b)^2 * d2y + (r+1)*(x+b) * dy - k^2*r^2 * (x+b)^r*exp(k*(x+b)^r) # = 0
(x+b)^2 * d2y + (r+1)*(x+b) * dy -((x+b)*dy + r*y)^2 / y # = 0

### ODE:
(x+b)^2 * y*d2y - (x+b)^2 * dy^2 - (r-1)*(x+b) * y*dy - r^2 * y^2 # = 0


######################
######################

### Type: Exp of 2 Radicals
# ODE-Type: y^2 * d2y^2

### y = x^p * exp((x+1)^(1/3) + k*(x+1)^(2/3))

# Check:
p = 1/3; # p = 0;
k = sqrt(2);
x = sqrt(3);
params = list(x=x, k=k, p=p);
e = expression(x^p * exp((x+1)^(1/3) + k/2*(x+1)^(2/3)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);
r1  = (x+1)^(1/3); r2 = r1^2;

# D =>
3*x*(x+1)*dy - 3*p*(x+1)*y - x*(r1 + k*r2) * y # = 0

# D2 =>
9*x*(x+1) * d2y + 9*(2*x+1) * dy - 9*p*(x+1) * dy - 9*p*y +
	- 3*x*(r1 + k*r2) * dy +
	- 3*(r1 + k*r2) * y - x*(r1 + 2*k*r2)/(x+1) * y # = 0
9*x*(x+1)^2 * d2y + 9*(x+1)*(2*x+1) * dy - 9*p*(x+1)^2 * dy +
	- (r1 + k*r2) * (3*p*(x+1) + x*(r1 + k*r2)) * y +
	- 3*(x+1)*(r1 + k*r2) * y - x*(r1 + 2*k*r2) * y - 9*p*(x+1) * y # = 0
9*x*(x+1)^2 * d2y + 9*(x+1)*(2*x+1) * dy - 9*p*(x+1)^2 * dy +
	- (5*k*x*r2 + 3*k*p*x*r2 + x*r2 + 3*k*p*r2 + 3*k*r2 +
		+ k^2*x*(x+1)*r1 + 4*x*r1 + 3*p*x*r1 + 3*p*r1 + 3*r1) * y +
	- (2*k*x*(x+1) + 9*p*x + 9*p) * y # = 0

# TODO

81*k^2*x^4*(x+1)^4 * y^2*d2y^2 +
 - 9*x*(x+1)^2 * (6*k*((6*k*p - k + 1)*x^2 + (12*k*p - k + 1)*x + 6*k*p)*x^2 * y*dy +
	+ (3*k^3*x^4 + (3*k^3 - 18*p^2*k^2 - 12*p*k^2 - 6*p*k + k + 1)*x^3 +
		- 6*p*k*(6*p*k + 5*k + 1)*x^2 - 18*p*(p+1)*k^2*x) * y^2) * y*d2y +
 + (9*x^4 - 18*k*x^4 + 9*k^2*x^4 + 108*k*p*x^4 - 108*k^2*p*x^4 + 324*k^2*p^2*x^4 +
	+ 18*x^3 - 36*k*x^3 + 18*k^2*x^3 + 324*k*p*x^3 - 324*k^2*p*x^3 + 1296*k^2*p^2*x^3 +
	+ 9*x^2 - 18*k*x^2 + 9*k^2*x^2 + 324*k*p*x^2 - 324*k^2*p*x^2 + 1944*k^2*p^2*x^2 +
	+ 108*k*p*x - 108*k^2*p*x + 1296*k^2*p^2*x + 324*k^2*p^2)*x^2 * y^2*dy^2 +

(- 6*y^3*x^4 - 6*y^3*x^5 - 54*y^3*p*x^2*k - 162*y^3*p^2*x^2*k - 108*y^3*p*x^3*k - 486*y^3*p^2*x^3*k +
 - 6*y^3*x^4*k - 54*y^3*p*x^4*k - 486*y^3*p^2*x^4*k - 6*y^3*x^5*k - 162*y^3*p^2*x^5*k +
 - 324*y^3*p^2*x*k^2 - 324*y^3*p^3*x*k^2 + 54*y^3*p*x^2*k^2 - 1134*y^3*p^2*x^2*k^2 - 1296*y^3*p^3*x^2*k^2 +
 + 144*y^3*p*x^3*k^2 - 1458*y^3*p^2*x^3*k^2 - 1944*y^3*p^3*x^3*k^2 + 15*y^3*x^4*k^2 + 126*y^3*p*x^4*k^2 +
 - 810*y^3*p^2*x^4*k^2 - 1296*y^3*p^3*x^4*k^2 + 30*y^3*x^5*k^2 + 36*y^3*p*x^5*k^2 - 162*y^3*p^2*x^5*k^2 +
 - 324*y^3*p^3*x^5*k^2 + 15*y^3*x^6*k^2 + 54*y^3*p*x^3*k^3 - 3*y^3*x^4*k^3 + 162*y^3*p*x^4*k^3 +
 - 6*y^3*x^5*k^3 + 162*y^3*p*x^5*k^3 - 3*y^3*x^6*k^3 + 54*y^3*p*x^6*k^3 - 3*y^3*x^4*k^5 - 9*y^3*x^5*k^5 +
 - 9*y^3*x^6*k^5 - 3*y^3*x^7*k^5) * dy +
- 9*y^4*p*x^2 - 12*y^4*p*x^3 - 3*y^4*p*x^4 + 54*y^4*p^2*x*k + 54*y^4*p^3*x*k - 9*y^4*p*x^2*k +
 + 135*y^4*p^2*x^2*k + 162*y^4*p^3*x^2*k - 12*y^4*p*x^3*k + 108*y^4*p^2*x^3*k + 162*y^4*p^3*x^3*k +
 + 2*y^4*x^4*k - 3*y^4*p*x^4*k + 27*y^4*p^2*x^4*k + 54*y^4*p^3*x^4*k + 2*y^4*x^5*k + 81*y^4*p^2*k^2 +
 + 162*y^4*p^3*k^2 + 81*y^4*p^4*k^2 + 270*y^4*p^2*x*k^2 + 594*y^4*p^3*x*k^2 + 324*y^4*p^4*x*k^2 +
 + 333*y^4*p^2*x^2*k^2 + 810*y^4*p^3*x^2*k^2 + 486*y^4*p^4*x^2*k^2 - 15*y^4*p*x^3*k^2 + 180*y^4*p^2*x^3*k^2 +
 + 486*y^4*p^3*x^3*k^2 + 324*y^4*p^4*x^3*k^2 + 2*y^4*x^4*k^2 - 30*y^4*p*x^4*k^2 + 36*y^4*p^2*x^4*k^2 +
 + 108*y^4*p^3*x^4*k^2 + 81*y^4*p^4*x^4*k^2 + 2*y^4*x^5*k^2 - 15*y^4*p*x^5*k^2 - 27*y^4*p*x^2*k^3 +
 - 27*y^4*p^2*x^2*k^3 - 78*y^4*p*x^3*k^3 - 81*y^4*p^2*x^3*k^3 - 75*y^4*p*x^4*k^3 - 81*y^4*p^2*x^4*k^3 +
 - 24*y^4*p*x^5*k^3 - 27*y^4*p^2*x^5*k^3 + 2*y^4*x^4*k^4 + 4*y^4*x^5*k^4 + 2*y^4*x^6*k^4 + 3*y^4*p*x^3*k^5 +
 + 9*y^4*p*x^4*k^5 + 9*y^4*p*x^5*k^5 + 3*y^4*p*x^6*k^5 # = 0

