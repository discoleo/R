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


######################

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
((- 6*k*x^4 + 6*k^2*x^4 - 36*k^2*p*x^4 - 6*k*x^3 + 6*k^2*x^3 - 72*k^2*p*x^3 - 36*k^2*p*x^2) * y*dy +
 - (3*k^3*x^4 + x^3 + k*x^3 - 6*p*k*x^3 - 12*p*k^2*x^3 - 18*p^2*k^2*x^3 + 3*k^3*x^3 - 6*p*k*x^2 +
 - 30*p*k^2*x^2 - 36*p^2*k^2*x^2 - 18*p*k^2*x - 18*p^2*k^2*x) * y^2) * 9*x*(x+1)^2 * y*d2y +
 
(9*y^2*x^4 + 18*y^2*x^5 + 9*y^2*x^6 + 108*y^2*p*x^3*k - 18*y^2*x^4*k + 324*y^2*p*x^4*k - 36*y^2*x^5*k +
 + 324*y^2*p*x^5*k - 18*y^2*x^6*k + 108*y^2*p*x^6*k + 324*y^2*p^2*x^2*k^2 - 108*y^2*p*x^3*k^2 +
 + 1296*y^2*p^2*x^3*k^2 + 9*y^2*x^4*k^2 - 324*y^2*p*x^4*k^2 + 1944*y^2*p^2*x^4*k^2 + 18*y^2*x^5*k^2 +
 - 324*y^2*p*x^5*k^2 + 1296*y^2*p^2*x^5*k^2 + 9*y^2*x^6*k^2 - 108*y^2*p*x^6*k^2 +
 + 324*y^2*p^2*x^6*k^2) * dy^2 +
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

