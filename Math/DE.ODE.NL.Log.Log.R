########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log of Log
##
## draft v.0.1d


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


#######################

### Variant: + extra Log-term;
### y = x * log(x + k*log(x)) + c1*x*log(x)

# Check:
c1 = 1/sqrt(5); # c1 = -1; # c1 = -1/2;
k = sqrt(2);
x = sqrt(3);
params = list(x=x, k=k);
e = expression(x * log(x + k*log(x)) + c1*x*log(x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - y - x * (x+k) / (x + k*log(x)) - c1*x # = 0

# D2 =>
x*d2y + (x+k)^2 / (x + k*log(x))^2 +
	- (2*x+k) / (x + k*log(x)) - c1 # = 0
x^3*(x+k) * d2y + (x+k)*(x*dy - y - c1*x)^2 - x*(2*x+k) * (x*dy - y - c1*x) +
	- c1*x^2*(x+k) # = 0

### ODE:
x^3*(x+k) * d2y + (x+k) * (x*dy - y)^2 +
	- x*(2*(c1+1)*x + k*(2*c1+1)) * (x*dy - y) +
	+ c1^2*x^2*(x+k) + c1*x^3 # = 0

### Special Cases:

### c1 = -1;
x^3*(x+k) * d2y + (x+k) * (x*dy - y)^2 + k*x * (x*dy - y) + k*x^2 # = 0

### c1 = -1/2;
4*x^3*(x+k) * d2y + 4*(x+k) * (x*dy - y)^2 - 4*x^2 * (x*dy - y) - x^2*(x-k) # = 0


###########################
###########################

### y = log(P(x)) * log(log(P(x))) + F0(x)

### Example:
### y = log(x^2 + k) * log(log(x^2 + k)) + f

### D =>
(x^2+k)*log(x^2 + k)*dy - 2*x*(y - f) - 2*x*log(x^2 + k) - (x^2+k)*df*log(x^2 + k) # = 0
# log(x^2 + k) = 2*x*(y - f) / ((x^2+k)*(dy - df) - 2*x);

### D2 =>
(x^2+k)*log(x^2 + k)*d2y + 2*x*dy + 2*x/(x^2+k) * log(x^2 + k)*dy +
	- 2*x*(dy - df) - 2*(y - f) - 4*x^2/(x^2 + k) - 2*log(x^2 + k) +
	- 2*x*df - 2*x*df*log(x^2 + k) - (x^2+k)*d2f*log(x^2 + k) # = 0

### ODE:
(x^5 + 2*k*x^3 + k^2*x)*y*d2y - (x^5 + 2*k*x^3 + k^2*x)*f*d2y +
	- (x^4 - 2*x^2 + 2*k*x^2 + k^2)*y*dy +
	- (2*x^4 - f*x^4 + 2*k*x^2 + 2*f*x^2 - 2*k*f*x^2 - k^2*f)*dy +
	- (d2f*x^5 + df*x^4 + 2*k*d2f*x^3 + k^2*d2f*x - df*k^2)*y +
	+ d2f*f*x^5 + 2*df*x^4 + df*f*x^4 + 4*x^3 + 2*k*d2f*f*x^3 + 2*df*k*x^2 + k^2*d2f*f*x - df*k^2*f

### Special Case:
# f = ct;
(x^5 + 2*k*x^3 + k^2*x)*y*d2y - f*(x^5 + 2*k*x^3 + k^2*x)*d2y +
	- (x^4 - 2*x^2 + 2*k*x^2 + k^2)*y*dy +
	- (2*x^4 - f*x^4 + 2*k*x^2 + 2*f*x^2 - 2*k*f*x^2 - k^2*f)*dy + 4*x^3 # = 0
# f = 2; k = -2;
x*(x^2-2)^2*y*d2y - 2*x*(x^2-2)^2*d2y - (x^2-2)^2*y*dy + 8*(x^2 - 1)*dy + 4*x^3 # = 0

# TODO: check;

### Derivation:
pxk = toPoly.pm("x^2 + k")
pD = toPoly.pm("pxk*L*d2y + 2*x*dy + 2*x*pxk.inv * L*dy +
	- 2*x*dy + 2*x*df - 2*y + 2*f - 4*x^2*pxk.inv - 2*L +
	- 2*x*df - 2*x*df*L - pxk*d2f*L")
pL = toPoly.pm("2*x*y - 2*x*f")
pLDiv = toPoly.pm("pxk*dy - pxk*df - 2*x")
#
pR = replace.fr.pm(pD, toPoly.pm(1), pxk, x="pxk.inv", pow=1)
pR = replace.fr.pm(pR, pL, pLDiv, x="L", pow=1)
pR = replace.pm(pR, pxk, x="pxk", pow=1)
pR = simplify.spm(pR, do.gcd=TRUE)
print.dpm(pR)


#########################

### y = log(P1(x)) * log(log(P1(x))) + log(P2(x)) * log(log(P2(x))) + F0(x)
### Persistence of: log(P1) & log(P2)

### D =>
dy - LL1*dp1/p1 - LL2*dp2/p2 - dp1/p1 - dp2/p2 - df0 # = 0
p1*p2*dy - p2*dp1*LL1 - p1*dp2*LL2 - p2*dp1 - p1*dp2 - df0 # = 0

### D2 =>

### Linear System:
# TODO
# Note: one level of log() will persist;

