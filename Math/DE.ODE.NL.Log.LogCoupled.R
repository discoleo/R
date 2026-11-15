########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Coupled Log of Log
##
## draft v.0.1a


### Coupled Logs:
# y = log(P(x)) * log(log(P(x))) + F0(x)


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### Examples:

### y = log(x^2 + k) * log(log(x^2 + k)) + f

# Check:
k = 1/sqrt(5);
x = sqrt(3); params = list(x=x, k=k);
e = expression(log(x^2 + k) * log(log(x^2 + k)) + exp(2/x))[[1]];
f = exp(2/x); df = -2/x^2 * f; d2f = 4*(x+1)/x^4 * f;
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### D =>
(x^2+k)*log(x^2 + k)*dy - 2*x*(y - f) - 2*x*log(x^2 + k) - (x^2+k)*df*log(x^2 + k) # = 0
# log(x^2 + k) = 2*x*(y - f) / ((x^2+k)*(dy - df) - 2*x);

### D2 =>
(x^2+k)*log(x^2 + k)*d2y + 2*x*dy + 2*x * log(x^2 + k)*dy +
	- 2*x*(dy - df) - 2*(y - f) - 4*x^2/(x^2 + k) - 2*log(x^2 + k) +
	- 2*x*df - 2*x*df*log(x^2 + k) - (x^2+k)*d2f*log(x^2 + k) # = 0

### ODE:
x*(x^2+k)^2 * y*d2y - x*f*(x^2+k)^2 * d2y +
	+ (x^4 - k^2) * y*dy - (x^2+k)*(2*x^2 + f*(x^2 - k)) * dy +
	- (x^2+k)*(d2f*x^3 + df*x^2 + k*d2f*x - k*df) * y +
 + d2f*f*x^5 + 2*df*x^4 + df*f*x^4 + 4*x^3 + 2*k*d2f*f*x^3 + 2*df*k*x^2 + k^2*d2f*f*x - df*k^2*f # = 0


### Special Case:
# Case: f = ct;
f = 2/3; df = d2f = 0; # f = -2; df = d2f = 0;
e = expression(log(x^2 + k) * log(log(x^2 + k)) + f)[[1]];
### ODE:
x*(x^2+k)^2 * y*d2y - x*f*(x^2+k)^2 * d2y +
	+ (x^4 - k^2) * y*dy - (x^2+k)*(2*x^2 + f*(x^2 - k)) * dy + 4*x^3 # = 0

# Case: f = -2;
x*(x^2+k)^2 * y*d2y + 2*x*(x^2+k)^2 * d2y +
	+ (x^4 - k^2) * y*dy - 2*k*(x^2+k) * dy + 4*x^3 # = 0


### Derivation:
pxk = toPoly.pm("x^2 + k")
pD = toPoly.pm("pxk*L*d2y + 2*x*dy + 2*x * L*dy +
	- 2*x*dy + 2*x*df - 2*y + 2*f - 4*x^2*pxk.inv - 2*L +
	- 2*x*df - 2*x*df*L - pxk*d2f*L")
pL = toPoly.pm("2*x*y - 2*x*f")
pLDiv = toPoly.pm("pxk*dy - pxk*df - 2*x")
#
pR = replace.fr.pm(pD, toPoly.pm(1), pxk, x="pxk.inv", pow=1)
pR = replace.fr.pm(pR, pL, pLDiv, x="L", pow=1)
pR = replace.pm(pR, pxk, x="pxk", pow=1)
pR = simplify.spm(pR, do.gcd=TRUE)
sort.dpm(pR)


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

