#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Exp of Exp
##
## draft v.0.1b


### Exp of Atan to Power:
# y = Exp( Exp(P(x)) )
# y = Exp( Exp(Atan(P(x))) )


### Examples:

# x*y*d2y - x * dy^2 - (k2*n*x^n + n-1) * y*dy = 0;
# (x^2+k^2) * y*d2y - (x^2+k^2) * dy^2 + (2*x - k2*k) * y*dy = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y = Exp(k1 * Exp(k2*x^n) )
# Note: k1 plays NO role in the ODE;

# Check:
k1 = exp(-2/3); k2 = exp(-1/3);
n = 1/sqrt(5);
x = sqrt(3); params = list(x=x, n=n, k1=k1, k2=k2);
e = expression(exp(k1*exp(k2*x^n)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
dy - k1*k2*n*x^(n-1) * exp(k2*x^n) * y # = 0

# D2 =>
d2y - k1*k2^2*n^2*x^(2*n-2) * exp(k2*x^n) * y +
	- k1*k2*n*(n-1)*x^(n-2) * exp(k2*x^n) * y +
	- k1*k2*n*x^(n-1) * exp(k2*x^n) * dy # = 0

### ODE:
x*y*d2y - x * dy^2 - (k2*n*x^n + n-1) * y*dy # = 0


#########################

### y = Exp(k1 * Exp(k2 * Atan(x/k)) )
# Note: k1 plays NO role in the ODE;

# Check:
k = 1/sqrt(5);
k1 = exp(-2/3); k2 = exp(-1/3);
x = sqrt(3); params = list(x=x, k=k, k1=k1, k2=k2);
e = expression(exp(k1*exp(k2*atan(x/k))))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2+k^2)*dy - k1*k2*k * exp(k2*atan(x/k)) * exp(k1*exp(k2*atan(x/k))) # = 0
(x^2+k^2)*dy - k1*k2*k * exp(k2*atan(x/k)) * y # = 0

# D2 =>
(x^2+k^2)^2 * d2y + 2*x*(x^2+k^2) * dy +
	- k1*k2^2*k^2 * exp(k2*atan(x/k)) * y +
	- k1*k2*k*(x^2+k^2) * exp(k2*atan(x/k)) * dy # = 0

### ODE:
(x^2+k^2) * y*d2y - (x^2+k^2) * dy^2 + (2*x - k2*k) * y*dy # = 0

