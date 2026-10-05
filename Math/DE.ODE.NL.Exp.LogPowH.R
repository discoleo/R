#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Exp of Log^Power
##
## draft v.0.1a


### Exp of Log to Power:
# y = Exp( Log(P(x))^r )


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### Power = 1/2

### y = Log(x)^(5/2) * Exp(k * Log(x)^(1/2) )

# Check:
k = 1/sqrt(5);
x = sqrt(3); params = list(x=x, k=k);
e = expression(log(x)^(5/2) * exp(k*log(x)^0.5))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*x*dy - 5*log(x)^(3/2) * exp(k*log(x)^0.5) +
	- k*log(x)^2 * exp(k*log(x)^0.5) # = 0

# D2 =>
4*x^2*d2y + 4*x*dy +
	- 15*log(x)^(1/2) * exp(k*log(x)^0.5) +
	- 9*k*log(x) * exp(k*log(x)^0.5) +
	- k^2*log(x)^(3/2) * exp(k*log(x)^0.5) # = 0

### ODE:
500*x^3 * y^2*d2y^2 +
	- (600*x^3 * y*dy^2 - 50*(k^2 + 20)*x^2 * y^2*dy - k^4*x * y^3) * d2y +
	+ 180*x^3 * dy^4 - 24*(2*k^2 + 25)*x^2 * y*dy^3 + 50*(k^2 + 10)*x * y^2*dy^2 +
	- k^4*x * y^2*dy^2 + k^4 * y^3*dy # = 0


# Derivation:
EXP = exp(k*log(x)^0.5);
L3 = log(x)^(3/2) * EXP;
L5 = L3 * log(x);
L2 = sqrt(L3*L5);
#
y - L5 # = 0
L2^2 - y*L3 # = 0
5*L3 + k*L2 - 2*x*dy # = 0
k^2*L3 + 9*k*L3^2/L2 + 15*L3^2/y - 4*x^2*d2y - 4*x*dy # = 0
# =>
5*L2^2 + k*y*L2 - 2*x*y*dy # = 0
k^2*L2^2 - 6*k*x*dy*L2 + 20*x^2*y*d2y - 12*x^2*dy^2 + 20*x*y*dy # = 0


#########################
#########################

### Arbitrary Powers:

### y = Exp(k * Log(x)^n )

# Check:
n = exp(1);
k = 1/sqrt(5);
x = sqrt(3); params = list(x=x, k=k, n=n);
e = expression(exp(k*log(x)^n))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - n*k*log(x)^(n-1) * exp(k*log(x)^n) # = 0
x*dy - n*k*log(x)^(n-1) * y # = 0

# D2 =>
x^2 * d2y + x*dy +
	- n*(n-1)*k*log(x)^(n-2) * y +
	- n*k*x * log(x)^(n-1) * dy # = 0
x^2 * y*d2y - x^2 * dy^2 + x * y*dy +
	- n*(n-1)*k*log(x)^(n-2) * y^2 # = 0

# TODO: uglier NL ODE;


######################

### y = Log(x)^(n-1) * Exp(k * Log(x)^n )

# Check:
n = exp(1);
k = 1/sqrt(5);
x = sqrt(3); params = list(x=x, k=k, n=n);
e = expression(log(x)^(n-1) * exp(k*log(x)^n))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - (n-1) * log(x)^(n-2) * exp(k*log(x)^n) +
	- k*n * log(x)^(2*n-2) * exp(k*log(x)^n) # = 0

# D2 =>
x^2*d2y + x*dy - (n-1)*(n-2) * log(x)^(n-3) * exp(k*log(x)^n) +
	- 3*k*n*(n-1) * log(x)^(2*n-3) * exp(k*log(x)^n) +
	- k^2*n^2 * log(x)^(3*n-3) * exp(k*log(x)^n) # = 0

# System:
Lb = log(x); EXP = exp(k*log(x)^n);
L1 = log(x)^(n-1); L2 = log(x)^(n-2);
#
L1 * EXP - y # = 0
L2 * Lb - L1 # = 0
k*n * y * Lb*L1 + (n-1) * y - x*dy * Lb # = 0
k^2*n^2 * y * L1^3 + 3*k*n*(n-1) * y * L1*L2 +
	+ (n-1)*(n-2) * L2^2 * EXP - (x^2*d2y + x*dy) * L1 # = 0

# =>
k*n * y * L1^2 + (n-1) * y * L2 - x*dy * L1 # = 0
k^2*n^2 * y * L1^4 + 3*k*n*(n-1) * y * L1^2*L2 +
	+ (n-1)*(n-2) * y * L2^2 - (x^2*d2y + x*dy) * L1^2 # = 0
# =>
k^2*(n-1)*n^2 * y^2 * L1^2 - 3*k*n*(n-1) * y * L1 * (k*n * y * L1 - x*dy) +
	+ (n-2) * (k*n * y * L1 - x*dy)^2 - (n-1)*(x^2*d2y + x*dy) * y # = 0
# =>
(n-1)*(x^2*d2y + x*dy) * y - (n-2)*x^2 * dy^2 - k*n*(n+1)*x * y*dy * L1 +
	+ k^2*n^3 * y^2 * L1^2 # = 0

# TODO: NO easy way to substitute L1;

