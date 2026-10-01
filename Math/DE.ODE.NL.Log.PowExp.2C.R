#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs: Log * Exp
## w. 2 Coupled Components
##
## draft v.0.1c

### Log to Power: 2 Coupled Components

# Base: y0 = B1(x) * Log(P1(x))^2 * Exp(P2(x)});
# y = y0 + B2(x) * Exp(P2(x)});
# y = y0 + B2(x) * Log(P1(x)) * Exp(P2(x)});

# - For Base-Case, see file:
#   DE.ODE.NL.Log.PowExp.R;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### Type: 2 Components w. Exp
# y = Log(x)^2 * Exp(k*x) + c1*Exp(k*x)

# Check:
k  = 1/sqrt(5);
c1 = 2^(1/3);
x = sqrt(3); params = list(x=x, k=k, c1=c1);
e = expression(log(x)^2 * exp(k*x) + c1*exp(k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - k*x*y - 2*log(x)*exp(k*x) # = 0

# D2 =>
x^2*d2y - x*(2*k*x - 1)*dy + k*x*(k*x - 1)*y - 2*exp(k*x) # = 0

# =>
2*log(x)^2 * exp(k*x) # ==
- c1*x^2*d2y + c1*x*(2*k*x - 1)*dy - c1*k*x*(k*x - 1)*y + 2*y;
# =>
x^2*d2y - x*(2*k*x - 1)*dy + k*x*(k*x - 1)*y +
	+ (x*dy - k*x*y)^2 / (c1*x^2*d2y - c1*x*(2*k*x - 1)*dy + c1*k*x*(k*x - 1)*y - 2*y) # = 0

# TODO: simplify;


#########################

### y = Log(x)^3 * Exp(k*x) + c1*Exp(k*x)

# Check:
k  = 1/sqrt(5);
c1 = 2^(1/3);
x = sqrt(3); params = list(x=x, k=k, c1=c1);
e = expression(log(x)^3 * exp(k*x) + c1*exp(k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - k*x*y - 3*log(x)^2*exp(k*x) # = 0

# D2 =>
x^2*d2y - x*(k*x-1)*dy - k*x*y +
	- 3*k*x*log(x)^2*exp(k*x) - 6*log(x)*exp(k*x) # = 0
x^2*d2y - x*(2*k*x-1)*dy + k*x*(k*x-1)*y +
	- 6*log(x)*exp(k*x) # = 0

# System:
Lb = log(x); EXP = exp(k*x);
L1 = Lb * EXP; L2 = L1 * Lb; L3 = L2 * Lb;
#
L3 + c1*EXP - y # = 0
3*L2 - x*dy + k*x*y # = 0
6*L1 - x^2*d2y + x*(2*k*x-1)*dy - k*x*(k*x-1)*y # = 0
L2^2 - L3*L1 # = 0
L1^2 - L2*EXP # = 0

### ODE:
c1*x^4*(dy - k*y)^4 * d2y^3 +
 - 6*k*c1*x^4*dy^5*d2y^2 + 3*c1*x^3*dy^5*d2y^2 + 27*k^2*c1*x^4*y*dy^4*d2y^2 +
 - 15*k*c1*x^3*y*dy^4*d2y^2 - 48*k^3*c1*x^4*y^2*dy^3*d2y^2 + 30*k^2*c1*x^3*y^2*dy^3*d2y^2 +
 + 42*k^4*c1*x^4*y^3*dy^2*d2y^2 - 30*k^3*c1*x^3*y^3*dy^2*d2y^2 - 18*k^5*c1*x^4*y^4*dy*d2y^2 +
 + 15*k^4*c1*x^3*y^4*dy*d2y^2 + 3*k^6*c1*x^4*y^5*d2y^2 - 3*k^5*c1*x^3*y^5*d2y^2 +
 + 12*k^2*c1*x^4*dy^6*d2y - 12*k*c1*x^3*dy^6*d2y + 3*c1*x^2*dy^6*d2y +
 - 60*k^3*c1*x^4*y*dy^5*d2y + 66*k^2*c1*x^3*y*dy^5*d2y - 18*k*c1*x^2*y*dy^5*d2y - 12*x*y*dy^5*d2y +
 + 123*k^4*c1*x^4*y^2*dy^4*d2y - 150*k^3*c1*x^3*y^2*dy^4*d2y +
 + 45*k^2*c1*x^2*y^2*dy^4*d2y + 60*k*x*y^2*dy^4*d2y - 132*k^5*c1*x^4*y^3*dy^3*d2y + 180*k^4*c1*x^3*y^3*dy^3*d2y +
 - 60*k^3*c1*x^2*y^3*dy^3*d2y - 120*k^2*x*y^3*dy^3*d2y + 78*k^6*c1*x^4*y^4*dy^2*d2y +
 - 120*k^5*c1*x^3*y^4*dy^2*d2y + 45*k^4*c1*x^2*y^4*dy^2*d2y + 120*k^3*x*y^4*dy^2*d2y - 24*k^7*c1*x^4*y^5*dy*d2y +
 + 42*k^6*c1*x^3*y^5*dy*d2y - 18*k^5*c1*x^2*y^5*dy*d2y - 60*k^4*x*y^5*dy*d2y + 3*k^8*c1*x^4*y^6*d2y +
 - 6*k^7*c1*x^3*y^6*d2y + 3*k^6*c1*x^2*y^6*d2y + 12*k^5*x*y^6*d2y - 8*k^3*c1*x^4*dy^7 + 12*k^2*c1*x^3*dy^7 +
 - 6*k*c1*x^2*dy^7 + 8*x*dy^7 + c1*x*dy^7 + 44*k^4*c1*x^4*y*dy^6 - 72*k^3*c1*x^3*y*dy^6 +
 + 39*k^2*c1*x^2*y*dy^6 - 32*k*x*y*dy^6 - 7*k*c1*x*y*dy^6 - 12*y*dy^6 - 102*k^5*c1*x^4*y^2*dy^5 +
 + 183*k^4*c1*x^3*y^2*dy^5 - 108*k^3*c1*x^2*y^2*dy^5 + 36*k^2*x*y^2*dy^5 + 21*k^2*c1*x*y^2*dy^5 + 72*k*y^2*dy^5 +
 + 129*k^6*c1*x^4*y^3*dy^4 - 255*k^5*c1*x^3*y^3*dy^4 + 165*k^4*c1*x^2*y^3*dy^4 + 20*k^3*x*y^3*dy^4 +
 - 35*k^3*c1*x*y^3*dy^4 - 180*k^2*y^3*dy^4 - 96*k^7*c1*x^4*y^4*dy^3 + 210*k^6*c1*x^3*y^4*dy^3 +
 - 150*k^5*c1*x^2*y^4*dy^3 - 80*k^4*x*y^4*dy^3 + 35*k^4*c1*x*y^4*dy^3 + 240*k^3*y^4*dy^3 + 42*k^8*c1*x^4*y^5*dy^2 +
 - 102*k^7*c1*x^3*y^5*dy^2 + 81*k^6*c1*x^2*y^5*dy^2 + 72*k^5*x*y^5*dy^2 - 21*k^5*c1*x*y^5*dy^2 +
 - 180*k^4*y^5*dy^2 - 10*k^9*c1*x^4*y^6*dy + 27*k^8*c1*x^3*y^6*dy - 24*k^7*c1*x^2*y^6*dy - 28*k^6*x*y^6*dy +
 + 7*k^6*c1*x*y^6*dy + 72*k^5*y^6*dy + k^10*c1*x^4*y^7 - 3*k^9*c1*x^3*y^7 + 3*k^8*c1*x^2*y^7 + 4*k^7*x*y^7 +
 - k^7*c1*x*y^7 - 12*k^6*y^7 # = 0

# TODO: format;


#########################
#########################

### Type: 2 Components w. Log * Exp
# y = Log(x)^2 * Exp(k*x) + c1*Log(x)*Exp(k*x)

# Check:
k  = 1/sqrt(5);
c1 = 2^(1/3); # c1 = 2i;
x = sqrt(3); params = list(x=x, k=k, c1=c1);
e = expression(log(x)^2 * exp(k*x) + c1*log(x)*exp(k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - k*x*y - 2*log(x)*exp(k*x) - c1*exp(k*x) # = 0

# D2 =>
x^2*d2y - x*(k*x - 1)*dy - k*x*y +
	- 2*k*x*log(x)*exp(k*x) - (c1*k*x + 2)*exp(k*x) # = 0
x^2*d2y - x*(k*x - 1)*dy - k*x*y +
	- k*x* (x*dy - k*x*y - c1*exp(k*x)) - (c1*k*x + 2)*exp(k*x) # = 0
x^2*d2y - x*(2*k*x - 1)*dy + k*x*(k*x - 1)*y - 2*exp(k*x) # = 0

# Substitution Eq. D in initial Eq:
2*y - 2*log(x)^2 * exp(k*x) - c1*(x*dy - k*x*y - c1*exp(k*x)) # = 0
2*log(x)^2 * exp(k*x) + c1*x*dy - (c1*k*x + 2)*y - c1^2*exp(k*x) # = 0
4*log(x)^2 * exp(k*x) + 2*c1*x*dy - 2*(c1*k*x + 2)*y +
	- c1^2* (x^2*d2y - x*(2*k*x - 1)*dy + k*x*(k*x - 1)*y) # = 0

# from D2 =>
c1*x^2*d2y - c1*x*(k*x - 1)*dy - c1*k*x*y +
	- 2*c1*k*x*log(x)*exp(k*x) +
	- (c1*k*x + 2) * (x*dy - k*x*y - 2*log(x)*exp(k*x)) # = 0
c1*x^2 * d2y - x*(2*c1*k*x - (c1-2)) * dy +
	+ (c1*k*x - (c1-2)) * k*x*y + 4*log(x)*exp(k*x) # = 0

# =>
2*x^2*d2y - 2*x*(2*k*x - 1)*dy + 2*k*x*(k*x - 1)*y +
	+ (c1*x^2 * d2y - x*(2*c1*k*x - (c1-2)) * dy +
		+ (c1*k*x - (c1-2)) * k*x*y)^2 /
	(2*c1*x*dy - 2*(c1*k*x + 2)*y +
		- c1^2* (x^2*d2y - x*(2*k*x - 1)*dy + k*x*(k*x - 1)*y)) # = 0

### ODE:
c1^2*x^3*d2y^2 - 2*c1^2*x^2*(2*k*x - 1) * dy*d2y +
	+ 2*x*(c1^2*k^2*x^2 - c1^2*k*x + 4) * y*d2y +
	+ x*(4*c1^2*k^2*x^2 - 4*c1^2*k*x + c1^2-4) * dy^2 +
	- (4*c1^2*k^3*x^3 - 6*c1^2*k^2*x^2 + 2*(c1^2+4)*k*x - 8) * y*dy +
	+ (c1^2*k^4*x^3 - 2*c1^2*k^3*x^2 + (c1^2+4)*k^2*x - 8*k) * y^2 # = 0

### Special Cases:

### Case: c1 = 2i;
x^3*d2y^2 - 2*x^2*(2*k*x - 1) * dy*d2y +
	+ 2*x*(k^2*x^2 - k*x - 1) * y*d2y +
	+ 2*x*(2*k^2*x^2 - 2*k*x + 1) * dy^2 +
	- 2*(2*k^3*x^3 - 3*k^2*x^2 + 1) * y*dy +
	+ k*(k^3*x^3 - 2*k^2*x^2 + 2) * y^2 # = 0


# Derivation:

(2*x^2*d2y - 2*x*(2*k*x - 1)*dy + 2*k*x*(k*x - 1)*y) *
	(2*c1*x*dy - 2*(c1*k*x + 2)*y +
		- c1^2* (x^2*d2y - x*(2*k*x - 1)*dy + k*x*(k*x - 1)*y)) +
	+ (c1*x^2 * d2y - x*(2*c1*k*x - (c1-2)) * dy + (c1*k*x - (c1-2)) * k*x*y)^2 # = 0

# as.pm(pp) |> mult.pm(-1) |> sort.dpm()

