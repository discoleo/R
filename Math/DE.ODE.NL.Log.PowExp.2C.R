#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs: Log * Exp
## w. 2 Coupled Components
##
## draft v.0.1b

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

