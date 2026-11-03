########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log( SQRT(Log) )
##
## draft v.0.1b


### NL ODEs: Log-Type
#   Subtype: Log( SQRT(Log) )


### Examples:

# x^4*(x^2*d2y - x*dy + y)^2 + b*(x*dy - y)^6 - x^2 * (x*dy - y)^4 = 0;


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


########################

### Simple Component

### y = x * log(sqrt(log(x)^2 + b) - log(x))

# Check:
e = expression(x * log(sqrt(log(x)^2 + b) - log(x)));
b = sqrt(5);
x = sqrt(3);
params = list(x=x, b=b);
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*dy - y + x / sqrt(log(x)^2 + b) # = 0
#
log(x)^2 + b - x^2 / (x*dy - y)^2 # = 0

# D2 =>
x*d2y - (x*dy - y) / x - log(x) / sqrt(log(x)^2 + b)^3 # = 0
x^4*d2y - x^3*dy + x^2*y + log(x) * (x*dy - y)^3 # = 0
(x^4*d2y - x^3*dy + x^2*y)^2 - (x^2 - b*(x*dy - y)^2) * (x*dy - y)^4 # = 0
(x^4*d2y - x^3*dy + x^2*y)^2 + b*(x*dy - y)^6 - x^2 * (x*dy - y)^4 # = 0

### ODE:
x^4*(x^2*d2y - x*dy + y)^2 + b*(x*dy - y)^6 - x^2 * (x*dy - y)^4 # = 0


########################

### Coupled Component

### y = x * log(sqrt(log(x)^2 + b) - log(x)) * sqrt(log(x)^2 + b)

# Check:
e = expression(x * log(sqrt(log(x)^2 + b) - log(x)) * sqrt(log(x)^2 + b));
b = sqrt(5);
x = sqrt(3);
params = list(x=x, b=b);
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*dy - y - log(x) / (log(x)^2 + b) * y + x # = 0
(x*dy - y + x)*(log(x)^2 + b) - log(x) * y # = 0

# D2 =>
x^2*d2y - x*log(x) / (log(x)^2 + b) * dy +
	- 1 / (log(x)^2 + b) * y +
	+ 2*log(x)^2 / (log(x)^2 + b)^2 * y + x # = 0
x^2 * y*d2y - x*(x*dy - y + x) * dy +
	- (x*dy - y + x) / log(x) * y +
	+ 2*(x*dy - y + x)^2 + x*y # = 0


### ODE:
b*x^4*(x*dy - y + x) * y^2*d2y^2 +
	+ (2*b*x^5 * dy^3 - 8*b*x^4 * (y - x)*dy^2 +
		+ x^3*((10*b-1)*y^2 - 18*b*x*y + 10*b*x^2) * dy +
		- (4*b-1)*x^2*y^3 + (10*b-1)*x^3*y^2 - 10*b*x^4*y + 4*b*x^5) * y*d2y +
 + b*x^5*dy^5 - 7*b*x^4*y*dy^4 + 7*b*x^5*dy^4 + 19*b*x^3*y^2*dy^3 - 36*b*x^4*y*dy^3 + 19*b*x^5*dy^3 +
 + x^2*y^3*dy^2 - 25*b*x^2*y^3*dy^2 - x^3*y^2*dy^2 + 67*b*x^3*y^2*dy^2 - 67*b*x^4*y*dy^2 + 25*b*x^5*dy^2 +
 + ((16*b-2)*x*y^4 - (54*b-3)*x^2*y^3 + (77*b-2)*x^3*y^2 - 54*b*x^4*y + 16*b*x^5) * dy +
 - (4*b-1)*y^5 + (16*b-2)*x*y^4 - (29*b-2)*x^2*y^3 + (29*b-1)*x^3*y^2 - 16*b*x^4*y + 4*b*x^5 # = 0

# TODO: format;

