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

# TODO: Elimination;

