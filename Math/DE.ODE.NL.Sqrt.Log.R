########################$
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - SQRT( Log )
##
## draft v.0.1b


### NL ODEs: SQRT-Type
#   Subtype: SQRT( Log )


### Examples:

# Order 1:
# x^2*y^2*dy^2 - y^2 + b = 0;

# Order 2:
# x^2*y^3*d2y + x*y^3*dy - b = 0;


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


########################

### Simple Components
# (Non-Coupled)

### y = sqrt(log(x)^2 + b)

# Check:
e = expression(sqrt(log(x)^2 + b));
b = sqrt(5);
x = sqrt(3);
params = list(x=x, b=b);
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*dy - log(x) / y # = 0
x*y*dy - log(x) # = 0
#
x^2*y^2*dy^2 - y^2 + b # = 0

# D2 =>
x^2*y*d2y + x^2*dy^2 + x*y*dy - 1 # = 0


### Variant 1:
x^2*y^3*d2y + x*y^3*dy - b # = 0

### Variant 2:
x^2*y^2*dy^2 - y^2 + b # = 0
x^2*y^2*dy*d2y + x^2*y*dy^3 + x*y^2*dy^2 - y*dy # = 0
(y^2-b)*d2y + x^2*y*dy^4 + x*y^2*dy^3 - y*dy^2 # = 0
x*(y^2-b)*d2y + x^3*y*dy^4 - x*y*dy^2 + (y^2-b)*dy # = 0


#########################

### y = log(x) * sqrt(log(x)^2 + b)

# Check:
e = expression(log(x) * sqrt(log(x)^2 + b));
b = sqrt(5);
x = sqrt(3);
params = list(x=x, b=b);
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*dy - y / log(x) - log(x)^2 / sqrt(log(x)^2 + b) # = 0
x*dy - y / log(x) - log(x)^3 / y # = 0
x*log(x) * y*dy - y^2 - log(x)^4 # = 0
x*log(x) * y*dy - 2*y^2 + b*log(x)^2 # = 0

### ODE Order 1:
# Solve:
x^4*y^2*dy^4 + 4*b*x^2*y^2*dy^2 + b^3*x^2*dy^2 - 16*y^4 - 8*b^2*y^2 - b^4 # = 0

# D2 =>
x^2*log(x) * y*d2y + x^2*log(x) * dy^2 +
	+ x*(log(x) - 3) * y*dy + 2*b*log(x) # = 0
# Substitute:
x^2*(2*x^2 * y^2*dy^2 + 4*b*y^2 + b^3) * y*d2y - x^4*y^2*dy^4 + 2*x^3*y^3*dy^3 +
 - 4*b*x^2*y^2*dy^2 - 2*b^3*x^2*dy^2 + 4*b*x*y^3*dy + b^3*x*y*dy + 8*b^2*y^2 + 2*b^4 # = 0

### ODE:
x^2*(2*x^2 * y^2*dy^2 + 4*b*y^2 + b^3) * y*d2y + 2*x^3 * y^3*dy^3 +
	- b^3*x^2*dy^2 + 4*b*x*y^3*dy + b^3*x*y*dy - 16*y^4 + b^4 # = 0
