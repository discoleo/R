########################$
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - SQRT( Log )
##
## draft v.0.1a


### NL ODEs: SQRT-Type
#   Subtype: SQRT( Log )


### Examples:

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

