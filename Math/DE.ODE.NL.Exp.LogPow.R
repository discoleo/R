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


### Examples:

# x^2 * y*d2y - x^2*dy^2 + x * y*dy - 2*k * y^2 = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y = Exp(k * Log(x)^2 )
# - For a Generalization of the Power,
#   see the section on Higher Powers;

# Check:
k = 1/sqrt(5);
x = sqrt(3); params = list(x=x, k=k);
e = expression(exp(k*log(x)^2))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - 2*k*log(x) * exp(k*log(x)^2) # = 0

# D2 =>
x^2 * d2y + x*dy - 2*k*y - 4*k^2 * log(x)^2 * exp(k*log(x)^2) # = 0
x^2 * d2y + x*dy - 2*k*y - x^2*dy^2 / y # = 0

### ODE:
x^2 * y*d2y - x^2*dy^2 + x * y*dy - 2*k * y^2 # = 0

