#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Power of Log
##
## draft v.0.1a


### Exp of Atan to Power:
# y = Log(P(x))^n


### Examples:

# x*y*d2y - n*x*dy^2 + y*dy = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y = Log(P(x))^n

# Check:
n = sqrt(5)
x = sqrt(3); params = list(x=x, n=n);
e = expression(log(x)^n)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - n*log(x)^(n-1) # = 0

# D2 =>
x^2*d2y + x*dy - n*(n-1)*log(x)^(n-2) # = 0

### ODE:
n*x * y*d2y - (n-1)*x * dy^2 + n*y*dy # = 0

# Alternative:
m = 1 - 1/n; # n = 1/(1-m);
x*y*d2y - m*x*dy^2 + y*dy # = 0

