#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log * Exp
##
## draft v.0.1a

### Log to Power:
# y = Log(P1(x))^2 * Exp(P2(x)})


### Examples:

# 2*x*y*d2y - x*dy^2 - 2*(k*x - 1) * y*dy + k*(k*x - 2) * y^2 = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y = log(x)^2 * Exp(k*x)

# Check:
k = 1/sqrt(5);
x = sqrt(3); params = list(x=x);
e = expression(log(x)^2 * exp(k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - k*x*y - 2*log(x)*exp(k*x) # = 0

# D2 =>
x*d2y + dy - k*x*dy - k*y +
	- 2*k*log(x)*exp(k*x) - 2*exp(k*x)/x # = 0
x^2*d2y - x*(2*k*x - 1)*dy + k*x*(k*x - 1)*y - 2*exp(k*x) # = 0
2*x^2*d2y - 2*x*(2*k*x - 1)*dy + 2*k*x*(k*x - 1)*y - (x*dy - k*x*y)^2 / y # = 0
2*x^2 * y*d2y - 2*x*(2*k*x - 1) * y*dy + 2*k*x*(k*x - 1) * y^2 +
	- (x*dy - k*x*y)^2 # = 0

### ODE:
2*x*y*d2y - x*dy^2 - 2*(k*x - 1) * y*dy + k*(k*x - 2) * y^2 # = 0

