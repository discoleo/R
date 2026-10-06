#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Exp of Atan^2
##
## draft v.0.1a


### Exp of Atan to Power:
# y = Exp( Atan(P(x))^2 )

# Note: NO Generalization yet to higher powers;


### Examples:

(x^2+k^2)^2 * y*d2y - (x^2+k^2)^2 * dy^2 + 2*x*(x^2+k^2) * y*dy - 2*a*k^2 * y^2 # = 0


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y = Exp(k * Atan(x)^2 )

# Check:
k = 1/sqrt(5); a = exp(1/3);
x = sqrt(3); params = list(x=x, k=k);
e = expression(exp(a*atan(x/k)^2))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2+k^2)*dy - 2*a*k * atan(x/k) * exp(a*atan(x/k)^2) # = 0

# D2 =>
(x^2+k^2)^2 * d2y + 2*x*(x^2+k^2) * dy +
	- 2*a*k^2 * exp(a*atan(x/k)^2) +
	- 4*a^2*k^2 * atan(x/k)^2 * exp(a*atan(x/k)^2) # = 0

### ODE:
(x^2+k^2)^2 * y*d2y - (x^2+k^2)^2 * dy^2 + 2*x*(x^2+k^2) * y*dy - 2*a*k^2 * y^2 # = 0

