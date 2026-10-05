#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Exp of Log^Power
##
## draft v.0.1b


### Exp of Log to Power:
# y = Exp( Log(P(x))^2 )

# Note: NO Generalization to higher powers;


### Examples:

# x^2 * y*d2y - x^2*dy^2 + x * y*dy - 2*k * y^2 = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y = Exp(k * Log(x)^2 )

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


#########################

### y = x^p * Exp(k * Log(x)^2 )
# Note: p has NO impact on ODE;

# Check:
k = 1/sqrt(5);
p = 1/sqrt(2);
x = sqrt(3); params = list(x=x, k=k, p=p);
e = expression(x^p * exp(k*log(x)^2))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - p*y - 2*k*x^p * log(x) * exp(k*log(x)^2) # = 0

# D2 =>
x^2 * d2y - (p-1)*x*dy - 2*k*y +
	- 2*p*k*x^p * log(x) * exp(k*log(x)^2) +
	- 4*k^2*x^p * log(x)^2 * exp(k*log(x)^2) # = 0
x^2 * d2y - (2*p-1)*x*dy + (p^2 - 2*k) * y +
	- 4*k^2*x^p * log(x)^2 * exp(k*log(x)^2) # = 0
x^2 * y*d2y - (2*p-1)*x * y*dy + (p^2 - 2*k) * y^2 +
	- (x*dy - p*y)^2 # = 0

### ODE:
x^2 * y*d2y - x^2 * dy^2 + x * y*dy - 2*k * y^2 # = 0

