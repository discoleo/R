#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Power of Log
##
## draft v.0.1b


### Simple Power of Log:
# y = Log(P(x))^n


### Examples:

# x*y*d2y - n*x*dy^2 + y*dy = 0;
# n*x^2 * y*d2y - (n-1)*x^2 * dy^2 + (n-2*p)*x * y*dy + p^2 * y^2 = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y = Log(P(x))^n

# Check:
n = sqrt(5);
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


#########################

### y = x^p * Log(P(x))^n

# Check:
n = sqrt(5);
p = sqrt(2); # n = 2*p;
x = sqrt(3); params = list(x=x, n=n);
e = expression(x^p * log(x)^n)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - p*y - n*x^p * log(x)^(n-1) # = 0

# D2 =>
x^2 * d2y - (2*p-1)*x * dy + p^2 * y +
	- n*(n-1)*x^p * log(x)^(n-2) # = 0
n*x^2 * d2y - n*(2*p-1)*x * dy + n*p^2 * y +
	- (n-1)* (x*dy - p*y)^2 / y # = 0

### ODE:
n*x^2 * y*d2y - (n-1)*x^2 * dy^2 + (n-2*p)*x * y*dy + p^2 * y^2 # = 0


### Special Cases:

### Case: n = 2*p;
2*p*x^2 * y*d2y - (2*p-1)*x^2 * dy^2 + p^2 * y^2 # = 0
# Alternative:
m = 1/(2*p);
x^2 * y*d2y + (m-1)*x^2 * dy^2 + 1/(4*m) * y^2 # = 0

