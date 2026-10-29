########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log Compositions
## Log( Radical(Exp) )
##
## draft v.0.1c


### Examples:

# x*d2y + k*dy^3 - 1/k * (x^2 + k)*dy = 0;

# x^5*d2y + (x*dy - y)^3 - (x + 2)*x^3 * (x*dy - y) = 0;
# x^5*d2y + k*(x*dy - y)^3 - 1/k * (x + 2*k)*x^3 * (x*dy - y) = 0;
# x^2*d2y - (x*dy - p*y)^3 + (x^(-2*p) + p) * (x*dy - p*y) - (p-1)*x * dy = 0;


####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### y = x * log(sqrt(exp(2*x) + b0) - exp(x))
# Note: Very simple variant;

# Check:
b0 = 2/3;
# Note: b0 has no impact on ODE;
x = sqrt(3); params = list(x=x, b0=b0);
e = expression(x * log(sqrt(exp(2*x) + b0) - exp(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*dy - y + x^2 * exp(x) / sqrt(exp(2*x) + b0) # = 0

# D2 =>
x*d2y + (x^2 + 2*x) * exp(x) / sqrt(exp(2*x) + b0) +
	- x^2 * exp(3*x) / sqrt(exp(2*x) + b0)^3 # = 0

### ODE:
x^5*d2y + (x*dy - y)^3 - (x + 2)*x^3 * (x*dy - y) # = 0


#######################

### y = x * log(sqrt(exp(2*k*x) + b0) - exp(k*x))

# Check:
# Note: b0 has no impact on ODE;
b0 = 2/3;
k = 1/sqrt(3); # k = 2;
x = sqrt(3); params = list(x=x, k=k, b0=b0);
e = expression(x * log(sqrt(exp(2*k*x) + b0) - exp(k*x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*dy - y + k*x^2 * exp(k*x) / sqrt(exp(2*k*x) + b0) # = 0

# D2 =>
x*d2y + k*(k*x^2 + 2*x) * exp(k*x) / sqrt(exp(2*k*x) + b0) +
	- k^2*x^2 * exp(3*k*x) / sqrt(exp(2*k*x) + b0)^3 # = 0

### ODE:
k*x^5*d2y + (x*dy - y)^3 - k*(k*x + 2)*x^3 * (x*dy - y) # = 0
# Re-formulated:
k = 1/k;
x^5*d2y + k*(x*dy - y)^3 - 1/k * (x + 2*k)*x^3 * (x*dy - y) # = 0

### Special Cases:

### Case: k = 2;
2*x^5*d2y + (x*dy - y)^3 - 4*(x+1)*x^3 * (x*dy - y) # = 0


#######################

### Gen: y = x^p * log(sqrt(exp(2*k*x^n) + b0) - exp(k*x^n))

# Check:
k = 1/sqrt(3); # k = 2;
n = 1/sqrt(5); p = -sqrt(2);
# p = 0; # p = 0; n = 2;  # n = -2*p; # n = -2*p; k = -1/n;
# p = 1; n = -2; k = 1/2; # p = -1; n = 2; k = -1/2;
b0 = 2/3;
# Note: b0 has no impact on ODE;
x = sqrt(3); params = list(x=x, k=k, b0=b0);
e = expression(x^p * log(sqrt(exp(2*k*x^n) + b0) - exp(k*x^n)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*dy - p*y + n*k*x^(n+p) * exp(k*x^n) / sqrt(exp(2*k*x^n) + b0) # = 0

# D2 =>
x*d2y - (p-1)*dy + n*k*x^(n+p-1)*(n*k*x^n + n+p) * exp(k*x^n) / sqrt(exp(2*k*x^n) + b0) +
	- n^2*k^2*x^(2*n+p-1) * exp(3*k*x^n) / sqrt(exp(2*k*x^n) + b0)^3 # = 0
x^2*d2y - (p-1)*x*dy - (n*k*x^n + n+p) * (x*dy - p*y) +
	- n^2*k^2*x^(2*n+p) * exp(3*k*x^n) / sqrt(exp(2*k*x^n) + b0)^3 # = 0

### ODE:
n*k*x^(n+2*p+2) * d2y + (x*dy - p*y)^3 - n*k*(p-1)*x^(n+2*p+1) * dy +
	- n*k*(n*k*x^n + n+p)*x^(n+2*p) * (x*dy - p*y) # = 0

### Special Cases:

### Case: p = 0;
n*k*x^(n+2) * d2y + x^3 * dy^3 - n*k*(n*k*x^n + n-1)*x^(n+1) * dy # = 0
### Case: p = 0; n = 2;
2*k*x*d2y + dy^3 - 2*k*(2*k*x^2 + 1)*dy # = 0
Ki = 1/(2*k);
x*d2y + Ki*dy^3 - 1/Ki * (x^2 + Ki)*dy # = 0

### Case: n = -2*p;
2*p*k*x^2 * d2y - (x*dy - p*y)^3 - 2*p*k*(p-1)*x * dy +
	+ 2*p^2*k*(2*k*x^n + 1) * (x*dy - p*y) # = 0
### Case: n = -2*p; k = 1/(2*p);
x^2 * d2y - (x*dy - p*y)^3 - (p-1)*x * dy + (x^(-2*p) + p) * (x*dy - p*y) # = 0
### Case: p = 1; n = -2; k = 1/2;
x^4 * d2y - x^2 * (x*dy - y)^3 + (x^2+1) * (x*dy - y) # = 0
### Case: p = -1; n = 2; k = -1/2;
x^2 * d2y - (x*dy + y)^3 + (x^2-1) * (x*dy + y) + 2*x*dy # = 0

