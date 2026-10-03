#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Power of Log(SQRT)
##
## draft v.0.1c

### Log to Power:
# y = Sqrt(P(x)) * Log(Sqrt(P(x)^2 + b) - P(x))^n


### Examples:

# 2*(x^2+1)^2 * y*d2y - (x^2+1)^2 * dy^2 + (x^2-2) * y^2 = 0;
# 3*(x^2+1)^2 * y*d2y - 2*(x^2+1)^2 * dy^2 + x*(x^2+1) * y*dy + (x^2-3) * y^2 = 0;
# n*(x^2+k)^2 * y*d2y - (n-1)*(x^2+k)^2 * dy^2 + (n-2)*x*(x^2+k) * y*dy + (x^2-k*n) * y^2 = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y = sqrt(x^2+1) * log(sqrt(x^2+1) - x)^2

# Check:
k = 1/sqrt(5);
x = sqrt(3); params = list(x=x, k=k);
e = expression(sqrt(x^2+1) * log(sqrt(x^2+1) - x)^2)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
dy - x*y/(x^2+1) + 2*log(sqrt(x^2+1) - x) # = 0
(x^2+1)*dy - x*y + 2*(x^2+1)*log(sqrt(x^2+1) - x) # = 0

# D2 =>
(x^2+1)*d2y + x*dy - y +
	+ 4*x*log(sqrt(x^2+1) - x) - 2*sqrt(x^2+1) # = 0
(x^2+1)^2*d2y + x*(x^2+1)*dy - (x^2+1)*y +
	- 2*x*((x^2+1)*dy - x*y) - 2*(x^2+1)*sqrt(x^2+1) # = 0
2*(x^2+1)^2*d2y - 2*x*(x^2+1)*dy + 2*(x^2-1)*y +
	- ((x^2+1)*dy - x*y)^2 / y # = 0

### ODE:
2*(x^2+1)^2 * y*d2y - (x^2+1)^2 * dy^2 + (x^2-2)*y^2 # = 0


#########################

### y = sqrt(x^2+1) * log(sqrt(x^2+1) - x)^3

# Check:
x = sqrt(3); params = list(x=x);
e = expression(sqrt(x^2+1) * log(sqrt(x^2+1) - x)^3)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2+1)*dy - x*y + 3*(x^2+1)*log(sqrt(x^2+1) - x)^2 # = 0

# D2 =>
(x^2+1)*d2y + x*dy - y +
	+ 6*x*log(sqrt(x^2+1) - x)^2 - 6*sqrt(x^2+1)*log(sqrt(x^2+1) - x) # = 0
(x^2+1)^2*d2y - x*(x^2+1)*dy + (x^2-1)*y +
	- 6*(x^2+1)*sqrt(x^2+1)*log(sqrt(x^2+1) - x) # = 0
3*(x^2+1)^2*d2y - 3*x*(x^2+1)*dy + 3*(x^2-1)*y +
	- 2*((x^2+1)*dy - x*y)^2 / y # = 0

### ODE:
3*(x^2+1)^2 * y*d2y - 2*(x^2+1)^2 * dy^2 + x*(x^2+1) * y*dy + (x^2-3) * y^2 # = 0


#####################

### Gen: y = sqrt(x^2+1) * log(sqrt(x^2+1) - x)^n

# Check:
k = sqrt(11); # k = sqrt(13);
n = exp(1);
x = 1/sqrt(3); params = list(x=x, n=n);
e = expression(sqrt(x^2+k) * log(sqrt(x^2+k) - x)^n)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2+k)*dy - x*y + n*(x^2+k)*log(sqrt(x^2+k) - x)^(n-1) # = 0

# D2 =>
(x^2+k)*d2y + x*dy - y +
	+ 2*n*x*log(sqrt(x^2+k) - x)^(n-1) +
	- n*(n-1)*sqrt(x^2+k)*log(sqrt(x^2+k) - x)^(n-2) # = 0

### ODE:
n*(x^2+k)^2 * y*d2y - (n-1)*(x^2+k)^2 * dy^2 + (n-2)*x*(x^2+k) * y*dy + (x^2-k*n) * y^2 # = 0

