#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Atan: Power
##
## draft v.0.1b

### Sum of Radicals to Power:
# y = B(x) * Atan(P(x))^n


### Examples:

# 3*x*(x^2+k^2) * y*d2y - 2*x*(x^2+k^2) * dy^2 + 6*(2*x^2 + k^2) * y*dy + 18*x * y^2 = 0;
# n*x*(x^2+k^2) * y*d2y - (n-1)*x*(x^2+k^2) * dy^2 + 2*n*(2*x^2 + k^2) * y*dy + 2*n^2*x * y^2 = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y = x^p * atan(x/k)^2
# Note: Linear ODE;

# Check:
k = 1/sqrt(5);
p = 1/sqrt(3);
x = sqrt(3); params = list(x=x, k=k, p=p);
e = expression(x^p * atan(x/k)^2)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x^2+k^2)*dy - p*(x^2+k^2)*y - 2*k*x^(p+1) * atan(x/k) # = 0

# D2 =>
x*(x^2+k^2)*d2y - ((p-3)*x^2 + (p-1)*k^2)*dy - 2*p*x*y +
	- 2*k*(p+1)*x^p * atan(x/k) - 2*k^2*x^(p+1) / (x^2+k^2) # = 0

### ODE:
x^2*(x^2+k^2)^2 * d2y - x*(x^2+k^2)*((2*p-2)*x^2 + 2*p*k^2)*dy +
	+ p*(x^2+k^2)*((p-1)*x^2 + (p+1)*k^2)*y - 2*k^2*x^(p+2) # = 0


#########################

### y = x^p * atan(x/k)^3
# Note: Linear ODE;

# Check:
k = 1/sqrt(5);
p = 1/sqrt(3); # p = 1; # p = -3;
x = sqrt(3); params = list(x=x, k=k, p=p);
e = expression(x^p * atan(x/k)^3)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x^2+k^2)*dy - p*(x^2+k^2)*y - 3*k*x^(p+1) * atan(x/k)^2 # = 0

# D2 =>
x*(x^2+k^2)*d2y - ((p-3)*x^2 + (p-1)*k^2)*dy - 2*p*x*y +
	- 3*k*(p+1)*x^p * atan(x/k)^2 - 6*k^2*x^(p+1) / (x^2+k^2) * atan(x/k) # = 0
x^2*(x^2+k^2)^2 * d2y - x*(x^2+k^2)*((2*p-2)*x^2 + 2*p*k^2) * dy +
	+ p*(x^2+k^2)*((p-1)*x^2 + (p+1)*k^2) * y - 6*k^2*x^(p+2) * atan(x/k) # = 0
3*x^2*(x^2+k^2)^2 * d2y - 3*x*(x^2+k^2)*((2*p-2)*x^2 + 2*p*k^2) * dy +
	+ 3*p*(x^2+k^2)*((p-1)*x^2 + (p+1)*k^2) * y +
	- 2 * (x*(x^2+k^2)*dy - p*(x^2+k^2)*y)^2 / y # = 0

### ODE:
3*x^2*(x^2+k^2) * y*d2y - 2*x^2*(x^2+k^2) * dy^2 +
	- 2*x*((p-3)*x^2 + p*k^2) * y*dy +
	+ p*((p-3)*x^2 + (p+3)*k^2) * y^2 # = 0

### Special Cases:

### Case: p = 1;
3*x^2*(x^2+k^2) * y*d2y - 2*x^2*(x^2+k^2) * dy^2 +
	+ 2*x*(2*x^2 - k^2) * y*dy - 2*(x^2 - 2*k^2) * y^2 # = 0

### Case: p = -3;
3*x*(x^2+k^2) * y*d2y - 2*x*(x^2+k^2) * dy^2 +
	+ 6*(2*x^2 + k^2) * y*dy + 18*x * y^2 # = 0


#########################

### y = x^p * atan(x/k)^n
# Note: Linear ODE;

# Check:
n = exp(1) - 1;
k = 1/sqrt(5);
p = 1/sqrt(3); # p = 1; # p = - n;
x = sqrt(3); params = list(x=x, k=k, n=n, p=p);
e = expression(x^p * atan(x/k)^n)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x^2+k^2)*dy - p*(x^2+k^2)*y - n*k*x^(p+1) * atan(x/k)^(n-1) # = 0

# D2 =>
x*(x^2+k^2)*d2y - ((p-3)*x^2 + (p-1)*k^2)*dy - 2*p*x*y +
	- n*k*(p+1)*x^p * atan(x/k)^(n-1) - n*(n-1)*k^2*x^(p+1) / (x^2+k^2) * atan(x/k)^(n-2) # = 0
x^2*(x^2+k^2)^2 * d2y - x*(x^2+k^2)*((2*p-2)*x^2 + 2*p*k^2) * dy +
	+ p*(x^2+k^2)*((p-1)*x^2 + (p+1)*k^2) * y - n*(n-1)*k^2*x^(p+2) * atan(x/k)^(n-2) # = 0
n*x^2*(x^2+k^2)^2 * d2y - n*x*(x^2+k^2)*((2*p-2)*x^2 + 2*p*k^2) * dy +
	+ n*p*(x^2+k^2)*((p-1)*x^2 + (p+1)*k^2) * y +
	- (n-1) * (x*(x^2+k^2)*dy - p*(x^2+k^2)*y)^2 / y # = 0

### ODE:
n*x^2*(x^2+k^2) * y*d2y - (n-1)*x^2*(x^2+k^2) * dy^2 +
	- 2*x*((p-n)*x^2 + p*k^2) * y*dy +
	+ p*((p-n)*x^2 + (p+n)*k^2) * y^2 # = 0

### Special Cases:

### Case: p = - n;
n*x*(x^2+k^2) * y*d2y - (n-1)*x*(x^2+k^2) * dy^2 +
	+ 2*n*(2*x^2 + k^2) * y*dy + 2*n^2*x * y^2 # = 0

