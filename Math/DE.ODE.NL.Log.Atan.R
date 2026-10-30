########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log of Atan
##
## draft v.0.1c


### Examples:

# (x^2+k^2) * d2y + (x^2+k^2) * dy^2 + 2*x * dy = 0;

# ODE-Type: d2y + dy^2 + y*dy
# - Usually with large coefficients;

# 4*(x-b)*x^2 * d2y + 4*(x-b)*x^3 * dy^2 + 8*(x-b)*x^2 * y*dy +
#   + 2*x*(7*x-6*b) * dy + 4*x*(x-b) * y^2 + 2*(3*x-2*b) * y = 0;

# Note:
# - For NL ODEs from y = Log(Log), see file:
#   DE.ODE.NL.Log.Log.R;

####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### y = x * log(atan(x/k))

# Check:
k = 1/5;
x = sqrt(3); params = list(x=x, k=k);
e = expression(x * log(atan(x/k)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*(x^2+k^2)*dy - (x^2+k^2)*y - k*x^2 / atan(x/k) # = 0

# D2 =>
x*(x^2+k^2)*d2y + 2*x^2*dy - 2*x*y +
	- 2*k*x / atan(x/k) + k^2*x^2 / (x^2+k^2) / atan(x/k)^2 # = 0
x^2*(x^2+k^2)*d2y - 2*k^2*x*dy + 2*k^2*y +
	+ k^2*x^3 / (x^2+k^2) / atan(x/k)^2 # = 0
x^3*(x^2+k^2) * d2y - 2*k^2*x^2 * dy + 2*k^2*x * y +
	+ (x^2+k^2)*(x*dy - y)^2 # = 0

### ODE:
x^3*(x^2+k^2) * d2y + x^2*(x^2+k^2) * dy^2 - 2*x*(x^2+k^2) * y*dy +
	+ (x^2+k^2) * y^2 - 2*k^2*x^2 * dy + 2*k^2*x * y # = 0


#######################

### Gen: y = x^p * log(atan(x^n/k))

# Check:
k = 1/5;
n = 1/sqrt(2); p = sqrt(3); # p = n;
# p = n = 1/2; # p = 0; # p = 0; n = 1; # n = 1/2;
x = sqrt(3); params = list(x=x, k=k);
e = expression(x^p * log(atan(x^n/k)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*(x^(2*n)+k^2)*dy - p*(x^(2*n)+k^2)*y - n*k*x^(n+p) / atan(x^n/k) # = 0

# D2 =>
x^2*(x^(2*n)+k^2) * d2y +
	+ x*((n-2*p+1)*x^(2*n) - (n+2*p-1)*k^2) * dy +
	- p*((n-p)*x^(2*n) - (n+p)*k^2) * y +
	+ n^2*k^2*x^(2*n+p) / (x^(2*n)+k^2) / atan(x^n/k)^2 # = 0
x^(p+2)*(x^(2*n)+k^2) * d2y +
	+ x^(p+1)*((n-2*p+1)*x^(2*n) - (n+2*p-1)*k^2) * dy +
	- p*x^p*((n-p)*x^(2*n) - (n+p)*k^2) * y +
	+ (x^(2*n)+k^2) * (x*dy - p*y)^2 # = 0

### ODE:
x^(p+2)*(x^(2*n)+k^2) * d2y + x^2*(x^(2*n)+k^2) * dy^2 +
	- 2*p*x*(x^(2*n)+k^2) * y*dy +
	+ x^(p+1)*((n-2*p+1)*x^(2*n) - (n+2*p-1)*k^2) * dy +
	+ p^2*(x^(2*n)+k^2) * y^2 - p*x^p*((n-p)*x^(2*n) - (n+p)*k^2) * y # = 0


### Special Cases:

### Case: p = 0
x*(x^(2*n)+k^2) * d2y + x*(x^(2*n)+k^2) * dy^2 +
	+ ((n+1)*x^(2*n) - (n-1)*k^2) * dy # = 0
### Case: p = 0; n = 1;
(x^2+k^2) * d2y + (x^2+k^2) * dy^2 + 2*x * dy # = 0


### Case: p = n;
x^(n+2)*(x^(2*n)+k^2) * d2y + x^2*(x^(2*n)+k^2) * dy^2 +
	- 2*n*x*(x^(2*n)+k^2) * y*dy +
	- x^(n+1)*((n-1)*x^(2*n) + (3*n-1)*k^2) * dy +
	+ n^2*(x^(2*n)+k^2) * y^2 + 2*k^2*n^2*x^n * y # = 0
### Case: p = n = 1/2;
4*x^(5/2)*(x+k^2) * d2y + 4*x^2*(x+k^2) * dy^2 +
	- 4*x*(x+k^2) * y*dy + 2*x^(3/2)*(x - k^2) * dy +
	+ (x+k^2) * y^2 + 2*k^2*x^(1/2) * y # = 0

### Case: n = 1/2;
x^(p+2)*(x+k^2) * d2y + x^2*(x+k^2) * dy^2 - 2*p*x*(x+k^2) * y*dy +
	- x^(p+1)*((2*p-3/2)*x + (2*p-1/2)*k^2) * dy +
	+ p^2*(x+k^2) * y^2 + p*x^p*((p-1/2)*x + (p+1/2)*k^2) * y # = 0


#########################
#########################

# Type: Log(Atan( SQRT ))

### y = x^p * log(atan(sqrt(k*x - 1)))
# Note: should have been sqrt(x/k - 1);

# Check:
k = sqrt(5); p = 1/sqrt(2); # p = -1;
x = sqrt(3); params = list(x=x, k=k); b = 1/k;
e = expression(x^p * log(atan(sqrt(k*x-1))))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
2*x*dy - 2*p*y - x^p / sqrt(k*x - 1) / atan(sqrt(k*x - 1))

# D2 =>
4*x^2*(k*x-1)*d2y - 4*(p-1)*x*(k*x-1)*dy +
	- x^p * (2*p*(k*x-1) - k*x) / sqrt(k*x - 1) / atan(sqrt(k*x - 1)) +
	+ (k*x-1)*x^p / sqrt(k*x - 1)^2 / atan(sqrt(k*x - 1))^2 # = 0
4*x^2*(k*x-1)*d2y - 2*x*((4*p-3)*k*x - 4*p+2) * dy +
	+ 2*p*((2*p-1)*k*x - 2*p) * y +
	+ (k*x-1)*x^p / sqrt(k*x - 1)^2 / atan(sqrt(k*x - 1))^2 # = 0
4*x^(p+2)*(k*x-1)*d2y - 2*x^(p+1)*((4*p-3)*k*x - 4*p+2) * dy +
	+ 2*p*((2*p-1)*k*x - 2*p)*x^p * y +
	+ 4*(k*x-1) * (x*dy - p*y)^2 # = 0

### ODE:
4*(k*x-1)*x^(p+2) * d2y + 4*(k*x-1)*x^2 * dy^2 +
	- 8*p*x*(k*x-1) * y*dy - 2*x^(p+1)*((4*p-3)*k*x - 4*p+2) * dy +
	+ 4*p^2*(k*x-1) * y^2 + 2*p*((2*p-1)*k*x - 2*p)*x^p * y # = 0

# Reformulated:
b = 1/k;
4*(x-b)*x^(p+2) * d2y + 4*(x-b)*x^2 * dy^2 +
	- 8*p*x*(x-b) * y*dy - 2*((4*p-3)*x - (4*p-2)*b)*x^(p+1) * dy +
	+ 4*p^2*(x-b) * y^2 + 2*p*((2*p-1)*x - 2*p*b)*x^p * y # = 0


#### Special Cases

### Case: p = -1;
4*(x-b)*x^2 * d2y + 4*(x-b)*x^3 * dy^2 +
	+ 8*(x-b)*x^2 * y*dy + 2*x*(7*x - 6*b) * dy +
	+ 4*(x-b)*x * y^2 + 2*(3*x - 2*b) * y # = 0

