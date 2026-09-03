########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log in Fractions
##
## draft v.0.1c


####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### y = 1 / (x + k*log(x))

# Check:
k = sqrt(2);
x = sqrt(3);
params = list(x=x, k=k);
e = expression(1 / (x + k*log(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy + (x+k) / (x + k*log(x))^2 # = 0
x*dy + (x+k)*y^2 # = 0

# D2 =>
x*d2y + 2*(x+k)*y*dy + y^2 + dy # = 0
# Substitution y^2 =>

### ODE:
x*(x+k) * d2y + 2*(x+k)^2 * y*dy + k*dy # = 0

##########################

### Partial Generalization
# y = x^p / (x + k*log(x))

# Check:
p = 1/3; # p = -1; # p = 1;
k = sqrt(2);
x = sqrt(3);
params = list(x=x, k=k, p=p);
e = expression(x^p / (x + k*log(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x^(p+1) * dy - p*x^p * y + (x+k)*y^2 # = 0

# D2 & Substitution =>

### ODE:
x^(p+1)*(x+k) * d2y + 2*(x+k)^2 * y*dy + k*x^p * dy +
	- (p^2-p)*x^p * y - k*p^2*x^(p-1) * y # = 0

### Special Cases:

# Note: for case p = 0 see a previous section;

# Case: p = 1;
x^2*(x+k) * d2y + 2*(x+k)^2 * y*dy + k*x * dy - k*y # = 0

# Case: p = -1;
x^2*(x+k) * d2y + 2*x^2*(x+k)^2 * y*dy + k*x * dy - (2*x + k) * y # = 0


##########################

### Partial Generalization
# y = x^p / (x^n + k*log(x))

# Check:
p = 1/3; # p = -1; # p = 1; # p = 0;
n = 7/5;
k = sqrt(2);    # k = n; # k = n = -1; p = 0;
# k = n; p = 1; # # k = n; p = 0;
# k = n*(n-1); # k = n*(n-1); p = 1;
x = sqrt(3);
params = list(x=x, p=p, n=n, k=k);
e = expression(x^p / (x^n + k*log(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x^(p+1) * dy - p*x^p * y + (n*x^n + k) * y^2 # = 0

# D2 =>
x^(p+1) * d2y + x^p * dy +
	+ 2*(n*x^n + k) * y*dy + n^2*x^(n-1) * y^2 - p^2*x^(p-1) * y # = 0

### ODE:
x^(p+1)*(n*x^n + k) * d2y + 2*(n*x^n + k)^2 * y*dy +
	- x^p*(n*(n-1)*x^n - k) * dy +
	+ p*x^(p-1)*(n*(n-p)*x^n - p*k) * y # = 0

### Special Cases:

# Case: k = n;
x^(p+1)*(x^n + 1) * d2y + 2*n*(x^n + 1)^2 * y*dy +
	- x^p*((n-1)*x^n - 1) * dy +
	+ p*x^(p-1)*((n-p)*x^n - p) * y # = 0

# Case: k = n, p = 1;
x^2*(x^n + 1) * d2y + 2*n*(x^n + 1)^2 * y*dy +
	- x*((n-1)*x^n - 1) * dy + ((n-1)*x^n - 1) * y # = 0

# Case: k = n*(n-1);
x^(p+1)*(x^n + n-1) * d2y + 2*n*(x^n + n-1)^2 * y*dy +
	- (n-1)*x^p*(x^n - 1) * dy +
	+ p*x^(p-1)*((n-p)*x^n - p*(n-1)) * y # = 0

# Case: k = n*(n-1); p = 1;
x^2*(x^n + n-1) * d2y + 2*n*(x^n + n-1)^2 * y*dy +
	- (n-1)*x*(x^n - 1) * dy + (n-1)*(x^n - 1) * y # = 0

# Case: p = 0;
x^2*(n*x^n + k) * d2y + 2*x*(n*x^n + k)^2 * y*dy +
	- x*(n*(n-1)*x^n - k) * dy # = 0

# Case: p = 0; k = n;
x^2*(x^n + 1) * d2y + 2*n*x*(x^n + 1)^2 * y*dy +
	- x*((n-1)*x^n - 1) * dy # = 0
# p = 0; k = n = -1;
x^2*(x+1) * d2y - 2*(x+1)^2 * y*dy + x*(x+2) * dy # = 0


##########################

# y = x^p / (x^2+1) / (x^n + k*log(x^2+1))

# Check:
p = 1/3; # p = -1; # p = 1; # p = 0;
n = 4/5; # n = 0;
k = sqrt(2);    # k = n;
# n = 0; p = 0; # n = -2; p = 0; # n = 0; p = 2;
x = sqrt(3);
params = list(x=x, p=p, n=n, k=k);
e = expression(x^p / (x^2+1) / (x^n + k*log(x^2+1)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x^2+1)*dy - p*(x^2+1)*y + 2*x^2 * y +
	+ (x^2+1) * (n*x^n*(x^2+1) + 2*k*x^2) / x^p * y^2 # = 0
x*(x^2+1) * dy - ((p-2)*x^2 + p) * y +
	+ (x^2+1) * (n*x^(n-p+2) + n*x^(n-p) + 2*k*x^(2-p)) * y^2 # = 0

# D2 =>
x*(x^2+1) * d2y - ((p-5)*x^2 + p-1) * dy +
	+ 2*(x^2+1) * (n*x^(n+2) + n*x^n + 2*k*x^2) / x^p * y*dy +
	- 2*(p-2)*x * y +
	+ x*(n*(n-p+4)*x^(n+2) + 2*n*(n-p+2)*x^n + n*(n-p)*x^(n-2) +
		+ 2*(4-p)*k*x^2 + 2*k*(2-p)) / x^p * y^2 # = 0
# Substitution for y^2 =>

### ODE:
x*(x^2+1)^2 * (n*x^(n+2) + n*x^n + 2*k*x^2) * d2y +
	+ 2*(x^2+1)^2 * (n*x^(n+2) + n*x^n + 2*k*x^2)^2 / x^p * y*dy +
	- (n*(n-1)*(x^2+1)^3 * x^n - 2*k*x^2*(x^4-1)) * dy +
	+ x * (n*(p-2)*(n-p+2)*x^(n+4) - n*(3*p^2 - 3*n*p - 8*p + 4*n)*x^(n+2) +
		- n*(3*p^2 - 4*p - 3*n*p - 4 + 2*n)*x^n + n*p*(n-p)*x^(n-2) +
		- 2*k*(p-2)^2*x^4 - 4*k*p*(p-3)*x^2 - 2*k*p*(p-2)) * y # = 0

### Special Cases:

# Case: p = 0;
x*(x^2+1)^2 * (n*x^(n+2) + n*x^n + 2*k*x^2) * d2y +
	+ 2*(x^2+1)^2 * (n*x^(n+2) + n*x^n + 2*k*x^2)^2 * y*dy +
	- (n*(n-1)*(x^2+1)^3 * x^n - 2*k*x^2*(x^4-1)) * dy +
	- 2*x * (n*(n+2)*x^(n+4) + 2*n^2*x^(n+2) + n*(n-2)*x^n + 4*k*x^4) * y # = 0

# Case: p = 0; n = 0;
x*(x^2+1)^2 * d2y + 4*k*x^2*(x^2+1)^2 * y*dy + (x^4-1) * dy - 4*x^3 * y # = 0

# Case: p = 0; n = -2;
x^3*(x^2+1)^2 * (k*x^4 - x^2 - 1) * d2y +
	+ 4*(x^2+1)^2 * (k*x^4 - x^2 - 1)^2 * y*dy +
	+ x^2*(k*x^4*(x^4-1) - 3*(x^2+1)^3) * dy +
	- 4*x^3 * (k*x^6 + 2*x^2 + 2) * y # = 0

# Case: n = 0;
x^2*(x^2+1)^2 * d2y +
	+ 4*k*(x^2+1)^2 * x^3 / x^p * y*dy +
	+ x*(x^4-1) * dy +
	- ((p-2)^2*x^4 + 2*p*(p-3)*x^2 + p*(p-2)) * y # = 0

# Case: n = 0; p = 2;
x*(x^2+1)^2 * d2y + 4*k*(x^2+1)^2 * y*dy +
	+ (x^4-1) * dy + 4*x * y # = 0

