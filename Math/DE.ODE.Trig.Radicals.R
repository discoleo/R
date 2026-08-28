########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Linear ODEs - Trig w. Radicals
##
## draft v.0.1c

### Linear ODE: Trig w. Radicals
# y = B1(x) * sin(SQRT(R(x))) + B2(x) * SQRT(R(x)) * cos(R(x));

### Theory
# D(Sin(SQRT(...))) => SQRT * Cos(SQRT(...));
# D(SQRT * Cos(SQRT(...))) => Sin(SQRT(...)) + SQRT * Cos(SQRT(...));

# ODE Type: Homogeneous;


####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### y = x * sin(k*SQRT(x)) + c1 * SQRT(x) * cos(k*SQRT(x));

# Check:
k = sqrt(2); c1 = -1/3;
# k = sqrt(2); c1 = 1/k; # k = 1i * sqrt(3); c1 = 1/k;
x = sqrt(3);
params = list(x=x, k=k, c1=c1);
e = expression(x * sin(k*sqrt(x)) + c1*sqrt(x)*cos(k*sqrt(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*x*dy + (c1*k-2) * x * sin(k*sqrt(x)) +
	- (k*x + c1) * sqrt(x) * cos(k*sqrt(x)) # = 0

# System:
(k*x + c1*(c1*k-1)) * x * sin(k*sqrt(x)) # ==
- 2*c1*x * dy + (k*x + c1) * y;
#
(k*x + c1*(c1*k-1)) * sqrt(x) * cos(k*sqrt(x)) # ==
2*x * dy + (c1*k-2) * y;

# D2 =>
4*x^2*d2y + 4*x*dy +
	+ (k*(k*x + c1) + 2*(c1*k-2)) * x * sin(k*sqrt(x)) +
	+ ((c1*k-5)*k*x - c1) * sqrt(x) * cos(k*sqrt(x)) # = 0
# Substitution =>

### ODE:
4*x^2*(k*x + c1*(c1*k-1)) * d2y +
	- 2*x * (3*k*x + c1*(c1*k-1)) * dy +
	+ (k^3*x^2 + k*(c1^2*k^2 - 3*c1*k + 6)*x + 2*c1*(c1*k - 1)) * y # = 0

### Special Cases:

# Case: c1*k = 1;
4*x^3 * d2y - 6*x^2 * dy + x*(k^2*x + 4) * y # = 0


#####################

### y = x * cos(k*SQRT(x)) + c1 * SQRT(x) * sin(k*SQRT(x));
# Basically: k => -k;

# Check:
k = sqrt(2); c1 = -1/3;
# k = sqrt(2); c1 = -1/k; # k = 1i * sqrt(3); c1 = -1/k;
x = sqrt(3);
params = list(x=x, k=k, c1=c1);
e = expression(x * cos(k*sqrt(x)) + c1*sqrt(x)*sin(k*sqrt(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*x*dy - (c1*k+2) * x * cos(k*sqrt(x)) +
	+ (k*x - c1) * sqrt(x) * sin(k*sqrt(x)) # = 0

### ODE:
4*x^2*(k*x + c1*(c1*k+1)) * d2y +
	- 2*x * (3*k*x + c1*(c1*k+1)) * dy +
	+ (k^3*x^2 + k*(c1^2*k^2 + 3*c1*k + 6)*x + 2*c1*(c1*k + 1)) * y # = 0

### Special Cases:

# Case: c1*k = -1;
4*x^3 * d2y - 6*x^2 * dy + x*(k^2*x + 4) * y # = 0


#######################

### y = x * sin(k*SQRT(x^2+b)) + c1 * SQRT(x^2+b) * cos(k*SQRT(x^2+b));

# Note:
# - the "x"-Coefficient does NOT couple well with "(x^2+b)";
# - should probably try "x^0" or "x^2";

# Check:
c1 = -1/3; # c1 = 0;
k = sqrt(2);
x = sqrt(3); b = -1/sqrt(7);
params = list(x=x, k=k, c1=c1);
e = expression(x * sin(k*sqrt(x^2+b)) + c1*sqrt(x^2+b)*cos(k*sqrt(x^2+b)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x^2+b) * dy + (k*c1*x - 1)*(x^2+b) * x * sin(k*sqrt(x^2+b)) +
	- (k*x^3 + c1*x^2) * sqrt(x^2+b) * cos(k*sqrt(x^2+b)) # = 0
c1*x*(x^2+b) * dy - (k*x^3 + c1*x^2) * y +
	+ (k*c1^2*x*(x^2+b) + k*x^3 - c1*b) * x * sin(k*sqrt(x^2+b)) # = 0

# System:
(c1*(k*c1*x - 1)*(x^2+b) + (k*x^3 + c1*x^2)) * x * sin(k*sqrt(x^2+b)) # ==
- c1*x*(x^2+b) * dy + (k*x^3 + c1*x^2) * y;
#
(c1*(k*c1*x - 1)*(x^2+b) + (k*x^3 + c1*x^2)) * sqrt(x^2+b) * cos(k*sqrt(x^2+b)) # ==
x*(x^2+b) * dy + (k*c1*x - 1)*(x^2+b) * y;
# Note: Div = (c1^2+1)*k*x^3 + c1^2*k*b*x - c1*b;

# D2 =>
c1*x^2*(x^2+b)^2 * d2y - x*(x^2+b)*(k*x^3 - 2*c1*x^2 - b*c1) * dy +
	- x*(x^2+b)*(3*k*x^2 + 2*c1*x) * y +
	+ (x^2+b)*(k*c1^2*x*(4*x^2+2*b) + 4*k*x^3 - c1*b) * x * sin(k*sqrt(x^2+b)) +
	+ k*x^3*(k*c1^2*x*(x^2+b) + k*x^3 - c1*b) * sqrt(x^2+b) * cos(k*sqrt(x^2+b)) # = 0
c1^2*x^2*(x^2+b)^2 * d2y +
	- c1*x*(x^2+b)*(k*x^3 - 2*c1*x^2 - b*c1) * dy +
	+ x*(k^2*(c1^2+1)*x^5 - 3*c1*k*x^4 - 2*c1^2*x^3 + b*c1^2*k^2*x^3 +
		- 4*b*c1*k*x^2 - 2*b*c1^2*x) * y +
	+ (c1*(x^2+b)*(k*c1^2*x*(4*x^2+2*b) + 4*k*x^3 - c1*b) +
		- k*x^3*(k*c1^2*x*(x^2+b) + k*x^3 - c1*b)) * x * sin(k*sqrt(x^2+b)) # = 0
# =>

### ODE:
(x^2+b) * ((c1^2+1)*k*x^3 + c1^2*k*b*x - c1*b) * d2y +
	- (2*k*(c1^2+1)*x^4 + 3*k*b*(c1^2+1)*x^2 + b*c1*x + k*b^2*c1^2) * dy +
	+ (k^3*(c1^2+1)*x^5 + 2*k*(c1^2+1)*x^3 + k^3*b*c1^2*x^3 +
		- 3*c1*b*k^2*x^2 + 3*b*k*x + c1*b) * y # = 0

### Special Cases:

# Case: c1 = 0;
x^2*(x^2+b) * d2y - x*(2*x^2 + 3*b) * dy +
	+ (k^2*x^4 + 2*x^2 + 3*b) * y # = 0


#######################

### y = sin(k*SQRT(x^2+b)) + c1 * SQRT(x^2+b) * cos(k*SQRT(x^2+b));

# Check:
c1 = -1/3; # c1 = 0;
k = sqrt(2); b = -1/sqrt(7);
# c1 = k = b = 1; # c1 = -1; k = 1; b = 1/3;
x = sqrt(3);
params = list(x=x, k=k, c1=c1);
e = expression(sin(k*sqrt(x^2+b)) + c1*sqrt(x^2+b)*cos(k*sqrt(x^2+b)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
dy - (k+c1)*x * sqrt(x^2+b)*cos(k*sqrt(x^2+b)) / (x^2+b) +
	+ c1*k*x * sin(k*sqrt(x^2+b)) # = 0
c1*(x^2+b) * dy - (k+c1)*x * y +
	+ (c1^2*k*x*(x^2+b) + (k+c1)*x) * sin(k*sqrt(x^2+b)) # = 0

# D2 =>
c1^2*(x^2+b)^2 * d2y +
	+ c1*(c1-k)*x*(x^2+b) * dy +
	+ (c1^2*k^2*(x^2+b)*x^2 + (k^2-c1^2)*x^2 - b*c1*(k+c1)) * y +
	+ (3*c1^3*k*x^4 - c1^2*k^2*x^4 + c1^2*x^2 - k^2*x^2 + 4*c1^3*k*b*x^2 - c1^2*k^2*b*x^2 +
		+ c1^2*b + c1*k*b + c1^3*k*b^2) * sin(k*sqrt(x^2+b)) # = 0

### ODE:
x*(x^2+b)^2 * (c1^2*k*(x^2+b) + k+c1) * d2y +
	- (x^2+b) * (2*k*c1^2*x^4 + 3*b*k*c1^2*x^2 + b^2*k*c1^2 + b*(k+c1)) * dy +
	+ k*x^3 * (c1^2*k^2*x^4 + 2*c1^2*x^2 + 3*c1*k*x^2 + k^2*x^2 + 2*c1^2*k^2*b*x^2 +
		+ 2*c1^2*b + 3*c1*k*b + k^2*b + c1^2*k^2*b^2) * y # = 0

### Special Cases:

# Case: c1 = 0;
x*(x^2+b) * d2y - b*dy + k^2*x^3 * y # = 0

# Case: c1 = k = b = 1;
x*(x^2+1)*(x^2+3) * d2y +
	- (2*x^4 + 3*x^2 + 3) * dy + x^3 * (x^2+7) * y # = 0

# Case: c1 = -1; k = 1;
x*(x^2+b) * d2y - (2*x^2+b) * dy + x^3 * y # = 0

