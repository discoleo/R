#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODE: Coupled SQRT * Log(SQRT)
##
## draft v.0.1b


### Coupled SQRT * Log(SQRT)
# y = SQRT(P(x)^2 + b) * Log(SQRT(P(x)^2 + b) - P(x))


### Examples:

# 4*(x^2+b)^2 * d2y + (x^2 - 2*b) * y = 0;
# (x^2-1)^2 * d2y - (x^2-1)*(2*x^2 + x - 2) * dy + (x^4 + x^3 - x^2 - x + 2) * y = 0;
# (x^2-1)*(x^2 - x - 1) * d2y - (2*x^4 - x^3 - 4*x^2 + x + 3) * dy + (x^4 - x^2 + 3) * y = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### Mixed w. Exp

### y = sqrt(x^2 + b) * log(sqrt(x^2 + b) - x) * exp(k*x)

# Check:
k = 1/sqrt(5); b = exp(1/3); # b = -1 + 0i; k = 1;
x = sqrt(3); params = list(x=x, k=k, b=b);
e = expression(sqrt(x^2 + b) * log(sqrt(x^2 + b) - x) * exp(k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2 + b) * dy - k*(x^2 + b) * y - x*y + (x^2 + b)*exp(k*x) # = 0

# D2 =>
(x^2 + b) * d2y - (k*x^2 - x + b*k) * dy - (2*k*x + 1) * y +
	+ k*(x^2 + b)*exp(k*x) + 2*x*exp(k*x) # = 0
(x^2 + b) * d2y - (2*k*x^2 - x + 2*b*k) * dy +
	+ (k^2*x^2 - k*x + b*k^2 - 1) * y + 2*x*exp(k*x) # = 0

### ODE:
(x^2+b)^2 * d2y - (x^2+b)*(2*k*x^2 + x + 2*b*k) * dy +
	+ ((x^2+b)*(k^2*x^2 + k*x + b*k^2 + 1) - 2*b) * y # = 0

### Special Cases:

### Case: b = -1; k = 1;
(x^2-1)^2 * d2y - (x^2-1)*(2*x^2 + x - 2) * dy +
	+ (x^4 + x^3 - x^2 - x + 2) * y # = 0


#########################

### Mixed: Sum w. Exp

### y = ... + C(x) * exp(k*x)

# Check:
c1 = exp(-1); # c1 = -1;
k = 1/sqrt(5); b = exp(1/3);
# c1 = -1; k = 1; # c1 = -1; k = -1;
# c1 = -1; k = 1; b = 1; # c1 = -1; k = 1; b = -1 + 0i;
x = sqrt(3); params = list(x=x, k=k, b=b);
e = expression(sqrt(x^2 + b) * log(sqrt(x^2 + b) - x) * exp(k*x) + c1*exp(k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2 + b) * dy - (k*x^2 + x + b*k) * (y - c1*exp(k*x)) +
	- (c1*k - 1)*(x^2 + b)*exp(k*x) # = 0
(x^2 + b) * dy - (k*x^2 + x + b*k) * y +
	+ (x^2 + c1*x + b)*exp(k*x) # = 0

# D2 =>
(x^2 + b) * d2y - (k*x^2 - x + b*k) * dy +
	- (2*k*x + 1) * y +
	+ (k*x^2 + (c1*k+2)*x + b*k + c1)*exp(k*x) # = 0
(x^2 + b) * d2y - (2*k*x^2 - x + 2*b*k) * dy +
	+ (k^2*x^2 - k*x + b*k^2 - 1) * y +
	+ (2*x + c1)*exp(k*x) # = 0

### ODE:
(x^2 + b)*(x^2 + c1*x + b) * d2y +
	- (2*k*x^4 + (2*c1*k+1)*x^3 + 4*b*k*x^2 + b*(2*c1*k+1)*x + b*c1 + 2*b^2*k) * dy +
	+ (k^2*x^4 + k*(c1*k+1)*x^3 + (2*k^2*b+1)*x^2 + k*b*(k*c1+1)*x + k^2*b^2 + k*b*c1 - b) * y # = 0


### Special Cases:

### Case: c1 = -1;
(x^2 + b)*(x^2 - x + b) * d2y +
	- (2*k*x^4 - (2*k-1)*x^3 + 4*b*k*x^2 - b*(2*k-1)*x - b + 2*b^2*k) * dy +
	+ (k^2*x^4 - k*(k-1)*x^3 + (2*k^2*b+1)*x^2 - k*b*(k-1)*x + k^2*b^2 - k*b - b) * y # = 0
# Case: c1 = -1; k = 1;
(x^2 + b)*(x^2 - x + b) * d2y +
	- (2*x^4 - x^3 + 4*b*x^2 - b*x + 2*b^2 - b) * dy +
	+ (x^4 + (2*b+1)*x^2 + b^2 - 2*b) * y # = 0
# Case: c1 = -1; k = -1;
(x^2+b)*(x^2 - x + b) * d2y +
	+ ((x^2+b)*(2*x^2 - 3*x + 2*b) + b) * dy +
	+ (x^4 - 2*x^3 + (2*b+1)*x^2 - 2*b*x + b^2) * y # = 0
# Case: c1 = -1; k = 1; b = 1;
(x^2 + 1)*(x^2 - x + 1) * d2y - (2*x^4 - x^3 + 4*x^2 - x + 1) * dy +
	+ (x^4 + 3*x^2 - 1) * y # = 0
# Case: c1 = -1; k = 1; b = -1;
(x^2 - 1)*(x^2 - x - 1) * d2y +
	- (2*x^4 - x^3 - 4*x^2 + x + 3) * dy +
	+ (x^4 - x^2 + 3) * y # = 0


############################

### Mixed w. Another Radical

### y = sqrt(x^2 + b) * log(sqrt(x^2 + b) - x) * P2(x)^r
# - for simplicity: P2(x) = (x^2+b);

# Check:
r = exp(-1); # r = -1/4;
b = exp(1/3);
x = sqrt(3); params = list(x=x, b=b, r=r);
e = expression(sqrt(x^2 + b) * log(sqrt(x^2 + b) - x) * (x^2+b)^r)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2+b)*dy - (2*r+1)*x*y + (x^2+b) * (x^2+b)^r # = 0

# D2 =>
(x^2+b)*d2y - (2*r-1)*x*dy - (2*r+1)*y + 2*(r+1)*x * (x^2+b)^r # = 0

### ODE:
(x^2+b)^2 * d2y - (4*r+1)*x*(x^2+b) * dy + (2*r+1)  * ((2*r+1)*x^2 - b) * y # = 0


### Special Cases:

### Case: r = -1/4;
4*(x^2+b)^2 * d2y + (x^2 - 2*b) * y # = 0

