########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs - Coupled SQRT
##
## draft v.0.1i


### Examples:

# 2*x*(x-1)*(x^2-1) * d2y - x*(x-1)*(x+3) * dy + x*(x+2) * y - (x-1)^2 = 0;
# 9*x*(x^2-1)^2 * d2y - 3*(x^2-1) * (2*x^2 - 2*x + 3) * dy + (4*x^3 - 8*x^2 + 25*x - 3) * y = 0;



####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### Type: 1 Component

### Type: Atan * Exp
# y = B(x) * Sqrt(P1(x)) * Atan(K(x) * Sqrt(P1(x))) * exp(P2(x));
# Note: Homogeneous ODE;


### y = sqrt(x^2-1) * atan(sqrt(x^2-1)) * exp(k*x)
# - Simple example;

# Check:
k = sqrt(3); # k = 1;
x = 3^(3/5);
params = list(x=x, k=k);
e = expression(sqrt(x^2-1) * atan(sqrt(x^2-1)) * exp(k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x^2-1)*dy - k*x*(x^2-1)*y - x^2 * y - (x^2-1)*exp(k*x) # = 0

# D2 =>
x*(x^2-1) * d2y +
	- (k*x*(x^2-1) - 2*x^2 + 1) * dy +
	- (k*(3*x^2-1) + 2*x) * y +
	- k*(x^2-1)*exp(k*x) - 2*x*exp(k*x) # = 0
x*(x^2-1) * d2y +
	- (2*k*x*(x^2-1) - 2*x^2 + 1) * dy +
	+ (k^2*x^3 - 2*k*x^2 - (k^2+2)*x + k) * y - 2*x*exp(k*x) # = 0

### ODE:
x*(x^2-1)^2 * d2y +
	- (x^2-1)*(2*k*x*(x^2-1) + 1) * dy +
	+ (k^2*x^5 - 2*k^2*x^3 + k*x^2 + (k^2+2)*x - k) * y # = 0


### Special Cases:

### Case: k = 1;
x*(x^2-1)^2 * d2y +
	- (x^2-1)*(2*x*(x^2-1) + 1) * dy +
	+ (x^5 - 2*x^3 + x^2 + 3*x - 1) * y # = 0

##############

### Extension: w. coupled Exp(k*x)
# y = ... + c1*exp(k*x);

# Check:
k  = sqrt(3); # k = 1;
c1 = 1/5^(1/3); # c1 = 1; # c1 = 2;
# c1 = 1; k = -2;
x = 3^(3/5);
params = list(x=x, k=k, c1=c1);
e = expression(sqrt(x^2-1) * atan(sqrt(x^2-1)) * exp(k*x) + c1*exp(k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x^2-1)*dy - (k*x*(x^2-1) + x^2) * (y - c1*exp(k*x)) +
	- (x^2-1)*(c1*k*x + 1) * exp(k*x) # = 0
x*(x^2-1) * dy - (k*x*(x^2-1) + x^2) * y +
	+ ((c1-1)*x^2 + 1) * exp(k*x) # = 0

# D2 =>
x*(x^2-1) * d2y +
	- (k*x*(x^2-1) - 2*x^2 + 1) * dy +
	- (k*(3*x^2-1) + 2*x) * y +
	+ (k*(c1-1)*x^2 + k + 2*(c1-1)*x) * exp(k*x) # = 0
x*(x^2-1) * d2y +
	- (2*k*x*(x^2-1) - 2*x^2 + 1) * dy +
	+ (k^2*x*(x^2-1) - 2*k*x^2 - 2*x + k) * y +
	+ 2*(c1-1)*x * exp(k*x) # = 0
x^2*(x^2-1) * d2y +
	- x*(2*k*x*(x^2-1) - 1) * dy +
	+ x*(k^2*x*(x^2-1) - k) * y +
	- 2*exp(k*x) # = 0

### ODE:
x^2*(x^2-1)*((c1-1)*x^2 + 1) * d2y +
	- x*(2*(c1-1)*k*x^5 - 2*(c1-2)*k*x^3 +
		- (c1+1)*x^2 - 2*k*x + 1) * dy +
	+ x*(k^2*(c1-1)*x^5 - k^2*(c1-2)*x^3 +
		- k*(c1+1)*x^2 - (k^2+2)*x + k) * y # = 0


### Special Cases:

### Case: c1 = 1;
x*(x^2-1) * d2y +
	- (2*k*x*(x^2-1) - 2*x^2 + 1) * dy +
	+ (k^2*x*(x^2-1) - 2*k*x^2 - 2*x + k) * y # = 0
### Case: c1 = 1; k = -2;
x*(x^2-1) * d2y +
	+ (4*x*(x^2-1) + 2*x^2 - 1) * dy +
	+ 2*(x-1)*(2*x^2 + 4*x + 1) * y # = 0

### Case: c1 = 2;
x*(x^2-1)*(x^2 + 1) * d2y +
	- (2*k*x^5 - 3*x^2 - 2*k*x + 1) * dy +
	+ (k^2*x^5 - 3*k*x^2 - (k^2+2)*x + k) * y # = 0


###################
###################

### y = sqrt(x^2-1) * atan(sqrt(x^2-1)) * exp(k/x)

# Check:
k = sqrt(3); # k = 1;
x = 3^(3/5);
params = list(x=x, k=k);
e = expression(sqrt(x^2-1) * atan(sqrt(x^2-1)) * exp(k/x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x^2*(x^2-1) * dy - (x^3 - k*x^2 + k) * y - x*(x^2-1)*exp(k/x) # = 0

# D2 =>
x^3*(x^2-1) * d2y + x*(3*x^3 + k*x^2 - 2*x - k) * dy +
	- x*(3*x^2 - 2*k*x) * y +
	+ k*(x^2-1)*exp(k/x) - (3*x^3 - x)*exp(k/x) # = 0

### ODE:
x^4*(x^2-1)^2 * d2y + x^2*(x^2-1)*(2*k*x^2 - x - 2*k) * dy +
	- (2*k*x^5 - (k^2+2)*x^4 - 3*k*x^3 + 2*k^2*x^2 + k*x - k^2) * y # = 0

# Alternative notation:
x^4*(x^2-1)^2 * d2y + x^2*(x^2-1)*(2*k*x^2 - x - 2*k) * dy +
	- (k*x*(2*x^2-1)*(x^2-1) - k^2*(x^2-1)^2 - 2*x^4) * y # = 0


###################

### y = sqrt(x-1) * atan(sqrt(x-1)) * exp(k/x)

# Check:
k = sqrt(3); # k = 1;
x = 3^(3/5);
params = list(x=x, k=k);
e = expression(sqrt(x-1) * atan(sqrt(x-1)) * exp(k/x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*x^2*(x-1) * dy - (x^2 - 2*k*x + 2*k) * y - x*(x-1)*exp(k/x) # = 0

# D2 =>
2*x^3*(x-1) * d2y + x*(5*x^2 + 2*k*x - 4*x - 2*k) * dy +
	- x*(2*x - 2*k) * y +
	- (2*x^2 - (k+1)*x + k)*exp(k/x) # = 0

### ODE:
2*x^4*(x-1)^2 * d2y +
	+ x^2*(x-1) * (x^2 + (4*k-2)*x - 4*k) * dy +
	- ((3*k-1)*x^3 - k*(2*k+5)*x^2 + 2*k*(2*k+1)*x - 2*k^2) * y # = 0


##############
### Extension: Variant w. 2 Components

### y = sqrt(x-1) * atan(sqrt(x-1)) * exp(k/x) + B(x) * Exp(k/x)

# Check:
k  = sqrt(3); # k = 1;
c1 = - 1/4^(1/3); # c1 = 1; # c1 = 2; # c1 = -1;
x = 3^(3/5);
params = list(x=x, k=k);
e = expression(sqrt(x-1) * atan(sqrt(x-1)) * exp(k/x) + c1*exp(k/x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*x^2*(x-1) * dy - (x^2 - 2*k*x + 2*k) * (y - c1*exp(k/x)) +
	- (x*(x-1) - 2*c1*k*(x-1)) * exp(k/x) # = 0
2*x^2*(x-1) * dy - (x^2 - 2*k*x + 2*k) * y +
	+ ((c1-1)*x^2 + x) * exp(k/x) # = 0

# D2 =>
2*x^3*(x-1) * d2y + (5*x^3 + (2*k-4)*x^2 - 2*k*x) * dy +
	- x*(2*x - 2*k) * y +
	+ (2*(c1-1)*x^2 - k*(c1-1)*x + x - k) * exp(k/x) # = 0

### ODE:
2*x^4*(x-1)*((c1-1)*x + 1) * d2y +
	+ x^2 * ((c1-1)*x^3 + (4*(c1-1)*k + 3)*x^2 - (4*(c1-2)*k + 2)*x - 4*k) * dy +
	- ((3*k*(c1-1)+1)*x^3 - k*(2*k*(c1-1) + (4*c1 - 5))*x^2 +
		+ 2*k*(k*(c1-2) - 1)*x + 2*k^2) * y # = 0

### Special Cases:

### Case: c1 = 1;
2*x^4*(x-1) * d2y + x^2*(3*x^2 + (4*k-2)*x - 4*k) * dy +
	- (x-k)*(x^2 + 2*k*x - 2*k) * y # = 0

### Case: c1 = 2;
2*x^4*(x^2 - 1) * d2y +
	+ x^2 * (x^3 + (4*k + 3)*x^2 - 2*x - 4*k) * dy +
	- ((3*k+1)*x^3 - k*(2*k + 3)*x^2 - 2*k*x + 2*k^2) * y # = 0

### Case: c1 = -1;
2*x^4*(x-1)*(2*x-1) * d2y +
	+ x^2 * (2*x^3 + (8*k-3)*x^2 - (12*k-2)*x + 4*k) * dy +
	- ((6*k-1)*x^3 - k*(4*k+9)*x^2 + 2*k*(3*k+1)*x - 2*k^2) * y # = 0
# Alternative notation:
2*x^4*(x-1)*(2*x-1) * d2y +
	+ x^2 * ((x-1)*(2*x-1)*(x + 4*k) + x) * dy +
	- ((6*k-1)*x^3 - k*(4*k+9)*x^2 + 2*k*(3*k+1)*x - 2*k^2) * y # = 0


#########################
#########################

### Type: Atan * Radical
# y = B(x) * Sqrt(P1(x)) * Atan(K(x) * Sqrt(P1(x))) * P2(x)^p;
# Note: Homogeneous ODE;


### y = sqrt(x^2-1) * atan(sqrt(x^2-1)) * (x^3-1)^(1/3)
# - Simple example;

# Check:
x = 3^(3/5);
params = list(x=x);
e = expression(sqrt(x^2-1) * atan(sqrt(x^2-1)) * (x^3-1)^(1/3))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x^2-1)*(x^3-1)*dy - x^2*(x^3-1) * y - x^3*(x^2-1) * y +
	- (x^2-1)*(x^3-1) * (x^3-1)^(1/3) # = 0
x*(x^2-1)*(x^3-1)*dy - (2*x^5 - x^3 - x^2) * y +
	- (x^2-1)*(x^3-1) * (x^3-1)^(1/3) # = 0
x*(x+1)*(x^3-1)*dy - (2*x^4 + 2*x^3 + x^2) * y +
	- (x+1)*(x^3-1) * (x^3-1)^(1/3) # = 0

# D2 =>
x*(x^2-1)*(x^3-1) * d2y + (4*x^5 - 3*x^3 - 2*x^2 + 1) * dy +
	- (10*x^4 - 3*x^2 - 2*x) * y +
	- (x-1)*(6*x^3 + 6*x^2 + 2*x) * (x^3-1)^(1/3) # = 0

### ODE:
x*(x+1)^2*(x^3-1)^2 * d2y +
	- (x+1)*(x^3-1) * (2*x^4 + 2*x^3 + x^2 + x + 1) * dy +
	+ x*(2*x^6 + 4*x^5 + 5*x^4 + 8*x^3 + 12*x^2 + 7*x + 2) * y # = 0


#########################

### y = sqrt(x^2-1) * atan(sqrt(x^2-1)) * (x+1)^(1/3)
# - Slightly simpler example;

# Check:
x = 3^(3/5);
params = list(x=x);
e = expression(sqrt(x^2-1) * atan(sqrt(x^2-1)) * (x+1)^(1/3))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
3*x*(x+1)*(x^2-1) * dy - 3*x^2*(x+1) * y - x*(x^2-1) * y +
	- 3*(x+1)*(x^2-1) * (x+1)^(1/3) # = 0
3*x*(x^2-1) * dy - x*(4*x-1) * y +
	- 3*(x^2-1) * (x+1)^(1/3) # = 0

# D2 =>
3*x*(x^2-1) * d2y + (5*x^2 + x - 3) * dy - (8*x-1) * y +
	- (7*x-1) * (x+1)^(1/3) # = 0

### ODE:
9*x*(x^2-1)^2 * d2y - 3*(x^2-1) * (2*x^2 - 2*x + 3) * dy +
	+ (4*x^3 - 8*x^2 + 25*x - 3) * y # = 0


#########################

### y = sqrt(x^2-1) * atan(sqrt(x^2-1)) * (3*x+1)^(1/3)
# - Slightly simpler example;

# Check:
x = 3^(3/5);
params = list(x=x);
e = expression(sqrt(x^2-1) * atan(sqrt(x^2-1)) * (3*x+1)^(1/3))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(3*x+1)*(x^2-1) * dy - (4*x^3 + x^2 - x) * y +
	- (3*x+1)*(x^2-1) * (3*x+1)^(1/3) # = 0

# D2 =>
x*(3*x+1)*(x^2-1) * d2y + (8*x^3 + 2*x^2 - 5*x - 1) * dy +
	- (12*x^2 + 2*x - 1) * y +
	- (10*x^2 + 2*x - 4) * (3*x+1)^(1/3) # = 0

### ODE:
x*(3*x+1)^2*(x^2-1)^2 * d2y +
	- (3*x+1)*(x^2-1) * (2*x^3 + x + 1) * dy +
	+ (4*x^5 + 13*x^3 + 13*x^2 + 3*x - 1) * y # = 0


#########################

### y = sqrt(x^2-1) * atan(sqrt(x^2-1)) * (x^2+1)^(1/2)
# - Slightly simpler example;

# Check:
x = 3^(3/5);
params = list(x=x);
e = expression(sqrt(x^2-1) * atan(sqrt(x^2-1)) * (x^2+1)^(1/2))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x^2+1)*(x^2-1) * dy - x^2*(x^2+1) * y - x^2*(x^2-1) * y +
	- (x^2+1)*(x^2-1) * (x^2+1)^(1/2) # = 0
x*(x^4-1) * dy - 2*x^4 * y - (x^4-1) * (x^2+1)^(1/2) # = 0

# D2 =>
x*(x^4-1) * d2y + (3*x^4-1) * dy - 8*x^3 * y +
	- x*(5*x^2-1) * (x^2+1)^(1/2) # = 0

### ODE:
x*(x^4-1)^2 * d2y - (x^4-1)*(2*x^4 - x^2 + 1) * dy +
	+ 2*x^3*(x^4 - x^2 + 4) * y # = 0


#########################
#########################

### Type: 2 Components

### y = B1(x) * SQRT(P1(x)) * ATAN(SQRT(P1(x))) + B2(x) * SQRT(P2(x)) * ATAN(SQRT(P2(x)))
# Note: Non-Homogenous;

### y = sqrt(x^2-1) * atan(sqrt(x^2-1)) + c1*sqrt(x-1) * atan(sqrt(x-1))
# - Simple example;

# Check:
c1 = sqrt(2); # c1 = -1; # c1 = -1/2;
x = 3^(3/5);
params = list(x=x, c1=c1);
e = expression(sqrt(x^2-1) * atan(sqrt(x^2-1)) + 2*c1*sqrt(x-1) * atan(sqrt(x-1)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x^2-1)*dy +
	- x^2 * sqrt(x^2-1) * atan(sqrt(x^2-1)) +
	- c1*x*(x+1) * sqrt(x-1) * atan(sqrt(x-1)) +
	- (c1+1)*(x^2-1) # = 0

# Linear System:
div = x*(x-1); # 2*x^2 - x*(x+1);
sqrt(x^2-1) * atan(sqrt(x^2-1)) # ==
(2*x*(x^2-1)*dy - 2*(c1+1)*(x^2-1) - x*(x+1)*y) / div;
#
c1*sqrt(x-1) * atan(sqrt(x-1)) # ==
- (x*(x^2-1)*dy - (c1+1)*(x^2-1) - x^2*y) / div;

# D2 =>
2*x*(x^2-1)*d2y + 2*(3*x^2-1)*dy +
	- 2*x*(3*x^2-2) * sqrt(x^2-1) * atan(sqrt(x^2-1)) / (x^2-1) +
	- c1*(5*x^2 - x - 2) * sqrt(x-1) * atan(sqrt(x-1)) / (x-1) +
	- c1*x - 6*x - 4*c1*x - c1 # = 0

### ODE:
2*x^2*(x-1)*(x^2-1) * d2y - x^2*(x-1)*(x+3) * dy + x^2*(x+2) * y +
	+ (x-1)*((2*c1+1)*x^2 + (2*c1+3)*x - 2*(c1+1)) # = 0


### Special Cases:

### Case: c1 = -1;
2*x*(x-1)*(x^2-1) * d2y - x*(x-1)*(x+3) * dy + x*(x+2) * y - (x-1)^2 # = 0

### Case: c1 = -1/2;
2*x^2*(x-1)*(x^2-1) * d2y - x^2*(x-1)*(x+3) * dy + x^2*(x+2) * y +
	+ (x-1)*(2*x - 1) # = 0


#####################

### Example 2:
### y = sqrt(x^2-1) * atan(sqrt(x^2-1)) + c1*sqrt(x^4-1) * atan(sqrt(x^4-1))

# Check:
c1 = sqrt(2); # c1 = -1; # c1 = -1/2;
x = 3^(3/5);
params = list(x=x, c1=c1);
e = expression(sqrt(x^2-1) * atan(sqrt(x^2-1)) + c1*sqrt(x^4-1) * atan(sqrt(x^4-1)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x^2-1)*(x^4-1)*dy +
	- x^2*(x^4-1) * sqrt(x^2-1) * atan(sqrt(x^2-1)) +
	- 2*c1*x^4*(x^2-1) * sqrt(x^4-1) * atan(sqrt(x^4-1)) +
	- (2*c1+1)*(x^2-1)*(x^4-1) # = 0

# Linear System:
div  = 2*x^4*(x^2-1) - x^2*(x^4-1);
div1 = x^2*(x^2-1);
sqrt(x^2-1) * atan(sqrt(x^2-1)) # ==
- (x*(x^2-1)*(x^4-1)*dy - 2*x^4*(x^2-1)*y - (2*c1+1)*(x^2-1)*(x^4-1)) / div;
- (x*(x^4-1)*dy - 2*x^4*y - (2*c1+1)*(x^4-1)) / div1;
#
c1*sqrt(x^4-1) * atan(sqrt(x^4-1)) # ==
(x*(x^2-1)*(x^4-1)*dy - x^2*(x^4-1)*y - (2*c1+1)*(x^2-1)*(x^4-1)) / div;
(x*(x^4-1)*dy - x^2*(x^2+1)*y - (2*c1+1)*(x^4-1)) / div1;

# D2 =>
x*(x^2-1)*(x^4-1) * d2y + (7*x^6 - 5*x^4 - 3*x^2 + 1) * dy +
	- x^3*(x^4-1) * sqrt(x^2-1) * atan(sqrt(x^2-1)) / (x^2-1) +
	- (6*x^5-2*x) * sqrt(x^2-1) * atan(sqrt(x^2-1)) +
	- 4*c1*x^7*(x^2-1) * sqrt(x^4-1) * atan(sqrt(x^4-1)) / (x^4-1) +
	- 2*c1*(6*x^5-4*x^3) * sqrt(x^4-1) * atan(sqrt(x^4-1)) +
	- x*(x^4-1) - 4*c1*x^3*(x^2-1) - (2*c1+1)*(6*x^5 - 4*x^3 - 2*x) # = 0
x*(x^2-1)*(x^4-1) * d2y + (7*x^6 - 5*x^4 - 3*x^2 + 1) * dy +
	- (7*x^5+x^3-2*x) * sqrt(x^2-1) * atan(sqrt(x^2-1)) +
	- 4*c1*x^3 * (4*x^4 + x^2 - 2) * sqrt(x^4-1) * atan(sqrt(x^4-1)) / (x^2+1) +
	- (7*x^5 - 4*x^3 - 3*x) - 4*c1*(4*x^5 - 3*x^3 - x) # = 0
x*(x^2-1)*(x^4-1) * d2y + (7*x^6 - 5*x^4 - 3*x^2 + 1) * dy +
	+ (7*x^5+x^3-2*x) * (x*(x^4-1)*dy - 2*x^4*y - (2*c1+1)*(x^4-1)) / div1 +
	- 4*x^3 * (4*x^4 + x^2 - 2) * (x*(x^2-1)*dy - x^2*y - (2*c1+1)*(x^2-1)) / div1 +
	- (7*x^5 - 4*x^3 - 3*x) - 4*c1*(4*x^5 - 3*x^3 - x) # = 0

### ODE:
x^2*(x^2-1)*(x^4-1) * d2y +
	- x*(x^2-1)*(2*x^4 + 3*x^2 - 1) * dy +
	+ 2*x^4 * (x^2 + 2) * y +
	+ 2*c1*(x^6 + 2*x^4 - 5*x^2 + 2) + 2*(x^6 - 2*x^2 + 1) # = 0

### Special Cases:

### Case: c1 = -1;
x^2*(x^2-1)*(x^4-1) * d2y +
	- x*(x^2-1)*(2*x^4 + 3*x^2 - 1) * dy +
	+ 2*x^4 * (x^2 + 2) * y  - 2*(x^2-1)*(2*x^2-1) # = 0

### Case: c1 = -1/2;
x*(x^2-1)*(x^4-1) * d2y +
	- (x^2-1)*(2*x^4 + 3*x^2 - 1) * dy +
	+ 2*x^3 * (x^2 + 2) * y + x^5 - 2*x^3 + x # = 0


####################

### Mixed Components

