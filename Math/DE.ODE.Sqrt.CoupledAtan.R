########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs - Coupled SQRT
##
## draft v.0.1c


### Examples:

# 2*x*(x-1)*(x^2-1) * d2y - x*(x-1)*(x+3) * dy + x*(x+2) * y - (x-1)^2 = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### Type: 1 Component



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

