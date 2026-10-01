#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Radicals: Power
##
## draft v.0.1e

### Log to Power:
# y = P(x)^n + P(x)^(j*n)


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### ODE: Order 1

### y = (x^2 + k)^n + (x^2 + k)^(2*n)

# Check:
k = 1/sqrt(5);
n = 1/sqrt(7);
x = sqrt(3); params = list(x=x, k=k);
e = expression((x^2 + k)^n + (x^2 + k)^(2*n))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2 + k)*dy - 2*n*x*((x^2 + k)^n + 2*(x^2 + k)^(2*n)) # = 0

# System:
R = (x^2 + k)^n; # =>
R^2 + R - y # = 0
2*n*x*(2*R^2 + R) - (x^2 + k)*dy # = 0
# =>
2*n*x*R + (x^2 + k)*dy - 4*n*x*y  # = 0
# =>
4*n^2*x^2*R^2 + 4*n^2*x^2*R - 4*n^2*x^2*y # = 0
((x^2 + k)*dy - 4*n*x*y)^2 +
	- 2*n*x*((x^2 + k)*dy - 4*n*x*y) - 4*n^2*x^2*y # = 0

### ODE:
(x^2+k)^2 * dy^2 - 8*n*x*(x^2+k) * y*dy + 16*n^2*x^2 * y^2 +
	- 2*n*x*((x^2+k)*dy - 2*n*x*y) # = 0


################

### Variation of Coefficients

### y = c1*(x^2 + k)^n + (x^2 + k)^(2*n)

# Check:
c1 = 2^(2/3);
k = 1/sqrt(5);
n = 1/sqrt(7);
x = sqrt(3); params = list(x=x, k=k, n=n, c1=c1);
e = expression(c1*(x^2 + k)^n + (x^2 + k)^(2*n))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2 + k)*dy - 2*n*x*(c1*(x^2 + k)^n + 2*(x^2 + k)^(2*n)) # = 0

# System:
R = (x^2 + k)^n; # =>
R^2 + c1*R - y # = 0
2*n*x*(2*R^2 + c1*R) - (x^2 + k)*dy # = 0
# =>
2*c1*n*x*R + (x^2 + k)*dy - 4*n*x*y # = 0

### ODE:
(x^2+k)^2 * dy^2 - 8*n*x*(x^2+k) * y*dy - 2*n*c1^2*x*(x^2+k) * dy +
	+ 16*n^2*x^2 * y^2 + 4*n^2*c1^2*x^2 * y # = 0

# Alternative formulation:
((x^2+k)*dy - 4*n*x*y)^2 - 2*n*c1^2*x*((x^2+k)*dy - 2*n*x*y) # = 0


##################
##################

### 3 Radicals ###

### y = (x^2 + k)^n + (x^2 + k)^(2*n) + (x^2 + k)^(3*n)
# Note: Order 1 will be presented separately;

# Check:
k = 1/sqrt(5);
n = 1/sqrt(7); # n = 1/5; # n = k = 1/9;
x = sqrt(3); params = list(x=x, k=k);
e = expression((x^2 + k)^n + (x^2 + k)^(2*n) + (x^2 + k)^(3*n))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2 + k)*dy - 2*n*x*((x^2 + k)^n + 2*(x^2 + k)^(2*n) + 3*(x^2 + k)^(3*n)) # = 0

### System
R = (x^2 + k)^n;
R^3 + R^2 + R - y # = 0
2*n*x*(3*R^3 + 2*R^2 + R) - (x^2 + k)*dy # = 0

### ODE:
(x^2+k)^3 * dy^3 +
	- 18*n*x*(x^2+k)^2 * y*dy^2 - 6*n*x*(x^2+k)^2 * dy^2 +
	+ 108*n^2*x^2*(x^2+k) * y^2*dy + 56*n^2*x^2*(x^2+k) * y*dy +
	+ 12*n^2*x^2*(x^2+k) * dy +
	- 216*n^3*x^3*y^3 - 112*n^3*x^3*y^2 - 24*n^3*x^3*y # = 0

# Alternative formulation:
((x^2+k)*dy - 6*n*x*y)^3 +
	- 6*n*x*((x^2+k)*dy - 6*n*x*y)^2 +
	- 16*n^2*x^2*(x^2+k) * y*dy +
	+ 12*n^2*x^2*(x^2+k) * dy +
	+ 104*n^3*x^3*y^2 - 24*n^3*x^3*y # = 0


#########################
#########################

### ODE: Order 2

### y = (x^2 + k)^n + (x^2 + k)^(2*n) + (x^2 + k)^(3*n)
# Note: Order 1 will be presented separately;

# Check:
k = 1/sqrt(5);
n = 1/sqrt(7); # n = 1/5; # n = k = 1/9;
x = sqrt(3); params = list(x=x, k=k);
e = expression((x^2 + k)^n + (x^2 + k)^(2*n) + (x^2 + k)^(3*n))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2 + k)*dy - 2*n*x*((x^2 + k)^n + 2*(x^2 + k)^(2*n) + 3*(x^2 + k)^(3*n)) # = 0

# D2 =>
(x^2 + k)*d2y + 2*x*dy +
	- 4*n^2*x^2*((x^2 + k)^n + 4*(x^2 + k)^(2*n) + 9*(x^2 + k)^(3*n)) / (x^2 + k) +
	- 2*n*((x^2 + k)^n + 2*(x^2 + k)^(2*n) + 3*(x^2 + k)^(3*n)) # = 0

### System:
# - there are distinct options regarding dependent radicals;
R = (x^2 + k)^n; R3 = (x^2 + k)^(3*n); # =>
R^2 + R - y + R3 # = 0
2*n*x*(2*R^2 + R + 3*R3) - (x^2+k)*dy # = 0
(x^2+k)^2*d2y + 2*x*(x^2+k)*dy +
	- 4*n^2*x^2*(4*R^2 + R + 9*R3) +
	- 2*n*(x^2+k)*(2*R^2 + R + 3*R3) # = 0


### ODE:
x^2*(x^2+k)^4 * d2y^2 +
	- 2*x*((10*n-1)*x^2 + k)*(x^2+k)^3 * dy*d2y +
	+ 48*n^2*x^4*(x^2+k)^2 * y*d2y  +  16*n^2*x^4*(x^2+k)^2 * d2y +
	+ (x^2+k)^2*((10*n-1)*x^2 + k)^2 * dy^2 +
	- 48*n^2*x^3*(x^2+k)*((10*n-1)*x^2 + k) * y*dy +
	- 16*n^2*x^3*(x^2+k)*((8*n-1)*x^2 + k) * dy +
	+ 576*n^4*x^6*y^2 + 192*n^4*x^6*y # = 0


### Special Cases:

### Case: n = k = 1/9;
x^2*(9*x^2+1)^4 * d2y^2 +
	- 2*x*(x^2+1)*(9*x^2+1)^3 * dy*d2y +
	+ 48*x^4*(9*x^2+1)^2 * y*d2y  +  16*x^4*(9*x^2+1)^2 * d2y +
	+ (x^2+1)^2*(9*x^2+1)^2 * dy^2 +
	- 48*x^3*(x^2+1)*(9*x^2+1) * y*dy +
	+ 16*x^3*(x^2-1)*(9*x^2+1) * dy +
	+ 576*x^6*y^2 + 192*x^6*y # = 0

### Case: n = 1/5;
x^2*(x^2+k)^4 * d2y^2 - 2*x*(x^2+k)^4 * dy*d2y +
	+ 48/5^2*x^4*(x^2+k)^2 * y*d2y  +  16/5^2*x^4*(x^2+k)^2 * d2y +
	+ (x^2+k)^4 * dy^2 - 48/5^2*x^3*(x^2+k)^2 * y*dy +
	- 16/5^3*x^3*(x^2+k)*(3*x^2 + 5*k) * dy +
	+ 576/5^4*x^6*y^2 + 192/5^4*x^6*y # = 0

