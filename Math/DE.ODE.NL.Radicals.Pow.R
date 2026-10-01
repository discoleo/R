#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Radicals: Power
##
## draft v.0.1b

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


#########################
#########################

### ODE: Order 2

### y = (x^2 + k)^n + (x^2 + k)^(2*n) + (x^2 + k)^(3*n)
# Note: Order 1 will be presented separately;

# Check:
k = 1/sqrt(5);
n = 1/sqrt(7);
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

x^10*d2y^2 + 4*k*x^8*d2y^2 + 6*k^2*x^6*d2y^2 + 4*k^3*x^4*d2y^2 + k^4*x^2*d2y^2 +
	+ 2*x^9*dy*d2y - 20*n*x^9*dy*d2y + 4*k*x^7*dy*d2y - 60*k*n*x^7*dy*d2y - 60*k^2*n*x^5*dy*d2y +
	- 4*k^3*x^3*dy*d2y - 20*k^3*n*x^3*dy*d2y - 2*k^4*x*dy*d2y +
	+ 48*n^2*x^8*y*d2y + 96*k*n^2*x^6*y*d2y + 48*k^2*n^2*x^4*y*d2y +
	+ 16*n^2*x^8*d2y + 32*k*n^2*x^6*d2y + 16*k^2*n^2*x^4*d2y +
	+ x^8*dy^2 - 20*n*x^8*dy^2 + 100*n^2*x^8*dy^2 - 20*k*n*x^6*dy^2 + 200*k*n^2*x^6*dy^2 +
	- 2*k^2*x^4*dy^2 + 20*k^2*n*x^4*dy^2 + 100*k^2*n^2*x^4*dy^2 + 20*k^3*n*x^2*dy^2 +
	+ k^4*dy^2 + 48*n^2*x^7*y*dy - 480*n^3*x^7*y*dy - 480*k*n^3*x^5*y*dy - 48*k^2*n^2*x^3*y*dy +
	+ 16*n^2*x^7*dy - 128*n^3*x^7*dy - 128*k*n^3*x^5*dy - 16*k^2*n^2*x^3*dy +
	+ 576*n^4*x^6*y^2 + 192*n^4*x^6*y # = 0

# TODO: format;

