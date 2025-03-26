########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs - Exponential
## "Linear" Combinations
##
## draft v.0.1a


### History:

### draft v.0.1a:
# - moved to this file from:
#   DE.ODE.Gaussian.R (Section B.1: Mixed)


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### Linear Combinations

### y = a2*e^(-x^2) + a3*e^(-x^3)
# [not run]
dy  = -2*a2*x*exp(-x^2) - 3*a3*x^2*exp(-x^3)
d2y = 4*a2*x^2*exp(-x^2) - 2*a2*exp(-x^2) + 9*a3*x^4*exp(-x^3) - 6*a3*x*exp(-x^3)
# Quasi-Linear System:
# e^(-x^2) = (dy + 3*x^2*y) / (a2*(3*x^2 - 2*x))
# e^(-x^3) = - (dy + 2*x*y) / (a3*(3*x^2 - 2*x))
# =>

### ODE:
(3*x^2 - 2*x)*d2y + (9*x^4 - 4*x^2 - 6*x + 2)*dy + 6*x^2*(3*x^3 - 2*x^2 - 1)*y # = 0


### Plot:
y = function(x, a=c(1, 1)) {
	y = a[1] * exp(-x^2) + a[2] * exp(-x^3)
	y = round0(y)
	return(y)
}
dy = function(x, a=c(1, 1)) {
	y.x = y(x, a=a)
	dp = - (2*a[1]*x*exp(-x^2) + 3*a[2]*x^2*exp(-x^3))
	return(dp)
}
d2y = function(x, a=c(1, 1)) {
	y.x = y(x, a=a)
	dy.x = dy(x, a=a)
	dp = (9*x^4 - 4*x^2 - 6*x + 2)*dy.x + 6*x^2*(3*x^3 - 2*x^2 - 1)*y.x;
	dp = - dp;
	div = 3*x^2 - 2*x;
	dp = ifelse(div != 0, dp / div, -Inf); # TODO
	return(dp)
}
### Plot:
a = c(1, 1) # has NO effect on eq of D2;
curve(y(x, a=a), from= -1, to = 3, ylim=c(-2, 3))
line.tan(c(-3:2 * 2/5, 3/2, 2), dx=3, p=y, dp=dy, a=a)
# non-sigmoidal
curve(dy(x, a=a), add=T, col="green")
line.tan(c(-3:2 * 2/5, 3/2, 2), dx=3, p=dy, dp=d2y, a=a, col="orange")

###
a = c(1, -1) # a = c(2, -1/2)
# although has NO effect on eq. of D2, d2y depends indirectly;
curve(y(x, a=a), from= -1, to = 3, ylim=c(-2, 3))
line.tan(c(-3:2 * 2/5, 3/2, 2), dx=3, p=y, dp=dy, a=a)
# non-sigmoidal
curve(dy(x, a=a), add=T, col="green")
px = c(0.05 + -3:3 /5, 3/2, 2);
line.tan(px, dx=3, p=dy, dp=d2y, a=a, col="orange")
points(px, dy(px, a=a), col = "#FF3232B2")


######################
### Generalization ###
### (part)         ###

### y = a1*e^(-(x^n + b1*x)) + a2*e^(-(x^n + b2*x))

# Check:
z = sqrt(c(2,3,5,7)); n = 1/3; x = 1/sqrt(5);
a1 = z[1]; a2 = z[2]; b1 = z[3]; b2 = z[4];
#
params = list(x=x, a1=a1, a2=a2, b1=b1, b2=b2, n=n);
e   = expression(a1*exp(-(x^n + b1*x)) + a2*exp(-(x^n + b2*x)))[[1]];
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

dy = - a1*(n*x^(n-1) + b1)*exp(-(x^n + b1*x)) +
     - a2*(n*x^(n-1) + b2)*exp(-(x^n + b2*x))
d2y = a1*(n*x^(n-1) + b1)^2 * exp(-(x^n + b1*x)) +
      a2*(n*x^(n-1) + b2)^2 * exp(-(x^n + b2*x)) +
	- a1*n*(n-1)*x^(n-2) * exp(-(x^n + b1*x)) +
	- a2*n*(n-1)*x^(n-2) * exp(-(x^n + b2*x));
# e^(-(x^n + b1*x)) =   (dy + (n*x^(n-1) + b2)*y) / (a1*(b2 - b1))
# e^(-(x^n + b2*x)) = - (dy + (n*x^(n-1) + b1)*y) / (a2*(b2 - b1))
# =>
(b2 - b1)*d2y - (n*x^(n-1) + b1)^2 * (dy + (n*x^(n-1) + b2)*y) +
  + (n*x^(n-1) + b2)^2 * (dy + (n*x^(n-1) + b1)*y) +
  + n*(n-1)*x^(n-2) * (dy + (n*x^(n-1) + b2)*y) +
  - n*(n-1)*x^(n-2) * (dy + (n*x^(n-1) + b1)*y) # = 0
(b2 - b1)*d2y - (n*x^(n-1) + b1)^2 * (dy + (n*x^(n-1) + b2)*y) +
  + (n*x^(n-1) + b2)^2 * (dy + (n*x^(n-1) + b1)*y) +
  - n*(n-1)*(b1 - b2)*x^(n-2)*y # = 0

### ODE:
d2y + (2*n*x^(n-1) + b1 + b2)*dy +
  + (n^2*x^(2*n-2) + n*(b1 + b2)*x^(n-1) + n*(n-1)*x^(n-2) + b1*b2)*y # = 0

### Examples:
### n = 2
d2y + (4*x + b1 + b2)*dy + (4*x^2 + 2*(b1 + b2)*x + b1*b2 + 2)*y # = 0
### b1 + b2 = 0
d2y + 2*n*x^(n-1) * dy + (n^2*x^(2*n-2) + n*(n-1)*x^(n-2) - b1^2)*y # = 0
### n = 2; b1 = - b2 = sqrt(2);
d2y + 4*x*dy + 4*x^2*y # = 0


### Solution:
y = function(x, a=c(1, 1), b=c(-1, 1), n=2) {
	y = sapply(x, function(x) sum(a * exp(-x^n - b*x)))
	y = round0(y)
	return(y)
}
dy = function(x, a=c(1, 1), b=c(-1, 1), n=2) {
	# y.x = y(x, a=a, b=b, n=n)
	nxn = n*x^(n-1);
	dp = sapply(seq_along(x), function(id) {
		xx = x[id]
		x.exp = -a * exp(-xx^n - b*xx);
		sum((nxn[id] + b) * x.exp)
	})
	return(dp)
}
d2y = function(x, a=c(1, 1), b=c(-1, 1), n=2) {
	y.x = y(x, a=a, b=b, n=n)
	dy.x = dy(x, a=a, b=b, n=n)
	nxn1 = n*x^(n-2);
	nxn = x * nxn1; nxsq = nxn^2;
	dp = (2*nxn + b[1] + b[2])*dy.x +
		+ (nxsq + (b[1] + b[2])*nxn + (n-1)*nxn1 + b[1]*b[2])*y.x
	div = -1;
	dp = ifelse(div != 0, dp / div, -Inf); # TODO
	return(dp)
}
### Plot:
n = 2; a = c(1, 1) # a[] has NO effect on eq of D2;
curve(y(x, a=a, n=n), from= -2, to = 3, ylim=c(-2, 3))
sapply(c(-3:2 * 4/7, 2), line.tan, dx=3, p=y, dp=dy, a=a, n=n)
# wavelet
curve(dy(x, a=a, n=n), add=T, col="green")
sapply(c(-3:2 * 4/7, 3/2, 2), line.tan, dx=3, p=dy, dp=d2y, a=a, n=n, col="orange")

### n = 3
n = 3; a = c(1, 1) # a[] has NO effect on eq of D2;
curve(y(x, a=a, n=n), from= -1, to = 3, ylim=c(-4, 4))
sapply(c(-3:0 * 1/7, 0.6, 1, 3/2), line.tan, dx=3, p=y, dp=dy, a=a, n=n)
# wavelet
curve(dy(x, a=a, n=n), add=T, col="green")
sapply(c(-3:3 * 1/7, 1, 3/2, 2), line.tan, dx=3, p=dy, dp=d2y, a=a, n=n, col="orange")
### TODO: check d2y(1)!


#################

### Simple Variant:
### y = a1*e^(x^n) + a2*e^(-x^n)

# Check:
x  = 1/sqrt(3); n = sqrt(2); a1 = 5/4; a2 = -1/5;
y  = a1*exp(x^n) + a2*exp(-x^n);
# D =>
dy = a1*n*x^(n-1) * exp(x^n) - a2*n*x^(n-1) * exp(-x^n)
### D2 =>
d2y = a1*n^2*x^(2*n-2) * exp(x^n) + a1*n*(n-1)*x^(n-2) * exp(x^n) +
	+ a2*n^2*x^(2*n-2) * exp(-x^n) - a2*n*(n-1)*x^(n-2) * exp(-x^n)
### Linear system:
exp(x^n)  - ( dy + n*x^(n-1)*y) / (2*a1*n*x^(n-1)) # = 0
exp(-x^n) - (-dy + n*x^(n-1)*y) / (2*a2*n*x^(n-1)) # = 0
# =>
2*x*d2y # =
(n*x^n + (n-1))*(dy + n*x^(n-1)*y) - (n*x^n - (n-1))*(dy - n*x^(n-1)*y)

### ODE:
x*d2y - (n-1)*dy - n^2*x^(2*n-1)*y # = 0

### Solution:
y = function(x, n=2, a=c(1, 1), b=c(1, -1)) {
	y = sapply(x, function(x) sum(a * exp(b*x^n)))
	y = round0(y)
	return(y)
}
dy = function(x, n=2, a=c(1, 1), b=c(1, -1)) {
	# y.x = y(x, a=a, b=b, n=n)
	a[2] = -a[2];
	nxn = n*x^(n-1);
	dp = sapply(seq_along(x), function(id) {
		xx = x[id]
		sum(a * nxn[id] * exp(b*xx^n));
	})
	return(dp)
}
d2y = function(x, n=2, a=c(1, 1), b=c(1, -1)) {
	y.x = y(x, n=n, a=a, b=b)
	dy.x = dy(x, n=n, a=a, b=b)
	dp = (n-1)*dy.x + n^2*x^(2*n-1)*y.x
	div = x;
	dp = ifelse(div != 0, dp / div, 0); # TODO
	return(dp)
}
### Plot:
n = 2; a = c(1, 1) # a[] has NO effect on eq of D2;
curve(y(x, a=a, n=n), from= -2, to = 2, ylim=c(-2, 5))
sapply(c(-3:3 * 4/7, 2), line.tan, dx=3, p=y, dp=dy, a=a, n=n)
#
curve(dy(x, a=a, n=n), add=T, col="green")
sapply(c(-3:3 * 2/7), line.tan, dx=3, p=dy, dp=d2y, a=a, n=n, col="orange")


###########################
###########################

### Partial Extension
### y = x^2*e^(x^3) + x^3*e^(x^2) + F0(x)

### D(y)
(3*x^4 + 2*x)*e^(x^3) + (2*x^4 + 3*x^2)*e^(x^2) + df0

### Solve liniar system:
div = (2*x^4 + 3*x^2)*x^2 - x^3*(3*x^4 + 2*x)
div = -3*x^7 + 2*x^6 + x^4;
e^(x^3) = ((2*x^4 + 3*x^2)*(y - f0) - x^3*(dy - df0)) / div;
e^(x^2) = ((3*x^4 + 2*x)*(y - f0) - x^2*(dy - df0)) / - div;

### D2(y)
(3*x^2*(3*x^4 + 2*x) + 12*x^3 + 2)*e^(x^3) +
	+ (2*x*(2*x^4 + 3*x^2) + 8*x^3 + 6*x)*e^(x^2) + d2f0
(9*x^6 + 18*x^3 + 2)*e^(x^3) +
	+ (4*x^5 + 14*x^3 + 6*x)*e^(x^2) + d2f0
### div * d2y =
(9*x^6 + 18*x^3 + 2)*((2*x^4 + 3*x^2)*(y - f0) - x^3*(dy - df0)) +
	- (4*x^5 + 14*x^3 + 6*x)*((3*x^4 + 2*x)*(y - f0) - x^2*(dy - df0)) + div*d2f0

### Examples:
### f = x^3
div = -3*x^7 + 2*x^6 + x^4;
### div * d2y =
(9*x^6 + 18*x^3 + 2)*((2*x^4 + 3*x^2)*y - x^3*dy) +
	- (4*x^5 + 14*x^3 + 6*x)*((3*x^4 + 2*x)*y - x^2*dy) +
	- 2*x^7*(9*x^6 + 18*x^3 + 2) + x^3*(3*x^4 - x)*(4*x^5 + 14*x^3 + 6*x) + 6*x*(-3*x^7 + 2*x^6 + x^4)
(4*x^3 + 14*x^5 - 18*x^6 + 4*x^7 - 9*x^9)*dy +
	+ (-6*x^2 - 24*x^4 + 36*x^5 - 8*x^6 - 6*x^7 + 27*x^8 - 12*x^9 + 18*x^10)*y +
	- 6*x^7 - 4*x^9 + 6*x^10 + 12*x^12 - 18*x^13


### Solution & Plot:
y = function(x, a=c(1, 1)) {
	x2 = x^2; x3 = x^3;
	y = a[1]*x2*exp(x3) +
		a[2]*x3*exp(x2) + x3;
	return(y)
}
dy = function(x, a=c(1, 1)) {
	x2 = x^2; x4 = x2*x2;
	dy.v = a[1]*(3*x4 + 2*x)*exp(x^3) +
		a[2]*(2*x4 + 3*x2)*exp(x2) + 3*x2;
	return(dy.v)
}
d2y = function(x, a=c(1, 1)) {
	y.x = y(x, a=a)
	dy.x = dy(x, a=a)
	dp = (4*x^3 + 14*x^5 - 18*x^6 + 4*x^7 - 9*x^9)*dy.x +
		+ (-6*x^2 - 24*x^4 + 36*x^5 - 8*x^6 - 6*x^7 + 27*x^8 - 12*x^9 + 18*x^10)*y.x +
		- 6*x^7 - 4*x^9 + 6*x^10 + 12*x^12 - 18*x^13;
	div = -3*x^7 + 2*x^6 + x^4;
	dp = ifelse(div != 0, dp / div, d2y(x + 1E-5, a=a)); # TODO
	return(dp)
}
### Plot:
a = c(1, 1) # a[] has NO effect on eq of D2;
px = c(-3:3 * 2/7);
curve(y(x, a=a), from= -1.5, to = 1.5, ylim=c(-10, 10))
sapply(px * 1.2, line.tan, dx=3, p=y, dp=dy, a=a)
#
curve(dy(x, a=a), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, a=a, col="orange")


### Ex 2:
a = c(1.5, -2) # a[] has NO effect on eq of D2;
px = c(-3:3 * 2/7);
curve(y(x, a=a), from= -1.5, to = 1.5, ylim=c(-12, 8))
sapply(px * 1.2, line.tan, dx=2.2, p=y, dp=dy, a=a)
#
curve(dy(x, a=a), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, a=a, col="orange")


### Partial Extension
### y = m*x^m*e^(k*x^n) - n*x^n*e^(k*x^m) + F0(x)

### D(y)
m^2*x^(m-1)*e^(k*x^n) - n^2*x^(n-1)*e^(k*x^m) + df0

### Solve Liniar system:
e^(k*x^n) = (n*(y - f0) - x*(dy - df0)) / (m*(n-m)*x^m)
e^(k*x^m) = (m*(y - f0) - x*(dy - df0)) / (n*(n-m)*x^n)

### D2(y)
m^2*(m-1)*x^(m-2)*e^(k*x^n) - n^2*(n-1)*x^(n-2)*e^(k*x^m) +
	+ k*m^2*n*x^(m+n-2)*e^(k*x^n) - k*m*n^2*x^(m+n-2)*e^(k*x^m) + d2f0
### (n-m)*x^2*D2(y) =
m*(m-1)*(n*y - x*dy) - n*(n-1)*(m*y - x*dy) +
	+ k*m*n*x^n*(n*y - x*(dy - df0)) - k*m*n*x^m*(m*y - x*(dy - df0)) +
	- m*n*(m-n)*f0 - k*m*n*(n*x^n - m*x^m)*f0 +
	+ (m-n)*(m+n-1)*x*df0 +
	+ (n-m)*x^2*d2f0
k*m*n*x*(x^m - x^n)*dy + (n-m)*(n + m - 1)*x*dy +
	+ k*m*n*(n*x^n - m*x^m)*y + m*n*(m-n)*y +
	- k*m*n*(n*x^n - m*x^m)*f0 - m*n*(m-n)*f0 +
	+ k*m*n*x*(x^n - x^m)*df0 + (m-n)*(m + n - 1)*x*df0 + (n-m)*x^2*d2f0

### Examples:
### f0 = x^2
### (n-m)*x^2*D2(y) =
k*m*n*x*(x^m - x^n)*dy + (n-m)*(n + m - 1)*x*dy +
	+ k*m*n*(n*x^n - m*x^m)*y + m*n*(m-n)*y +
	- k*m*n*(n*x^n - m*x^m)*x^2 + 2*k*m*n*(x^n - x^m)*x^2 +
	+ 2*(m-n)*(m+n - 2)*x^2 - m*n*(m-n)*x^2;


### Solution & Plot:
y = function(x, m=c(2,3), k=1, a=c(1, 1)) {
	xm = x^m[1]; xn = x^m[2];
	y = a[1]*m[1]*xm*exp(k*xn) - a[2]*m[2]*xn*exp(k*xm) + x^2;
	return(y)
}
dy = function(x, m=c(2,3), k=1, a=c(1, 1)) {
	# TODO: a;
	xm = x^m[1]; xn = x^m[2];
	dy.v = m[1]^2*xm*exp(k*xn) - m[2]^2*xn*exp(k*xm) + 2*x^2; # x*df0
	div = x;
	dy.v = ifelse(div != 0, dy.v / div, dy(x + 1E-3, m=m, k=k, a=a)); # TODO
	return(dy.v)
}
d2y = function(x, m=c(2,3), k=1, a=c(1, 1)) {
	y.x = y(x, m=m, k=k, a=a)
	dy.x = dy(x, m=m, k=k, a=a)
	xm = x^m[1]; xn = x^m[2]; x2 = x^2;
	n = m[2]; m = m[1];
	mn = m*n; kmn = k*mn;
	dp = kmn*x*(xm - xn)*dy.x + (n-m)*(n + m - 1)*x*dy.x +
		+ kmn*(n*xn - m*xm)*y.x + mn*(m-n)*y.x +
		- kmn*(n*xn - m*xm)*x2 +
		+ 2*kmn*(xn - xm)*x2 + 2*(m-n)*(n + m - 2)*x2 - mn*(m-n)*x2;
	div = (n-m)*x2;
	dp = ifelse(div != 0, dp / div, d2y(x + 1E-2, m=c(m,n), k=k, a=a)); # TODO
	return(dp)
}
### Plot:
m = c(2, 3); k = 1;
a = c(1, 1) # a[] has NO effect on eq of D2;
px = c(-3:5 * 1/7);
curve(y(x, m=m, k=k, a=a), from= -1.5, to = 1.5, ylim=c(-10, 10))
sapply(px * 2, line.tan, dx=3, p=y, dp=dy, m=m, k=k, a=a)
#
curve(dy(x, m=m, k=k, a=a), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, m=m, k=k, a=a, col="orange")


### Ex 2:
m = c(2, 3); k = -1;
a = c(1, 1) # a[] has NO effect on eq of D2;
px = c(-2:1 * 2/7, 0.8, 1.3);
curve(y(x, m=m, k=k, a=a), from= -1.5, to = 1.8, ylim=c(-10, 10))
sapply(px * 1.2, line.tan, dx=3, p=y, dp=dy, m=m, k=k, a=a)
#
curve(dy(x, m=m, k=k, a=a), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, m=m, k=k, a=a, col="orange")


####################

# y = exp(1/(x+b)) + exp(-1/(x+b)) + f0(x)

# D =>
(x+b)^2 * dy + (exp(1/(x+b)) - exp(-1/(x+b))) - (x+b)^2 * df0 # = 0

# System: (but not needed)
# 2*exp(1/(x+b)) =
-(x+b)^2 * dy + y + (x+b)^2 * df0 - f0;
# 2*exp(-1/(x+b)) =
(x+b)^2 * dy + y - (x+b)^2 * df0 - f0;

### ODE:
(x+b)^4 * d2y + 2*(x+b)^3 * dy - y +
	- (x+b)^4 * d2f0 - 2*(x+b)^3 * df0 + f0 # = 0


### Solution & Plot:
y = function(x, b = sqrt(5), FUN = NULL) {
	xi = 1 / (x + b);
	y = exp(xi) + exp(- xi) + FUN(x);
	return(y)
}
dy = function(x, b = sqrt(5), FUN = NULL) {
	xb = x + b; xi = 1 / xb;
	div = xb^2;
	dy = exp(- xi) - exp(xi) + div * FUN(x);
	dy = dy / div;
	# dy = ifelse(div != 0, dy / div, dy(x + 1E-3, b=b, FUN=FUN));
	return(dy)
}
d2y = function(x, b = sqrt(5), FUN = NULL, DFUN = NULL, D2FUN) {
	# (x+b)^4 * d2y + 2*(x+b)^3 * dy - y +
	# - (x+b)^4 * d2f0 - 2*(x+b)^3 * df0 + f0 # = 0
	y  =  y(x, b=b, FUN=FUN);
	dy = dy(x, b=b, FUN=DFUN);
	xb = x + b; xb3 = xb^3;
	div = xb^4;
	dp  = - 2*xb3 * dy + y +
		+ div * D2FUN(x) + 2*xb3 * DFUN(x) - FUN(x);
	dp = dp / div;
	# dp = ifelse(div != 0, dp / div, d2y(x + 1E-2, ...);
	return(dp)
}
### Plot:
b = sqrt(5);
FUN   = \(x) x*exp(-x^2 / 3);
DFUN  = \(x) -1/3 * (2*x^2 - 3)*exp(-x^2 / 3);
D2FUN = \(x)  1/9 * (4*x^3 - 6*x - 12*x)*exp(-x^2 / 3);
yf   = \(...)   y(..., FUN=FUN);
dyf  = \(...)  dy(..., FUN=DFUN);
d2yf = \(...) d2y(..., FUN=FUN, DFUN=DFUN, D2FUN=D2FUN);
px = c(-0.75, 0:4);
curve(y(x, b=b, FUN=FUN), from= -1, to = 4.5, ylim=c(-1.5, 4.5))
sapply(px, line.tan, dx=3, p=yf, dp=dyf, b=b)
#
curve(dy(x, b=b, FUN = DFUN), add=T, col="green")
sapply(px, line.tan, dx=3, p=dyf, dp=d2yf, b=b, col="orange")

