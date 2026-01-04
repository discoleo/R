########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs - Logarithms
##
## draft v.0.3p


### ODEs Derived from Logarithms

### Note:
# Lambert W examples:
# - should be in a different file;

###############

###############
### History ###
###############

### draft v.0.3m - v.0.3p:
# - [refactor] moved ODEs based on Hidden Log
#   to new file: DE.ODE.NL.Log.Other.R;
# - [refactor] Prod( LOG ) moved to new file:
#   DE.ODE.NL.Log.Prod.R;
# - [refactor] Int( LOG ) moved to new file:
#   DE.ODE.NL.Log.Int.R;
### draft v.0.3i - v.0.3j:
# - Mixed variants: Sum(LOG, EXP)
#   y = B1(x)*log(P1(x)) + B2(x)*exp(P2(x)) + F0(x);
# - Sum(LOG, I( EXP )):
#   y = B1(x)*log(P1(x)) + B2(x)*exp(P2(x))*I( exp(-P2(x)) ) + F0(x); [v.0.3j]
### draft v.0.3h:
# - workout of case:
#   y = log(P) * log(log(P));
### draft v.0.3g - v.0.3g-ex2:
# - Automatic generation of simple types of ODEs;
# - more examples & checks; [v.0.3g-ex2]
### draft v.0.3f - v.0.3f-st1:
# - derived from: [TODO full derivation]
#   y = log(exp(P1(x)) + P2(x)) + F0(x); [started]
#   y = log(log(P(x))) + F0(x);
#   y = log(P(x)) * log(log(P(x))) + F0(x); [v.0.3f-bis]
### draft v.0.3e:
# - derived from:
#   y*log(x+k) = I(log(x+k) / x) + F0(x);
### draft v.0.3d:
# - cleanup;
### draft v.0.3c - v.0.3c-ex2:
# - Mixed Log-Exp:
#   y = log(P1(x)) * exp(P2(x)) + F(x);
#   [refactor] moved to file: DE.ODE.Mixed.Exp.Log.R;
# - more examples; [v.0.3c-ex2]
### draft v.0.3a - v.0.3b:
#   [refactor] moved to new file;
# - from:
#   y = (x + k1)^(x + k2) + F0(x);
#   y * (x + k1)^(x + k2) = F0(x);
# - example:
#   (x+k1)*(y - f0)*d2y = (x+k1)*(dy - df0)^2 + (y - f0)^2 + (x+k1)*d2f0*(y - f0);
### draft v.0.2c:
# - cleanup:
#   moved section: y = log(P1(x))*log(P2(x))
#   from DE.ODE.Fractions.Lambert.R to this file;
### draft v.0.2b:
# - derived from:
#   y * I(1/log(x + k)) dx = F0(x);
### draft v.0.2a:
# - derived from:
#   y = x * I(1/log(x + k)) dx + F0(x);
### draft v.0.1d - v.0.1g:
# - derived from: y = x^log(x)
#   x^2*y*d2y - x^2*dy^2 + x*y*dy - 2*y^2 = 0;
# - extension: y = (x+k)^log(x+k)
#  (x+k)^2*y*d2y - (x+k)^2*dy^2 + (x+k)*y*dy - 2*y^2 = 0;
# - generalization to higher orders:
#   y = (x^n+k)^log(x^n+k); [v.0.1f]
# - generalization:
#   y = P(x)^log(P(x));


#########################
#########################

### Helper Functions

# include: Polynomials.Helper.ODE.R;
# include: Polynomials.Helper.R; [automatically]
# include: DE.ODE.Helper.R;
source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################


### Section A: Linear ODEs

### Derived from:
### y = F1(x)*log(P1(x)) + F2(x)*log(P2(x))

### Automatic:

### Simple Example:
p1 = toPoly.pm("x^2 + c")
p2 = toPoly.pm("x^2 - c")
pL1 = toPoly.pm("x^2 + c")
pL2 = toPoly.pm("x^2 - c")
pDiv = toPoly.pm("c*x^8 - 2*x^4*c^3 + c^5")
pR = genODE.Log.pm(p1, p2, pL1, pL2, pDiv=pDiv, div.by="x");
print.dpm(pR, do.sort=FALSE)

### ODE:
x*(x^4 - c^2)*d2y - (x^4 - c^2)*dy - 8*x^5 # = 0


### Simple Example 2:
n = 3
p1 = toPoly.pm("x^n + c")
p2 = toPoly.pm("x^n - c")
pL1 = toPoly.pm("x^n + c")
pL2 = toPoly.pm("x^n - c")
pDiv = toPoly.pm("c*x^(4*n) - 2*x^(2*n)*c^3 + c^5");
pR = genODE.Log.pm(p1, p2, pL1, pL2, pDiv=pDiv, div.by="x");
print.dpm(pR, do.sort=FALSE)

### ODE:
x^2*(x^6 - c^2)*d2y - 2*x*(x^6 - c^2)*dy - 18*x^9 # = 0


### Example 3:
n = 2
p1 = toPoly.pm("x^n + c")
p2 = toPoly.pm("x^n - c")
pL1 = toPoly.pm("x^n - c")
pL2 = toPoly.pm("x^n + c")
pDiv = toPoly.pm("c*x^(2*n) - c^3")
pR = genODE.Log.pm(p1, p2, pL1, pL2, pDiv=pDiv, div.by="x");
print.dpm(pR, do.sort=FALSE)

### ODE:
x*(x^4 - c^2)^2*d2y - (x^4 - c^2)^2*dy - 8*x^5*(x^4 - 5*c^2) # = 0


### Example 4:
n = 2; m = 4; k = 4;
p1 = toPoly.pm("x^n + c")
p2 = toPoly.pm("x^n - c")
f0 = toPoly.pm("k()*x^m") # evaluate k();
pL1 = toPoly.pm("x^n + c")
pL2 = toPoly.pm("x^n - c")
pDiv = toPoly.pm("c*x^(4*n) - 2*c^3*x^(2*n) + c^5")
pR = genODE.Log.pm(p1, p2, pL1, pL2, f0=f0, pDiv=pDiv, div.by="x");
print.dpm(pR, do.sort=FALSE)

### ODE:
x*(x^4 - c^2)*d2y - (x^4 - c^2)*dy - 8*x^5 + k*m*(m-2)*(x^4 - c^2)*x^(m-1) # = 0


### Example 5:
# - non-correlated: a & c;
p1 = toPoly.pm("x^2 + c")
p2 = toPoly.pm("x^2 - c")
pL1 = toPoly.pm("x^2 + a")
pL2 = toPoly.pm("x^2 - a")
pDiv = toPoly.pm("c*x^4 - c*a^2")
pR = genODE.Log.pm(p1, p2, pL1, pL2, pDiv=pDiv, div.by="x");
print.dpm(pR, do.sort=FALSE)

### ODE:
(x^9 - 2*a^2*x^5 + a^4*x)*d2y - (x^8 - 2*a^2*x^4 + a^4)*dy +
	- 8*x^9 + 24*a^2*x^5 - 16*a*c*x^5 # = 0


### Example 6:
Log.ODE.gen = function(n, m, k, print=TRUE) {
	p1 = toPoly.pm("x^n")
	p2 = toPoly.pm("1")
	f0 = toPoly.pm("k()*x^m") # evaluate k();
	pL1 = toPoly.pm("x + b")
	pL2 = toPoly.pm("x + c")
	pDiv = toPoly.pm("x^2 + b*x + c*x + b*c")
	pR = genODE.Log.pm(p1, p2, pL1, pL2, f0=f0, pDiv=pDiv, div.by="x");
	if(print) print.dpm(pR, do.sort=FALSE)
	invisible(pR);
}
### Ex 6a:
n = 1; m = 0; k = 4;
Log.ODE.gen(n, m, k)

### ODE:
(x + b)^2*(x + c)^2*d2y +
	- x^3 + x^2 - 2*b*x^2 - 2*c*x^2 + 2*b*x - 4*b*c*x - c^2*x + b^2 - 2*b*c^2

### Ex 6b:
n = 1; m = 2; k = 3;
pR0 = Log.ODE.gen(n, m=0, k=k, print=FALSE);
pR  = diff.pm(pR0, Log.ODE.gen(n, m=m, k=k, print=FALSE))
print.dpm(pR)
# Diff(ODE)
# k*m*(m-1)*(x+b)^2*(x+c)^2*x^(m-2)
# Check:
div.pm(pR, pow.pm(toPoly.pm("x^2 + b*x + c*x + b*c"), 2), by="x")


################
### Explicit ###
################

### Examples:
### y = (x^m + c)*log(x^n + a) + (x^m - c)*log(x^n + b)
# - when b = -a: simplifies massively & eliminates y from D2(y);
# - does NOT simplify when: F1 - F2 != constant;

# TODO: check computations;

### D(y)
# [not run]
m*x^(m-1)*log(x^n + a) + m*x^(m-1)*log(x^n + b) +
	+ n*x^(n-1)*(x^m + c) / (x^n + a) +
	+ n*x^(n-1)*(x^m - c) / (x^n + b)
m*x^(m-1)*log(x^n + a) + m*x^(m-1)*log(x^n + b) +
	+ n*x^(n-1)*(2*x^(m+n) + (a+b)*x^m - c*(a-b)) / ((x^n + a) * (x^n + b))
# let: T(x) = x^(n-1)*(2*x^(m+n) + (a+b)*x^m - c*(a-b)) / ((x^n + a) * (x^n + b))
# Note: dy = ... + n * T(x); =>
m*x^(m-1)*log(x^n + a) + m*x^(m-1)*log(x^n + b) + n*T;

### Solve linearly for log(P(x))
# Check:
# m = 4; toPoly.pm("m*x^(m-1)*(x^m - c) - m*x^(m-1)*(x^m + c)");
log(x^n + a) =
	((x^m - c)*dy - m*x^(m-1)*y - n*T*(x^m - c)) / - (2*c*m*x^(m-1));
log(x^n + b) =
	((x^m + c)*dy - m*x^(m-1)*y - n*T*(x^m + c)) /   (2*c*m*x^(m-1));
log(x^n + a) + log(x^n + b) = (dy - n*T) / (m*x^(m-1));

### D(y):
(x^n + a)*(x^n + b) * dy +
	- m*x^(m-1)*(x^n + a)*(x^n + b)*(log(x^n + a) + log(x^n + b)) +
	- n*x^(n-1)*(2*x^(m+n) + (a+b)*x^m - c*(a-b));
T = x^(n-1)*(2*x^(m+n) + (a+b)*x^m - c*(a-b)) / ((x^n + a) * (x^n + b))

### D2(y)
(x^n + a)*(x^n + b) * d2y.x + n*x^(n-1)*(2*x^n + a + b) * dy.x +
	- (m*(m-1)*x^(m-2)*(x^n + a)*(x^n + b) + m*n*x^(m+n-2)*(2*x^n + a + b)) *
		(log(x^n + a) + log(x^n + b)) +
	- n*x^(m+n-2)*(2*m*x^n + m*(a + b)) +
	- n*x^(m+n-2)*(2*(m+n)*x^n + m*(a+b)) +
	- n*(n-1)*x^(n-2)*(2*x^(m+n) + (a+b)*x^m - c*(a-b))
(x^n + a)*(x^n + b) * d2y.x + n*x^(n-1)*(2*x^n + a + b) * dy.x +
	- (m*(m-1)*x^(m-2)*(x^n + a)*(x^n + b) + m*n*x^(m+n-2)*(2*x^n + a + b)) *
		(log(x^n + a) + log(x^n + b)) +
	- n*x^(m+n-2)*(2*(2*m+2*n-1)*x^n + (2*m+n-1)*(a+b)) +
	+ n*(n-1)*c*(a-b)*x^(n-2)
x*(x^n + a)*(x^n + b) * d2y + n*x^n*(2*x^n + a + b) * dy +
	- ((m-1)*(x^n + a)*(x^n + b) + n*x^n*(2*x^n + a + b)) *
		(dy - n*T) +
	- n*x^(m+n-1)*(2*(2*m+2*n-1)*x^n + (2*m+n-1)*(a+b)) +
	+ n*(n-1)*c*(a-b)*x^(n-1)
x*(x^n + a)*(x^n + b) * d2y +
	- (m-1)*x^n*(x^n + (a+b)) * dy  - (m-1)*a*b*dy + # - (m-1)*(x^n+a)*(x^n+b)*dy
	+ n*((m-1)*(x^n + a)*(x^n + b) + n*x^n*(2*x^n + a + b)) * T +
	- n*x^(m+n-1)*(2*(2*m+2*n-1)*x^n + (2*m+n-1)*(a+b)) +
	+ n*(n-1)*c*(a-b)*x^(n-1)
x*(x^n + a)*(x^n + b) * d2y +
	- (m-1)*(x^n + a)*(x^n + b)*dy +
	- n*x^(m+n-1)*(2*(m+2*n)*x^n + (m+n)*(a+b)) +
	+ n*(n-m)*c*(a-b)*x^(n-1) +
	+ n^2*x^n*(2*x^n + a + b) * T
# D(D(y) / x^(m-1)) = ... * x^(m-1) / (x*(x^n + a)*(x^n + b));

n = 3; m = 4;
p = toPoly.pm("(x^n + a)*(x^n + b) * dy +
	- m*x^(m-1)*(x^n + a)*(x^n + b)*(L1 + L2) +
	- n*x^(n-1)*(2*x^(m+n) + (a+b)*x^m - c*(a-b))");
dy = lapply(c("y", "L1", "L2"), function(y) dy.pm(p, y, x=NULL));
dx = dp.pm(p, "x");
pR = sum.lpm(c(dy, list(dx)));
dL = toPoly.pm("n*x^(n-1)");
pLa = toPoly.pm("x^n + a"); pLb = toPoly.pm("x^n + b");
pR = replace.fr.pm(pR, dL, pLa, "dL1");
pR = replace.fr.pm(pR, dL, pLb, "dL2");
# TODO

### Special Case 1:
m = 2;
b = -a; # a+b = 0;
T = 2*x^(n-1)*(x^(n+2) - c*a) / ((x^n + a) * (x^n - a))
### ODE:
x*(x^n + a)*(x^n - a) * d2y - (x^n + a)*(x^n - a)*dy +
	- 4*n*(n+1)*x^(2*n+1) + c*n*(n-2)*(a-b)*x^(n-1) +
	+ 2*n^2*x^(2*n) * T # = 0


### Solution & Plot:
y = function(x, a, b=-a, ct=0, n=2, m=2) {
	xm = x^m; xn = x^n;
	x.prod = (xn + a)*(xn +b); x.div = (xn + a) / (xn + b);
	val = xm*log(x.prod) + ct*log(x.div);
	return(val)
}
dy = function(x, a, b=-a, ct=0, n=2, m=2) {
	xm = x^m; xn = x^n;
	# y.x = y(x, a=a, b=b, ct=ct, n=n, m=m);
	x.prod = (xn + a)*(xn + b);
	div = x;
	Tnx = n*xn*(2*xm*xn + (a+b)*xm - ct*(a-b)) / (x.prod);
	dp = m*xm*log(x.prod) + Tnx;
	dp = ifelse(div != 0, dp/div, 1); # TODO
	return(dp)
}
d2y = function(x, a, b=-a, ct=0, n=2, m=2) {
	xm = x^m; xn = x^n;
	x.prod = (xn + a)*(xn + b);
	### T
	div = x * x.prod
	T = xn*(2*xm*xn + (a+b)*xm - ct*(a-b)) / div;
	### d2y
	# y.x  =  y(x, a=a, b=b, ct=ct, n=n, m=m);
	dy.x = dy(x, a=a, b=b, ct=ct, n=n, m=m);
	d2p =
	(m-1) * x.prod * dy.x +
		+ n*(xm*(2*(m+2*n)*xn + (m+n)*(a+b)) - (n-m)*ct*(a-b)) * xn / x +
		- n^2*xn*(2*xn + a + b) * T;
	d2p / div;
}
### Plot:
n = 2; m = 2;
a = 1; b = -a;
ct = 2;
px = 6/7 + (1:4)*2/7
curve(y(x, a=a, b=b, ct=ct, n=n, m=m), from= 1, to = 2.5, ylim=c(-2,20))
line.tan(px, dx=3, p=y, dp=dy, a=a, b=b, ct=ct, n=n, m=m)
# x*log(x)-like:
curve(dy(x, a=a, b=b, ct=ct, n=n, m=m), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, a=a, b=b, ct=ct, n=n, m=m, col="orange")


### Ex 2:
n = 2; m = 1/2;
a = 1; b = 3/2;
ct = 2;
px = 6/7 + c(1, 4)*2/7
curve(y(x, a=a, b=b, ct=ct, n=n, m=m), from= 1, to = 2.5, ylim=c(1/2,7))
line.tan(px, dx=3, p=y, dp=dy, a=a, b=b, ct=ct, n=n, m=m)
#
curve(dy(x, a=a, b=b, ct=ct, n=n, m=m), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, a=a, b=b, ct=ct, n=n, m=m, col="orange")


############
### Example:
y = (x^2 + k*x)*log(x^n + a) + (x^2 - k*x)*log(x^n + b)

### D(y)
(2*x + k)*log(x^n + a) + (2*x - k)*log(x^n + b) + n*x^n*((x+k)/(x^n+a) + (x-k)/(x^n+b))
(2*x + k)*log(x^n + a) + (2*x - k)*log(x^n + b) +
	+ n*x^n*(2*x^(n+1) + (a+b)*x - k*(a-b))/((x^n+a)*(x^n+b))

### D2(y)
2*(log(x^n + a) + log(x^n + b)) +
	+ n*x^(n-1)*(2*x+k) / (x^n+a) + n*x^(n-1)*(2*x-k) / (x^n+b) +
	+ D( n*x^n*(2*x^(n+1) + (a+b)*x - k*(a-b))/((x^n+a)*(x^n+b)) );
2*(x*dy - y) / x^2 +
	- n*x^(n-1)*(2*x^(n+1) + (a+b)*x - k*(a-b))/((x^n+a)*(x^n+b)) + # from x*dy
	+ n*x^(n-1)*(2*x+k) / (x^n+a) + n*x^(n-1)*(2*x-k) / (x^n+b) +
	+ D( n*x^n*(2*x^(n+1) + (a+b)*x - k*(a-b))/((x^n+a)*(x^n+b)) );
# TODO: ...


############
### Example:
y = x^3*log(x^n + a) + x^2*log(x^n + b)

### D(y)
3*x^2*log(x^n + a) + 2*x*log(x^n + b) + n*x^(n+2)/(x^n+a) + n*x^(n+1)/(x^n+b)

### D2(y)
6*x*log(x^n + a) + 2*log(x^n + b) +
	+ 3*n*x^(n+1) / (x^n + a) + 2*n*x^n / (x^n + b) +
	+ D( n*x^(n+2)/(x^n+a) + n*x^(n+1)/(x^n+b) )
2*(2*dy - 3*y) / x^2 +
	+ 3*n*x^(n+1) / (x^n + a) + 2*n*x^n / (x^n + b) +
	+ D( n*x^(n+2)/(x^n+a) + n*x^(n+1)/(x^n+b) )


############
### Example:
y = sin(x)*log(p1) + cos(x)*log(p2)

### D(y)
cos(x)*log(p1) - sin(x)*log(p2) + sin(x) / p1 * dp1 + cos(x) / p2 * dp2

### D2(y)
- sin(x)*log(p1) - cos(x)*log(p2) +
	+ cos(x) / p1 * dp1 - sin(x) / p2 * dp2 +
	D( sin(x) / p1 * dp1 + cos(x) / p2 * dp2 )
- y +
	+ cos(x) / p1 * dp1 - sin(x) / p2 * dp2 +
	D( sin(x) / p1 * dp1 + cos(x) / p2 * dp2 )


#########################

###############
### Mixed:  ###
### Log-Exp ###
###############

### Product:
### y = log(P1(x)) * exp(P2(x)) + F(x);
# - [refactor] moved to file:
#   DE.ODE.Mixed.Exp.Log.R;

#################

#################
### Mixed Sum ###
#################

### y = B1(x) * log(P1(x)) + B2(x) * exp(P2(x)) + F0(x)

### Example:
# y = x*log(x) + x^2*exp(k/x) + f0

### D =>
# dy = log(x) + (2*x - k)*exp(k/x) + df0 + 1;
dy - log(x) - (2*x - k)*exp(k/x) - df0 - 1 # = 0

### Linear System =>
# log(x) =
(x^2*dy - (2*x-k)*y - x^2*df0 - x^2 + (2*x-k)*f0) / - (x^2 - k*x);
# exp(k/x) =
(x*dy - y + f0 - x*df0 - x) / (x^2 - k*x);

### D2 =>
x^2*d2y - 2*x^2*exp(x/k) + k*(2*x - k) * exp(k/x) - x^2*d2f0 - x # = 0
x^2*(x^2 - k*x)*d2y - (2*x^2 - 2*x*k + k^2)*(x*dy - y + f0 - x*df0 - x) +
	- x^2*(x^2 - k*x)*d2f0 - x*(x^2 - k*x) # = 0

### ODE:
x^2*(x^2 - k*x)*d2y - (2*x^2 - 2*x*k + k^2)*(x*dy - y) +
	- (2*x^2 - 2*x*k + k^2)*(f0 - x*df0 - x) - x^2*(x^2 - k*x)*d2f0 - x*(x^2 - k*x) # = 0


### Solution & Plot:
y = function(x, k=1, f=NULL) {
	val = x*log(x) + x^2*exp(k/x);
	if( ! is.null(f)) {
		fx = eval.vpm(f, x);
		val = val + fx;
	}
	return(val)
}
dy = function(x, k=1, f=NULL) {
	dyx = log(x) + (2*x - k)*exp(k/x) + 1;
	if( ! is.null(f)) {
		df0 = dp.pm(f, xn="x");
		dyx = dyx + eval.vpm(df0, x);
	}
	return(dyx)
}
d2y = function(x, k=1, f=NULL) {
	yx  =  y(x, k=k, f=f);
	dyx = dy(x, k=k, f=f);
	#
	px  = (2*x^2 - 2*x*k + k^2);
	x2k = (x^2 - k*x);
	d2p = px*(x*dyx - yx - x) + x*x2k;
	if( ! is.null(f)) {
		fx0  = eval.vpm(f, x);
		df0  = dp.pm(f, xn="x");
		d2f0 = dp.pm(df0, xn="x");
		dfx0 = eval.vpm(df0, x);
		d2fx = eval.vpm(d2f0, x);
		d2p = d2p + px*(fx0 - x*dfx0) + x^2*x2k*d2fx;
	}
	div = x^2*x2k;
	d2p = ifelse(div != 0, d2p/div, 1); # TODO: check;
	return(d2p);
}
### Plot:
k = 3
px = (1:3)*2/7; px = c(px + 0.25, px + 1);
curve(y(x, k=k), from = 0.5, to = 2, ylim = c(10, 40))
line.tan(px, dx=3, p=y, dp=dy, k=k)

#
curve(dy(x, k=k), from = 0.5, to = 2, ylim=c(-50, 10), col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, col="orange")


### Ex 2:
k = -2;
f = toPoly.pm("x^2 + 3*x")
px = (1:3)*1/7; px1 = c(2*px, 2*px + 1); px2 = c(px, 2*px + 1);
curve(y(x, k=k, f=f), from = 0, to = 2, ylim = c(0, 12))
line.tan(px1, dx=3, p=y, dp=dy, k=k, f=f)
#
curve(dy(x, k=k, f=f), add=TRUE, col="green")
line.tan(px2, dx=3, p=dy, dp=d2y, k=k, f=f, col="orange")

##################

##################
### Mixed Sum  ###
### LOG + I()  ###
##################

### y = B1(x) * log(P1(x)) + B2(x) * exp(P2(x)) * I( exp( - P2(x)) ) + F0(x)

### Example:
# y = x*log(x) + x^2 * exp(k/x) * I( exp(-k/x) ) + f0

### D =>
# x^2 * dy = x^2 * log(x) + (2*x - k)*(y - x*log(x) - f0) + x^2*df0 + x^4 + x^2;
# x^2 * dy = - (x^2 - k*x)*log(x) + (2*x - k)*(y - f0) + x^2*df0 + x^4 + x^2;

### Linear system:
# - x*(x - k) * log(x) =
x^2*dy - (2*x - k)*(y - f0) - x^2*df0 - x^4 - x^2;

### D2 =>
x^2*d2y + 2*x*dy # =
- (x^2 - k*x)/x - (2*x - k)*log(x) + (2*x - k)*(dy - df0) + 2*(y - f0) +
	+ x^2*d2f0 + 2*x*df0 + 4*x^3 + 2*x;
# =>
x^3*d2y + 2*x^2*dy + x*(2*x - k)*log(x) - x*(2*x - k)*(dy - df0) - 2*x*(y - f0) +
	- x^3*d2f0 - 2*x^2*df0 - 4*x^4 - x^2 - k*x # = 0
x^3*(x - k)*d2y + 2*x^2*(x - k)*dy +
	- (2*x - k)*(x^2*dy - (2*x - k)*(y - f0) - x^2*df0 - x^4 - x^2) +
	- x*(x - k)*(2*x - k)*(dy - df0) - 2*x*(x - k)*(y - f0) +
	- (x - k)*(x^3*d2f0 + 2*x^2*df0 + 4*x^4 + x^2 + k*x) # = 0

### ODE:
# P(d2y, dy, y) = P(d2f, df, f) + someTrP(x);
x^3*(x - k)*d2y - x*(2*x^2 - 2*k*x + k^2)*dy + (2*x^2 - 2*k*x + k^2)*y +
	- 2*x^5 + 3*k*x^4 + x^3 - k*x^2 + k^2*x +
	- x^3*(x - k)*d2f0 + x*(2*x^2 - 2*k*x + k^2)*df0 - (2*x^2 - 2*k*x + k^2)*f0 # = 0

# TODO: check;


###
# p = toPoly.pm("")
p = sort.pm(p, xn=c("d2y", "dy", "y"), xn2 = c("x", "d2f0", "df0", "f0"), sort.coeff=c(10:12,13:16))
p = as.pm.first(p, "k")
print.pm(p, lead="d2y", do.sort=FALSE)


#########################
#########################

#######################
### Section B:      ### 
### Non-Linear ODEs ###
#######################

####################
### Logarithmic  ###
### Higher Power ###
####################

### Derived from:
### y = (log(P1(x)))^2 + (log(P2(x)))^2

### Example:
y = (log(x^2 + a))^2 + (log(x^2 + b))^2

### D(y)
4*x*log(x^2 + a) / (x^2 + a) + 4*x*log(x^2 + b) / (x^2 + b)
### (x^2+a)*(x^2+b)*dy
4*x*(x^2+b)*log(x^2 + a) + 4*x*(x^2+a)*log(x^2 + b)

### (x^2+a)*(x^2+b)*D2(y) + D((x^2+a)*(x^2+b))*dy
4*(3*x^2 + b)*log(x^2 + a) + 4*(3*x^2 + a)*log(x^2 + b) +
	+ 8*x^2*(x^2+b) / (x^2 + a) + 8*x^2*(x^2+a) / (x^2 + b)
### Solve Linear system:
# T = 8*x^2*(x^2+b) / (x^2 + a) + 8*x^2*(x^2+a) / (x^2 + b) - D((x^2+a)*(x^2+b))*dy;
log(x^2 + a) =
	(x*(x^2+a)^2*(x^2+b)*d2y - (3*x^2 + a)*dy - x*(x^2+a)*T) / (8*(a-b)*x^3)
log(x^2 + b) =
	(x*(x^2+b)^2*(x^2+a)*d2y - (3*x^2 + b)*dy - x*(x^2+b)*T) / -(8*(a-b)*x^3)
# TODO: check!

### ODE:
y +
	(x*(x^2+a)^2*(x^2+b)*d2y - (3*x^2 + a)*dy - x*(x^2+a)*T) *
	(x*(x^2+b)^2*(x^2+a)*d2y - (3*x^2 + b)*dy - x*(x^2+b)*T) / (64*x^6) # = 0


######################
######################

### y = x^log(x)

### D(y)
2*log(x)/x * x^log(x)
2*log(x)/x * y

### D2(y)
2*log(x)/x * dy + 2*(1 - log(x))/x^2 * y
### x^2 * D2(y)
2*x*log(x) * dy + 2*(1 - log(x)) * y
x^2/y * dy^2 - x*dy + 2*y

### ODE:
x^2*y*d2y - x^2*dy^2 + x*y*dy - 2*y^2 # = 0


### Extension
### y = (x + k)^log(x + k)

### D(y)
2*log(x+k)*y / (x+k)

### ODE:
(x+k)^2*y*d2y - (x+k)^2*dy^2 + (x+k)*y*dy - 2*y^2 # = 0


### Solution & Plot:
y = function(x, k=0) {
	xk = if(k == 0) x else x+k;
	val = xk^log(xk);
	return(val)
}
dy = function(x, k=0) {
	xk = x + k;
	logx = log(xk);
	dp = 2*logx * (xk)^logx
	div = xk;
	dp = ifelse(div != 0, dp/div, 1); # TODO
	return(dp)
}
d2y = function(x, k=0) {
	### D()
	y.x  =  y(x, k=k);
	dy.x = dy(x, k=k);
	xk = if(k == 0) x else x+k;
	div = xk^2 * y.x;
	d2p = xk^2*dy.x^2 - xk*y.x*dy.x + 2*y.x^2
	d2p = ifelse(div != 0, d2p/div, 1); # TODO
	return(d2p)
}
### Plot:
px = 3/7 + (0:4)*2/7
curve(y(x), from= 0+1E-1, to = 2.5, ylim=c(-2,5))
line.tan(px, dx=3, p=y, dp=dy)
#
curve(dy(x), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, col="orange")


### Ex 2:
k = 2
px = 3/7 + (0:4)*2/7 - k;
curve(y(x, k=k), from= -k+1E-1, to = 2.5, ylim=c(-2,5))
line.tan(px, dx=3, p=y, dp=dy, k=k)
#
curve(dy(x, k=k), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, col="orange")


### Extensions: Higher Order
### y = (x^n + k)^log(x^n + k)

### D(y)
2*n*log(x^n+k)*x^(n-1) / (x^n+k) * y

### D2(y)
2*n*log(x^n+k)*x^(n-1) / (x^n+k) * dy +
	+ 2*n^2*x^(2*n-2) / (x^n+k)^2 * y +
	+ 2*n*log(x^n+k)*((n-1)*x^(n-2)*(x^n+k) - n*x^(2*n-2)) / (x^n+k)^2 * y
(dy)^2 / y +
	+ 2*n^2*x^(2*n-2) / (x^n+k)^2 * y +
	- (x^n + k - n*k) / (x^n+k) * dy / x

### ODE:
x*(x^n+k)^2*y*d2y - x*(x^n+k)^2*dy^2 + (x^n+k-n*k)*(x^n+k)*y*dy - 2*n^2*x^(2*n-1)*y^2 # = 0

### Case: n = 2
x*(x^2+k)^2*y*d2y - x*(x^2+k)^2*dy^2 + (x^4-k^2)*y*dy - 8*x^3*y^2 # = 0


### Solution & Plot:
y = function(x, k=0, n=2) {
	xk = if(k == 0) x^n else x^n + k;
	val = xk^log(xk);
	return(val)
}
dy = function(x, k=0, n=2) {
	xn = x^n;
	xk = if(k == 0) xn else xn + k;
	logx = log(xk);
	dp = 2*n * logx * xn * (xk)^logx;
	div = x * xk;
	dp = ifelse(div != 0, dp/div, 0); # TODO
	return(dp)
}
d2y = function(x, k=0, n=2) {
	### D()
	y.x  =  y(x, k=k, n=n);
	dy.x = dy(x, k=k, n=n);
	xn = x^n; x2 = x^2;
	xk = if(k == 0) xn else xn + k;
	div = x2 * xk^2 * y.x;
	d2p = x2*xk^2*dy.x^2 - x*(xk - n*k)*xk*y.x*dy.x + 2*n^2*xn^2*y.x^2
	d2p = ifelse(div != 0, d2p/div, d2y(1E-5, k=k, n=n)); # TODO
	return(d2p)
}
### Plot:
n = 2; k = 2;
px = (-3:3)*2/7
curve(y(x, k=k, n=n), from= -2.5, to = 2.5, ylim=c(-2,5))
line.tan(px, dx=3, p=y, dp=dy, k=k, n=n)
# TODO: px = 0;
curve(dy(x, k=k, n=n), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, n=n, col="orange")


### Ex 2:
n = 2; k = 0.9;
px = (-3:3)*2/7
curve(y(x, k=k, n=n), from= -2.5, to = 2.5, ylim=c(-2,5))
line.tan(px, dx=3, p=y, dp=dy, k=k, n=n)
# TODO: px = 0;
curve(dy(x, k=k, n=n), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, n=n, col="orange")


######################
### Generalization ###

### y = P(x)^log(P(x))

### D(y)
2*log(p)*dp / p * y

### D2(y)
# p*d2y + dp*dy
2*log(p)*dp*dy + 2*(dp)^2 / p * y + 2*log(p)*d2p*y
p / y * (dy)^2 + 2*(dp)^2 / p * y + p*d2p/dp * dy

### ODE:
dp*p^2*y*d2y - dp*p^2*(dy)^2 + dp^2*p*y*dy - d2p*p^2*y*dy - 2*(dp)^3*y^2

### Examples:
### Check: p = x + k; dp = 1;
(x+k)^2*y*d2y - (x+k)^2*(dy)^2 + (x+k)*y*dy - 2*y^2 # = 0
### p = exp(x) + k; dp = exp(x)
(e^x+k)^2*y*d2y - (e^x+k)^2*(dy)^2 + e^x*(e^x+k)*y*dy - (e^x+k)^2*y*dy - 2*e^(2*x)*y^2
### p = ln(k*x); dp = 1/x
x^2*ln(k*x)^2*y*d2y - x^2*ln(k*x)^2*(dy)^2 + x*ln(k*x)*y*dy + x*ln(k*x)^2*y*dy - 2*y^2


### Solution & Plot:
y = function(x, k=2) {
	# p = exp(x) + k;
	p = if(k == 0) exp(x) else exp(x) + k;
	val = p^log(p);
	return(val)
}
dy = function(x, k=2) {
	# 2*log(p)*dp / p * y
	x.exp = exp(x);
	p = if(k == 0) x.exp else x.exp + k;
	x.log = log(p);
	dp = 2 * x.log * x.exp * (p)^x.log;
	div = p;
	dp = ifelse(div != 0, dp/div, 0); # TODO
	return(dp)
}
d2y = function(x, k=2) {
	### D()
	y.x  =  y(x, k=k);
	dy.x = dy(x, k=k);
	x.exp = exp(x);
	xk = if(k == 0) x.exp else x.exp + k;
	div = xk^2 * y.x;
	d2p = xk^2*(dy.x)^2 - x.exp*xk*y.x*dy.x + xk^2*y.x*dy.x + 2*x.exp*x.exp*y.x^2
	d2p = ifelse(div != 0, d2p/div, 1); # TODO
	return(d2p)
}
### Plot:
k = 2;
px = (-4:2)*3/7
curve(y(x, k=k), from= -2.5, to = 2, ylim=c(-1,10))
line.tan(px, dx=3, p=y, dp=dy, k=k)
#
curve(dy(x, k=k), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, col="orange")


### Ex 2:
k = -1/2;
px = (-1:3)*3/7
curve(y(x, k=k), from= -1/2, to = 2, ylim=c(-2,10))
line.tan(px, dx=3, p=y, dp=dy, k=k)
#
curve(dy(x, k=k), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, col="orange")


###########################
###########################

###########################
### Compositions of Log ###
###########################

### y = log(exp(P1(x)) + P2(x)) + F0(x)

### y = log(exp(x^n) + k) + f

### D =>
(exp(x^n)+k)*dy - n*x^(n-1)*exp(x^n) - df*(exp(x^n)+k) # = 0
# exp(x^n) = - k*(dy - df) / (dy - n*x^(n-1) - df);

### D2 =>
(exp(x^n)+k)*d2y + n*x^(n-1)*exp(x^n)*dy +
	- n*(n-1)*x^(n-2)*exp(x^n) - n^2*x^(2*n-2)*exp(x^n) +
	- d2f*(exp(x^n)+k) - n*x^(n-1)*df*exp(x^n) # = 0
(d2y + n*x^(n-1)*dy +
		- n*(n-1)*x^(n-2) - n^2*x^(2*n-2) - n*x^(n-1)*df - d2f) * exp(x^n) +
	+ k*d2y - d2f*k # = 0
# Subst =>
k*(d2y + n*x^(n-1)*dy +
		- n*(n-1)*x^(n-2) - n^2*x^(2*n-2) - n*x^(n-1)*df - d2f) * (dy - df) +
	- (k*d2y - k*d2f)*(dy - n*x^(n-1) - df) # = 0
k*n*x^(n-1)*d2y + k*n*x^(n-1)*dy^2 - 2*k*n*x^(n-1)*df*dy +
	- k*(n*(n-1)*x^(n-2) + n^2*x^(2*n-2)) * (dy - df) +
	+ k*n*x^(n-1)*df^2 - k*n*x^(n-1)*d2f # = 0
### ODE:
x*d2y + x*dy^2 - (2*x*df + n*x^n + (n-1))*dy +
	+ (n*x^n + (n-1))*df + x*df^2 - x*d2f # = 0

# TODO: check!


#######################
#######################

### y = log(log(P(x))) + F0(x)

### y = log(log(x^2 + k)) + f

### D =>
(x^2+k)*log(x^2 + k)*dy - (x^2+k)*df*log(x^2 + k) - 2*x # = 0
# (x^2+k)*log(x^2 + k) = 2*x / (dy - df);

### D2 =>
(x^2+k)*log(x^2 + k)*d2y + 2*x*dy + 2*x*log(x^2 + k)*dy +
	- 2*x*(x^2+k)*df - 2*x*df*log(x^2 + k) - (x^2+k)*d2f*log(x^2 + k) - 2 # = 0
2*x*(x^2+k)/(dy - df) * d2y + 2*x*(x^2+k)*dy + 4*x^2/(dy - df) * dy +
	- 2*x*(x^2+k)^2*df - 4*x^2*df/(dy - df) - 2*x*(x^2+k)*d2f/(dy - df) - 2*(x^2+k) # = 0

# TODO


#######################
#######################

### y = log(P(x)) * log(log(P(x))) + F0(x)

### Example:
### y = log(x^2 + k) * log(log(x^2 + k)) + f

### D =>
(x^2+k)*log(x^2 + k)*dy - 2*x*(y - f) - 2*x*log(x^2 + k) - (x^2+k)*df*log(x^2 + k) # = 0
# log(x^2 + k) = 2*x*(y - f) / ((x^2+k)*(dy - df) - 2*x);

### D2 =>
(x^2+k)*log(x^2 + k)*d2y + 2*x*dy + 2*x/(x^2+k) * log(x^2 + k)*dy +
	- 2*x*(dy - df) - 2*(y - f) - 4*x^2/(x^2 + k) - 2*log(x^2 + k) +
	- 2*x*df - 2*x*df*log(x^2 + k) - (x^2+k)*d2f*log(x^2 + k) # = 0

### ODE:
(x^5 + 2*k*x^3 + k^2*x)*y*d2y - (x^5 + 2*k*x^3 + k^2*x)*f*d2y +
	- (x^4 - 2*x^2 + 2*k*x^2 + k^2)*y*dy +
	- (2*x^4 - f*x^4 + 2*k*x^2 + 2*f*x^2 - 2*k*f*x^2 - k^2*f)*dy +
	- (d2f*x^5 + df*x^4 + 2*k*d2f*x^3 + k^2*d2f*x - df*k^2)*y +
	+ d2f*f*x^5 + 2*df*x^4 + df*f*x^4 + 4*x^3 + 2*k*d2f*f*x^3 + 2*df*k*x^2 + k^2*d2f*f*x - df*k^2*f

### Special Case:
# f = ct;
(x^5 + 2*k*x^3 + k^2*x)*y*d2y - f*(x^5 + 2*k*x^3 + k^2*x)*d2y +
	- (x^4 - 2*x^2 + 2*k*x^2 + k^2)*y*dy +
	- (2*x^4 - f*x^4 + 2*k*x^2 + 2*f*x^2 - 2*k*f*x^2 - k^2*f)*dy + 4*x^3 # = 0
# f = 2; k = -2;
x*(x^2-2)^2*y*d2y - 2*x*(x^2-2)^2*d2y - (x^2-2)^2*y*dy + 8*(x^2 - 1)*dy + 4*x^3 # = 0

# TODO: check;

### Derivation:
pxk = toPoly.pm("x^2 + k")
pD = toPoly.pm("pxk*L*d2y + 2*x*dy + 2*x*pxk.inv * L*dy +
	- 2*x*dy + 2*x*df - 2*y + 2*f - 4*x^2*pxk.inv - 2*L +
	- 2*x*df - 2*x*df*L - pxk*d2f*L")
pL = toPoly.pm("2*x*y - 2*x*f")
pLDiv = toPoly.pm("pxk*dy - pxk*df - 2*x")
#
pR = replace.fr.pm(pD, toPoly.pm(1), pxk, x="pxk.inv", pow=1)
pR = replace.fr.pm(pR, pL, pLDiv, x="L", pow=1)
pR = replace.pm(pR, pxk, x="pxk", pow=1)
pR = simplify.spm(pR, do.gcd=TRUE)
print.dpm(pR)


#########################

### y = log(P1(x)) * log(log(P1(x))) + log(P2(x)) * log(log(P2(x))) + F0(x)
### Persistence of: log(P1) & log(P2)

### D =>
dy - LL1*dp1/p1 - LL2*dp2/p2 - dp1/p1 - dp2/p2 - df0 # = 0
p1*p2*dy - p2*dp1*LL1 - p1*dp2*LL2 - p2*dp1 - p1*dp2 - df0 # = 0

### D2 =>

### Linear System:
# TODO
# Note: one level of log() will persist;


#######################
#######################

##############
### Log(y) ###
##############

### y*log(k1*y + f) = f - k2*y

### Note:
# y*log(exp(k2)*(k1*y + f)) = f

### D =>
log(k1*y + f)*dy + (k1*y/(k1*y+f))*dy + k2*dy - df # = 0 # * y =>
(f - k2*y)*dy + (k1*y^2/(k1*y+f))*dy + k2*y*dy - df*y # = 0
(f^2 + (k1-k2)*f*y - k1*k2*y^2)*dy + k1*y^2*dy + k2*(k1*y+f)*y*dy - df*y*(k1*y+f) # = 0
### ODE:
k1*y^2*dy + k1*f*y*dy + f^2*dy - k1*df*y^2 - f*df*y # = 0

### Special Cases:
# f = x; df = 1;
k1*y^2*dy + k1*x*y*dy + x^2*dy - k1*y^2 - x*y # = 0

