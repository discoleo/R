########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs - Logarithms
##
## draft v.0.3t


### ODEs Derived from Logarithms

### Note:
# Lambert W examples:
# - should be in a different file;

###############
### History ###

# - Moved to file:
#   DE.ODE.Log.History.md;


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

