
########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Logarithms
###
### draft v.0.3f


### ODEs Derived from Logarithms

### Note:
# Lambert W examples:
# - should be in a different file;

###############

###############
### History ###
###############

### draft v.0.3f:
# - derived from: [TODO full derivation]
#   y = log(exp(P1(x)) + P2(x)) + F0(x);
#   y = log(log(P(x))) + F0(x);
### draft v.0.3e:
# - derived from:
#   y*log(x+k) = I(log(x+k) / x) + F0(x);
### draft v.0.3d:
# - cleanup;
### draft v.0.3c - v.0.3c-ex2:
# - Mixed Log-Exp:
#   y = log(P1(x)) * exp(P2(x)) + F(x);
# - more examples; [v.0.3c-ex2]
### draft v.0.3a - v.0.3b:
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

### Helper functions

library(pracma)
# - may be needed to solve various equations;

# include: DE.ODE.Helper.R;
source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#########################
#########################


### Section A: Linear ODEs

### Derived from:
### y = F1(x)*log(P1(x)) + F2(x)*log(P2(x))

### Examples:
### y = (x^m + c)*log(x^n + a) + (x^m - c)*log(x^n + b)
# - simplifies massively & eliminates y from D2(y);
# - does NOT simplify when: F1 - F2 != constant;

### D(y)
# [not run]
m*x^(m-1)*log(x^n + a) + m*x^(m-1)*log(x^n + b) +
	+ n*x^(n-1)*(x^m + c) / (x^n + a) +
	+ n*x^(n-1)*(x^m - c) / (x^n + b)
m*x^(m-1)*log(x^n + a) + m*x^(m-1)*log(x^n + b) +
	+ n*x^(n-1)*(2*x^(m+n) + (a+b)*x^m - c*(a-b)) / ((x^n + a) * (x^n + b))
# let: T(x) = x^(n-1)*(2*x^(m+n) + (a+b)*x^m - c*(a-b)) / ((x^n + a) * (x^n + b))
# Note: + n * T(x);
m*x^(m-1)*log(x^n + a) + m*x^(m-1)*log(x^n + b) + n*T
### Solve linearly for log(P(x))
log(x^n + a) =
	((x^m - c)*dy - m*x^(m-1)*y - n*T*(x^m - c)) / - (2*c*m*x^(m-1));
log(x^n + b) =
	((x^m + c)*dy - m*x^(m-1)*y - n*T*(x^m + c)) /   (2*c*m*x^(m-1));
log(x^n + a) + log(x^n + b) = (dy - n*T) / (m*x^(m-1));

###
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

### y = log(P1(x)) * exp(P2(x)) + F(x);

### y = log(x^2 + b0) * exp(x^2)
### D =>
dy - 2*x*y - 2*x*exp(x^2)/(x^2 + b0) # = 0
(x^2 + b0)*dy - 2*x*(x^2 + b0)*y - 2*x*exp(x^2) # = 0

### D2 =>
(x^2 + b0)*d2y - 2*x*(x^2 + b0 - 1)*dy - (6*x^2 + 2*b0)*y +
	- 2*(2*x^2 + 1)*exp(x^2) # = 0
x*(x^2 + b0)*d2y - 2*x^2*(x^2 + b0 - 1)*dy - x*(6*x^2 + 2*b0)*y +
	- (2*x^2 + 1)*((x^2 + b0)*dy - 2*x*(x^2 + b0)*y) # = 0
### ODE:
x*(x^2 + b0)*d2y - (4*x^4 + (4*b0 - 1)*x^2 + b0)*dy +
	+ 4*x^3*(x^2 + (b0 - 1))*y # = 0


### Solution & Plot:
y = function(x, b=1) {
	x2 = x^2;
	val = log(x2 + b[1]) * exp(x2);
	return(val)
}
dy = function(x, b=1) {
	yx = y(x, b=b);
	x2 = x^2; x2b = x2 + b[1];
	dyx = 2*x*(x2b*yx + exp(x2));
	div = x2b;
	dyx = ifelse(div != 0, dyx/div, 1); # TODO: check;
	return(dyx)
}
d2y = function(x, b=1) {
	yx  =  y(x, b=b);
	dyx = dy(x, b=b);
	#
	x2 = x^2; x2b = x2 + b[1];
	d2p = (4*x2*x2b - x2 + b[1])*dyx - 4*x2*x*(x2b - 1)*yx;
	div = x*x2b;
	d2p = ifelse(div != 0, d2p/div, 1); # TODO: check;
	return(d2p);
}
### Plot:
b = 1
px = (1:3)*2/7 + 1; px = c(0, -px, px);
curve(y(x, b=b), from = -2, to = 2, ylim=c(-30, 80))
line.tan(px, dx=3, p=y, dp=dy, b=b)
# global minimum:
curve(dy(x, b=b), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, b=b, col="orange")

### Ex 2:
b = -1
px = (1:4)*2/7 + 1;
curve(y(x, b=b), from = 1, to = 2.5, ylim=c(-10, 200))
line.tan(px, dx=3, p=y, dp=dy, b=b)
# global minimum:
curve(dy(x, b=b), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, b=b, col="orange")


#########################
#########################

#######################
### Section B:      ### 
### Non-Linear ODEs ###
#######################

###################
### Logarithmic ###
###################

### y^n = log(f(x)) * log(g(x))

### D =>
n*y^(n-1)*dy = df * log(g)/f + dg * log(f)/g # * f*g
n*f*g*y^(n-1)*dy = g*df * log(g) + f*dg * log(f)

### D2 =>
n*f*g*y^(n-1)*d2y + n*(n-1)*f*g*y^(n-2)*dy^2 + n*df*g*y^(n-1)*dy + n*f*dg*y^(n-1)*dy =
	(g*d2f + df*dg) * log(g) + (f*d2g + df*dg) * log(f) + df + dg
### Solve Linear =>
# log(f) = ...
# log(g) = ...

### Special Cases:

### Order: n = 1
f*g*dy = g*df * log(g) + f*dg * log(f)
f*g*d2y + df*g*dy + f*dg*dy - df - dg =
	(g*d2f + df*dg) * log(g) + (f*d2g + df*dg) * log(f)

### Solve Liniar =>
(f*dg*(g*d2f + df*dg) - g*df*(f*d2g + df*dg))*log(g) = ...
(f*g*dg*d2f + f*df*dg^2 - f*g*df*d2g - g*df^2*dg)*log(g) = ...
#
log(f) = ...


### Order: n = 1/2
1/2*f*g*y^(-1/2)*dy = g*df * log(g) + f*dg * log(f)
1/2*f*g*y^(-1/2)*d2y - 1/4*f*g*y^(-3/2)*dy^2 + 1/2*df*g*y^(-1/2)*dy + 1/2*f*dg*y^(-1/2)*dy -df - dg =
	(g*d2f + df*dg) * log(g) + (f*d2g + df*dg) * log(f)


### Examples:

### y = log(x + k) * log(x - k)
dy = log(x-k)/(x+k) + log(x+k)/(x-k)
(x^2 - k^2)*dy = (x-k)*log(x-k) + (x+k)*log(x+k)
### D2 =>
(x^2 - k^2)*d2y + 2*x*dy - 2 = log(x-k) + log(x+k)
### Solve Liniar =>
2*k*log(x-k) = (x+k)*((x^2 - k^2)*d2y + (x + k)*dy - 2)
2*k*log(x+k) = -(x-k)*((x^2 - k^2)*d2y + (x - k)*dy - 2)

### ODE:
(x^2 - k^2)*((x^2 - k^2)*d2y + (x + k)*dy - 2) * ((x^2 - k^2)*d2y + (x - k)*dy - 2) + 4*k^2*y = 0

### Ex 2:
### y^(1/2) = log(x + k) * log(x - k)
1/2*y^(-1/2)*dy = log(x-k)/(x+k) + log(x+k)/(x-k)
(x^2 - k^2)*y^(-1/2)*dy = 2*(x-k)*log(x-k) + 2*(x+k)*log(x+k)
### D2 =>
(x^2 - k^2)*y^(-1/2)*d2y + 2*x*y^(-1/2)*dy - 1/2*(x^2 - k^2)*y^(-3/2)*dy^2 - 4 = 2*log(x-k) + 2*log(x+k)
### Solve Liniar =>
4*k*log(x-k) =  (x+k)*y^(-1/2)*((x^2 - k^2)*d2y + (x + k)*dy - 1/2*(x^2 - k^2)*dy^2/y - 4)
4*k*log(x+k) = -(x-k)*y^(-1/2)*((x^2 - k^2)*d2y + (x - k)*dy - 1/2*(x^2 - k^2)*dy^2/y - 4)

### ODE:
(x^2 - k^2)*(...)*(...) + 4*k^2*y^(3/2) = 0


### Solution & Plot
y = function(x, k=1, n=1, v.dy, v.d2y) {
	if(missing(v.dy)) v.dy = dy(x, k=k, n=n)
	if(missing(v.d2y)) v.d2y = d2y(x, k=k, n=n, v.dy=v.dy)
	x2 = x^2 - k^2
	y = - x2*(x2*v.d2y + (x + k)*v.dy - 2)*(x2*v.d2y + (x - k)*v.dy - 2)
	return(y / 4 / k^2)
}
dy = function(x, k=1, n=1) {
	dp = log(x-k)/(x+k) + log(x+k)/(x-k);
	return(dp)
}
d2y = function(x, k=1, n=1, v.dy) {
	if(missing(v.dy)) v.dy = dy(x, k=k, n=n)
	dp =log(x-k) + log(x+k) + 2 - 2*x*v.dy;
	dp = dp / (x^2 - k^2)
	return(dp)
}
### Plot:
k = 1; n = 1;
px = c(3/5 + (1:3)*3/5);
curve(y(x, k=k, n=n), from= 1.01, to = 3, ylim=c(-3, 3))
# global "minimum" / horn;
line.tan(px, dx=3, p=y, dp=dy, k=k, n=n)
#
curve(dy(x, k=k, n=n), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, n=n, col="orange")


### Ex 2:
k = 3; n = 1;
px = c(3 + (1:3)*3/5);
curve(y(x, k=k, n=n), from= 3.01, to = 6, ylim=c(-3, 3))
# global "minimum" / horn;
line.tan(px, dx=3, p=y, dp=dy, k=k, n=n)
#
curve(dy(x, k=k, n=n), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, n=n, col="orange")


#######################
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


##################

##################
### Log: Other ###
##################

### y = (x + k1)^(x + k2) + F0(x)

# "a problem which stumped Photomath!"

### D(y)
(log(x+k1) + (x+k2)/(x+k1))*(y - f0) + df0
# (x+k1)*dy =
((x+k1)*log(x+k1) + (x+k2))*(y - f0) + (x+k1)*df0

# (x+k1)*log(x+k1) = (x+k1)*(dy - df0) / (y - f0) - (x + k2);

### D2(y)
# (x+k1)*d2y + dy =
((x+k1)*log(x+k1) + (x+k2))*(dy - df0) +
	+ (log(x+k1) + 2)*(y - f0) +
	+ (x+k1)*d2f0 + df0
(x+k1)*(dy - df0)^2 / (y - f0) +
	+ (log(x+k1) + 2)*(y - f0) +
	+ (x+k1)*d2f0 + df0
# (x+k1)*(y - f0)*d2y + (y - f0)*dy =
(x+k1)*(dy - df0)^2 +
	+ log(x+k1)*(y - f0)^2 + 2*(y - f0)^2 +
	+ (x+k1)*d2f0*(y - f0) + df0*(y - f0)
# (x+k1)^2*(y - f0)*d2y + (x+k1)*(y - f0)*dy =
(x+k1)^2*(dy - df0)^2 +
	+ (x+k1)*(y - f0)*(dy - df0) - (x + k2)*(y - f0)^2 + 2*(x+k1)*(y - f0)^2 +
	+ (x+k1)^2*d2f0*(y - f0) + (x+k1)*df0*(y - f0)

### ODE:
(x+k1)^2*(y - f0)*d2y # =
(x+k1)^2*(dy - df0)^2 + (x+k1)*(y - f0)^2 + (k1 - k2)*(y - f0)^2 +
	+ (x+k1)^2*d2f0*(y - f0);

### Examples:
### k1 == k2
(x+k1)*(y - f0)*d2y # =
(x+k1)*(dy - df0)^2 + (y - f0)^2 + (x+k1)*d2f0*(y - f0);


### Solution & Plot:
y = function(x, k=c(0,0), b=1, asList=FALSE) {
	f0 = eval.pol(x, b=b);
	y.x = (x + k[1])^(x + k[2]) + f0;
	if(asList) return(list(y=y.x, f0=f0));
	return(y.x)
}
dy = function(x, k=c(0,0), b=1, asList=FALSE) {
	y.all = y(x, k=k, b=b, asList=TRUE);
	y.x = y.all$y; f0 = y.all$f0;
	df0 = deriv.pol(x, b=b, dn=1);
	xk = x + k[1]; x.log = log(xk);
	dp = (xk*x.log + (x+k[2]))*(y.x - f0) + xk*df0
	div = xk;
	dp = ifelse(div != 0, dp/div, 0); # TODO
	if(asList) return(c(y.all, dy=dp, df0=df0));
	return(dp)
}
d2y = function(x, k=c(0,0), b=1) {
	y.all = dy(x, k=k, b=b, asList=TRUE);
	y.x = y.all$y; dy.x = y.all$dy;
	f0 = y.all$f0; df0 = y.all$df0;
	d2f0 = deriv.pol(x, b=b, dn=2);
	Tx = y.x - f0; xk = x + k[1]; Txk = xk*Tx;
	dTx = dy.x - df0;
	div = xk * Txk;
	d2p = xk^2*dTx^2 +
		+ (k[1] - k[2])*Tx^2 + Txk*Tx + xk*d2f0*Txk;
	d2p = ifelse(div != 0, d2p/div, 1); # TODO
	return(d2p)
}
### Plot:
k = c(0,0); b = c(0, 1)
px = (0:4)*4/7 + 0.1;
curve(y(x, k=k, b=b), from = 1E-3, to = 3, ylim=c(-2,8))
line.tan(px, dx=3, p=y, dp=dy, k=k, b=b)
#
curve(dy(x, k=k, b=b), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, b=b, col="orange")


### Ex 2:
k = c(0,0); b = c(-1, 1)
px = (0:4)*4/7 + 0.1;
curve(y(x, k=k, b=b), from = 1E-3, to = 3, ylim=c(-2,8))
line.tan(px, dx=3, p=y, dp=dy, k=k, b=b)
#
curve(dy(x, k=k, b=b), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, b=b, col="orange")


### Ex 3:
k = c(1,-1); b = c(-1, 1)
px = (0:4)*4/7 + 0.1; dx = c(1.5, 2, 2, 3, 3)
curve(y(x, k=k, b=b), from = 1E-3, to = 3, ylim=c(-0.5,8))
line.tan(px, dx=dx, p=y, dp=dy, k=k, b=b)
#
curve(dy(x, k=k, b=b), add=T, col="green")
line.tan(px, dx=dx, p=dy, dp=d2y, k=k, b=b, col="orange")


#######################
#######################

### y * (x + k1)^(x + k2) = F0(x)

# "another problem which stumped Photomath!" ;-)

### D(y)
f0*dy/y + (log(x+k1) + (x+k2)/(x+k1))*y*f0/y - df0 # = 0
f0*dy + (log(x+k1) + (x+k2)/(x+k1))*f0*y - df0*y # = 0
(x+k1)*f0*dy + ((x+k1)*log(x+k1) + (x+k2))*f0*y - (x+k1)*df0*y # = 0

# (x+k1)*log(x+k1)*f0*y = -(x+k1)*(f0*dy - df0*y) - (x + k2)*f0*y;

### D2(y)
(x+k1)*f0*d2y + f0*dy - df0*y - (x+k1)*d2f0*y +
	+ ((x+k1)*log(x+k1) + (x+k2))*(f0*dy + df0*y) +
	+ (log(x+k1) + 2)*f0*y # = 0
(x+k1)*f0*d2y + f0*dy - df0*y - (x+k1)*d2f0*y +
	- (x+k1)*(f0*dy - df0*y)*(dy/y + df0/f0) +
	+ log(x+k1)*f0*y + 2*f0*y # = 0
(x+k1)^2*f0^2*y*d2y +
	- (x+k1)^2*(f0*dy + df0*y)*(f0*dy - df0*y) +
	- (x+k1)^2*f0*d2f0*y^2 + (x+k1)*f0^2*y^2 + (k1 - k2)*f0^2*y^2 # = 0

### Examples:
### k1 == k2
(x+k1)*f0^2*y*d2y +
	- (x+k1)*(f0*dy + df0*y)*(f0*dy - df0*y) +
	- (x+k1)*f0*d2f0*y^2 + f0^2*y^2 # = 0


### Solution & Plot:
y = function(x, k=c(0,0), b=1, asList=FALSE) {
	f0 = eval.pol(x, b=b);
	div = (x + k[1])^(x + k[2]);
	y.x = ifelse(div != 0, f0/div, 0); # TODO
	if(asList) return(list(y=y.x, f0=f0));
	return(y.x)
}
dy = function(x, k=c(0,0), b=1, asList=FALSE) {
	y.all = y(x, k=k, b=b, asList=TRUE);
	y.x = y.all$y; f0 = y.all$f0;
	df0 = deriv.pol(x, b=b, dn=1);
	xk = x + k[1]; x.log = log(xk);
	dp = -((xk*x.log + x + k[2])*f0 - xk*df0)*y.x;
	div = xk * f0;
	dp = ifelse(div != 0, dp/div,
		dy(x + 1E-3, k=k, b=b)); # TODO
	if(asList) return(c(y.all, dy=dp, df0=df0));
	return(dp)
}
d2y = function(x, k=c(0,0), b=1) {
	y.all = dy(x, k=k, b=b, asList=TRUE);
	y.x = y.all$y; dy.x = y.all$dy;
	f0 = y.all$f0; df0 = y.all$df0;
	d2f0 = deriv.pol(x, b=b, dn=2);
	Tx  = f0*dy.x - df0*y.x;
	Txp = f0*dy.x + df0*y.x;
	xk = x + k[1]; xk2 = xk*xk; y2 = y.x*y.x;
	div = xk2*f0^2*y.x;
	d2p = xk2*Txp*Tx +
		+ xk2*f0*d2f0*y2 - xk*f0^2*y2 - (k[1] - k[2])*f0^2*y2;
	d2p = ifelse(div != 0, d2p/div,
		d2y(x + 1E-3, k=k, b=b)); # TODO
	return(d2p)
}
### Plot:
k = c(0,0); b = c(0, 1)
px = (0:4)*4/7 + 0.1;
curve(y(x, k=k, b=b), from = 1E-3, to = 3, ylim=c(-1.5,2))
line.tan(px, dx=3, p=y, dp=dy, k=k, b=b)
#
curve(dy(x, k=k, b=b), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, b=b, col="orange")


### Ex 2:
k = c(1, -2); b = c(0, 1)
px = (0:4)*4/7 + 0.1;
curve(y(x, k=k, b=b), from = 0, to = 3, ylim=c(-1.5, 3))
line.tan(px, dx=3, p=y, dp=dy, k=k, b=b)
#
curve(dy(x, k=k, b=b), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, b=b, col="orange")


#######################
#######################

### y = log(exp(P1(x)) + P2(x)) + F0(x)

### y = log(exp(x^n) + k) + f

### D =>
(exp(x^n)+k)*dy - n*x^(n-1)*exp(x^n) - df*(exp(x^n)+k) # = 0
# exp(x^n) = - k*(dy - df) / (dy - n*x^(n-1) - df);

### D2 =>

# TODO


#######################
#######################

### y = log(log(P(x))) + F0(x)

### y = log(log(x^2 + k)) + f

### D =>
(x^2+k)*log(x^2 + k)*dy - (x^2+k)*df*log(x^2 + k) - 2*x # = 0
# (x^2+k)*log(x^2 + k) = 2*x / (dy - df);

### D2 =>

# TODO


#######################
#######################

#################
### Integrals ###
#################

### y = x * I(1/log(x + k)) dx + F0(x)

### D(y)
I + x/log(x + k) + df0
# x*dy =
y - f0 + x^2/log(x + k) + x*df0

### D2(y)
# x*d2y + dy =
dy + x*d2f0 + 2*x/log(x + k) - x^2/((x+k)*log(x+k)^2)
dy + x*d2f0 + 2*(x*dy - y + f0 - x*df0)/x - (x*dy - y + f0 - x*df0)^2/(x^2*(x+k))
# x^3*(x+k)*d2y =
2*x*(x+k)*(x*dy - y + f0 - x*df0) - (x*dy - y + f0 - x*df0)^2 + x^3*(x+k)*d2f0;

### Examples:
### k = 0
x^4*d2y - 2*x^2*(x*dy - y + f0 - x*df0) + (x*dy - y + f0 - x*df0)^2 - x^4*d2f0 # = 0
### f0 = x
x^4*d2y + x^2*dy^2 - 2*x*y*dy - 2*x^3*dy + y^2 + 2*x^2*y # = 0


### Solution & Plot:
y.I = function(x, k=0, n=1, lower=1+1E-3) {
	sapply(x, function(upper)
		integrate(function(x) 1/log(x^n + k), lower=lower, upper=upper)$value)
}
y = function(x, b=0, k=0, n=1, lower=1+1E-3) {
	# x * I(1/log(x^n + k)) + F0(x)
	I.v = y.I(x, k=k, n=n, lower=lower)
	r = x*I.v + eval.pol(x, b);
	return(r)
}
dy = function(x, b=0, k=0, n=1, lower=1+1E-3) {
	I.v = y.I(x, k=k, n=n, lower=lower)
	# (y - f0) + x^2/log(x^n + k) + x*df0
	# y - f0 = x*I.v;
	r = x*I.v + x^2 / log(x^n + k) + deriv.pol(x, b, x.mult=1); # x*df0
	div = x;
	r = ifelse(div != 0, r / div,
		dy(x + 1E-3, b=b, k=k, n=n, lower=lower)) # TODO
	return(r)
}
d2y = function(x, b=0, k=0, n=1, lower=1+1E-3) {
	y.x = y(x, b=b, k=k, n=n, lower=lower)
	dy.x = dy(x, b=b, k=k, n=n, lower=lower)
	# 2*x*(x+k)*(x*dy - y + f0 - x*df0) - (x*dy - y + f0 - x*df0)^2 + x^3*(x+k)*d2f0
	f0 = eval.pol(x, b); xdf0 = deriv.pol(x, b, dn=1, x.mult=1);
	d2f0 = deriv.pol(x, b, dn=2, x.mult=0);
	x3 = x^3; x3k = x3 * (x+k);
	Tx = x*dy.x - y.x + f0 - xdf0;
	dp = 2*x*(x+k)*Tx - Tx^2 + x3k*d2f0;
	div = x3k;
	dp = ifelse(div != 0, dp / div,
		d2y(x + 1E-3, b=b, k=k, n=n, lower=lower)); # TODO
	return(dp)
}
### Plot:
b = c(0, 1); k = 0;
x.px = c(0.01, 1:3 * 2/13) + (1-k); xlim = c(1-k + 1E-3, 1-k + 3);
curve(y(x, b=b, k=k), from= xlim[1], to = xlim[2], ylim=c(0, 30))
line.tan(x.px, dx=1.6, p=y, dp=dy, b=b, k=k)
#
curve(dy(x, b=b, k=k), add=T, col="green")
line.tan(x.px, dx=1.5, p=dy, dp=d2y, b=b, k=k, col="orange")


### Ex 2:
b = c(0, 1); k = 2;
x.px = c(0.01, 1:4 * 3/7) + (1-k); xlim = c(1-k + 1E-3, 1-k + 3);
curve(y(x, b=b, k=k), from= xlim[1], to = xlim[2], ylim=c(-5, 10))
line.tan(x.px, dx=1.6, p=y, dp=dy, b=b, k=k)
#
curve(dy(x, b=b, k=k), add=T, col="green")
line.tan(x.px*1.5 + 0.8, dx=1.5, p=dy, dp=d2y, b=b, k=k, col="orange")


#########################
#########################

#########################
### Inverse Integrals ###
#########################

### y * I(1/log(x + k)) dx = F0(x)

### D(y)
I*dy + 1/log(x+k) * y - df0 # = 0
f0*dy + 1/log(x+k) * y^2 - df0*y # = 0

### D2(y)
f0*d2y + 2/log(x+k) * y*dy - y^2/((x+k)*log(x+k)^2) - d2f0*y # = 0
(x+k)*f0*d2y + 2*(x+k)/log(x+k) * y*dy - y^2/log(x+k)^2 - (x+k)*d2f0*y # = 0
# 1/log(x+k) = - (f0*dy - df0*y) / y^2;
(x+k)*f0*d2y - 2*(x+k)*(f0*dy - df0*y)/y * dy +
	- (f0*dy - df0*y)^2/y^2 - (x+k)*d2f0*y # = 0
(x+k)*f0*y^2*d2y - 2*(x+k)*(f0*dy - df0*y)*y*dy +
	- (f0*dy - df0*y)^2 - (x+k)*d2f0*y^3 # = 0

### ODE:
(x+k)*f0*y^2*d2y - 2*(x+k)*(f0*dy - df0*y)*y*dy - (f0*dy - df0*y)^2 - (x+k)*d2f0*y^3 # = 0

### Examples:
### k = 0
x*f0*y^2*d2y - 2*x*(f0*dy - df0*y)*y*dy - (f0*dy - df0*y)^2 - x*d2f0*y^3 # = 0
### f0 = 1
(x+k)*y^2*d2y - 2*(x+k)*y*dy^2 - dy^2 # = 0


### Solution & Plot:
y.I = function(x, k=0, n=1, lower=1+1E-3) {
	sapply(x, function(upper)
		integrate(function(x) 1/log(x^n + k), lower=lower, upper=upper)$value)
}
y = function(x, b=0, k=0, n=1, lower=1+1E-3) {
	# y * I(1/log(x + k)) = F0(x)
	I.v = y.I(x, k=k, n=n, lower=lower)
	r = eval.pol(x, b);
	div = I.v;
	r = ifelse(I.v != 0, r / div, Inf);
	return(r)
}
dy = function(x, b=0, k=0, n=1, lower=1+1E-3) {
	I.v = y.I(x, k=k, n=n, lower=lower)
	# f0*log(x+k)*dy = - y^2 + log(x+k)*df0*y # = 0
	logx = log(x^n + k); f0 = eval.pol(x, b);
	r = - f0 + I.v * logx * deriv.pol(x, b);
	div = I.v^2 * logx;
	r = ifelse(div != 0, r / div,
		dy(x + 1E-3, b=b, k=k, n=n, lower=lower)) # TODO
	return(r)
}
d2y = function(x, b=0, k=0, n=1, lower=1+1E-3) {
	y.x = y(x, b=b, k=k, n=n, lower=lower)
	dy.x = dy(x, b=b, k=k, n=n, lower=lower)
	# (x+k)*f0*y^2*d2y - 2*(x+k)*(f0*dy - df0*y)*y*dy + (f0*dy - df0*y)^2 - (x+k)*d2f0*y^3
	f0 = eval.pol(x, b); df0 = deriv.pol(x, b, dn=1);
	d2f0 = deriv.pol(x, b, dn=2, x.mult=0);
	xk = x + k; xkf = xk * f0;
	Tx = f0*dy.x - df0*y.x;
	dp = 2*xk*Tx*y.x*dy.x + Tx^2 + xk*d2f0*y.x^3;
	div = xkf * y.x^2;
	dp = ifelse(div != 0, dp / div,
		d2y(x + 1E-3, b=b, k=k, n=n, lower=lower)); # TODO
	return(dp)
}
### Plot:
b = c(0, 1); k = 0;
x.px = c(0.1, 1:3 * 2/13) + (1-k); xlim = c(1-k + 1E-3, 1-k + 3);
#
curve(y(x, b=b, k=k), from= xlim[1], to = xlim[2], ylim=c(-0.5, 0.45))
line.tan(x.px, dx=1.6, p=y, dp=dy, b=b, k=k)
#
curve(dy(x, b=b, k=k), add=T, col="green")
line.tan(x.px, dx=1.5, p=dy, dp=d2y, b=b, k=k, col="orange")


### Ex 2:
b = c(0, 1); k = 2;
x.px = c(0.1, 1:3 * 2/13) + (1-k); xlim = c(1-k + 1E-3, 1-k + 1);
#
curve(y(x, b=b, k=k), from= xlim[1], to = xlim[2], ylim=c(-0.5, 1))
line.tan(x.px, dx=1.6, p=y, dp=dy, b=b, k=k)
#
curve(dy(x, b=b, k=k), add=T, col="green")
line.tan(x.px, dx=1.5, p=dy, dp=d2y, b=b, k=k, col="orange")

### Ex 2 (cont):
b = c(0, 1); k = 2;
x.px = c(0.1, 1:4 * 2/5) + (1-k); xlim = c(1-k + 1E-3, 1-k + 3);
#
curve(y(x, b=b, k=k), from= xlim[1], to = xlim[2], ylim=c(-5, 1))
line.tan(x.px, dx=1.6, p=y, dp=dy, b=b, k=k)
#
curve(dy(x, b=b, k=k), add=T, col="green")
line.tan(x.px, dx=0.75, p=dy, dp=d2y, b=b, k=k, col="orange")


#####################
#####################

### y*log(x+k) = I(log(x+k) / x) + F0(x)

### D =>
x*(x+k)*log(x+k)*dy + x*y - (x+k)*log(x+k) - x*(x+k)*df # = 0
# (x+k)*log(x+k)*(x*dy - 1) = - x*y + x*(x+k)*df;

### D2 =>
x*(x+k)*log(x+k)*d2y + 2*x*dy + ((2*x+k)*log(x+k))*dy + y +
	- log(x+k) - 1 - (2*x+k)*df - x*(x+k)*d2f # = 0
(x*(x+k)*d2y + (2*x+k)*dy - 1)*log(x+k) +
	+ 2*x*dy + y - (2*x+k)*df - x*(x+k)*d2f - 1 # = 0
(x^2*(x+k)*d2y + x*(2*x+k)*dy - x)*(y - (x+k)*df) +
	- (x+k)*(x*dy - 1)*(2*x*dy + y - (2*x+k)*df - x*(x+k)*d2f - 1) # = 0
### ODE:
x^2*(x+k)*y*d2y - x^2*(x+k)^2*df*d2y +
	- 2*x^2*(x+k)*dy^2 + x^2*y*dy +
	+ x^2*(x+k)^2*d2f*dy + 3*x*(x+k)*dy + k*y +
	- (x+k)*((x+k)*df + x*(x+k)*d2f + 1) # = 0

### Special Case:
# f = ct
x^2*(x+k)*y*d2y - 2*x^2*(x+k)*dy^2 + x^2*y*dy +
	+ 3*x*(x+k)*dy + k*y - (x+k) # = 0


### Solution & Plot:
y.I = function(x, k=1, n=1, lower=1) {
	sapply(x, function(upper)
		integrate(function(x) log(x^n + k) / x, lower=lower, upper=upper)$value)
}
y = function(x, k=1, n=1, f=NULL, lower=1) {
	yx = y.I(x, k=k, n=n, lower=lower);
	if( ! is.null(f)) {
		f0 = sapply(x, function(x) eval.pm(f, x));
		yx = yx + f0;
	}
	div = log(x^n + k);
	r = ifelse(div != 0, yx / div, Inf); # TODO
	return(r)
}
dy = function(x, k=1, n=1, f=NULL, lower=1) {
	# x*(x+k)*log(x+k)*dy + x*y - (x+k)*log(x+k) - x*(x+k)*df
	yx = y(x, k=k, n=n, f=f, lower=lower);
	xk = x+k; logxk = xk*log(xk);
	r  = logxk - x*yx;
	if( ! is.null(f)) {
		df = dp.pm(f, xn="x");
		dfx = sapply(x, function(x) eval.pm(df, x));
		r = r + x*(x+k)*dfx;
	}
	div = x*logxk;
	r = ifelse(div != 0, r / div, 0) # TODO
	return(r)
}
d2y = function(x, k=1, n=1, f=NULL, lower=1) {
	yx = y(x, k=k, n=n, f=f, lower=lower);
	dyx = dy(x, k=k, n=n, f=f, lower=lower);
	xk = x + k; x2 = x^2;
	dp = 2*x2*xk*dyx^2 - x2*yx*dyx +
		- 3*x*xk*dyx + xk - k*yx;
	div = x2*xk*yx;
	if( ! is.null(f)) {
		df = dp.pm(f, xn="x");
		dfx = sapply(x, function(x) eval.pm(df, x));
		d2f = dp.pm(df, xn="x");
		d2fx = sapply(x, function(x) eval.pm(d2f, x));
		dp = dp + xk*(xk*dfx + x*xk*d2fx) - x2*xk^2*d2fx*dyx;
		div = div - x2*xk^2*dfx;
	}
	dp = ifelse(div != 0, dp / div, 0); # TODO
	return(dp)
}
### Plot:
k = 1;
f = toPoly.pm("x^2 - 3*x + 4")
x.px = c(0.1, 1:3 * 5/13) + 1; xlim = c(1, 3);
#
curve(y(x, k=k, f=f), from= xlim[1], to = xlim[2], ylim=c(-1, 3.6))
line.tan(x.px, dx=1.6, p=y, dp=dy, k=k, f=f)
# global minimum
curve(dy(x, k=k, f=f), add=T, col="green")
line.tan(x.px, dx=1.5, p=dy, dp=d2y, k=k, f=f, col="orange")


### Example 2:
k = 3;
f = toPoly.pm("x^3 - 3*x - 4")
x.px = c(0.1, 1:3 * 5/13) + 1; xlim = c(1, 3);
#
curve(y(x, k=k, f=f), from= xlim[1], to = xlim[2], ylim=c(-5, 10))
line.tan(x.px, dx=1.6, p=y, dp=dy, k=k, f=f)
#
curve(dy(x, k=k, f=f), add=T, col="green")
line.tan(x.px, dx=1.5, p=dy, dp=d2y, k=k, f=f, col="orange")

