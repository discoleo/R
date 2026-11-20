#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals:
## Derived from I( Trig-Fractions )
## Type: Fresnel-Integrals
##
## v.0.2c


### Examples:

# x*d2y - (n-1)*dy + n^2*x^(2*n-1)*y - n*x^n = 0;
# x*d2y - (n-1)*dy + k^2*n^2*x^(2*n-1)*y - k*n*x^n = 0;
# x^2*d2y - (n+1)*x*dy + (K^2*x^(2*n) + (n+1))*y - K*x^(n+2) = 0;


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#########################

### y = sin(x^n) * I(cos(x^n)) - cos(x^n) * I(sin(x^n))

# Check:
x  = 4/5 * pi; n = sqrt(3); # x = x + 1E-5; # Numeric D();
Ic = integrate(\(x) cos(x^n), 0, x)$value;
Is = integrate(\(x) sin(x^n), 0, x)$value;
y  = sin(x^n) * Ic - cos(x^n) * Is;

### D(y)
dy = n*x^(n-1)*cos(x^n) * Ic + n*x^(n-1)*sin(x^n) * Is;

### D2(y)
n*x^(n-1) - n^2*x^(2*n-2)*(sin(x^n) * Ic - cos(x^n) * Is) + (n-1) * dy / x;
n*x^(n-1) - n^2*x^(2*n-2)*y + (n-1) * dy / x;
d2y = (x*d2y - (n-1)*dy + n^2*x^(2*n-1)*y - n*x^n) / x;

### ODE:
x*d2y - (n-1)*dy + n^2*x^(2*n-1)*y - n*x^n # = 0


### Solution & Plot:
y.sin = function(x, n=2, lower=0) {
	sapply(x, function(upper)
		integrate(\(x) sin(x^n), lower=lower, upper=upper)$value )
}
y.cos = function(x, n=2, lower=0) {
	sapply(x, function(upper)
		integrate(\(x) cos(x^n), lower=lower, upper=upper)$value )
}
y = function(x, n=2, lower=0) {
	xn = x^n
	I.s = y.sin(x, n=n, lower=lower)
	I.c = y.cos(x, n=n, lower=lower)
	r = sin(xn)*I.c - cos(xn)*I.s;
	return(r)
}
dy = function(x, n=2, lower=0) {
	I.s = y.sin(x, n=n, lower=lower)
	I.c = y.cos(x, n=n, lower=lower)
	xn = x^n;
	r = n*xn*cos(xn) * I.c + n*xn*sin(xn) * I.s
	div = x;
	r = ifelse(div != 0, r / div, 0) # TODO: check!
	return(r)
}
d2y = function(x, n=2, lower=0) {
	y.x = y(x, n=n, lower=lower)
	dy.x = dy(x, n=n, lower=lower)
	xn = x^n;
	#
	dp  = (n-1)*dy.x - n^2*x^(2*n-1)*y.x + n*xn;
	div = x;
	dp = ifelse(div != 0, dp / div, 0); # TODO: needs correction!
	return(dp)
}
### Plot:
n = 2; lim=3.5;
x.px = c(-5:7 * 3/7)
curve(y(x, n=n), from= -3, to = 3, ylim=c(-lim, lim))
sapply(x.px, line.tan, dx=1.5, p=y, dp=dy, n=n)
#
points(x.px, dy(x.px, n=n), col = "#328632")
curve(dy(x, n=n), add=T, col="green")
sapply(x.px, line.tan, dx= 1.5, p=dy, dp=d2y, n=n, col="#409232A0")
sapply(x.px, line.tan, dx=-0.5, p=dy, dp=d2y, n=n, col="#409232A0")


#########
### Ex 2:
n = 3; lim=3.5;
x.px = c(-5:5 * 3/7)
curve(y(x, n=n), from= -2.5, to = 2.5, ylim=c(-lim, lim))
sapply(x.px, line.tan, dx=1.5, p=y, dp=dy, n=n)
#
curve(dy(x, n=n), add=T, col="green")
sapply(x.px, line.tan, dx=1.5, p=dy, dp=d2y, n=n, col="orange")


#########
### Ex 3:
n = 1.5; lim=2;
x.px = c(0:5 * 3/7)
curve(y(x, n=n), from= 0, to = 2.5, ylim=c(-lim, lim))
sapply(x.px, line.tan, dx=1.5, p=y, dp=dy, n=n)
#
curve(dy(x, n=n), add=T, col="green")
sapply(x.px, line.tan, dx=1.5, p=dy, dp=d2y, n=n, col="orange")


###################

### Extension 1:
# y = sin(k*x^n) * I(cos(k*x^n)) - cos(k*x^n) * I(sin(k*x^n))

# Check:
n = sqrt(3);  k = sqrt(2/3);
x  = 4/5 * pi; # x = x + 1E-5; # Numeric D();
Ic = integrate(\(x) cos(k*x^n), 0, x)$value;
Is = integrate(\(x) sin(k*x^n), 0, x)$value;
y  = sin(k*x^n) * Ic - cos(k*x^n) * Is;

### D(y)
dy = k*n*x^(n-1)*cos(k*x^n) * Ic + k*n*x^(n-1)*sin(k*x^n) * Is;

### D2(y)
k*n*x^(n-1) - k^2*n^2*x^(2*n-2)*(sin(k*x^n) * Ic - cos(k*x^n) * Is) + (n-1) * dy / x;
d2y = k*n*x^(n-1) - k^2*n^2*x^(2*n-2)*y + (n-1) * dy / x;

### ODE:
x*d2y - (n-1)*dy + k^2*n^2*x^(2*n-1)*y - k*n*x^n # = 0

# Reformulated:
Kn = k*n;
x*d2y - (n-1)*dy + Kn^2*x^(2*n-1)*y - Kn*x^n # = 0


###################

### Extension 2: Poly-Coefficient
# y = x * (sin(k*x^n) * I(cos(k*x^n)) - cos(k*x^n) * I(sin(k*x^n)))

# Check:
n = sqrt(3);  k = sqrt(2/3);
x  = 4/5 * pi; # x = x + 1E-5; # Numeric D();
Ic = integrate(\(x) cos(k*x^n), 0, x)$value;
Is = integrate(\(x) sin(k*x^n), 0, x)$value;
y  = x * (sin(k*x^n) * Ic - cos(k*x^n) * Is);

### D(y)
dy = k*n*x^n*cos(k*x^n) * Ic + k*n*x^n*sin(k*x^n) * Is + y/x;

### D2(y)
k*n*x^n - k^2*n^2*x^(2*n-1)*(sin(k*x^n) * Ic - cos(k*x^n) * Is) + n * (dy - y/x) / x + dy/x - y/x^2;
d2y = k*n*x^n - k^2*n^2*x^(2*n-2)*y + (n+1) * dy / x - (n+1)/x^2 * y;

### ODE:
x^2*d2y - (n+1)*x*dy + (k^2*n^2*x^(2*n) + (n+1))*y - k*n*x^(n+2) # = 0

# Reformulated:
Kn = k*n;
x^2*d2y - (n+1)*x*dy + (Kn^2*x^(2*n) + (n+1))*y - Kn*x^(n+2) # = 0


###################
###################

### Extension: Poly Argument
### y = sin(x^n + k*x) * I(cos(x^n + k*x)) - cos(x^n + k*x) * I(sin(x^n + k*x))

### D(y)
(n*x^(n-1) + k) * (cos(x^n + k*x) * Ic + sin(x^n + k*x) * Is)

### D2(y)
(n*x^(n-1) + k) - (n*x^(n-1) + k)^2 * (sin(x^n + k*x) * Ic - cos(x^n + k*x) * Is) +
	+ n*(n-1)*x^(n-2) / (n*x^(n-1) + k) * dy
n*x^(n-1) + k - (n*x^(n-1) + k)^2*y + n*(n-1)*x^(n-2) / (n*x^(n-1) + k) * dy

### ODE:
(n*x^(n-1) + k)*d2y - n*(n-1)*x^(n-2)*dy + (n*x^(n-1) + k)^3*y +
	- (n*x^(n-1) + k)^2 # = 0

### Example:
### n = 2
(2*x + k)*d2y - 2*dy + (2*x + k)^3*y + (2*x + k)^2 # = 0


### Solution & Plot:
y.sin = function(x, n=2, k=0, lower=0) {
	sapply(x, function(upper)
		integrate(function(x) sin(x^n + k*x), lower=lower, upper=upper)$value )
}
y.cos = function(x, n=2, k=0, lower=0) {
	sapply(x, function(upper)
		integrate(function(x) cos(x^n + k*x), lower=lower, upper=upper)$value )
}
y = function(x, n=2, k=0, lower=0) {
	xn = x^n + k*x;
	I.s = y.sin(x, n=n, k=k, lower=lower)
	I.c = y.cos(x, n=n, k=k, lower=lower)
	r = sin(xn)*I.c - cos(xn)*I.s;
	return(r)
}
dy = function(x, n=2, k=0, lower=0) {
	I.s = y.sin(x, n=n, k=k, lower=lower)
	I.c = y.cos(x, n=n, k=k, lower=lower)
	xn = x^n; xnk = xn + k*x; nxn = n*xn + k*x;
	r = nxn*cos(xnk) * I.c + nxn*sin(xnk) * I.s
	div = x;
	r = ifelse(div != 0, r / div, 0) # TODO: check!
	return(r)
}
d2y = function(x, n=2, k=0, lower=0) {
	y.x = y(x, n=n, k=k, lower=lower)
	dy.x = dy(x, n=n, k=k, lower=lower)
	xn = x^n; xnk = xn + k*x; nxn = n*xn + k*x;
	#
	dp  = n*(n-1)*x*xn*dy.x - nxn^3*y.x + x*nxn^2
	div = x^2 * nxn;
	dp = ifelse(div != 0, dp / div, n*k/2); # TODO: needs correction!
	return(dp)
}
### Plot:
n = 2; k = 1; lim=3.5;
x.px = c(-5:7 * 3/7)
curve(y(x, n=n, k=k), from= -3, to = 3, ylim=c(-lim, lim))
line.tan(x.px, dx=1.5, p=y, dp=dy, n=n, k=k)
#
curve(dy(x, n=n, k=k), add=T, col="green")
line.tan(x.px, dx=1.5, p=dy, dp=d2y, n=n, k=k, col="orange")


#########
### Ex 2:
n = 1/2; k = 3; lim=2;
x.px = c(0:7 * 3/7)
curve(y(x, n=n, k=k), from= 0, to = 3, ylim=c(-1, lim))
line.tan(x.px, dx=1.5, p=y, dp=dy, n=n, k=k)
# x.px = 0: needs correction;
curve(dy(x, n=n, k=k), add=T, col="green")
line.tan(x.px, dx=1.5, p=dy, dp=d2y, n=n, k=k, col="orange")


############################
############################

### y = x * I(tan(k*x^2)) + F0(x)

### D(y)
(y - f0)/x + x*tan(k*x^2) + df0
### x*dy =
(y - f0) + x^2*tan(k*x^2) + x*df0

### D2(y)
# x*d2y + dy =
dy + 2*x*tan(k*x^2) + 2*k*x^3*(tan(k*x^2)^2 + 1) + x*d2f0
# x*d2y =
2*x*tan(k*x^2) + 2*k*x^3*tan(k*x^2)^2 + 2*k*x^3 + x*d2f0
# x^2*d2y =
2*(x*dy - y + f0 - x*df0) + 2*k*(x*dy - y + f0 - x*df0)^2 + 2*k*x^4 + x^2*d2f0


### Solution & Plot:
y.I = function(x, k=1, n=2, lower=0) {
	sapply(x, function(upper)
		integrate(function(x) tan(k*x^n), lower=lower, upper=upper)$value)
}
y = function(x, b=0, k=1, n=2, lower=0) {
	# x * I(tan(k*x^2)) + F0(x)
	I.v = y.I(x, k=k, n=n, lower=lower)
	r = x*I.v + eval.pol(x, b);
	return(r)
}
dy = function(x, b=0, k=1, n=2, lower=0) {
	I.v = y.I(x, k=k, n=n, lower=lower)
	# (y - f0) + x^2*tan(k*x^2) + x*df0
	# y - f0 = x*I.v;
	r = x*I.v + x^2 * tan(k*x^2) + deriv.pol(x, b, x.mult=1); # x*df0
	div = x;
	r = ifelse(div != 0, r / div,
		dy(x + 1E-3, b=b, k=k, n=n, lower=lower)) # TODO
	return(r)
}
d2y = function(x, b=0, k=1, n=2, lower=0) {
	y.x = y(x, b=b, k=k, n=n, lower=lower)
	dy.x = dy(x, b=b, k=k, n=n, lower=lower)
	# 2*(x*dy - y + f0 - x*df0) + 2*k*(x*dy - y + f0 - x*df0)^2 + 2*k*x^4 + x^2*d2f0
	f0 = eval.pol(x, b); xdf0 = deriv.pol(x, b, dn=1, x.mult=1);
	d2f0 = deriv.pol(x, b, dn=2, x.mult=0);
	x2 = x^2;
	Tx = x*dy.x - y.x + f0 - xdf0;
	dp = 2*Tx + 2*k*Tx^2 + 2*k*x2*x2 + x2*d2f0;
	div = x2;
	dp = ifelse(div != 0, dp / div, d2y(x + 1E-3, b=b, n=n, lower=lower)); # TODO
	return(dp)
}
### Plot:
b = c(0, 1); k = 1;
x.px = c(-2:2 * 2/7); xlim = sqrt(pi/2)/k - 1E-2;
curve(y(x, b=b, k=k), from= -xlim, to = xlim, ylim=c(-2, 3.5))
line.tan(x.px*2.08, dx=1.6, p=y, dp=dy, b=b, k=k)
#
curve(dy(x, b=b, k=k), add=T, col="green")
line.tan(x.px*1.6, dx=1.5, p=dy, dp=d2y, b=b, k=k, col="orange")


### Ex 2:
b = c(0, 1); k = 2;
x.px = c(-2:2 * 1/7); xlim = sqrt(pi/2)/k - 1E-2;
curve(y(x, b=b, k=k), from= -xlim, to = xlim, ylim=c(-1, 2))
line.tan(x.px*2.08, dx=1.6, p=y, dp=dy, b=b, k=k)
#
curve(dy(x, b=b, k=k), add=T, col="green")
line.tan(x.px*1.7, dx=1.5, p=dy, dp=d2y, b=b, k=k, col="orange")

