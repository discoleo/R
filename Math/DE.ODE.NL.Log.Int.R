########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Integral( LOG )
##
## draft v.0.1b



####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

#################
### Integrals ###
#################

### y = x * I(1/log(x + k)) dx + F0(x)

# Check:
x = sqrt(3); k = 4/3;
# x = 1/sqrt(3); k = 0;
params = list(x=x, k=k);
I = integrate(\(x) 1/log(x+k), 0, x)$value;
e = expression(x * I + x^2/2 - 3*x)[[1]];
f0 = x^2/2 - 3*x; df0 = x - 3; d2f0 = 1;
#
y   = eval(e, params);
dI  = 1/log(x+k);
dy  = eval(D(e, "x"), params) + x*dI;
d2y = eval(D(D(e, "x"), "x"), params) + 2*dI - x/(x+k)/log(x+k)^2; # 2*dI!

### D(y)
dy   # =
I + x/log(x + k) + df0
x*dy # =
y - f0 + x^2/log(x + k) + x*df0

### D2(y)
x*d2y + dy # =
dy + x*d2f0 + 2*x/log(x + k) - x^2/((x+k)*log(x+k)^2);
dy + x*d2f0 + 2*(x*dy - y + f0 - x*df0)/x - (x*dy - y + f0 - x*df0)^2/(x^2*(x+k));
#
x^3*(x+k)*d2y # =
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

