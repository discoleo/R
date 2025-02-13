########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs - Trigonometric
##
## draft v.0.1a


### Non-Linear: Trigonometric Integrals
# - derived from I(Trig);


###############
### History ###
###############

### draft v.0.1a:
# - moved I(Trig) variants to this file
#   from file: DE.ODE.Trigonometric.R;

### [previous file]

### draft v.0.4a:
# - derived from Inverse Integrals:
#   y * I(sin(x^n)) = F0(x);
### draft v.0.3h:
# - minor extension to: [v.0.3f];
### draft v.0.3g:
# - derived from:
#   y = x * I(tan(k*x^2)) + F0(x);
### draft v.0.3f - v.0.3f-test:
# - derived from:
#   y = x * I(atan(k*x)^2) + F0(x); [extended in v.0.3h]
### draft v.0.3c-pre:
# - derived from:
#   y = I(sin(x^n)) * I(cos(x^n));


#################

#################
### Integrals ###
#################

### y = I(sin(x^n)) * I(cos(x^n))

### D(y)
sin(x^n) * Ic + cos(x^n) * Is

### D2(y)
n*x^(n-1)*cos(x^n) * Ic - n*x^(n-1)*sin(x^n) * Is +
	+ sin(2*x^n)

### Solve Linear system:
# TODO:
Is = (n*x^(n-1)*cos(x^n)*dy - sin(x^n)*d2y + sin(x^n)*sin(2*x^n));
Ic = (n*x^(n-1)*sin(x^n)*dy + cos(x^n)*d2y - cos(x^n)*sin(2*x^n));
# TODO: check;


#########################
#########################

### y = sin(x^n) * I(cos(x^n)) - cos(x^n) * I(sin(x^n))

### D(y)
n*x^(n-1)*cos(x^n) * Ic + n*x^(n-1)*sin(x^n) * Is

### D2(y)
n*x^(n-1) - n^2*x^(2*n-2)*(sin(x^n) * Ic - cos(x^n) * Is) + (n-1) * dy / x
n*x^(n-1) - n^2*x^(2*n-2)*y + (n-1) * dy / x

### ODE:
x*d2y - (n-1)*dy + n^2*x^(2*n-1)*y - n*x^n # = 0


### Solution & Plot:
y.sin = function(x, n=2, lower=0) {
	sapply(x, function(upper)
		integrate(function(x) sin(x^n), lower=lower, upper=upper)$value )
}
y.cos = function(x, n=2, lower=0) {
	sapply(x, function(upper)
		integrate(function(x) cos(x^n), lower=lower, upper=upper)$value )
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
curve(dy(x, n=n), add=T, col="green")
sapply(x.px, line.tan, dx=1.5, p=dy, dp=d2y, n=n, col="orange")


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

### Extension
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


############################
############################

### y = x * I(atan(k*x)^2) + F0(x)

### D(y)
(y - f0)/x + x*atan(k*x)^2 + df0
### x*dy =
(y - f0) + x^2*atan(k*x)^2 + x*df0

### D2(y)
# x*d2y + dy =
dy + 2*x*atan(k*x)^2 + 2*k*x^2 / (k^2*x^2+1) * atan(x) + x*d2f0
dy + 2*(x*dy - y + f0 - x*df0)/x +
	+ 2*k*x / (k^2*x^2+1) * sqrt(x*dy - y + f0 - x*df0) + x*d2f0
# x^2*(k^2*x^2+1)*d2y =
2*(k^2*x^2+1)*(x*dy - y + f0 - x*df0) +
	+ 2*k*x^2*sqrt(x*dy - y + f0 - x*df0) + x^2*(k^2*x^2+1)*d2f0

### ODE:
(x^2*(k^2*x^2+1)*d2y - 2*(k^2*x^2+1)*(x*dy - y + f0 - x*df0) - x^2*(k^2*x^2+1)*d2f0)^2 +
	- 4*k^2*x^4*(x*dy - y + f0 - x*df0) # = 0

### Examples:
### f0 = x; df0 = 1;
(x^2*(k^2*x^2+1)*d2y - 2*(k^2*x^2+1)*(x*dy - y))^2 - 4*k^2*x^4*(x*dy - y) # = 0


### Solution & Plot:
y.I = function(x, k=1, n=2, lower=0) {
	sapply(x, function(upper)
		integrate(function(x) atan(k*x)^n, lower=lower, upper=upper)$value)
}
deriv.pol = function(x, b, dn=1, x.mult=0) {
	coeff = tail(b, -1);
	pow = seq(length(coeff));
	coeff = coeff * pow;
	if(dn > 1) return(deriv.pol(x, coeff, dn=dn-1, x.mult=x.mult))
	if(x.mult != 1) pow = pow + x.mult - 1; # x * df0;
	sapply(x, function(x) sum(coeff * x^pow));
}
eval.pol = function(x, b) {
	pow = seq(0, length(b) - 1)
	sapply(x, function(x) sum(b * x^pow))
}
y = function(x, b=0, k=1, n=2, lower=0) {
	# x * I(atan(x)^2) + F0(x)
	I.v = y.I(x, k=k, n=n, lower=lower)
	r = x*I.v + eval.pol(x, b);
	return(r)
}
dy = function(x, b=0, k=1, n=2, lower=0) {
	I.v = y.I(x, k=k, n=n, lower=lower)
	# (y - f0) + x^2*atan(k*x)^2 + x*df0
	# y - f0 = x*I.v;
	r = x*I.v + x^2 * atan(k*x)^2 + deriv.pol(x, b, x.mult=1); # x*df0
	div = x;
	r = ifelse(div != 0, r / div,
		dy(x + 1E-3, b=b, k=k, n=n, lower=lower)) # TODO
	return(r)
}
d2y = function(x, b=0, k=1, n=2, lower=0) {
	y.x = y(x, b=b, k=k, n=n, lower=lower)
	dy.x = dy(x, b=b, k=k, n=n, lower=lower)
	# (x^2*(k^2*x^2+1)*d2y - 2*(k^2*x^2+1)*(x*dy - y + f0 - x*df0) - x^2*(k^2*x^2+1)*d2f0)^2 +
	#  - 4*k^2*x^4*(x*dy - y + f0 - x*df0) # = 0
	x2 = x^2; x21 = k^2 * x2 + 1;
	f0 = eval.pol(x, b); xdf0 = deriv.pol(x, b, dn=1, x.mult=1);
	d2f0 = deriv.pol(x, b, dn=2, x.mult=0);
	Tx = x*dy.x - y.x + f0 - xdf0;
	dp  = 4*k^2*x2*x2*Tx;
	dp = sqrt(dp);
	dp = dp + 2*x21*Tx + x2*x21*d2f0;
	div = x2 * x21;
	dp = ifelse(div != 0, dp / div, d2y(x + 1E-3, b=b, n=n, lower=lower)); # TODO
	return(dp)
}
### Plot:
b = c(0, 1);
x.px = c(-5:7 * 3/7)
curve(y(x, b=b), from= -3, to = 2, ylim=c(-2, 3.5))
line.tan(x.px, dx=1.6, p=y, dp=dy, b=b)
#
curve(dy(x, b=b), add=T, col="green")
line.tan(x.px, dx=1.5, p=dy, dp=d2y, b=b, col="orange")


### Ex 2:
b = c(0, -2, 1);
x.px = c(-5:7 * 3/7)
curve(y(x, b=b), from= -1.5, to = 2.5, ylim=c(-3, 4))
line.tan(x.px, dx=1.6, p=y, dp=dy, b=b)
#
curve(dy(x, b=b), add=T, col="green")
line.tan(x.px, dx=1.5, p=dy, dp=d2y, b=b, col="orange")


### Ex 3:
b = c(0, -2, 1); k =1/3;
x.scale = 1; xlim = c(-1.5, 2.5) * x.scale;
x.px = c(-5:7 * 3/7) * x.scale;
curve(y(x, b=b, k=k), from= xlim[1], to = xlim[2], ylim=c(-3, 4))
line.tan(x.px, dx=1.6, p=y, dp=dy, b=b, k=k)
#
curve(dy(x, b=b, k=k), add=T, col="green")
line.tan(x.px, dx=1.5, p=dy, dp=d2y, b=b, k=k, col="orange")


#########################
#########################

#########################
### Inverse Integrals ###
#########################

### y * I(sin(x^n)) = F0(x)

### D(y)
I*dy + y*sin(x^n) - df0 # = 0 # * y =>
f0*dy + sin(x^n)*y^2 - y*df0 # = 0

### D2(y)
# Note: basic equation contains cos(x^n);
f0*d2y + 2*sin(x^n)*y*dy + n*x^(n-1)*cos(x^n)*y^2 - y*d2f0 # = 0
f0*y*d2y - 2*(f0*dy - y*df0)*dy + n*x^(n-1)*cos(x^n)*y^3 - y^2*d2f0 # = 0
f0*y*d2y - 2*(f0*dy - y*df0)*dy - y^2*d2f0 + n*x^(n-1)*sqrt(y^4 - (f0*dy - y*df0)^2)*y # = 0
(f0*y*d2y - 2*(f0*dy - y*df0)*dy - y^2*d2f0)^2 +
	- n^2*x^(2*n-2)*(y^4 - (f0*dy - y*df0)^2)*y^2 # = 0

### TODO: check;

