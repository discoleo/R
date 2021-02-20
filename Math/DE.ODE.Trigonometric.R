########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Trigonometric
###
### draft v.0.3d-pre


### History

### Order 1 Non-Liniar & Linear:
### Trigonometric Variants

### Order 1 & 2 Linear:
### Trigonometric Variants
###
### draft v.0.2e - v.0.2f:
# - minor fixes in older formulas;
# - better comments & formatting;
### draft v.0.2c: [06-12-2020]
# - re-organizing sections;
### draft v.0.2b - v.0.2b-ex:
# - based on solving for sin/cos(log(P(x))):
#   (x+b)^2 * d2y + (x+b)*dy + y = 0;
# - more examples (based on generalization); [v.0.2b-ex]
#   [are a special case of v.0.2a]
### draft v.0.2a - v.0.2a-form: [2020-12-04]
# - based on solving for sin/cos:
#   (2*x+b) * d2y - 2*dy + (2*x+b)^3 * y = 0;
#   [including partial generalization]
# - added a more formal approach & generalization; [v.0.2a-form]
#   dP*d2y - d2P*dy + dP^3 * y = 0; [full eq. in v.0.2b-gen]

### Order 1 & 2 Non-Liniar:
### Trigonometric Variants
###
### draft v.0.3d-pre:
# - ODE: d2y + n^2*x^(2*n-2)*y = n*x^(n-1);
### draft v.0.3c-pre:
# - derived from:
#   y = I(sin(x^n)) * I(cos(x^n));
### draft v.0.3b:
# - linearly combined Trig functions of y:
#   f1*sin(y) + f2*cos(y) = P(x);
### draft v.0.3a: [27-01-2021]
# - derived from: y^2 = x^p * sin(x^m);
#   2*x^4*y*d2y + 2*x^4*dy^2 + y^2 = 0;
### draft v.0.2d: [06-12-2020]
# - generalization:
#   P(y)*sin(P(y)) + cos(P(y)) = f(x);
### ... [see Section on Liniar Fx]
### draft v.0.1e - v.0.1g:
# - solved:
#   x^2*dy*d2y + x/(x+1) * dy^2 + (x+1)^2 * y*dy = 0;
#   [includes generalization]
#   x*y^2*d2y + x*y*dy^2 - (k^2+1)/2 * x^2*dy^3 - y^2*dy = 0; [v.0.1f]
# - another simple trigonometric example: (type: sin(y^2) = sqrt(x))
#   2*x*(1 - x)*y^3*d2y + (1 - 2*x)*y^3*dy + 1/8 = 0; [v.0.1g]
# - fixed comments [minor fix];
### draft v.0.1d: [14-11-2020]
# - integration by parts:
#   x*d4z - 3/2 * d3z + 9/4*a^2*x^2*d2z - 9/2*a^2*x*dz + 9/2*a^2*z = 0;
### draft v.0.1c:
# - moved Section on Trigonometric Functions
#   from DE.ODE.Polynomial.R;
# - type: sin(x*y + a) = x^p;
### draft v.0.1b:
# - added types of form: y*sin(y) + cos(y) = f(x);
#   y*d2y + 2*(dy)^2 - (x + b)*y*(dy)^3 = 0;
### draft v.0.1a: [08-11-2020]
# - examples for type: sin(y^2) = f(x)
#   x*(1-x^4)*y^3*d2y - (x^4+1)*y^3*dy + x^3 = 0;
#   (x^2 - 2)*y^3*d2y + x*y^3*dy = 1;
# - moved section Trigonometric Variants
#   to this new file;
# - renamed to Trigonometric (variants);

### [old file] DE.ODE.Fractions.Lambert.R
### draft v.0.1h:
# - solved a trigonometric type:
#   x*(1-x^4)*y^3*d2y - (x^4+1)*y^3*dy + x^3 = 0;
# - [TODO] move to separate file; [DONE]


#########################

### Helper functions

library(pracma)
# needed for Lambert W;


# include: DE.ODE.Helper.R;
source("DE.ODE.Helper.R")

#########################
#########################

### Trigonometric Functions:
### Simple


#########################
### Section A: Non-Linear
### Higher Power

### y^2 = x^p * sin(x^m)

### D =>
# 2*x*y*dy = p*y^2 + m*x^(p+m)*cos(x^m)

### D2 =>
2*x^2*y*d2y + 2*x^2*dy^2 - 2*(2*p + m - 1)*x*y*dy +
	+ m^2*x^(2*m)*y^2 + p*(p+m)*y^2 # = 0

### Special cases:

### m = 1
2*x^2*y*d2y + 2*x^2*dy^2 - 4*p*x*y*dy + x^2*y^2 + p*(p+1)*y^2 # = 0

### p = - m
2*x^2*y*d2y + 2*x^2*dy^2 + 2*(m + 1)*x*y*dy + m^2*x^(2*m)*y^2 # = 0
### p = - 1; m = 1;
2*x*y*d2y + 2*x*dy^2 + 4*y*dy + x*y^2 # = 0
### p = 1; m = -1;
2*x^4*y*d2y + 2*x^4*dy^2 + y^2 # = 0

### Plot:
y = function(x, pp=1, m=1, posRoot=TRUE) {
	# root
	y = x^pp * sin(x^m)
	y = sqrt(y)
	if( ! posRoot) y = -y;
	return(y)
}
dy = function(x, pp=1, m=1, posRoot=TRUE, y.x) {
	if(missing(y.x)) y.x = y(x, pp=pp, m=m, posRoot=posRoot);
	dp = pp*y.x^2 + m*x^(pp+m)*cos(x^m)
	div = 2*x*y.x;
	dp = ifelse(div != 0, dp / div, 0); # TODO: may need correction
	return(dp)
}
d2y = function(x, pp=1, m=1, posRoot=TRUE) {
	y.x  = y(x, pp=pp, m=m, posRoot=posRoot);
	dy.x = dy(x, pp=pp, m=m, posRoot=posRoot, y.x=y.x);
	dp = 2*x^2*dy.x^2 - 2*(2*pp + m - 1)*x*y.x*dy.x + m^2*x^(2*m)*y.x^2 + pp*(pp+m)*y.x^2;
	div = - 2*x^2*y.x;
	dp = ifelse(div != 0, dp / div, (2*pp + m - 1) / m); # TODO: needs correction!
	return(dp)
}
### Test

pp = 2; m = 1;
curve(y(x, pp=pp, m=m), from= 0, to = 2)
# global minimum;
sapply(c((0:4)/2.2), line.tan, dx=1.5, p=y, dp=dy, pp=pp, m=m)
# pseudo-sigmoidal
curve(dy(x, pp=pp, m=m), add=T, col="green")
sapply(c((0:4)/2.2), line.tan, dx=1/5, p=dy, dp=d2y, pp=pp, m=m, col="orange")


### m = 2
pp = 2; m = 2;
curve(y(x, pp=pp, m=m), from= -5/3, to = 5/3, ylim=c(-2, 2))
# global minimum;
sapply(c((-3:3)/2.2), line.tan, dx=1.5, p=y, dp=dy, pp=pp, m=m)
# pseudo-sigmoidal
curve(dy(x, pp=pp, m=m), add=T, col="green")
sapply(c((-3:3)/2.2), line.tan, dx=1/5, p=dy, dp=d2y, pp=pp, m=m, col="orange")


#########################
### Section B: Non-Linear
### Inverse Trigonometric

### sin(y^2) = f(x);
(1 - f^2)*df*y^3*d2y + (d2f*f^2 - df^2*f - d2f)*y^3*dy + 1/4 * df^3 # = 0

### Generalization:
### sin(y^n) = f(x)
### D =>
# n*y^(n-1)*cos(y^n)*dy = df;
### D2 =>
n*y^(n-1)*cos(y^n)*d2y - n^2*y^(2*n-2)*sin(y^n)*dy^2 + n*(n-1)*y^(n-2)*cos(y^n)*dy^2 - d2f # = 0
df/dy * d2y - n^2*y^(2*n-2)*f*dy^2 + (n-1)*df/y * dy - d2f # = 0
### Eq:
df*y*d2y - n^2*f*y^(2*n-1)*dy^3 + (n-1)*df*dy^2 - d2f*y*dy # = 0
### n^2*y^(2*n-2)*cos(y^n)^2*dy^2 = df^2
# dy^2 = 1/n^2 * y^(2-2*n) * df^2 / (1 - f^2)
### Alternative Eq:
(1-f^2)*df*y*d2y + (f^2*d2f - f*df^2 - d2f)*y*dy + (n-1)/n^2*df^3*y^(2-2*n) # = 0 # * y^(2*n-2)
(1-f^2)*df*y^(2*n-1)*d2y + (f^2*d2f - f*df^2 - d2f)*y^(2*n-1)*dy + (n-1)/n^2*df^3 # = 0

### Examples:
# [not run]
# f = x^2;
x*(1-x^4)*y^3*d2y - (x^4+1)*y^3*dy + x^3 # = 0
# f = x^2 + b;
x*(1-(x^2+b)^2)*y^3*d2y - (x^4 - b^2 + 1)*y^3*dy + x^3 # = 0
# f = sqrt(x); df = 1/2 / sqrt(x); d2f = -1/4 * x^(-3/2);
# df*(1 - x)*y^3*d2y - d2f*(1 - 2*x)*y^3*dy + 1/4 * 1/4 * 1/x * df # = 0 # * 4 * x^(3/2)
# 2*x*(1 - x)*y^3*d2y + (1 - 2*x)*y^3*dy + 1/8 # = 0

### Plot:
y = function(x, b=0) {
	# root
	y = asin(x^2 + b)
	y[y < 0] = y[y < 0] + 2*pi;
	y = sqrt(y)
	y = sapply(y, round0)
	return(y)
}
dy = function(x, b=0, y.x) {
	if(missing(y.x)) y.x = y(x, b=b);
	div = y.x * sqrt(1 - (x^2+b)^2)
	dp = x;
	dp = ifelse(div != 0, dp / div, -1 - b); # may need correction
	return(dp)
}
d2y = function(x, b=0) {
	y.x = y(x, b=b)
	z = dy(x, b=b, y.x=y.x)
	x2 = x*x; x4 = x2*x2;
	div = x*(1 - (x2 + b)^2)*y.x^3;
	dp = (x4 - b^2 + 1)*y.x^3*z - x^3;
	dp = ifelse(div != 0, dp / div, -1); # TODO: needs correction!
	return(dp)
}

### Test

### b = 0;
curve(y(x), from= -1, to = 1)
# global minimum;
sapply(c((-2:2)/2.2), line.tan, dx=3, p=y, dp=dy)
# pseudo-sigmoidal
curve(dy(x), add=T, col="green")
sapply(c(-2, -1.5, 1.5, 2)/2.1, line.tan, dx=1/5, p=dy, dp=d2y, col="orange")

# check full pseudo-sigmoidal
curve(dy(x), from= -1, to = 1, col="green")
sapply(c(-2, -1.5, 1.5, 2)/2.1, line.tan, dx=1/5, p=dy, dp=d2y, col="orange")


### b = -1;
# x^3*(x^2 - 2)*y^3*d2y + x^4*y^3*dy - x^3 = 0
# (x^2 - 2)*y^3*d2y + x*y^3*dy = 1
b = -1;
#
curve(y(x, b=b), from= -sqrt(2) + 1E-10, to = sqrt(2) - 1E-10)
# global minimum;
sapply(c((-2:2)/2.2), line.tan, dx=3, p=y, dp=dy, b=b)
# pseudo-sigmoidal
curve(dy(x, b=b), add=T, col="green")
sapply(c(-2, -1.5, 1.5, 2)/2.1, line.tan, dx=1/5, p=dy, dp=d2y, b=b, col="orange")

# check full pseudo-sigmoidal
curve(dy(x, b=b), from= -sqrt(2) + 1E-10, to = sqrt(2) - 1E-10, col="green", ylim=c(-4,4))
sapply(c(-2.8, -2.4, -2, -1.5, 1.5, 2, 2.4, 2.8)/2.1, line.tan, dx=1/5, p=dy, dp=d2y, b=b, col="orange")


### b = -1/2;
b = -1/2
curve(dy(x, b=b), from= -sqrt(5/4), to = sqrt(5/4), col="green", ylim=c(-4,4))
sapply(c(-2.4, -2, -1.5, 1.5, 2, 2.4)/2.1, line.tan, dx=1/5, p=dy, dp=d2y, b=b, col="orange")


### f = sqrt(x);
# df = 1/2 / sqrt(x); d2f = -1/4 * x^(-3/2);
# df*(1 - x)*y^3*d2y - d2f*(1 - 2*x)*y^3*dy + 1/4 * 1/4 * 1/x * df # = 0 # * 4 * x^(3/2)
2*x*(1 - x)*y^3*d2y + (1 - 2*x)*y^3*dy + 1/8 # = 0

### Plot:
y = function(x, b=0) {
	# root: TODO: compute also with b;
	y = asin(sqrt(x) + b)
	y[y < 0] = y[y < 0] + 2*pi;
	y = sqrt(y)
	y = sapply(y, round0)
	return(y)
}
dy = function(x, b=0, y.x) {
	if(missing(y.x)) y.x = y(x, b=b);
	div = y.x * sqrt(1 - (x+b)^2) * sqrt(x) # cos(y^2) = sqrt(...);
	dp = 1/4;
	dp = ifelse(div != 0, dp / div, 1E+3); # may need correction
	return(dp)
}
d2y = function(x, b=0) {
	y.x  = y(x, b=b)
	dy.x = dy(x, b=b, y.x=y.x)
	div = - 2*x*(1 - x)*y.x^3;
	dp  = (1 - 2*x)*y.x^3*dy.x + 1/8
	dp = ifelse(div != 0, dp / div, 1E+3); # may need correction
	return(dp)
}
### b = 0;
curve(y(x), from= 0, to= 1, ylim=c(0, 1.5))
# quasi/inverted sigmoidal;
sapply(c(0, 1/5, 2/3, 0.94), line.tan, dx=3, p=y, dp=dy)
# global minimum
curve(dy(x), add=T, col="green")
sapply(c(0:4/5, 0.9), line.tan, dx=1/5, p=dy, dp=d2y, col="orange")


# check full D(y) curve:
curve(dy(x), from= 0, to = 1, col="green", ylim=c(0, 2))
sapply(c(0:4/5, 0.9), line.tan, dx=1/5, p=dy, dp=d2y, col="orange")


############################
############################

### y*sin(y) + cos(y) = f(x)

### D =>
y*dy*cos(y) - df # = 0
# D2 =>
df*y*d2y - f*y*(dy)^3 + 2*df*(dy)^2 - d2f*y*dy # = 0
### Alternative Eq:
# y*dy = df / cos(y)
# needs also D3;
df*y*d2y - f*y*(dy)^3 + 2*df^3/(y^2 - (f-cos(y))^2) - d2f*y*dy # = 0

### Generalization:
### P(y)*sin(P(y)) + cos(P(y)) = f(x)
### D =>
# P*cos(P)*dP = df
# P*sin(P) = f - df / (P*dP)
### D2 =>
P*cos(P)*d2P + cos(P)*dP^2 - P*sin(P)*dP^2 - d2f # = 0
df/dP*d2P + df/P*dP - (f*dP - df/P)*dP - d2f # = 0
df*P*d2P - f*P*dP^3 + 2*df*dP^2 - d2f*P*dP # = 0

### Examples:

### Generalized Case:
### P(y) = y^2
2*df*y^2*(y*d2y + dy^2) - 8*f*y^5*dy^3 + 8*df*y^2*dy^2 - 2*d2f*y^3*dy # = 0
df*(y*d2y + dy^2) - 4*f*y^3*dy^3 + 4*df*dy^2 - d2f*y*dy # = 0
df*y*d2y - 4*f*y^3*dy^3 + 5*df*dy^2 - d2f*y*dy # = 0
### P(y) = y^2; f(x) = x + b;
y*d2y - 4*(x+b)*y^3*dy^3 + 5*dy^2 # = 0


### Simple Case:
# f = x + b, where b = constant;
y*d2y - (x+b)*y*(dy)^3 + 2*(dy)^2 # = 0
### Solution:
y = function(x, b) {
	# root
	y.f = function(x, v) x*sin(x) + cos(x) - v - b;
	dy.f = function(x, v) x*cos(x);
	x0.f = function(x) {
		xb = x + b
		x0 = if(xb >= 1 & xb <= pi/2) 1 else if(xb < 1) 3 else 7;
		return(x0);
	}
	y = sapply(x, function(x) newtonRaphson(y.f, x0.f(x), dfun=dy.f, v=x)[[1]])
	y = sapply(y, round0)
	return(y)
}
dy = function(x, b, y.x) {
	if(missing(y.x)) y.x = y(x, b=b);
	div = y.x * cos(y.x)
	dp = 1;
	dp = ifelse(div != 0, dp / div, -1); # may need correction
	return(dp)
}
d2y = function(x, b) {
	y.x = y(x, b=b)
	z = dy(x, b=b, y.x=y.x)
	div = y.x
	dp = - 2*z^2 + (x+b)*y.x*z^3;
	dp = ifelse(div != 0, dp / div, -1); # TODO: needs correction!
	return(dp)
}
###
b = -1/2
curve(y(x, b=b), from= -3, to = 3, ylim=c(-2,7))
# oscillating function with local minimum;
sapply(c(-1, 1.3, 1.6, 1.8, 1.95), line.tan, dx=3, p=y, dp=dy, b=b)
# spikes
curve(dy(x, b=b), from= -3, to = 3, add=T, col="green")
sapply(c(-1, 1.3, 1.6, 1.8, 1.95), line.tan, dx=3, p=dy, dp=d2y, b=b, col="orange")


#####################
### Generalized Case:
# P(y) = y^2
# f = x + b, where b = constant;
y*d2y - 4*(x+b)*y^3*dy^3 + 5*dy^2 # = 0
### Solution:
y = function(x, b, n=2) {
	# root
	y.f = function(x, v) {
		x = if(n == 1) x else x^n;
		x*sin(x) + cos(x) - v - b;
	}
	dy.f = function(x, v) {
		if(n == 1) return(x*cos(x));
		xn = if(n == 1) x else x^n;
		n*xn*xn/x*cos(xn);
	}
	x0.f = function(x) {
		xb = x + b
		x0 = if(xb >= 1 & xb <= pi/2) 3/4 else if(xb < 1) 1.5 else 2.5;
		return(x0);
	}
	y = sapply(x, function(x) newtonRaphson(y.f, x0.f(x), dfun=dy.f, v=x)[[1]])
	# if(n != 1) y = rootn(y, n)
	y = sapply(y, round0)
	return(y)
}
dy = function(x, b, y.x, n=2, dF=1) {
	if(missing(y.x)) y.x = y(x, b=b, n=n);
	dp = if(n == 1) dF else dF * y.x;
	y.x = if(n == 1) y.x else y.x^n;
	div = n * y.x * y.x * cos(y.x)
	dp = ifelse(div != 0, dp / div, 1E+3); # TODO: needs correction
	return(dp)
}
d2y = function(x, b, n=2) {
	y.x = y(x, b=b, n=n)
	dy.x = dy(x, b=b, y.x=y.x, n=n)
	div = y.x
	dp = 4*(x+b)*y.x^3*dy.x^3 - 5*dy.x^2;
	dp = ifelse(div != 0, dp / div, -1); # TODO: needs correction!
	return(dp)
}
###
b = -1/2
curve(y(x, b=b), from= -3, to = 3, ylim=c(-2,7))
# oscillating/spikes (ECG-like) function with local minima;
sapply(c(seq(1.4, 2.2, by=0.2) - 0.05), line.tan, dx=1.5, p=y, dp=dy, b=b)
# spikes
curve(dy(x, b=b), from= -3, to = 3, add=T, col="green")
sapply(c(seq(1.4, 2.2, by=0.2) - 0.05), line.tan, dx=1.4, p=dy, dp=d2y, b=b, col="orange")


#####################

#####################
### Trigonometric ###
###   Functions   ###

### Type: sin(f(x) * g(y))

### sin(x*y) = x^2
x*(1 - x^4)*d2y - (4*x^4 - 2)*dy - 2*x^3*y = 2*sqrt(1 - x^4)
### Solution
# x*y = t => dy = dt/x - t/x^2; d2y = d2t/x - 2*dt/x^2 + 2*t/x^3;
y = function(x, n=2) {
	r = asin(x^n)/x
	r[x == 0] = 0
	return(r)
}
dy = function(x, n=2) {
	y.x = y(x, n=n)
	dp = (n*x^(n-1) / sqrt(1- x^(2*n))) - y.x
	dp = dp/x
	zero = if(n == 2) 1 else 0; # TODO: dependent on n;
	dp[x == 0] = zero
	return(dp)
}
d2y = function(x, n=2) {
	y.x = y(x, n=n)
	dy.x = dy(x, n=n)
	x4 = x^(2*n)
	dp = if(n == 2) 1 else (n-1)*x^(n-2);
	dp = (dp * sqrt(1 - x4) + x^(2*n-1)*y.x)*n + ((n+2)*x4 - 2)*dy.x
	dp = dp / x / (1-x4)
	return(dp)
}
### Plot
curve(y, from=-1, to=1)
div = 23
sapply(c(-1 + (1:4)/div, 1 - (1:4)/div), line.tan, dx=0.5, p=y, dp=dy)
### D2(y):
curve(dy, from=-1, to=1, col="green", ylim=c(-1, 6))
curve(y, from=-1, to=1, col="grey", add=T)
div = 23
sapply(c(-1 + (1:4)/div, 1 - (1:4)/div), line.tan, dx=0.5, p=dy, dp=d2y)


### sin(x*y + a) = x^p
# Note: parameter [a] does NOT seem to have any impact;
(x*dy + y)*cos(x*y + a) = p*x^(p-1)
(x*dy + y) * sqrt(1 - x^(2*p)) = p*x^(p-1)
# D2 =>
(x*d2y + 2*dy) * sqrt(1 - x^(2*p)) - p*x^(2*p-1)*(x*dy + y)/sqrt(1 - x^(2*p))  = p*(p-1)*x^(p-2)
(x*d2y + 2*dy) * (1 - x^(2*p)) - p*x^(2*p-1)*(x*dy + y)  = p*(p-1)*x^(p-2)*sqrt(1 - x^(2*p))
x*(1 - x^(2*p))*d2y - ((p+2)*x^(2*p) - 2)*dy - p*x^(2*p-1)*y  = p*(p-1)*x^(p-2)*sqrt(1 - x^(2*p))
### p =3
x*(1 - x^6)*d2y - (5*x^6 - 2)*dy - 3*x^5*y  = 6*x*sqrt(1 - x^6)
### Plot:
### D2(y):
curve(dy(x, n=3), from=-1, to=1, col="green", ylim=c(-4, 6))
curve(y(x, n=3), from=-1, to=1, col="grey", add=T)
div = 23
sapply(c(-1 + (1:4)/div, 1 - (1:4)/div), line.tan, dx=0.5, p=dy, dp=d2y, n=3)

### Combinations:
# [not run]
x*(1 - x^6)*d2y - (5*x^6 - 2)*dy - 3*x^2*y*sin(x*y + a) - 6*x*sqrt(1 - x^6) # = 0
### Solution & Plot:
y = function(x, n=2, a=1/3) {
	r = (asin(x^n) - a)/x
	r[x == 0] = if(a == 0) 0 else Inf;
	return(r)
}
dy = function(x, n=2, a=1/3) {
	y.x = y(x, n=n, a=a)
	dp = (n*x^(n-1) / sqrt(1- x^(2*n))) - y.x
	dp = dp/x
	zero = if(n == 2) 1 else 0; # TODO: dependent on n;
	zero = if(a == 0) zero else Inf;
	dp[x == 0] = zero
	return(dp)
}
d2y = function(x, n=2, a=1/3) {
	y.x = y(x, n=n, a=a)
	dy.x = dy(x, n=n, a=a)
	x4 = x^(2*n)
	dp = if(n == 2) 1 else (n-1)*x^(n-2);
	dp = (dp * sqrt(1 - x4) + x^(n-1)*y.x*sin(x*y.x + a))*n + ((n+2)*x4 - 2)*dy.x
	dp = dp / x / (1-x4)
	return(dp)
}
### Plot
curve(y, from=-1, to=1)
div = 4.5
sapply(c(-1 + (1:4)/div, 1 - (1:4)/div), line.tan, dx=0.5, p=y, dp=dy)
### D2(y):
curve(dy, from=-1, to=1, col="green", ylim=c(-1, 8))
curve(y, from=-1, to=1, col="grey", add=T)
div = 4.5
sapply(c(-1 + (1:4)/div, 1 - (1:4)/div), line.tan, dx=0.5, p=dy, dp=d2y)


### TODO: tan, ln;

########################

########################
### Integration by parts

# based on Linear simple (Basics);
# y = sin(a*x^n);

### n = 1
x*d2y + a^2*x*y # = 0
### I() =>
x*dy - y + a^2*x*I(y) - a^2*I(I(y)) # = 0
# z = I(I(y))
x*d3z - d2z + a^2*x*dz - a^2*z # = 0
# D(d2z/x) + a^2*(z/x) = 0;

### n = 3/2
# x*d2y - (n-1)*dy + n^2*a^2*x^(2*n-1)*y # = 0
x*d2y - 1/2 * dy + 9/4*a^2*x^2*y # = 0
### I() =>
x*dy - y - 1/2 * y + 9/4*a^2*x^2*I(y) - 9/2*a^2*x*I(I(y)) + 9/2*a^2*I(I(I(y))) # = 0
# z = I(I(I(y)))
x*d4z - 3/2 * d3z + 9/4*a^2*x^2*d2z - 9/2*a^2*x*dz + 9/2*a^2*z # = 0
### TODO: test;


### Test
# only n == 1
y.base = function(x, a=1, n=1) {
	r = sin(a*x^n)
	return(r)
}
dy.base = function(x, a=1, n=1) {
	xn1 = if(n == 2) x else x^(n-1);
	xn  = xn1 * x;
	dp = n*a*xn1*cos(a*xn)
	return(dp)
}
y = function(x, a=1, n=1) {
	d3z = dy.base(x, a=a, n=n)
	d2z = y.base(x, a=a, n=n)
	dz  = dy(x, a=a, n=n)
	z = (x*d3z - d2z + a^2*x*dz) / a^2
	return(z)
}
dy = function(x, a=1, n=1, lower=0) {
	r = sapply(x, function(x) integrate(y.base, lower=lower, upper=x, a=a, n=n)$value)
	return(r)
}
d2y = function(x, a=1, n=1) {
	return(y.base(x, a=a, n=n))
}
### Plot:
a = 1; n = 1;
curve(y(x, a=a, n=n), from= -3, to= 3, ylim=c(-2, 1.5))
# sinus wave;
sapply(c((-5:5)/2.2), line.tan, dx=3, p=y, dp=dy, a=a, n=n)
# wave
curve(dy(x, a=a, n=n), add=T, col="green")
sapply(c((-5:5)/2.2), line.tan, dx=1/5, p=dy, dp=d2y, a=a, n=n, col="orange")

### Test separately:
curve(dy(x, a=a, n=n), from= -3, to= 3, col="green")
sapply(c((-5:5)/2.2), line.tan, dx=1/5, p=dy, dp=d2y, a=a, n=n, col="orange")


### Test
# only n == 3/2
# x*d4z - 3/2 * d3z + 9/4*a^2*x^2*d2z - 9/2*a^2*x*dz + 9/2*a^2*z # = 0
y.base = function(x, a=1, n=3/2) {
	r = sin(a*x^n)
	return(r)
}
dy.base = function(x, a=1, n=3/2) {
	xn1 = if(n == 2) x else x^(n-1);
	xn  = xn1 * x;
	dp = n*a*xn1*cos(a*xn)
	return(dp)
}
y = function(x, a=1, n=3/2) {
	d4z = dy.base(x, a=a, n=n)
	d3z = y.base(x, a=a, n=n)
	d2z = d2y(x, a=a, n=n)
	dz  = dy(x, a=a, n=n)
	z = -(x*d4z - 3/2 * d3z + 9/4*a^2*x^2*d2z - 9/2*a^2*x*dz) * 2/9 / a^2
	return(z)
}
d2y = function(x, a=1, n=3/2, lower=0) {
	r = sapply(x, function(x) integrate(y.base, lower=lower, upper=x, a=a, n=n)$value)
	return(r)
}
dy = function(x, a=1, n=3/2, lower=0) {
	r = sapply(x, function(x) integrate(d2y, lower=lower, upper=x, a=a, n=n)$value)
	return(r)
}
### Plot
a = 1; n = 3/2;
# quasi-exponential;
curve(y(x, a=a, n=n), from= 0, to= 3)
sapply(c((0:5)/2.2), line.tan, dx=3, p=y, dp=dy, a=a, n=n)
# wave
curve(dy(x, a=a, n=n), add=T, col="green")
sapply(c((0:5)/2.2), line.tan, dx=1/5, p=dy, dp=d2y, a=a, n=n, col="orange")


### Test separately:
curve(dy(x, a=a, n=n), from= 0, to= 3, col="green")
sapply(c((0:5)*2/3.2), line.tan, dx=1/5, p=dy, dp=d2y, a=a, n=n, col="orange")


############################
############################

### Special / Combined Functions

### y*sin(k * log(y)) = f(x)

### y*sin(k * log(y)) = x^2
# k*dy * cos(k * log(y)) = (2*x*y - x^2*dy) / y;
x*y^2*d2y + x*y*dy^2 - (k^2+1)/2 * x^2*dy^3 - y^2*dy # = 0;
### Solution:
y = function(x, k, start=1) {
	# start: default = 1 functions for 0 < x <= sqrt(7) & k == 1;
	# start = 500: for x <= 9.5 (k == 1/2);
	# root
	y.f = function(x, v) {
		r = if(x == 0) -v^2 else if(x < 0) -Inf else x*sin(k * log(x)) - v^2;
		return(r)
	}
	dy.f = function(x, v) {
		if(x == 0) return(sign(k));
		if(x < 0) return(1E+3);
		x.log = k * log(x)
		k*cos(x.log) + sin(x.log);
	}
	# damped waves;
	x0.f = function(x) {
		x2 = x^2
		x0 = if(x2 >= 0 & x2 <= 7.46) start else if(x2 < 0 & x2 > -0.32) 1/4 else 500; # for k == 1;
		return(x0);
	}
	y = sapply(x, function(x) newtonRaphson(y.f, x0.f(x), dfun=dy.f, v=x)[[1]])
	y = sapply(y, round0)
	return(y)
}
dy = function(x, k, y.x) {
	if(missing(y.x)) y.x = y(x, k=k);
	y.log = ifelse(x == 0, -1, log(y.x))
	div = x^2 + k * y.x * cos(k * y.log)
	dp = 2*x*y.x;
	dp = ifelse(div != 0, dp / div, 0); # may need correction
	return(dp)
}
d2y = function(x, k) {
	y.x = y(x, k=k)
	dy.x = dy(x, k=k, y.x=y.x)
	div = - x*y.x^2
	dp  = x*y.x*dy.x^2 - (k^2+1)/2 * x^2*dy.x^3 - y.x^2*dy.x
	dp = ifelse(div != 0, dp / div, 0); # TODO: needs correction!
	return(dp)
}
### Plot:
k = 1;
curve(y(x, k=k), from= 0+1E-3, to = sqrt(7))
# oscillating function with local minimum;
sapply(c(1/3, 1.1, 1.6, 1.8, 1.95), line.tan, dx=3, p=y, dp=dy, k=k)
# log-like
curve(dy(x, k=k), add=T, col="green")
sapply(c(1/3, 1.1, 1.7, 2, 2.5), line.tan, dx=3, p=dy, dp=d2y, k=k, col="orange")

# separately D2:
curve(dy(x, k=k), from= 0+1E-3, to = sqrt(7), col="green")
sapply(c(1/3, 1.1, 1.7, 2, 2.5), line.tan, dx=3, p=dy, dp=d2y, k=k, col="orange")


### k == 1/2; # TODO: may need correcting y0();
k = 1/2;
curve(y(x, k=k), from= 0+1E-3, to = 2.7)
# for x > 2.7 up to 9.5: needs start=500;
# oscillating function with local minimum;
sapply(c(1/3, 1.1, 1.6, 1.8, 1.95), line.tan, dx=3, p=y, dp=dy, k=k)
# log-like
curve(dy(x, k=k), add=T, col="green")
sapply(c(1/3, 1/2, 1, 2.2), line.tan, dx=3, p=dy, dp=d2y, k=k, col="orange")


# separately D2:
curve(dy(x, k=k), from= 0+1E-3, to = 2.7, col="green")
sapply(c(1/3, 1/2, 1, 2.2), line.tan, dx=3, p=dy, dp=d2y, k=k, col="orange")


############################
############################
############################

#####################
### Linear Simple ###
###   (Basics)    ###
#####################

### Simple / Power
# - generalization in the next sections;

### y = sin(a*x^n)
dy = n*a*x^(n-1)*cos(a*x^n)
x*d2y - (n-1)*dy + n^2*a^2*x^(2*n-1)*y # = 0

### Test & Plot:
y = function(x, a=1, n=2) {
	r = sin(a*x^n)
	return(r)
}
dy = function(x, a=1, n=2) {
	xn1 = if(n == 2) x else x^(n-1);
	xn  = xn1 * x;
	dp = n*a*xn1*cos(a*xn)
	return(dp)
}
d2y = function(x, a=1, n=2) {
	y.x = y(x, a=a, n=n)
	dy.x = dy(x, a=a, n=n)
	div = x;
	dp = (n-1)*dy.x - n^2*a^2*x^(2*n-1)*y.x;
	dp = ifelse(div != 0, dp / div, n*(n-1)*a); # TODO: needs correction!
	return(dp)
}
### Plot
a = 1; n = 2;
curve(y(x, a=a, n=n), from= -3, to= 3, ylim=c(-2, 1.5))
# sinus wave;
sapply(c((-5:5)/2.2), line.tan, dx=3, p=y, dp=dy, a=a, n=n)
# wave
curve(dy(x, a=a, n=n), add=T, col="green")
sapply(c((-5:5)/2.2), line.tan, dx=1/5, p=dy, dp=d2y, a=a, n=n, col="orange")

### Test wave
curve(dy(x, a=a, n=n), from=-3, to=3, col="green")
sapply(c((-5:5)/2.2), line.tan, dx=1/5, p=dy, dp=d2y, a=a, n=n, col="orange")


### Example 2:
a = 1/3; n = 2;
curve(y(x, a=a, n=n), from= -3, to= 3, ylim=c(-1, 1.5))
# sinus wave;
sapply(c((-5:5)/2.2), line.tan, dx=3, p=y, dp=dy, a=a, n=n)
# wave
curve(dy(x, a=a, n=n), add=T, col="green")
sapply(c((-5:5)/2.2), line.tan, dx=1/5, p=dy, dp=d2y, a=a, n=n, col="orange")

#####################

##################
### Non-Linear ###
###   Simple   ###
##################

### Generalization:

### y = sin(f(x))
# [but still basics: 1 trigonometric term]
# [not run]
dy = df * cos(f)
dy^2 = df^2 * cos(f)^2
dy^2 = df^2 * (1 - y^2)
### Non-linear ODE:
dy^2 + df^2 * y^2 - df^2 # = 0


### Alternative Linear ODE:
dP*d2y - d2P*dy + dP^3 * y # = 0
# [see next Section]

### Examples:

### f = x + ln(x)
dy^2 + ((x+1)/x)^2 * y^2 - ((x+1)/x)^2 # = 0
x^2 * dy^2 + (x+1)^2 * y^2 - (x+1)^2 # = 0
### Solution & Plot:
y = function(x) {
	f = x + log(x)
	r = sin(f)
	return(r)
}
dy = function(x) {
	# Test formula
	y.x = y(x)
	x1 = (x+1)^2
	dp = x1 * (1 - y.x^2)
	div = x^2
	dp = ifelse(div != 0, dp / div, 1E+3); # TODO: check!
	f.sign = ifelse(x + log(x) <= pi/2, FALSE, TRUE)
	dp = sqrt(dp); dp[f.sign] = - dp[f.sign];
	return(dp)
}
### Plot:
curve(y(x), from= 0+1E-3, to=3, ylim=c(-1.5, 2))
sapply(c((1:8)/3.2), line.tan, dx=3, p=y, dp=dy)


### D2:
# dy^2 + df^2 * y^2 - df^2 # = 0
2*dy*d2y + 2*df^2 * y*dy + 2*df*d2f * y^2 - 2*df*d2f # = 0
dy*d2y + df^2 * y*dy + df*d2f * y^2 - df*d2f # = 0
# alternative:
dy*d2y + df^2 * y*dy + d2f * (df^2 - dy^2)/df - df*d2f # = 0
dy*d2y - d2f/df * dy^2 + df^2 * y*dy # = 0

### Examples:

### f = x + ln(x)
dy*d2y + (1/x^2)/((x+1)/x) * dy^2 + ((x+1)/x)^2 * y*dy # = 0
dy*d2y + 1/((x+1)*x) * dy^2 + ((x+1)/x)^2 * y*dy # = 0 # * x^2
x^2*dy*d2y + x/(x+1) * dy^2 + (x+1)^2 * y*dy # = 0
### Solution & Plot:
# reuses functions y(x) & dy(x) from above;
d2y = function(x) {
	y.x = y(x)
	dy.x = dy(x)
	dp = - (x/(x+1) * dy.x^2 + (x+1)^2 * y.x*dy.x);
	div = x^2 * dy.x;
	dp = ifelse(div != 0, dp / div, 1E+3); # TODO: check!
	return(dp)
}
### Plot:
curve(y(x), from= 0+1E-3, to=3, ylim=c(-1.5, 3))
sapply(c((1:5)/2.2), line.tan, dx=3, p=y, dp=dy)
# wave
curve(dy(x), add=T, col="green")
sapply(c((1:5)/2.2), line.tan, dx=1/5, p=dy, dp=d2y, col="orange")


#####################
#####################

#####################
### Linear Simple ###
###  (Polynomial) ###
#####################

### y = a1*sin(P(x)) + a2*cos(P(x))
# - both sin(P(x)) & cos(P(x)) are solutions;

### D =>
# dy = dP * (a1*cos(P(x)) - a2*sin(P(x)))
### D2 =>
# d2y = d2P*(a1*cos(P(x)) - a2*sin(P(x))) - dP^2 * y
### Solve for sin(P(x)), cos(P(x)) using y & dy;
# dP * sin(P(x)) = (-a2*dy + a1*dP*y) / (a1^2 + a2^2)
# dP * cos(P(x)) = ( a1*dy + a2*dP*y) / (a1^2 + a2^2)
### =>
# (a1^2 + a2^2)*dP*d2y = d2P*(a1*(a1*dy + a2*dP*y) + a2*(a2*dy - a1*dP*y)) - (a1^2 + a2^2)*dP^3 * y
dP*d2y - d2P*dy + dP^3 * y # = 0

### Examples:
### P(x) = x^2
x*d2y - dy + 4*x^2 * y # = 0
### P(x) = x^n => * x^(2-n)
x*d2y - (n-1)*dy + n^2*x^(2*n-1) * y # = 0
### P(x) = x^n + b*x
(n*x^(n-1)+b)*d2y - n*(n-1)*x^(n-2)*dy + (n*x^(n-1)+b)^3 * y # = 0
# [is detailed below]

### Example: P(x) = x^n + b*x;

### y = a1*sin(x^n + b*x) + a2*cos(x^n + b*x)
# dy = - (n*x^(n-1) + b) * (a2*sin(x^n + b*x) - a1*cos(x^n + b*x))
# sin(x^n + b*x) = - (a2*dy - a1*(n*x^(n-1) + b)*y) / (a1^2 + a2^2) / (n*x^(n-1) + b)
# cos(x^n + b*x) =   (a1*dy + a2*(n*x^(n-1) + b)*y) / (a1^2 + a2^2) / (n*x^(n-1) + b)
d2y = - n*(n-1)*x^(n-2) * (a2*sin(x^n + b*x) - a1*cos(x^n + b*x)) +
  - (n*x^(n-1) + b)^2 * (a1*sin(x^n + b*x) + a2*cos(x^n + b*x))
= - n*(n-1)*x^(n-2) * (a[2]*sin(x^n + b*x) - a[1]*cos(x^n + b*x)) - (n*x^(n-1) + b)^2 * y
(a1^2 + a2^2)*(n*x^(n-1) + b) * d2y +
  + n*(n-1)*x^(n-2)*( - a[2]*(a[2]*dy(x) - a[1]*(n*x^(n-1) + b)*y(x)) - a[1]*(a[1]*dy(x) + a[2]*(n*x^(n-1) + b)*y(x))) +
  + (a[1]^2 + a[2]^2)*(n*x^(n-1) + b)^3 * y(x)

### ODE:
(n*x^(n-1) + b) * d2y - n*(n-1)*x^(n-2)*dy(x) + (n*x^(n-1) + b)^3 * y(x) # = 0
### Example: n = 2
(2*x + b) * d2y(x, n=2) - 2*dy(x, n=2) + (2*x + b)^3 * y(x, n=2) # = 0
### Example: n = 3
(3*x^2 + b) * d2y(x, n=3) - 6*x*dy(x, n=3) + (3*x^2 + b)^3 * y(x, n=3) # = 0

### Solution:
y = function(x, a=c(1, 1), b=1, n=2) {
	xn = x^n + b*x
	r = a[1]*sin(xn) + a[2]*cos(xn)
	return(r)
}
dy = function(x, a=c(1, 1), b=1, n=2) {
	xn1 = x^(n-1)
	xn = (xn1 + b)*x
	r = - (n*xn1 + b) * (a[2]*sin(xn) - a[1]*cos(xn))
	r = round0(r)
	return(r)
}
d2y = function(x, a=c(1, 1), b=1, n=2) {
	y.x = y(x, a=a, b=b, n=n)
	dy.x = dy(x, a=a, b=b, n=n)
	xn2 = n * x^(n-2)
	xnb = xn2 * x + b
	#
	div = xnb;
	dp  = (n-1)*xn2*dy.x - xnb^3*y.x
	dp = ifelse(div != 0, dp / div, n*(n-1)); # TODO: needs correction!
	return(dp)
}
### Plot:
n = 2; b = 1;
curve(y(x, b=b, n=n), from= -3, to = 3, ylim=c(-3, 3))
# oscillating function with local minima;
# slightly shifted: + 1/2;
# dy == 0 also for: x^2 + b*x - pi/4 = 0;
sapply(c(-5:7 * 3/7 - 1/2), line.tan, dx=1.5, p=y, dp=dy, b=b, n=n)
# also sinusoidal:
curve(dy(x, b=b, n=n), add=T, col="green")
sapply(c(-5:7 * 3/7 - 1/2), line.tan, dx=1.5, p=dy, dp=d2y, b=b, n=n, col="orange")


### Ex 2:
n = 2; b = -2;
curve(y(x, b=b, n=n), from= -3, to = 3, ylim=c(-4, 4))
# oscillating function with local minima;
# slightly shifted: + 1/2
sapply(c(-5:7 * 3/7 - 1/2), line.tan, dx=1.5, p=y, dp=dy, b=b, n=n)
# also sinusoidal:
curve(dy(x, b=b, n=n), add=T, col="green")
sapply(c(-5:7 * 3/7 - 1/2), line.tan, dx=1.5, p=dy, dp=d2y, b=b, n=n, col="orange")


########################

########################
###  Linear Complex  ###
### (Non-Polynomial) ###
########################

###########
### LOG ###

### y = a1*sin(log(P(x))) + a2*cos(log(P(x)))

### D =>
# P(x) * dy = dP * (a1*cos(log(P(x))) - a2*sin(log(P(x))))
### D2 * P(x) =>
# P(x)^2 * d2y + P(x) * dP*dy = P(x)*d2P(x)*(a1*cos(log(P(x))) - a2*sin(log(P(x)))) - dP(x)^2 * y
### Solve for sin(log(P(x))), cos(log(P(x))) using y & dy;
# sin(log(P(x))) = (-a2*P*dy + a1*dP*y) / ((a1^2 + a2^2)*dP)
# cos(log(P(x))) = ( a1*P*dy + a2*dP*y) / ((a1^2 + a2^2)*dP)
(a1^2 + a2^2)*dP*(P(x)^2 * d2y + P(x)*dP*dy) =
	P(x)*d2P(x)*(a[1]*(a[1]*P(x)*dy + a[2]*dP(x)*y) + a[2]*(a[2]*P(x)*dy - a[1]*dP(x)*y)) - (a1^2 + a2^2)*dP(x)^3 * y
### Eq:
dP*P(x)^2 * d2y + P(x)*dP^2*dy - P(x)^2*d2P*dy + dP^3*y # = 0

### Examples:
### P(x) = x + b
(x+b)^2 * d2y + (x+b)*dy + y # = 0
### P(x) = x^2 + b1*x + b0
(2*x+b1)*(x^2+b1*x+b0)^2 * d2y + (x^2+b1*x+b0)*(2*x+b1)^2*dy - 2*(x^2+b1*x+b0)^2*dy + (2*x+b1)^3*y # = 0
### P(x) = x^n + b0
n*x^(n-1)*(x^n+b0)^2 * d2y + n^2*(x^n+b0)*x^(2*n-2)*dy - n*(n-1)*(x^n+b0)^2*x^(n-2)*dy + n^3*x^(3*n-3)*y # = 0 # * x^(2-n) / n
x*(x^n+b0)^2 * d2y + n*x^n*(x^n+b0)*dy - (n-1)*(x^n+b0)^2*dy + n^2*x^(2*n-1)*y # = 0
x*(x^n+b0)^2 * d2y + (x^(2*n) - (n-2)*b0*x^n - (n-1)*b0^2)*dy + n^2*x^(2*n-1)*y # = 0


### Solution:
y = function(x, a=c(1, 1), FUN.list, ...) {
	xlog = log(FUN.list[[1]](x, ...));
	r = a[1]*sin(xlog) + a[2]*cos(xlog)
	return(r)
}
dy = function(x, a=c(1, 1), FUN.list, ...) {
	div = FUN.list[[1]](x, ...)
	xlog = log(div);
	dp = (a[1]*cos(xlog) - a[2]*sin(xlog)) * FUN.list[[2]](x, ...)
	dp = ifelse(div != 0, dp / div, 0); # TODO: check;
	dp = round0(dp)
	return(dp)
}
d2y = function(x, a=c(1, 1), FUN.list, ...) {
	y.x  =  y(x, a=a, FUN.list=FUN.list, ...)
	dy.x = dy(x, a=a, FUN.list=FUN.list, ...)
	Px = FUN.list[[1]](x, ...)
	dP = FUN.list[[2]](x, ...); dPsq = dP^2;
	d2P = FUN.list[[3]](x, ...)
	#
	div = - dP * Px^2;
	dp  = Px*(dPsq - Px*d2P)*dy.x + dPsq*dP*y.x
	dp = ifelse(div != 0, dp / div, d2P/Px*(a[1]*cos(log(Px)) - a[2]*sin(log(Px)))); # TODO: check!
	return(dp)
}
P = function(x, b=c(1, 1)) {
	x^2 + b[2]*x + b[1]
}
dP = function(x, b=c(1, 1)) {
	2*x + b[2]
}
d2P = function(x, b=c(1, 1)) {
	2
}
p.list = list(P, dP, d2P)
### Plot:
b = c(1, 1);
curve(y(x, b=b, FUN.list=p.list), from= -3, to = 3, ylim=c(-2, 2))
# oscillating function with local minima;
# slightly shifted: + 1/2
sapply(c(-5:7 * 3/7 - 1/2), FUN=line.tan, dx=2, p=y, dp=dy, FUN.list=p.list, b=b)
# also sinusoidal:
curve(dy(x, b=b, FUN.list=p.list), add=T, col="green")
sapply(c(-5:7 * 3/7 - 1/2), line.tan, dx=1.5, p=dy, dp=d2y, FUN.list=p.list, b=b, col="orange")

###
b = c(3, -1);
curve(y(x, b=b, FUN.list=p.list), from= -3, to = 3, ylim=c(-1, 1.5))
# possible local maximum;
# slightly shifted: ???
sapply(c(-5:7 * 3/7 - 1/2), FUN=line.tan, dx=2, p=y, dp=dy, FUN.list=p.list, b=b)
#
curve(dy(x, b=b, FUN.list=p.list), add=T, col="green")
sapply(c(-5:7 * 3/7 - 1/2), line.tan, dx=1.5, p=dy, dp=d2y, FUN.list=p.list, b=b, col="orange")


############################
############################

### Combined Functions:
### Linearly combined


### f1(x)*sin(y) + f2(x)*cos(y) = P(x)

f2(x)*cos(y) + f1(x)*sin(y) # = p
### D =>
(f1*dy + df2)*cos(y) - (f2*dy - df1)*sin(y) - dp # = 0
### Solve linear =>
sin(y) = ((f1*dy + df2)*p - f2*dp) / ((f1*dy + df2)*f1 + (f2*dy - df1)*f2)
cos(y) = ((f2*dy - df1)*p + f1*dp) / ((f2*dy - df1)*f2 + (f1*dy + df2)*f1)

### D2 =>
(f1*d2y + df1*dy + d2f2)*cos(y) - (f1*dy + df2)*dy*sin(y) +
	- (f2*d2y + df2*dy - d2f1)*sin(y) - (f2*dy - df1)*dy*cos(y) - d2p # = 0
(f1*d2y - f2*dy^2 + 2*df1*dy + d2f2)*cos(y) +
	- (f2*d2y + f1*dy^2 + 2*df2*dy - d2f1)*sin(y) - d2p # = 0

### ODE:
(f1*d2y - f2*dy^2 + 2*df1*dy + d2f2)*((f2*dy - df1)*p + f1*dp) +
	- (f2*d2y + f1*dy^2 + 2*df2*dy - d2f1)*((f1*dy + df2)*p - f2*dp) +
	- d2p*((f1*dy + df2)*f1 + (f2*dy - df1)*f2) # = 0

### Examples:

### f1(x) = x; f2(x) = 1;
(x*d2y - dy^2 + 2*dy)*(p*dy - p + x*dp) - (d2y + x*dy^2)*(x*p*dy - dp) +
	- d2p*(x^2*dy + dy - 1) # = 0
(dp*x^2 - x*p + dp)*d2y - p*(x^2 + 1)*dy^3 + 3*p*dy^2 - 2*(x^2*d2p - x*dp + p + d2p)*dy + d2p # = 0
### P(x) = x^2
(x^3 + 2*x)*d2y - x^2*(x^2 + 1)*dy^3 + 3*x^2*dy^2 - 2*(x^2 + 2)*dy + 2 # = 0

### TODO: check!


#########################
#########################

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
n*x^(n-1) - n^2*x^(2*n-2)*(sin(x^n) * Ic - cos(x^n) * Is)
n*x^(n-1) - n^2*x^(2*n-2)*y
# TODO: check;

### ODE:
d2y + n^2*x^(2*n-2)*y - n*x^(n-1) # = 0