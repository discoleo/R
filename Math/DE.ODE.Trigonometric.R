########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Trigonometric
###
### draft v.0.1c


### History

### Order 1 Non-Liniar:
### Trigonometric Variants
###
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

### Simple Trigonometric Functions:


### sin(y^2) = f(x);
df*(1 - f^2)*y^3*d2y - (d2f - d2f*f^2 + df^2*f)*y^3*dy + 1/4 * df^3 # = 0

### Examples:
# [not run]
# f = x^2;
x*(1-x^4)*y^3*d2y - (x^4+1)*y^3*dy + x^3 # = 0
# f = x^2 + b;
x*(1-(x^2+b)^2)*y^3*d2y - (x^4 - b^2 + 1)*y^3*dy + x^3 # = 0
### Solution:
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
	dp =(x4 - b^2 + 1)*y.x^3*z - x^3;
	dp = ifelse(div != 0, dp / div, -1); # TODO: needs correction!
	return(dp)
}
### b = 0;
curve(y(x), from= -1, to = 1)
# global minimum;
sapply(c((-2:2)/2.2), line.tan, dx=3, p=y, dp=dy)
# pseudo-sigmoidal
curve(dy(x), from= -1, to = 1, add=T, col="green")
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
curve(dy(x, b=b), from= -1, to = 1, add=T, col="green")
sapply(c(-2, -1.5, 1.5, 2)/2.1, line.tan, dx=1/5, p=dy, dp=d2y, b=b, col="orange")

# check full pseudo-sigmoidal
curve(dy(x, b=b), from= -sqrt(2) + 1E-10, to = sqrt(2) - 1E-10, col="green", ylim=c(-4,4))
sapply(c(-2.8, -2.4, -2, -1.5, 1.5, 2, 2.4, 2.8)/2.1, line.tan, dx=1/5, p=dy, dp=d2y, b=b, col="orange")


### b = -1/2;
b = -1/2
curve(dy(x, b=b), from= -sqrt(5/4), to = sqrt(5/4), col="green", ylim=c(-4,4))
sapply(c(-2.4, -2, -1.5, 1.5, 2, 2.4)/2.1, line.tan, dx=1/5, p=dy, dp=d2y, b=b, col="orange")



############################

### y*sin(y) + cos(y) = f(x)

### D =>
y*dy*cos(y) - df # = 0
# D2 =>
df*y*d2y + 2*df*(dy)^2 - f*y*(dy)^3 - d2f*y*dy # = 0

### Examples:
# f = x + b, where b = constant;
y*d2y + 2*(dy)^2 - (x + b)*y*(dy)^3 # = 0
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

