########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Trigonometric
###
### draft v.0.1f


### History

### Order 1 Non-Liniar:
### Trigonometric Variants
###
### draft v.0.1e - v.0.1f:
# - solved:
#   x^2*dy*d2y + x/(x+1) * dy^2 + (x+1)^2 * y*dy = 0;
#   [includes generalization]
#   x*y^2*d2y + x*y*dy^2 - (k^2+1)/2 * x^2*dy^3 - y^2*dy = 0;
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
curve(dy(x, b=b), add=T, col="green")
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


######################

### Simple / Power

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


##########################
### Simple: Generalization

### y = sin(f(x))
dy = df * cos(f)
dy^2 = df^2 * cos(f)^2
dy^2 = df^2 * (1 - y^2)
dy^2 + df^2 * y^2 - df^2 # = 0

### Examples:

### f = x + ln(x)
dy^2 + ((x+1)/x)^2 * y^2 -  ((x+1)/x)^2 # = 0
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




########################
### Integration by parts

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

### y*sin(k * log(y)) = f(x)

### y*sin(k * log(y)) = x^2
# k*dy * cos(k * log(y)) = (2*x*y - x^2*dy) / y;
x*y^2*d2y + x*y*dy^2 - (k^2+1)/2 * x^2*dy^3 - y^2*dy # = 0;
### Solution & Plot:
y = function(x, k) {
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
		x0 = if(x2 >= 0 & x2 <= 7.46) 1 else if(x2 < 0 & x2 > -0.32) 1/4 else 500; # for k == 1;
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
# spikes
curve(dy(x, k=k), add=T, col="green")
sapply(c(1/3, 1.1, 1.7, 2, 2.5), line.tan, dx=3, p=dy, dp=d2y, k=k, col="orange")

# separately D2:
curve(dy(x, k=k), from= 0+1E-3, to = sqrt(7), col="green")
sapply(c(1/3, 1.1, 1.7, 2, 2.5), line.tan, dx=3, p=dy, dp=d2y, k=k, col="orange")


### k == 1/2; # TODO: may need correcting y0();
k = 1/2;
curve(y(x, k=k), from= 0+1E-3, to = 2.7)
# oscillating function with local minimum;
sapply(c(1/3, 1.1, 1.6, 1.8, 1.95), line.tan, dx=3, p=y, dp=dy, k=k)
# spikes
curve(dy(x, k=k), add=T, col="green")
sapply(c(1/3, 1/2, 1, 2.2), line.tan, dx=3, p=dy, dp=d2y, k=k, col="orange")

# separately D2:
curve(dy(x, k=k), from= 0+1E-3, to = 2.7, col="green")
sapply(c(1/3, 1/2, 1, 2.2), line.tan, dx=3, p=dy, dp=d2y, k=k, col="orange")


