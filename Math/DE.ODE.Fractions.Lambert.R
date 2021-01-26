
########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Fractions: Lambert
###
### draft v.0.3e


### History

### Order 1 Non-Liniar
###
### draft v.0.3e:
# - ODEs derived from:
#   y^n = log(f(x)) * log(g(x));
# - simple example from: y = log(x + k) * log(x - k);
### draft v.0.3d:
# - ODEs derived from e^e^f(y, x);
### draft v.0.3b - v.0.3c-ex:
# - [Technique] Integration by parts:
#   dz^2 + 2*(a*x+1)*dz - 2*a*z + 2*a*x = 0;
#   I(y^2) = f(y, x); [v.0.3c & v.0.3c-ex]
#   a/2*I(y^2) = 1/3*y^3 + 1/2*a*x*y^2 + y^2 + (a*x+1)*y + a*x;
#   TODO:
#   y*dy + (a*x + 1)*dy + 2*a*x*y + 2*a^2*x^2 + a = 0;
### draft v.0.3a: [16-11-2020]
# - added some theoretical aspects:
#   f*g*Dg + f*Dg - df*g = 0;
#   where g = g(y), f = f(x);
# - examples (trivial, but can be extended as in ex.2.):
#   n*f*y^n*dy + n*f*dy - df*y = 0;
#   x*y*dy + x*(a*x + 1)*dy - 2*y - a*x = 0;
### draft v.0.2b: [12-11-2020]
# - cleanup: moved Gaussian type to
#   new file: DE.ODE.Gaussian.R;
### draft v.0.2a - v.0.2a-xTerm: [11-11-2020]
# - solved [& moved to new file]:
#   d2z + 2*x*dz - 2*z = 0;
#   d2z + 2*x*dz - 4*z = 0; [v.0.2a-x2Lev3]
#   d2z + 2*x*dz - 6*z = 0; [v.0.2a-x2Lev4]
#   d2z + 2*x*dz - 2*z = b0*x^2; [v.0.2a-xTerm]
#   d3z + 3*x^2*d2z - 6*x*dz + 6*z = 0; [v.0.2a-pow3]
### draft v.0.2-pre-a:
# - moved Trigonometric variants to new file:
#   DE.ODE.Trigonometric.R;
### draft v.0.1h:
# - solved a trigonometric type:
#   x*(1-x^4)*y^3*d2y - (x^4+1)*y^3*dy + x^3 = 0;
# - [TODO] move to separate file; [DONE]
### draft v.0.1g - v.0.1g-ln: [07-11-2020]
# - solved:
#   y^2*dy + c*y*dy - c*x^2*dy - 2*x*y^2 = 0; [v.0.1g]
#   x^2*dy + n*y^2 - x*y = 0; [v.0.1g-ln]
#   where c, n = constants;
### draft v.0.1f: [05-11-2020]
# - solved: x*dz - x^2*z^3 + 2*x*z^2 - z = 0;
### draft v.0.1e - v.0.1e-more: [05-11-2020]
# - started to explore various y powers:
#   x*y^(1/2)*dy + x*dy - 2*y = 0;
#   x*sqrt(y + x + 1)*dy + x*dy - 2*y + x*sqrt(y + x + 1) = x + 2; (v.0.1e-more)
# - find bug; [FIXED] (v.0.1e-fix)
### draft v.0.1d: [02-11-2020]
# - added variants based on logarythms (equivalent to y^y):
#   y*dy + (x+b)*dy - y = 0;
### draft v.0.1c:
# - 3-parameter generalization of ODEs of type:
#   y*dy + dy + y + f(x) = 0;
### draft v.0.1b - v.0.1b-2:
# - trigonometric coefficients;
# - more examples (but with same structure of trig coeffs);
# - minor improvements: using ifelse() (v.0.1b-2);
### draft v.0.1a:
# - moved section Exponential/Lambert
#   to this new file;
# - renamed to Fractions: Lambert;

### [old file] DE.ODE.Polynomial.R
### draft v.0.1f - v.0.1f-2:
# - solved: x*y*dy + x*(x+1)*dy - 2*y = x;
#   (a Lambert snack)
# - various generalizations, e.g.:
#   (x^2+b)*y*dy + (x^2+b)*(x+1)*dy - 2*x*y = x^2 - b; (v.0.1f-2)


#########################

### Helper functions

library(pracma)
# needed for Lambert W;


# include: DE.ODE.Helper.R;
source("DE.ODE.Helper.R")

#########################

#########################
### W Exponentials =>
### Fractions & Lambert W

### Theory

### g(y) * e^(g(y)) = f(x)

### D =>
g*e^g*Dg + e^g*Dg - df # = 0 # * g
g^2*e^g*Dg + g*e^g*Dg - df*g # = 0
# g*e^g = f =>
f*g*Dg + f*Dg - df*g # = 0

### Examples:
# g = y:
f*y*dy + f*dy - df*y # = 0
# g = y^2
2*f*y^3*dy + 2*f*y*dy - df*y^2 # = 0
2*f*y^2*dy + 2*f*dy - df*y # = 0
# g = y^n
n*f*y^(2*n-1)*dy + n*f*y^(n-1)*dy - df*y^n # = 0
n*f*y^n*dy + n*f*dy - df*y # = 0

# g = y + a*x; a = constant;
f*(y + a*x)*(dy + a) + f*(dy + a) - df*(y + a*x) # = 0
f*y*dy + f*(a*x + 1)*dy + (a*f - df)*y + a^2*f*x - a*x*df + a*f # = 0
# for f = x^2 * e^(a*x):
# e^(a*x) can be factored both from f & df, as df = (a*x^2 + 2*x)*e^(ax):
x^2*y*dy + x^2*(a*x + 1)*dy + (a*x^2 - a*x^2 - 2*x)*y + a^2*x^3 - a*x*(a*x^2 + 2*x) + a*x^2 # = 0
x*y*dy + x*(a*x + 1)*dy - 2*y - a*x # = 0


################
### Examples ###

### (x + y)*e^y = x^2
x*y*dy + x*(x+1)*dy - 2*y - x # = 0
y*dy + (x+1)*dy - 2/x * y - 1 # = 0
### Solution:
y = function(x, a=1) {
	# root
	y = lambertWp(x^2 * exp(a*x)) - a*x
	y = sapply(y, round0)
	return(y)
}
dy = function(x, a=1) {
	y.x = y(x, a=a)
	ax = a*x;
	div = x*y.x + x*(ax+1)
	dp = 2*y.x + ax;
	# TODO: find BUG: must be -1 (but why?);
	dp = ifelse(div != 0, dp / div, -a);
	return(dp)
}
###
curve(y(x), from=-1/5, to=3)
# a nice global minimum
sapply(c((0:4)/2), line.tan, dx=3, p=y, dp=dy)

###
curve(y(x), from=-1/5, to=5)
# a nice global minimum
sapply(c((0:3)*4/3), line.tan, dx=3, p=y, dp=dy)

### a == 2
a = 2;
curve(y(x, a=a), from=-1/5, to=3)
# a nice global minimum
sapply(c(0:4, 7.5, 10)/5, line.tan, dx=3, p=y, dp=dy, a=a)

### a == 1/2
a = 1/2;
curve(y(x, a=a), from=-1/5, to=3)
# a nice global minimum
sapply(c(0:4, 10)/5, line.tan, dx=3, p=y, dp=dy, a=a)



#####################
### (x + y)*e^y = x^2 + b
# [not run]
(x^2+b)*y*dy + (x^2+b)*(x+1)*dy - 2*x*y - x^2 + b # = 0
### Solution:
y = function(x, b) {
	# root
	y = lambertWp((x^2+b) * exp(x)) - x
	y = sapply(y, round0)
	return(y)
}
dy = function(x, b) {
	y.x = y(x, b)
	x2 = x^2 + b
	div = x2*(y.x + x + 1)
	dp = 2*x*y.x + x2 - 2*b # using x2 instead of x^2
	dp = ifelse(div != 0, dp / div, -1);
	return(dp)
}
curve(y(x, b=1/2), from=-1, to=3)
# a nice global minimum
sapply(c((0:4)/2), line.tan, dx=3, p=y, dp=dy, b=1/2)


#####################
### (x + a + y)*e^y = x^2 + b
# [not run]
(1 + dy)*e^y + (x + a + y)*e^y*dy = 2*x
(1 + dy)*(x^2 + b)/(y + x + a) + (x^2 + b)*dy - 2*x # = 0
(1 + dy)*(x^2 + b) + (x^2 + b)*(y + x + a)*dy - 2*x*(y + x + a) # = 0
(x^2 + b)*y*dy + (x + a + 1)*(x^2 + b)*dy - 2*x*y - x^2 - 2*a*x + b # = 0
###
(x^2 + b)*y*dy + (x + a + 1)*(x^2 + b)*dy - 2*x*y - x^2 - 2*a*x + b # = 0
### Solution:
y = function(x, a, b) {
	# root
	y = lambertWp((x^2+b) * exp(x+a)) - x - a
	y = sapply(y, round0)
	return(y)
}
dy = function(x, a, b) {
	y.x = y(x, a, b)
	x2 = x^2 + b
	div = x2 * (y.x + x + a + 1)
	dp = 2*x*y.x + x^2 + 2*a*x - b
	dp = ifelse(div != 0, dp / div, -1); # may need correction
	return(dp)
}
curve(y(x, a=-1, b=1/2), from=-2, to=3)
# a nice global minimum
sapply(c(-1, (0:4)/3), line.tan, dx=3, p=y, dp=dy, a=-1, b=1/2)


###################
### Generalisations

### (x + y)*e^y = x^n + b*x
(x^n+b*x)*y*dy + (x^n+b*x)*(x+1)*dy - (n*x^(n-1) + b)*y - (n-1)*x^n # = 0
(x^3+x)*y*dy + (x^3+x)*(x+1)*dy - (3*x^2 + 1)*y - 2*x^3 # for: n = 3; b = 1;
### Solution:
y = function(x, n, b) {
	# root
	y = lambertWp((x^n + b*x) * exp(x)) - x
	y = sapply(y, round0)
	return(y)
}
dy = function(x, n, b) {
	y.x = y(x, n, b)
	dp = (n*x^(n-1) + b)*y.x + (n-1)*x^n;
	xn = x^n + b*x
	div = xn*y.x + xn*(x+1)
	# TODO: correct Limit;
	dp = ifelse(div != 0, dp / div, 0);
	return(dp)
}
curve(y(x, n=3, b=1), from=-1/5, to=3)
# a nice global minimum
sapply(c((0:4)/2), line.tan, dx=3, p=y, dp=dy, n=3, b=1)


### Generalisations

### (y + ax)*e^(y + ex) = fx
(dy + dax)*e^(y + ex) + (dy + dex)*(y + ax)*e^(y + ex) - dfx # = 0
(dy + dax)*fx/(y + ax) + (dy + dex)*fx - dfx # = 0 # * (y + ax)
(y + ax)*(dy + dex)*fx + (dy + dax)*fx - dfx*(y + ax) # = 0
fx*(y + ax)*(dy + dex) + fx*(dy + dax) - dfx*(y + ax) # = 0
fx*y*dy + fx*(ax + 1)*dy + (fx*dex - dfx)*y + ax*(fx*dex - dfx) + fx*dax # = 0

### Example 1:
# fx = x^2
# ax = x^2 + x - 1
# ex = -2*x
x^2*y*dy + x^2*(x^2 + x)*dy - 2*(x^2 + x)*y - 2*(x^2 + x - 1)*(x^2 + x) + x^2*(2*x+1) # = 0
y*dy + (x^2 + x)*dy - 2*(x + 1)/x * y - 2*(x^2 + x - 1)*(x + 1)/x + 2*x+1 # = 0
### Solution:
y = function(x) {
	# root
	ax = x^2 + x - 1
	y = lambertWp(x^2 * exp(ax + 2*x)) - ax
	y = sapply(y, round0)
	return(y)
}
dy = function(x) {
	y.x = y(x)
	ax = x^2 + x - 1
	dp = 2*(y.x + ax)*(x+1)/ x - 2*x - 1
	div = y.x + ax + 1;
	# TODO: correct Limit;
	dp = ifelse(div != 0, dp / div, 0);
	return(dp)
}
curve(y(x), from=-1, to=1)
# a nice global minimum
sapply(c((-4:4)/5), line.tan, dx=3, p=y, dp=dy)


#######################

#######################
### Different Power ###

### sqrt(y)*e^(sqrt(y)) = f(x)
# [not run]
1/2 * y^(-1/2)*e^(sqrt(y))*dy + 1/2 * y^(-1/2)*y^(1/2)*e^(sqrt(y))*dy - df # = 0
f/y * dy + f*y^(-1/2)*dy - 2*df # = 0 # * y
f*y^(1/2)*dy + f*dy - 2*df*y # = 0

### Example:
# f = x;
x*y^(1/2)*dy + x*dy - 2*y # = 0
### Solution:
y = function(x, useNeg=FALSE) {
	# root
	y = lambertWp(x)^2
	if(useNeg) {
		isNeg = (x < 0)
		if(any(isNeg)) {
			y[isNeg] = lambertWn(x[isNeg])^2
		}
	}
	y = sapply(y, round0)
	return(y)
}
dy = function(x, useNeg=FALSE) {
	y.x = y(x, useNeg=useNeg)
	dp = 2*y.x
	y.sq = sqrt(y.x)
	isNeg = lambertWp(x) < 0
	if(any(isNeg)) {
		y.sq[isNeg] = -y.sq[isNeg]
	}
	div = x*(y.sq + 1)
	dp = ifelse(div != 0, dp / div, 0);
	return(dp)
}
curve(y(x), from= -0.5, to=3)
# a nice global minimum
sapply(c(-1/3, -1/5, (0:4)/3), line.tan, dx=3, p=y, dp=dy)


##################

### sqrt(y + g(x))*e^(sqrt(y + g(x))) = f(x)
# [not run]
(dy + dg)/sqrt(y + g)*e^(sqrt(y + g)) + (dy + dg)/sqrt(y + g) * sqrt(y + g)*e^(sqrt(y + g)) - 2*df # = 0
f*(dy + dg)/(y + g) + f*(dy + dg)/sqrt(y + g) - 2*df # = 0 # *(y+g)
f*(dy + dg)*sqrt(y + g) + f*(dy + dg) - 2*df*(y + g) # = 0
f*(dy + dg)*sqrt(y + g) + f*dy - 2*df*y + f*dg - 2*df*g # = 0

### Example:
# f = x;
# g = x + 1
x*sqrt(y + x + 1)*dy + x*dy - 2*y + x*sqrt(y + x + 1) - x - 2 # = 0
### Solution:
y = function(x, useNeg=FALSE) {
	# root
	y = lambertWp(x)^2 - x - 1
	if(useNeg) {
		isNeg = (x < 0)
		if(any(isNeg)) {
			y[isNeg] = lambertWn(x[isNeg])^2 - x - 1
		}
	}
	y = sapply(y, round0)
	return(y)
}
dy = function(x, useNeg=FALSE) {
	y.x = y(x, useNeg=useNeg)
	y.sq = sqrt(y.x + x + 1)
	isNeg = lambertWp(x) < 0
	if(any(isNeg)) {
		y.sq[isNeg] = -y.sq[isNeg]
	}
	dp = 2*y.x - x*y.sq + x + 2
	div = x*(y.sq + 1)
	dp = ifelse(div != 0, dp / div, -1); # TODO: correct limit!
	return(dp)
}
curve(y(x), from= -0.5, to=3)
sapply(c(-1/3, -1/5, (0:4)/3), line.tan, dx=3, p=y, dp=dy)


##################

### c * e^(1/y) + y = f(x)
# where c = constant;
# [not run]
y^2*dy + c*y*dy - c*f*dy - df*y^2 # = 0

### Example:
# f = x^2
y^2*dy + c*y*dy - c*x^2*dy - 2*x*y^2 # = 0
### Solution
y = function(x, b) {
	# root
	y.f = function(x, v) b*exp(1/x) + x - v^2;
	dy.f = function(x, v) 1 - b*exp(1/x)/x^2;
	y = sapply(x, function(x) newtonRaphson(y.f, ifelse(x < 0, x, x - 1), dfun=dy.f, v=x)[[1]])
	y = sapply(y, round0)
	return(y)
}
dy = function(x, b) {
	y.x = y(x, b)
	dp = 2*x*y.x^2
	div = y.x^2 + b*y.x - b*x^2
	dp = ifelse(div != 0, dp / div, -1); # TODO: correct limit!
	return(dp)
}
###
c = 1;
curve(y(x, b=c), from= 1.9, to=4)
sapply(c(1.9, (4:6)/2), line.tan, dx=3, p=y, dp=dy, b=c)

##############
### Simple ###

### e^(x/y) = x^n
x^2*dy + n*y^2 - x*y # = 0
### Solution:
y = function(x, n) {
	# root
	y = x
	if(is.integer(n) && n %% 2 == 0) x = abs(x)
	y = y / (n*log(x))
	y = sapply(y, round0)
	return(y)
}
dy = function(x, n) {
	y.x = y(x, n)
	dp = - y.x*(n*y.x - x)
	div = x^2
	dp = ifelse(div != 0, dp / div, 0); # TODO: correct limit!
	return(dp)
}
###
n = 4L;
curve(y(x, n=n), from= -2, to=2)
sapply(c((-5:5)/2.8), line.tan, dx=3, p=y, dp=dy, n=n)



##################################

##################################
### Trigonometric Coefficients ###

#####################
### (x + a + y)*e^y = sin(x)^2 + b
# [not run]
(1 + dy)*e^y + (x + a + y)*e^y*dy = 2*sin(x)*cos(x)
(1 + dy)*(sin(x)^2 + b)/(y + x + a) + (sin(x)^2 + b)*dy - sin(2*x) # = 0
(1 + dy)*(sin(x)^2 + b) + (sin(x)^2 + b)*(y + x + a)*dy - sin(2*x)*(y + x + a) # = 0
(sin(x)^2 + b)*y*dy + (x + a + 1)*(sin(x)^2 + b)*dy - sin(2*x)*y - x*sin(2*x) + sin(x)^2 - a*sin(2*x) + b # = 0
###
(sin(x)^2 + b)*y*dy + (x + a + 1)*(sin(x)^2 + b)*dy - sin(2*x)*y - x*sin(2*x) - a*sin(2*x) + sin(x)^2 + b # = 0
### Solution:
y = function(x, a, b) {
	# root
	y = lambertWp((sin(x)^2+b) * exp(x+a)) - x - a
	isNA = is.na(y)
	if(any(isNA)) {
		x = x[isNA] # may not play a role ???
		y[isNA] = lambertWn((sin(x)^2+b) * exp(x+a)) - x - a
	}
	y = sapply(y, round0)
	return(y)
}
dy = function(x, a, b) {
	y.x = y(x, a, b)
	x2 = sin(x)^2 + b; x2a = sin(2*x)
	div = x2 * (y.x + x + a + 1)
	dp = x2a*y.x + x2a*(x + a) - x2;
	dp = ifelse(div != 0, dp / div, -1); # may need correction
	return(dp)
}
###
a = -1; b = 1;
curve(y(x, a=a, b=b), from=-2, to=3)
# NO minimum: only inflexion;
sapply(c(-1, (0:3)), line.tan, dx=3, p=y, dp=dy, a=a, b=b)

### y*dy + x*dy + 2*tan(x)*y + 2*x*tan(x) - 2*tan(x) + 1 = 0
a = -1; b = -1;
curve(y(x, a=a, b=b), from=-3, to=2)
# a "virtual" minimum
sapply(c(-1, 1+(0:3)/3), line.tan, dx=3, p=y, dp=dy, a=a, b=b)

###
a = -1; b = -1/4;
curve(y(x, a=a, b=b), from=-3, to=2)
# NO minimum: only inflexion;
sapply(c(-1, (0:3)/2), line.tan, dx=3, p=y, dp=dy, a=a, b=b)


#########################

###################
### Logarythmic ###
###################

### y * log(y) = x + b
# [not run]
y*dy + (x+b)*dy - y # = 0
### Solution:
y = function(x, b) {
	# root
	y = exp(lambertWp(x + b))
	isNA = is.na(y)
	if(any(isNA)) {
		x = x[isNA] # may not play a role ???
		y[isNA] = exp(lambertWn(x + b))
	}
	y = sapply(y, round0)
	return(y)
}
dy = function(x, b) {
	y.x = y(x, b)
	div = y.x + x + b
	dp = y.x;
	dp = ifelse(div != 0, dp / div, -1); # may need correction
	return(dp)
}
###
b = 1;
curve(y(x, b=b), from= - b - 2/3, to = b + 2)
# NO minimum: only inflexion;
sapply(c(-1, (0:3)), line.tan, dx=3, p=y, dp=dy, b=b)

###
b = 2;
curve(y(x, b=b), from= - b - 2/3, to = b + 2)
# NO minimum: only inflexion;
sapply(c(-2, -1, 1, 3), line.tan, dx=3, p=y, dp=dy, b=b)


### y*log(y + a) = x + b
y^2*dy + (x+b)*y*dy + a*(x+b)*dy - y^2 - a*y # = 0
### Solution:
# TODO:
# - solve the equation;
# - check if derivation is correct;


#########################

### Exponentials

### y*e^y + h*e^y = f
# e^y*dy + y*e^y*dy + h*e^y*dy + dh*e^y = df # *y
# y*e^y*dy + y*y*e^y*dy + h*y*e^y*dy + dh*y*e^y = df*y
# (f - h*e^y)*y*dy + (f - h*e^y)*dy + h*(f - h*e^y)*dy + (f - h*e^y)*dh = df*y
# f*y*dy - h*e^y*y*dy + f*dy - h*e^y*dy + h*(f - h*e^y)*dy + (f - h*e^y)*dh = df*y
# f*y*dy - h*(f - h*e^y)*dy + f*dy - h*e^y*dy + h*(f - h*e^y)*dy + (f - h*e^y)*dh = df*y
# f*y*dy - 2*h*e^y*dy + 2*f*dy - df*y = 0

### Examples:

###
# f = x
# h = x - 1
# ...
# TODO: check result;
# TODO: implement snack;


##############################

### e^y * sin(y) = x^2
# z = dy; dz = d2y;
x*dz - x^2*z^3 + 2*x*z^2 - z # = 0
### Solution
y = function(x) {
	# root
	y.f = function(x, v) exp(x)*sin(x) - v^2;
	dy.f = function(x, v) exp(x)*(sin(x) + cos(x))
	y = sapply(x, function(x) newtonRaphson(y.f, 0, dfun=dy.f, v=x)[[1]])
	y = sapply(y, round0)
	return(y)
}
dy = function(x, y.x) {
	if(missing(y.x)) y.x = y(x);
	div = (x^2 + sqrt(exp(2*y.x) - x^4));
	dp = 2*x;
	dp = ifelse(div != 0, dp / div, -1); # may need correction
	return(dp)
}
d2y = function(x) {
	y.x = y(x)
	z = dy(x, y.x=y.x)
	div = x
	dp = x^2*z^3 - 2*x*z^2 + z;
	dp = ifelse(div != 0, dp / div, -1); # TODO: needs correction!
	return(dp)
}
curve(y(x), from= -1.4, to = 1.4)
# global minimum;
sapply(c((-2:2)*2/3), line.tan, dx=3, p=y, dp=dy)
# sigmoidal
curve(dy(x), from= -1.4, to = 1.4, add=T, col="green")
sapply(c((0:5)/5 + 0.01), line.tan, dx=3, p=dy, dp=d2y, col="orange")

# check full sigmoidal
curve(dy(x), from= -1.4, to = 1.4, col="green")
sapply(c((0:5)/5 + 0.01), line.tan, dx=3, p=dy, dp=d2y, col="orange")


#########################
#########################

### Integration by Parts:

# g = y + a*x; a = constant;
f*y*dy + f*(a*x + 1)*dy + (a*f - df)*y + a^2*f*x - a*x*df + a*f # = 0

### f = e^(a*x); df = a*f;
y*dy + (a*x + 1)*dy + a # = 0
# I() =>
1/2*y^2 + (a*x+1)*y - a*I(y) + a*x # = 0
# z = I(y); dz = y;
dz^2 + 2*(a*x+1)*dz - 2*a*z + 2*a*x # = 0

### Solution:
y = function(x, a=1) {
	dz = dy(x, a=a)
	ax = a*x;
	div = 2*a
	dp = dz^2 + 2*(ax+1)*dz + 2*ax;
	dp = if(div != 0) dp / div else  0; # a != 0
	return(dp)
}
dy = function(x, a=1) {
	# root
	y = lambertWp(exp(a*x)) - a*x
	y = sapply(y, round0)
	return(y)
}

### Examples:

### a == 1
curve(y(x), from=-3, to=3)
# a nice global minimum
sapply(c((-2:4)/2), line.tan, dx=3, p=y, dp=dy)

### a == -1
a = -1
curve(y(x, a=a), from=-3, to=3)
# a nice global minimum
sapply(c((-2:4)/2), line.tan, dx=3, p=y, dp=dy, a=a)

### a == 1/2
a = 1/2
curve(y(x, a=a), from=-3, to=3)
# a nice global minimum
sapply(c((-2:4)/2), line.tan, dx=3, p=y, dp=dy, a=a)


### Higher Powers:
# (y + a*x)*e^y = 1;
y*dy + (a*x + 1)*dy + a # = 0
# I() =>
# z = I(y); dz = y;
dz^2 + 2*(a*x+1)*dz - 2*a*z + 2*a*x # = 0
# Original * y =>
y^2*dy + (a*x + 1)*y*dy + a*y # = 0
# I() =>
1/3*y^3 + 1/2*(a*x + 1)*y^2 - a/2*I(y^2) + a*I(y) # = 0
# a*I(y) = a*z = 1/2 * (y^2 + 2*(a*x+1)*y + 2*a*x);
1/3*y^3 + 1/2*(a*x + 1)*y^2 - a/2*I(y^2) + 1/2 * (y^2 + 2*(a*x+1)*y + 2*a*x) # = 0
1/3*y^3 + 1/2*a*x*y^2 + y^2 + (a*x+1)*y - a/2*I(y^2) + a*x # = 0

### Solution:
# y = actual dy from above;
y2.f = function(x, a=1) {
	yx = dy(x, a=a)
	ax = a*x
	Iy2 = 1/3*yx^3 + 1/2*ax*yx^2 + yx^2 + (ax+1)*yx + ax;
	Iy2 = Iy2 * 2/a;
	return(Iy2)
}
y2.I = function(x, a=1, lower=0) {
	Iy2 = sapply(x, function(x) integrate(function(x, a) dy(x, a)^2, lower=lower, upper=x, a=a)$value)
	return(Iy2)
}
dy2 = function(x, a) dy(x, a=a)^2
###
a = 1
curve(y2.f(x, a=a), from=-2, to=3)
# + 0.1 to separate curves;
curve(y2.I(x, a=a, lower=-1.473744 + 0.1), add=T, col="green")
sapply(c((-3:2)/2, 2), line.tan, dx=3, p=y2.f, dp=dy2, a=a)
# TODO: compute exact value for -1.473744;

# separately: only the integration: I(y^2)
curve(y2.I(x, a=a), from=-2, to=3)
sapply(c((-3:2)/2, 2), line.tan, dx=3, p=y2.I, dp=dy2, a=a)

###
a = 3/2
curve(y2.f(x, a=a), from=-2, to=3)
# + 0.1 to separate curves;
curve(y2.I(x, a=a, lower=-0.98 + 0), add=T, col="green")
sapply(c((-3:2)/2, 2), line.tan, dx=3, p=y2.f, dp=dy2, a=a)



### Example 2:
### f = e^(-a*x^2); df = -2*a*x*f;
# complete f = e^(-a*x^2 + a*x);
y*dy + (a*x + 1)*dy + 2*a*x*y + 2*a^2*x^2 + a # = 0

### Solution:
y = function(x, a=1) {
	# root
	y = lambertWp(exp(-a*x^2 + a*x)) - a*x
	y = sapply(y, round0)
	return(y)
}
dy = function(x, a=1) {
	yx = y(x, a=a)
	ax = a*x;
	div = -(yx + ax + 1)
	dp = 2*ax*yx + 2*ax^2 + a;
	dp = ifelse(div != 0, dp / div,  -1/2); # TODO: check!
	return(dp)
}
### a == 1
lim = 3
curve(y(x), from=-lim, to=lim)
# quasi-bi-sigmoidal
sapply(c((-3:3)/1.5, -0.50404905), line.tan, dx=3, p=y, dp=dy)

### a == 1/2
a = 1/2
curve(y(x, a=a), from=-3, to=3)
# quasi-bi-sigmoidal
sapply(c((-3:3)/1.5), line.tan, dx=3, p=y, dp=dy, a=a)

### a == 1/3
a = 1/3
curve(y(x, a=a), from=-3, to=3)
# quasi-bi-sigmoidal
sapply(c((-3:3)/1.5), line.tan, dx=3, p=y, dp=dy, a=a)

### a == -1/3
a = -1/3
curve(y(x, a=a), from=-3, to=3)
# global minimum
sapply(c((-3:3)/1.5), line.tan, dx=3, p=y, dp=dy, a=a)

### TODO: Integration by parts;


###########################
###########################

### e^(y + a1*x) * e^(e^(y + a1*x)) = x^n

### D =>
(dy + a1)*(1 + e^(y + a1*x)) - n/x # = 0
### D2 =>
n*x*d2y - x^2*(dy + a1)^3 + n*x*(dy + a1)^2 + n*(dy + a1) # = 0;
# w = dy:
# n*x*dw - x^2*(w + a1)^3 + n*x*(w + a1)^2 + n*(w + a1) = 0;
### Solution & Plot
y = function(x, a=1, n=2) {
	# root
	y = log(lambertWp(x^n)) - a[1]*x
	y = sapply(y, round0)
	return(y)
}
dy = function(x, a=1, n=2) {
	yx = y(x, a=a, n=n)
	fx = x*(1 + exp(yx + a[1]*x));
	div = fx;
	dp = n - a[1]*fx;
	dp = ifelse(div != 0, dp / div,  1E+3); # TODO: check!
	return(dp)
}
d2y = function(x, a=1, n=2) {
	yx = y(x, a=a, n=n)
	dyx = dy(x, a=a, n=n)
	f = (dyx + a[1]);
	div = n*x;
	dp = x^2*f^3 - n*x*f^2 - n*f;
	dp = ifelse(div != 0, dp / div,  0); # TODO: check!
	return(dp)
}
### Plot:
a = 1; n = 2;
curve(y(x, a=a, n=n), from= -2, to = 2)
# global "minimum" / horn;
sapply(c((-3:3)*2/5), line.tan, dx=3, p=y, dp=dy, a=a, n=n)
#
curve(dy(x, a=a, n=n), add=T, col="green")
sapply(c(-(3:1)/5, 1:5/5), line.tan, dx=3, p=dy, dp=d2y, a=a, n=n, col="orange")


### Example 2:
a = 1/3; n = 2;
curve(y(x, a=a, n=n), from= -2, to = 2, ylim=c(-7, 3))
# global "minimum" / horn;
sapply(c((-3:3)*2/5), line.tan, dx=3, p=y, dp=dy, a=a, n=n)
#
curve(dy(x, a=a, n=n), add=T, col="green")
sapply(c(-(3:1)/5, 1:5/5), line.tan, dx=3, p=dy, dp=d2y, a=a, n=n, col="orange")


##########################
##########################

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
### Solve Liniar =>
# ...

### Special Cases:

### Order: n = 1
f*g*dy = g*df * log(g) + f*dg * log(f)
f*g*d2y + df*g*dy + f*dg*dy - df - dg =
	(g*d2f + df*dg) * log(g) + (f*d2g + df*dg) * log(f)

### Solve Liniar =>
log(g) = ...
log(f) = ...

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
curve(y(x, k=k, n=n), from= 1.01, to = 3, ylim=c(-3, 3))
# global "minimum" / horn;
sapply(c(3/5 + (1:3)*3/5), line.tan, dx=3, p=y, dp=dy, k=k, n=n)
#
curve(dy(x, k=k, n=n), add=T, col="green")
sapply(c(3/5 + (1:3)*3/5), line.tan, dx=3, p=dy, dp=d2y, k=k, n=n, col="orange")


### Ex 2:
k = 3; n = 1;
curve(y(x, k=k, n=n), from= 3.01, to = 6, ylim=c(-3, 3))
# global "minimum" / horn;
sapply(c(3 + (1:3)*3/5), line.tan, dx=3, p=y, dp=dy, k=k, n=n)
#
curve(dy(x, k=k, n=n), add=T, col="green")
sapply(c(3 + (1:3)*3/5), line.tan, dx=3, p=dy, dp=d2y, k=k, n=n, col="orange")


