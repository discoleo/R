
########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Fractions: Lambert
###
### draft v.0.2b


### History

### Order 1 Non-Liniar
###
### draft v.0.2b: [12-11-2020]
# - cleanup: moved Gaussian type to
#   new file: DE.ODE.Gaussian.R;
### draft v.0.2a - v.0.2a-xTerm: [11-11-2020]
# - solved:
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

### Examples

### (x + y)*e^y = x^2
x*y*dy + x*(x+1)*dy - 2*y - x # = 0
y*dy + (x+1)*dy - 2/x * y - 1 # = 0
### Solution:
y = function(x) {
	# root
	y = lambertWp(x^2 * exp(x)) - x
	y = sapply(y, round0)
	return(y)
}
dy = function(x) {
	y.x = y(x)
	div = x*y.x + x*(x+1)
	dp = 2*y.x + x
	# TODO: find BUG: must be -1 (but why?);
	dp = ifelse(div != 0, dp / div, -1);
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

