
########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Fractions: Lambert
###
### draft v.0.1e-fix


### History

### Order 1 Non-Liniar
###
### draft v.0.1e - v.0.1e-fix: [05-11-2020]
# - started to explore various y powers:
#   x*y^(1/2)*dy + x*dy - 2*y = 0;
# - TODO: find bug; [FIXED]
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
