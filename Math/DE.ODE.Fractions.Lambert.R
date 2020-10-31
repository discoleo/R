
########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Fractions: Lambert
###
### draft v.0.1a


### History

### Order 1 Non-Liniar
###
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
	dp = if(div != 0) dp / div else -1;
	return(dp)
}
curve(y(x), from=-1/5, to=3)
# a nice global minimum
sapply(c((0:4)/2), line.tan, dx=3, p=y, dp=dy)

#####################
### (x + y)*e^y = x^2 + b
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
	div = x2*y.x + x2*(x+1)
	dp = 2*x*y.x + x2 - 2*b
	# TODO: find BUG: must be -1 (but why?);
	dp = if(div != 0) dp / div else -1;
	return(dp)
}
curve(y(x, b=1/2), from=-1, to=3)
# a nice global minimum
sapply(c((0:4)/2), line.tan, dx=3, p=y, dp=dy, b=1/2)


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
	dp = if(div != 0) dp / div else 0;
	return(dp)
}
curve(y(x, n=3, b=1), from=-1/5, to=3)
# a nice global minimum
sapply(c((0:4)/2), line.tan, dx=3, p=y, dp=dy, n=3, b=1)


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
