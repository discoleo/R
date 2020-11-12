########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Gaussian
###
### draft v.0.1a


### History

### Liniar / Non-Liniar Gaussian-type
###
### draft v.0.1a: [12-11-2020]
# - moved section on Gaussian-types
#   to this new file;

### [old file]
### Exponential/Lambert: DE.ODE.Fractions.Lambert.R
### draft v.0.2a - v.0.2a-xTerm: [11-11-2020]
# - solved:
#   d2z + 2*x*dz - 2*z = 0;
#   d2z + 2*x*dz - 4*z = 0; [v.0.2a-x2Lev3]
#   d2z + 2*x*dz - 6*z = 0; [v.0.2a-x2Lev4]
#   d2z + 2*x*dz - 2*z = b0*x^2; [v.0.2a-xTerm]
#   d3z + 3*x^2*d2z - 6*x*dz + 6*z = 0; [v.0.2a-pow3]


#########################

### Helper functions

library(pracma)
# needed for Lambert W;


# include: DE.ODE.Helper.R;
source("DE.ODE.Helper.R")

#########################
#########################
#########################

##############
### Theory ###
##############

### by Integration

### Base:
# y = e^f(x) + b(x);
### D() =>
# dy = df * y + db - df*b;
### Integration by parts =>
# y = df * I(y) + I(d2f * I(y)) + I(df*b);
### z = I(y):
# dz = df*z + I(d2f * z) + I(df*b);


####################
####################

### Examples:

### y = e^(-x^2)
# [not run]
# dy = -2*x*y;
d2z + 2*x*dz - 2*z # = 0
# where dz = sqrt(pi)/2 * erf(x); d2z = e^(-x^2);
### Solution:
y = function(x, b0=0, lower = if(b0 == 0) -Inf else 0) {
	dy.x = dy(x, b0=b0, lower=lower)
	d2y.x = d2y(x, b0=b0)
	y = 1/2 * d2y.x + x*dy.x - b0/2 *x^2;
	y = sapply(y, round0)
	return(y)
}
dy = function(x, b0=0, lower = if(b0 == 0) -Inf else 0) {
	if(b0 == 0) {
		dp = pnorm(x * sqrt(2)) * sqrt(pi)
	} else {
		dp = sapply(x, function(x) integrate(function(x) exp(-x^2) + b0, lower=lower, upper=x)$value)
	}
	return(dp)
}
d2y = function(x, b0=0) {
	dp = exp(-x^2) + b0;
	return(dp)
}
# b0 == 0;
curve(y(x), from= -3, to = 3)
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=y, dp=dy)
# sigmoidal
curve(dy(x), add=T, col="green")
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=dy, dp=d2y, col="orange")


### y = e^(-x^2) + b0
# [not run]
# dy = -2*x*y + 2*b0*x;
# d2z = y;
d2z + 2*x*dz - 2*z - b0*x^2 # = 0
### Plot:
# using functions defined above;
b0 = -1;
curve(y(x, b0=b0), from= -3, to = 3, ylim=c(-3, 2))
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=y, dp=dy, b0=b0)
# sigmoidal
curve(dy(x, b0=b0), add=T, col="green")
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=dy, dp=d2y, b0=b0, col="orange")


#####################
### Higher Levels ###
#####################

### y = e^(-x^2)
# d3z = y;
# [not run]
# Level 3:
d2z + 2*x*dz - 4*z # = 0
# Level 4:
d2z + 2*x*dz - 6*z # = 0
### Solution:
dny = function(x, n=2) {
	dp = exp(-x^n)
	return(dp)
}
y = function(x, n=2, level=3) {
	dz = if(level == 3) dp2y(x, n=n) else dp3y(x, n=n);
	d2z = if(level == 3) dp1y(x, n=n) else dp2y(x, n=n);
	z = (1/2 * d2z + x*dz) / (level - 1)
	return(z);
}
dp3y = function(x, n=2, ...) {
	# could also use (recursively) the formula in this+previous sections;
	integrate.y(x, FUN=dp2y, n=n, lower=-Inf)
}
dp2y = function(x, n=2) {
	# could also use (recursively) the formula in previous section;
	integrate.y(x, FUN=dp1y, n=n, lower=-Inf)
}
dp1y = function(x, n=2) {
	integrate.y(x, FUN=dny, n=n, lower=-Inf)
}
integrate.y = function(x, FUN=dny, n=2, lower=-Inf) {
	dp = sapply(x, function(x) integrate(FUN, lower=lower, upper=x, n=n)$value)
	return(dp)
}
### Levels: +3, +2, +1;
curve(y(x), from= -3, to = 2.2)
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=y, dp=dp2y)
# *NON*-sigmoidal / seems increasing
curve(dp2y(x), add=T, col="green")
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=dp2y, dp=dp1y, col="orange")


### Levels: +4, +3, +2;
# Note: takes very long !!!
# - should be implemented using the compact/derived formulas!
curve(y(x, level=4), from= -3, to = 2.2)
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=y, dp=dp3y, level=4)
# *NON*-sigmoidal / seems increasing
curve(dp3y(x), add=T, col="green")
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=dp3y, dp=dp2y, col="orange")


#####################
#####################

#####################
### Higher Powers ###
#####################

### dy = e^(-x^3)
# [not run]
# d2y = -3*x^2*dy;
d3z + 3*x^2*d2z - 6*x*dz + 6*z # = 0
# where d3z = e^(-x^3);
### Solution:
y = function(x) {
	dz = dy(x)
	d2z = d2y(x)
	d3z = d3y(x)
	y = - 1/6 * (d3z + 3*x^2*d2z - 6*x*dz);
	y = sapply(y, round0)
	return(y)
}
dy = function(x) {
	dp = sapply(x, function(x) integrate(d2y, lower=0, upper=x)$value)
	return(dp)
}
d2y = function(x) {
	dp = sapply(x, function(x) integrate(d3y, lower=0, upper=x)$value)
	return(dp)
}
d3y = function(x) {
	dp = exp(-x^3)
	return(dp)
}
curve(y(x), from= -0.01, to = 3)
#
sapply(c(0:3 * 3/4), line.tan, dx=3, p=y, dp=dy)
# sigmoidal
curve(dy(x), add=T, col="green")
sapply(c(0:3 * 3/4), line.tan, dx=3, p=dy, dp=d2y, col="orange")


