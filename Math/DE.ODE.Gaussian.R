########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Gaussian
###
### draft v.0.3g-test

#############
### Types ###
#############

### Simple:
# d2z + 2*x*dz - 2*n*z = 0;
### e^f(x)-type, f = order 2:
# Level 1: d2z + (2*k2*x + k1)*dz - 2*k2*z = 0;
# Level n: TODO;
### Higher powers:
# Level 1: d3z + 3*x^2*d2z - 6*x*dz + 6*z = 0
# Level n: TODO;
### Others:
# - and many more: non-homogenous variants;


###############
### History ###
###############

### Liniar / Non-Liniar Gaussian-type
###
### draft v.0.3g - v.0.3g-test:
# - derived from:
#   y = x^m * I(e^(x^n)) * I(e^(-x^n));
### draft v.0.3d - v.0.3f:
# - derived from: y = exp(I),
#   where I = I(exp(P(x))) dx;
#   y*d2y - dy^2 - dp*y*dy = 0;
# - derived from: y = exp(I^2),
#   where I = I(exp(P(x))) dx; [v.0.3e]
#   y*d2y - dy^2 - dp*y*dy - 2*e^(2*p) * y^2 = 0;
# - derived from: y^2 + b*y = e^I; [v.0.3f]
#   (y^2 + b*y)*(2*y + b)*d2y - ... = 0;
### draft v.0.3c:
# - derived from:
#   y = sin(x^n) * exp(-x^n);
### draft v.0.3b:
# - more variants & examples;
### draft v.0.3a:
# - y = I(e^f(x)) * I(e^(-f(x))), e.g. f(x) = x^n:
#   (d2y + n*x^(n-1)*dy - 2)*(d2y - n*x^(n-1)*dy - 2) + 4*n^2*x^(2*n-2)*y = 0;
### draft v.0.2d - v.0.2e:
# - derived from Combinations of Exponentials:
#   (3*x^2 - 2*x)*d2y + (9*x^4 - 4*x^2 - 6*x + 2)*dy + 6*x^2*(3*x^3 - 2*x^2 - 1)*y = 0;
# - slight generalization; [v.0.2e]
### draft v.0.2c:
# - derived from Fractions:
#   x*d2y - k*x*dy + 2*dy - k*y = 0;
#   x*d2y + 2*dy - k^2*x*y = 0;
### draft v.0.2b - v.0.2b-gen:
# - sinh-type variants, including generalization:
#   x*d2y + k*x^2*sin(y)*dy - dy = 0;
#   x*d2y + k*x^n*sin(y)*dy - (n-1)*dy = 0; [v.0.2b-gen]
### draft v.0.2a: [12-11-2020]
# - solved: d2z + (4*x + 1)*dz - 4*z = 0;
#   including generalisations of type:
#   d2z + (2*k2*x + k1)*dz - 2*k2*z = 0;
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

### Section A: by Integration

### Base:
# y = e^f(x) + b(x);
### D() =>
# dy = df * y + db - df*b;
### Integration by parts =>
# y = df * I(y) - I(d2f * I(y)) + b - I(df*b);
### z = I(y):
# dz = df*z - I(d2f * z) + b - I(df*b);


### Section B: Mixt

### B.1. Mixt Liniar
### y = e^(F1(x)) + e^(F2(x))

### B.2. Mixt Exponential-Trigonometric
### y = sin(F1(x)) * exp(-F1(x))

### B.3. Mixt Non-Liniar
### y = Integral(exp(F1(x))) * Integral(exp(-F1(x)))

### B.4. Non-Liniar Double Exp
### y = e^I, where I = I(e^P(x)) dx;


####################
####################

### Examples:

### y = e^(-x^2)
# [not run]
# dy = -2*x*y;
# d2y + 2*x*dy + 2*y = 0
### I(D1) by parts =>
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
curve(y(x), from= -3, to = 3, ylim=c(-1/3, 4))
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


##################
### y = e^f(x) ###

### y = e^-(k2*x^2 + k1*x)
# dy = - (2*k2*x + k1)*y
### Integration by parts =>
# y = -(2*k2*x + k1)*I(y) + 2*k2*I(I(y));
# [not run]
d2z + (2*k2*x + k1)*dz - 2*k2*z # = 0

### Solution & Plot:
y = function(x, k, lower = -Inf) {
	dz = dy(x, k=k, lower=lower)
	d2z = d2y(x, k=k)
	# Note: rev(k)
	y = 1/2 * (d2z + (2*k[1]*x + k[2])*dz) / k[1];
	y = sapply(y, round0)
	return(y)
}
dy = function(x, k, lower = -Inf) {
	dp = sapply(x, function(x) integrate(d2y, lower=lower, upper=x, k=k)$value)
	return(dp)
}
d2y = function(x, k) {
	if(length(k) == 1) k = c(k, 0)
	len = length(k)
	dp = sapply(x, function(x) exp(-sum(x^rev(seq(len)) * k)));
	return(dp)
}
# d2z + (4*x + 1)*dz - 4*z = 0
k = c(2, 1)
curve(y(x, k), from= -3, to = 3)
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=y, dp=dy, k=k)
# sigmoidal
curve(dy(x, k), add=T, col="green")
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=dy, dp=d2y, k=k, col="orange")



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

### y = I(e^(-x^n))

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
### Plot:
curve(y(x), from= -0.01, to = 3)
sapply(c(0:3 * 3/4), line.tan, dx=3, p=y, dp=dy)
# sigmoidal
curve(dy(x), add=T, col="green")
sapply(c(0:3 * 3/4), line.tan, dx=3, p=dy, dp=d2y, col="orange")


#######################
#######################

#######################
### Section B: Mixt ###
#######################

################
### Section B.1:

### Liniar Combinations

### y = a2*e^(-x^2) + a3*e^(-x^3)
# [not run]
dy  = -2*a2*x*exp(-x^2) - 3*a3*x^2*exp(-x^3)
d2y = 4*a2*x^2*exp(-x^2) - 2*a2*exp(-x^2) + 9*a3*x^4*exp(-x^3) - 6*a3*x*exp(-x^3)
# e^(-x^2) = (dy + 3*x^2*y) / (a2*(3*x^2 - 2*x))
# e^(-x^3) = - (dy + 2*x*y) / (a3*(3*x^2 - 2*x))
# =>

### ODE:
(3*x^2 - 2*x)*d2y + (9*x^4 - 4*x^2 - 6*x + 2)*dy + 6*x^2*(3*x^3 - 2*x^2 - 1)*y # = 0


### Plot:
y = function(x, a=c(1, 1)) {
	y = a[1] * exp(-x^2) + a[2] * exp(-x^3)
	y = round0(y)
	return(y)
}
dy = function(x, a=c(1, 1)) {
	y.x = y(x, a=a)
	dp = - (2*a[1]*x*exp(-x^2) + 3*a[2]*x^2*exp(-x^3))
	return(dp)
}
d2y = function(x, a=c(1, 1)) {
	y.x = y(x, a=a)
	dy.x = dy(x, a=a)
	dp = (9*x^4 - 4*x^2 - 6*x + 2)*dy.x + 6*x^2*(3*x^3 - 2*x^2 - 1)*y.x;
	dp = - dp;
	div = 3*x^2 - 2*x;
	dp = ifelse(div != 0, dp / div, -Inf); # TODO
	return(dp)
}
### Plot:
a = c(1, 1) # has NO effect on eq of D2;
curve(y(x, a=a), from= -1, to = 3, ylim=c(-2, 3))
sapply(c(-3:2 * 2/5, 3/2, 2), line.tan, dx=3, p=y, dp=dy, a=a)
# non-sigmoidal
curve(dy(x, a=a), add=T, col="green")
sapply(c(-3:2 * 2/5, 3/2, 2), line.tan, dx=3, p=dy, dp=d2y, a=a, col="orange")

###
a = c(1, -1) # a = c(2, -1/2)
# although has NO effect on eq. of D2, d2y depends indirectly;
curve(y(x, a=a), from= -1, to = 3, a=a, ylim=c(-2, 3))
sapply(c(-3:2 * 2/5, 3/2, 2), line.tan, dx=3, p=y, dp=dy, a=a)
# non-sigmoidal
curve(dy(x, a=a), add=T, col="green")
sapply(c(0.05 + -3:3 /5, 3/2, 2), line.tan, dx=3, p=dy, dp=d2y, a=a, col="orange")


######################
### Generalization ###
### (part)         ###

### y = a1*e^(-(x^n + b1*x)) + a2*e^(-(x^n + b2*x))
# [not run]
dy = - a1*(n*x^(n-1) + b1)*exp(-(x^n + b1*x)) +
     - a2*(n*x^(n-1) + b2)*exp(-(x^n + b2*x))
d2y = a1*(n*x^(n-1) + b1)^2 * exp(-(x^n + b1*x)) +
      a2*(n*x^(n-1) + b2)^2 * exp(-(x^n + b2*x)) +
	- a1*n*(n-1)*x^(n-2) * exp(-(x^n + b1*x)) +
	- a2*n*(n-1)*x^(n-2) * exp(-(x^n + b2*x));
# e^(-(x^n + b1*x)) =   (dy + (n*x^(n-1) + b2)*y) / (a1*(b2 - b1))
# e^(-(x^n + b2*x)) = - (dy + (n*x^(n-1) + b1)*y) / (a2*(b2 - b1))
# =>
(b2 - b1)*d2y - (n*x^(n-1) + b1)^2 * (dy + (n*x^(n-1) + b2)*y) +
  + (n*x^(n-1) + b2)^2 * (dy + (n*x^(n-1) + b1)*y) +
  + n*(n-1)*x^(n-2) * (dy + (n*x^(n-1) + b2)*y) +
  - n*(n-1)*x^(n-2) * (dy + (n*x^(n-1) + b1)*y) # = 0
(b2 - b1)*d2y - (n*x^(n-1) + b1)^2 * (dy + (n*x^(n-1) + b2)*y) +
  + (n*x^(n-1) + b2)^2 * (dy + (n*x^(n-1) + b1)*y) +
  - n*(n-1)*(b1 - b2)*x^(n-2)*y # = 0
### ODE:
d2y + (2*n*x^(n-1) + b1 + b2)*dy +
  + (n^2*x^(2*n-2) + n*(b1 + b2)*x^(n-1) + n*(n-1)*x^(n-2) + b1*b2)*y # = 0

### Examples:
### n = 2
d2y + (4*x + b1 + b2)*dy + (4*x^2 + 2*(b1 + b2)*x + b1*b2 + 2)*y # = 0
### b1 + b2 = 0
d2y + 2*n*x^(n-1) * dy + (n^2*x^(2*n-2) + n*(n-1)*x^(n-2) - b1^2)*y # = 0
### n = 2; b1 = - b2 = sqrt(2);
d2y + 4*x*dy + 4*x^2*y # = 0


### Solution:
y = function(x, a=c(1, 1), b=c(-1, 1), n=2) {
	y = sapply(x, function(x) sum(a * exp(-x^n - b*x)))
	y = round0(y)
	return(y)
}
dy = function(x, a=c(1, 1), b=c(-1, 1), n=2) {
	# y.x = y(x, a=a, b=b, n=n)
	nxn = n*x^(n-1);
	dp = sapply(seq_along(x), function(id) {
		xx = x[id]
		x.exp = -a * exp(-xx^n - b*xx);
		sum((nxn[id] + b) * x.exp)
	})
	return(dp)
}
d2y = function(x, a=c(1, 1), b=c(-1, 1), n=2) {
	y.x = y(x, a=a, b=b, n=n)
	dy.x = dy(x, a=a, b=b, n=n)
	nxn1 = n*x^(n-2);
	nxn = x * nxn1; nxsq = nxn^2;
	dp = (2*nxn + b[1] + b[2])*dy.x +
		+ (nxsq + (b[1] + b[2])*nxn + (n-1)*nxn1 + b[1]*b[2])*y.x
	div = -1;
	dp = ifelse(div != 0, dp / div, -Inf); # TODO
	return(dp)
}
### Plot:
n = 2; a = c(1, 1) # a[] has NO effect on eq of D2;
curve(y(x, a=a, n=n), from= -2, to = 3, ylim=c(-2, 3))
sapply(c(-3:2 * 4/7, 2), line.tan, dx=3, p=y, dp=dy, a=a, n=n)
# wavelet
curve(dy(x, a=a, n=n), add=T, col="green")
sapply(c(-3:2 * 4/7, 3/2, 2), line.tan, dx=3, p=dy, dp=d2y, a=a, n=n, col="orange")

### n = 3
n = 3; a = c(1, 1) # a[] has NO effect on eq of D2;
curve(y(x, a=a, n=n), from= -1, to = 3, ylim=c(-4, 4))
sapply(c(-3:0 * 1/7, 0.6, 1, 3/2), line.tan, dx=3, p=y, dp=dy, a=a, n=n)
# wavelet
curve(dy(x, a=a, n=n), add=T, col="green")
sapply(c(-3:3 * 1/7, 1, 3/2, 2), line.tan, dx=3, p=dy, dp=d2y, a=a, n=n, col="orange")
### TODO: check d2y(1)!


#################

### Variant:
### y = a1*e^(x^n) + a2*e^(-x^n)
# [not run]
dy = a1*n*x^(n-1) * exp(x^n) - a2*n*x^(n-1) * exp(-x^n)
### D2 =>
d2y = a1*n^2*x^(2*n-2) * exp(x^n) + a1*n*(n-1)*x^(n-2) * exp(x^n) +
	+ a2*n^2*x^(2*n-2) * exp(-x^n) - a2*n*(n-1)*x^(n-2) * exp(-x^n)
### Linear system:
exp(x^n)  = ( dy + n*x^(n-1)*y) / (2*a1*n*x^(n-1))
exp(-x^n) = (-dy + n*x^(n-1)*y) / (2*a2*n*x^(n-1))
# =>
2*x*d2y = (n*x^n + (n-1))*(dy + n*x^(n-1)*y) - (n*x^n - (n-1))*(dy - n*x^(n-1)*y)

### ODE:
x*d2y - (n-1)*dy - n^2*x^(2*n-1)*y # = 0

### Solution:
y = function(x, n=2, a=c(1, 1), b=c(1, -1)) {
	y = sapply(x, function(x) sum(a * exp(b*x^n)))
	y = round0(y)
	return(y)
}
dy = function(x, n=2, a=c(1, 1), b=c(1, -1)) {
	# y.x = y(x, a=a, b=b, n=n)
	a[2] = -a[2];
	nxn = n*x^(n-1);
	dp = sapply(seq_along(x), function(id) {
		xx = x[id]
		sum(a * nxn[id] * exp(b*xx^n));
	})
	return(dp)
}
d2y = function(x, n=2, a=c(1, 1), b=c(1, -1)) {
	y.x = y(x, n=n, a=a, b=b)
	dy.x = dy(x, n=n, a=a, b=b)
	dp = (n-1)*dy.x + n^2*x^(2*n-1)*y.x
	div = x;
	dp = ifelse(div != 0, dp / div, 0); # TODO
	return(dp)
}
### Plot:
n = 2; a = c(1, 1) # a[] has NO effect on eq of D2;
curve(y(x, a=a, n=n), from= -2, to = 2, ylim=c(-2, 5))
sapply(c(-3:3 * 4/7, 2), line.tan, dx=3, p=y, dp=dy, a=a, n=n)
#
curve(dy(x, a=a, n=n), add=T, col="green")
sapply(c(-3:3 * 2/7), line.tan, dx=3, p=dy, dp=d2y, a=a, n=n, col="orange")


###########################
###########################

### Type: sinh(x) / cosh(x)

### tan(y) = 1/2 * (e^(k/2 * x^2) - e^(-k/2 * x^2))
# [not run]
dy = k*x*cos(y)
x*d2y + k*x^2*sin(y)*dy - dy # = 0
# alternative: x*d2y - dy + 1/2 * k^2*x^3*sin(2*y) # = 0
### Solution:
y = function(x, k=2, n=2) {
	val = k*x^n / n;
	y = 1/2 * (exp(val) - exp(-val))
	y = atan(y)
	y = sapply(y, round0)
	return(y)
}
dy = function(x, k=2, n=2) {
	y.x = y(x, k=k, n=n)
	xn = if(n == 2) x else x^(n-1);
	dp = k * xn * cos(y.x)
	return(dp)
}
d2y = function(x, k=2, n=2) {
	y.x = y(x, k=k, n=n)
	dy.x = dy(x, k=k, n=n)
	dp = (n-1)*dy.x - k*x^n * sin(y.x)*dy.x;
	div = x;
	lim = if(n > 2) 0 else if(n == 2) (n-1)*k else 1E+3;
	dp = ifelse(div != 0, dp / div, lim); # probably correct
	return(dp)
}
### k/2 == 1
curve(y(x), from= -3, to = 3, ylim=c(-2, 2))
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=y, dp=dy)
# quasi-sigmoidal
curve(dy(x), add=T, col="green")
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=dy, dp=d2y, col="orange")


### k/2 == 1/2
k = 1
curve(y(x, k=k), from= -3, to = 3, ylim=c(-1, 2))
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=y, dp=dy, k=k)
# quasi-sigmoidal
curve(dy(x, k=k), add=T, col="green")
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=dy, dp=d2y, k=k, col="orange")


### tan(y) = 1/2 * (e^(k/n*x^n) - e^(-k/n*x^n))
# [not run]
dy = k*x^(n-1)*cos(y)
x*d2y + k*x^n*sin(y)*dy - (n-1)*dy # = 0
x*d2y - (n-1)*dy + 1/2 * k^2*x^(2*n-1)*sin(2*y) # = 0

### Plot:
### n = 3; k/n == 1/3
n = 3; k = 1;
curve(y(x, k=k, n=n), from= -3, to = 3, ylim=c(-2, 2))
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=y, dp=dy, k=k, n=n)
# non-sigmoidal
curve(dy(x, k=k, n=n), add=T, col="green")
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=dy, dp=d2y, k=k, n=n, col="orange")

### n = 3; k/n == 1
n = 3; k = 3;
curve(y(x, k=k, n=n), from= -3, to = 3, ylim=c(-2, 2))
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=y, dp=dy, k=k, n=n)
# non-sigmoidal
curve(dy(x, k=k, n=n), add=T, col="green")
sapply(c(-3:3 * 3/4), line.tan, dx=3, p=dy, dp=d2y, k=k, n=n, col="orange")


### Integration by parts
# n = 2
x*d2y + k*x^2*sin(y)*dy - dy # = 0
### I() =>
x*dy - k*x^2*cos(y) + 2*k*I(x*cos(y)) - 2*y # = 0
# x*dy = k*x^2*cos(y)
k*x*I(cos(y)) - k*I(I(cos(y))) - y # = 0

### I() =>
x*dy - k*x^2*cos(y) + 2*k*x*I(cos(y)) - 2*k*I(I(cos(y))) - 2*y # = 0


##################
### Extensions ###

### tan(y) = 1/2 * (e^(k/n*x^n + y) - e^(-k/n*x^n - y))
# [not run]
dy = (k*x^(n-1) + dy)*cos(y)
# =>
sin(y)^2*dy^2 = k*x^(n-1)*(k*x^(n-1) + 2*dy)
#
x*d2y - x*d2y*cos(y) + x*sin(y)*dy^2 + k*x^n*sin(y)*dy - (n-1)*dy # = 0
### n = 2
dy = (k*x + dy)*cos(y) # * y



###########################
###########################

### Simple Fractions

### y = e^(k*x) / x
# [not run]
x*dy = k*x*y - y
### D2 =>
x*d2y - k*x*dy + 2*dy - k*y # = 0
### Variant 1:
x*d2y + 2*dy - k^2*x*y # = 0
### Variant 2:
x^2*d2y + 2*x*dy - k^2*x^2*y # = 0
x^2*d2y - k^2*x^2*y + 2*k*x*y - 2*y # = 0


### Solution:
y = function(x, k=1, n=1) {
	xn = if(n == 1) x else x^n;
	y = exp(k*x) / xn
	# y = sapply(y, round0)
	y = round0(y)
	return(y)
}
dy = function(x, k=1, n=1, variant=0) {
	y.x = y(x, k=k, n=n)
	dp = (k*x - 1)*y.x;
	div = x;
	dp = ifelse(div != 0, dp / div, -Inf); # probably correct
	return(dp)
}
d2y = function(x, k=1, n=1, variant=0) {
	# variant = test equations for variants;
	y.x = y(x, k=k, n=n)
	dy.x = dy(x, k=k, n=n)
	dp = if(variant == 0) { (k*x - 2)*dy.x + k*y.x; }
	else if(variant == 1) { -(2*dy.x - k^2*x*y.x); }
	else if(variant == 2) { (k^2*x^2 - 2*k*x + 2)*y.x; }
	div = if(variant == 2) x^2 else x;
	dp = ifelse(div != 0, dp / div, -Inf); # probably correct (check sign)
	return(dp)
}
### Plot:
k = 1;
curve(y(x, k=k), from= -3, to = 3)
sapply(c(-3:3 * 2/7), line.tan, dx=3, p=y, dp=dy, k=k)
# non-sigmoidal
curve(dy(x, k=k), add=T, col="green")
sapply(c(-3:3 * 2/7), line.tan, dx=3, p=dy, dp=d2y, k=k, variant=2, col="orange")

### k = 1/2
k = 1/2;
curve(y(x, k=k), from= -3, to = 3)
sapply(c(-3:3 * 2/7), line.tan, dx=3, p=y, dp=dy, k=k)
# non-sigmoidal
curve(dy(x, k=k), add=T, col="green")
sapply(c(-3:3 * 2/7), line.tan, dx=3, p=dy, dp=d2y, k=k, variant=1, col="orange")

### k = -1
k = -1;
curve(y(x, k=k), from= -3, to = 3)
sapply(c(-3:3 * 2/7), line.tan, dx=3, p=y, dp=dy, k=k)
# non-sigmoidal
curve(dy(x, k=k), add=T, col="green")
sapply(c(-3:3 * 2/7), line.tan, dx=3, p=dy, dp=d2y, k=k, col="orange")

### k = -3/2
k = -3/2;
curve(y(x, k=k), from= -3, to = 3)
sapply(c(-3:3 * 2/7), line.tan, dx=3, p=y, dp=dy, k=k)
# non-sigmoidal
curve(dy(x, k=k), add=T, col="green")
sapply(c(-3:3 * 2/7), line.tan, dx=3, p=dy, dp=d2y, k=k, variant=2, col="orange")


####################
####################

####################
### Section B.2. ###
####################

### Mixt Exp-Trig

### y = sin(F(x)) * exp(-F(x))

### Examples:

### y = sin(x^n) * exp(-x^n)

### D(y)
n*x^(n-1)*cos(x^n)*exp(-x^n) - n*x^(n-1)*sin(x^n)*exp(-x^n)
n*x^(n-1)*cos(x^n)*exp(-x^n) - n*x^(n-1) * y
# n*x^(n-1)*cos(x^n)*exp(-x^n) =
dy + n*x^(n-1) * y

### D2(y)
n*(n-1)*x^(n-2)*cos(x^n)*exp(-x^n) +
	- n^2*x^(2*n-2)*sin(x^n)*exp(-x^n) +
	- n^2*x^(2*n-2)*cos(x^n)*exp(-x^n) +
	- n*x^(n-1) * dy - n*(n-1)*x^(n-2) * y
(n-1)*(dy + n*x^(n-1) * y) / x +
	- n*x^(n-1) * (dy + n*x^(n-1) * y) +
	- n*x^(n-1) * dy +
	- n^2*x^(2*n-2) * y - n*(n-1)*x^(n-2) * y

### ODE:
x*d2y - (n-1)*dy + 2*n*x^n*dy + 2*n^2*x^(2*n-1)*y

### Examples
### n = 2
x*d2y - dy + 4*x^2*dy + 8*x^3*y # = 0


### Solution & Plot:
y = function(x, n=2) {
	xn = if(n == 1) x else x^n;
	y = sin(xn) * exp(-xn);
	return(y)
}
dy = function(x, n=2) {
	xn = if(n == 1) x else x^n;
	x.e = exp(-xn);
	dp = n*xn*x.e*(cos(xn) - sin(xn))
	div = x;
	dp = ifelse(div != 0, dp / div,
		if(n > 1) 0 else if(n == 1) 1 else -Inf); # probably correct
	return(dp)
}
d2y = function(x, n=2) {
	# variant = test equations for variants;
	y.x = y(x, n=n)
	dy.x = dy(x, n=n)
	xn = if(n == 1) x else x^n;
	dp = (n-1)*x*dy.x - 2*n*x*xn*dy.x - 2*n^2*xn*xn*y.x;
	div = x*x;
	dp = ifelse(div != 0, dp / div,
		if(n > 1) 1 else if(n == 1) 1 else -Inf); # TODO: correct
	return(dp)
}
### Plot:
n = 2;
x.px = c(-4:4 * 4/7)
curve(y(x, n=n), from= -3, to = 3, ylim=c(-0.5, 0.5))
sapply(x.px, line.tan, dx=2, p=y, dp=dy, n=n)
# wave
curve(dy(x, n=n), add=T, col="green")
sapply(x.px, line.tan, dx=1.25, p=dy, dp=d2y, n=n, col="orange")



####################
####################

####################
### Section B.3. ###
####################

################
### Combined ###
################

### y = I(e^(k*x^n)) * I(e^(-k*x^n))

### ODE:
(d2y + n*k*x^(n-1)*dy - 2)*(d2y - n*k*x^(n-1)*dy - 2) + 4*n^2*k^2*x^(2*n-2)*y = 0


### Solution & Plot:
base.I = function(x, n=2, k=1, low=0) {
	Ip = sapply(x, function(upper) integrate(function(x) exp( k*x^n), lower=low, upper=upper)$value)
	In = sapply(x, function(upper) integrate(function(x) exp(-k*x^n), lower=low, upper=upper)$value)
	return(list(Ip=Ip, In=In))
}
y = function(x, n=2, k=1, I, low=0) {
	if(missing(I)) I = base.I(x, n=n, k=k, low=low);
	# return(I$Ip * I$In)
	# by formula
	dy = dy(x, n=n, k=k, I=I, low=low)
	d2y = d2y(x, n=n, k=k, I=I, low=low)
	# TODO: check k;
	- (d2y + n*k*x^(n-1)*dy - 2)*(d2y - n*k*x^(n-1)*dy - 2) / (4*n^2*k^2*x^(2*n-2))
}
dy = function(x, n=2, k=1, I, low=0) {
	if(missing(I)) I = base.I(x, n=n, k=k, low=low);
	e = exp(k * x^n);
	s = I$In * e + I$Ip / e;
	return(s)
}
d2y = function(x, n=2, k=1, I, low=0) {
	if(missing(I)) I = base.I(x, n=n, k=k, low=low);
	e = exp(k * x^n);
	dx = n*k*x^(n-1);
	s = dx * (I$In * e - I$Ip / e) + 2;
	return(s)
}
### Plot:
n = 2; k = 1;
px = c(-3:3 * 3/5);
curve(y(x, n=n, k=k), from= -2, to = 2, ylim=c(-10,10))
sapply(px, line.tan, dx=3, p=y, dp=dy, n=n, k=k)
# non-sigmoidal
curve(dy(x, n=n, k=k), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, n=n, k=k, col="orange")


### Ex 2:
n = 3.5; k = 1;
px = c(0:5 * 3/8);
curve(y(x, n=n, k=k), from= 0, to = 2, ylim=c(0,10))
sapply(px, line.tan, dx=3, p=y, dp=dy, n=n, k=k)
# non-sigmoidal
curve(dy(x, n=n, k=k), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, n=n, k=k, col="orange")


### Ex 3:
n = 2; k = 3;
px = c(-3:3 * 4/13);
curve(y(x, n=n, k=k), from= -2, to = 2, ylim=c(-10,10))
sapply(px, line.tan, dx=3, p=y, dp=dy, n=n, k=k)
# non-sigmoidal / tan-like
curve(dy(x, n=n, k=k), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, n=n, k=k, col="orange")


#######################

### Combined:
### Extensions

### y = x^m * I(e^(x^n)) * I(e^(-x^n))

### D(y)
m*y / x + x^m*e^(x^n)*In + x^m*e^(-x^n)*Ip
# x*dy =
m*y + x^(m+1)*(e^(x^n)*In + e^(-x^n)*Ip)

### D2(y)
# x*d2y + dy =
m*dy + (m+1)*(x*dy - m*y) / x + 2*x^(m+1) +
	+ n*x^(m+n)*(e^(x^n)*In - e^(-x^n)*Ip)
# x^2*d2y + x*dy =
(2*m+1)*x*dy - m*(m+1)*y + 2*x^(m+2) +
	+ n*x^(m+n+1)*(e^(x^n)*In - e^(-x^n)*Ip)
# x^2*d2y =
2*m*x*dy - m*(m+1)*y + 2*x^(m+2) +
	+ n*x^(m+n+1)*(e^(x^n)*In - e^(-x^n)*Ip)

### Solve liniar system:
T = 2*m*x*dy - m*(m+1)*y + 2*x^(m+2)
e^(x^n)*In  = (n*x^n*(x*dy - m*y) + x^2*d2y - T) / (2*n*x^(m+n+1))
e^(x^-n)*Ip = (n*x^n*(x*dy - m*y) - x^2*d2y + T) / (2*n*x^(m+n+1))

### ODE:
(n*x^n*(x*dy - m*y) + x^2*d2y - T) * (n*x^n*(x*dy - m*y) - x^2*d2y + T) +
	- 4*n^2*x^(m+2*n+2)*y # = 0

### TODO: verify;

### Solution & Plot:
base.I = function(x, n=2, k=1, low=0) {
	Ip = sapply(x, function(x) integrate(function(v) exp(v^n/k), lower=low, upper=x)$value)
	In = sapply(x, function(x) integrate(function(v) exp(-v^n/k), lower=low, upper=x)$value)
	return(list(Ip=Ip, In=In))
}
y = function(x, m=1, n=2, k=1, I, low=0) {
	if(missing(I)) I = base.I(x, n=n, k=k, low=low);
	return(x^m * I$Ip * I$In)
}
dy = function(x, m=1, n=2, k=1, I, low=0) {
	if(missing(I)) I = base.I(x, n=n, k=k, low=low);
	e = exp(x^n / k);
	s = I$In * e + I$Ip / e;
	xm = x^m;
	s = m*xm*I$Ip*I$In + xm*x*s;
	div = x;
	s = ifelse(div != 0, s/div, 0) # TODO
	return(s)
}
d2y = function(x, m=1, n=2, k=1, I, low=0) {
	# TODO: k;
	if(missing(I)) I = base.I(x, n=n, k=k, low=low);
	y.x  =  y(x, m=m, n=n, k=k, I=I, low=low);
	dy.x = dy(x, m=m, n=n, k=k, I=I, low=low);
	xnn = n*x^n; xm = x^(m+2);
	s.sq = (xnn*(x*dy.x - m*y.x))^2 - 4*xnn*xnn*xm*y.x;
	s = sqrt(s.sq) * sign(x); # TODO: sign?
	T = 2*m*x*dy.x - m*(m+1)*y.x + 2*xm;
	s = s + T; div = x^2;
	s = ifelse(div != 0, s/div, 0); # TODO
	return(s)
}
### Plot:
m = 1; n = 2; k = 1;
px = -3:3 * 3/7;
curve(y(x, m=m, n=n, k=k), from= -2, to = 2, ylim=c(-10,10))
sapply(px, line.tan, dx=3, p=y, dp=dy, m=m, n=n, k=k)
# non-sigmoidal
curve(dy(x, m=m, n=n, k=k), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, m=m, n=n, k=k, col="orange")



#######################

########################
### B.4. Non-Liniar: ###
###      Double Exp  ###
########################

### y = e^I, where I = I(e^P(x)) dx;

### Note:
# - for similar variants see file:
#   DE.ODE.Log.R;
#   [I is replaced by log(P(x))^2]

### D(y)
exp(p) * y

### D2(y)
exp(p)*dy + dp*exp(p)*y
(dy)^2 / y + dp*dy

### ODE:
y*d2y - dy^2 - dp*y*dy # = 0

### Examples:
### Case: P(x) = -k*x^n;
y*d2y - dy^2 + k*n*x^(n-1)*y*dy # = 0


### Solution & Plot:
base.I = function(x, n=2, k=1, low=0) {
	I.v = sapply(x, function(x) integrate(function(v) exp(-k * v^n), lower=low, upper=x)$value)
	return(I.v)
}
y = function(x, FUN=base.I, I, ..., low=0) {
	if(missing(I)) I = FUN(x, ..., low=low);
	return(exp(I));
}
dy = function(x, FUN=base.I, I, ..., low=0) {
	if(missing(I)) I = FUN(x, ..., low=low);
	e = exp(-k * x^n);
	return(e * exp(I))
}
d2y = function(x, FUN=base.I, n=2, k=1, low=0) {
	y.x = y(x, FUN=FUN, n=n, k=k, low=low)
	dy.x = dy(x, FUN=FUN, n=n, k=k, low=low)
	v = dy.x^2 - k*n*x^(n-1)*y.x*dy.x;
	div = y.x;
	r = ifelse(div != 0, v/div, 1) # TODO
	return(r)
}
### Plot:
n = 2; k = 1;
px = c(-3:3 * 3/5);
curve(y(x, n=n, k=k), from= -2, to = 2, ylim=c(-0.5,3.5))
sapply(px, line.tan, dx=3, p=y, dp=dy, n=n, k=k)
# gaussian
curve(dy(x, n=n, k=k), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, n=n, k=k, col="orange")


###################
###################

### y = e^(I^2)
# where I = I(e^P(x)) dx;

### Note:
# - for other variants, see file:
#   DE.ODE.Log.R;

### D(y)
2*e^p * I*y

### D2(y)
2*e^p * I*dy + 2*e^(2*p) * y + 2*dp*e^p * I*y
(dy)^2 / y + 2*e^(2*p) * y + dp*dy

### ODE:
y*d2y - (dy)^2 - dp*y*dy - 2*e^(2*p) * y^2 # = 0

### Examples:
### Case: P(x) = - k * x^2;
y*d2y - (dy)^2 + 2*k*x*y*dy - 2*e^(-2*k*x^2) * y^2 # = 0


### Solution & Plot:
base.I = function(x, n=2, k=1, low=0) {
	I.v = sapply(x, function(lim)
		integrate(function(x) exp(-k * x^n), lower=low, upper=lim)$value)
	return(I.v)
}
y = function(x, FUN=base.I, I, ..., low=0) {
	if(missing(I)) I = FUN(x, ..., low=low);
	return(exp(I^2));
}
dy = function(x, FUN=base.I, I, ..., low=0) {
	if(missing(I)) I = FUN(x, ..., low=low);
	e = 2*exp(-k * x^n);
	return(e * I * exp(I^2))
}
d2y = function(x, FUN=base.I, n=2, k=1, low=0) {
	y.x = y(x, FUN=FUN, n=n, k=k, low=low)
	dy.x = dy(x, FUN=FUN, n=n, k=k, low=low)
	v = dy.x^2 - k*n*x^(n-1)*y.x*dy.x + 2*exp(-2*k*x^n) * y.x^2;
	div = y.x;
	r = ifelse(div != 0, v/div, 1) # TODO
	return(r)
}
### Plot:
n = 2; k = 1;
px = c(-3:3 * 3/5);
# inverse gaussian
curve(y(x, n=n, k=k), from= -2, to = 2, ylim=c(-2,3))
sapply(px, line.tan, dx=2, p=y, dp=dy, n=n, k=k)
# 2-wave
curve(dy(x, n=n, k=k), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, n=n, k=k, col="orange")


######################
######################

### Polynomial Extensions
### Y(y) = e^I, where I = I(e^P(x)) dx;

### Example:
### y^2 + b*y = e^I

### Note:
# - for similar variants see file:
#   DE.ODE.Log.R;
#   [I is replaced by log(P(x))^2]

### D(Y(y))
# (2*y + b)*dy
exp(p) * (y^2 + b*y)

### D2(Y(y))
# (2*y + b)*d2y + 2*(dy)^2
exp(p) * (2*y + b)*dy + dp*exp(p)*(y^2 + b*y)
exp(p) * (2*y + b)*dy + dp*(2*y + b)*dy
(2*y + b)^2*dy^2 / (y^2 + b*y) + dp*(2*y + b)*dy
# =>
(y^2 + b*y)*(2*y + b)*d2y + 2*(y^2 + b*y)*dy^2 +
	- (2*y + b)^2*dy^2 - dp*(y^2 + b*y)*(2*y + b)*dy

### ODE:
(y^2 + b*y)*(2*y + b)*d2y - (2*y^2 + 2*b*y + b^2)*dy^2 +
	- dp*(y^2 + b*y)*(2*y + b)*dy # = 0

### Examples:
### Case: P(x) = -k*x^n;
(y^2 + b*y)*(2*y + b)*d2y - (2*y^2 + 2*b*y + b^2)*dy^2 +
	+ k*n*x^(n-1)*(y^2 + b*y)*(2*y + b)*dy # = 0

### Solution & Plot:
base.I = function(x, n=2, k=1, low=0) {
	I.v = sapply(x, function(x) integrate(function(v) exp(-k * v^n), lower=low, upper=x)$value)
	return(I.v)
}
y = function(x, b=2, FUN=base.I, I, ..., low=0, only.y=TRUE) {
	if(missing(I)) I = FUN(x, ..., low=low);
	I = exp(I);
	# y^2 + b*y - I = 0;
	y.v = -b/2 + sqrt(b^2 + 4*I)/2;
	if(only.y) {
		return(y.v);
	}
	return(list(y=y.v, Iexp=I));
}
dy = function(x, b=2, FUN=base.I, I, ..., low=0) {
	if(missing(I)) {
		sol = y(x, b=b, FUN=FUN, ..., low=low, only.y=FALSE);
	} else sol = I;
	# f.args = as.list(match.call());
	f.args = as.list(match.call(expand.dots = FALSE)$'...');
	k = eval(f.args$k, list2env(list(...), parent = parent.frame()));
	n = eval(f.args$n, list2env(list(...), parent = parent.frame()));
	# exp(p) * (y^2 + b*y)
	e = exp(-k * x^n);
	dy.v = e * sol$Iexp;
	div = (2*sol$y + b);
	dy.v = ifelse(div != 0, dy.v/div, Inf) # TODO
	return(dy.v)
}
d2y = function(x, b=2, FUN=base.I, n=2, k=1, low=0) {
	sol = y(x, b=b, FUN=FUN, n=n, k=k, low=low, only.y=FALSE)
	dy.x = dy(x, b=b, FUN=FUN, I=sol, n=n, k=k, low=low)
	y.x = sol$y; yf.x = sol$Iexp; y.2 = 2*y.x + b;
	v = (2*yf.x + b^2)*dy.x^2 +
		- k*n*x^(n-1)*yf.x*y.2*dy.x;
	div = yf.x*y.2;
	r = ifelse(div != 0, v/div, 1) # TODO
	return(r)
}
### Plot:
n = 2; k = 1; b = 3;
px = c(-2:2 * 4/7);
curve(y(x, b=b, n=n, k=k), from= -2, to = 2, ylim=c(-0.3, 1))
sapply(px, line.tan, dx=3, p=y, dp=dy, b=b, n=n, k=k)
# gaussian
curve(dy(x, n=n, k=k), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, n=n, k=k, col="orange")

