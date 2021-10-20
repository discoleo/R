########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Gaussian
###
### draft v.0.4c

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
# - and many more: non-homogeneous variants;


###############
### History ###
###############

### Linear / Non-Linear Gaussian-type

### draft v.0.4c:
# - derived from:
#   y = k * exp(x^n) * I(exp(-x^n)) + F0(x);
### draft v.0.4a - v.0.4b-fix:
# - automatic generation of exponential type ODEs;
# - preparation for extension; [v.0.4a-ext0/ext1]
# - extension & more examples; [v.0.4b & v.0.4b-fix]
### draft v.0.3n - v.0.3n-chk0:
# - derived from: I(exp(y^n)) = P1(x) * P2(y);
#   x*d2y - n*x*y^(n-1)*dy^2 - n*y^n*dy + 2*dy = 0;
# - pseudo-check based on pre-defined y;
### draft v.0.3m - v.0.3m-check:
# - derived from: y = f(x) / (k + I(exp(x^n) dx)):
#   dy + y^2 - n*x^(n-1)*y = 0; [+ check]
# - more variants & partial generalization; [v.0.3m-bis & v.0.3m-gen]
### draft v.0.3j - v.0.3k:
# - derived from:
#   y = x^2*e^(x^3) + x^3*e^(x^2) + F0(x);
#   y = m*x^m*e^(k*x^n) - n*x^n*e^(k*x^m) + F0(x);
### draft v.0.3i:
# - extension to: v.0.3g;
### draft v.0.3h:
# - generalization:
#   y = F1(x) * I(e^(P(x))) * I(e^(-P(x))) + F0(x);
# - complicated example: see specific section;
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
source("Polynomials.Helper.R")


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


### Section B: Mixed

### B.1. Mixed Linear
### y = e^(F1(x)) + e^(F2(x))

### B.2. Mixed Exponential-Trigonometric
### y = sin(F1(x)) * exp(-F1(x))

### B.3. Mixed Non-Linear
### y = Integral(exp(F1(x))) * Integral(exp(-F1(x)))

### B.4. Non-Linear Double Exp
### y = e^I, where I = I(e^P(x)) dx;


####################
####################

### Examples:

### y = e^(-x^2)
# [not run]
# dy = -2*x*y;
### ODE:
d2y + 2*x*dy + 2*y # = 0
### I(D1) by parts =>
# y + 2*x*I(y) - 2*I(I(y)) = 0
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
### Plot:
# b0 == 0;
curve(y(x), from= -3, to = 3, ylim=c(-1/3, 4))
line.tan(c(-3:3 * 3/4), dx=3, p=y, dp=dy)
# sigmoidal
curve(dy(x), add=T, col="green")
line.tan(c(-3:3 * 3/4), dx=3, p=dy, dp=d2y, col="orange")


### y = e^(-x^2) + b0
# [not run]
# dy = -2*x*y + 2*b0*x;
# d2z = y;
d2z + 2*x*dz - 2*z - b0*x^2 # = 0
### Plot:
# using functions defined above;
b0 = -1;
curve(y(x, b0=b0), from= -3, to = 3, ylim=c(-3, 2))
line.tan(c(-3:3 * 3/4), dx=3, p=y, dp=dy, b0=b0)
# sigmoidal
curve(dy(x, b0=b0), add=T, col="green")
line.tan(c(-3:3 * 3/4), dx=3, p=dy, dp=d2y, b0=b0, col="orange")


##################
### y = e^f(x) ###

### y = e^-(k2*x^2 + k1*x)
# dy = - (2*k2*x + k1)*y
### Integration by parts =>
# y = -(2*k2*x + k1)*I(y) + 2*k2*I(I(y));

### ODE:
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
### Plot:
# d2z + (4*x + 1)*dz - 4*z = 0
k = c(2, 1)
curve(y(x, k), from= -3, to = 3)
line.tan(c(-3:3 * 3/4), dx=3, p=y, dp=dy, k=k)
# sigmoidal
curve(dy(x, k), add=T, col="green")
line.tan(c(-3:3 * 3/4), dx=3, p=dy, dp=d2y, k=k, col="orange")



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
### Plot:

### Levels: +3, +2, +1;
curve(y(x), from= -3, to = 2.2)
line.tan(c(-3:3 * 3/4), dx=3, p=y, dp=dp2y)
# *NON*-sigmoidal / seems increasing
curve(dp2y(x), add=T, col="green")
line.tan(c(-3:3 * 3/4), dx=3, p=dp2y, dp=dp1y, col="orange")


### Levels: +4, +3, +2;
# Note: takes very long !!!
# - should be implemented using the compact/derived formulas!
curve(y(x, level=4), from= -3, to = 2.2)
line.tan(c(-3:3 * 3/4), dx=3, p=y, dp=dp3y, level=4)
# *NON*-sigmoidal / seems increasing
curve(dp3y(x), add=T, col="green")
line.tan(c(-3:3 * 3/4), dx=3, p=dp3y, dp=dp2y, col="orange")


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

### Solution & Plot:
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
line.tan(c(0:3 * 3/4), dx=3, p=y, dp=dy)
# sigmoidal
curve(dy(x), add=T, col="green")
line.tan(c(0:3 * 3/4), dx=3, p=dy, dp=d2y, col="orange")


################

### y = k * exp(x^n) * I(exp(-x^n)) + F0(x)

### D =>
dy - n*x^(n-1)*y + n*x^(n-1)*f - k - df # = 0

### D2 =>
d2y - n*x^(n-1)*dy - n*(n-1)*x^(n-2)*y +
	+ n*(n-1)*x^(n-2)*f + n*x^(n-1)*df - d2f # = 0
# Variant 1:
d2y - n*x^(n-1)*(n*x^(n-1)*y - n*x^(n-1)*f + df + k) - n*(n-1)*x^(n-2)*y +
	+ n*(n-1)*x^(n-2)*f + n*x^(n-1)*df - d2f # = 0
d2y - n^2*x^(2*n-2)*y - n*(n-1)*x^(n-2)*y +
	+ n^2*x^(2*n-2)*f + n*(n-1)*x^(n-2)*f - k*n*x^(n-1) - d2f # = 0

### Special Cases:
### n = 2
d2y - (4*x^2 + 2)*y + 4*x^2*f + 2*f - 2*k*x - d2f# = 0

### Solution & Plot:
y = function(x, n=2, k=1, f=NULL) {
	fy = function(x) exp(-x^n);
	y = sapply(x, function(up) integrate(fy, lower=0, upper=up)$value);
	y = k * exp(x^n) * y;
	if( ! is.null(f)) {
		fx = sapply(x, function(x) eval.pm(f, x));
		y  = y + fx;
	}
	return(y)
}
dy = function(x, n=2, k=1, f=NULL, y=NULL) {
	yx = if( ! is.null(y)) y else y(x, n=n, k=k, f=f);
	xn1 = n*x^(n-1);
	dp  = xn1*yx + k; # + df;
	if( ! is.null(f)) {
		fx = sapply(x, function(x) eval.pm(f, x));
		dp = dp - xn1*fx;
		df = dp.pm(f);
		dfx = sapply(x, function(x) eval.pm(df, x));
		dp = dp + dfx;
	}
	return(dp)
}
d2y = function(x, n=2, k=1, f=NULL) {
	yx = y(x, n=n, k=k, f=f);
	# dyx = dy(x, n=n, k=k, f=f, y=yx);
	xn2 = x^(n-2); xn1 = x*xn2;
	d2p = n^2*xn1^2*yx + n*(n-1)*xn2*yx;
	if( ! is.null(f)) {
		fx = sapply(x, function(x) eval.pm(f, x));
		d2p = d2p - n*(n-1)*xn2*fx - n^2*xn1^2*fx + k*n*xn1;
		df  = dp.pm(f);
		d2f = dp.pm(df);
		d2fx = sapply(x, function(x) eval.pm(d2f, x));
		d2p = d2p + d2fx;
	}
	return(d2p)
}
### Plot:
n = 2; k = 2;
f = toPoly.pm("x^2 - 3*x - 3")
px = c(0, rep(1,5)) +  (0:5)*1/5;
curve(y(x, n=n, k=k, f=f), from= 0, to = 2, ylim=c(0, 100))
line.tan(px, dx=3, p=y, dp=dy, n=n, k=k, f=f)
#
curve(dy(x, n=n, k=k, f=f), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, n=n, k=k, f=f, col="orange")


### Ex 2:
n = 2; k = -2;
f = toPoly.pm("x^2 - 3*x - 3")
px = c(0, rep(1,5)) +  (0:5)*1/5;
curve(y(x, n=n, k=k, f=f), from= 0, to = 2)
line.tan(px, dx=3, p=y, dp=dy, n=n, k=k, f=f)
#
curve(dy(x, n=n, k=k, f=f), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, n=n, k=k, f=f, col="orange")


################
################

### y = e^(x^n) / (k + I(e^(x^n) dx));

### D =>
dy + y^2 - n*x^(n-1)*y # = 0

### D2 =>
d2y + 2*y*dy - n*x^(n-1)*dy - n*(n-1)*x^(n-2)*y # = 0
# variant 1:
x*d2y + 2*x*y*dy - n*x^n*dy - (n-1)*dy - (n-1)*y^2 # = 0
# variant 2:
y*d2y + y^2*dy - dy^2 - n*(n-1)*x^(n-2)*y^2 # = 0
# variant 3:
x*y^2*dy + x*dy^2 - n*x^n*y*dy - (n-1)*y*dy  - (n-1)*y^3 + n*(n-1)*x^(n-1)*y^2 # = 0

### Example:
n = 2
### V1:
x*d2y + 2*x*y*dy - 2*x^2*dy - dy  - y^2 # = 0
### V3:
x*y^2*dy + x*dy^2 - 2*x^2*y*dy - y*dy  - y^3 + 2*x*y^2 # = 0

### Solution:
y = function(x, n=2, k=1, from=0) {
	y = exp(x^n);
	div = sapply(x, function(xu)
		integrate(function(x) exp(x^n), lower=from, upper=xu)$value);
	y = y / (k + div);
	return(y);
}
dy = function(x, n=2, k=1, from=0, y=NULL) {
	yx = if( ! is.null(y)) y else yx = y(x, n=n, k=k, from=from);
	xn = if(n != 1) x^(n-1) else 1;
	dy = - yx^2 + n*xn*yx;
	return(dy);
}
d2y = function(x, n=2, k=1, from=0) {
	yx  = y(x, n=n, k=k, from=from);
	dyx = dy(x, n=n, k=k, from=from, y=yx);
	d2y = (2*x*yx - n*x^n - (n-1))*dyx  - (n-1)*yx^2;
	d2y = - d2y / x;
	d2y[x == 0] = 0; # TODO
	return(d2y);
}
### Plot:
n = 2
px = c(-0.92, -0.875, -0.75, -0.65, -0.35);
curve(y(x, n=n), from= -1.25, to = 0.75, ylim=c(-150, 75))
line.tan(px, dx=3, p=y, dp=dy, n=n)
# sigmoidal
curve(dy(x, n=n), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, n=n, col="orange")



### Generalisation
### (y^2 + b1*y + b0) * (k + I(e^(x^n) dx)) = e^(x^n);

### D =>
(2*y + b1)*dy + (y^2 + b1*y + b0)^2 - n*x^(n-1)*(y^2 + b1*y + b0) # = 0

### D2 =>
(2*y + b1)*d2y + 2*dy^2 + 2*(y^2 + b1*y + b0)*(2*y + b1)*dy +
	- n*x^(n-1)*(2*y + b1)*dy - n*(n-1)*x^(n-2)*(y^2 + b1*y + b0) # = 0
# Variant 1:
x*(2*y + b1)*d2y + 2*x*dy^2 + 2*x*(y^2 + b1*y + b0)*(2*y + b1)*dy +
	- n*x^n*(2*y + b1)*dy - (n-1)*(2*y + b1)*dy - (n-1)*(y^2 + b1*y + b0)^2 # = 0

# TODO


#######################
#######################

########################
### Section B: Mixed ###
########################

################
### Section B.1:

### Linear Combinations

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

### Partial Extension
### y = x^2*e^(x^3) + x^3*e^(x^2) + F0(x)

### D(y)
(3*x^4 + 2*x)*e^(x^3) + (2*x^4 + 3*x^2)*e^(x^2) + df0

### Solve liniar system:
div = (2*x^4 + 3*x^2)*x^2 - x^3*(3*x^4 + 2*x)
div = -3*x^7 + 2*x^6 + x^4;
e^(x^3) = ((2*x^4 + 3*x^2)*(y - f0) - x^3*(dy - df0)) / div;
e^(x^2) = ((3*x^4 + 2*x)*(y - f0) - x^2*(dy - df0)) / - div;

### D2(y)
(3*x^2*(3*x^4 + 2*x) + 12*x^3 + 2)*e^(x^3) +
	+ (2*x*(2*x^4 + 3*x^2) + 8*x^3 + 6*x)*e^(x^2) + d2f0
(9*x^6 + 18*x^3 + 2)*e^(x^3) +
	+ (4*x^5 + 14*x^3 + 6*x)*e^(x^2) + d2f0
### div * d2y =
(9*x^6 + 18*x^3 + 2)*((2*x^4 + 3*x^2)*(y - f0) - x^3*(dy - df0)) +
	- (4*x^5 + 14*x^3 + 6*x)*((3*x^4 + 2*x)*(y - f0) - x^2*(dy - df0)) + div*d2f0

### Examples:
### f = x^3
div = -3*x^7 + 2*x^6 + x^4;
### div * d2y =
(9*x^6 + 18*x^3 + 2)*((2*x^4 + 3*x^2)*y - x^3*dy) +
	- (4*x^5 + 14*x^3 + 6*x)*((3*x^4 + 2*x)*y - x^2*dy) +
	- 2*x^7*(9*x^6 + 18*x^3 + 2) + x^3*(3*x^4 - x)*(4*x^5 + 14*x^3 + 6*x) + 6*x*(-3*x^7 + 2*x^6 + x^4)
(4*x^3 + 14*x^5 - 18*x^6 + 4*x^7 - 9*x^9)*dy +
	+ (-6*x^2 - 24*x^4 + 36*x^5 - 8*x^6 - 6*x^7 + 27*x^8 - 12*x^9 + 18*x^10)*y +
	- 6*x^7 - 4*x^9 + 6*x^10 + 12*x^12 - 18*x^13


### Solution & Plot:
y = function(x, a=c(1, 1)) {
	x2 = x^2; x3 = x^3;
	y = a[1]*x2*exp(x3) +
		a[2]*x3*exp(x2) + x3;
	return(y)
}
dy = function(x, a=c(1, 1)) {
	x2 = x^2; x4 = x2*x2;
	dy.v = a[1]*(3*x4 + 2*x)*exp(x^3) +
		a[2]*(2*x4 + 3*x2)*exp(x2) + 3*x2;
	return(dy.v)
}
d2y = function(x, a=c(1, 1)) {
	y.x = y(x, a=a)
	dy.x = dy(x, a=a)
	dp = (4*x^3 + 14*x^5 - 18*x^6 + 4*x^7 - 9*x^9)*dy.x +
		+ (-6*x^2 - 24*x^4 + 36*x^5 - 8*x^6 - 6*x^7 + 27*x^8 - 12*x^9 + 18*x^10)*y.x +
		- 6*x^7 - 4*x^9 + 6*x^10 + 12*x^12 - 18*x^13;
	div = -3*x^7 + 2*x^6 + x^4;
	dp = ifelse(div != 0, dp / div, d2y(x + 1E-5, a=a)); # TODO
	return(dp)
}
### Plot:
a = c(1, 1) # a[] has NO effect on eq of D2;
px = c(-3:3 * 2/7);
curve(y(x, a=a), from= -1.5, to = 1.5, ylim=c(-10, 10))
sapply(px * 1.2, line.tan, dx=3, p=y, dp=dy, a=a)
#
curve(dy(x, a=a), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, a=a, col="orange")


### Ex 2:
a = c(1.5, -2) # a[] has NO effect on eq of D2;
px = c(-3:3 * 2/7);
curve(y(x, a=a), from= -1.5, to = 1.5, ylim=c(-12, 8))
sapply(px * 1.2, line.tan, dx=2.2, p=y, dp=dy, a=a)
#
curve(dy(x, a=a), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, a=a, col="orange")


### Partial Extension
### y = m*x^m*e^(k*x^n) - n*x^n*e^(k*x^m) + F0(x)

### D(y)
m^2*x^(m-1)*e^(k*x^n) - n^2*x^(n-1)*e^(k*x^m) + df0

### Solve Liniar system:
e^(k*x^n) = (n*(y - f0) - x*(dy - df0)) / (m*(n-m)*x^m)
e^(k*x^m) = (m*(y - f0) - x*(dy - df0)) / (n*(n-m)*x^n)

### D2(y)
m^2*(m-1)*x^(m-2)*e^(k*x^n) - n^2*(n-1)*x^(n-2)*e^(k*x^m) +
	+ k*m^2*n*x^(m+n-2)*e^(k*x^n) - k*m*n^2*x^(m+n-2)*e^(k*x^m) + d2f0
### (n-m)*x^2*D2(y) =
m*(m-1)*(n*y - x*dy) - n*(n-1)*(m*y - x*dy) +
	+ k*m*n*x^n*(n*y - x*(dy - df0)) - k*m*n*x^m*(m*y - x*(dy - df0)) +
	- m*n*(m-n)*f0 - k*m*n*(n*x^n - m*x^m)*f0 +
	+ (m-n)*(m+n-1)*x*df0 +
	+ (n-m)*x^2*d2f0
k*m*n*x*(x^m - x^n)*dy + (n-m)*(n + m - 1)*x*dy +
	+ k*m*n*(n*x^n - m*x^m)*y + m*n*(m-n)*y +
	- k*m*n*(n*x^n - m*x^m)*f0 - m*n*(m-n)*f0 +
	+ k*m*n*x*(x^n - x^m)*df0 + (m-n)*(m + n - 1)*x*df0 + (n-m)*x^2*d2f0

### Examples:
### f0 = x^2
### (n-m)*x^2*D2(y) =
k*m*n*x*(x^m - x^n)*dy + (n-m)*(n + m - 1)*x*dy +
	+ k*m*n*(n*x^n - m*x^m)*y + m*n*(m-n)*y +
	- k*m*n*(n*x^n - m*x^m)*x^2 + 2*k*m*n*(x^n - x^m)*x^2 +
	+ 2*(m-n)*(m+n - 2)*x^2 - m*n*(m-n)*x^2;


### Solution & Plot:
y = function(x, m=c(2,3), k=1, a=c(1, 1)) {
	xm = x^m[1]; xn = x^m[2];
	y = a[1]*m[1]*xm*exp(k*xn) - a[2]*m[2]*xn*exp(k*xm) + x^2;
	return(y)
}
dy = function(x, m=c(2,3), k=1, a=c(1, 1)) {
	# TODO: a;
	xm = x^m[1]; xn = x^m[2];
	dy.v = m[1]^2*xm*exp(k*xn) - m[2]^2*xn*exp(k*xm) + 2*x^2; # x*df0
	div = x;
	dy.v = ifelse(div != 0, dy.v / div, dy(x + 1E-3, m=m, k=k, a=a)); # TODO
	return(dy.v)
}
d2y = function(x, m=c(2,3), k=1, a=c(1, 1)) {
	y.x = y(x, m=m, k=k, a=a)
	dy.x = dy(x, m=m, k=k, a=a)
	xm = x^m[1]; xn = x^m[2]; x2 = x^2;
	n = m[2]; m = m[1];
	mn = m*n; kmn = k*mn;
	dp = kmn*x*(xm - xn)*dy.x + (n-m)*(n + m - 1)*x*dy.x +
		+ kmn*(n*xn - m*xm)*y.x + mn*(m-n)*y.x +
		- kmn*(n*xn - m*xm)*x2 +
		+ 2*kmn*(xn - xm)*x2 + 2*(m-n)*(n + m - 2)*x2 - mn*(m-n)*x2;
	div = (n-m)*x2;
	dp = ifelse(div != 0, dp / div, d2y(x + 1E-2, m=c(m,n), k=k, a=a)); # TODO
	return(dp)
}
### Plot:
m = c(2, 3); k = 1;
a = c(1, 1) # a[] has NO effect on eq of D2;
px = c(-3:5 * 1/7);
curve(y(x, m=m, k=k, a=a), from= -1.5, to = 1.5, ylim=c(-10, 10))
sapply(px * 2, line.tan, dx=3, p=y, dp=dy, m=m, k=k, a=a)
#
curve(dy(x, m=m, k=k, a=a), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, m=m, k=k, a=a, col="orange")


### Ex 2:
m = c(2, 3); k = -1;
a = c(1, 1) # a[] has NO effect on eq of D2;
px = c(-2:1 * 2/7, 0.8, 1.3);
curve(y(x, m=m, k=k, a=a), from= -1.5, to = 1.8, ylim=c(-10, 10))
sapply(px * 1.2, line.tan, dx=3, p=y, dp=dy, m=m, k=k, a=a)
#
curve(dy(x, m=m, k=k, a=a), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, m=m, k=k, a=a, col="orange")



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

### y = x^m * I(e^(k*x^n)) * I(e^(-k*x^n))

### D(y)
m*y / x + x^m*e^(k*x^n)*In + x^m*e^(-k*x^n)*Ip
# x*dy =
m*y + x^(m+1)*(e^(k*x^n)*In + e^(-k*x^n)*Ip)

### D2(y)
# x*d2y + dy =
m*dy + (m+1)*(x*dy - m*y) / x + 2*x^(m+1) +
	+ n*k*x^(m+n)*(e^(k*x^n)*In - e^(-k*x^n)*Ip)
# x^2*d2y + x*dy =
(2*m+1)*x*dy - m*(m+1)*y + 2*x^(m+2) +
	+ n*k*x^(m+n+1)*(e^(k*x^n)*In - e^(-k*x^n)*Ip)
# x^2*d2y =
2*m*x*dy - m*(m+1)*y + 2*x^(m+2) +
	+ n*k*x^(m+n+1)*(e^(k*x^n)*In - e^(-k*x^n)*Ip)

### Solve liniar system:
T = 2*m*x*dy - m*(m+1)*y + 2*x^(m+2)
e^( k*x^n)*In = (n*k*x^n*(x*dy - m*y) + x^2*d2y - T) / (2*n*k*x^(m+n+1))
e^(-k*x^n)*Ip = (n*k*x^n*(x*dy - m*y) - x^2*d2y + T) / (2*n*k*x^(m+n+1))

### ODE:
(n*k*x^n*(x*dy - m*y) + x^2*d2y - T) * (n*k*x^n*(x*dy - m*y) - x^2*d2y + T) +
	- 4*n^2*k^2*x^(m+2*n+2)*y # = 0


### Solution & Plot:
base.I = function(x, n=2, k=1, low=0) {
	Ip = sapply(x, function(upper) integrate(function(x) exp( k*x^n), lower=low, upper=upper)$value)
	In = sapply(x, function(upper) integrate(function(x) exp(-k*x^n), lower=low, upper=upper)$value)
	return(list(Ip=Ip, In=In))
}
y = function(x, m=1, n=2, k=1, I, low=0) {
	if(missing(I)) I = base.I(x, n=n, k=k, low=low);
	return(x^m * I$Ip * I$In)
}
dy = function(x, m=1, n=2, k=1, I, low=0) {
	if(missing(I)) I = base.I(x, n=n, k=k, low=low);
	e = exp(k * x^n);
	s = I$In * e + I$Ip / e;
	xm = x^m;
	s = m*xm*I$Ip*I$In + xm*x*s;
	div = x;
	s = ifelse(div != 0, s/div, 0) # TODO
	return(s)
}
d2y = function(x, m=1, n=2, k=1, I, low=0) {
	if(missing(I)) I = base.I(x, n=n, k=k, low=low);
	y.x  =  y(x, m=m, n=n, k=k, I=I, low=low);
	dy.x = dy(x, m=m, n=n, k=k, I=I, low=low);
	xnn = n*k*x^n; xm = x^(m+2);
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


### Ex 2:
m = 1; n = 2; k = 3;
px = -3:3 * 2/7;
curve(y(x, m=m, n=n, k=k), from= -2, to = 2, ylim=c(-10,10))
sapply(px * 1.3, line.tan, dx=3, p=y, dp=dy, m=m, n=n, k=k)
# non-sigmoidal
curve(dy(x, m=m, n=n, k=k), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, m=m, n=n, k=k, col="orange")


#######################

### Combined:
### Fully Generalized

### y = F1(x) * I(e^(P(x))) * I(e^(-P(x))) + F0(x)

### D(y)
df1 / f1 * (y - f0) + df0 + f1*(e^p*In + e^(-p)*Ip)
# f1*dy =
df1*y - f0*df1 + f1*df0 + f1^2*(e^p*In + e^(-p)*Ip)

### D2(y)
# f1*d2y =
d2f1*y - f0*d2f1 + f1*d2f0 + 2*f1^2 +
	+ 2*f1*df1 / f1^2*(f1*dy - df1*y + f0*df1 - f1*df0) +
	+ f1^2*dp*(e^p*In - e^(-p)*Ip);

### Solve liniar system:
f1^2*(e^p*In + e^(-p)*Ip) # =
	f1*dy - df1*y + f0*df1 - f1*df0;
f1^3*dp*(e^p*In - e^(-p)*Ip) # =
	f1^2*d2y - 2*f1*df1*dy + 2*df1^2*y - f1*d2f1*y +
	+ f0*f1*d2f1 - f1^2*d2f0 - 2*f1^3 - 2*f0*df1^2 + 2*f1*df0*df1;


### Examples:
# f0 = x + b0; df0 = 1;
# f1 = 1/x; df1 = -1/x^2; d2f1 = 2/x^3;
(e^p*In + e^(-p)*Ip) # =
	x*dy + y - 2*x - b0;
x*dp*(e^p*In - e^(-p)*Ip) # =
	x^2*d2y + 2*x*dy - 4*x;
###
2*x*dp*e^(-p)*Ip = x*dp*(x*dy + y - 2*x - b0) - x^2*d2y - 2*x*dy + 4*x;
2*x*dp*e^( p)*In = x*dp*(x*dy + y - 2*x - b0) + x^2*d2y + 2*x*dy - 4*x;

### ODE:
4*x^3*dp^2*y # =
	(x*dp*(x*dy + y - 2*x - b0) - x^2*d2y - 2*x*dy + 4*x) *
	(x*dp*(x*dy + y - 2*x - b0) + x^2*d2y + 2*x*dy - 4*x) + 4*x^3*(x + b0)*dp^2;
### p = k*x^n
4*k^2*n^2*x^(2*n+1)*y # =
	(k*n*x^n*(x*dy + y - 2*x - b0) - x^2*d2y - 2*x*dy + 4*x) *
	(k*n*x^n*(x*dy + y - 2*x - b0) + x^2*d2y + 2*x*dy - 4*x) + 4*k^2*n^2*x^(2*n+1)*(x + b0);


### Solution & Plot:
base.I = function(x, n=2, k=1, low=0) {
	Ip = sapply(x, function(upper) integrate(function(x) exp( k*x^n), lower=low, upper=upper)$value)
	In = sapply(x, function(upper) integrate(function(x) exp(-k*x^n), lower=low, upper=upper)$value)
	return(list(Ip=Ip, In=In))
}
y = function(x, m=-1, n=2, k=1, b0=1, I, low=0) {
	if(missing(I)) I = base.I(x, n=n, k=k, low=low);
	if(m > 0) {
		r = x^m * I$Ip * I$In;
	} else r = ifelse(x != 0, x^m * I$Ip * I$In, 0);
	r = r + x + b0;
	return(r)
}
dy = function(x, m=-1, n=2, k=1, b0=1, I, low=0) {
	if(missing(I)) I = base.I(x, n=n, k=k, low=low);
	e = exp(k * x^n);
	s = I$In * e + I$Ip / e;
	s = 2*x + b0 + s - y(x, m=m, n=n, k=k, b0=b0, I=I, low=low);
	div = x; # m = -1;
	s = ifelse(div != 0, s/div, 2) # TODO
	return(s)
}
d2y = function(x, m=-1, n=2, k=1, b0=1, I, low=0) {
	# TODO: k;
	if(missing(I)) I = base.I(x, n=n, k=k, low=low);
	y.x  =  y(x, m=m, n=n, k=k, b0=b0, I=I, low=low);
	dy.x = dy(x, m=m, n=n, k=k, b0=b0, I=I, low=low);
	xnn = n*x^n; xm = x^(m+2);
	# 4*k^2*n^2*x^(2*n+1)*y = (k*n*x^n*(x*dy + y - 2*x - b0) - x^2*d2y - 2*x*dy + 4*x) *
	#  (k*n*x^n*(x*dy + y - 2*x - b0) + x^2*d2y + 2*x*dy - 4*x) + 4*k^2*n^2*x^(2*n+1)*(x + b0)
	s.sq = (k*xnn*(x*dy.x + y.x - 2*x - b0))^2 +
		- 4*k^2*xnn*xnn*x * (y.x - x - b0);
	s = sqrt(s.sq) * sign(x); # TODO: sign?
	T = 2*x*dy.x - 4*x;
	s = s - T; div = x^2;
	s = ifelse(div != 0, s/div, 0); # TODO
	return(s)
}
### Plot:
n = 2; k = 1; b0 = 1; # m = -1; is fixed;
px = -3:3 * 4/7;
curve(y(x, n=n, k=k, b0=b0), from= -2, to = 2, ylim=c(-6,10))
sapply(px, line.tan, dx=3, p=y, dp=dy, n=n, k=k, b0=b0)
# non-sigmoidal
curve(dy(x, n=n, k=k, b=b0), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, n=n, k=k, b=b0, col="orange")


### Ex 2:
n = 2; k = 1; b0 = 3; # m = -1; is fixed;
px = -3:3 * 4/7;
curve(y(x, n=n, k=k, b0=b0), from= -2, to = 2, ylim=c(-6,10))
sapply(px, line.tan, dx=3, p=y, dp=dy, n=n, k=k, b0=b0)
# non-sigmoidal
curve(dy(x, n=n, k=k, b=b0), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, n=n, k=k, b=b0, col="orange")



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
line.tan(px, dx=3, p=y, dp=dy, b=b, n=n, k=k)
# gaussian
curve(dy(x, n=n, k=k), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, n=n, k=k, col="orange")


##################
##################

### I(exp(y^n) dy) = P1(x) * P2(y)

### Example:
### I(exp(y^n) dy) = x*y

### D =>
exp(y^n) - x*dy - y # = 0

### D2 =>
n*y^(n-1)*exp(y^n)*dy - x*d2y - 2*dy # = 0
n*y^(n-1)*(x*dy + y)*dy - x*d2y - 2*dy # = 0
# =>
x*d2y - n*x*y^(n-1)*dy^2 - n*y^n*dy + 2*dy # = 0

### Example:
# n = 2
x*d2y - 2*x*y*dy^2 - 2*y^2*dy + 2*dy # = 0

### Variant:
### Integration by parts =>
x*dy - 2*I(x*dy*y*dy) - 2/3*y^3 + 2*y - y # = 0
x*dy - x*y^2*dy + I((x*d2y + dy) * y^2) - 2/3*y^3 + y # = 0
x*dy - x*y^2*dy + I(x*y^2*d2y) - 1/3*y^3 + y # = 0
# [...]

# TODO: check
# Q: How ???


### Special Example:
# - based on a particular case:
#   n = 2; y = sqrt(x);
# => I(exp(y^n) dy) = exp(x) / sqrt(x) * y;

### D =>
x*sqrt(x)*exp(y^n) - x*exp(x)*dy - (x - 1/2)*exp(x)*y # = 0

### Solution:
y  = function(x, n=2) rootn(x, n=n);
dy = function(x, n=2) {
	yx = y(x, n=n);
	xn = rootn(x, n=n); ex = exp(x);
	dy = x*xn*exp(yx^n) - (x - 1/n)*ex*yx;
	div = x*ex;
	dy = dy / div;
	return(dy);
}
### Plot:
n = 2;
px = c(1E-2, (1:5) / 5);
curve(y(x, n=n), from = 0, to = 2)
line.tan(px, dx=3, p=y, dp=dy, n=n)


######################
######################

solve.Dexp = function(pl1, pl2, p0=NULL) {
	dp1 = dp.exp.pm(pl1, xn="x");
	dp2 = dp.exp.pm(pl2, xn="x");
	### Substitution
	pE1  = dp2$Poly; pE1$y = 1;
	pE1b = pl2$Poly; pE1b$dy = 1;
	pE1  = diff.pm(pE1, pE1b);
	pE2  = dp1$Poly; pE2$y = 1;
	pE2b = pl1$Poly; pE2b$dy = 1;
	pE2  = diff.pm(pE2, pE2b);
	#
	div  = diff.pm(
		mult.pm(pl1$Poly, dp2$Poly),
		mult.pm(pl2$Poly, dp1$Poly));
	if( ! is.null(p0)) {
		pE1 = diff.pm(pE1, mult.pm(p0, dp2$Poly));
		pE2 = diff.pm(pE2, mult.pm(p0, dp1$Poly));
		# dp0
		dp0 = dp.pm(p0, xn="x");
		pE1 = sum.pm(pE1, mult.pm(dp0, pl2$Poly));
		pE2 = sum.pm(pE2, mult.pm(dp0, pl1$Poly));
		# d2p0
		d2p0 = dp.pm(dp0, xn="x");
		if(is.data.frame(d2p0)) {
			d2p0 = mult.pm(div, d2p0);
		} else d2p0 = NULL;
	} else d2p0 = NULL;
	### D2
	d2p1 = dp.exp.pm(dp1);
	d2p2 = dp.exp.pm(dp2);
	# Subst =>
	div$d2y = 1;
	pE1 = mult.pm(d2p1$Poly, pE1);
	pE2 = mult.pm(d2p2$Poly, pE2);
	pR = diff.pm(pE1, pE2);
	pR = diff.pm(div, pR);
	if( ! is.null(d2p0)) {
		pR = diff.pm(pR, d2p0);
	}
	return(pR)
}
solveP.Dexp = function(pl1, pl2, p0=NULL) {
	pR = solve.Dexp(pl1, pl2, p0);
	pR = sort.pm(pR, xn=c("d2y", "dy", "y", "x"), 10:13)
	if(pR$coeff[1] < 0) pR$coeff = - pR$coeff;
	pp = print.pm(pR, leading=c("y", "dy", "d2y"), do.sort=FALSE);
	print(pp);
	invisible(pR);
}

### Generic Solution:
eval.dexp = function(x, pl1, pl2, p0=NULL) {
	y1  = sapply(x, function(x) eval.pm(pl1$Poly, x));
	y2  = sapply(x, function(x) eval.pm(pl2$Poly, x));
	yE1 = sapply(x, function(x) eval.pm(pl1$Exp, x));
	yE2 = sapply(x, function(x) eval.pm(pl2$Exp, x));
	y = y1 * exp(yE1) + y2 * exp(yE2);
	if( ! is.null(p0)) {
		y0 = sapply(x, function(x) eval.pm(p0, x));
		y  = y + y0;
	}
	return(y);
}
y = function(x) {
	return(eval.dexp(x, pl1, pl2, p0));
}
dy = function(x) {
	dp1 = dp.exp.pm(pl1)
	dp2 = dp.exp.pm(pl2)
	dp0 = if(is.null(p0)) NULL else dp.pm(p0);
	return(eval.dexp(x, dp1, dp2, p0=dp0));
}
# d2y = is specific;


### Examples:

### Ex 1:
n = 2
m = 2
# too long!
# p1 = toPoly.pm("x^n + b11*x + b10")
# p2 = toPoly.pm("x^n + b21*x + b20")
p0 = NULL
p1 = toPoly.pm("x^n + 3*x + 3")
p2 = toPoly.pm("x^n - 2*x + 3")
pE1 = toPoly.pm("x^m + 2*x")
pE2 = toPoly.pm("x^m - 2*x")
pl1 = list(Exp=pE1, Poly=p1)
pl2 = list(Exp=pE2, Poly=p2)
#
solveP.Dexp(pl1, pl2, p0=p0)

### ODE:
(4*x^4 + 4*x^3 - 5*x^2 + 12*x + 51)*d2y +
	- (16*x^5 + 16*x^4 - 4*x^3 + 60*x^2 + 194*x + 12)*dy +
	+ (16*x^6 + 16*x^5 - 12*x^4 + 48*x^3 + 278*x^2 - 36*x - 508)*y


### Solution:
# pl1, pl2, p0: see above;
d2y = function(x) {
	yx = y(x);
	dyx = dy(x);
	div = (4*x^4 + 4*x^3 - 5*x^2 + 12*x + 51);
	d2y = (16*x^5 + 16*x^4 - 4*x^3 + 60*x^2 + 194*x + 12)*dyx +
		- (16*x^6 + 16*x^5 - 12*x^4 + 48*x^3 + 278*x^2 - 36*x - 508)*yx;
	d2y = ifelse(div != 0, d2y/div, 0); # TODO: check!
	return(d2y)
}
### Plot:
p0 = NULL;
px = c((-5:5) / 5);
curve(y(x), from = -1, to = 1, ylim=c(-60, 150))
line.tan(px, dx=3, p=y, dp=dy)
# global minimum
curve(dy(x), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, col="orange")

# TODO: more checks;

#########
### Ex 2:
p0 = NULL
p1 = toPoly.pm("x")
p2 = toPoly.pm("x")
pE1 = toPoly.pm("x^2 + 2*x")
pE2 = toPoly.pm("x^2 - 2*x")
pl1 = list(Exp=pE1, Poly=p1)
pl2 = list(Exp=pE2, Poly=p2)
#
solveP.Dexp(pl1, pl2, p0=p0)

### ODE:
x^2*d2y - 2*x*(2*x^2 + 1)*dy + 2*(2*x^4 - x^2 + 1)*y # = 0

### Plot:
# pl1, pl2, p0: see above;
d2y = function(x) {
	yx = y(x);
	dyx = dy(x);
	x2 = x^2;
	div = x2;
	d2y = 2*x*(2*x2 + 1)*dyx - 2*(2*x2^2 - x2 + 1)*yx;
	d2y = ifelse(div != 0, d2y/div, 0); # TODO: check!
	return(d2y)
}
### Plot:
p0 = NULL;
px = c((-3:3) * 3/9);
curve(y(x), from = -1, to = 1, ylim=c(-30, 30))
line.tan(px, dx=3, p=y, dp=dy)
# Inflexion
curve(dy(x), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, col="orange")


#########
### Ex 3:
p0 = toPoly.pm("x^2")
p1 = toPoly.pm("x")
p2 = toPoly.pm("x")
pE1 = toPoly.pm("x^2 + 2*x")
pE2 = toPoly.pm("x^2 - 2*x")
pl1 = list(Exp=pE1, Poly=p1)
pl2 = list(Exp=pE2, Poly=p2)
#
solveP.Dexp(pl1, pl2, p0=p0)

### ODE:
x^2*d2y - 2*x*(2*x^2 + 1)*dy + 2*(2*x^4 - x^2 + 1)*y - 4*x^6 + 10*x^4 # = 0

### Plot:
# pl1, pl2, p0: see above;
d2y = function(x) {
	yx = y(x);
	dyx = dy(x);
	x2 = x^2;
	div = x2;
	d2y = 2*x*(2*x2 + 1)*dyx - 2*(2*x2^2 - x2 + 1)*yx + 4*x2^3 - 10*x2^2;
	d2y = ifelse(div != 0, d2y/div, 3.5); # TODO: correct!
	return(d2y)
}
### Plot:
px = c((-3:3) * 3/9);
curve(y(x), from = -1, to = 1, ylim=c(-20, 30))
line.tan(px, dx=3, p=y, dp=dy)
# Inflexion
curve(dy(x), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, col="orange")


#########
### Ex 4:
p0 = toPoly.pm("x^2")
p1 = toPoly.pm("x + 1")
p2 = toPoly.pm("x - 1")
pE1 = toPoly.pm("x^2 + 2*x")
pE2 = toPoly.pm("x^2 - 2*x")
pl1 = list(Exp=pE1, Poly=p1)
pl2 = list(Exp=pE2, Poly=p2)
#
solveP.Dexp(pl1, pl2, p0=p0)

### ODE:
(2*x^2 - 3)*d2y - 8*x*(x^2 - 1)*dy + 2*(4*x^4 - 8*x^2 + 15)*y +
	- 8*x^6 + 32*x^4 - 50*x^2 + 6 # = 0

### Plot:
# pl1, pl2, p0: see above;
d2y = function(x) {
	yx = y(x);
	dyx = dy(x);
	x2 = x^2;
	div = 2*x2 - 3;
	d2y = 8*x*(x2 - 1)*dyx - 2*(4*x2^2 - 8*x2 + 15)*yx +
		+ 8*x2^3 - 32*x2^2 + 50*x2 - 6;
	d2y = ifelse(div != 0, d2y/div, 0); # TODO: check!
	return(d2y)
}
### Plot:
px = c((-3:3) * 2/9);
curve(y(x), from = -0.8, to = 0.8, ylim=c(-20, 30))
line.tan(px, dx=3, p=y, dp=dy)
# Inflexion
curve(dy(x), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, col="orange")



#########
### Ex 5:
p0 = toPoly.pm("x^3 - 1")
p1 = toPoly.pm("x")
p2 = toPoly.pm("x^2")
pE1 = toPoly.pm("x^2 + 2*x")
pE2 = toPoly.pm("x^2 - 2*x")
pl1 = list(Exp=pE1, Poly=p1)
pl2 = list(Exp=pE2, Poly=p2)
#
solveP.Dexp(pl1, pl2, p0=p0)

### ODE:
x^2*(4*x - 1)*d2y - 2*(8*x^4 - 2*x^3 + 6*x^2 - x)*dy +
	+ 2*(8*x^5 - 2*x^4 + 5*x^2 + 6*x - 1)*y +
	- (16*x^8 - 4*x^7 - 48*x^6 + 6*x^5 + 4*x^4 - 2*x^3 - 10*x^2 - 12*x + 2) # = 0


### Plot:
# pl1, pl2, p0: see above;
d2y = function(x) {
	yx = y(x);
	dyx = dy(x);
	x2 = x^2;
	div = x2*(4*x - 1);
	d2y = 2*(8*x2^2 + 6*x2 - x*(2*x2 + 1))*dyx +
		- 2*(x*(8*x2^2 + 6) - 2*x2^2 + 5*x2 - 1)*yx +
		+ (16*x^8 - 4*x^7 - 48*x^6 + 6*x^5 + 4*x^4 - 2*x^3 - 10*x^2 - 12*x + 2)
	d2y = ifelse(div != 0, d2y/div, 6); # TODO: correct!
	return(d2y)
}
### Plot:
px = c((-3:3) * 2/9);
curve(y(x), from = -0.8, to = 0.8, ylim=c(-20, 30))
line.tan(px, dx=3, p=y, dp=dy)
# Inflexion
curve(dy(x), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, col="orange")

###
px = 1 + (1:5)/13;
curve(y(x), from = 1, to = 1.5, ylim=c(30, 500))
line.tan(px, dx=3, p=y, dp=dy)
# exponential sub-domain
curve(dy(x), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, col="orange")

