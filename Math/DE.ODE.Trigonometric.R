########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Trigonometric
###
### draft v.0.4e


### Non-Linear & Linear:
### Trigonometric Variants


###############
### History ###
###############

### draft v.0.4e: [refactor]
# - moved NL ODEs of type Trig(y) to new file:
#   DE.ODE.NL.TrigY.R;
### draft v.0.4d: [refactor]
# - moved Section on I( Trig ) to new file:
#   DE.ODE.NL.Trig.Int.R;
### draft v.0.4b - v.0.4b-clean4:
# - moved Section on Automatic generation
#   & Basic types to a new file:
#   DE.ODE.Trigonometric.Basic.R;

### draft v.0.4c:
# - concept: sin(exp(y)) - y = F(x);

### Order 1 & 2 Linear:
### Trigonometric Variants

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

### draft v.0.3e:
# - extension of [v.0.3d] to cases derived from:
#   sin(x^n + k*x);
### draft v.0.3d-pre - v.0.3d:
# - ODE: x*d2y - (n-1)*dy + n^2*x^(2*n-1)*y = n*x^n; [fixed]
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


####################

### Helper Functions

# library(pracma)
# needed for Lambert W;

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#########################
#########################

### Trigonometric Functions

### Section A: Simple variants
#   G(y) = P1(x) * sin(T(x)) + F(x)
#   where G(y) = polynomial(y);
# - moved to file:
#   DE.ODE.Trigonometric.Basic.R;


#########################
### Section B: Non-Linear
### Inverse Trigonometric
### [Trig(y)]

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
line.tan(c((-2:2)/2.2), dx=3, p=y, dp=dy)
# pseudo-sigmoidal
curve(dy(x), add=T, col="green")
line.tan(c(-2, -1.5, 1.5, 2)/2.1, dx=1/5, p=dy, dp=d2y, col="orange")

# check full pseudo-sigmoidal
curve(dy(x), from= -1, to = 1, col="green")
line.tan(c(-2, -1.5, 1.5, 2)/2.1, dx=1/5, p=dy, dp=d2y, col="orange")


### b = -1;
# x^3*(x^2 - 2)*y^3*d2y + x^4*y^3*dy - x^3 = 0
# (x^2 - 2)*y^3*d2y + x*y^3*dy = 1
b = -1;
#
curve(y(x, b=b), from= -sqrt(2) + 1E-10, to = sqrt(2) - 1E-10)
# global minimum;
line.tan(c((-2:2)/2.2), dx=3, p=y, dp=dy, b=b)
# pseudo-sigmoidal
curve(dy(x, b=b), add=T, col="green")
line.tan(c(-2, -1.5, 1.5, 2)/2.1, dx=1/5, p=dy, dp=d2y, b=b, col="orange")

# check full pseudo-sigmoidal
curve(dy(x, b=b), from= -sqrt(2) + 1E-10, to = sqrt(2) - 1E-10, col="green", ylim=c(-4,4))
line.tan(c(-2.8, -2.4, -2, -1.5, 1.5, 2, 2.4, 2.8)/2.1, dx=1/5, p=dy, dp=d2y, b=b, col="orange")


### b = -1/2;
b = -1/2
curve(dy(x, b=b), from= -sqrt(5/4), to = sqrt(5/4), col="green", ylim=c(-4,4))
line.tan(c(-2.4, -2, -1.5, 1.5, 2, 2.4)/2.1, dx=1/5, p=dy, dp=d2y, b=b, col="orange")


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
### Plot:

### b = 0;
curve(y(x), from= 0, to= 1, ylim=c(0, 1.5))
# quasi/inverted sigmoidal;
line.tan(c(0, 1/5, 2/3, 0.94), dx=3, p=y, dp=dy)
# global minimum
curve(dy(x), add=T, col="green")
line.tan(c(0:4/5, 0.9), dx=1/5, p=dy, dp=d2y, col="orange")


# check full D(y) curve:
curve(dy(x), from= 0, to = 1, col="green", ylim=c(0, 2))
line.tan(c(0:4/5, 0.9), dx=1/5, p=dy, dp=d2y, col="orange")


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

### ODE:
### Combinations:
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
div = 4.55
px  = c(-1 + (1:4)/div, 1 - (1:4)/div);
#
curve(y, from=-1, to=1)
line.tan(px, dx=1.2, p=y, dp=dy)
### D2(y):
curve(dy, from=-1, to=1, col="green", ylim=c(-1, 8))
curve(y, from=-1, to=1, col="grey", add=T)
line.tan(px, dx=0.5, p=dy, dp=d2y)


### TODO: tan, ln;

########################
########################

#################
### Section C ###
#################

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
line.tan(c((-5:5)/2.2), dx=3, p=y, dp=dy, a=a, n=n)
# wave
curve(dy(x, a=a, n=n), add=T, col="green")
line.tan(c((-5:5)/2.2), dx=1/5, p=dy, dp=d2y, a=a, n=n, col="orange")

### Test separately:
curve(dy(x, a=a, n=n), from= -3, to= 3, col="green")
line.tan(c((-5:5)/2.2), dx=1/5, p=dy, dp=d2y, a=a, n=n, col="orange")


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
line.tan(c((0:5)/2.2), dx=3, p=y, dp=dy, a=a, n=n)
# wave
curve(dy(x, a=a, n=n), add=T, col="green")
line.tan(c((0:5)/2.2), dx=1/5, p=dy, dp=d2y, a=a, n=n, col="orange")


### Test separately:
curve(dy(x, a=a, n=n), from= 0, to= 3, col="green")
line.tan(c((0:5)*2/3.2), dx=1/5, p=dy, dp=d2y, a=a, n=n, col="orange")



############################
############################
############################

##############
### Basics ###
##############

### Linear Simple
### (Polynomial)
###       &
### Linear Complex
### (Non-Polynomial)

# y = P(x)*sin(T(x))
# where T(x) = Polynomial or Non-Polynomial;

# - moved to file:
#   DE.ODE.Trigonometric.Basic.R;


############################
############################

### Note:
# - ODE derived from Integrals:
#   moved to file: DE.ODE.NL.Trig.Int.R;

