########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Fractions: Lambert
###
### draft v.0.4i


###############
### History ###
###############

### Order 1 Non-Linear
###
### draft v.0.4g - v.0.4i:
# - derived from:
#   y^n * exp(y^m) = P1(x)*y + F(x);
#   y^n * exp(G(x)/y^m) = P1(x)*y + F(x); [v.0.4h]
# - e.g. for m = -1 & p1 = f:
#   f*dy + df*y^3 + df*y^2 = 0; [& solution in v.0.4i]
#   p1*g^2*dy + dp1*y^3 + dp1*g*y^2 - p1*g*dg*y = 0; [v.0.4h]
### draft v.0.4f - v.0.4f-ex2:
# - trivial Lambert example:
#   W(y+x) = F(x);
# - more examples for:
#   W(y+x) = y + b*x;
### draft v.0.4e - v.0.4e-various:
# - based on: y^m * exp(y^n) = F(x);
#   e.g.: (m*y^m + n)*dy - (k+1)*x^k * y = 0;
# - more cleanup; [v.0.4e-clean]
### draft v.0.4d: pre & test:
# - generalization:
#   (y + F1(x)) * (W(x+k) + F2(x)) = F0(x);
### draft v.0.4c:
# - derived from:
#   y = F1(x) * W(x+k) + F0(x);
# - examples:
#   (x+k)*f1*y*dy - (x+k)*df1*y^2 - f1^2*y + f1^3 = 0;
### draft v.0.4b:
# - cleanup:
#   moved section with y = log(P1(x))*log(P2(x))
#   to file: DE.ODE.Log.R;
### draft v.0.4a - v.0.4a-ext:
# - derived from:
#   y * W(x) = F0(x);
#   y * W(x + k) = F0(x); [v.0.4a-ext]
# - examples:
#   x*y*dy + b0*x*dy + y^2 = 0;
#   (x + k)*y*dy + b0*(x + k)*dy + y^2 = 0;
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

# library(pracma)
# needed for Lambert W;


# include: Polynomials.Helper.R;
# include: DE.ODE.Helper.R;
source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")

#########################

#####################
### G(y) = W(...) ###
### W = Lambert W ###
#####################

### Theory

### g(y) * e^(g(y)) = f(x)

### D =>
g*e^g*Dg + e^g*Dg - df # = 0 # * g
g^2*e^g*Dg + g*e^g*Dg - df*g # = 0
# g*e^g = f =>
f*g*Dg + f*Dg - df*g # = 0

#############
### Examples:

### Trivial:
# g = y:
f*y*dy + f*dy - df*y # = 0

# g = y^2
2*f*y^3*dy + 2*f*y*dy - df*y^2 # = 0
2*f*y^2*dy + 2*f*dy - df*y # = 0

# g = y^n
n*f*y^(2*n-1)*dy + n*f*y^(n-1)*dy - df*y^n # = 0
n*f*y^n*dy + n*f*dy - df*y # = 0
n*y^n*dy + n*dy - (df/f)*y # = 0


### Non-Trivial:
# g = y + a*x; a = constant;
f*(y + a*x)*(dy + a) + f*(dy + a) - df*(y + a*x) # = 0
f*y*dy + f*(a*x + 1)*dy + (a*f - df)*y + a^2*f*x - a*x*df + a*f # = 0
# for f = x^2 * e^(a*x):
# e^(a*x) can be factored both from f & df, as df = (a*x^2 + 2*x)*e^(a*x):
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
### Plot:
px = c(0, 1, 2, 4.5) / 2;
curve(y(x), from=-1/5, to=3)
# a nice global minimum
line.tan(px, dx=3, p=y, dp=dy)

###
curve(y(x), from=-1/5, to=5)
# a nice global minimum
line.tan(c((0:3)*4/3), dx=3, p=y, dp=dy)

### a == 2
a = 2;
curve(y(x, a=a), from=-1/5, to=3)
# a nice global minimum
line.tan(c(0:4, 7.5, 10)/5, dx=3, p=y, dp=dy, a=a)

### a == 1/2
a = 1/2;
curve(y(x, a=a), from=-1/5, to=3)
# a nice global minimum
line.tan(c(0:4, 10)/5, dx=3, p=y, dp=dy, a=a)


#####################

##################
### Asymmetric ###
##################

### (x + y)*e^y = x^2 + b
# ODE:
(x^2+b)*y*dy + (x^2+b)*(x+1)*dy - 2*x*y - x^2 + b # = 0
(x^n+b)*y*dy + (x^n+b)*(x+1)*dy - n*x^(n-1)*y - (n-1)*x^n + b # = 0

### Solution:
y = function(x, b, n=2) {
	# root
	y = lambertWp((x^n+b) * exp(x)) - x
	y = sapply(y, round0)
	return(y)
}
dy = function(x, b, n=2) {
	y.x = y(x, b, n=n)
	xn = x^n; xnb = xn + b;
	div = xnb*(y.x + x + 1)
	xn1 = if(n == 2) x else x^(n-1);
	dp = n * xn1 * y.x + (n-1)*xn - b;
	dp = ifelse(div != 0, dp / div, -1);
	return(dp)
}
### Plot:
n = 2;
curve(y(x, b=1/2, n=n), from=-1, to=3)
# a nice global minimum
line.tan(c((0:4)/2), dx=3, p=y, dp=dy, b=1/2, n=n)

### Ex 2:
n = 3
curve(y(x, b=1/2, n=n), from=-1, to=3)
# a nice minimum (pseudo-global)
line.tan(c((0:4)/2), dx=3, p=y, dp=dy, b=1/2, n=n)


#####################
### (x + a + y)*e^y = x^2 + b
# [not run]
(1 + dy)*e^y + (x + a + y)*e^y*dy = 2*x
(1 + dy)*(x^2 + b)/(y + x + a) + (x^2 + b)*dy - 2*x # = 0
(1 + dy)*(x^2 + b) + (x^2 + b)*(y + x + a)*dy - 2*x*(y + x + a) # = 0
(x^2 + b)*y*dy + (x + a + 1)*(x^2 + b)*dy - 2*x*y - x^2 - 2*a*x + b # = 0
### ODE:
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
### Plot:
curve(y(x, a=-1, b=1/2), from=-2, to=3)
# a nice global minimum
line.tan(c(-1, (0:4)/3), dx=3, p=y, dp=dy, a=-1, b=1/2)


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
### Plot:
px = c(-0.42, (0:2)/2, 1.8)
curve(y(x, n=3, b=1), from=-1/2, to=3)
# possible inflexion
line.tan(px, dx=3, p=y, dp=dy, n=3, b=1)


### Generalisations

### (y + Ax)*e^(y + Ex) = Fx
# Ax, Ex, Fx = f(x);
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
x*y*dy + x*(x^2 + x)*dy - 2*(x + 1)*y - (2*x^3 + 2*x^2 - x - 2) # = 0
y*dy + (x^2 + x)*dy - 2*(x + 1)/x * y - (2*x^3 + 2*x^2 - x - 2)/x # = 0

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
### Plot:
curve(y(x), from=-1, to=1)
# a nice global minimum
line.tan(c((-4:4)/5), dx=3, p=y, dp=dy)


#######################

#######################
### Different Power ###
#######################

### y^n * exp(y^m) = F(x)
# - the simple cases are trivially solvable;
# - see other sections for variants of type:
#   G1(y, x) * exp(G2(y, x)) = F(x);

### D =>
(m*y^(m+n-1) + n*y^(n-1)) * exp(y^m)*dy - df # = 0 # * y =>
### ODE:
f*(m*y^m + n)*dy - df*y # = 0

### Example:
# m = 1:
f*(y + n)*dy - df*y # = 0
# m = n:
n*f*(y^n + 1)*dy - df*y # = 0
# f = exp(x^(k+1))
(m*y^m + n)*dy - (k+1)*x^k * y # = 0

### Solution & Plot:
y = function(x, m=c(2,2), FUN, DFUN=NULL, useNeg=FALSE) {
	# m = m[1]; n = m[2];
	fr = m[1] / m[2];
	# root
	fx = fr * FUN(x)^fr;
	y = lambertWp(fx);
	if(useNeg) {
		isNeg = (fx < 0)
		if(any(isNeg)) {
			y[isNeg] = lambertWn(fx);
		}
	}
	y = rootn(y, m[1]);
	y = sapply(y, round0)
	return(y)
}
dy = function(x, m=c(2,2), FUN, DFUN, useNeg=FALSE) {
	y.x = y(x, m=m, FUN=FUN, useNeg=useNeg)
	# f*(m*y^m + n)*dy - df*y
	div = FUN(x) * (m[1]*y.x^m[1] + m[2]);
	dp = DFUN(x) * y.x;
	dp = ifelse(div != 0, dp / div, 0); # TODO: check;
	return(dp)
}
### Plot:
b = 3;
FUN  = function(x) x^2 + b[1];
DFUN = function(x) 2*x;
#
curve(y(x, FUN=FUN), from= -1, to=3)
# global minimum
line.tan(c(-1/2, -1/5, 0:4 / 3), dx=3, p=y, dp=dy, FUN=FUN, DFUN=DFUN)


##################
### Polynomial ###
### Terms      ###
##################

### y^n * exp(y^m) = P1(x)*y + F(x)

### D =>
(m*y^(m+n-1) + n*y^(n-1)) * exp(y^m)*dy - p1*dy - dp1*y - df # = 0 # * y =>
(m*y^m + n)*(p1*y + f)*dy - p1*y*dy - dp1*y^2 - df*y # = 0

### Special Cases:
# m = 1
(y + n)*(p1*y + f)*dy - p1*y*dy - dp1*y^2 - df*y # = 0
(p1*y^2 + ((n-1)*p1+f)*y + n*f)*dy - dp1*y^2 - df*y # = 0
# p1 = f = exp(k*x)
(m*y^m + n)*(y + 1)*dy - y*dy - k*y^2 - k*y # = 0
(m*y^(m+1) + m*y^m + (n-1)*y + n)*dy - k*y^2 - k*y # = 0
# & m = -1:
((n-1)*y^2 + (n-1)*y - 1)*dy - k*y^3 - k*y^2 # = 0
# & m = -1, n = 1:
dy + k*y^3 + k*y^2 # = 0

### Special Case: [check]
### y * exp(1/y) = p1*y + f

### D =>
(1 - 1/y)*exp(1/y)*dy - p1*dy - dp1*y - df # = 0 # * y^2
(y - 1)*(p1*y + f)*dy - p1*y^2*dy - dp1*y^3 - df*y^2 # = 0
((f-p1)*y - f)*dy - dp1*y^3 - df*y^2 # = 0
### Case p1 = f
f*dy + df*y^3 + df*y^2 # = 0


### Solution:

### Step 1:
# A*dy + B*y^3 + B*y^2 = 0
# D(A*h) = B*h =>
A*dh + (dA - B)*h # = 0
# log(h) = I( (B-dA) / A ) dx;
### Step 1: [short]
# h = exp(I(B/A)); # where p1 = h!

### Step 2:
# y * exp(1/y) = h*y + h
# y = 1/z =>
# exp(z) = h*z + h # EXP() =>
# exp(exp(z)) = exp(h*z) * exp(h) # POW(1/h) =>
# exp(exp(z) / h) = exp(z) * exp(1) =>
# - exp(z)/h * exp( - exp(z)/h) = - exp(-1)/h;
# - exp(z)/h = W(- exp(-1)/h);
# exp(z) = - h * W(- exp(-1)/h);
z = log(- h * W(- exp(-1)/h));

### Solution & Plot:
eval.H = function(x, FA, FB, lower=0) {
	f = if(is.null(FA)) {
		function(x) eval.FUN(x, FB);
	} else {
		function(x) eval.FUN(x, FB) / eval.FUN(x, FA);
	}
	# H(x)
	h = I.FUN(x, F=f, lower=lower);
	h = exp(h);
	return(h);
}
y = function(x, FB, FA=NULL, lower=0, valH=FALSE) {
	# H(x)
	h = eval.H(x, FA=FA, FB=FB, lower=lower)
	#
	z = lambertWn(- exp(-1)/h);
	z = log(-h * z);
	y = 1 / z;
	y = round0(y);
	if(valH) attr(y, "H") = h;
	return(y)
}
dy = function(x, FB, FA=NULL, lower=0) {
	y.x = y(x, FA=FA, FB=FB, lower=lower);
	dfx = eval.FUN(x, FB);
	fx  = if(is.null(FA)) 1 else eval.FUN(x, FA);
	# f*dy + df*y^3 + df*y^2 = 0
	dp  = dfx*y.x^3 + dfx*y.x^2;
	div = - fx;
	dp = ifelse(div != 0, dp / div, 0); # TODO: check;
	return(dp)
}
### Plot:
FA  = toPoly.pm("x^3 - 2*x + 5")
FB = toPoly.pm("x^2 + 5*x + 1")
px = 1 + (1:6) * 1/3
#
curve(y(x, FB=FB, FA=FA), from = 1, to = 3)
line.tan(px, dx=3, p=y, dp=dy, FA=FA, FB=FB)


### Ex 2:
FA  = toPoly.pm("5 - x^-1")
FB = toPoly.pm("x^2 + 5*x + 1")
px = 1 + (1:6) * 1/3
#
curve(y(x, FB=FB, FA=FA, lower=1), from = 1, to = 3)
line.tan(px, dx=3, p=y, dp=dy, FA=FA, FB=FB, lower=1)


##########
### Slight Generalization
### y * exp(G(x)/y) = P1(x)*y + F(x)

### D =>
((1 - g/y)*dy + dg)*exp(g/y) - p1*dy - dp1*y - df # = 0 # * y^2 =>
((y - g)*dy + dg*y)*(p1*y + f) - p1*y^2*dy - dp1*y^3 - df*y^2 # = 0
(y - g)*(p1*y + f)*dy + dg*(p1*y + f)*y - p1*y^2*dy - dp1*y^3 - df*y^2 # = 0
(f*y - p1*g*y - f*g)*dy - dp1*y^3 - df*y^2 + dg*(p1*y + f)*y # = 0

### Special Cases:
# f = p1*g =>
p1*g^2*dy + dp1*y^3 + dp1*g*y^2 - p1*g*dg*y # = 0
# p1 = x^2;
x*g^2*dy + 2*y^3 + 2*g*y^2 - x*g*dg*y # = 0


################

################
### Radicals ###
################

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
line.tan(c(-1/3, -1/5, (0:4)/3), dx=3, p=y, dp=dy)


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
### Plot:
curve(y(x), from= -0.5, to=3)
line.tan(c(-1/3, -1/5, (0:4)/3), dx=3, p=y, dp=dy)


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
### Plot:
c = 1;
curve(y(x, b=c), from= 1.9, to=4)
line.tan(c(1.9, (4:6)/2), dx=3, p=y, dp=dy, b=c)


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
### Plot:
n = 4L;
curve(y(x, n=n), from= -2, to=2)
line.tan(c((-5:5)/2.8), dx=3, p=y, dp=dy, n=n)


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
### ODE:
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
### Plot:
a = -1; b = 1;
curve(y(x, a=a, b=b), from=-2, to=3)
# NO minimum: only inflexion;
line.tan(c(-1, (0:3)), dx=3, p=y, dp=dy, a=a, b=b)

### y*dy + x*dy + 2*tan(x)*y + 2*x*tan(x) - 2*tan(x) + 1 = 0
a = -1; b = -1;
curve(y(x, a=a, b=b), from=-3, to=2)
# a "virtual" minimum
line.tan(c(-1, 1+(0:3)/3), dx=3, p=y, dp=dy, a=a, b=b)

###
a = -1; b = -1/4;
curve(y(x, a=a, b=b), from=-3, to=2)
# NO minimum: only inflexion;
line.tan(c(-1, (0:3)/2), dx=3, p=y, dp=dy, a=a, b=b)


#########################

###################
### Logarithmic ###
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
line.tan(c((-2:4)/2), dx=3, p=y, dp=dy)

### a == -1
a = -1
curve(y(x, a=a), from=-3, to=3)
# a nice global minimum
line.tan(c((-2:4)/2), dx=3, p=y, dp=dy, a=a)

### a == 1/2
a = 1/2
curve(y(x, a=a), from=-3, to=3)
# a nice global minimum
line.tan(c((-2:4)/2), dx=3, p=y, dp=dy, a=a)


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
line.tan(c((-3:2)/2, 2), dx=3, p=y2.f, dp=dy2, a=a)
# TODO: compute exact value for -1.473744;

# separately: only the integration: I(y^2)
curve(y2.I(x, a=a), from=-2, to=3)
line.tan(c((-3:2)/2, 2), dx=3, p=y2.I, dp=dy2, a=a)

###
a = 3/2
curve(y2.f(x, a=a), from=-2, to=3)
# + 0.1 to separate curves;
curve(y2.I(x, a=a, lower=-0.98 + 0), add=T, col="green")
line.tan(c((-3:2)/2, 2), dx=3, p=y2.f, dp=dy2, a=a)


### Example 2:
### f = e^(-a*x^2); df = -2*a*x*f;
# complete f = e^(-a*x^2 + a*x);
y*dy + (a*x + 1)*dy + 2*a*x*y + 2*a^2*x^2 + a # = 0


### Solution & Plot:
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
line.tan(c((-3:3)/1.5, -0.50404905), dx=3, p=y, dp=dy)

### a == 1/2
a = 1/2
curve(y(x, a=a), from=-3, to=3)
# quasi-bi-sigmoidal
line.tan(c((-3:3)/1.5), dx=3, p=y, dp=dy, a=a)

### a == 1/3
a = 1/3
curve(y(x, a=a), from=-3, to=3)
# quasi-bi-sigmoidal
line.tan(c((-3:3)/1.5), dx=3, p=y, dp=dy, a=a)

### a == -1/3
a = -1/3
curve(y(x, a=a), from=-3, to=3)
# global minimum
line.tan(c((-3:3)/1.5), dx=3, p=y, dp=dy, a=a)

### TODO: Integration by parts;


###########################
###########################

### e^(y + a1*x) * e^(e^(y + a1*x)) = x^n

# TODO: generalize F0(x);

### D =>
(dy + a1)*(1 + e^(y + a1*x)) - n/x # = 0
### D2 =>
n*x*d2y - x^2*(dy + a1)^3 + n*x*(dy + a1)^2 + n*(dy + a1) # = 0;
# w = dy:
n*x*dw - x^2*(w + a1)^3 + n*x*(w + a1)^2 + n*(w + a1) # = 0;


### Solution & Plot
y = function(x, a=1, n=2) {
	# root
	y = log(lambertWp(x^n)) - a[1]*x
	y = round0(y);
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
# global "minimum" / asymmetric horn;
line.tan(c((-3:3)*2/5), dx=3, p=y, dp=dy, a=a, n=n)
# +/- Inf discontinuity
curve(dy(x, a=a, n=n), add=T, col="green")
line.tan(c((-4:4)/5), dx=3, p=dy, dp=d2y, a=a, n=n, col="orange")


### Example 2:
a = 1/3; n = 2;
curve(y(x, a=a, n=n), from= -2, to = 2, ylim=c(-7, 3))
# global "minimum" / horn;
line.tan(c((-3:3)*2/5), dx=3, p=y, dp=dy, a=a, n=n)
#
curve(dy(x, a=a, n=n), add=T, col="green")
line.tan(c((-4:4)/5), dx=3, p=dy, dp=d2y, a=a, n=n, col="orange")


#####################
#####################

#################
### Lambert W ###
#################

### y * W(x) = F0(x)

### D(y)
W*dy + W/(W*x + x) * y - df0 # = 0
y*W*dy + W/(x*W + x) * y^2 - df0*y # = 0
f0*dy + f0/(x*y*W + x*y) * y^2 - df0*y # = 0
f0*dy + f0/(x*f0 + x*y) * y^2 - df0*y # = 0
f0*(x*f0 + x*y)*dy + f0*y^2 - df0*(x*f0 + x*y)*y # = 0
### ODE:
x*f0*y*dy + x*f0^2*dy - x*df0*y^2 + f0*y^2 - x*f0*df0*y # = 0

### Extensions:
### y * W(x + k) = F0(x)

### D(y)
W*dy + W/((x + k)*(W + 1)) * y - df0 # = 0
f0*(x+k)*(f0 + y)*dy + f0*y^2 - df0*(x+k)*(f0 + y)*y # = 0


### Examples:
### f0 = b0
x*y*dy + b0*x*dy + y^2 # = 0
### Ext:
(x + k)*y*dy + b0*(x + k)*dy + y^2 # = 0


### Solution & Plot
y = function(x, b=0, k=0, n=1, pos.br=TRUE) {
	W.x = if(pos.br) lambertWp(x^n + k) else lambertWn(x^n + k);
	f0 = eval.pol(x, b);
	div = W.x;
	r = ifelse(div != 0, f0/div, Inf) # TODO
	return(r)
}
dy = function(x, b=0, k=0, n=1, pos.br=TRUE) {
	y.x = y(x, b=b, k=k, n=n, pos.br=pos.br);
	# x*f0*(f0 + y)*dy + f0*y^2 - df0*x*(f0 + y)*y
	f0 = eval.pol(x, b); df0 = deriv.pol(x, b, dn=1);
	xf0 = (x + k)*(f0 + y.x);
	dp = - f0*y.x^2 + df0*xf0*y.x;
	div = xf0 * f0;
	dp = ifelse(div != 0, dp/div, 1) # TODO
	return(dp)
}
### Plot:
b = c(1, 1)
lim = -exp(-1)
px = c(-(4:1) * 1/13)
curve(y(x, b=b), from= lim[1], to = 0, ylim=c(-15, 2))
#
line.tan(px, dx=3, p=y, dp=dy, b=b)


### Ex 2:
b = c(1, 1)
lim = -exp(-1)
px = - c(9,7, 5, 2,1) * 1/23 + 0.03
curve(y(x, b=b, pos.br=FALSE), from= lim[1], to = 0, ylim=c(-0.8, 0))
#
line.tan(px, dx=c(0.32, 0.37, 0.37, 1/4, 0.2), p=y, dp=dy, b=b, pos.br=FALSE)


### Ex 3:
b = c(1, 1); k = 2;
px = c(-(4:1) * 6/13) + 0.1
curve(y(x, b=b, k=k), from= -k, to = -k + 2, ylim=c(-5, 2))
#
line.tan(px, dx=3, p=y, dp=dy, b=b, k=k)


#######################
#######################

### y = F1(x) * W(x+k) + F0(x)

### D(y)
df1*W + f1*W / ((x+k)*(W+1)) + df0
### f1*dy =
df1*(y - f0) + f1*(y - f0) / ((x+k)*(W+1)) + df0*f1
### (x+k)*f1*dy =
(x+k)*df1*(y - f0) + f1^2*(y - f0) / (f1*W+f1) + (x+k)*df0*f1
(x+k)*df1*(y - f0) + f1^2*(y - f0) / (y - f0 + f1) + (x+k)*df0*f1
(x+k)*df1*y + f1^2*(y - f0) / (y - f0 + f1) + (x+k)*(df0*f1 - f0*df1)

### ODE:
### (x+k)*f1*(y - f0 + f1)*dy =
(x+k)*df1*(y - f0 + f1)*y + f1^2*(y - f0) + (x+k)*(df0*f1 - f0*df1)*(y - f0 + f1)

### Examples:
### f0 = f1
(x+k)*f1*y*dy - (x+k)*df1*y^2 - f1^2*y + f1^3 # = 0


### Solution & Plot
y = function(x, b0=1, b1=1, k=0, n=1, pos.br=TRUE) {
	W.x = if(pos.br) lambertWp(x^n + k) else lambertWn(x^n + k);
	f0 = eval.pol(x, b0);
	f1 = eval.pol(x, b1);
	r = f1 * W.x + f0;
	return(r)
}
dy = function(x, b0=1, b1=1, k=0, n=1, pos.br=TRUE) {
	y.x = y(x, b0=b0, b1=b1, k=k, n=n, pos.br=pos.br);
	f0 = eval.pol(x, b0); df0 = deriv.pol(x, b0, dn=1);
	f1 = eval.pol(x, b1); df1 = deriv.pol(x, b1, dn=1);
	# (x+k)*df1*(y - f0 + f1)*y + f1^2*(y - f0) + (x+k)*(df0*f1 - f0*df1)*(y - f0 + f1)
	fs = (y.x - f0 + f1); xf = (x + k)*fs;
	dp = xf*df1*y.x + f1^2*fs + xf*(df0*f1 - f0*df1) - f1^3;
	div = xf*f1;
	dp = ifelse(div != 0, dp/div,
		dy(x + 1E-3, b0=b0, b1=b1, k=k, n=n, pos.br=pos.br)) # TODO
	return(dp)
}
### Plot:
b0 = c(1, 1); b1 = c(1, 1)
k = 0;
lim = -exp(-1) - k;
px = c(-4,-3, -0.3) * 1/12
curve(y(x, b0=b0, b1=b1, k=k), from= lim[1], to = 0, ylim=c(-0.5, 1))
#
line.tan(px, dx=3, p=y, dp=dy, b0=b0, b1=b1, k=k)


### Ex 2:
b0 = c(1, 1); b1 = c(1, 1)
k = 2;
lim = -exp(-1) - k;
px = c(-7,-6, -5, -2) * 1/3
curve(y(x, b0=b0, b1=b1, k=k), from= lim[1], to = 0, ylim=c(-3, 1))
#
line.tan(px, dx=3, p=y, dp=dy, b0=b0, b1=b1, k=k)


### Ex 3:
b0 = c(1, 1); b1 = c(1, -2, -3)
k = 2;
lim = -exp(-1) - k;
px = c(-10,-9, -8, -6, -4, -3) * 1/5
curve(y(x, b0=b0, b1=b1, k=k), from= lim[1], to = 0, ylim=c(-3, 1))
# slight concavity between -1.3 to 0;
line.tan(px, dx=1.4, p=y, dp=dy, b0=b0, b1=b1, k=k)


########################
########################

### (y + F1(x)) * (W(x+k) + F2(x)) = F0(x)

### D(y)
(dy + df1)*(W + f2) + (y + f1)*(W/((x+k)*(W+1)) + df2) - df0 # = 0
(dy + df1)*f0 / (y + f1) + (y + f1)*(W/((x+k)*(W+1)) + df2) - df0 # = 0
(dy + df1)*f0 + (y + f1)^2*(W/((x+k)*(W+1)) + df2) - df0*(y + f1) # = 0
(x+k)*f0*(dy + df1) + (y + f1)^2*(W/(W+1) + (x+k)*df2) - (x+k)*df0*(y + f1) # = 0
(x+k)*f0*(dy + df1) + (y + f1)^2*W/(W+1) + (x+k)*df2*(y + f1)^2 - (x+k)*df0*(y + f1) # = 0
(x+k)*f0*(W+1)*(dy + df1) + (f0 - f2*(y + f1))*(y + f1) +
	+ (x+k)*df2*(W+1)*(y + f1)^2 - (x+k)*df0*(W+1)*(y + f1) # = 0
(x+k)*f0*(W+1)*(dy + df1) + (f0 - f2*(y + f1))*(y + f1) +
	+ (x+k)*df2*(f0 - f2*(y + f1) + y + f1)*(y + f1) - (x+k)*df0*(f0 - f2*(y + f1) + y + f1) # = 0
(x+k)*f0*(W+1)*(dy + df1) + (f0 - f2*(y + f1))*(y + f1) +
	+ (x+k)*(df2*(y + f1) - df0)*(f0 - f2*(y + f1) + y + f1) # = 0
(x+k)*f0*(f0 - f2*(y + f1) + y + f1)*(dy + df1) + (f0 - f2*(y + f1))*(y + f1)^2 +
	+ (x+k)*(df2*(y + f1) - df0)*(f0 - f2*(y + f1) + y + f1)*(y + f1) # = 0
# T = f0 - f2*(y + f1) + y + f1;
(x+k)*f0*T*(dy + df1) + (T - (y + f1))*(y + f1)^2 +
	+ (x+k)*(df2*(y + f1) - df0)*T*(y + f1) # = 0
(x+k)*f0*T*dy - (y + f1)^3 + (x+k)*df2*T*(y + f1)^2 + T*(y + f1)^2 +
	- (x+k)*df0*T*(y + f1) + (x+k)*f0*df1*T # = 0

### ODE:
T = f0 - (f2 - 1)*(y + f1);
(x+k)*f0*T*dy - f2*(y + f1)^3 + (x+k)*df2*T*(y + f1)^2 + f0*(y + f1)^2 +
	- (x+k)*df0*T*(y + f1) + (x+k)*f0*df1*T # = 0


### Solution & Plot
y = function(x, b, k=0, n=1, pos.br=TRUE, all=FALSE) {
	W.x = if(pos.br) lambertWp(x^n + k) else lambertWn(x^n + k);
	f0 = eval.pol(x, b[[1]]);
	f1 = eval.pol(x, b[[2]]);
	f2 = eval.pol(x, b[[3]]);
	div = f2 + W.x;
	# (y + F1(x)) * (W(x+k) + F2(x)) = F0(x)
	r = ifelse(div != 0, f0/div,
		y(x - 1E-3, b=b, k=k, n=n, pos.br=pos.br)); # TODO
	r = r - f1;
	if(all) return(list(y=r, f0=f0, f1=f1, f2=f2));
	return(r)
}
dy = function(x, b, k=0, n=1, pos.br=TRUE) {
	y.all = y(x, b=b, k=k, n=n, pos.br=pos.br, all=T);
	f0 = y.all$f0; f1 = y.all$f1; f2 = y.all$f2;
	y.x = y.all$y;
	df.all = lapply(b, deriv.pol, x=x);
	df0 = df.all[[1]]; df1 = df.all[[2]]; df2 = df.all[[3]];
	#
	yf = y.x + f1;
	T = f0 - (f2 - 1)*yf; xk = x + k; Tk = xk*T;
	div = f0*Tk;
	dp = f2*yf^3 - df2*Tk*yf^2 - f0*yf^2 + df0*Tk*yf - div*df1;
	dp = ifelse(div != 0, dp/div,
		dy(x + 1E-3, b=b, k=k, n=n, pos.br=pos.br)) # TODO
	return(dp)
}
### Plot:
b = list(b0 = c(1, 1), b1 = c(-1, 2, 3), b2 = c(0, 1, -1));
k = 0;
lim = -exp(-1) - k;
px = c(-4,-3, -2, -1, -0.75) * 1/12
curve(y(x, b=b, k=k), from= lim[1], to = 0, ylim=c(-15, 1))
#
line.tan(px, dx=3, p=y, dp=dy, b=b, k=k)


### Ex 2:
b = list(b0 = c(1, 1), b1 = c(-1, 2, 3), b2 = c(0, 1, -1));
k = 0;
px = c(1:5) * 3/11
curve(y(x, b=b, k=k), from = 1E-3, to = 1.5, ylim=c(-1, 10))
#
line.tan(px, dx=3, p=y, dp=dy, b=b, k=k)


### Ex 3:
b = list(b0 = c(1, 1), b1 = c(-1, 2, 3), b2 = c(0, 1, -1));
k = 2;
px = c(1,2,4,5, 5.3) * 3/11
curve(y(x, b=b, k=k), from = 1E-3, to = 1.5, ylim=c(-5, 2.5))
#
line.tan(px, dx=c(1,1,1.5, 2,2), p=y, dp=dy, b=b, k=k)
curve(dy(x, b=b, k=k), add=T, col="green")


#####################
#####################

### W(y+x) = F0(x)
# - trivial example;

### D =>
dW = W*(dy + 1) / ((W + 1)*(y + x)); # !! dW(y+x) !!
dW = f0 * (dy + 1) / ((f0 + 1)*(y + x))
# =>
f0 * (dy + 1) - df0*(f0 + 1)*(y + x) # = 0
### ODE:
f0*dy - df0*(f0 + 1)*y + f0 - x*df0*(f0 + 1) # = 0

### Check:
# f0 = x - 1
(x-1)*dy - x*y - x^2 + x - 1 # = 0
# (x-1)*exp(x-1) = y + x
# D => (x-1)*dy - x*y - x^2 + x - 1 = 0;


####################

### W(y+x) = y + b*x

### D =>
dW = W*(dy + 1) / ((W + 1)*(y + x)); # !! dW(y+x) !!
(y + b*x)*(dy + 1) - (y + b*x + 1)*(y + x)*(dy + b) # = 0
(y + b*x)*dy - (y^2 + b*x*y + y + x*y + b*x^2 + x)*dy +
	- b*(y^2 + b*x*y + y + x*y + b*x^2 + x) + (y + b*x) # = 0
(y^2 + (b+1)*x*y + b*x^2 - (b-1)*x)*dy +
	+ b*y^2 + b*(b+1)*x*y + (b-1)*y + b^2*x^2 # = 0

### Special Cases:
# b = 0
(y^2 + x*y + x)*dy - y # = 0
# b = -1
(y^2 - x^2 + 2*x)*dy - y^2 - 2*y + x^2 # = 0


### Check
# b = 1 =>
(y + x)^2*(dy + 1) # = 0 # OK

