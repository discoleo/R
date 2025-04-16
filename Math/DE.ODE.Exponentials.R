########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Exponentials
###
### draft v.0.2f


### ODEs Derived from Exponentials

### TODO:
# - move relevant sections from
#   DE.ODE.Fractions.Lambert.R
#   & DE.ODE.Gaussian.R to this file;


###############
### History ###
###############


### draft v.0.2f:
# - [Refactor] moved Trig(EXP) to new file:
#   DE.ODE.NL.Trig.Exp.R;
# - [Refactor] moved log(y + g1) * log(y + g2)
#   to new file: DE.ODE.NL.LogY.R;
### draft v.0.2d - v.0.2e:
# - Exponential extensions:
#   exp(k*x) * y^n = n*p1*y + n*f0;
#   e^(e^y) = y + P(x);
### draft v.0.2c:
# - from I( exp(exp(x)) ):
#   k*y*d2y - k*dy^2 + y^2*dy - k*y*dy - y^3  = 0;
### draft v.0.2b - v.0.2b-check:
# - from Mixed Exp-Trig:
#   y*d2y - dy^2 - y*dy - y^2 + 1 = 0;
#   y*d2y - dy^2 - y*dy - y^2 + k^2 = 0; [v.0.2b-gen & v.0.2b-check]
### draft v.0.2a:
# - moved Section with ODEs based on another ODE
#   to file: DE.ODE.FromODEs.R;
### v.0.1f - v.0.1f-PowFr:
# - derived from ODEs of Higher Power:
#   dy - G(x)*y^2 = F(x);
# - slight generalization:
#   dy - G2(x)*y^2 - G1(x)*y = F(x); [v.0.1f-Gen & v.0.1f-PowFr]
### v.0.1e - v.0.1e-var + v.0.1f-P1Var4:
# - derived from: dy - G(x)*y = F(x)
#   y*d2y - dy^2 + f*dy - dg*y^2 - df*y = 0;
# - more variants; [v.0.1e-var & v.0.1f-P1Var4]
### v.0.1d:
# - [started] derived from another ODE;


####################

### Helper Functions

# library(pracma)
# - may be needed to solve various equations;
#   (e.g. Lambert W)


# include: DE.ODE.Helper.R;
# include: Polynomials.Helper.R;
source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")

### Other
# moved to DE.ODE.Helper.R;
eval.FUN = function(x, F, ...) {
	xl = c(...);
	r = sapply(x, function(x) eval.pm(F, c(x, xl)));
	return(r);
}
eval.DFUN = function(x, F, ...) {
	xl = c(...);
	r = sapply(x, function(x) eval.pm(F, c(x, xl)));
	DF = dp.pm(F, xn="x");
	nms = c("x", names(xl), "coeff");
	DF  = DF[, nms];
	r1 = sapply(x, function(x) eval.pm(DF, c(x, xl)));
	return(data.frame(f=r, df=r1));
}
eval.d2y = function(F, vals) {
	sapply(seq(nrow(vals)), function(nr) eval.pm(F, vals[nr,]));
}
I.FUN = function(x, F, lower=0, ...) {
	sapply(x, function(up) integrate(F, lower=lower, upper=up, ...)$value);
}

### Roots
filter.root = function(r, n=1, tol=1E-8) {
	# only real roots
	r = round0(r, tol=tol);
	r = r[Im(r) == 0];
	r = head(r, n);
	return(r);
}

#########################
#########################


### Section B: Simple
### Non-Linear ODEs

##############
### e^y = P(x)

### D =>
e^y * dy = dp
### D2 =>
e^y * d2y + e^y * dy^2 - d2p # = 0

### ODE:
p*d2y + p*dy^2 - d2p # = 0
# - can be converted trivially to order 1 (but still non-liniar);

### Examples:
### P(x) = x^2
x^2*d2y + x^2*dy^2 - 2 # = 0
### P(x) = sin(a*x)
d2y + dy^2 + a^2 # = 0
### P(x) = cos(x)^3
d2y + dy^2 - 6*tan(x)^2 + 3 # = 0


##################

### Section C: High-Power
### Non-Linear ODEs

##################
### e^(y^2) = P(x)

### D =>
2*e^(y^2) * y*dy = dp
### D2 =>
2*e^(y^2) * (y*d2y + dy^2) + 4*e^(y^2) * y^2*dy^2 - d2p # = 0
2*p*y*d2y + 4*p*y^2*dy^2 + 2*p*dy^2 - d2p # = 0

### ODE:
2*p*y*d2y + 4*p*y^2*dy^2 + 2*p*dy^2 - d2p # = 0

### Examples:
### P(x) = x
y*d2y + 2*y^2*dy^2 + dy^2 # = 0
### P(x) = x^2
x^2*y*d2y + 2*x^2*y^2*dy^2 + x^2*dy^2 - 1 # = 0
### P(x) = sqrt(x^3)
8*x^2*y*d2y + 16*x^2*y^2*dy^2 + 8*x^2*dy^2 - 3 # = 0
### P(x) = sqrt(x^5)
8*x^2*y*d2y + 16*x^2*y^2*dy^2 + 8*x^2*dy^2 - 15 # = 0
### P(x) = sin(a*x)
2*y*d2y + 4*y^2*dy^2 + 2*dy^2 + a^2 # = 0
### P(x) = sin(sqrt(3)/2 * x)
8*y*d2y + 16*y^2*dy^2 + 8*dy^2 + 3 # = 0


### Solution & Plot:
y.gen = function(n) {
	n1=n;
	y1.f = function(x, n=n1) { x^n }
	y1.d1.df = function(x, n=n1) { n*x^(n-1) }
	y1.d2.df = function(x, n=n1) { if(n == 2) 2  else n*(n-1)*x^(n-2) }
	list(y1.f, y1.d1.df, y1.d2.df)
}
y1.lst = y.gen(2)
y = function(x, PFUN, posRoot=TRUE) {
	val = sqrt(log(PFUN[[1]](x)));
	if( ! posRoot) val = - val;
	return(val)
}
dy = function(x, PFUN, posRoot=TRUE) {
	y.x = y(x, PFUN=PFUN, posRoot=posRoot)
	dp = PFUN[[2]](x) / y.x / exp(y.x^2) / 2;
	return(dp)
}
d2y = function(x, PFUN, posRoot=TRUE) {
	p.v = PFUN[[1]](x);
	y.x  = y(x, PFUN=PFUN, posRoot=posRoot);
	dy.x = dy(x, PFUN=PFUN, posRoot=posRoot);
	dp = PFUN[[3]](x) - 4*p.v*y.x^2*dy.x^2 - 2*p.v*dy.x^2
	div = 2 * p.v * y.x;
	dp = ifelse(div != 0, dp/ div, 1) # TODO
	return(dp)
}
### Plot:
px = c(3/5 + (1:3)*3/5);
curve(y(x, PFUN=y1.lst), from = 1, to = 3)
line.tan(px, dx=3, p=y, dp=dy, PFUN=y1.lst)
# inverse-exp-like:
curve(dy(x, PFUN=y1.lst), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, PFUN=y1.lst, col="orange")


### Ex 2:
n = 3/2;
px = c(3/5 + (1:3)*3/5);
y1.lst = y.gen(n=n)
curve(y(x, PFUN=y1.lst), from= 1, to = 3)
line.tan(px, dx=3, p=y, dp=dy, PFUN=y1.lst)
# inverse-exp-like:
curve(dy(x, PFUN=y1.lst), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, PFUN=y1.lst, col="orange")


#####################

### Section D: Higher-Exponentials
### Non-Linear ODEs

######################
### e^(e^y) - y = P(x)

### D =>
exp(y)*exp(exp(y))*dy - dy - dp # = 0
# =>
exp(y)*(y + p)*dy - dy - dp # = 0

# D2 =>
# TODO


########################
########################

### Mixed Exp(x) & y^n

### exp(k*x)*y^n = n*p1*y + n*f0

### D =>
(k*y^n + n*y^(n-1)*dy)*exp(k*x) - n*p1*dy - n*dp1*y - n*df0 # = 0 # * y =>
(k*y + n*dy)*(p1*y + f0) - p1*y*dy - dp1*y^2 - df0*y # = 0

### ODE:
(n-1)*p1*y*dy + n*f0*dy + (k*p1 - dp1)*y^2 + (k*f0 - df0)*y # = 0

### Solution & Plot:
y = function(x, n=3, k=1, FF0, FP1) {
	f0 = if(inherits(FF0, "pm")) eval.FUN(x, FF0) else FF0;
	p1 = if(inherits(FP1, "pm")) eval.FUN(x, FP1) else FP1;
	xe = exp(k*x);
	if(n == 1) {
		yx = f0 / (xe - p1);
	} else {
		coeff0 = rep(0, n-2);
		yx = sapply(seq_along(x), function(id) {
			coeff = c(xe[id], coeff0, - n*p1[id], - n*f0[id]);
			yx = roots(coeff);
			yx = filter.root(yx, n=1);
			return(yx);
		})
	}
	return(yx);
}
dy = function(x, n=3, k=1, FF0, FP1) {
	f.l = eval.DFUN(x, FF0);
	p.l = eval.DFUN(x, FP1);
	yx  = y(x, n=n, k=k, FF0=f.l$f, FP1=p.l$f);
	#
	dp  = (k*p.l$f - p.l$df)*yx^2 + (k*f.l$f - f.l$df)*yx;
	div = (n-1)*p.l$f*yx + n*f.l$f;
	dp  = ifelse(div != 0, - dp/div, 0); # TODO: check;
	return(dp)
}
### Plot:
F0 = toPoly.pm("x^2 + 3*x + 3");
P1 = toPoly.pm("x^3 - 3*x + 5");
n = 3; k = 2;
px = c(-2.25, -2.05, (0:4)*3/5 - 1.5);
curve(y(x, n=n, k=k, FF0=F0, FP1=P1), from = -3, to = 3, n=256)
line.tan(px, dx=2, p=y, dp=dy, n=n, k=k, FF0=F0, FP1=P1)


########################
########################

########################
### Section D:
### Non-Linear ODEs
### Product-Type
########################

######################
### Integral-Based ###
######################

### y * I( exp(exp(x)) ) = exp(exp(x)) * F(x)

### D =>
I*dy + exp(exp(x))*y - (df + f*exp(x))*exp(exp(x)) # = 0
exp(exp(x))*f*dy + exp(exp(x))*y^2 - (df + f*exp(x))*exp(exp(x))*y # = 0
f*dy + y^2 - (df + f*exp(x))*y # = 0
# f*exp(x)*y = f*dy + y^2 - df*y;

### D2 =>
f*d2y + df*dy + 2*y*dy - (df + f*exp(x))*dy - (d2f + (f+df)*exp(x))*y # = 0
f*d2y + 2*y*dy - f*exp(x)*dy - d2f*y - (f+df)*exp(x)*y # = 0 # * f*y =>
f^2*y*d2y + 2*f*y^2*dy - f*(f*dy + y^2 - df*y)*dy - f*d2f*y^2 - (f+df)*(f*dy + y^2 - df*y)*y # = 0
### ODE:
f^2*y*d2y - f^2*dy^2 + f*y^2*dy - f^2*y*dy - (f+df)*y^3 - f*d2f*y^2 + (f+df)*df*y^2 # = 0

### Special Cases:
# f = k;
k*y*d2y - k*dy^2 + y^2*dy - k*y*dy - y^3 # = 0
# f = x;
x^2*y*d2y - x^2*dy^2 + x*y^2*dy - x^2*y*dy - (x+1)*y^3 + (x+1)*y^2 # = 0

### Solution & Plot:
y = function(x, PFUN, n=1, lower=0) {
	xe = exp(x^n); xee = exp(xe);
	fx = eval.FUN(x, PFUN);
	val = fx * xee;
	I = I.FUN(x, function(x) exp(exp(x^n)), lower=lower);
	return(val / I);
}
dy = function(x, PFUN, n=1, lower=0) {
	yx = y(x, PFUN=PFUN, n=n, lower=lower);
	# f*dy + y^2 - (df + f*exp(x))*y = 0;
	xe = exp(x^n);
	fx = eval.FUN(x, PFUN);
	df  = dp.pm(PFUN, xn="x");
	dfx = eval.FUN(x, df);
	dp  = - yx^2 + (dfx + fx*xe)*yx;
	div = fx;
	dp  = ifelse(div != 0, dp/div, 1) # TODO;
	return(dp)
}
d2y = function(x, PFUN, n=1, lower=0) {
	# F(x):
	fx  = eval.FUN(x, PFUN);
	dp  = dp.pm(PFUN, xn="x");
	dfx = eval.FUN(x, dp);
	d2p = dp.pm(dp, xn="x");
	d2fx = eval.FUN(x, d2p);
	# y & dy:
	yx  = y(x, PFUN=PFUN, n=n, lower=lower);
	dyx = dy(x, PFUN=PFUN, n=n, lower=lower);
	F = toPoly.pm("- f^2*dy^2 + f*y^2*dy - f^2*y*dy +
		- f*y^3 - df*y^3 - f*d2f*y^2 + f*df*y^2 + df^2*y^2");
	F = F[, c("dy", "y", "f", "df", "d2f", "coeff")];
	vals = data.frame("dy"=dyx, "y"=yx, "f"=fx, "df"=dfx, "d2f"=d2fx);
	d2y = eval.d2y(F, vals);
	div = - fx^2*yx;
	d2y = ifelse(div != 0, d2y / div, 1) # TODO
	return(d2y)
}
### Plot:
f = toPoly.pm("x^2 - 3*x + 5")
px = c(0.4, 0.65, 0.75, -0.35 - (1:5)*5/17);
curve(y(x, PFUN=f), from = -3, to = 1, n=512, ylim=c(-15, 15))
line.tan(px, dx=3, p=y, dp=dy, PFUN=f)
# discontinuous
curve(dy(x, PFUN=f), add=T, col="green")
line.tan(px, dx=1.6, p=dy, dp=d2y, PFUN=f, col="orange")


########################
### Integral-Product ###
########################

### I(e^y) dx * I(e^(-y)) dx = P(x)

# - trivial: y^1;
#   [but interesting approach to some integrals]

### D =>
e^y * In + e^(-y) * Ip - dp # = 0
### D2 =>
e^y * dy * In - e^(-y) * dy * Ip - d2p + 2 # = 0
### Solve liniar system =>
e^y * In = (dp*dy + d2p - 2) / (2*dy)
e^(-y) * Ip = (dp*dy - d2p + 2) / (2*dy)

### ODE:
(dp*dy + d2p - 2)*(dp*dy - d2p + 2) - 4*p*dy^2 # = 0
(dp^2 - 4*p)*dy^2 - (d2p - 2)^2 # = 0

### Examples:
### P(x) = x
(1 - 4*x)*dy^2 - 4 # = 0
# y = sqrt(1 - 4*x)

### Plot
# TODO: proper sign pairing;
y = function(x, posRoot=TRUE) {
	val = sqrt(1 - 4*x)
	if( ! posRoot) val = - val;
	return(val)
}
y.exp = function(x, posInt=TRUE) {
	val = sqrt(1 - 4*x)
	if( ! posInt) val = - val;
	val = exp(val);
	return(val)
}
dy = function(x, posRoot=TRUE) {
	dp = 2 / sqrt(1 - 4*x)
	if( ! posRoot) dp = - dp;
	return(dp)
}
Ip.f = function(x) {
	dy.x = dy(x)
	(dy.x + 2) * y.exp(x) / (2*dy.x)
}
Ip.int = function(x, upper=0, diff=3/4) {
	val = sapply(x - diff,
		function(lower) {
		integrate(y.exp, lower=lower, upper=upper)$value
		})
}
y.int = function(lower, upper=0) {
	val = sapply(lower,
		function(lower) {
		integrate(y.exp, lower=lower, upper=upper)$value *
		integrate(y.exp, lower=lower, upper=upper, posInt=FALSE)$value # - upper
	} )
}
### Plot:
lim = c(-7, 1/4)
# interesting integral
# D ( 1/2*(1-sqrt(1-4*x))*exp(sqrt(1-4*x))) = exp(sqrt(1-4*x));
# Note: switched signs; [TODO]
curve(Ip.f(x), from=lim[1], to=lim[2])
curve(Ip.int(x, upper=1/4, diff=3/4), add=T, col="green")
# TODO:
# curve(y.int(x), from=lim[1], to=lim[2])


#########################
#########################

####################
### Derived from ###
###  other ODEs  ###
####################

# - moved to file: DE.ODE.FromODEs.R;

### Basic Example:
### dy - G(x)*y = F(x)

### D =>
d2y - g*dy - dg*y - df # = 0

### Variant 1:
d2y - g*(g*y + f) - dg*y - df # = 0
d2y - (g^2 + dg)*y - g*f - df # = 0

### Variant 2:
y*d2y - (dy - f)*dy - dg*y^2 - df*y # = 0
y*d2y - dy^2 + f*dy - dg*y^2 - df*y # = 0

### Variant 3:
g*d2y - g^2*(g*y + f) - dg*(dy - f) - g*df # = 0
g*d2y - dg*dy - g^3*y - g^2*f + dg*f - g*df # = 0

### Variant 4:
y*d2y - (dy - f)*(g*y + f) - dg*y^2 - df*y # = 0
y*d2y - g*y*dy - f*dy - dg*y^2 + (g*f - df)*y + f^2 # = 0


### Examples:
# G(x) = exp(x^2) =>
y*d2y - dy^2 + f*dy - 2*x*exp(x^2)*y^2 - df*y # = 0
### Extra-Variant:
y*d2y - dy^2 + f*dy - 2*x*(dy - f)*y - df*y # = 0
y*d2y - dy^2 - 2*x*y*dy + f*dy + 2*x*f*y - df*y # = 0


### Example 1:
# F(x) = exp(x^3) + b0
# G(x) = -3*x^2
d2y - (9*x^4 - 6*x)*y + 3*b0*x^2 # = 0

### Example 2:
# G(x) = -1/x;
x^2*d2y + x*dy + y - x^2*df # = 0
# F(x) = ln(x^n + b0);
x^2*d2y + x*dy + y - n*x^(n+1)/(x^n + b0) # = 0


########################
########################

### Higher Powers & y-Exponentials:
# 1.) dy - G(x)*y^2 = F(x)
# 2.) dy - G2(x)*y^2 - G1(x)*y = F(x)
# 3.) x*y*dy + y^2 - exp(y) # = 0

# - moved to file: DE.ODE.FromODEs.R;


#########################
#########################
