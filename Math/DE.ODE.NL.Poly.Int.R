########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs: Polynomial types
## Derived from Integrals of Poly-Fractions
##
## draft v.0.1a


# Note:
# - Derivation of the ODE does NOT require
#   the explicit computation of the Integral;


####################
####################

### Helper Functions

# include: DE.ODE.Helper.R;
source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")

dya = function(x, FUN, eps = 1E-4) {
	(FUN(x + eps) - FUN(x)) / eps;
}


#################

#################
### Integrals ###
#################

### y = x * I(1 / (x^n + k) dx) + F0(x)

### D(y)
# dy =
(y - f0)/x + x / (x^n + k) + df0
# x*dy =
(y - f0) + x^2 / (x^n + k) + x*df0
# x^n =
x^2 / (x*dy - y + f0 - x*df0) - k
# or: 1/(x^n + k) = (x*dy - y + f0 - x*df0) / x^2;

### D2(y)
# x*d2y + dy =
dy + 2*x / (x^n + k) - n*x^(n + 1) / (x^n + k)^2 + x*d2f0;
# Double Entanglement:
# x^2*d2y =   # various variants;
2*(x*dy - y + f0 - x*df0) - n*x^(n + 2) / (x^n + k)^2 + x^2*d2f0;
2*(x*dy - y + f0 - x*df0) - n*x^n*(x*dy - y + f0 - x*df0) / (x^n + k) + x^2*d2f0;
2*(x*dy - y + f0 - x*df0) - n*x^(n - 2)*(x*dy - y + f0 - x*df0)^2 + x^2*d2f0;


### ODE:
x^2*d2y +
	- 2*(x*dy - y + f0 - x*df0) + n*x^(n - 2)*(x*dy - y + f0 - x*df0)^2 - x^2*d2f0;

### Note:
# - it is possible to leave n*x^(n + 1) / (x^n + k)^2 as such;
# - it is possible to express x^n as a function of: dy, y;
# - this latter case becomes interesting when n = sqrt(n0):
x^4*d2y - n*k*(x*dy - y + f0 - x*df0)^2 +
	+ (n - 2)*x^2*(x*dy - y + f0 - x*df0) - x^4*d2f0;


### Solution & Plot:
y.Int = function(x, n, k, lower=0) {
	sapply(x, function(upper)
		integrate(function(x) 1/(x^n + k), lower=lower, upper=upper)$value)
}
y = function(x, b0=1, n=2, k=1, lower=0, asList=FALSE) {
	I.x = y.Int(x, n=n, k=k, lower=lower);
	f0 = eval.pol(x, b0);
	y.x = x*I.x + f0;
	if(asList) return(list(y=y.x, f0=f0));
	return(y.x)
}
dy = function(x, b0=1, n=2, k=1, lower=0, asList=FALSE) {
	y.all = y(x, b0=b0, n=n, k=k, lower=lower, asList=TRUE)
	y.x = y.all$y; f0 = y.all$f0;
	df0 = deriv.pol(x, b0, dn=1);
	div = x;
	dp = (y.x - f0) + x^2 / (x^n + k) + x*df0
	dp = ifelse(div != 0, dp / div,
		dy(x + 1E-3, b0=b0, n=n, k=k, lower=lower)); # TODO: check
	if(asList) return(list(y=y.x, dy=dp, f0=f0, df0=df0));
	return(dp)
}
d2y = function(x, b0=1, n=2, k=1, lower=0) {
	y.all = dy(x, b0=b0, n=n, k=k, lower=lower, asList=TRUE);
	y.x = y.all$y; dy.x = y.all$dy;
	f0 = y.all$f0; df0 = y.all$df0; d2f0 = deriv.pol(x, b0, dn=2);
	Tx = (x*dy.x - y.x + f0 - x*df0);
	dp = n*k*Tx^2 - (n - 2)*x^2*Tx + x^4*d2f0;
	div = x^4;
	dp = ifelse(div != 0, dp / div,
		d2y(x + 1E-3, b0=b0, n=n, k=k, lower=lower)); # TODO
	return(dp)
}
### Plot
b0 = c(1, 1); n = 2; k = 1;
px = c(0:3, 5) * 2/5 - 1/2;
curve(y(x, b0=b0, n=n, k=k), from=-0.5, to=3, ylim=c(-0.5, 6.5))
line.tan(px, dx=2.5, p=y, dp=dy, b0=b0, n=n, k=k)
# D(y)
curve(dy(x, b0=b0, n=n, k=k), add=T, col="green")
line.tan(px, dx=2.5, p=dy, dp=d2y, b0=b0, n=n, k=k, col="orange")


### Ex 2:
b0 = c(1, -1); n = 2; k = 3;
px = c(0:3, 5) * 2/5 - 1/2;
curve(y(x, b0=b0, n=n, k=k), from=-0.5, to=3, ylim=c(-2, 2))
line.tan(px, dx=2.5, p=y, dp=dy, b0=b0, n=n, k=k)
# D(y)
curve(dy(x, b0=b0, n=n, k=k), add=T, col="green")
line.tan(px, dx=2.5, p=dy, dp=d2y, b0=b0, n=n, k=k, col="orange")


### Ex 3:
b0 = c(1, -1); n = 5^(1/5); k = 3;
px = c(0:3, 5) * 2/5;
curve(y(x, b0=b0, n=n, k=k), from=0, to=3, ylim=c(-1.2, 1.2))
line.tan(px, dx=2.5, p=y, dp=dy, b0=b0, n=n, k=k)
# D(y)
curve(dy(x, b0=b0, n=n, k=k), add=T, col="green")
line.tan(px, dx=2.5, p=dy, dp=d2y, b0=b0, n=n, k=k, col="orange")

### [cont]
# nice global minimum
px = c(0.5, 1:4)*2 - 1/2;
curve(y(x, b0=b0, n=n, k=k), from=0, to=10, ylim=c(-1.2, 2.5))
line.tan(px, dx=2.5, p=y, dp=dy, b0=b0, n=n, k=k)
# D(y)
curve(dy(x, b0=b0, n=n, k=k), add=T, col="green")
line.tan(px, dx=2.5, p=dy, dp=d2y, b0=b0, n=n, k=k, col="orange")


#########################
#########################

##################
### Integrals: ###
### M-Inverse  ###
##################

### y * I(1 / (x^n + k) dx) = F0(x)


### D(y)
I*dy + y / (x^n + k) - df0 # = 0 # * y =>
f0*dy + y^2 / (x^n + k) - df0*y # = 0

# y^2 / (x^n + k) = - (f0*dy - df0*y);
# x^n = - y^2 / (f0*dy - df0*y) - k;

### D2(y)
f0*d2y + 2*y*dy / (x^n + k) - n*x^(n-1)*y^2 / (x^n + k)^2 - d2f0*y # = 0
f0*y*d2y - 2*(f0*dy - df0*y)*dy - n*x^(n-1)*y^3 / (x^n + k)^2 - d2f0*y^2 # = 0
x*f0*y^2*d2y - 2*x*(f0*dy - df0*y)*y*dy - n*x^n*(f0*dy - df0*y)^2 - x*d2f0*y^3 # = 0
x*f0*y^2*d2y - 2*x*(f0*dy - df0*y)*y*dy +
	+ n*(y^2 / (f0*dy - df0*y) + k)*(f0*dy - df0*y)^2 - x*d2f0*y^3 # = 0

### ODE:
x*f0*y^2*d2y - 2*x*(f0*dy - df0*y)*y*dy +
	+ n*y^2*(f0*dy - df0*y) + n*k*(f0*dy - df0*y)^2 - x*d2f0*y^3 # = 0


### Solution & Plot:
y.Int = function(x, n, k, lower=0) {
	sapply(x, function(upper)
		integrate(function(x) 1/(x^n + k), lower=lower, upper=upper)$value)
}
y = function(x, b0=1, n=2, k=1, lower=0, asList=FALSE) {
	I.x = y.Int(x, n=n, k=k, lower=lower);
	f0 = eval.pol(x, b0);
	y.x = f0; div = I.x;
	y.x = ifelse(div != 0, y.x/div, 100) # TODO
	if(asList) return(list(y=y.x, f0=f0));
	return(y.x)
}
dy = function(x, b0=1, n=2, k=1, lower=0, asList=FALSE) {
	y.all = y(x, b0=b0, n=n, k=k, lower=lower, asList=TRUE)
	y.x = y.all$y; f0 = y.all$f0;
	df0 = deriv.pol(x, b0, dn=1);
	div = f0;
	dp =  - y.x^2 / (x^n + k) + df0*y.x;
	dp = ifelse(div != 0, dp / div,
		dy(x + 1E-3, b0=b0, n=n, k=k, lower=lower)); # TODO: check
	if(asList) return(list(y=y.x, dy=dp, f0=f0, df0=df0));
	return(dp)
}
d2y = function(x, b0=1, n=2, k=1, lower=0) {
	y.all = dy(x, b0=b0, n=n, k=k, lower=lower, asList=TRUE);
	y.x = y.all$y; dy.x = y.all$dy;
	f0 = y.all$f0; df0 = y.all$df0; d2f0 = deriv.pol(x, b0, dn=2);
	Tx = (f0*dy.x - df0*y.x);
	dp = 2*x*Tx*y.x*dy.x - n*Tx*y.x^2 - n*k*Tx^2 + x*d2f0*y.x^3;
	div = x*f0*y.x^2;
	dp = ifelse(div != 0, dp / div,
		d2y(x + 1E-3, b0=b0, n=n, k=k, lower=lower)); # TODO
	return(dp)
}
### Plot
b0 = c(1, 1); n = 2; k = 1;
px = - c(4:0) * 2/5 - 1/2;
curve(y(x, b0=b0, n=n, k=k), from=-2, to=0, ylim=c(-5, 3))
line.tan(px, dx=2.5, p=y, dp=dy, b0=b0, n=n, k=k)
# D(y)
curve(dy(x, b0=b0, n=n, k=k), add=T, col="green")
line.tan(px, dx=2.5, p=dy, dp=d2y, b0=b0, n=n, k=k, col="orange")


### Ex 2:
b0 = c(1, 1); n = 2; k = 1;
px = c(0:3, 5) * 2/5 + 1/2;
curve(y(x, b0=b0, n=n, k=k), from=0, to=3, ylim=c(-2.5, 5.5))
line.tan(px, dx=2.5, p=y, dp=dy, b0=b0, n=n, k=k)
# D(y)
curve(dy(x, b0=b0, n=n, k=k), add=T, col="green")
line.tan(px, dx=2.5, p=dy, dp=d2y, b0=b0, n=n, k=k, col="orange")


### Ex 3:
b0 = c(-1/2, 2); n = 5^(1/5); k = 1;
px = c(0:3, 5) * 2/5 + 1/4;
curve(y(x, b0=b0, n=n, k=k), from=0+1E-3, to=3, ylim=c(-2.5, 5.5))
line.tan(px, dx=2.5, p=y, dp=dy, b0=b0, n=n, k=k)
# D(y)
curve(dy(x, b0=b0, n=n, k=k), add=T, col="green")
line.tan(px, dx=2.5, p=dy, dp=d2y, b0=b0, n=n, k=k, col="orange")

