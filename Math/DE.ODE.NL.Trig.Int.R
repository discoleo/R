########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs - Trigonometric
##
## draft v.0.2a


### Non-Linear: Trigonometric Integrals
# - derived from I(Trig);


###############
### History ###
###############

### draft v.0.2a:
# - [refactor] moved ODEs derived from Fresnel integrals
#   to file: DE.ODE.Trig.Int.Fresnel.R;
### draft v.0.1a:
# - moved I(Trig) variants to this file
#   from file: DE.ODE.Trigonometric.R;

### [previous file]

### draft v.0.4a:
# - derived from Inverse Integrals:
#   y * I(sin(x^n)) = F0(x);
### draft v.0.3h:
# - minor extension to: [v.0.3f];
### draft v.0.3g:
# - derived from:
#   y = x * I(tan(k*x^2)) + F0(x);
### draft v.0.3f - v.0.3f-test:
# - derived from:
#   y = x * I(atan(k*x)^2) + F0(x); [extended in v.0.3h]
### draft v.0.3c-pre:
# - derived from:
#   y = I(sin(x^n)) * I(cos(x^n));


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#################

#################
### Integrals ###
#################

### y = I(sin(x^n)) * I(cos(x^n))

# Check:
n = sqrt(3); x = 2*pi;
Ic = integrate(\(x) cos(x^n), 0, x)$value;
Is = integrate(\(x) sin(x^n), 0, x)$value;
y  = Ic * Is;

### D(y)
dy = sin(x^n) * Ic + cos(x^n) * Is;

### D2(y)
d2y = n*x^(n-1)*cos(x^n) * Ic - n*x^(n-1)*sin(x^n) * Is +
	+ sin(2*x^n);

### Solve Linear system:
div = n*x^(n-1);
Is # ==
(n*x^(n-1)*cos(x^n)*dy - sin(x^n)*d2y + sin(x^n)*sin(2*x^n)) / div;
Ic # ==
(n*x^(n-1)*sin(x^n)*dy + cos(x^n)*d2y - cos(x^n)*sin(2*x^n)) / div;

### ODE:
y - (n*x^(n-1)*cos(x^n)*dy - sin(x^n)*d2y + sin(x^n)*sin(2*x^n)) *
	(n*x^(n-1)*sin(x^n)*dy + cos(x^n)*d2y - cos(x^n)*sin(2*x^n)) / x^(2*n-2) / n^2;

# TODO: check;


############################
############################

### y = x * I(atan(k*x)^2) + F0(x)

### D(y)
(y - f0)/x + x*atan(k*x)^2 + df0
### x*dy =
(y - f0) + x^2*atan(k*x)^2 + x*df0

### D2(y)
# x*d2y + dy =
dy + 2*x*atan(k*x)^2 + 2*k*x^2 / (k^2*x^2+1) * atan(x) + x*d2f0
dy + 2*(x*dy - y + f0 - x*df0)/x +
	+ 2*k*x / (k^2*x^2+1) * sqrt(x*dy - y + f0 - x*df0) + x*d2f0
# x^2*(k^2*x^2+1)*d2y =
2*(k^2*x^2+1)*(x*dy - y + f0 - x*df0) +
	+ 2*k*x^2*sqrt(x*dy - y + f0 - x*df0) + x^2*(k^2*x^2+1)*d2f0

### ODE:
(x^2*(k^2*x^2+1)*d2y - 2*(k^2*x^2+1)*(x*dy - y + f0 - x*df0) - x^2*(k^2*x^2+1)*d2f0)^2 +
	- 4*k^2*x^4*(x*dy - y + f0 - x*df0) # = 0

### Examples:
### f0 = x; df0 = 1;
(x^2*(k^2*x^2+1)*d2y - 2*(k^2*x^2+1)*(x*dy - y))^2 - 4*k^2*x^4*(x*dy - y) # = 0


### Solution & Plot:
y.I = function(x, k=1, n=2, lower=0) {
	sapply(x, function(upper)
		integrate(function(x) atan(k*x)^n, lower=lower, upper=upper)$value)
}
deriv.pol = function(x, b, dn=1, x.mult=0) {
	coeff = tail(b, -1);
	pow = seq(length(coeff));
	coeff = coeff * pow;
	if(dn > 1) return(deriv.pol(x, coeff, dn=dn-1, x.mult=x.mult))
	if(x.mult != 1) pow = pow + x.mult - 1; # x * df0;
	sapply(x, function(x) sum(coeff * x^pow));
}
eval.pol = function(x, b) {
	pow = seq(0, length(b) - 1)
	sapply(x, function(x) sum(b * x^pow))
}
y = function(x, b=0, k=1, n=2, lower=0) {
	# x * I(atan(x)^2) + F0(x)
	I.v = y.I(x, k=k, n=n, lower=lower)
	r = x*I.v + eval.pol(x, b);
	return(r)
}
dy = function(x, b=0, k=1, n=2, lower=0) {
	I.v = y.I(x, k=k, n=n, lower=lower)
	# (y - f0) + x^2*atan(k*x)^2 + x*df0
	# y - f0 = x*I.v;
	r = x*I.v + x^2 * atan(k*x)^2 + deriv.pol(x, b, x.mult=1); # x*df0
	div = x;
	r = ifelse(div != 0, r / div,
		dy(x + 1E-3, b=b, k=k, n=n, lower=lower)) # TODO
	return(r)
}
d2y = function(x, b=0, k=1, n=2, lower=0) {
	y.x = y(x, b=b, k=k, n=n, lower=lower)
	dy.x = dy(x, b=b, k=k, n=n, lower=lower)
	# (x^2*(k^2*x^2+1)*d2y - 2*(k^2*x^2+1)*(x*dy - y + f0 - x*df0) - x^2*(k^2*x^2+1)*d2f0)^2 +
	#  - 4*k^2*x^4*(x*dy - y + f0 - x*df0) # = 0
	x2 = x^2; x21 = k^2 * x2 + 1;
	f0 = eval.pol(x, b); xdf0 = deriv.pol(x, b, dn=1, x.mult=1);
	d2f0 = deriv.pol(x, b, dn=2, x.mult=0);
	Tx = x*dy.x - y.x + f0 - xdf0;
	dp  = 4*k^2*x2*x2*Tx;
	dp = sqrt(dp);
	dp = dp + 2*x21*Tx + x2*x21*d2f0;
	div = x2 * x21;
	dp = ifelse(div != 0, dp / div, d2y(x + 1E-3, b=b, n=n, lower=lower)); # TODO
	return(dp)
}
### Plot:
b = c(0, 1);
x.px = c(-5:7 * 3/7)
curve(y(x, b=b), from= -3, to = 2, ylim=c(-2, 3.5))
line.tan(x.px, dx=1.6, p=y, dp=dy, b=b)
#
curve(dy(x, b=b), add=T, col="green")
line.tan(x.px, dx=1.5, p=dy, dp=d2y, b=b, col="orange")


### Ex 2:
b = c(0, -2, 1);
x.px = c(-5:7 * 3/7)
curve(y(x, b=b), from= -1.5, to = 2.5, ylim=c(-3, 4))
line.tan(x.px, dx=1.6, p=y, dp=dy, b=b)
#
curve(dy(x, b=b), add=T, col="green")
line.tan(x.px, dx=1.5, p=dy, dp=d2y, b=b, col="orange")


### Ex 3:
b = c(0, -2, 1); k =1/3;
x.scale = 1; xlim = c(-1.5, 2.5) * x.scale;
x.px = c(-5:7 * 3/7) * x.scale;
curve(y(x, b=b, k=k), from= xlim[1], to = xlim[2], ylim=c(-3, 4))
line.tan(x.px, dx=1.6, p=y, dp=dy, b=b, k=k)
#
curve(dy(x, b=b, k=k), add=T, col="green")
line.tan(x.px, dx=1.5, p=dy, dp=d2y, b=b, k=k, col="orange")


#########################
#########################

#########################
### Inverse Integrals ###
#########################

### y * I(sin(x^n)) = F0(x)

### D(y)
I*dy + y*sin(x^n) - df0 # = 0 # * y =>
f0*dy + sin(x^n)*y^2 - y*df0 # = 0

### D2(y)
# Note: basic equation contains cos(x^n);
f0*d2y + 2*sin(x^n)*y*dy + n*x^(n-1)*cos(x^n)*y^2 - y*d2f0 # = 0
f0*y*d2y - 2*(f0*dy - y*df0)*dy + n*x^(n-1)*cos(x^n)*y^3 - y^2*d2f0 # = 0
f0*y*d2y - 2*(f0*dy - y*df0)*dy - y^2*d2f0 + n*x^(n-1)*sqrt(y^4 - (f0*dy - y*df0)^2)*y # = 0
(f0*y*d2y - 2*(f0*dy - y*df0)*dy - y^2*d2f0)^2 +
	- n^2*x^(2*n-2)*(y^4 - (f0*dy - y*df0)^2)*y^2 # = 0

### TODO: check;

