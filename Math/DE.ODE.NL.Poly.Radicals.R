########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs: Polynomial types
## with Radicals
##
## draft v.0.1a



### Examples:

# 2*x*y*dy - n*y^2 + n = 0;
# x*y^2*d2y - 2*dy # = 0
# (x^6+1)*dy^2 - 4*x^6 * y = 8*x^8;


####################
####################

### Helper Functions

# include: DE.ODE.Helper.R;
source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


################
################

################
### Radicals ###
################

### Simple Radicals

### y = sqrt(x^n + 1)
# dy = n/2 * x^(n-1) / y # * x
# 2*x*dy = n * x^n / y
# 2*x*y*dy = n * (y^2 - 1)

### D() =>
# 2*x*y*d2y + 2*x*(dy)^2 + 2*y*dy = 2*n*y*dy
### Variant 1:
# x*y*d2y + x*(dy)^2 - (n-1)*y*dy # = 0

### Variant 2:
# Eq 1 * y =>
# x*y^2*d2y + x*y*(dy)^2 - (n-1)*y^2*dy # = 0
# x*y^2*d2y + n/2 * (y^2 - 1)*dy - (n-1)*y^2*dy # = 0
# x*y^2*d2y - (n-2)/2 * y^2*dy - n/2 * dy # = 0
# 2*x*y^2*d2y - (n-2)*y^2*dy - n*dy # = 0

### Variant 3:
# Eq[2] + D1:
# 2*x*y^2*d2y - (n-2)*y^2*dy + 2*x*y*dy - n*dy - n*y^2 + n # = 0

### Examples:

### n = 2
x*y^2*d2y - 2*dy # = 0
### Variant 3:
x*y^2*d2y + x*y*dy - dy - y^2 + 1 # = 0

### Solution:
y = function(x, n=4, b0=1) {
	f.x = x^n + b0
	y.f = sqrt(f.x)
	return( y.f )
}
dy = function(x, n=4, b0=1, variant=FALSE) {
	y.x = y(x, n=n, b0=b0)
	div = 2*x*y.x
	dp = n * (y.x^2 - b0)
	dp = ifelse(div != 0, dp / div, 0); # TODO: check
	return(dp)
}
d2y = function(x, n=4, b0=1, variant=FALSE) {
	# TODO: b0
	# uses only equation for b0=1;
	y.x = y(x, n=n, b0=b0);
	dy.x = dy(x, n=n, b0=b0);
	div = 2*x*y.x^2;
	if(variant) {
		dp = (n-2)*y.x^2*dy.x - 2*x*y.x*dy.x + n*dy.x + n*y.x^2 - n
	} else {
		dp = ((n-2)*y.x^2 + n)*dy.x;
	}
	lim = if(n > 2) 0 else 1; # TODO
	dp = ifelse(div != 0, dp / div, lim);
	return(dp)
}
### Plot
n = 4;
# nice global minimum
curve(y(x, n=n), from=-3, to=3, ylim=c(-2, 8))
line.tan(c(-2:2), dx=3, p=y, dp=dy, n=n)
# D(y)
curve(dy(x, n=n), add=T, col="green")
line.tan(c(-3:3)/2.7, dx=3, p=dy, dp=d2y, n=n, col="orange")


### Variant
n = 2;
# nice global minimum
curve(y(x, n=n), from=-3, to=3, ylim=c(-2, 4))
line.tan(c(-2:2), dx=3, p=y, dp=dy, n=n)
# D(y)
curve(dy(x, n=n), add=T, col="green")
line.tan(c(-3:3)/2.4, dx=3, p=dy, dp=d2y, n=n, variant=T, col="orange")


######################
######################

######################

### Cardano-type Roots

### Pow = 6
f0 = 0; df0 = 0; # for simplicity
x = sqrt(3); # Test

y = (sqrt(x^6 + 1) + 1)^(2/3) + (sqrt(x^6 + 1) - 1)^(2/3) + f0;
dy = 2 * x^5 / sqrt(x^6 + 1) * 
	((sqrt(x^6 + 1) + 1)^(-1/3) + (sqrt(x^6 + 1) - 1)^(-1/3)) + df0;

### ODE:
dy - 2*x^3 / sqrt(x^6 + 1) * sqrt(y + 2*x^2 - f0) - df0 # = 0

# Alternative Eq:
(x^6+1)*dy^2 - 2*(x^6+1)*df0*dy - 4*x^6 * y +
	+ (x^6+1)*df0^2 + 4*x^6*f0 - 8*x^8 # = 0


### Special Cases:

###
f0 = x^2; df0 = 2*x;
(x^6+1)*dy^2 - 4*x*(x^6+1)*dy - 4*x^6 * y + 4*x^2 # = 0

###
f0 = 0; df0 = 0;
(x^6+1)*dy^2 - 4*x^6 * y - 8*x^8 # = 0


### Test:
f0  = \(x) x^2; df0 = \(x) 2*x;
yf  = \(x) (sqrt(x^6 + 1) + 1)^(2/3) + (sqrt(x^6 + 1) - 1)^(2/3) + f0(x);
dyf = \(x) 2 * x^5 / sqrt(x^6 + 1) * 
	((sqrt(x^6 + 1) + 1)^(-1/3) + (sqrt(x^6 + 1) - 1)^(-1/3)) + df0(x);
dya = \(x, eps = 1E-4) (yf(x + eps) - yf(x)) / eps;

x = c(1/3, 1/2, 1, 3/2);
dyf(x); dya(x);


