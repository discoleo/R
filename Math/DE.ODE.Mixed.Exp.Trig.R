#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Mixed Exp & Trig
##
## draft v.0.1b


### History:

### draft v.0.1a:
# - moved to this file from:
#   DE.ODE.Gaussian.R (Section B.2: Mixed)


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


######################

######################
### Mixed Exp-Trig ###

### Examples:

### y = sin(x^n) * exp(-x^n)

### Check:
ye = expression(sin(x^n) * exp(-x^n))[[1]]
x = sqrt(3); n = sqrt(2);
params = list(x=x, n=n);
#
y = eval(ye, params); dy = eval(D(ye, "x"), params);
d2y = eval(D(D(ye, "x"), "x"), params);


### ODE:
x*d2y - (n-1)*dy + 2*n*x^n*dy + 2*n^2*x^(2*n-1)*y # = 0


### D(y)
dy - n*x^(n-1)*cos(x^n)*exp(-x^n) + n*x^(n-1)*sin(x^n)*exp(-x^n) # = 0
dy - n*x^(n-1)*cos(x^n)*exp(-x^n) + n*x^(n-1) * y # = 0
# =>
n*x^(n-1)*cos(x^n)*exp(-x^n) # =
dy + n*x^(n-1) * y

### D2(y)
d2y # =
n*(n-1)*x^(n-2)*cos(x^n)*exp(-x^n) +
	- n^2*x^(2*n-2)*sin(x^n)*exp(-x^n) +
	- n^2*x^(2*n-2)*cos(x^n)*exp(-x^n) +
	- n*x^(n-1) * dy - n*(n-1)*x^(n-2) * y
(n-1)*(dy + n*x^(n-1) * y) / x +
	- n*x^(n-1) * (dy + n*x^(n-1) * y) +
	- n*x^(n-1) * dy +
	- n^2*x^(2*n-2) * y - n*(n-1)*x^(n-2) * y


### Examples
### n = 2
x*d2y + (4*x^2 - 1)*dy + 8*x^3*y # = 0


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

### y = x * sin(x^n) * exp(-k*x^n)

### Check:
ye = expression(x * sin(x^n) * exp(-k*x^n))[[1]]
x = sqrt(3); n = sqrt(2); k = 1/sqrt(5);
params = list(x=x, n=n, k=k);
#
y = eval(ye, params); dy = eval(D(ye, "x"), params);
d2y = eval(D(D(ye, "x"), "x"), params);


### ODE:
x^2*d2y + x*(2*k*n*x^n - (n+1))*dy +
	+ ((k^2+1)*n^2*x^(2*n) - 2*k*n*x^n + n + 1)*y # = 0


### Derivation:

### D(y):
x*dy - y - n*x^(n+1)*cos(x^n)*exp(-k*x^n) + k*n*x^n*y # = 0

### D2(y):
x*d2y +
	- n*(n+1)*x^n*cos(x^n)*exp(-k*x^n) +
	+ n^2*x^(2*n)*sin(x^n)*exp(-k*x^n) +
	+ k*n^2*x^(2*n)*cos(x^n)*exp(-k*x^n) +
	+ k*n*x^n*dy + k*n^2*x^(n-1)*y # = 0
x^2*d2y + x*(2*k*n*x^n - (n+1))*dy +
	+ ((k^2+1)*n^2*x^(2*n) - 2*k*n*x^n + n + 1)*y # = 0


### Special Cases:

### k = 1
x^2*d2y + x*(2*n*x^n - (n+1))*dy +
	+ (2*n^2*x^(2*n) - 2*n*x^n + n + 1)*y # = 0

### n = -1
x^4*d2y - 2*k*x^2*dy + (2*k*x + k^2 + 1)*y # = 0

