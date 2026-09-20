########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Power-Log
##
## draft v.0.1b

### Power with Log
# y = P1(x)^log(P2(x))

### Examples:

# x^2*y*d2y - x^2*dy^2 + x*y*dy - 2*y^2 = 0
# x^2*y*d2y - x^2*dy^2 + x*y*dy - k*y^2 = 0 # Gen


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y = x^log(x)

# Check:
x = sqrt(3); params = list(x=x);
e = expression(x^log(x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*y*d2y - x^2*dy^2 + x*y*dy - 2*y^2 # = 0


### D(y)
dy - 2*log(x)/x * x^log(x) # = 0
dy - 2*log(x)/x * y # = 0

### D2(y)
d2y # ==
2*log(x)/x * dy + 2*(1 - log(x))/x^2 * y
### x^2 * D2(y)
x^2 * d2y # ==
2*x*log(x) * dy + 2*(1 - log(x)) * y
x^2/y * dy^2 - x*dy + 2*y


##################
### Generalization

### y = x^p * x^(k*log(x))
# Note: p has no impact on ODE;
# (except for initial conditions)

# Check:
p = 2/5; # p = 7/3;
k = sqrt(5/3);
x = sqrt(3); params = list(x=x, k=k, p=p);
e = expression(x^p * x^(k*log(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2 * y*d2y - x^2 * dy^2 + x*y*dy - 2*k*y^2 # = 0

# D =>
x*dy - p*y - 2*k*log(x)*y # = 0

# D2 =>
x^2 * d2y - (p-1)*x*dy +
	- 2*k*x*log(x)*dy - 2*k*y # = 0
x^2 * y*d2y - (p-1)*x * y*dy +
	- x*(x*dy - p*y) * dy - 2*k*y^2 # = 0


#############
### Extension
### y = (x + k)^log(x + k)

### D(y)
2*log(x+k)*y / (x+k)

### ODE:
(x+k)^2*y*d2y - (x+k)^2*dy^2 + (x+k)*y*dy - 2*y^2 # = 0


### Solution & Plot:
y = function(x, k=0) {
	xk = if(k == 0) x else x+k;
	val = xk^log(xk);
	return(val)
}
dy = function(x, k=0) {
	xk = x + k;
	logx = log(xk);
	dp = 2*logx * (xk)^logx
	div = xk;
	dp = ifelse(div != 0, dp/div, 1); # TODO
	return(dp)
}
d2y = function(x, k=0) {
	### D()
	y.x  =  y(x, k=k);
	dy.x = dy(x, k=k);
	xk = if(k == 0) x else x+k;
	div = xk^2 * y.x;
	d2p = xk^2*dy.x^2 - xk*y.x*dy.x + 2*y.x^2
	d2p = ifelse(div != 0, d2p/div, 1); # TODO
	return(d2p)
}
### Plot:
px = 3/7 + (0:4)*2/7
curve(y(x), from= 0+1E-1, to = 2.5, ylim=c(-2,5))
line.tan(px, dx=3, p=y, dp=dy)
#
curve(dy(x), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, col="orange")


### Ex 2:
k = 2
px = 3/7 + (0:4)*2/7 - k;
curve(y(x, k=k), from= -k+1E-1, to = 2.5, ylim=c(-2,5))
line.tan(px, dx=3, p=y, dp=dy, k=k)
#
curve(dy(x, k=k), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, col="orange")


### Extensions: Higher Order
### y = (x^n + k)^log(x^n + k)

### D(y)
2*n*log(x^n+k)*x^(n-1) / (x^n+k) * y

### D2(y)
2*n*log(x^n+k)*x^(n-1) / (x^n+k) * dy +
	+ 2*n^2*x^(2*n-2) / (x^n+k)^2 * y +
	+ 2*n*log(x^n+k)*((n-1)*x^(n-2)*(x^n+k) - n*x^(2*n-2)) / (x^n+k)^2 * y
(dy)^2 / y +
	+ 2*n^2*x^(2*n-2) / (x^n+k)^2 * y +
	- (x^n + k - n*k) / (x^n+k) * dy / x

### ODE:
x*(x^n+k)^2*y*d2y - x*(x^n+k)^2*dy^2 + (x^n+k-n*k)*(x^n+k)*y*dy - 2*n^2*x^(2*n-1)*y^2 # = 0

### Case: n = 2
x*(x^2+k)^2*y*d2y - x*(x^2+k)^2*dy^2 + (x^4-k^2)*y*dy - 8*x^3*y^2 # = 0


### Solution & Plot:
y = function(x, k=0, n=2) {
	xk = if(k == 0) x^n else x^n + k;
	val = xk^log(xk);
	return(val)
}
dy = function(x, k=0, n=2) {
	xn = x^n;
	xk = if(k == 0) xn else xn + k;
	logx = log(xk);
	dp = 2*n * logx * xn * (xk)^logx;
	div = x * xk;
	dp = ifelse(div != 0, dp/div, 0); # TODO
	return(dp)
}
d2y = function(x, k=0, n=2) {
	### D()
	y.x  =  y(x, k=k, n=n);
	dy.x = dy(x, k=k, n=n);
	xn = x^n; x2 = x^2;
	xk = if(k == 0) xn else xn + k;
	div = x2 * xk^2 * y.x;
	d2p = x2*xk^2*dy.x^2 - x*(xk - n*k)*xk*y.x*dy.x + 2*n^2*xn^2*y.x^2
	d2p = ifelse(div != 0, d2p/div, d2y(1E-5, k=k, n=n)); # TODO
	return(d2p)
}
### Plot:
n = 2; k = 2;
px = (-3:3)*2/7
curve(y(x, k=k, n=n), from= -2.5, to = 2.5, ylim=c(-2,5))
line.tan(px, dx=3, p=y, dp=dy, k=k, n=n)
# TODO: px = 0;
curve(dy(x, k=k, n=n), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, n=n, col="orange")


### Ex 2:
n = 2; k = 0.9;
px = (-3:3)*2/7
curve(y(x, k=k, n=n), from= -2.5, to = 2.5, ylim=c(-2,5))
line.tan(px, dx=3, p=y, dp=dy, k=k, n=n)
# TODO: px = 0;
curve(dy(x, k=k, n=n), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, n=n, col="orange")


######################
### Generalization ###

### y = P(x)^log(P(x))

### D(y)
2*log(p)*dp / p * y

### D2(y)
# p*d2y + dp*dy
2*log(p)*dp*dy + 2*(dp)^2 / p * y + 2*log(p)*d2p*y
p / y * (dy)^2 + 2*(dp)^2 / p * y + p*d2p/dp * dy

### ODE:
dp*p^2 * y*d2y - dp*p^2 * (dy)^2 + (dp^2*p - d2p*p^2) * y*dy - 2*(dp)^3 * y^2 # = 0

### Examples:
### Check: p = x + k; dp = 1;
(x+k)^2*y*d2y - (x+k)^2*(dy)^2 + (x+k)*y*dy - 2*y^2 # = 0
### p = exp(x) + k; dp = exp(x)
(e^x+k)^2*y*d2y - (e^x+k)^2*(dy)^2 + e^x*(e^x+k)*y*dy - (e^x+k)^2*y*dy - 2*e^(2*x)*y^2
### p = ln(k*x); dp = 1/x
x^2*ln(k*x)^2*y*d2y - x^2*ln(k*x)^2*(dy)^2 + x*ln(k*x)*y*dy + x*ln(k*x)^2*y*dy - 2*y^2


### Solution & Plot:
y = function(x, k=2) {
	# p = exp(x) + k;
	p = if(k == 0) exp(x) else exp(x) + k;
	val = p^log(p);
	return(val)
}
dy = function(x, k=2) {
	# 2*log(p)*dp / p * y
	x.exp = exp(x);
	p = if(k == 0) x.exp else x.exp + k;
	x.log = log(p);
	dp = 2 * x.log * x.exp * (p)^x.log;
	div = p;
	dp = ifelse(div != 0, dp/div, 0); # TODO
	return(dp)
}
d2y = function(x, k=2) {
	### D()
	y.x  =  y(x, k=k);
	dy.x = dy(x, k=k);
	x.exp = exp(x);
	xk = if(k == 0) x.exp else x.exp + k;
	div = xk^2 * y.x;
	d2p = xk^2*(dy.x)^2 - x.exp*xk*y.x*dy.x + xk^2*y.x*dy.x + 2*x.exp*x.exp*y.x^2
	d2p = ifelse(div != 0, d2p/div, 1); # TODO
	return(d2p)
}
### Plot:
k = 2;
px = (-4:2)*3/7
curve(y(x, k=k), from= -2.5, to = 2, ylim=c(-1,10))
line.tan(px, dx=3, p=y, dp=dy, k=k)
#
curve(dy(x, k=k), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, col="orange")


### Ex 2:
k = -1/2;
px = (-1:3)*3/7
curve(y(x, k=k), from= -1/2, to = 2, ylim=c(-2,10))
line.tan(px, dx=3, p=y, dp=dy, k=k)
#
curve(dy(x, k=k), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, col="orange")

