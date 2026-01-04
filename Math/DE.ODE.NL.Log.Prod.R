########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs - Logarithms: Prod
##
## draft v.0.1a


### NL ODEs Derived from Logarithms
# y = PROD( LOG )


###################
### Logarithmic ###
###################

### y^n = log(f(x)) * log(g(x)) + F0(x)

# Check:
# Homogenous: F0 = 0;
n = 2; # n = 1;
x = sqrt(3); b0 = 1/2; c1 = 1/2; c0 = 1;
params = list(x=x, b0=b0, c1=c1, c0=c0);
e = expression((log(x+b0) * log(x^2+c1*x+c0))^(1/n))[[1]];
f = x+b0; g = x^2+c1*x+c0;
df = 1; d2f = 0; dg = 2*x + c1; d2g = 2;
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### D =>
n*y^(n-1)*dy - df * log(g)/f - dg * log(f)/g # = 0 # * f*g
n*f*g*y^(n-1)*dy - g*df * log(g) - f*dg * log(f) # = 0

### D2 =>
n*f*g*y^(n-1)*d2y + n*(n-1)*f*g*y^(n-2)*dy^2 + n*df*g*y^(n-1)*dy + n*f*dg*y^(n-1)*dy # =
	(g*d2f + df*dg) * log(g) + (f*d2g + df*dg) * log(f) + df*dg + dg*df;

### Solve Linear: D1 & D2 =>
# log(f) = ...
# log(g) = ...


### Special Cases:

### Order: n = 1
f*g*dy - g*df * log(g) - f*dg * log(f) # == 0
f*g*d2y + df*g*dy + f*dg*dy - df*dg - dg*df # =
	(g*d2f + df*dg) * log(g) + (f*d2g + df*dg) * log(f);


### Solve Linear: D1 & D2 =>
# Log(f):
(f*dg*(g*d2f + df*dg) - g*df*(f*d2g + df*dg)) * log(f) # =
(g*d2f + df*dg) * f*g*dy - g*df * (f*g*d2y + df*g*dy + f*dg*dy - df*dg - dg*df);

# Log(g):
(g*df*(f*d2g + df*dg) - f*dg*(g*d2f + df*dg)) * log(g) # =
(f*d2g + df*dg) * f*g*dy - dg*f * (f*g*d2y + df*g*dy + f*dg*dy - df*dg - dg*df);

# OR

### Solve Quadratic =>
(g*d2f + df*dg) * log(g)^2 +
	- (f*g*d2y + df*g*dy + f*dg*dy - df*dg - dg*df) * log(g) +
	+ (f*d2g + df*dg) * y # == 0
#
log(f) = ...
log(g) = ...


### Order: n = 1/2
1/2*f*g*y^(-1/2)*dy = g*df * log(g) + f*dg * log(f)
1/2*f*g*y^(-1/2)*d2y - 1/4*f*g*y^(-3/2)*dy^2 + 1/2*df*g*y^(-1/2)*dy + 1/2*f*dg*y^(-1/2)*dy -df - dg =
	(g*d2f + df*dg) * log(g) + (f*d2g + df*dg) * log(f)


### Examples:

### y = log(x + k) * log(x - k)
dy = log(x-k)/(x+k) + log(x+k)/(x-k)
(x^2 - k^2)*dy = (x-k)*log(x-k) + (x+k)*log(x+k)
### D2 =>
(x^2 - k^2)*d2y + 2*x*dy - 2 = log(x-k) + log(x+k)
### Solve Linear: D1 & D2 =>
2*k*log(x-k) = (x+k)*((x^2 - k^2)*d2y + (x + k)*dy - 2)
2*k*log(x+k) = -(x-k)*((x^2 - k^2)*d2y + (x - k)*dy - 2)

### ODE:
(x^2 - k^2)*((x^2 - k^2)*d2y + (x + k)*dy - 2) * ((x^2 - k^2)*d2y + (x - k)*dy - 2) + 4*k^2*y = 0

### Ex 2:
### y^(1/2) = log(x + k) * log(x - k)
1/2*y^(-1/2)*dy = log(x-k)/(x+k) + log(x+k)/(x-k)
(x^2 - k^2)*y^(-1/2)*dy = 2*(x-k)*log(x-k) + 2*(x+k)*log(x+k)
### D2 =>
(x^2 - k^2)*y^(-1/2)*d2y + 2*x*y^(-1/2)*dy - 1/2*(x^2 - k^2)*y^(-3/2)*dy^2 - 4 = 2*log(x-k) + 2*log(x+k)
### Solve Liniar =>
4*k*log(x-k) =  (x+k)*y^(-1/2)*((x^2 - k^2)*d2y + (x + k)*dy - 1/2*(x^2 - k^2)*dy^2/y - 4)
4*k*log(x+k) = -(x-k)*y^(-1/2)*((x^2 - k^2)*d2y + (x - k)*dy - 1/2*(x^2 - k^2)*dy^2/y - 4)

### ODE:
(x^2 - k^2)*(...)*(...) + 4*k^2*y^(3/2) = 0


### Solution & Plot
y = function(x, k=1, n=1, v.dy, v.d2y) {
	if(missing(v.dy)) v.dy = dy(x, k=k, n=n)
	if(missing(v.d2y)) v.d2y = d2y(x, k=k, n=n, v.dy=v.dy)
	x2 = x^2 - k^2
	y = - x2*(x2*v.d2y + (x + k)*v.dy - 2)*(x2*v.d2y + (x - k)*v.dy - 2)
	return(y / 4 / k^2)
}
dy = function(x, k=1, n=1) {
	dp = log(x-k)/(x+k) + log(x+k)/(x-k);
	return(dp)
}
d2y = function(x, k=1, n=1, v.dy) {
	if(missing(v.dy)) v.dy = dy(x, k=k, n=n)
	dp =log(x-k) + log(x+k) + 2 - 2*x*v.dy;
	dp = dp / (x^2 - k^2)
	return(dp)
}
### Plot:
k = 1; n = 1;
px = c(3/5 + (1:3)*3/5);
curve(y(x, k=k, n=n), from= 1.01, to = 3, ylim=c(-3, 3))
# global "minimum" / horn;
line.tan(px, dx=3, p=y, dp=dy, k=k, n=n)
#
curve(dy(x, k=k, n=n), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, n=n, col="orange")


### Ex 2:
k = 3; n = 1;
px = c(3 + (1:3)*3/5);
curve(y(x, k=k, n=n), from= 3.01, to = 6, ylim=c(-3, 3))
# global "minimum" / horn;
line.tan(px, dx=3, p=y, dp=dy, k=k, n=n)
#
curve(dy(x, k=k, n=n), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, n=n, col="orange")

