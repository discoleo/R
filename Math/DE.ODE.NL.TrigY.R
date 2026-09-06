#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Trigonometric Variants
##
## draft v.0.2a


### Non-Linear ODEs:
### Trigonometric Variants of y

### Theory:
# f1(x)*sin(y) + f2(x)*cos(y) = F0(x)


####################

### Helper Functions


source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#########################
#########################

###############
### Trig(y) ###
###############

### Combined Functions:
### Linearly combined

### f1(x)*sin(y) + f2(x)*cos(y) = F0(x)

# Check:
ye = expression(asin(1/sqrt(x^2+2)))[1]
fe = expression(x/sqrt(x^2+2) + x^2*sqrt(x^2+1)/sqrt(x^2+2))
x = sqrt(3);
params = list(x=x);
#
f1 = x;   df1 = 1; d2f1 = 0;
f2 = x^2; df2 = 2*x; d2f2 = 2;
#
y = eval(ye, params); dy = eval(D(ye, "x"), params);
f = eval(fe, params); df = eval(D(fe, "x"), params);
d2y = eval(D(D(ye, "x"), "x"), params);
d2f = eval(D(D(fe, "x"), "x"), params);
p = f; dp = df; d2p = d2f;


# f2(x)*cos(y) + f1(x)*sin(y) # = F0(x)
f2*cos(y) + f1*sin(y) - f # = 0

### D =>
(f1*dy + df2)*cos(y) - (f2*dy - df1)*sin(y) - dp # = 0
### Solve linear =>
sin(y) - ((f1*dy + df2)*p - f2*dp) / ((f1*dy + df2)*f1 + (f2*dy - df1)*f2) # = 0
cos(y) - ((f2*dy - df1)*p + f1*dp) / ((f2*dy - df1)*f2 + (f1*dy + df2)*f1) # = 0

### D2 =>
(f1*d2y + df1*dy + d2f2)*cos(y) - (f1*dy + df2)*dy*sin(y) +
	- (f2*d2y + df2*dy - d2f1)*sin(y) - (f2*dy - df1)*dy*cos(y) - d2p # = 0
(f1*d2y - f2*dy^2 + 2*df1*dy + d2f2)*cos(y) +
	- (f2*d2y + f1*dy^2 + 2*df2*dy - d2f1)*sin(y) - d2p # = 0
(f1*d2y - f2*dy^2 + 2*df1*dy + d2f2)*((f2*dy - df1)*p + f1*dp) +
	- (f2*d2y + f1*dy^2 + 2*df2*dy - d2f1)*((f1*dy + df2)*p - f2*dp) +
	- d2p*((f1*dy + df2)*f1 + (f2*dy - df1)*f2) # = 0

### ODE:
# - Riccati-type NL ODE - but w. Power = 3;
#   Note: NO y, as y was only present as Trig(y);
((f1^2+f2^2)*dp - (df1*f1 + f2*df2)*p) * d2y +
	- p*(f1^2 + f2^2) * dy^3 + 3*p*(f2*df1 - f1*df2) * dy^2 +
	+ (p*(f1*d2f1 + d2f2*f2) - 2*p*(df1^2 + df2^2) + 2*dp*(df1*f1 + f2*df2) +
		- (f1^2+f2^2)*d2p) * dy +
	+ p*(d2f1*df2 - d2f2*df1) + dp*(d2f2*f1 - d2f1*f2) + d2p*(f2*df1 - f1*df2) # = 0


### Examples:

### Ex: f1(x) = x; f2(x) = 1;
# Check:
fe = expression(x/sqrt(x^2+2) + sqrt(x^2+1)/sqrt(x^2+2))
f1 = x; df1 = 1; d2f1 = 0;
f2 = 1; df2 = d2f2 = 0;
#
(x*d2y - dy^2 + 2*dy)*(p*dy - p + x*dp) - (d2y + x*dy^2)*(x*p*dy - dp) +
	- d2p*(x^2*dy + dy - 1) # = 0
### ODE:
(dp*x^2 - p*x + dp)*d2y - p*(x^2 + 1)*dy^3 + 3*p*dy^2 - (x^2*d2p - 2*x*dp + 2*p + d2p)*dy + d2p # = 0
# Case: F0(x) = x^2;
(x^3 + 2*x)*d2y - x^2*(x^2 + 1)*dy^3 + 3*x^2*dy^2 - 2*dy + 2 # = 0



### Ex: f1 = x; f2 = x^2;
((x^4 + x^2)*dp - (2*x^3 + x)*p) * d2y +
	- p*(x^4 + x^2) * dy^3 - 3*p*x^2 * dy^2 +
	- (2*p*(3*x^2 + 1) - 2*dp*(2*x^3 + x) + d2p*(x^4 + x^2)) * dy +
	- d2p*x^2 + 2*dp*x - 2*p # = 0
# Case: F0(x) = x^2;
d2y - x*(x^2+1)*dy^3 - 3*x*dy^2 # = 0
# Case: F0(x) = x;
x*d2y + (x^2+1)*dy^3 + 3*dy^2 + 2*dy # = 0


##########################
##########################

### Non-Linear Combination

### y*sin(y) + cos(y) = f(x)

### D =>
y*dy*cos(y) - df # = 0

### D2 =>
### ODE:
df*y*d2y - f*y*(dy)^3 + 2*df*(dy)^2 - d2f*y*dy # = 0

### Alternative Eq:
# y*dy = df / cos(y)
# needs also D3;
df*y*d2y - f*y*(dy)^3 + 2*df^3/(y^2 - (f-cos(y))^2) - d2f*y*dy # = 0


### Generalization:
### P(y)*sin(P(y)) + cos(P(y)) = f(x)
### D =>
# P*cos(P)*dP = df
# P*sin(P) = f - df / (P*dP)
### D2 =>
P*cos(P)*d2P + cos(P)*dP^2 - P*sin(P)*dP^2 - d2f # = 0
df/dP*d2P + df/P*dP - (f*dP - df/P)*dP - d2f # = 0
df*P*d2P - f*P*dP^3 + 2*df*dP^2 - d2f*P*dP # = 0

### Examples:

### Generalized Case:
### P(y) = y^2
2*df*y^2*(y*d2y + dy^2) - 8*f*y^5*dy^3 + 8*df*y^2*dy^2 - 2*d2f*y^3*dy # = 0
df*(y*d2y + dy^2) - 4*f*y^3*dy^3 + 4*df*dy^2 - d2f*y*dy # = 0
df*y*d2y - 4*f*y^3*dy^3 + 5*df*dy^2 - d2f*y*dy # = 0
### P(y) = y^2; f(x) = x + b;
y*d2y - 4*(x+b)*y^3*dy^3 + 5*dy^2 # = 0


### Simple Case:
# f = x + b, where b = constant;
y*d2y - (x+b)*y*(dy)^3 + 2*(dy)^2 # = 0

### Solution & Plot:
y = function(x, b) {
	# root
	y.f = function(x, v) x*sin(x) + cos(x) - v - b;
	dy.f = function(x, v) x*cos(x);
	x0.f = function(x) {
		xb = x + b
		x0 = if(xb >= 1 & xb <= pi/2) 1 else if(xb < 1) 3 else 7;
		return(x0);
	}
	y = sapply(x, function(x) newtonRaphson(y.f, x0.f(x), dfun=dy.f, v=x)[[1]])
	y = sapply(y, round0)
	return(y)
}
dy = function(x, b, y.x) {
	if(missing(y.x)) y.x = y(x, b=b);
	div = y.x * cos(y.x)
	dp = 1;
	dp = ifelse(div != 0, dp / div, -1); # may need correction
	return(dp)
}
d2y = function(x, b) {
	y.x = y(x, b=b)
	z = dy(x, b=b, y.x=y.x)
	div = y.x
	dp = - 2*z^2 + (x+b)*y.x*z^3;
	dp = ifelse(div != 0, dp / div, -1); # TODO: needs correction!
	return(dp)
}
### Plot:
b = -1/2
curve(y(x, b=b), from= -3, to = 3, ylim=c(-2,7))
# oscillating function with local minimum;
line.tan(c(-1, 1.3, 1.6, 1.8, 1.95), dx=3, p=y, dp=dy, b=b)
# spikes
curve(dy(x, b=b), from= -3, to = 3, add=T, col="green")
line.tan(c(-1, 1.3, 1.6, 1.8, 1.95), dx=3, p=dy, dp=d2y, b=b, col="orange")


#####################
### Generalized Case:
# P(y) = y^2
# f = x + b, where b = constant;
y*d2y - 4*(x+b)*y^3*dy^3 + 5*dy^2 # = 0

### Solution:
y = function(x, b, n=2) {
	# root
	y.f = function(x, v) {
		x = if(n == 1) x else x^n;
		x*sin(x) + cos(x) - v - b;
	}
	dy.f = function(x, v) {
		if(n == 1) return(x*cos(x));
		xn = if(n == 1) x else x^n;
		n*xn*xn/x*cos(xn);
	}
	x0.f = function(x) {
		xb = x + b
		x0 = if(xb >= 1 & xb <= pi/2) 3/4 else if(xb < 1) 1.5 else 2.5;
		return(x0);
	}
	y = sapply(x, function(x) newtonRaphson(y.f, x0.f(x), dfun=dy.f, v=x)[[1]])
	# if(n != 1) y = rootn(y, n)
	y = sapply(y, round0)
	return(y)
}
dy = function(x, b, y.x, n=2, dF=1) {
	if(missing(y.x)) y.x = y(x, b=b, n=n);
	dp = if(n == 1) dF else dF * y.x;
	y.x = if(n == 1) y.x else y.x^n;
	div = n * y.x * y.x * cos(y.x)
	dp = ifelse(div != 0, dp / div, 1E+3); # TODO: needs correction
	return(dp)
}
d2y = function(x, b, n=2) {
	y.x = y(x, b=b, n=n)
	dy.x = dy(x, b=b, y.x=y.x, n=n)
	div = y.x
	dp = 4*(x+b)*y.x^3*dy.x^3 - 5*dy.x^2;
	dp = ifelse(div != 0, dp / div, -1); # TODO: needs correction!
	return(dp)
}
###
b = -1/2
curve(y(x, b=b), from= -3, to = 3, ylim=c(-2,7))
# oscillating/spikes (ECG-like) function with local minima;
line.tan(c(seq(1.4, 2.2, by=0.2) - 0.05), dx=1.5, p=y, dp=dy, b=b)
# spikes
curve(dy(x, b=b), from= -3, to = 3, add=T, col="green")
line.tan(c(seq(1.4, 2.2, by=0.2) - 0.05), dx=1.4, p=dy, dp=d2y, b=b, col="orange")


##########################

##########################
### Composed Functions ###
##########################

### Trig(LOG(y))

### y * sin(k * log(y)) = f(x)

### y * sin(k * log(y)) = x^2
# k*dy * cos(k * log(y)) = (2*x*y - x^2*dy) / y;
x*y^2*d2y + x*y*dy^2 - (k^2+1)/2 * x^2*dy^3 - y^2*dy # = 0;
### Solution:
y = function(x, k, start=1) {
	# start: default = 1 functions for 0 < x <= sqrt(7) & k == 1;
	# start = 500: for x <= 9.5 (k == 1/2);
	# root
	y.f = function(x, v) {
		r = if(x == 0) -v^2 else if(x < 0) -Inf else x*sin(k * log(x)) - v^2;
		return(r)
	}
	dy.f = function(x, v) {
		if(x == 0) return(sign(k));
		if(x < 0) return(1E+3);
		x.log = k * log(x)
		k*cos(x.log) + sin(x.log);
	}
	# damped waves;
	x0.f = function(x) {
		x2 = x^2
		x0 = if(x2 >= 0 & x2 <= 7.46) start else if(x2 < 0 & x2 > -0.32) 1/4 else 500; # for k == 1;
		return(x0);
	}
	y = sapply(x, function(x) newtonRaphson(y.f, x0.f(x), dfun=dy.f, v=x)[[1]])
	y = sapply(y, round0)
	return(y)
}
dy = function(x, k, y.x) {
	if(missing(y.x)) y.x = y(x, k=k);
	y.log = ifelse(x == 0, -1, log(y.x))
	div = x^2 + k * y.x * cos(k * y.log)
	dp = 2*x*y.x;
	dp = ifelse(div != 0, dp / div, 0); # may need correction
	return(dp)
}
d2y = function(x, k) {
	y.x = y(x, k=k)
	dy.x = dy(x, k=k, y.x=y.x)
	div = - x*y.x^2
	dp  = x*y.x*dy.x^2 - (k^2+1)/2 * x^2*dy.x^3 - y.x^2*dy.x
	dp = ifelse(div != 0, dp / div, 0); # TODO: needs correction!
	return(dp)
}
### Plot:
k = 1;
curve(y(x, k=k), from= 0+1E-3, to = sqrt(7))
# oscillating function with local minimum;
line.tan(c(1/3, 1.1, 1.6, 1.8, 1.95), dx=3, p=y, dp=dy, k=k)
# log-like
curve(dy(x, k=k), add=T, col="green")
line.tan(c(1/3, 1.1, 1.7, 2, 2.5), dx=3, p=dy, dp=d2y, k=k, col="orange")

# separately D2:
curve(dy(x, k=k), from= 0+1E-3, to = sqrt(7), col="green")
line.tan(c(1/3, 1.1, 1.7, 2, 2.5), dx=3, p=dy, dp=d2y, k=k, col="orange")


### k == 1/2; # TODO: may need correcting y0();
k = 1/2;
curve(y(x, k=k), from= 0+1E-3, to = 2.7)
# for x > 2.7 up to 9.5: needs start=500;
# oscillating function with local minimum;
line.tan(c(1/3, 1.1, 1.6, 1.8, 1.95), dx=3, p=y, dp=dy, k=k)
# log-like
curve(dy(x, k=k), add=T, col="green")
line.tan(c(1/3, 1/2, 1, 2.2), dx=3, p=dy, dp=d2y, k=k, col="orange")


# separately D2:
curve(dy(x, k=k), from= 0+1E-3, to = 2.7, col="green")
sapply(c(1/3, 1/2, 1, 2.2), line.tan, dx=3, p=dy, dp=d2y, k=k, col="orange")


##########################

### Trig(EXP(y))

### sin(exp(y)) - y = F(x)

### D(y)
exp(y)*cos(exp(y))*dy - dy - df # = 0

### D2(y)
exp(y)*cos(exp(y))*d2y + exp(y)*cos(exp(y))*dy^2 +
	- exp(2*y)*sin(exp(y))*dy^2 - d2y - d2f # = 0
# =>
(dy + df)*d2y + (dy + df)*dy^2 +
	- (y + f)*exp(2*y)*dy^3 - dy*d2y - d2f*dy # = 0

### ODE:
df*d2y + dy^3 + df*dy^2 - d2f*dy - (y + f)*exp(2*y)*dy^3 # = 0

# TODO: check;

