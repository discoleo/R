#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs: Polynomial Types
## Cardano-Polynomial-System
##
## draft v.0.2a


### Cardano Polynomials
# & corresponding System;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")

#########################


###########################
### Cardano-Polynomials ###
###########################

### System:
# p^n + q^n = 2*f(x)
# p*q = h(x)

### Full P3: y = p + q;
# - see section in file: DE.ODE.NL.Polynomial.R;

### P6-partial:
# => q = h / p
# =>
# n*p^(n-1)*dp + n*q^(n-1)*dq = 2 * df
# p^(n-1)*dp + q^(n-1)*dq = 2/n * df
# p^(n-1)*dp + (h/p)^(n-1) * (dh/p - h*dp/p^2) = 2/n * df
# p^(2*n)*dp + h^(n-1) * (p*dh - h*dp) - 2/n * p^(n+1) * df = 0
#
# p^(2*n)*dp - h^n*dp + h^(n-1)*p*dh - 2/n * p^(n+1) * df = 0

### Solutions
# p = (f + sqrt(f^2 - h^n))^(1/n)
# q = (f - sqrt(f^2 - h^n))^(1/n)
# Note:
# - 2 basic solutions are possible;
# - it is possible to rotate these solutions using the roots of unity;


### Transformations:

### Order 3: [redundant]
# same as regular solution;
p^2*dp + q^2*dq = 2/3*df # * p^2
p^4*dp + p^2*q^2*dq - 2/3*df*p^2 # = 0
p^4*dp + h^2*dq - 2/3*df*p^2 # = 0
p^4*dp + h^2*(dh - q*dp)/p - 2/3*df*p^2 # = 0 # * p
p^5*dp + h^2*(dh - q*dp) - 2/3*df*p^3 # = 0 # *p
p^6*dp + h^2*(dh*p - h*dp) - 2/3*df*p^4 # = 0

### Variant:
p^2*dp + q^2*dq = 2/3*df
p^4*dp^2 + q^4*dq^2 + 2*p^2*q^2*dp*dq = 4/9*df^2
p^4*dp^2 + q^4*((dh - q*dp)/p)^2 + 2*h^2*dp*(dh - q*dp)/p = 4/9*df^2 # *p^2
p^6*dp^2 + q^4*(dh - q*dp)^2 + 2*h^2*dp*(p*dh - p*q*dp) = 4/9*df^2
p^6*dp^2 + q^4*(dh - q*dp)^2 + 2*h^2*dp*(p*dh - h*dp) = 4/9*df^2
p^6*dp^2 + q^4*(dh - q*dp)^2 + 2*h^2*dh*p*dp - 2*h^3*dp^2 = 4/9*df^2
# TODO: ...


################
### Examples ###

#########
### n = 2
p^4*dp - h^2*dp + h*p*dh - p^3 * df = 0

###
# h(x) = x
# f(x) = 1
p^4*dp - x^2*dp + x*p = 0
# p = sqrt(f + sqrt(f^2 - h^2))
p = function(x) {
	sqrt(1 + sqrt(1 - x^2))
}
dp = function(x) {
	p.x = p(x)
	div = (p.x^4 - x^2)
	dp = - x*p.x
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=-1, to=1)
sapply(c((0:5)/6), line.tan, dx=2, p=p, dp=dp)


###
# h(x) = x
# f(x) = 2*x
p^4*dp - x^2*dp - 2*p^3 + x*p = 0
# p = sqrt(f + sqrt(f^2 - h^2))
p = function(x) {
	sqrt(2*x + x*sqrt(3))
}
dp = function(x) {
	p.x = p(x)
	div = (p.x^4 - x^2)
	dp = 2*p.x^3 - x*p.x
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=0, to=3)
sapply(c((0:5)/4.5), line.tan, dx=2, p=p, dp=dp)

#########

#########
### n = 3
p^6*dp - h^3*dp - 2/3 * p^4 * df + h^2*p*dh = 0
# where df, dh = given;

###
# h(x) = x
# f(x) = 1
p^6*dp - x^3*dp + x^2*p = 0
# Solution & Plot:
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (1 + sqrt(1 - x^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
dp = function(x) {
	div = (p(x)^6 - x^3)
	dp = - x^2*p(x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=-3, to=1)
sapply(c((0:5)/6), line.tan, dx=2, p=p, dp=dp)


###
# h(x) = x
# f(x) = 3*x
p^6*dp - x^3*dp - 2*p^4 + x^2*p = 0
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (3*x + sqrt(9*x^2 - x^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
dp = function(x) {
	div = (p(x)^6 - x^3)
	dp = 2*p(x)^4 - x^2*p(x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=-3, to=9)
sapply(c(0.001, 2*(1:3), -1), line.tan, dx=3)

###
# h(x) = x
# f(x) = x^3 + 3*x
p^6*dp - x^3*dp - 2*(x^2+1)*p^4 + x^2*p = 0
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (x^3 + 3*x + sqrt((x^3 + 3*x)^2 - x^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
dp = function(x) {
	div = (p(x)^6 - x^3)
	dp = 2*(x^2+1)*p(x)^4 - x^2*p(x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=-3, to=9)
sapply(c(0.001, 2*(1:3), -1), line.tan, dx=2)


###
# h(x) = x
# f(x) = 3/2 * log(x)
p^6*dp - h^3*dp - 2/3 * p^4 * df + h^2*p*dh = 0
p^6*dp - x^3*dp - 1/x*p^4 + x^2*p = 0
x*p^6*dp - x^4*dp - p^4 + x^3*p = 0
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (3/2 * log(x) + sqrt((3/2 * log(x))^2 - x^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
dp = function(x) {
	div = x*(p(x)^6 - x^3)
	dp = p(x)^4 - x^3*p(x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=0, to=0.8)
sapply(c(0.001, (1:3)/6), line.tan, dx=2)


###
# h(x) = 1/(x^2 + 1)
# f(x) = 3/2 * log(x^2 + 1)
p^6*dp - h^3*dp - 2/3 * p^4 * df + h^2*p*dh = 0
(x^2 + 1)^4 * p^6 * dp - (x^2 + 1)*dp - 2*x*(x^2 + 1)^3 * p^4 - 2*x*p = 0
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (3/2 * log(x^2+1) + sqrt(9/4 * log(x^2+1)^2 - 1/(x^2+1)^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
dp = function(x) {
	x.mult = (x^2+1)^3
	p.x = p(x)
	div = (x^2+1)*(x.mult*p.x^6 - 1)
	dp = 2*x*(x.mult*p.x^4 - p.x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=-5, to=5)
sapply(c(-(1:4), 1:4), line.tan, dx=3)




####################
####################
####################


####################

### F(D1): dz + f(z)

### y^3 - 3*c*y - 2*f = 0
# [not run]
c*y*d2y + f*d2y + c*dy^2 + 2/3*df*dy - 1/2*d2c*y^2 - 1/3*d2f*y # = 0

### Special Cases:
# c = constant; f = x^3 - 1
x^2*(x^3-1)*d2y - c^2*dy^3 - c*x^2*dy^2 + 2*x*dy # = 0
# z = dy; dz = d2y;
x^2*(x^3-1)*dz - c^2*z^3 - c*x^2*z^2 + 2*x*z # = 0

### Solution & Plot:
# y = (f + sqrt(f^2 - c^3))^(1/3) + (f - sqrt(f^2 - c^3))^(1/3)
y = function(x, h, n=3) {
	f.x = x^3 - 1
	det = sqrt(f.x^2 - h^3 + 0i)
	r1 = round0(rootn(f.x + det, n=n))
	r2 = round0(rootn(f.x - det, n=n))
	return( round0(r1 + r2) )
}
dy = function(x, h, y.x) {
	if(missing(y.x)) y.x = y(x, h);
	div = (h*y.x + x^3 - 1)
	dp = x^2*y.x;
	dp = ifelse(div != 0, dp / div, Inf); # TODO
	return(dp)
}
d2y = function(x, h) {
	dy.x = dy(x, h);
	div = x^2*(x^3 - 1)
	dp = h^2*dy.x^3 + h*x^2*dy.x^2 - 2*x*dy.x;
	dp = ifelse(div != 0, dp / div, Inf); # TODO
	return(dp)
}
### Plot:
c = 1;
# a small local minimum
curve(y(x, c), from=-3, to=3)
line.tan(c(-2, c(-(2:1), 1:2)/1.7, 2), dx=3, p=y, dp=dy, h=c)
# D(y)
curve(dy(x, c), from=-3, to=3, add=T, col="green")
line.tan(c(-2, c(-(2:1), 1:2)/1.7, 2), dx=3, p=dy, dp=d2y, h=c, col="orange")


#######################

### (x^n - 1)*dy = f(x);
# [simple Model]
# x*(x^n - 1)*d2y + n*dy = x*df - n*f;
x*(x^n - 1)*d2y + b*(x^n - 1)*y*dy + n*dy - b*f*y - x*df + n*f # = 0
# where b = a constant or a function;
# e.g. b = x; f = b0 (constant);
x*(x^n - 1)*d2y + x*(x^n - 1)*y*dy + n*dy - b0*x*y + n*b0 # = 0

### Solution & Plot:
y = function(x, n=3, b0=1) {
	f.x = x^n - 1
	y.f = function(x) b0/(x^n - 1)
	y = sapply(x, function(x) integrate(y.f, lower=1 + 1E-5, upper=x)$value);
	return( y )
}
dy = function(x, n, b0=1) {
	div = x^n - 1
	dp = b0; # f = b0;
	dp = ifelse(div != 0, dp / div, Inf); # TODO
	return(dp)
}
d2y = function(x, n, b0=1) {
	# uses simple equation;
	dy.x = dy(x, n=n, b0=b0);
	div = x*(x^n - 1);
	dp = - n * (dy.x + b0);
	dp = ifelse(div != 0, dp / div, Inf); # TODO
	return(dp)
}
d2ymixt = function(x, n, b0=1) {
	# checking the mixt equation;
	y.x = y(x, n, b0)
	dy.x = dy(x, n=n, b0=b0);
	xn = x^n - 1
	dp = - (x*xn*y.x*dy.x + n*dy.x - b0*x*y.x + n*b0);
	div = x*xn;
	dp = ifelse(div != 0, dp / div, Inf); # TODO
}
#
c = 1; n = 5;
#
curve(y(x, n, c), from=1 + 1E-5, to=3)
line.tan(c(1.1, 2:5/1.7), dx=3, p=y, dp=dy, n=n, b0=c)
# D(y)
curve(dy(x, n, c), add=T, col="green")
# line.tan(c(1.1, 2:5/1.7), dx=3, p=dy, dp=d2y, n=n, b0=c, col="orange")
line.tan(c(1.1, 2:5/1.7), dx=3, p=dy, dp=d2ymixt, n=n, b0=c, col="orange")


#########################

### Integration by parts:

### Base ODE:
h*y*dy + f*dy - 1/2*dh*y^2 - 1/3*df*y # = 0

### Case: h = c (constant);
c*y*dy + f*dy - 1/3*df*y # = 0

### Ex 1:
# f = b1*x + b0;
c*y*dy + (b1*x+b0)*dy - b1/3 * y # = 0
# I() =>
c/2*y^2 + (b1*x+b0)*y - 4/3*b1*I(y) # = 0
# z = I(y); dz = y;
# c/2*dz^2 + (3*x+b0)*dz - 4/3*b1*z = 0;
### Variant:
c/3*y^3 + 1/2*(b1*x+b0)*y^2 - 2/3*b1*I(y^2) # = 0

### Solution & Plot:
# y = (f + sqrt(f^2 - h^3))^(1/3) + (f - sqrt(f^2 - h^3))^(1/3)
y = function(x, b, h, n=3) {
	d = (b[1]*x + b[2]);
	det = sqrt(d^2 - h^n + 0i)
	r1 = (d + det)
	r2 = (d - det)
	r = round0(rootn(r1, n=n) + rootn(r2, n=n))
	return(r)
}
Iy = function(x, b, h, n=3) {
	yx = y(x, b=b, h=h, n=n)
	yI = (h/2*yx^2 + (b[1]*x+b[2])*yx) / 4 * 3/b[1]
	return(yI)
}
Iy2 = function(x, b, h, n=3) {
	yx = y(x, b=b, h=h, n=n)
	yI = (h/3*yx^3 + 1/2*(b[1]*x+b[2])*yx^2) / 2 * 3/b[1]
	return(yI)
}
dI2 = function(x, b, h, n=3) {
	yx = y(x, b=b, h=h, n=n)
	return(yx^2)
}
### Plot
b = c(3, 1); h = 1;
#
curve(Iy(x, b=b, h=h), from=-3, to=3)
line.tan(c(-2:2), dx=3, p=Iy, dp=y, b=b, h=h)

#
curve(Iy2(x, b=b, h=h), from=-3, to=3)
line.tan(c(-2:2), dx=3, p=Iy2, dp=dI2, b=b, h=h)

###
b = c(3, -1/3); h = 1;
#
curve(Iy(x, b=b, h=h), from=-3, to=3)
line.tan(c(-2:2), dx=3, p=Iy, dp=y, b=b, h=h)

###
b = c(3, -2); h = 2;
#
curve(Iy(x, b=b, h=h), from=-3, to=3)
line.tan(c(-2:2), dx=3, p=Iy, dp=y, b=b, h=h)

###
b = c(2, 2); h = -3;
#
curve(Iy(x, b=b, h=h), from=-3, to=3)
line.tan(c(-2:2), dx=3, p=Iy, dp=y, b=b, h=h)

###
b = c(-2, 2); h = -3;
#
curve(Iy(x, b=b, h=h), from=-3, to=3)
line.tan(c(-2:2), dx=3, p=Iy, dp=y, b=b, h=h)


########################
########################

##################
### y-Radicals ###
##################

### Initial "Polynomial"
# y - G(x)*y^(1/n) = F(x)

### D =>
dy - 1/n*g*y^(1/n-1)*dy - dg*y^(1/n) - df # = 0
g*y*dy - 1/n*g*(y - f)*dy - dg*y*(y - f) - g*df*y # = 0
### ODE:
(n-1)*g*y*dy + g*f*dy - n*dg*y^2 + n*(dg*f - g*df)*y # = 0

### Example:
# g = x^2;
# f = x^2;
(n-1)*x*y*dy + x^3*dy - 2*n*y^2 # = 0

### Solution & Plot:
y = function(x, n=3) {
	coeff0 = c(1, rep(0, n-2));
	r = sapply(x, function(x) {
		coeff = c(coeff0, -x^2, -x^2);
		r = roots(coeff);
		r = r[Im(r) == 0][1];
		return(r);
	})
	r = r^n;
	return(r)
}
dy = function(x, n=3) {
	# (n-1)*x*y*dy + x^3*dy - 2*n*y^2 # = 0
	yx = y(x, n=n); x3 = x^3;
	dy = 2*n*yx^2;
	div = (n-1)*x*yx + x3;
	dy = dy / div;
	return(dy)
}
### Plot
n = 3;
px = -3:3 * 3/5;
curve(y(x, n=n), from=-2, to=2)
line.tan(px, p=y, dp=dy, n=n)

