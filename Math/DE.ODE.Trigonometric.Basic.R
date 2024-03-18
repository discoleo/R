########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Trigonometric: Basic
###
### draft v.0.2a


### Trigonometric ODEs
### Basic Variants:
### Order 1 & 2 Linear


###############
### History ###
###############


### draft v.0.2a:
# - moved Helper Functions (Generators) to new file:
#   Polynomials.Helper.ODE.R;
### draft v.0.1h - v.0.1h-fix2:
# - derived from:
#   y = P(x) * sin(T1(x))^2 + F0(x);
#   x^2*d2y - 3*x*dy + (16*x^4 + 3)*y - 8*x^5 = 0;
### draft v.0.1g:
# - [refactoring] use sort.dpm();
### draft v.0.1f - v.0.1f-ex2:
# - automatic generation of mixed Trig-Log ODEs:
#   y = P1(x) * sin(T0(x) + log(T1(x)));
# - more examples & checks; [v.0.1f-ex2]
### draft v.0.1e:
# - moved another Section with Basic Variants
#   from file: DE.ODE.Trigonometric.R;
# - type: y = P1(x)*sin(log(T(x)));
### draft v.0.1d - v.0.1d-check:
# - slight generalization:
#   y = ... + F0(x);
### draft v.0.1a - v.0.1c:
# - moved Section with Basic Variants
#   from file: DE.ODE.Trigonometric.R;



#########################

### Helper Functions

# library(pracma)
# needed for Lambert W;


# include: Polynomials.Helper.ODE.R;
# include: Polynomials.Helper.R; [automatically]
# include: DE.ODE.Helper.R;
source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")

### Other
# ...


#######################
#######################

#####################
### Linear Simple ###
###  (Polynomial) ###
#####################

### y = P1(x)*sin(T(x)) + P2(x)*cos(T(x)) + F0(x);
# - where T(x) is the same;


#################
### Automatic ###
#################

### Ex 1:
pT = toPoly.pm("x^3 + b*x")
p1 = toPoly.pm("a1");
p2 = toPoly.pm("a2");
pDiv = toPoly.pm("a1^2 + a2^2");
# Note: a1 & a2 do NOT contribute;
pR = genODE.Trig.pm(p1, p2, pT, pDiv=pDiv, div.by="a1");
print.dpm(pR, do.sort=FALSE)
### ODE:
(3*x^2 + b)*d2y - 6*x*dy + (27*x^6 + 27*b*x^4 + 9*b^2*x^2 + b^3)*y # = 0


### Ex 2:
pT = toPoly.pm("x^3 + b*x")
p1 = toPoly.pm("a1");
p2 = toPoly.pm("a2");
p0 = toPoly.pm("a3");
pDiv = toPoly.pm("a1^2 + a2^2");
# Note: a1 & a2 do NOT contribute;
pR = genODE.Trig.pm(p1, p2, pT, f0=p0, pDiv=pDiv, div.by="a1");
print.dpm(pR, do.sort=FALSE)
### ODE:
(3*x^2 + b)*d2y - 6*x*dy + (27*x^6 + 27*b*x^4 + 9*b^2*x^2 + b^3)*y +
	- 27*a3*x^6 - 27*b*a3*x^4 - 9*b^2*a3*x^2 - b^3*a3 # = 0


### Ex 3:
pT = toPoly.pm("x^2 + b*x")
p1 = toPoly.pm("a1*x");
p2 = toPoly.pm("a2");
p0 = toPoly.pm("a3*x");
#
pR = genODE.Trig.pm(p1, p2, pT, f0=p0, pDiv=NULL, div.by="a1");
print.dpm(pR, do.sort=FALSE)
### ODE:
(2*a1^2*x^3 + a1^2*b*x^2 + 2*a2^2*x - a1*a2 + b*a2^2)*d2y +
	- 2*(3*a1^2*x^2 + a1^2*b*x + a2^2)*dy +
	+ (8*a1^2*x^5 + 12*a1^2*b*x^4 + 6*a1^2*b^2*x^3 + 8*a2^2*x^3 + a1^2*b^3*x^2 - 12*a1*a2*x^2 +
		+ 12*b*a2^2*x^2 + 6*a1^2*x - 12*a1*b*a2*x + 6*b^2*a2^2*x + 2*a1^2*b - 3*a1*b^2*a2 + b^3*a2^2)*y +
	- 8*a1^2*a3*x^6 - 12*a1^2*b*a3*x^5 - 6*a1^2*b^2*a3*x^4 - 8*a2^2*a3*x^4 - a1^2*b^3*a3*x^3 + 12*a1*a2*a3*x^3 +
		- 12*b*a2^2*a3*x^3 + 12*a1*b*a2*a3*x^2 - 6*b^2*a2^2*a3*x^2 + 3*a1*b^2*a2*a3*x - b^3*a2^2*a3*x + 2*a2^2*a3


### Ex 4:
pT = toPoly.pm("x^3 - 2*x - 1")
p1 = toPoly.pm("x + 2")
p2 = toPoly.pm("x^2")
#
pR = genODE.Trig.pm(p1, p2, pT);
print.dpm(pR, do.sort=FALSE)

# TODO: check;

######################

################
### Explicit ###

### y = a1*sin(P(x)) + a2*cos(P(x))
# - both sin(P(x)) & cos(P(x)) are solutions;

### D =>
# dy = dP * (a1*cos(P(x)) - a2*sin(P(x)))
### D2 =>
# d2y = d2P*(a1*cos(P(x)) - a2*sin(P(x))) - dP^2 * y
### Solve for sin(P(x)), cos(P(x)) using y & dy;
# dP * sin(P(x)) = (-a2*dy + a1*dP*y) / (a1^2 + a2^2)
# dP * cos(P(x)) = ( a1*dy + a2*dP*y) / (a1^2 + a2^2)
### =>
# (a1^2 + a2^2)*dP*d2y = d2P*(a1*(a1*dy + a2*dP*y) + a2*(a2*dy - a1*dP*y)) - (a1^2 + a2^2)*dP^3 * y

### ODE:
dP*d2y - d2P*dy + dP^3 * y # = 0

### Examples:
### P(x) = x^2
x*d2y - dy + 4*x^2 * y # = 0
### P(x) = x^n => * x^(2-n)
x*d2y - (n-1)*dy + n^2*x^(2*n-1) * y # = 0
### P(x) = x^n + b*x
(n*x^(n-1)+b)*d2y - n*(n-1)*x^(n-2)*dy + (n*x^(n-1)+b)^3 * y # = 0
# [is detailed below]

### Example: P(x) = x^n + b*x;

### y = a1*sin(x^n + b*x) + a2*cos(x^n + b*x)
# dy = - (n*x^(n-1) + b) * (a2*sin(x^n + b*x) - a1*cos(x^n + b*x))
# sin(x^n + b*x) = - (a2*dy - a1*(n*x^(n-1) + b)*y) / (a1^2 + a2^2) / (n*x^(n-1) + b)
# cos(x^n + b*x) =   (a1*dy + a2*(n*x^(n-1) + b)*y) / (a1^2 + a2^2) / (n*x^(n-1) + b)
d2y = - n*(n-1)*x^(n-2) * (a2*sin(x^n + b*x) - a1*cos(x^n + b*x)) +
  - (n*x^(n-1) + b)^2 * (a1*sin(x^n + b*x) + a2*cos(x^n + b*x))
= - n*(n-1)*x^(n-2) * (a[2]*sin(x^n + b*x) - a[1]*cos(x^n + b*x)) - (n*x^(n-1) + b)^2 * y
(a1^2 + a2^2)*(n*x^(n-1) + b) * d2y +
  + n*(n-1)*x^(n-2)*( - a[2]*(a[2]*dy(x) - a[1]*(n*x^(n-1) + b)*y(x)) - a[1]*(a[1]*dy(x) + a[2]*(n*x^(n-1) + b)*y(x))) +
  + (a[1]^2 + a[2]^2)*(n*x^(n-1) + b)^3 * y(x)

### ODE:
(n*x^(n-1) + b) * d2y - n*(n-1)*x^(n-2)*dy(x) + (n*x^(n-1) + b)^3 * y(x) # = 0
### Example: n = 2
(2*x + b) * d2y(x, n=2) - 2*dy(x, n=2) + (2*x + b)^3 * y(x, n=2) # = 0
### Example: n = 3
(3*x^2 + b) * d2y(x, n=3) - 6*x*dy(x, n=3) + (3*x^2 + b)^3 * y(x, n=3) # = 0

### Solution & Plot:
y = function(x, b=1, n=2, a=c(1, 1), a0=0) {
	xn = x^n + b*x
	r = a[1]*sin(xn) + a[2]*cos(xn) + a0;
	return(r)
}
dy = function(x, b=1, n=2, a=c(1, 1), a0=0) {
	xn1 = x^(n-1)
	xn = (xn1 + b)*x
	r = - (n*xn1 + b) * (a[2]*sin(xn) - a[1]*cos(xn))
	r = round0(r)
	return(r)
}
d2y = function(x, b=1, n=2, a=c(1, 1), a0=0) {
	y.x = y(x, b=b, n=n, a=a, a0=a0)
	dy.x = dy(x, b=b, n=n, a=a, a0=a0)
	xn2 = n * x^(n-2)
	xnb = xn2 * x + b
	#
	div = xnb; xnb3 = xnb^3;
	dp  = (n-1)*xn2*dy.x - xnb3*y.x + xnb3*a0;
	dp = ifelse(div != 0, dp / div, n*(n-1)); # TODO: needs correction!
	return(dp)
}
### Plot:
n = 2; b = 1;
curve(y(x, b=b, n=n), from= -3, to = 3, ylim=c(-3, 3))
# oscillating function with local minima;
# slightly shifted: + 1/2;
# dy == 0 also for: x^2 + b*x - pi/4 = 0;
line.tan(c(-5:7 * 3/7 - 1/2), dx=1.5, p=y, dp=dy, b=b, n=n)
# also sinusoidal:
curve(dy(x, b=b, n=n), add=T, col="green")
line.tan(c(-5:7 * 3/7 - 1/2), dx=1.5, p=dy, dp=d2y, b=b, n=n, col="orange")


### Ex 2:
n = 2; b = -2;
curve(y(x, b=b, n=n), from= -3, to = 3, ylim=c(-4, 4))
# oscillating function with local minima;
# slightly shifted: + 1/2
line.tan(c(-5:7 * 3/7 - 1/2), dx=1.5, p=y, dp=dy, b=b, n=n)
# also sinusoidal:
curve(dy(x, b=b, n=n), add=T, col="green")
line.tan(c(-5:7 * 3/7 - 1/2), dx=1.5, p=dy, dp=d2y, b=b, n=n, col="orange")


### Ex 3:
n = 2; b = 1.5; a0=3; # a0 = -4;
px = c(-5:7 * 3/7 - 1/2);
curve(y(x, b=b, n=n, a0=a0), from= -3, to = 3, ylim=c(-4, 6))
# oscillating function with local minima;
# slightly shifted: + 1/2
line.tan(px, dx=1.5, p=y, dp=dy, b=b, n=n, a0=a0)
# also sinusoidal:
curve(dy(x, b=b, n=n, a0=a0), add=T, col="green")
line.tan(px, dx=1.4, p=dy, dp=d2y, b=b, n=n, a0=a0, col="orange")


#####################

### Basic Example: Simple / Power
# - generalization in the (next)/previous sections;

### y = sin(a*x^n)
dy = n*a*x^(n-1)*cos(a*x^n)

### ODE:
x*d2y - (n-1)*dy + n^2*a^2*x^(2*n-1)*y # = 0

### Test & Plot:
y = function(x, a=1, n=2) {
	r = sin(a*x^n)
	return(r)
}
dy = function(x, a=1, n=2) {
	xn1 = if(n == 2) x else x^(n-1);
	xn  = xn1 * x;
	dp = n*a*xn1*cos(a*xn)
	return(dp)
}
d2y = function(x, a=1, n=2) {
	y.x = y(x, a=a, n=n)
	dy.x = dy(x, a=a, n=n)
	div = x;
	dp = (n-1)*dy.x - n^2*a^2*x^(2*n-1)*y.x;
	dp = ifelse(div != 0, dp / div, n*(n-1)*a); # TODO: needs correction!
	return(dp)
}
### Plot
a = 1; n = 2;
curve(y(x, a=a, n=n), from= -3, to= 3, ylim=c(-2, 1.5))
# sinus wave;
line.tan(c((-5:5)/2.2), dx=3, p=y, dp=dy, a=a, n=n)
# wave
curve(dy(x, a=a, n=n), add=T, col="green")
line.tan(c((-5:5)/2.2), dx=1/5, p=dy, dp=d2y, a=a, n=n, col="orange")

### Test wave
curve(dy(x, a=a, n=n), from=-3, to=3, col="green")
line.tan(c((-5:5)/2.2), dx=1/5, p=dy, dp=d2y, a=a, n=n, col="orange")


### Example 2:
a = 1/3; n = 2;
curve(y(x, a=a, n=n), from= -3, to= 3, ylim=c(-1, 1.5))
# sinus wave;
line.tan(c((-5:5)/2.2), dx=3, p=y, dp=dy, a=a, n=n)
# wave
curve(dy(x, a=a, n=n), add=T, col="green")
line.tan(c((-5:5)/2.2), dx=1/5, p=dy, dp=d2y, a=a, n=n, col="orange")


#########################

### Simple Generalization

### y = sin(f(x))
# [but still basics: 1 trigonometric term]
# [not run]
dy = df * cos(f)
dy^2 = df^2 * cos(f)^2
dy^2 = df^2 * (1 - y^2)

### Variant 1: Non-linear ODE
dy^2 + df^2 * y^2 - df^2 # = 0

### Alternative: Linear ODE
dP*d2y - d2P*dy + dP^3 * y # = 0
# [see next Section]

### Examples:

### f = x + ln(x)
dy^2 + ((x+1)/x)^2 * y^2 - ((x+1)/x)^2 # = 0
x^2 * dy^2 + (x+1)^2 * y^2 - (x+1)^2 # = 0
### Solution & Plot:
y = function(x) {
	f = x + log(x)
	r = sin(f)
	return(r)
}
dy = function(x) {
	# Test formula
	y.x = y(x)
	x1 = (x+1)^2
	dp = x1 * (1 - y.x^2)
	div = x^2
	dp = ifelse(div != 0, dp / div, 1E+3); # TODO: check!
	f.sign = ifelse(x + log(x) <= pi/2, FALSE, TRUE)
	dp = sqrt(dp); dp[f.sign] = - dp[f.sign];
	return(dp)
}
### Plot:
curve(y(x), from= 0+1E-3, to=3, ylim=c(-1.5, 2))
line.tan(c((1:8)/3.2), dx=3, p=y, dp=dy)


### D2:
# dy^2 + df^2 * y^2 - df^2 # = 0
2*dy*d2y + 2*df^2 * y*dy + 2*df*d2f * y^2 - 2*df*d2f # = 0
dy*d2y + df^2 * y*dy + df*d2f * y^2 - df*d2f # = 0
# alternative:
dy*d2y + df^2 * y*dy + d2f * (df^2 - dy^2)/df - df*d2f # = 0
dy*d2y - d2f/df * dy^2 + df^2 * y*dy # = 0

### Examples:

### f = x + ln(x)
dy*d2y + (1/x^2)/((x+1)/x) * dy^2 + ((x+1)/x)^2 * y*dy # = 0
dy*d2y + 1/((x+1)*x) * dy^2 + ((x+1)/x)^2 * y*dy # = 0 # * x^2
x^2*dy*d2y + x/(x+1) * dy^2 + (x+1)^2 * y*dy # = 0
### Solution & Plot:
# reuses functions y(x) & dy(x) from above;
d2y = function(x) {
	y.x = y(x)
	dy.x = dy(x)
	dp = - (x/(x+1) * dy.x^2 + (x+1)^2 * y.x*dy.x);
	div = x^2 * dy.x;
	dp = ifelse(div != 0, dp / div, 1E+3); # TODO: check!
	return(dp)
}
### Plot:
curve(y(x), from= 0+1E-3, to=3, ylim=c(-1.5, 3))
line.tan(c((1:5)/2.2), dx=3, p=y, dp=dy)
# wave
curve(dy(x), add=T, col="green")
line.tan(c((1:5)/2.2), dx=1/5, p=dy, dp=d2y, col="orange")


######################
######################

######################
### Generalization ###
######################

### G(y) = P1(x) * sin(T(x)) + P2(x) * cos(T(x) + F(x)

### Simple: Higher Powers
### y^2 = x^p * sin(x^m)

### D =>
# 2*x*y*dy = p*y^2 + m*x^(p+m)*cos(x^m)

### D2 =>
2*x^2*y*d2y + 2*x^2*dy^2 - 2*(2*p + m - 1)*x*y*dy +
	+ m^2*x^(2*m)*y^2 + p*(p+m)*y^2 # = 0

### Special cases:

### m = 1
2*x^2*y*d2y + 2*x^2*dy^2 - 4*p*x*y*dy + x^2*y^2 + p*(p+1)*y^2 # = 0

### p = - m
2*x^2*y*d2y + 2*x^2*dy^2 + 2*(m + 1)*x*y*dy + m^2*x^(2*m)*y^2 # = 0
### p = - 1; m = 1;
2*x*y*d2y + 2*x*dy^2 + 4*y*dy + x*y^2 # = 0
### p = 1; m = -1;
2*x^4*y*d2y + 2*x^4*dy^2 + y^2 # = 0

### Plot:
y = function(x, pp=1, m=1, posRoot=TRUE) {
	# root
	y = x^pp * sin(x^m)
	y = sqrt(y)
	if( ! posRoot) y = -y;
	return(y)
}
dy = function(x, pp=1, m=1, posRoot=TRUE, y.x) {
	if(missing(y.x)) y.x = y(x, pp=pp, m=m, posRoot=posRoot);
	dp = pp*y.x^2 + m*x^(pp+m)*cos(x^m)
	div = 2*x*y.x;
	dp = ifelse(div != 0, dp / div, 0); # TODO: may need correction
	return(dp)
}
d2y = function(x, pp=1, m=1, posRoot=TRUE) {
	y.x  = y(x, pp=pp, m=m, posRoot=posRoot);
	dy.x = dy(x, pp=pp, m=m, posRoot=posRoot, y.x=y.x);
	dp = 2*x^2*dy.x^2 - 2*(2*pp + m - 1)*x*y.x*dy.x + m^2*x^(2*m)*y.x^2 + pp*(pp+m)*y.x^2;
	div = - 2*x^2*y.x;
	dp = ifelse(div != 0, dp / div, (2*pp + m - 1) / m); # TODO: needs correction!
	return(dp)
}
### Test

pp = 2; m = 1;
curve(y(x, pp=pp, m=m), from= 0, to = 2)
# global minimum;
line.tan(c((0:4)/2.2), dx=1.5, p=y, dp=dy, pp=pp, m=m)
# pseudo-sigmoidal
curve(dy(x, pp=pp, m=m), add=T, col="green")
line.tan(c((0:4)/2.2), dx=1/5, p=dy, dp=d2y, pp=pp, m=m, col="orange")


### m = 2
pp = 2; m = 2;
curve(y(x, pp=pp, m=m), from= -5/3, to = 5/3, ylim=c(-2, 2))
# global minimum;
line.tan(c((-3:3)/2.2), dx=1.5, p=y, dp=dy, pp=pp, m=m)
# pseudo-sigmoidal
curve(dy(x, pp=pp, m=m), add=T, col="green")
line.tan(c((-3:3)/2.2), dx=1/5, p=dy, dp=d2y, pp=pp, m=m, col="orange")


########################

########################
###  Linear Complex  ###
### (Non-Polynomial) ###
########################

###########
### LOG ###

### y = a1*sin(log(P(x))) + a2*cos(log(P(x)))


#################
### Automatic ###
#################

### genODE.TrigLog.pm():
# - Eq: y = p1(x) * sin(log(T(x))) + p2 * cos(log(T(x))) + F0(x);
# - pDiv used for simplification;

### Ex 1:
pT = toPoly.pm("x + b")
p1 = toPoly.pm("a1");
p2 = toPoly.pm("a2");
pDiv = toPoly.pm("a1^2 + a2^2");
# Note: a1 & a2 do NOT contribute;
pR = genODE.TrigLog.pm(p1, p2, pT, pDiv=pDiv, div.by="a1");
print.dpm(pR, do.sort=FALSE)
### ODE:
(x^2 + 2*b*x + b^2)*d2y + (x + b)*dy + y # = 0


### Ex 2:
pT1 = toPoly.pm("x + b")
pT0 = toPoly.pm("k*x")
p1 = toPoly.pm("a1");
p2 = toPoly.pm("a2");
pDiv = toPoly.pm("a1^2 + a2^2");
# Note: a1 & a2 do NOT contribute;
pR = genODE.TrigLog.pm(p1, p2, list(pT1, pT0), pDiv=pDiv, div.by="a1");
print.dpm(pR, do.sort=FALSE)
### ODE:
(x + b)*(k*x^2 + 2*k*b*x + x + k*b^2 + b)*d2y + (x + b)*dy +
	+ (k^3*x^3 + 3*k^2*x^2 + 3*b*k^3*x^2 + 3*k*x + 6*b*k^2*x + 3*b^2*k^3*x + 1 + 3*b*k + 3*b^2*k^2 + b^3*k^3)*y # = 0
(x + b)^2 * (k*x + k*b + 1)*d2y + (x + b)*dy +
	+ (k*x + k*b + 1)^3 * y # = 0


### Ex 3:
pT1 = toPoly.pm("x + b")
p1 = toPoly.pm("a1");
p2 = toPoly.pm("a2");
p0 = toPoly.pm("a3");
pDiv = toPoly.pm("a1^2 + a2^2");
# Note: a1 & a2 do NOT contribute;
pR = genODE.TrigLog.pm(p1, p2, pT1, f0=p0, pDiv=pDiv, div.by="a1");
print.dpm(pR, do.sort=FALSE)
### ODE:
(x + b)^2*d2y + (x + b)*dy + y - a3 # = 0


### Ex 4:
pT1 = toPoly.pm("x^2 + b1*x + b0")
p1 = toPoly.pm("a1");
p2 = toPoly.pm("a2");
p0 = toPoly.pm("a3");
pDiv = toPoly.pm("a1^2 + a2^2");
# Note: a1 & a2 do NOT contribute;
pR = genODE.TrigLog.pm(p1, p2, pT1, f0=p0, pDiv=pDiv, div.by="a1");
print.dpm(pR, do.sort=FALSE)
### ODE:
(2*x^5 + 5*b1*x^4 + 4*b0*x^3 + 4*b1^2*x^3 + 6*b0*b1*x^2 + b1^3*x^2 + 2*b0^2*x + 2*b0*b1^2*x + b0^2*b1)*d2y +
	+ (2*x^4 + 4*b1*x^3 + 3*b1^2*x^2 + b1^3*x - 2*b0^2 + b0*b1^2)*dy +
	+ (2*x + b1)^3*y - a3*(2*x + b1)^3 # = 0


### Ex 5: Simple
# x^p * y = sin(log(x)) + F0(x);

# D =>
x^p * dy + p*x^(p-1) * y - cos(log(x)) / x - df0 # = 0
x^(p+1) * dy + p*x^p * y - cos(log(x)) - x*df0 # = 0

# D2 =>
x^(p+1) * d2y + (p+1)*x^p * dy + p*x^p * dy + p^2*x^(p-1) * y +
	+ sin(log(x)) / x - x*d2f0 - df0 # = 0
x^(p+2) * d2y + (2*p + 1)*x^(p+1) * dy + (p^2 + 1)*x^p * y +
	- x^2*d2f0 - x*df0 - f0 # = 0

### Ex 5a: F0(x) = x^3
x^(p+2) * d2y + (2*p + 1)*x^(p+1) * dy + (p^2 + 1)*x^p * y - 10*x^3 # = 0

# Case: p = 1
genODE.TrigLog.pm(data.frame(x=-1, coeff=1), 0, pT=as.pm("x"), as.pm("x^2"))
# ODE:
x^2*d2y + 3*x*dy + 2*y - 10*x^2 # = 0


#############
### Explicit:

### D =>
# P(x) * dy = dP * (a1*cos(log(P(x))) - a2*sin(log(P(x))))
### D2 * P(x) =>
# P(x)^2 * d2y + P(x) * dP*dy = P(x)*d2P(x)*(a1*cos(log(P(x))) - a2*sin(log(P(x)))) - dP(x)^2 * y
### Solve for sin(log(P(x))), cos(log(P(x))) using y & dy;
# sin(log(P(x))) = (-a2*P*dy + a1*dP*y) / ((a1^2 + a2^2)*dP)
# cos(log(P(x))) = ( a1*P*dy + a2*dP*y) / ((a1^2 + a2^2)*dP)
(a1^2 + a2^2)*dP*(P(x)^2 * d2y + P(x)*dP*dy) =
	P(x)*d2P(x)*(a[1]*(a[1]*P(x)*dy + a[2]*dP(x)*y) + a[2]*(a[2]*P(x)*dy - a[1]*dP(x)*y)) - (a1^2 + a2^2)*dP(x)^3 * y
### Eq:
dP*P(x)^2 * d2y + P(x)*dP^2*dy - P(x)^2*d2P*dy + dP^3*y # = 0

### Examples:
### P(x) = x + b
(x+b)^2 * d2y + (x+b)*dy + y # = 0
### P(x) = x^2 + b1*x + b0
(2*x+b1)*(x^2+b1*x+b0)^2 * d2y + (x^2+b1*x+b0)*(2*x+b1)^2*dy - 2*(x^2+b1*x+b0)^2*dy + (2*x+b1)^3*y # = 0
### P(x) = x^n + b0
n*x^(n-1)*(x^n+b0)^2 * d2y + n^2*(x^n+b0)*x^(2*n-2)*dy - n*(n-1)*(x^n+b0)^2*x^(n-2)*dy + n^3*x^(3*n-3)*y # = 0 # * x^(2-n) / n
x*(x^n+b0)^2 * d2y + n*x^n*(x^n+b0)*dy - (n-1)*(x^n+b0)^2*dy + n^2*x^(2*n-1)*y # = 0
x*(x^n+b0)^2 * d2y + (x^(2*n) - (n-2)*b0*x^n - (n-1)*b0^2)*dy + n^2*x^(2*n-1)*y # = 0


### Solution & Plot:
y = function(x, a=c(1, 1, 0), FUN.list, ...) {
	if(length(a) < 3) a = c(a, 0);
	xlog = log(FUN.list[[1]](x, ...));
	r = a[1]*sin(xlog) + a[2]*cos(xlog) + a[3];
	return(r)
}
dy = function(x, a=c(1, 1, 0), FUN.list, ...) {
	div = FUN.list[[1]](x, ...)
	xlog = log(div);
	dp = (a[1]*cos(xlog) - a[2]*sin(xlog)) * FUN.list[[2]](x, ...)
	dp = ifelse(div != 0, dp / div, 0); # TODO: check;
	dp = round0(dp)
	return(dp)
}
d2y = function(x, a=c(1, 1, 0), FUN.list, ...) {
	if(length(a) < 3) a = c(a, 0);
	y.x  =  y(x, a=a, FUN.list=FUN.list, ...)
	dy.x = dy(x, a=a, FUN.list=FUN.list, ...)
	Px = FUN.list[[1]](x, ...)
	dP = FUN.list[[2]](x, ...); dPsq = dP^2;
	d2P = FUN.list[[3]](x, ...)
	#
	div = - dP * Px^2; # this is "-";
	dp  = Px*(dPsq - Px*d2P)*dy.x + dPsq*dP*y.x - a[3]*dP*dPsq;
	dp = ifelse(div != 0, dp / div,
		d2P/Px*(a[1]*cos(log(Px)) - a[2]*sin(log(Px)))); # TODO: check & a[3]!
	return(dp)
}
P = function(x, b=c(1, 1)) {
	x^2 + b[2]*x + b[1]
}
dP = function(x, b=c(1, 1)) {
	2*x + b[2]
}
d2P = function(x, b=c(1, 1)) {
	2
}
p.list = list(P, dP, d2P)

### Plot:
b = c(1, 1);
# slightly shifted: + 1/2
px = c(-5:7 * 3/7 - 1/2);
# oscillating function with local minima;
curve(y(x, b=b, FUN.list=p.list), from= -3, to = 3, ylim=c(-1.5, 2))
line.tan(px, dx=2, p=y, dp=dy, FUN.list=p.list, b=b)
# also sinusoidal:
curve(dy(x, b=b, FUN.list=p.list), add=T, col="green")
line.tan(px, dx=1.5, p=dy, dp=d2y, FUN.list=p.list, b=b, col="orange")


### Ex 2:
b = c(3, -1);
# slightly shifted: ???
px = c(-5:7 * 3/7 - 1/2);
# possible local/global maximum;
curve(y(x, b=b, FUN.list=p.list), from= -3, to = 3, ylim=c(-1, 1.5))
line.tan(px, dx=2, p=y, dp=dy, FUN.list=p.list, b=b)
#
curve(dy(x, b=b, FUN.list=p.list), add=T, col="green")
line.tan(px, dx=1.5, p=dy, dp=d2y, FUN.list=p.list, b=b, col="orange")


### Ex 3:
b = c(3, -1); a = c(1,1,3);
px = c(-5:7 * 3/7 - 1/2);
# possible local/global maximum;
curve(y(x, b=b, FUN.list=p.list, a=a), from= -3, to = 3, ylim=c(-1, 4))
line.tan(px, dx=2, p=y, dp=dy, FUN.list=p.list, b=b, a=a)
#
curve(dy(x, b=b, FUN.list=p.list, a=a), add=T, col="green")
line.tan(px, dx=1.5, p=dy, dp=d2y, FUN.list=p.list, b=b, a=a, col="orange")


#####################
#####################

#####################
### Higher Powers ###
#####################

### y = P(x) * sin(T1(x))^2 + F0(x);

# Note: sin(T(x))^2 = (1 - cos(2 * T(x))) / 2;

### D =>
dy - dp*sin(t1)^2 - 2*p*dt1*sin(t1)*cos(t1) - df0 # = 0 # * p =>
p*dy - dp*y - p^2*dt1*sin(2*t1) + dp*f0 - p*df0 # = 0

### D2 =>
p*d2y - d2p*y - (2*dp*dt1 + p*d2t1)*p*sin(2*t1) +
	- 2*p^2*dt1^2*cos(2*t1) + d2p*f0 - p*d2f0 # = 0 # * p =>
p^2*d2y - p*d2p*y - (2*dp*dt1 + p*d2t1)*(p*dy - dp*y + dp*f0 - p*df0)/dt1 +
	- 2*p^3*dt1^2*(1 - 2*sin(t1)^2) + p*d2p*f0 - p^2*d2f0 # = 0
p^2*dt1*d2y - p*d2p*dt1*y - (2*dp*dt1 + p*d2t1)*(p*dy - dp*y + dp*f0 - p*df0) +
	- 2*p^2*dt1^3*(p - 2*y + 2*f0) + p*d2p*dt1*f0 - p^2*dt1*d2f0 # = 0

### ODE:
p^2*dt1*d2y - p*(2*dp*dt1 + p*d2t1)*dy +
	+ 4*p^2*dt1^3*y - p*d2p*dt1*y + dp*(2*dp*dt1 + p*d2t1)*y +
	- (2*dp*dt1 + p*d2t1)*(dp*f0 - p*df0) +
	- 2*p^2*dt1^3*(p + 2*f0) + p*d2p*dt1*f0 - p^2*dt1*d2f0 # = 0

### Examples:
### p = x; dp = 1; f0 = 0;
x^2*dt1*d2y - x*(2*dt1 + x*d2t1)*dy +
	 + 4*x^2*dt1^3*y + (2*dt1 + x*d2t1)*y - 2*x^3*dt1^3 # = 0
# t1 = x^2; dt1 = 2*x;
x^2*d2y - 3*x*dy + (16*x^4 + 3)*y - 8*x^5 # = 0

# TODO: correct formulas;

### Solution & Plot:
y = function(x, p1, pT1, f0=NULL) {
	T1 = eval.FUN(x, pT1);
	p1 = eval.FUN(x, p1);
	yx = p1 * sin(T1)^2;
	if( ! is.null(f0)) yx = yx + eval.FUN(x, f0);
	return(yx)
}
dy = function(x, p1, pT1, f0=NULL) {
	# p*dy - dp*y - p^2*dt1*sin(2*t1) - p*df0
	t1.all = eval.DFUN(x, pT1);
	p1.all = eval.DFUN(x, p1);
	t1x = t1.all$f; dt1 = t1.all$df;
	p1x = p1.all$f; dp1 = p1.all$df;
	# y(x):
	yx  = p1x * sin(t1x)^2;
	if( ! is.null(f0)) yx = yx + eval.FUN(x, f0);
	dp  = dp1*yx + p1x^2*dt1*sin(2*t1x);
	if( ! is.null(f0)) {
		f0.all = eval.DFUN(x, f0);
		f0x = f0.all$f; df0 = f0.all$df;
		dp = dp - dp1*f0x + p1x*df0;
	}
	div = p1x;
	dp = ifelse(div != 0, dp / div, 0); # TODO: check;
	dp = round0(dp)
	return(dp)
}
d2y = function(x, p1, pT1, f0=NULL) {
	t1.all = eval.DFUN(x, dp.pm(pT1, xn="x"));
	p1.all = eval.DFUN(x, p1);
	dt1  = t1.all$f; d2t1 = t1.all$df;
	p    = p1.all$f; dp = p1.all$df;
	d2pf = dnp.pm(p1, n=2, xn="x");
	d2p = if(is.numeric(d2pf)) d2pf else eval.FUN(x, d2pf);
	yx  = y(x, p1=p1, pT1=pT1, f0=f0);
	dyx = dy(x, p1=p1, pT1=pT1, f0=f0);
	#
	div = - p^2 * dt1; # this is "-";
	ppr = (2*dp*dt1 + p*d2t1);
	dpR = - p*ppr*dyx +
		+ (4*p^2*dt1^3 - p*d2p*dt1 + dp*ppr)*yx +
		- 2*p^3*dt1^3;
	if( ! is.null(f0)) {
		f0.all = eval.DFUN(x, f0);
		d2f0   = eval.FUN(x, dnp.pm(f0, n=2, xn="x"));
		f0 = f0.all$f; df0 = f0.all$df;
		dpR = dpR - ppr*(dp*f0 - p*df0) +
			- 4*p^2*dt1^3*f0 + p*d2p*dt1*f0 - p^2*dt1*d2f0;
	}
	dpR = ifelse(div != 0, dpR / div, 0); # TODO: correct;
	return(dpR)
}

### Plot:
p1  = toPoly.pm("x");
pT1 = toPoly.pm("x^2 - 3*x - 3")
# slightly shifted:
px = c(-2.8, -2.65, -5:7 * 3/7 - 0.95/1.5, 3 - 0.95/1.5 + (1:3)/20);
# oscillating function;
curve(y(x, p1=p1, pT1=pT1), from= -3, to = 3)
line.tan(px, dx=1.5, p=y, dp=dy, p1=p1, pT1=pT1)
# also sinusoidal:
curve(dy(x, p1=p1, pT1=pT1), add=T, col="green")
line.tan(px, dx=1.5, p=dy, dp=d2y, p1=p1, pT1=pT1, col="orange")

### only D2:
curve(dy(x, p1=p1, pT1=pT1), from= -3, to = 3, n=256, col="green")
curve(y(x, p1=p1, pT1=pT1), add=T, col="gray")
line.tan(px, dx=1.4, p=dy, dp=d2y, p1=p1, pT1=pT1, col="orange")


# Test:
p1  = toPoly.pm("x")
pT1 = toPoly.pm("x^2")
x = seq(-3, 3, by=0.5)
yx = y(x, p1, pT1)
dyx = dy(x, p1, pT1)
d2yx = d2y(x, p1, pT1)
err = x^2*d2yx - 3*x*dyx + (16*x^4 + 3)*yx - 8*x^5;
round0(err)

pT1dbl = pT1; pT1dbl$coeff = 2*pT1dbl$coeff;
p1$coeff = - 1/2; # affects only f0!
f0 = p1; f0$coeff = - f0$coeff;
pDiv = toPoly.pm("x")
pR = genODE.Trig.pm(p1, NULL, pT1dbl, f0=f0, pDiv=pDiv, div.by="x", trig.order="cos")
print.dpm(pR, do.sort=FALSE)
### ODE:
x^2*d2y - 3*x*dy + 16*x^4*y + 3*y - 8*x^5 # = 0


### Ex 2:
p1  = toPoly.pm("x");
pT1 = toPoly.pm("x^2 - 3*x - 1");
f0  = toPoly.pm("x^3 + 1")
# slightly shifted:
px = c(-2.56, c(-5:0, 2, 4:7) * 3/7 - 0.95/1.5);
# oscillating function;
curve(y(x, p1=p1, pT1=pT1, f0=f0), from = -3, to = 3)
line.tan(px, dx=1.5, p=y, dp=dy, p1=p1, pT1=pT1, f0=f0)
# also sinusoidal:
curve(dy(x, p1=p1, pT1=pT1, f0=f0), add=T, col="green")
line.tan(px + 0.1, dx=1.4, p=dy, dp=d2y, p1=p1, pT1=pT1, f0=f0, col="orange")

# only D2:
curve(dy(x, p1=p1, pT1=pT1, f0=f0), from = -3, to = 3, col="green")
curve(y(x, p1=p1, pT1=pT1, f0=f0), add=T, col="grey")
line.tan(px + 0.1, dx=1.4, p=dy, dp=d2y, p1=p1, pT1=pT1, f0=f0, col="orange")

