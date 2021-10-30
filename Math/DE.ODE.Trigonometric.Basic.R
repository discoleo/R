########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Trigonometric: Basic
###
### draft v.0.1f


### Trigonometric ODEs
### Basic Variants:
### Order 1 & 2 Linear


###############
### History ###
###############


### draft v.0.1f:
# - automatic generation of mixed Trig-Log ODEs:
#   y = P1(x) * sin(T0(x) + log(T1(x)));
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

### Helper functions

# library(pracma)
# needed for Lambert W;


# include: Polynomials.Helper.R;
# include: DE.ODE.Helper.R;
source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")

### Other

isNZ.pm = function(p) {
	is.data.frame(p) && (nrow(p) > 0);
}

# p1*sin(pT) + p2*cos(pT)
genODE.Trig.pm = function(p1, p2, pT, f0=NULL, print=FALSE, pDiv=NULL, div.by=NULL) {
	pC = list(p1, p2);
	pD = dp.trig.pm(pC, pT);
	# Linear System
	pR  = list(toPoly.pm("y"), toPoly.pm("dy"));
	d2f = NULL;
	if( ! is.null(f0)) {
		pR[[1]] = diff.pm(pR[[1]], f0);
		df0 = dp.pm(f0, xn="x");
		hasD = isNZ.pm(df0);
		if(hasD) {
			pR[[2]] = diff.pm(pR[[2]], df0);
			d2f = dp.pm(df0, xn="x");
			if( ! isNZ.pm(d2f)) d2f = NULL;
		}
	}
	# lapply(pR, print.pm);
	pR = solve.LD.pm(c(pC, pD[c("C1", "C2")]), pR);
	# D2 =>
	pD2 = dp.trig.pm(pD);
	pD2R = mult.pm(pD2$C1, pR$C1);
	pD2R = sum.pm(pD2R, mult.pm(pD2$C2, pR$C2));
	pD2y = pR$div; pD2y$d2y = 1;
	if(! is.null(d2f)) {
		pD2y = diff.pm(pD2y, mult.pm(pR$Div, d2f));
	}
	pD2R = diff.pm(pD2y, pD2R);
	if( ! is.null(pDiv)) pD2R = div.pm(pD2R, pDiv, by=div.by)$Rez;
	nms = c("d2y", "dy", "y", "x");
	idSort = nms %in% names(pD2R);
	nms  = nms[idSort];
	pD2R = sort.pm(pD2R, nms, sort.coeff=seq(10, length.out=length(nms)));
	nms = c(names(pD2R)[ ! names(pD2R) %in% nms], rev(nms));
	pD2R = pD2R[, nms];
	if(print) print.pm(pD2R, do.sort=FALSE, leading=NA);
	return(pD2R);
}

# p1*sin(pT) + p2*cos(pT)
# where pT = pT0 + log(pT1)
genODE.TrigLog.pm = function(p1, p2, pT, f0=NULL, print=FALSE, pDiv=NULL, div.by=NULL) {
	pC = list(p1, p2);
	pD = dp.trigLog.pm(pC, pT);
	# Linear System
	pR  = list(toPoly.pm("y"), toPoly.pm("dy"));
	d2f = NULL;
	if( ! is.null(f0)) {
		pR[[1]] = diff.pm(pR[[1]], f0);
		df0 = dp.pm(f0, xn="x");
		hasD = isNZ.pm(df0);
		if(hasD) {
			pR[[2]] = diff.pm(pR[[2]], df0);
			d2f = dp.pm(df0, xn="x");
			if( ! isNZ.pm(d2f)) d2f = NULL;
		}
	}
	# convert Fractions from: D(log(...))
	pD2y = mult.pm(pR[[2]], pD$Div);
	pR[[2]] = pD2y;
	pR = solve.LD.pm(c(pC, pD[c("C1", "C2")]), pR);
	# lapply(pR, print.data.frame);
	# D2 =>
	pD$Div = NULL; # reset DIV; (could be useful in the future)
	pD2  = dp.trigLog.pm(pD);
	pD2R = mult.pm(pD2$C1, pR$C1); # pD2RC1
	pD2R = sum.pm(pD2R, mult.pm(pD2$C2, pR$C2)); # pD2RC2
	# d2y:
	pD2y = dy.pm(pD2y, yn="y", xn="x");
	if(! is.null(d2f)) {
		pD2y = diff.pm(pD2y, d2f);
	}
	pD2y = mult.pm(pD2y, mult.pm(pD2$Div, pR$Div));
	pD2R = diff.pm(pD2y, pD2R);
	if( ! is.null(pDiv)) pD2R = div.pm(pD2R, pDiv, by=div.by)$Rez;
	nms = c("d2y", "dy", "y", "x");
	idSort = nms %in% names(pD2R);
	nms  = nms[idSort];
	pD2R = sort.pm(pD2R, nms, sort.coeff=seq(10, length.out=length(nms)));
	nms = c(names(pD2R)[ ! names(pD2R) %in% nms], rev(nms));
	pD2R = pD2R[, nms];
	if(print) print.pm(pD2R, do.sort=FALSE, leading=NA);
	return(pD2R);
}


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
p1 = toPoly.pm("a1"); # TODO: fix;
p2 = toPoly.pm("a2");
pDiv = toPoly.pm("a1^2 + a2^2");
# Note: a1 & a2 do NOT contribute;
pR = genODE.Trig.pm(p1, p2, pT, pDiv=pDiv, div.by="a1");
print.pm(pR, do.sort=FALSE, leading=NA)
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
print.pm(pR, do.sort=FALSE, leading=NA)
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
print.pm(pR, do.sort=FALSE, leading=NA)
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
print.pm(pR, do.sort=FALSE, leading=NA)

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

### Ex 1:
pT = toPoly.pm("x + b")
p1 = toPoly.pm("a1");
p2 = toPoly.pm("a2");
pDiv = toPoly.pm("a1^2 + a2^2");
# Note: a1 & a2 do NOT contribute;
pR = genODE.TrigLog.pm(p1, p2, pT, pDiv=pDiv, div.by="a1");
print.pm(pR, do.sort=FALSE, leading=NA)
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
print.pm(pR, do.sort=FALSE, leading=NA)
### ODE:
(x + b)*(k*x^2 + 2*k*b*x + x + k*b^2 + b)*d2y + (x + b)*dy +
	+ (k^3*x^3 + 3*k^2*x^2 + 3*b*k^3*x^2 + 3*k*x + 6*b*k^2*x + 3*b^2*k^3*x + 1 + 3*b*k + 3*b^2*k^2 + b^3*k^3)*y # = 0


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


### Solution:
y = function(x, a=c(1, 1), FUN.list, ...) {
	xlog = log(FUN.list[[1]](x, ...));
	r = a[1]*sin(xlog) + a[2]*cos(xlog)
	return(r)
}
dy = function(x, a=c(1, 1), FUN.list, ...) {
	div = FUN.list[[1]](x, ...)
	xlog = log(div);
	dp = (a[1]*cos(xlog) - a[2]*sin(xlog)) * FUN.list[[2]](x, ...)
	dp = ifelse(div != 0, dp / div, 0); # TODO: check;
	dp = round0(dp)
	return(dp)
}
d2y = function(x, a=c(1, 1), FUN.list, ...) {
	y.x  =  y(x, a=a, FUN.list=FUN.list, ...)
	dy.x = dy(x, a=a, FUN.list=FUN.list, ...)
	Px = FUN.list[[1]](x, ...)
	dP = FUN.list[[2]](x, ...); dPsq = dP^2;
	d2P = FUN.list[[3]](x, ...)
	#
	div = - dP * Px^2;
	dp  = Px*(dPsq - Px*d2P)*dy.x + dPsq*dP*y.x
	dp = ifelse(div != 0, dp / div, d2P/Px*(a[1]*cos(log(Px)) - a[2]*sin(log(Px)))); # TODO: check!
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
curve(y(x, b=b, FUN.list=p.list), from= -3, to = 3, ylim=c(-2, 2))
# oscillating function with local minima;
# slightly shifted: + 1/2
line.tan(c(-5:7 * 3/7 - 1/2), dx=2, p=y, dp=dy, FUN.list=p.list, b=b)
# also sinusoidal:
curve(dy(x, b=b, FUN.list=p.list), add=T, col="green")
line.tan(c(-5:7 * 3/7 - 1/2), dx=1.5, p=dy, dp=d2y, FUN.list=p.list, b=b, col="orange")

###
b = c(3, -1);
curve(y(x, b=b, FUN.list=p.list), from= -3, to = 3, ylim=c(-1, 1.5))
# possible local maximum;
# slightly shifted: ???
line.tan(c(-5:7 * 3/7 - 1/2), dx=2, p=y, dp=dy, FUN.list=p.list, b=b)
#
curve(dy(x, b=b, FUN.list=p.list), add=T, col="green")
line.tan(c(-5:7 * 3/7 - 1/2), dx=1.5, p=dy, dp=d2y, FUN.list=p.list, b=b, col="orange")

