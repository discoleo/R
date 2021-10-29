########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Trigonometric: Basic
###
### draft v.0.1a


### Trigonometric ODEs
### Basic Variants:
### Order 1 & 2 Linear


###############
### History ###
###############

### draft v.0.1a:
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

# p1*sin(pT) + p2*cos(pT)
genODE.Trig.pm = function(p1, p2, pT, print=FALSE, pDiv=NULL, div.by=NULL) {
	pC = list(p1, p2);
	pD = dp.trig.pm(pC, pT);
	# Linear System
	pR = list(toPoly.pm("y"), toPoly.pm("dy"));
	# lapply(pR, print.pm);
	pR = solve.LD.pm(c(pC, pD[c("C1", "C2")]), pR);
	# D2 =>
	pD2 = dp.trig.pm(pD);
	pD2R = mult.pm(pD2$C1, pR$C1);
	pD2R = sum.pm(pD2R, mult.pm(pD2$C2, pR$C2));
	pD2d2y = pR$div; pD2d2y$d2y = 1;
	pD2R = diff.pm(pD2d2y, pD2R);
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
p1 = toPoly.pm("a1"); p1$x = 0;
p2 = toPoly.pm("a2"); p2$x = 0;
pDiv = toPoly.pm("a1^2 + a2^2");
#
pR = genODE.Trig.pm(p1, p2, pT, pDiv=pDiv, div.by="a1");
print.pm(pR, do.sort=FALSE, leading=NA)


### Ex 2:
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

### Solution:
y = function(x, a=c(1, 1), b=1, n=2) {
	xn = x^n + b*x
	r = a[1]*sin(xn) + a[2]*cos(xn)
	return(r)
}
dy = function(x, a=c(1, 1), b=1, n=2) {
	xn1 = x^(n-1)
	xn = (xn1 + b)*x
	r = - (n*xn1 + b) * (a[2]*sin(xn) - a[1]*cos(xn))
	r = round0(r)
	return(r)
}
d2y = function(x, a=c(1, 1), b=1, n=2) {
	y.x = y(x, a=a, b=b, n=n)
	dy.x = dy(x, a=a, b=b, n=n)
	xn2 = n * x^(n-2)
	xnb = xn2 * x + b
	#
	div = xnb;
	dp  = (n-1)*xn2*dy.x - xnb^3*y.x
	dp = ifelse(div != 0, dp / div, n*(n-1)); # TODO: needs correction!
	return(dp)
}
### Plot:
n = 2; b = 1;
curve(y(x, b=b, n=n), from= -3, to = 3, ylim=c(-3, 3))
# oscillating function with local minima;
# slightly shifted: + 1/2;
# dy == 0 also for: x^2 + b*x - pi/4 = 0;
sapply(c(-5:7 * 3/7 - 1/2), line.tan, dx=1.5, p=y, dp=dy, b=b, n=n)
# also sinusoidal:
curve(dy(x, b=b, n=n), add=T, col="green")
sapply(c(-5:7 * 3/7 - 1/2), line.tan, dx=1.5, p=dy, dp=d2y, b=b, n=n, col="orange")


### Ex 2:
n = 2; b = -2;
curve(y(x, b=b, n=n), from= -3, to = 3, ylim=c(-4, 4))
# oscillating function with local minima;
# slightly shifted: + 1/2
sapply(c(-5:7 * 3/7 - 1/2), line.tan, dx=1.5, p=y, dp=dy, b=b, n=n)
# also sinusoidal:
curve(dy(x, b=b, n=n), add=T, col="green")
sapply(c(-5:7 * 3/7 - 1/2), line.tan, dx=1.5, p=dy, dp=d2y, b=b, n=n, col="orange")


