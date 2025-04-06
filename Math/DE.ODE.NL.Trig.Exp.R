########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Trig( Exp )
##
## draft v.0.1b


### NL ODEs: Trigonometric-Type
#   Subtype: Trig( EXP )

### Simple:
# P1(x) * exp(x) * y = sin(P2(x) * exp(x))
### Higher Poly:
# P1(x) * exp(Q(x)) * y = sin(P2(x) * exp(Q(x)))


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


########################

### P1(x)*exp(x)*y = sin(P2(x)*exp(x))

### Simple Homogeneous:
### exp(x)*y = sin(exp(x))

# Check:
e = expression(sin(exp(x)) * exp(-x));
x = sqrt(3);
params = list(x=x);
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


### D =>
dy + y - cos(exp(x)) # = 0

### D2 =>
d2y + dy + exp(x)*sin(exp(x)) # = 0
y*d2y + y*dy + sin(exp(x))^2 # = 0
y*d2y + y*dy + (1 - cos(exp(x))^2) # = 0
y*d2y + y*dy - (dy + y)^2 + 1 # = 0

### ODE:
y*d2y - dy^2 - y*dy - y^2 + 1 # = 0


### Extended Homogeneous:
### exp(x)*y = sin(p*exp(x))

# Check:
e = expression(sin(k * x^2 * exp(x)) * exp(-x));
x = sqrt(3); k = sqrt(5);
p = k*x^2; dp = 2*k*x; d2p = 2*k;
params = list(x=x, p=p, dp=dp, d2p=d2p, k=k);
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


### D =>
dy + y - (p+dp)*cos(p*exp(x)) # = 0

### D2 =>
d2y + dy + (p+dp)^2*exp(x)*sin(p*exp(x)) +
	- (dp+d2p)*cos(p*exp(x)) # = 0
(p+dp)*d2y + (p+dp)*dy + (p+dp)^3*exp(x)*sin(p*exp(x)) +
	- (dp+d2p)*(dy + y) # = 0
(p+dp)*d2y + (p-d2p)*dy - (dp+d2p)*y +
	+ (p+dp)^3*exp(x)*sin(p*exp(x)) # = 0
(p+dp)*y*d2y + (p-d2p)*y*dy - (dp+d2p)*y^2 +
	+ (p+dp)^3*sin(p*exp(x))^2 # = 0
(p+dp)*y*d2y + (p-d2p)*y*dy - (dp+d2p)*y^2 +
	+ (p+dp)^3*(1 - cos(p*exp(x))^2) # = 0
(p+dp)*y*d2y + (p-d2p)*y*dy - (dp+d2p)*y^2 +
	- (p+dp)*(dy + y)^2 + (p+dp)^3 # = 0

### ODE:
(p+dp)*y*d2y - (p+dp)*dy^2 - (p + 2*dp + d2p)*y*dy +
	- (p + 2*dp + d2p)*y^2 + (p+dp)^3 # = 0

### Special Case:
# p = k; dp = 0;
y*d2y - dy^2 - y*dy - y^2 + k^2 # = 0
# p = x + k - 1; dp = 1; p + dp = x + k;
(x+k)*y*d2y - (x+k)*dy^2 - (x+k+1)*y*dy - (x+k+1)*y^2 + (x+k)^3 # = 0


### Solution & Plot:
y = function(x, PFUN, n=1) {
	xe = exp(x^n);
	fx = eval.FUN(x, PFUN);
	val = sin(fx*xe) / xe;
	return(val)
}
dy = function(x, PFUN, n=1) {
	yx = y(x, PFUN=PFUN, n=n);
	xe = exp(x^n);
	fx = eval.FUN(x, PFUN);
	df  = dp.pm(PFUN, xn="x");
	dfx = eval.FUN(x, df);
	dp =  - yx + (fx+dfx)*cos(fx*exp(x));
	return(dp)
}
d2y = function(x, PFUN, n=1) {
	px  = eval.FUN(x, PFUN);
	dp  = dp.pm(PFUN, xn="x");
	dpx = eval.FUN(x, dp); ps = px + dpx;
	d2p = dp.pm(dp, xn="x");
	d2px = eval.FUN(x, d2p);
	dps = dpx + d2px; sall = ps + dps;
	yx  = y(x, PFUN=PFUN, n=n);
	dyx = dy(x, PFUN=PFUN, n=n);
	d2y = ps*dyx^2 + sall*yx*dyx + sall*yx^2 - ps^3;
	div = ps*yx;
	d2y = ifelse(div != 0, d2y / div, 1) # TODO
	return(d2y)
}
### Plot:
f = as.pm("x^2 - 3*x + 5")
px = c(0.2, 0.55, 0.7, (5:10)*3/17) + 1/13;
curve(y(x, PFUN=f), from = 0, to = 2.5, n=512, ylim = c(-2,1))
line.tan(px, dx=3, p=y, dp=dy, PFUN=f)
# anti-damped sinusoidal
curve(dy(x, PFUN=f), add=T, col="green")
line.tan(px, dx=1.6, p=dy, dp=d2y, PFUN=f, col="orange")

###
px = c(0.2, 0.55, 0.7, (5:6)*3/17) + 1/13; px = c(px, 1.5 + (1:10)/13);
xlim = c(0, 2.5); # xlim = c(1, 2.5)
curve(y(x, PFUN=f), from = xlim[1], to = xlim[2], n=512, ylim=c(-2, 4))
line.tan(px, dx=1.3, p=y, dp=dy, PFUN=f)
# anti-damped sinusoidal
curve(dy(x, PFUN=f), add=T, col="green")
line.tan(px, dx=1.6, p=dy, dp=d2y, PFUN=f, col="orange")


#########################

### Extended Homogeneous:
### exp(k1*x) * y = sin(p * exp(k1*x))

# Check:
e = expression(sin(k2 * x^2 * exp(k1*x)) * exp(-k1*x));
x = sqrt(3); k1 = sqrt(2); k2 = sqrt(5);
p = k2*x^2; dp = 2*k2*x; d2p = 2*k2;
params = list(x=x, p=p, dp=dp, d2p=d2p, k1=k1, k2=k2);
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
(k1*p + dp)*y*d2y - (k1*p + dp)*dy^2 +
	- (k1^2*p + 2*k1*dp + d2p)*y*dy +
	- k1*(k1^2*p + 2*k1*dp + d2p)*y^2 + (k1*p + dp)^3 # = 0

### Special Cases:

### k1*p + dp = x^n
x*y*d2y - x*dy^2 - (k1*x + n)*y*dy +
	- k1*(k1*x + n)*y^2 + x^(2*n+1) # = 0

# TODO


### Derivation:

### D =>
dy + k1*y - (k1*p + dp)*cos(p*exp(k1*x)) # = 0
# =>
(k1*p + dp)^2 * exp(2*k1*x) * y^2 + (dy + k1*y)^2 - (k1*p + dp)^2 # = 0

### D2 =>
d2y + k1*dy +
	+ (k1*p + dp)^2*sin(p*exp(k1*x))*exp(k1*x) +
	- (k1*dp + d2p)*cos(p*exp(k1*x)) # = 0
(k1*p + dp)*d2y + k1*(k1*p + dp)*dy +
	- (k1*dp + d2p)*(dy + k1*y) +
	+ (k1*p + dp)^3*exp(2*k1*x)*y # = 0
(k1*p + dp)*y*d2y + k1*(k1*p + dp)*y*dy +
	- (k1*dp + d2p)*(dy + k1*y)*y +
	- (k1*p + dp) * (dy + k1*y)^2 + (k1*p + dp)^3 # = 0


#######################
#######################

### Gen: Fully-Extended Homogeneous:
### p1*exp(x)*y = sin(p2*exp(x))

# Check:
e = expression(sin(k1 * x^2 * exp(x)) * exp(-x) / (x+b0));
x = sqrt(3); k1 = sqrt(5); b0 = sqrt(2);
p1 = x + b0; dp1 = 1; d2p1 = 0;
p2 = k1*x^2; dp2 = 2*k1*x; d2p2 = 2*k1;
parf = list(p1=p1, dp1=dp1, d2p=d2p);
parf = c(parf, list(p2=p2, dp2=dp2, d2p2=d2p2));
params = c(parf, list(x=x, k1=k1, b0=b0));
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### D =>
p1*dy + (p1+dp1)*y - (p2+dp2)*cos(p2*exp(x)) # = 0
#
(p1*dy + (p1+dp1)*y)^2 + p1^2*(p2+dp2)^2*exp(2*x)*y^2 - (p2+dp2)^2 # = 0


### D2 =>
p1*d2y + (p1 + 2*dp1)*dy + (dp1+d2p1)*y +
	- (dp2+d2p2)*cos(p2*exp(x)) +
	+ (p2+dp2)*(dp2 + p2)*exp(x)*sin(p2*exp(x)) # = 0
p1*(p2+dp2)*d2y +
	+ ((p1 + 2*dp1)*(p2+dp2) - p1*(dp2+d2p2))*dy +
	+ ((dp1+d2p1)*(p2+dp2) - (dp2+d2p2)*(p1+dp1))*y +
	+ (p2+dp2)^2*(dp2 + p2)*exp(x)*sin(p2*exp(x)) # = 0
p1*(p2+dp2)*d2y +
	+ ((p1 + 2*dp1)*(p2+dp2) - p1*(dp2+d2p2))*dy +
	+ ((dp1+d2p1)*(p2+dp2) - (dp2+d2p2)*(p1+dp1))*y +
	+ p1*(p2+dp2)^2*(p2+dp2)*exp(2*x)*y # = 0
p1^2*(p2+dp2)*y*d2y +
	+ p1*((p1 + 2*dp1)*(p2+dp2) - p1*(dp2+d2p2))*y*dy +
	+ p1*((dp1+d2p1)*(p2+dp2) - (dp2+d2p2)*(p1+dp1))*y^2 +
	- (p2+dp2) * ((p1*dy + (p1+dp1)*y)^2) +
	+ (p2+dp2)^3 # = 0

# TODO: simplify;

