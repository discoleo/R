#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL Log: Hidden Log
##
## draft v.0.1a


### History:

### draft v.0.1a:
# - moved to this file from:
#   DE.ODE.Log.R;


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#####################

##################
### Log: Other ###
##################

### y = (x + k1)^(x + k2) + F0(x)

# "a problem which stumped Photomath!"

### Check:
fe = expression(sin(x^2/3))[[1]];
ye = substitute(expression((x + k1)^(x + k2) + f), list(f=fe))[[2]]
x = sqrt(3); k1 = sqrt(2); k2 = 1/sqrt(5);
params = list(x=x, k1=k1, k2=k2);
#
y = eval(ye, params); dy = eval(D(ye, "x"), params);
d2y = eval(D(D(ye, "x"), "x"), params);
f0  = eval(fe, params); df0 = eval(D(fe, "x"), params);
d2f0 = eval(D(D(fe, "x"), "x"), params);


### Derivation:

### D(y)
dy - (log(x+k1) + (x+k2)/(x+k1))*(y - f0) - df0 # = 0
(x+k1)*dy # =
((x+k1)*log(x+k1) + (x+k2))*(y - f0) + (x+k1)*df0;
# =>
(x+k1)*log(x+k1) # =
(x+k1)*(dy - df0) / (y - f0) - (x + k2);

### D2(y)
(x+k1)*d2y + dy # =
((x+k1)*log(x+k1) + (x+k2))*(dy - df0) +
	+ (log(x+k1) + 2)*(y - f0) +
	+ (x+k1)*d2f0 + df0;
(x+k1)*(dy - df0)^2 / (y - f0) +
	+ (log(x+k1) + 2)*(y - f0) +
	+ (x+k1)*d2f0 + df0;
# =>
(x+k1)*(y - f0)*d2y + (y - f0)*dy # =
(x+k1)*(dy - df0)^2 +
	+ log(x+k1)*(y - f0)^2 + 2*(y - f0)^2 +
	+ (x+k1)*d2f0*(y - f0) + df0*(y - f0);

#
(x+k1)^2*(y - f0)*d2y + (x+k1)*(y - f0)*dy # =
(x+k1)^2*(dy - df0)^2 +
	+ (x+k1)*(y - f0)*(dy - df0) - (x + k2)*(y - f0)^2 + 2*(x+k1)*(y - f0)^2 +
	+ (x+k1)^2*d2f0*(y - f0) + (x+k1)*df0*(y - f0);

### ODE:
(x+k1)^2*(y - f0)*d2y # =
(x+k1)^2*(dy - df0)^2 + (x+k1)*(y - f0)^2 + (k1 - k2)*(y - f0)^2 +
	+ (x+k1)^2*d2f0*(y - f0);

### Examples:
### k1 == k2
(x+k1)*(y - f0)*d2y # =
(x+k1)*(dy - df0)^2 + (y - f0)^2 + (x+k1)*d2f0*(y - f0);


### Solution & Plot:
y = function(x, k=c(0,0), b=1, asList=FALSE) {
	f0 = eval.pol(x, b=b);
	y.x = (x + k[1])^(x + k[2]) + f0;
	if(asList) return(list(y=y.x, f0=f0));
	return(y.x)
}
dy = function(x, k=c(0,0), b=1, asList=FALSE) {
	y.all = y(x, k=k, b=b, asList=TRUE);
	y.x = y.all$y; f0 = y.all$f0;
	df0 = deriv.pol(x, b=b, dn=1);
	xk = x + k[1]; x.log = log(xk);
	dp = (xk*x.log + (x+k[2]))*(y.x - f0) + xk*df0
	div = xk;
	dp = ifelse(div != 0, dp/div, 0); # TODO
	if(asList) return(c(y.all, dy=dp, df0=df0));
	return(dp)
}
d2y = function(x, k=c(0,0), b=1) {
	y.all = dy(x, k=k, b=b, asList=TRUE);
	y.x = y.all$y; dy.x = y.all$dy;
	f0 = y.all$f0; df0 = y.all$df0;
	d2f0 = deriv.pol(x, b=b, dn=2);
	Tx = y.x - f0; xk = x + k[1]; Txk = xk*Tx;
	dTx = dy.x - df0;
	div = xk * Txk;
	d2p = xk^2*dTx^2 +
		+ (k[1] - k[2])*Tx^2 + Txk*Tx + xk*d2f0*Txk;
	d2p = ifelse(div != 0, d2p/div, 1); # TODO
	return(d2p)
}
### Plot:
k = c(0,0); b = c(0, 1)
px = (0:4)*4/7 + 0.1;
curve(y(x, k=k, b=b), from = 1E-3, to = 3, ylim=c(-2,8))
line.tan(px, dx=3, p=y, dp=dy, k=k, b=b)
#
curve(dy(x, k=k, b=b), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, b=b, col="orange")


### Ex 2:
k = c(0,0); b = c(-1, 1)
px = (0:4)*4/7 + 0.1;
curve(y(x, k=k, b=b), from = 1E-3, to = 3, ylim=c(-2,8))
line.tan(px, dx=3, p=y, dp=dy, k=k, b=b)
#
curve(dy(x, k=k, b=b), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, b=b, col="orange")


### Ex 3:
k = c(1,-1); b = c(-1, 1)
px = (0:4)*4/7 + 0.1; dx = c(1.5, 2, 2, 3, 3)
curve(y(x, k=k, b=b), from = 1E-3, to = 3, ylim=c(-0.5,8))
line.tan(px, dx=dx, p=y, dp=dy, k=k, b=b)
#
curve(dy(x, k=k, b=b), add=T, col="green")
line.tan(px, dx=dx, p=dy, dp=d2y, k=k, b=b, col="orange")


#########################

### y * (x + k1)^(x + k2) = F0(x)

# "another problem which stumped Photomath!" ;-)

### D(y)
f0*dy/y + (log(x+k1) + (x+k2)/(x+k1))*y*f0/y - df0 # = 0
f0*dy + (log(x+k1) + (x+k2)/(x+k1))*f0*y - df0*y # = 0
(x+k1)*f0*dy + ((x+k1)*log(x+k1) + (x+k2))*f0*y - (x+k1)*df0*y # = 0

# (x+k1)*log(x+k1)*f0*y = -(x+k1)*(f0*dy - df0*y) - (x + k2)*f0*y;

### D2(y)
(x+k1)*f0*d2y + f0*dy - df0*y - (x+k1)*d2f0*y +
	+ ((x+k1)*log(x+k1) + (x+k2))*(f0*dy + df0*y) +
	+ (log(x+k1) + 2)*f0*y # = 0
(x+k1)*f0*d2y + f0*dy - df0*y - (x+k1)*d2f0*y +
	- (x+k1)*(f0*dy - df0*y)*(dy/y + df0/f0) +
	+ log(x+k1)*f0*y + 2*f0*y # = 0
(x+k1)^2*f0^2*y*d2y +
	- (x+k1)^2*(f0*dy + df0*y)*(f0*dy - df0*y) +
	- (x+k1)^2*f0*d2f0*y^2 + (x+k1)*f0^2*y^2 + (k1 - k2)*f0^2*y^2 # = 0

### Examples:
### k1 == k2
(x+k1)*f0^2*y*d2y +
	- (x+k1)*(f0*dy + df0*y)*(f0*dy - df0*y) +
	- (x+k1)*f0*d2f0*y^2 + f0^2*y^2 # = 0


### Solution & Plot:
y = function(x, k=c(0,0), b=1, asList=FALSE) {
	f0 = eval.pol(x, b=b);
	div = (x + k[1])^(x + k[2]);
	y.x = ifelse(div != 0, f0/div, 0); # TODO
	if(asList) return(list(y=y.x, f0=f0));
	return(y.x)
}
dy = function(x, k=c(0,0), b=1, asList=FALSE) {
	y.all = y(x, k=k, b=b, asList=TRUE);
	y.x = y.all$y; f0 = y.all$f0;
	df0 = deriv.pol(x, b=b, dn=1);
	xk = x + k[1]; x.log = log(xk);
	dp = -((xk*x.log + x + k[2])*f0 - xk*df0)*y.x;
	div = xk * f0;
	dp = ifelse(div != 0, dp/div,
		dy(x + 1E-3, k=k, b=b)); # TODO
	if(asList) return(c(y.all, dy=dp, df0=df0));
	return(dp)
}
d2y = function(x, k=c(0,0), b=1) {
	y.all = dy(x, k=k, b=b, asList=TRUE);
	y.x = y.all$y; dy.x = y.all$dy;
	f0 = y.all$f0; df0 = y.all$df0;
	d2f0 = deriv.pol(x, b=b, dn=2);
	Tx  = f0*dy.x - df0*y.x;
	Txp = f0*dy.x + df0*y.x;
	xk = x + k[1]; xk2 = xk*xk; y2 = y.x*y.x;
	div = xk2*f0^2*y.x;
	d2p = xk2*Txp*Tx +
		+ xk2*f0*d2f0*y2 - xk*f0^2*y2 - (k[1] - k[2])*f0^2*y2;
	d2p = ifelse(div != 0, d2p/div,
		d2y(x + 1E-3, k=k, b=b)); # TODO
	return(d2p)
}
### Plot:
k = c(0,0); b = c(0, 1)
px = (0:4)*4/7 + 0.1;
curve(y(x, k=k, b=b), from = 1E-3, to = 3, ylim=c(-1.5,2))
line.tan(px, dx=3, p=y, dp=dy, k=k, b=b)
#
curve(dy(x, k=k, b=b), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, b=b, col="orange")


### Ex 2:
k = c(1, -2); b = c(0, 1)
px = (0:4)*4/7 + 0.1;
curve(y(x, k=k, b=b), from = 0, to = 3, ylim=c(-1.5, 3))
line.tan(px, dx=3, p=y, dp=dy, k=k, b=b)
#
curve(dy(x, k=k, b=b), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, k=k, b=b, col="orange")

