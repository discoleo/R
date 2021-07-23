########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### DE Systems: Polynomial
###
### draft v.0.1e


#############
### Types ###
#############

### Simple:
# - TODO;
### Hetero-Symmetric:
# - Simple, derived from:
#   S.1.) y1^n + c1*y2 = R;
#   S.2.) (y1*y2)^n + c1*y2 = R;
### Others:
# TODO


###############
### History ###
###############


### draft v.0.1e:
# - plot/test for system derived from:
#   (y1*y2)^n + c1*y2 = R;
### draft v.0.1d:
# - improved document structure;
### draft v.0.1a - v.0.1c:
# - system:
#   n*(R - c0*y2)*dy1 + c0*y1*dy2 = y1*dR; [plot in v.0.1c]
# - system:
#   n*(R1 - c0*g)*(g*df + f*dg) + c0*f*g*dg = f*g*dR1; [v.0.1b]


#########################

### Helper functions

library(pracma)
# needed for Lambert W;


# include: DE.ODE.Helper.R;
source("DE.ODE.Helper.R")


### Real roots
isRe.f = function(x) {
	lapply(x, function(x) (Im(x) == 0));
}
Re.f = function(x, isRe=NULL) {
	if( ! is.null(isRe) && is.list(isRe) && ! is.list(x)) {
		x = lapply(seq_along(isRe), function(id) x);
	}
	if(is.null(isRe)) isRe = isRe.f(x);
	x = lapply(seq_along(x), function(id) Re(x[[id]][isRe[[id]]]));
	return(x);
}
range.c = function(x, isRe=NULL) {
	if(is.null(isRe)) isRe = isRe.f(x);
	x1 = Re(x[[1]][isRe[[1]]]);
	x2 = Re(x[[2]][isRe[[2]]]);
	xmax = max(x1, x2)
	xmin = min(x1, x2)
	return(list(rg=c(xmin, xmax), isRe=isRe));
}

########################
########################
########################

########################
### Hetero-Symmetric ###
########################

### Derived from:
# y1^n + c1*y2 = R1
# y2^n + c1*y1 = R2
# where (y1, y2) = functions in x;
# R  = R(x) is a given parameter function;
# c1 = numeric constant;
# Note:
# - the original polynomial system is much easier to solve
#   when R1 = R2;

### D & Mult(y1*...) =>
n*(R1 - c1*y2)*dy1 + c1*y1*dy2 - y1*dR1 # = 0
n*(R2 - c1*y1)*dy2 + c1*y2*dy1 - y2*dR2 # = 0

### Examples:

# R = x + b0; [R1 = R2 = R]
n*(x + b0 - c1*y2)*dy1 + c1*y1*dy2 - y1 # = 0
n*(x + b0 - c1*y1)*dy2 + c1*y2*dy1 - y2 # = 0
# =>
c1*y1*dy2 - n*c1*y2*dy1 + n*(x + b0)*dy1 - y1 # = 0
c1*y2*dy1 - n*c1*y1*dy2 + n*(x + b0)*dy2 - y2 # = 0

n*(x + b0 - c1*y2)*dy1 + c1*y1*dy2 - y1 # = 0
c1*y2*dy1 + n*(x + b0 - c1*y1)*dy2 - y2 # = 0

### Solution:
y.f = function(x, b0=1, c1=1, n=3, asMatrix=FALSE) {
	### for n = 3!
	if(n != 3) stop("Implemented only for n = 3!")
	# y1*y2 = S^2 - c1;
	# S^3 - 2*c1*S + R = 0;
	S = sapply(x, function(x) roots(c(1,0, -2*c1, x + b0)));
	# x = rep(x, each=3); # NOT needed;
	ypr = S^2 - c1;
	yd  = sqrt(S^2 - 4*ypr + 0i);
	y1 = (S + yd)/2; y2 = (S - yd)/2;
	if( ! asMatrix) {
		sol = list(y1=as.vector(y1), y2=as.vector(y2));
	} else {
		sol = list(y1=y1, y2=y2); # ??? also Matrix ???
	}
	return(sol);
}
dy.f = function(x, b0=1, c1=1, n=3) {
	y.all = y.f(x, b0=b0, c1=c1, n=n);
	y1 = y.all[[1]]; y2 = y.all[[2]];
	x = rep(x, each=(n^2 - n)/2); # only y1 != y2
	nxb = n*(x + b0);
	div = n^2*(x + b0 - c1*y2)*(x + b0 - c1*y1) - c1^2*y1*y2;
	dy1 = (nxb*y1 - n*c1*y1^2 - c1*y1*y2) / div;
	dy2 = (nxb*y2 - n*c1*y2^2 - c1*y1*y2) / div;
	return(list(dy1=dy1, dy2=dy2, y=list(y1=y1, y2=y2)));
}
# TODO:
# - handle multiple values for functions;

### Test
b0 = 2
c1 = 3
x = seq(-5, 5, by=0.1)
y = y.f(x, b0=b0, c1=c1)
### Plot
xr = rep(x, each=3);
ylim = range.c(y); isRe = ylim$isRe;
plot(xr[isRe[[1]]], y[[1]][isRe[[1]]], type="l", ylim=ylim$rg);
lines(xr[isRe[[2]]], y[[2]][isRe[[2]]], col="darkgreen");
### dy
x = c(-4, -2, 0, 1.5)
xr = rep(x, each=3)
dy = dy.f(x, b0=b0, c1=c1);
y = dy$y; dy = dy[-3];
isRe = isRe.f(dy);
dy = Re.f(dy, isRe);
y  = Re.f(y, isRe);
xr = Re.f(xr, isRe);
line.tan(xr[[1]], dx=3, p=y[[1]], dp=dy[[1]], col="green")
line.tan(xr[[2]], dx=3, p=y[[2]], dp=dy[[2]], col="red")


#########################
#########################

### Derived from:
# (y1*y2)^n + c1*y2 = R1
# (y1*y2)^n + c1*y1 = R2
# where (y1, y2) = functions in x;
# R1, R2 = are given parameters / functions;
# c1 = numeric constant;

### D & Mult(y1*y2*...) =>
n*(R1 - c1*y2)*(y2*dy1 + y1*dy2) + c1*y1*y2*dy2 - y1*y2*dR1 # = 0
n*(R2 - c1*y1)*(y2*dy1 + y1*dy2) + c1*y1*y2*dy1 - y1*y2*dR2 # = 0

### Examples:
# R1 = R2 = x + b0;
n*(x + b0 - c1*y2)*(y2*dy1 + y1*dy2) + c1*y1*y2*dy2 - y1*y2 # = 0
n*(x + b0 - c1*y1)*(y2*dy1 + y1*dy2) + c1*y1*y2*dy1 - y1*y2 # = 0
# =>
n*c1*y2^2*dy1 + (n-1)*c1*y1*y2*dy2 - n*(x + b0)*(y2*dy1 + y1*dy2) + y1*y2 # = 0
n*c1*y1^2*dy2 + (n-1)*c1*y1*y2*dy1 - n*(x + b0)*(y2*dy1 + y1*dy2) + y1*y2 # = 0


y2*(n*c1*y2 - n*(x + b0))*dy1 + y1*((n-1)*c1*y2 - n*(x + b0))*dy2 + y1*y2 # = 0
y2*((n-1)*c1*y1 - n*(x + b0))*dy1 + y1*(n*c1*y1 - n*(x + b0))*dy2 + y1*y2 # = 0


### Solution:
y.f = function(x, b0=c(1,1), c1=c(1,1), n=3, asMatrix=FALSE) {
	if(length(b0) == 1) b0 = rep(b0, 2);
	if(length(c1) == 1) c1 = rep(c1, 2);
	b0s = sum(b0);
	solve.ypr = function(x) {
		coeff = c(1, rep(0, n-1), - (2*x + b0s), rep(0, n-1), prod(x + b0));
		coeff[2*n] = - prod(c1);
		return(roots(coeff));
	}
	ypr = sapply(x, solve.ypr);
	yn = round0(ypr^n);
	x = rep(x, each=2*n);
	y1 = (x + b0[2] - yn)/c1[2]; # Eq 2: + c1[2]*y1;
	y2 = (x + b0[1] - yn)/c1[1];
	if( ! asMatrix) {
		sol = list(y1=as.vector(y1), y2=as.vector(y2));
	} else {
		sol = list(y1=y1, y2=y2); # ??? also Matrix ???
	}
	return(sol);
}
dy.f = function(x, b0=c(1,1), c1=c(1,1), n=3) {
	y.all = y.f(x, b0=b0, c1=c1, n=n);
	y1 = y.all[[1]]; y2 = y.all[[2]];
	x = rep(x, each=2*n);
	nxb1 = n*(x + b0[1]); nxb2 = n*(x + b0[2]); dnxb = nxb2 - nxb1;
	ypr = y1*y2;
	div = ypr*(n*c1[1]*y2 - nxb1)*(n*c1[2]*y1 - nxb2) +
		- ypr*((n-1)*c1[1]*y2 - nxb1)*((n-1)*c1[2]*y1 - nxb2);
	dy1 = - y1*ypr*(n*c1[2]*y1 - (n-1)*c1[1]*y2 - dnxb) / div;
	dy2 = - y2*ypr*(n*c1[1]*y2 - (n-1)*c1[2]*y1 + dnxb) / div;
	return(list(dy1=dy1, dy2=dy2, y=list(y1=y1, y2=y2)));
}
# TODO:
# - handle multiple values for functions;


### Test
b0 = c(1,1)
c1 = c(-1,-1)
n = 3
x = seq(-4, 4, by=0.1)
y = y.f(x, b0=b0, c1=c1, n=n)
### Plot
xr = rep(x, each=2*n);
ylim = range.c(y); isRe = ylim$isRe;
plot(xr[isRe[[1]]], y[[1]][isRe[[1]]], type="l", ylim=ylim$rg);
lines(xr[isRe[[2]]], y[[2]][isRe[[2]]], col="darkgreen");
### dy
x = c(-0.75, 0, 1, 2)
xr = rep(x, each=2*n)
dy = dy.f(x, b0=b0, c1=c1, n=n);
y = dy$y; dy = dy[-3];
isRe = isRe.f(dy);
dy = Re.f(dy, isRe);
y  = Re.f(y, isRe);
xr = Re.f(xr, isRe);
line.tan(xr[[1]], dx=3, p=y[[1]], dp=dy[[1]], col="green")
line.tan(xr[[2]], dx=3, p=y[[2]], dp=dy[[2]], col="red")


#########
### Ex 2:
b0 = c(1, 2)
c1 = c(-1, -2)
n = 3
x = seq(-4, 4, by=0.1)
y = y.f(x, b0=b0, c1=c1, n=n)
### Plot
xr = rep(x, each=2*n);
ylim = range.c(y); isRe = ylim$isRe;
plot(xr[isRe[[1]]], y[[1]][isRe[[1]]], type="l", ylim=ylim$rg);
lines(xr[isRe[[2]]], y[[2]][isRe[[2]]], col="darkgreen");
### TODO:
# Note: also wrong separation (a few points);
plot(xr[isRe[[1]]][c(T,F)], y[[1]][isRe[[1]]][c(T,F)], type="l", ylim=ylim$rg);
lines(xr[isRe[[1]]][c(F,T)], y[[1]][isRe[[1]]][c(F,T)], type="l", ylim=ylim$rg);
lines(xr[isRe[[2]]][c(T,F)], y[[2]][isRe[[2]]][c(T,F)], col="darkgreen");
lines(xr[isRe[[2]]][c(F,T)], y[[2]][isRe[[2]]][c(F,T)], col="darkgreen");

### dy
x = c(-1.25, 0, 1, 2.5)
xr = rep(x, each=2*n)
dy = dy.f(x, b0=b0, c1=c1, n=n);
y = dy$y; dy = dy[-3];
isRe = isRe.f(dy);
dy = Re.f(dy, isRe);
y  = Re.f(y, isRe);
xr = Re.f(xr, isRe);
line.tan(xr[[1]], dx=3, p=y[[1]], dp=dy[[1]], col="green")
line.tan(xr[[2]], dx=3, p=y[[2]], dp=dy[[2]], col="red")


