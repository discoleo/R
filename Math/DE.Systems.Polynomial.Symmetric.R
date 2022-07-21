########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### DE Systems: Polynomial
###
### draft v.0.1c


#############
### Types ###
#############

### Symmetric:
# - Simple, derived from:
#   y1^n + y2^n = R1;


####################

### Helper Functions

# include: DE.ODE.Helper.R;
source("DE.ODE.Helper.R")

# TODO:
# see also helper functions in DE.Systems.Polynomial.R;

unity = function(n=3, all=TRUE) {
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	if(all) {
		m = m^(0:(n-1))
	}
	return(m)
}


#################

#################
### Symmetric ###
#################

### Derived from:
# y1^n + y2^n = R1
# y1*y2 = cx
# where (y1, y2) = functions in x;
# R, c are given parameters / functions;


### D(Eq 1) =>
y1^(n-1)*dy1 + y2^(n-1)*dy2 - dR1 / n # = 0
# * y1 or * y2 =>
y1^n*dy1 + cx*y2^(n-2)*dy2 - dR1*y1 / n # = 0
y2^n*dy2 + cx*y1^(n-2)*dy1 - dR1*y2 / n # = 0

### Variant:
# - substituting y[i]^n:
(R1 - y2^n)*dy1 + cx*y2^(n-2)*dy2 - dR1*y1 / n # = 0
(R1 - y1^n)*dy2 + cx*y1^(n-2)*dy1 - dR1*y2 / n # = 0
# =>
y2^n*dy1 - cx*y2^(n-2)*dy2 - R1*dy1 + dR1*y1 / n # = 0
y1^n*dy2 - cx*y1^(n-2)*dy1 - R1*dy2 + dR1*y2 / n # = 0

### Note:
# - resulting system is NOT symmetric, but Hetero-Symmetric;


###############
### Order 3 ###
###############

### y1^3 + y2^3 = R1

### D =>
y1^2*dy1 + y2^2*dy2 - dR1/3 # = 0
### *y1 OR *y2 =>
y1^3*dy1 + cx*y2*dy2 - y1*dR1/3 # = 0
cx*y1*dy1 + y2^3*dy2 - y2*dR1/3 # = 0
# =>
# System variant:
y2^3*dy1 - R1*dy1 - cx*y2*dy2 + y1*dR1/3 # = 0
y1^3*dy2 - R1*dy2 - cx*y1*dy1 + y2*dR1/3 # = 0
# Another variant:
y2^2*(dcx - y1*dy2) - R1*dy1 - cx*y2*dy2 + y1*dR1/3 # = 0
y1^2*(dcx - y2*dy1) - R1*dy2 - cx*y1*dy1 + y2*dR1/3 # = 0
# =>
2*cx*y2*dy2 + R1*dy1 - y2^2*dcx - y1*dR1/3 # = 0
2*cx*y1*dy1 + R1*dy2 - y1^2*dcx - y2*dR1/3 # = 0

### D2 =>
3*y2^2*dy1*dy2 + y2^3*d2y1 - R1*d2y1 - dR1*dy1 +
	- cx*y2*d2y2 - cx*(dy2)^2 - dcx*y2*dy2 + dy1*dR1/3 + y1*d2R1/3 # = 0
# TODO


### Example:
# R1 = b01; cx = x + b02;
2*(x + b02)*y2*dy2 + b01*dy1 - y2^2 # = 0
2*(x + b02)*y1*dy1 + b01*dy2 - y1^2 # = 0


### Solution:
# y.f():
y.f = function(x, b0=c(1,1), n=3, asMatrix=FALSE) {
	if(length(b0) == 1) b0 = rep(b0, 2);
	m = unity(n, all=TRUE);
	d = b0[1]/2;
	solve.exact = function(x) {
		cx = x + b0[2];
		det = d^2 - cx^n;
		det = if(det >= 0) sqrt(det) else complex(re=0, im=sqrt(-det));
		y1 = rootn(d + det, n);
		y1 = y1 * m;
		y2 = rootn(d - det, n);
		y2 = y2 / m;
		if(n %% 2 == 0 && cx < 0) y2 = -y2;
		return(cbind(y1, y2));
	}
	sol = round0(sapply(x, solve.exact));
	y1 = sol[1:n,]; y2 = sol[seq(n+1, 2*n),];
	if( ! asMatrix) {
		sol = list(y1=as.vector(y1), y2=as.vector(y2));
	} else {
		sol = list(y1=y1, y2=y2); # ??? also Matrix ???
	}
	return(sol);
}
# dy.f():
dy.f = function(x, b0=c(1,1), n=3) {
	y.all = y.f(x, b0=b0, n=n);
	y1 = y.all[[1]]; y2 = y.all[[2]];
	x = rep(x, each=n);
	xb = x + b0[2]; b01 = b0[1];
	ypr = (y1*y2); yprn2 = ypr^(n-2);
	div = b01^2 - 4*xb^2*yprn2;
	#
	dy1 = (b01*y2^(n-1) - 2*xb*yprn2*y1) / div;
	dy2 = (b01*y1^(n-1) - 2*xb*yprn2*y2) / div;
	dy1[div == 0] = 0; dy2[div == 0] = 0; # TODO
	#
	return(list(dy1=dy1, dy2=dy2, y=list(y1=y1, y2=y2)));
}

### Test
b0 = c(1, 2)
x = seq(-4, -1, by=0.025)
y = y.f(x, b0=b0)
### Plot
xr = rep(x, each=3);
ylim = range.c(y); isRe = ylim$isRe;
plot(xr[isRe[[1]]], y[[1]][isRe[[1]]], type="l", ylim=ylim$rg);
lines(xr[isRe[[2]]], y[[2]][isRe[[2]]], col="darkgreen");

### dy
x = c(-3.5, -3, -1.75, -1.5)
xr = rep(x, each=3)
dy = dy.f(x, b0=b0);
y = dy$y; dy = dy[-3];
isRe = isRe.f(dy);
dy = Re.f(dy, isRe);
y  = Re.f(y, isRe);
xr = Re.f(xr, isRe);
line.tan(xr[[1]], dx=1.5, p=y[[1]], dp=dy[[1]], col="green")
line.tan(xr[[2]], dx=1.5, p=y[[2]], dp=dy[[2]], col="red")

### old code:
# plot(xr[isRe[[1]]][c(T,F)], y[[1]][isRe[[1]]][c(T,F)], type="l", ylim=ylim$rg);
# lines(xr[isRe[[1]]][c(F,T)], y[[1]][isRe[[1]]][c(F,T)], type="l", ylim=ylim$rg);
# lines(xr[isRe[[2]]][c(T,F)], y[[2]][isRe[[2]]][c(T,F)], col="darkgreen");
# lines(xr[isRe[[2]]][c(F,T)], y[[2]][isRe[[2]]][c(F,T)], col="darkgreen");


###################
### Other Variants:

# [old]

### Eq =>
(y1 + y2)^3 - 3*cx*(y1 + y2) - R1 # = 0
### D =>
3*(y1 + y2)^2*(dy1 + dy2) - 3*dcx*(y1 + y2) - 3*cx*(dy1 + dy2) - dR1 # = 0
### Mult((y1+y2)*) =>
3*(R1 + 3*cx*(y1 + y2))*(dy1 + dy2) - 3*dcx*(y1 + y2)^2 +
	- 3*cx*(y1 + y2)*(dy1 + dy2) - (y1 + y2)*dR1 # = 0
2*cx*(y1 + y2)*(dy1 + dy2) + R1*(dy1 + dy2) - dcx*(y1 + y2)^2 - (y1 + y2)*dR1/3 # = 0

### Eq 2 =>
y1*dy2 + y2*dy1 - dcx # = 0

### (Eq-2bis)*y1 =>
2*cx*(y1^2 - y2^2)*dy1 + R1*(y1 - y2)*dy1 +
	- dcx*y1*(y1 + y2)^2 - y1*(y1 + y2)*dR1/3 + 2*cx*dcx*(y1 + y2) + R1*dcx # = 0
### (Eq-2bis)*y2 =>
2*cx*(y2^2 - y1^2)*dy2 + R1*(y2 - y1)*dy2 +
	- dcx*y2*(y1 + y2)^2 - y2*(y1 + y2)*dR1/3 + 2*cx*dcx*(y1 + y2) + R1*dcx # = 0

### System:
2*cx*(y1^2 - y2^2)*dy1 + R1*(y1 - y2)*dy1 +
	- dcx*y1*(y1 + y2)^2 - y1*(y1 + y2)*dR1/3 + 2*cx*dcx*(y1 + y2) + R1*dcx # = 0
2*cx*(y2^2 - y1^2)*dy2 + R1*(y2 - y1)*dy2 +
	- dcx*y2*(y1 + y2)^2 - y2*(y1 + y2)*dR1/3 + 2*cx*dcx*(y1 + y2) + R1*dcx # = 0

### Variant
# D(Eq 1) - 2*cx*(y1*dy2 + y2*dy1) + 2*dcx*y1*y2 =>
2*cx*(y1*dy1 + y2*dy2) + R1*(dy1 + dy2) - dcx*(y1^2 + y2^2) - 2*dcx*y1*y2 +
	- (y1 + y2)*dR1/3 + 2*cx*dcx # = 0
# Eq 2 variant:
3*y1^3*dy1 + 3*y2^3*dy2 + 3*cx*(y1*dy1 + y2*dy2) - dR1*(y1 + y2) # = 0


### Example:
# R1 = b01; cx = x + b02;
2*(x+b02)*(y1^2 - y2^2)*dy1 + b01*(y1 - y2)*dy1 +
	- y1*(y1 + y2)^2 + 2*(x+b02)*(y1 + y2) + b01 # = 0
2*(x+b02)*(y2^2 - y1^2)*dy2 + b01*(y2 - y1)*dy2 +
	- y2*(y1 + y2)^2 + 2*(x+b02)*(y1 + y2) + b01 # = 0


### Solution:

# y.f(): see in previous section;
y.f.old = function(x, b0=c(1,1), asMatrix=FALSE) {
	if(length(b0) == 1) b0 = rep(b0, 2);
	solve.ypr = function(x) {
		coeff = c(1, 0, 0, - b0[1], 0, 0, (x + b0[2])^3);
		return(roots(coeff));
	}
	y1 = round0(sapply(x, solve.ypr));
	x  = rep(x, each=6);
	y2 = (x + b0[2]) / y1;
	y2[y1 == 0] = 0; # TODO: check;
	if( ! asMatrix) {
		sol = list(y1=as.vector(y1), y2=as.vector(y2));
	} else {
		sol = list(y1=y1, y2=y2); # ??? also Matrix ???
	}
	return(sol);
}
dy.f = function(x, b0=c(1,1), n=3) {
	# currently only n=3 !!
	y.all = y.f(x, b0=b0, n=n);
	y1 = y.all[[1]]; y2 = y.all[[2]];
	x = rep(x, each=n);
	xb = x + b0[2]; b01 = b0[1];
	ys = y1 + y2; nyb = 2*xb*ys + b01;
	div = 2*xb*(y1^2 - y2^2) + b01*(y1 - y2);
	#
	dy1 =   (y1*ys^2 - nyb) / div;
	dy2 = - (y2*ys^2 - nyb) / div;
	dy1[div == 0] = 0; dy2[div == 0] = 0; # TODO
	#
	return(list(dy1=dy1, dy2=dy2, y=list(y1=y1, y2=y2)));
}
# TODO:
# - handle multiple values for functions;


### Test
b0 = c(1, 2)
x = seq(-4, -1, by=0.025)
y = y.f(x, b0=b0)
### Plot
xr = rep(x, each=3);
ylim = range.c(y); isRe = ylim$isRe;
plot(xr[isRe[[1]]], y[[1]][isRe[[1]]], type="l", ylim=ylim$rg);
lines(xr[isRe[[2]]], y[[2]][isRe[[2]]], col="darkgreen");

### dy
x = c(-3.5, -3, -1.75, -1.5)
xr = rep(x, each=3)
dy = dy.f(x, b0=b0);
y = dy$y; dy = dy[-3];
isRe = isRe.f(dy);
dy = Re.f(dy, isRe);
y  = Re.f(y, isRe);
xr = Re.f(xr, isRe);
line.tan(xr[[1]], dx=1.5, p=y[[1]], dp=dy[[1]], col="green")
line.tan(xr[[2]], dx=1.5, p=y[[2]], dp=dy[[2]], col="red")

### with old code:
# plot(xr[isRe[[1]]][c(T,F)], y[[1]][isRe[[1]]][c(T,F)], type="l", ylim=ylim$rg);
# lines(xr[isRe[[1]]][c(F,T)], y[[1]][isRe[[1]]][c(F,T)], type="l", ylim=ylim$rg);
# lines(xr[isRe[[2]]][c(T,F)], y[[2]][isRe[[2]]][c(T,F)], col="darkgreen");
# lines(xr[isRe[[2]]][c(F,T)], y[[2]][isRe[[2]]][c(F,T)], col="darkgreen");


#########
### Ex 2:
b0 = c(-2, -3)
x = seq(-2, 4, by=0.025)
y = y.f(x, b0=b0)
### Plot
xr = rep(x, each=3);
ylim = range.c(y); isRe = ylim$isRe;
plot(xr[isRe[[1]]], y[[1]][isRe[[1]]], type="l", ylim=ylim$rg);
lines(xr[isRe[[2]]], y[[2]][isRe[[2]]], col="darkgreen");

### dy
x = c(-1.5, -0.5, 1, 2, 2.75, 3.75)
xr = rep(x, each=3)
dy = dy.f(x, b0=b0);
y = dy$y; dy = dy[-3];
isRe = isRe.f(dy);
dy = Re.f(dy, isRe);
y  = Re.f(y, isRe);
xr = Re.f(xr, isRe);
line.tan(xr[[1]], dx=1.5, p=y[[1]], dp=dy[[1]], col="green")
line.tan(xr[[2]], dx=1.5, p=y[[2]], dp=dy[[2]], col="red")

### with old code:
# plot(xr[isRe[[1]]][c(T,F)], y[[1]][isRe[[1]]][c(T,F)], type="l", ylim=ylim$rg);
# lines(xr[isRe[[1]]][c(F,T)], y[[1]][isRe[[1]]][c(F,T)], type="l", ylim=ylim$rg);
# lines(xr[isRe[[2]]][c(T,F)], y[[2]][isRe[[2]]][c(T,F)], col="darkgreen");
# lines(xr[isRe[[2]]][c(F,T)], y[[2]][isRe[[2]]][c(F,T)], col="darkgreen");


####################
####################

###############
### Order 4 ###
###############

# y1^4 + y2^4 = R1
# y1 * y2 = cx

### D =>
y1^3*dy1 + y2^3*dy2 - dR1/4 # = 0
### *y1 OR *y2 =>
y1^4*dy1 + y1*y2^3*dy2 - y1*dR1/4 # = 0
y2*y1^3*dy1 + y2^4*dy2 - y2*dR1/4 # = 0
# =>
(R1 - y2^4)*dy1 + y1*y2^3*dy2 - y1*dR1/4 # = 0
y2*y1^3*dy1 + (R1 - y1^4)*dy2 - y2*dR1/4 # = 0
# Substituting D(Eq 2) =>
y2^3*(y1*dy2 - dcx) + y1*y2^3*dy2 + R1*dy1 - y1*dR1/4 # = 0
y1^3*(y2*dy1 - dcx) + y2*y1^3*dy1 + R1*dy2 - y2*dR1/4 # = 0

### System:
2*cx*y2^2*dy2 + R1*dy1 - dcx*y2^3 - y1*dR1/4 # = 0
2*cx*y1^2*dy1 + R1*dy2 - dcx*y1^3 - y2*dR1/4 # = 0


### Example:
# R1 = b01; cx = x + b02;
2*(x + b02)*y2^2*dy2 + b01*dy1 - y2^3 # = 0
2*(x + b02)*y1^2*dy1 + b01*dy2 - y1^3 # = 0

### Solution:
# - uses the generalized functions from the previous section;


### Test
n = 4
b0 = c(1, 2)
x = seq(-4, -1, by=0.025)
y = y.f(x, b0=b0, n=n)
### Plot
xr = rep(x, each=n);
ylim = range.c(y); isRe = ylim$isRe;
# plot(xr[isRe[[1]]], y[[1]][isRe[[1]]], type="l", ylim=ylim$rg);
# lines(xr[isRe[[2]]], y[[2]][isRe[[2]]], col="darkgreen");
# TODO: still needed for n > 3
id = diag(T, n);
isT = isRe[[1]] & id[1,]; plot(xr[isT], y[[1]][isT], type="l", ylim=ylim$rg);
isT = isRe[[2]] & id[1,]; lines(xr[isT], y[[2]][isT], type="l", ylim=ylim$rg, col="darkgreen");
zero = sapply(seq(2, n), function(id1) {
	isT = isRe[[1]] & id[id1,]; lines(xr[isT], y[[1]][isT], type="l", ylim=ylim$rg);
	isT = isRe[[2]] & id[id1,]; lines(xr[isT], y[[2]][isT], type="l", ylim=ylim$rg, col="darkgreen");
})

### dy
x = c(-2.65, -2.25, -1, -1.35)
xr = rep(x, each=n)
dy = dy.f(x, b0=b0, n=n);
y = dy$y; dy = dy[-3];
isRe = isRe.f(dy);
dy = Re.f(dy, isRe);
y  = Re.f(y, isRe);
xr = Re.f(xr, isRe);
line.tan(xr[[1]], dx=1.5, p=y[[1]], dp=dy[[1]], col="green")
line.tan(xr[[2]], dx=1.5, p=y[[2]], dp=dy[[2]], col="red")


####################
####################

##################
### Extensions ###
### Type A     ###
##################

# y1^n + y2^n = R1
# y1 * y2  + c1*y1 + c2*y2 = cx

### Note:
# - the polynomial system is easier to solve
#   when c1 == c2;

### D =>
y1^(n-1)*dy1 + y2^(n-1)*dy2 - dR1/n # = 0
### *y1 OR *y2 =>
y1^n*dy1 + y1*y2^(n-1)*dy2 - y1*dR1/n # = 0
y2*y1^(n-1)*dy1 + y2^n*dy2 - y2*dR1/n # = 0
# =>
(R1 - y2^n)*dy1 + y1*y2^(n-1)*dy2 - y1*dR1/n # = 0
y2*y1^(n-1)*dy1 + (R1 - y1^n)*dy2 - y2*dR1/n # = 0
# Substituting D(Eq 2) =>
# - assuming (c1, c2) = numeric constants;
y2^(n-1)*(y1*dy2 + c1*dy1 + c2*dy2 - dcx) + y1*y2^(n-1)*dy2 + R1*dy1 - y1*dR1/n # = 0
y1^(n-1)*(y2*dy1 + c1*dy1 + c2*dy2 - dcx) + y2*y1^(n-1)*dy1 + R1*dy2 - y2*dR1/n # = 0

### System:
2*(cx - c1*y1 - c2*y2)*y2^(n-2)*dy2 + c1*y2^(n-1)*dy1 + c2*y2^(n-1)*dy2 + R1*dy1 - dcx*y2^(n-1) - y1*dR1/n # = 0
2*(cx - c1*y1 - c2*y2)*y1^(n-2)*dy1 + c1*y1^(n-1)*dy1 + c2*y1^(n-1)*dy2 + R1*dy2 - dcx*y1^(n-1) - y2*dR1/n # = 0

### Special Cases:

### n = 2:
  c1*y2*dy1 + R1*dy1 - (2*c1*y1 + c2*y2)*dy2 + 2*cx*dy2 - dcx*y2 - y1*dR1/2 # = 0
- (c1*y1 + 2*c2*y2)*dy1 + 2*cx*dy1 + c2*y1*dy2 + R1*dy2 - dcx*y1 - y2*dR1/2 # = 0
### n = 2; c1 = c2;
  c1*y2*dy1 + R1*dy1 - c1*(2*y1 + y2)*dy2 + 2*cx*dy2 - dcx*y2 - y1*dR1/2 # = 0
- c1*(y1 + 2*y2)*dy1 + 2*cx*dy1 + c2*y1*dy2 + R1*dy2 - dcx*y1 - y2*dR1/2 # = 0

