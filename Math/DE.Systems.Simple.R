########################
##
## Leonard Mada
## [the one and only]
##
## Systems of Differential Equations
## Conversion to 1 ODE
##
## draft v.0.1d



### Transformations
# A. Conversion of a linear system of ODEs to *one* ODE
# B. Conversion of a non-linear system of ODEs to *one* ODE


####################

### Helper Functions


library(deSolve)

# source("Polynomials.Helper.R")

subst.part = function(x, param) {
	if(! inherits(param, "list")) param = as.list(param);
	do.call(substitute, list(x, param));
}

# Note:
# - Functions to simplify the expressions are in file:
#   Differentiation.R;


##########################

### Simple Coupled Systems

# dy1 = b1*y1 + b2*y2 + f1;
# dy2 = c1*y1 + c2*y2 + f2;

# =>
c1*dy1 - b1*dy2 - c1*(b2*y2 + f1) + b1*(c2*y2 + f2) # = 0
c1*dy1 - b1*dy2 + (b1*c2 - b2*c1)*y2 - c1*f1 + b1*f2 # = 0


### D(Eq Sys 2) =>

### Case: b, c = constants;
d2y2 - c1*dy1 - c2*dy2 - df2 # = 0
d2y2 - (b1 + c2)*dy2 + (b1*c2 - b2*c1)*y2 - c1*f1 + b1*f2 - df2 # = 0

# => Solve ODE of Order 2;
# see Examples below;


### Case: b, c = f(x)
d2y2 - c1*dy1 - dc1*y1 - c2*dy2 - dc2*y2 - df2 # = 0
d2y2 - (b1 + c2)*dy2 - dc1*y1 +
	+ (b1*c2 - b2*c1 - dc2)*y2 - c1*f1 + b1*f2 - df2 # = 0
c1*d2y2 - c1*(b1 + c2)*dy2 - dc1*dy2 +
	+ c1*(b1*c2 - b2*c1 - dc2)*y2 + dc1*c2*y2 +
	- c1^2*f1 + b1*c1*f2 - c1*df2 + dc1*f2 # = 0

# TODO: check;

### Test:

b = c(1,-2,-3,2)
names(b) = c("b1","b2","c1","c2");

### ODE:
d2y2 - 3*dy2 - 4*y2 + 3*f1 + f2 - df2 # = 0


### Solve:
ode.sys = function(t, y, parms, ...) {
	y1 = y[1]; y2 = y[2]; x = t;
	with(as.list(parms), {
		dy1 = b1*y1 + b2*y2 + f1(x);
		dy2 = c1*y1 + c2*y2 + f2(x);
		return(list(c(dy1, dy2)));
	});
}
ode.f = function(t, y, parms, ...) {
	y2 = y[1]; dy2 = y[2]; x = t;
	with(as.list(parms), {
		d2y2 = - (b1 + c2)*dy2 + (b1*c2 - b2*c1)*y2 +
			- c1*f1(x) + b1*f2(x) - df2(x);
		d2y2 = - d2y2;
		return(list(c(dy2, d2y2)));
	});
}

### Test:
params.f = list(
	f1  = \(x)  5/3*x - 1/2,
	f2  = \(x) -1/5*x + 251/257,
	df2 = \(x) -1/5 );

b = c(-3,1,-5,2)
# b = c(1,-4,-3,1); # b = c(1,-2,-3,2);
names(b) = c("b1","b2","c1","c2");
params = c(params.f, as.list(b));

x  = seq(0, 4, by = 1/128);
y1 = 1/2; y2 = 1/5;
# y1 = -2/7; y2 = -1/5;
ys = c(y1, y2);
ye = c(y2 = y2, dy2 = b[3]*y1 + b[4]*y2 + params$f2(x[1]));

sols = ode(ys, x, func = ode.sys, parms = params)
sole = ode(ye, x, func = ode.f,   parms = params)
col = c("#2472E6A0", "#FF643290")
matplot(sols[,1], sols[, 3], type = "l", lwd=2, col=col[1])
matplot(sole[,1], sole[, 2], type = "l", lwd=3, col=col[2], lty=2, add=T)



#####################
#####################

#####################
### System: 3 Eqs ###

# Converting Sys 3 * Linear ODE
# to 1 Linear ODE

# dy1 = b1*y1 + b2*y2 + b3*y3 + f1;
# dy2 = c1*y1 + c2*y2 + c3*y3 + f2;
# dy3 = e1*y1 + e2*y2 + e3*y3 + f3;

### Case: Coefficients = constants;

# D(Eq 3):
d2y3 - e1*dy1 - e2*dy2 - e3*dy3 - df3 # = 0
# Substitute: Eq 1 & Eq 2 =>
d2y3 - (e1*b1 + e2*c1)*y1 - (e1*b2 + e2*c2)*y2 +
	- e1*b3*y3 - e2*c3*y3 - e3*dy3 +
	- e1*f1 - e2*f2 - df3 # = 0

# D2(Eq 3):
d3y3 - e1*d2y1 - e2*d2y2 - e3*d2y3 - d2f3 # = 0
# Substitute: D(Eq 1) & D(Eq 2) =>
d3y3 - e1*(b1*dy1 + b2*dy2 + b3*dy3) +
	- e2*(c1*dy1 + c2*dy2 + c3*dy3) +
	- e3*d2y3 - d2f3 - e1*df1 - e2*df2 # = 0
# Substitute: dy1 & dy2 =>
d3y3 - (e1*b1 + e2*c1)*(b1*y1 + b2*y2 + b3*y3 + f1) +
	- (e1*b2 + e2*c2)*(c1*y1 + c2*y2 + c3*y3 + f2) +
	- e3*d2y3 - (e1*b3 + e2*c3)*dy3 +
	- d2f3 - e1*df1 - e2*df2 # = 0
# => [see Eq 3t in Transformed System]

### Transformed System:

dy3 - e1*y1 - e2*y2 - e3*y3 - f3 # = 0;
d2y3 - (e1*b1 + e2*c1)*y1 - (e1*b2 + e2*c2)*y2 +
	- e1*b3*y3 - e2*c3*y3 - e3*dy3 +
	- e1*f1 - e2*f2 - df3 # = 0;
d3y3 - (b1*(e1*b1 + e2*c1) + c1*(e1*b2 + e2*c2))*y1 +
	- (b2*(e1*b1 + e2*c1) + c2*((e1*b2 + e2*c2)))*y2 +
	- e3*d2y3 - (e1*b3 + e2*c3)*dy3 +
		- b3*(e1*b1 + e2*c1)*y3 - c3*(e1*b2 + e2*c2)*y3 +
	- d2f3 - e1*df1 - e2*df2 +
		- (e1*b1 + e2*c1)*f1 - (e1*b2 + e2*c2)*f2 # = 0;

# Eliminate y1 & y2;
# =>
# ODE of Order 3

### Note:
# - the same steps can be performed when
#   the coefficients are NOT constant;
# - however, the expressions are much uglier;

##############
# [NOT needed]
# Eq1 - Eq3:
e1*dy1 - b1*dy3 + b1*(e2*y2 + e3*y3) +
	- e1*(b2*y2 + b3*y3) + b1*f3 - e1*f1 # = 0
# Eq2 - Eq3:
e1*dy2 - c1*dy3 + c1*(e2*y2 + e3*y3) +
	- e1*(c2*y2 + c3*y3) + c1*f3 - e1*f2 # = 0

# [NOT needed]
# D2(Eq 3b) =>
d3y3 - (b1 + e3)*d2y3 - e2*d2y2 +
	+ b1*(e2*dy2 + e3*dy3) - e1*(b2*dy2 + b3*dy3) +
	+ b1*df3 - e1*df1 - d2f3 # = 0

###############
### Example ###

bb = c(1,2,3)
cc = c(1,3,5)
ee = c(2,5,7)
names(bb) = paste0("b", 1:3)
names(cc) = paste0("c", 1:3)
names(ee) = paste0("e", 1:3)
params = c(bb, cc, ee)

eq3 = expression(e1*y1 + e2*y2 + e3*y3 + f3 - dy3)[[1]]
subst.part(eq3, params)


### Transformed System:
 2 * y1 +  5 * y2 +   7 * y3 + f3 - dy3 # = 0
 7 * y1 + 19 * y2 +  31 * y3 + 2 * f1 + 5 * f2 - df3 - d2y3 + 7 * dy3 # = 0
26 * y1 + 71 * y2 + 116 * y3 +
	+ d2f3 + 2 * df1 + 5 * df2 + 7 * f1 + 19 * f2 +
	- d3y3 + 7 * d2y3 + 31 * dy3 # = 0

### ODE:
d3y3 - 11*d2y3 - 2*dy3 + y3 - d2f3 - 4*df3 - f3 - 2*df1 + f1 - 5*df2 + f2 # = 0

pracma::roots(c(1,-11,-2,1))
# - for exact solution, see file:
#   Polynomials.CardanGeneralisation.R;


# TODO:
# Example with more sensible coefficients;


#########
### Ex 2:

bb = c( 1, 2,-3)
cc = c(-1, 3,-5)
ee = c( 2,-5, 2)
names(bb) = paste0("b", 1:3)
names(cc) = paste0("c", 1:3)
names(ee) = paste0("e", 1:3)
params = c(bb, cc, ee)

eq3 = expression(e1*y1 + e2*y2 + e3*y3 + f3 - dy3)[[1]]
subst.part(eq3, params)


### Transformed System:
 2 * y1 -  5 * y2 +  2 * y3 + f3 - dy3 # = 0
 7 * y1 - 11 * y2 + 19 * y3 +
	+ 2 * f1 - 5 * f2 + df3 - d2y3 + 2 * dy3 # = 0
18 * y1 - 19 * y2 + 34 * y3 +
	+ d2f3 + 2 * df1 + 7 * f1 - 5 * df2 - 11 * f2 +
	- d3y3 + 2 * d2y3 + 19 * dy3 # = 0


### ODE:
d3y3 - 6*d2y3 - 6*dy3 + 32*y3 +
	- d2f3 + 4*df3 - 5*f3 - 2*df1 + f1 + 5*df2 - 9*f2 # = 0

rr = pracma::roots(c(1,-6,-6,32))
eigen(matrix(params, ncol=3, byrow=TRUE))
print(rr)


#####################
#####################

##################
### Non-Linear ###
##################

### Example 1:
# dy1 = b2*y1^2 + b1*y1 + c1*y2 + f1;
# dy2 = d1*y1 + d2*y2 + d3*y2^2 + f2;

# Case: Coefficients = constant
# (for simplicity)

### D(Eq 2)
d2y2 - d1*dy1 - d2*dy2 - 2*d3*y2*dy2 - df2 # = 0
# Substitution =>
d2y2 - d1*(b2*y1^2 + b1*y1) +
	- 2*d3*y2*dy2 - d2*dy2 - c1*d1*y2 - df2 - d1*f1 # = 0

### Gaussian Elimination:
# Eq 2 & Eq 2b:

### ODE:
d1*d2y2 - b2*dy2^2 + 2*d3*b2*y2^2*dy2 + 2*(d2*b2 - d3*d1)*y2*dy2 +
	- (d2*d1 + d1*b1 + 2*f2*b2)*dy2 +
	- d3^2*b2*y2^4 - 2*d3*d2*b2*y2^3 +
	+ (d3*d1*b1 - d2^2*b2 - 2*f2*d3*b2)*y2^2 +
	+ (d2*d1*b1 - d1^2*c1 - 2*f2*d2*b2)*y2 +
	- d1^2*f1 - d1*df2 + d1*b1*f2 - b2*f2^2 # = 0


# p1 = as.pm("d1*y1 + d2*y2 + d3*y2^2 + f2 - dy2")
# p2 = as.pm(...)
# pR = solve.pm(p1, p2, by="y1")


#################

### Example 2:
# dy1 = b1*y1 + b2*y2^2;
# dy2 = c1*y2 + c2*y1^2;

# Case: Coefficients = constant
# (for simplicity)

# D(Eq 2) =>
d2y2 - c1*dy2 - 2*c2*y1*dy1 # = 0
# Subst =>
d2y2 - c1*dy2 - 2*c2*y1*(b1*y1 + b2*y2^2) # = 0

### Transformed System:
dy2  - c1*y2  - c2*y1^2 # = 0
d2y2 - c1*dy2 - 2*b1*c2*y1^2 - 2*b2*c2*y1*y2^2 # = 0


### ODE:
d2y2^2 - (4*b1 + 2*c1)*dy2*d2y2 + 4*b1*c1*y2*d2y2 +
	+ (4*b1^2 + 4*b1*c1 + c1^2)*dy2^2 +
	- 4*c2*b2^2*y2^4*dy2 - 4*b1*c1*(2*b1 + c1)*y2*dy2 +
	+ 4*b2^2*c1*c2*y2^5 + 4*b1^2*c1^2*y2^2 # = 0

# TODO: check;


# p1 = as.pm("c1*y2  + c2*y1^2 - dy2")
# p2 = as.pm(...)
# pR = solve.pm(p1, p2, by="y1")
# print.pm(sort.dpm(pR$Rez, y = "y2", x = NULL), lead = NA)

### Alternative Eq 2:
d2y2 - c1*(c1*y2  + c2*y1^2) - 2*b1*c2*y1^2 - 2*b2*c2*y1*y2^2 # = 0
# Note: generates same ODE;


### Solve:
ode.sys = function(t, y, parms, ...) {
	with(as.list(parms), {
		y1  = y[1]; y2 = y[2]; x = t;
		dy1 = b1*y1 + b2*y2^2;
		dy2 = c1*y2 + c2*y1^2;
		return(list(c(dy1, dy2)));
	});
}
ode.f = function(t, y, parms, ...) {
	with(as.list(parms), {
		x  = t;
		y2 = y[1]; dy2 = y[2]; d2y2 = y[3];
		# Based on D3:
		div  = 2*d2y2 - (4*b1 + 2*c1)*dy2 + 4*b1*c1*y2;
		d3y2 = - (4*b1 + 2*c1)*d2y2^2 +
			+ 4*b1*c1*dy2*d2y2 + 2*(4*b1^2 + 4*b1*c1 + c1^2)*dy2*d2y2 +
			- 4*c2*b2^2*y2^4*d2y2 - 4*b1*c1*(2*b1 + c1)*y2*d2y2 +
			- 16*c2*b2^2*y2^3*dy2^2 - 4*b1*c1*(2*b1 + c1)*dy2^2 +
			+ 20*b2^2*c1*c2*y2^4*dy2 + 8*b1^2*c1^2*y2*dy2;
		d3y2 = - d3y2 / div;
		return(list(c(dy2, d2y2, d3y2)));
	});
}
d2y2.f = function(x, y, params) {
	y1 = y[1]; y2 = y[2];
	with(params, {
		dy2  = c1*y2 + c2*y1^2;
		d2y2 = - c1*dy2 - 2*b1*c2*y1^2 - 2*b2*c2*y1*y2^2;
		d2y2 = - d2y2;
		return(c(dy2 = dy2, d2y2 = d2y2));
	});
}

### Example:
b = c(1/3,-2,-1/2,2)
# b = c(1,-4,-3,1); # b = c(1,-2,-3,2);
names(b) = c("b1","b2","c1","c2");
params = as.list(b);

x  = seq(0, 4, by = 1/512);
y1 = 1/2; y2 = 1/5;
# x  = seq(0, 2, by = 1/512); y1 = 4/5; y2 = 1/7;
ys = c(y1, y2);
ye = c(y2 = y2, d2y2.f(x[1], ys, params));

sols = ode(ys, x, func = ode.sys, parms = params)
sole = ode(ye, x, func = ode.f,   parms = params)
col = c("#2472E6A0", "#FF643290")
matplot(sole[,1], sole[, 2], type = "l", lwd=3, col=col[2], lty=2)
matplot(sols[,1], sols[, 3], type = "l", lwd=2, col=col[1], add=T)

