
### Transformations

# Conversion of a linear system of ODEs to *one* ODE


####################

### Helper Functions

# source("Polynomials.Helper.R")

subst.part = function(x, param) {
	if(! inherits(param, "list")) param = as.list(param);
	do.call(substitute, list(x, param));
}

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

# TODO: check;

### Case: b, c = f(x)
d2y2 - c1*dy1 - dc1*y1 - c2*dy2 - dc2*y2 - df2 # = 0
d2y2 - (b1 + c2)*dy2 - dc1*y1 +
	+ (b1*c2 - b2*c1 - dc2)*y2 - c1*f1 + b1*f2 - df2 # = 0
c1*d2y2 - c1*(b1 + c2)*dy2 - dc1*dy2 +
	+ c1*(b1*c2 - b2*c1 - dc2)*y2 + dc1*c2*y2 +
	- c1^2*f1 + b1*c1*f2 - c1*df2 + dc1*f2 # = 0

# TODO: check;


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

