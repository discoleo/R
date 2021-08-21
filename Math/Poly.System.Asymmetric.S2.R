########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Asymmetric S2:
### Base Types
###
### draft v.0.4k-robust


### Asymmetric Polynomial Systems: 2 Variables
### Base Types

### Reducible

### Example 1: (x*y)-Variants
# x^n + b1*x*y = R1
# y^n + b2*x*y = R2
# where it is possible to transform into P(x*y);
# => solve directly for (x*y).

### Example 2: Binomial Expansions
# x^n + P1(x, y) = R1
# y^n + P2(x, y) = R2
# where: x^n + y^n  + P1(x, y) + P2(x, y) = P(S) = R1 + R2;
# => solve directly for S = x+y;


### Example 3: Decomposable System
# x^n * P1(x, y) = R1 * P2(x, y)
# y^n * P2(x, y) = R2 * P1(x, y)

### Generalization 2:
# Q1(x, y) * P1(x, y) = R1 * P2(x, y)
# Q2(x, y) * P2(x, y) = R2 * P1(x, y)
# => Decomposition into 2 distinct Subsystems:

### Sub-Sys 1:
# P1(x, y) = 0
# P2(x, y) = 0

### Sub-Sys 2:
# Q1(x, y) * Q2(x, y) = R1*R2;


### Complicated

### Example 1:
# x^n + b*y = R1
# y^n + b*x = R2
# R1 != R2

### Example 2:
# x^n + b1*y = R
# y^n + b2*x = R
# b1 != b2

### Example 3: d given
# (x + d)^n + y^n = R1
# x^n + (y + d)^n = R2
# R1 != R2

### Example 4: d unknown
# (x + d)^n + y^n = R1
# x^n + (y + d)^n = R2
# x^n + y^n = R3
# R1 != R2

### TODO:
# - understand the advantages of the Dual system;


###############
### History ###
###############


### draft v.0.4k - v.0.4k-robust:
# - systems with Equalities of variables:
#   Order 1: x == y;
#   Order 3: x^3 = y^3; [v.0.4k-O3]
# - robust roots for Order 3; [v.0.4k-robust]
### draft v.0.4j - v.0.4j-sol:
# - decomposition using: (x^2+y^2) & x*y;
#   x^3*y + b1*y^2 = R1;
### draft v.0.4g - v.0.4i:
# - simple xy-Variant:
#   x^3*y^3 + b1*x^2*y + b2*x*y^2 = R1; (trivial P[6])
# - simple Cross-Product example:
#   x*(x^3 + y^3 + b1) = R1*(x*y + b2);
# - TODO: variants with roots {(x,y), m*(x,y), m^2*(x,y)};
### draft v.0.4e - v.0.4f:
# - more Cross-Products:
#   Ex 1: x^2*(x + b) = R1*(y + b); [v.0.4f]
#   Ex 2: x^2*(x^2 + b11*y + b10) = R1*(y^2 + b21*x + b20);
### draft v.0.4d:
# - better classification:
#   Reducible vs Complicated;
### draft v.0.4c:
# - special case: coupled asymmetry
#   x^2*y + b[1]*x^2 - R/b[2]*x = R;
### draft v.0.4b:
# - complex transformation (symmetry breaking) of:
#   x^3 + b*y = R;
### draft v.0.4a:
# - generalized:
#   a11*x^3*y + a12*x*y^3 + b12*(x*y)^2 + b11*(x*y) = R1;
### draft v.0.3n:
# - solved [Cross-Product] Order 3+1:
#   x^3*y + b3*(x*y)^2 + b1*x^2 = R;
### draft v.0.3m - v.0.3m-ext:
# - solved [Cross-Product] type:
#   x^2*y + b3*x*y + b1*x = R;
# - extension: + b4*(x*y)^2;
### draft v.0.3k - v.0.3l:
# - started work on:
#   x^3 + b1*x*y^2 = R;
# - solved special case: b1*b2 = 1 & extension;
# - simple system:
#   x^2*y + b1*x*y^2 = R;
### draft v.0.3h - v.0.3j:
# - more Binomial Expansions:
#   various derivatives of Order 3;
### draft v.0.3g:
# - solved composed system:
#   x^4*y^2 + b2*x^3*y^3 + b1*x*y^2 = R1;
### draft v.0.3f:
# - solved Binomial-Derived system:
#   x^3 + 3*x*y^2 + b2*y^2 + b2*x*y + b1*y = R1;
### draft v.0.3c - v.0.3e:
# - solved Mixed Leading Term:
#   x^2*y + b1*x*y = R1;
# - solved variants:
#   x^2*y + a*x*y^2 + b1*x*y = R1; [v.0.3d]
#   x^2*y + b*y + b1*x*y = R1; [v.0.3e]
### draft v.0.3a - v.0.3b:
# - solved:
#   x^n + b1*x*y = R1;
# - extension:
#   x^n + b2*(x*y)^2 + b1*x*y = R1;
### draft v.0.2g:
# - comments on basic transforms;
### draft v.0.2f - v.0.2f-ext:
# - generalized approach to Order 2 Asymmetric:
#   b1*x^2 + b2*y^2 = R2;
# - A-type extension to the Order 3 system;
### draft v.0.2c - v.0.2e:
# - generalized approach to Order 1 Asymmetric:
#   b1*x + b2*y = R2;
# - solved Order 3 & Order 4 systems; [v.0.2d, v.0.2e]
### draft v.0.2b:
# - moved Mixt Symmetric variant to new file:
#   Poly.System.MixtVar.Hetero.Asym.R;
# - solved: Simple Asymetric system:
#   x^2 + b*y^2 = R1;
#   y^2 + x*y = R2;
### draft v.0.2a-sym:
# - some initial work on:
#   (x+d)^3 + y^3 = R;
# - the symmetric variants (in v.0.2a-sym);
### draft v.0.1d:
# - asymmetric Coefficients: Order 3
#   x^3 + b1*y = R;
# - TODO: special Case b1 + b2 == 0;
### draft v.0.1c:
# - solved variant with asymmetric Coefficients:
#   x^2 + b1*y = R;
### draft v.0.1b - v.0.1b-dual:
# - solved Order 3:
#   x^3 + b*y = R1;
# - added explicitly the Dual system; [v.0.1b-dual]
### draft v.0.1a:
# - solved Order 2:
#   x^2 + b*y = R1;
# - special cases: when b1[S] == 0;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R;
# e.g. round0(), round0.p;

### Other

order.cm = function(x, digits=5) {
	if(is.null(dim(x)) || length(dim(x)) == 1) {
		id = order(round(Re(x), digits), round(Im(x), digits));
	} else {
		id = order(
			round(Re(x[,1]), digits), round(Im(x[,1]), digits),
			round(Re(x[,2]), digits), round(Im(x[,2]), digits));
	}
	return(id);
}


##########################

##########################
### Polynomial Systems ###
##########################

###############
### Simple  ###
###############

# x^n + b1*x*y = R1
# y^n + b2*x*y = R2

### Solution:

### Step 1:
# - solve for (x*y);
# x^n = R1 - b1*x*y
# y^n = R2 - b2*x*y
### Prod =>
(x*y)^n - b1*b2*(x*y)^2 + (b1*R2+b2*R1)*(x*y) - R1*R2 # = 0

### Step 2:
# - solve for (x+y); [robust]
x^n + y^n + (b1+b2)*x*y - R1 - R2 # = 0

### Examples:

###############
### Order 3 ###
###############

# x^3 + b1*x*y = R1
# y^3 + b2*x*y = R2

### Solution

### Prod =>
(x*y)^3 - b1*b2*(x*y)^2 + (b1*R2+b2*R1)*(x*y) - R1*R2 # = 0

### Step 2:
x^3 + y^3 + (b1+b2)*x*y - R1 - R2 # = 0
S^3 - 3*(x*y)*S + (b1+b2)*x*y - R1 - R2 # = 0

### Solver
solve.Simplxy.S2P3 = function(R, b, debug=TRUE) {
	coeff = c(1, -b[1]*b[2], (b[1]*R[2]+b[2]*R[1]), - R[1]*R[2])
	xy = roots(coeff);
	if(debug) print(xy);
	# S
	S = sapply(xy, function(xy) roots(c(1, 0, - 3*xy, (b[1]+b[2])*xy - (R[1]+R[2]))));
	xy = rep(xy, each=3);
	xy.diff = (R[1] - R[2] - (b[1] - b[2])*xy) / (S^2 - xy);
	x = (S + xy.diff) / 2;
	y = S - x;
	sol = cbind(x=as.vector(x), y=as.vector(y));
}

### Examples:
R = c(1, 2)
b = c(-1, 3)

sol = solve.Simplxy.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3 + b[1]*x*y # - R[1]
y^3 + b[2]*x*y # - R[2]


#########
### Ex 2:
R = c(3, -2)
b = c(-1, 2)

sol = solve.Simplxy.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3 + b[1]*x*y # - R[1]
y^3 + b[2]*x*y # - R[2]


################
### Extended ###
################

# x^3 + b12*(x*y)^2 + b11*x*y = R1
# y^3 + b22*(x*y)^2 + b21*x*y = R2

### Solution

### Prod =>
b12*b22*(x*y)^4 + (b12*b21+b11*b22 - 1)*(x*y)^3 + (b11*b21 - b12*R2 - b22*R1)*(x*y)^2 +
	- (b11*R2+b21*R1)*(x*y) + R1*R2 # = 0

### Step 2:
x^3 + y^3 + (b12+b22)*(x*y)^2 + (b11+b21)*x*y - R1 - R2 # = 0
S^3 - 3*(x*y)*S + (b12+b22)*(x*y)^2 + (b11+b21)*x*y - R1 - R2 # = 0


### Solver
solve.xy.S2P3 = function(R, b, debug=TRUE) {
	b11 = b[1]; b12 = b[2]; b21 = b[3]; b22 = b[4];
	coeff = c(b12*b22, (b12*b21+b11*b22 - 1), (b11*b21 - b12*R[2] - b22*R[1]),
		- (b11*R[2] + b21*R[1]), R[1]*R[2])
	xy = roots(coeff);
	if(debug) print(xy);
	# S
	S = sapply(xy, function(xy) roots(
		c(1, 0, - 3*xy, (b12+b22)*(xy)^2 + (b11+b21)*xy - (R[1]+R[2]))));
	xy = rep(xy, each=3);
	# x - y
	x3 = R[1] - b12*xy^2 - b11*xy;
	y3 = R[2] - b22*xy^2 - b21*xy;
	xy.diff = (x3 - y3) / (S^2 - xy);
	x = (S + xy.diff) / 2;
	y = S - x;
	sol = cbind(x=as.vector(x), y=as.vector(y));
}

### Examples:
R = c(1, 2)
b = c(-1, 3, 2, 2)

sol = solve.xy.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3 + b[2]*(x*y)^2 + b[1]*x*y # - R[1]
y^3 + b[4]*(x*y)^2 + b[3]*x*y # - R[2]


# degenerate polynomial:
round0.p(poly.calc(x)) * 3


##########################

##########################
### Mixed Leading Term ###
##########################

# x^2*y + b1*x*y = R1
# y^2*x + b2*x*y = R2

### Solution

### Prod =>
(x*y)^3 - b1*b2*(x*y)^2 + (b1*R2+b2*R1)*(x*y) - R1*R2 # = 0

### Step 2:
x^2*y + x*y^2 + (b1+b2)*x*y - R1 - R2 # = 0
x*y*S + (b1+b2)*x*y - R1 - R2 # = 0


### Solver
solve.MLxy.S2P3 = function(R, b, debug=TRUE) {
	coeff = c(1, -b[1]*b[2], (b[1]*R[2]+b[2]*R[1]), - R[1]*R[2])
	xy = roots(coeff);
	if(debug) print(xy);
	# S
	S = (R[1] + R[2] - (b[1]+b[2])*xy) / xy;
	xy.diff = (R[1] - R[2] - (b[1] - b[2])*xy) / xy;
	x = (S + xy.diff) / 2;
	y = S - x;
	sol = cbind(x=as.vector(x), y=as.vector(y));
}

### Examples:
R = c(-1, 2)
b = c(-1, 3)

sol = solve.MLxy.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2*y + b[1]*x*y # - R[1]
y^2*x + b[2]*x*y # - R[2]

# simple P[3];


##########################

### Extensions

# x^2*y + a*x*y^2 + b1*x*y = R1
# y^2*x + a*x^2*y + b2*x*y = R2

### Solution

### Sum:
(a+1)*x*y*(x+y) + (b1+b2)*x*y - R1 - R2 # = 0
(a+1)*x*y*S + (b1+b2)*x*y - R1 - R2 # = 0

### Prod =>
(a+1)*(x*y)^3 + a*(x*y)^2*(x^2+y^2) - b1*b2*(x*y)^2 + (b1*R2+b2*R1)*(x*y) - R1*R2 # = 0
(a+1)*(x*y)^3 + a*(x*y)^2*(S^2 - 2*x*y) - b1*b2*(x*y)^2 + (b1*R2+b2*R1)*(x*y) - R1*R2 # = 0
(1-a)*(x*y)^3 + a*(x*y)^2*S^2 - b1*b2*(x*y)^2 + (b1*R2+b2*R1)*(x*y) - R1*R2
(a+1)^2*(1-a)*(x*y)^3 + a*((b1+b2)*x*y - R1 - R2)^2 +
	- (a+1)^2*b1*b2*(x*y)^2 + (a+1)^2*(b1*R2+b2*R1)*(x*y) - (a+1)^2*R1*R2
(a+1)^2*(1-a)*(x*y)^3 + a*(b1+b2)^2*(x*y)^2 - (a+1)^2*b1*b2*(x*y)^2 +
	+ (a+1)^2*(b1*R2+b2*R1)*(x*y) - 2*a*(b1+b2)*(R1 + R2)*x*y + a*(R1+R2)^2 - (a+1)^2*R1*R2


### Solver
solve.MLxy.S2P3 = function(R, b, a, debug=TRUE) {
	bs = b[1] + b[2]; Rs = R[1] + R[2];
	coeff = c((a+1)^2*(1-a), a*(bs)^2 - (a+1)^2*b[1]*b[2],
		(a+1)^2*(b[1]*R[2]+b[2]*R[1]) - 2*a*bs*Rs, a*Rs^2 - (a+1)^2*R[1]*R[2])
	xy = roots(coeff);
	if(debug) print(xy);
	# S
	S = (R[1] + R[2] - (b[1]+b[2])*xy) / xy / (a+1); # TODO: a = -1;
	xy.diff = (R[1] - R[2] - (b[1] - b[2])*xy) / xy / (1-a);
	x = (S + xy.diff) / 2;
	y = S - x;
	sol = cbind(x=as.vector(x), y=as.vector(y));
}

### Examples:
R = c(-1, 2)
b = c(-1, 3)
a = 2

sol = solve.MLxy.S2P3(R, b, a=a)
x = sol[,1]; y = sol[,2];

### Test
x^2*y + a*x*y^2 + b[1]*x*y # - R[1]
y^2*x + a*x^2*y + b[2]*x*y # - R[2]

##########################

#################
### Extension ###

# x^2*y + b*y + b1*x*y = R1
# y^2*x + b*x + b2*x*y = R2

### Solution

### Sum:
x^2*y + x*y^2 + b*S + (b1+b2)*x*y - R1 - R2 # = 0
(x*y + b)*S + (b1+b2)*x*y - R1 - R2 # = 0
# (x*y + b)*S = R1 + R2 - (b1+b2)*x*y;

### Prod =>
(x*y)^3 + b*x*y*(x^2+y^2) + b^2*x*y - b1*b2*(x*y)^2 + (b1*R2+b2*R1)*(x*y) - R1*R2 # = 0
(x*y)^3 + b*x*y*(S^2 - 2*x*y) - b1*b2*(x*y)^2 + (b1*R2+b2*R1 + b^2)*(x*y) - R1*R2 # = 0
(x*y)^3 + b*x*y*S^2 - (2*b + b1*b2)*(x*y)^2 + (b1*R2+b2*R1 + b^2)*(x*y) - R1*R2
### Eq:
(- R1*R2*b^2) +
((R1*b2 + R2*b1)*b^2 + (R1^2 + R2^2)*b + b^4)*x*y +
- (R1*R2 + 2*(R1*b1 + R2*b2)*b + b^2*b1*b2)*(x*y)^2 +
(R1*b2 + R2*b1 + b*(b1^2 + b2^2) - 2*b^2)*(x*y)^3 +
(- b1*b2)*(x*y)^4 + (x*y)^5


### Solver
solve.MLxy.S2P3 = function(R, b, ba, debug=TRUE) {
	coeff = c(1, -ba[1]*ba[2], (ba[1]*R[2]+ba[2]*R[1] + b*(ba[1]^2 + ba[2]^2) - 2*b^2),
		- (R[1]*R[2] + 2*(R[1]*ba[1] + R[2]*ba[2])*b + b^2*ba[1]*ba[2]),
		((R[1]*ba[2] + R[2]*ba[1])*b^2 + (R[1]^2 + R[2]^2)*b + b^4),
		- R[1]*R[2]*b^2)
	xy = roots(coeff);
	if(debug) print(xy);
	# S
	S = (R[1] + R[2] - (ba[1]+ba[2])*xy) / (xy + b);
	xy.diff = (R[1] - R[2] - (ba[1] - ba[2])*xy) / (xy - b);
	x = (S + xy.diff) / 2;
	y = S - x;
	sol = cbind(x=as.vector(x), y=as.vector(y));
}

### Examples:
R = c(-1, 2)
b = 2;
ba = c(-1, 3)

sol = solve.MLxy.S2P3(R, b=b, ba=ba)
x = sol[,1]; y = sol[,2];

### Test
x^2*y + b*y + ba[1]*x*y # - R[1]
y^2*x + b*x + ba[2]*x*y # - R[2]

### simple P[5]
poly.calc(x)


##########################

### Reducible

# x^2*y + b11*x + b12*x*y = R1
# y^2*x + b21*y + b22*x*y = R2

# trivial system: P[3];

### Solution

### Prod: x^2*y + b11*x = R1 - b12*x*y
(x*y)^3 + (b11+b21-b12*b22)*(x*y)^2 + (b11*b21+b12*R2+b22*R1)*x*y - R1*R2 # = 0

### =>
x*(x*y + b11) + b12*x*y - R1 # = 0


############
### P3+1 ###
############

# x^3*y + b11*x^2 + b12*x*y = R1
# y^3*x + b21*y^2 + b22*x*y = R2

# still relatively trivial: (x, y), (-x, -y);

### Prod: x^3*y + b11*x^2 = R1 - b12*x*y
(x*y)^4 + (b11+b21)*(x*y)^3 + (b11*b21-b12*b22)*(x*y)^2 + (b12*R2+b22*R1)*x*y - R1*R2 # = 0


############
### P3+1 ###
############

# x^3*y + b12*x^2 + b*x = R
# y^3*x + b22*y^2 + b*y = R

### Solution:

### Prod: x^3*y + b12*x^2 = R - b*x
(x*y)^4 + (b12+b22)*(x*y)^3 + b12*b22*(x*y)^2 - b^2*x*y + b*R*S - R^2 # = 0

### Prod: x^3*y + b*x - R = - b12*x^2
(x*y)^4 + b*(x*y)^2*S - R*x*y*(S^2 - 2*x*y) - b12*b22*(x*y)^2 + b^2*x*y - b*R*S + R^2 # = 0

### Sum: Eq 1b + 2b =>
2*(x*y)^4 + (b12+b22)*(x*y)^3 + b*(x*y)^2*S - R*x*y*(S^2 - 2*x*y) # = 0
2*(x*y)^3 + (b12+b22)*(x*y)^2 + b*(x*y)*S - R*(S^2 - 2*x*y) # = 0
2*b^2*R*(x*y)^3 + b^2*R*(b12+b22)*(x*y)^2 + b^2*R*b*(x*y)*S - b^2*R^2*S^2 + 2*b^2*R^2*x*y
2*b^2*R*(x*y)^3 + b^2*R*(b12+b22)*(x*y)^2 +
	- b^2*(x*y)*((x*y)^4 + (b12+b22)*(x*y)^3 + b12*b22*(x*y)^2 - b^2*x*y - R^2) +
	- ((x*y)^4 + (b12+b22)*(x*y)^3 + b12*b22*(x*y)^2 - b^2*x*y - R^2)^2 + 2*b^2*R^2*x*y
b^2*(x*y)^5 + b^2*(b12+b22)*(x*y)^4 - b^2*(2*R - b12*b22)*(x*y)^3 +
	- b^2*R*(b12+b22)*(x*y)^2 - b^4*(x*y)^2 - 3*b^2*R^2*x*y +
	+ ((x*y)^4 + (b12+b22)*(x*y)^3 + b12*b22*(x*y)^2 - b^2*x*y - R^2)^2

### Diff(y^2*...) =>
(b12-b22)*x^2*y^2 - b*x*y*(x-y) + R*(x-y)*S # = 0

### Sum(y^2*...) =>
2*x^3*y^3 + (b12+b22)*x^2*y^2 + b*x*y*S - R*(S^2 - 2*x*y) # = 0 # same as Sum above;

### TODO!

###
R = 2
b = 3
bi = c(3, 2)

### Debug
x = -1.0901825396;
y = -1.3159531274;
S = x+y; b12 = bi[1]; b22 = bi[2];


###############
### P2+1    ###
### Special ###
###############

# x^2*y + b1*x^2 - R/b2*x = R
# y^2*x + b2*y^2 - R/b1*y = R

### Solution:

### Prod: - (x^2*y - R) = b1*x^2 - R/b2*x
(x*y)^3 - b1*b2*(x*y)^2 - R^2/(b1*b2)*x*y + R^2 # = 0
(x*y - b1*b2)*(b1*b2*(x*y)^2 - R^2) # = 0

### * y =>
# - b1*(x*y)*x + R*y =
(x*y)^2 - R/b2*(x*y)
# - b2*(x*y)*y + R*x =
(x*y)^2 - R/b1*(x*y)

### Solver:
solve.pAsym.S2P21 = function(R, b, debug=TRUE) {
	bpr = b[1]*b[2];
	coeff = c(bpr, 0, - R^2)
	xy2 = roots(coeff); xy = c(bpr, xy2);
	if(debug) print(xy);
	xy.s1 = (xy)^2 - R[1]/b[2]*(xy);
	xy.s2 = (xy)^2 - R[1]/b[1]*(xy);
	### Root 1:
	xy = xy[1];
	div = R[1]^2 - b[1]*b[2]*(xy)^2;
	x = (b[2]*(xy)*xy.s1 + R*xy.s2)[1] / div;
	y = (b[1]*(xy)*xy.s2 + R*xy.s1)[1] / div;
	### other Roots:
	xy = xy2; # xy.s1 = xy.s1[-1]; xy.s2 = xy.s2[-1];
	x2 = sapply(xy, function(xy) roots(c(b[1], xy - R[1]/b[2], - R[1])));
	xy = rep(xy, each=2); y2 = xy / x2;
	x = c(x, x2); y = c(y, y2);
	sol = cbind(x=as.vector(x), y=as.vector(y))
	return(sol);
}

### Examples:

R = 1
b = c(-1, 3)
sol = solve.pAsym.S2P21(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2*y + b[1]*x^2 - R/b[2]*x # - R
y^2*x + b[2]*y^2 - R/b[1]*y # - R

### Poly:
poly.calc(x[-1]) * 9


### Ex 2:
R = 3
b = c(-1, 3)
sol = solve.pAsym.S2P21(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2*y + b[1]*x^2 - R/b[2]*x # - R
y^2*x + b[2]*y^2 - R/b[1]*y # - R

### Poly:
poly.calc(x[-1])


##########################
##########################
##########################

###############
### Order 2 ###
### Simple  ###

# x^2 + y^2 = R1
# b1*x + b2*y = R2

### Eq 2:
(b1*x + b2*y)*(b2*x + b1*y) - R2*(b2*x + b1*y) # = 0
b1*b2*(x^2 + y^2) + (b1^2 + b2^2)*x*y - R2*(b2*x + b1*y) # = 0
### - R2 * Eq 2 =>
b1*b2*(x^2 + y^2) + (b1^2 + b2^2)*x*y - R2*(b1 + b2)*(x + y) + R2^2 # = 0
### generalized Eq:
b1*b2*(x^2 + y^2) + (b1^2 + b2^2)*x*y - R2*(b1 + b2)*S + R2^2 # = 0
### =>
b1*b2*R1 + (b1^2 + b2^2)*x*y - R2*(b1 + b2)*S + R2^2 # = 0
# (b1^2 + b2^2)*x*y = R2*(b1 + b2)*S - b1*b2*R1 - R2^2

### Eq 1 =>
S^2 - 2*x*y - R1 # = 0
(b1^2 + b2^2)*S^2 - 2*(b1^2 + b2^2)*x*y - (b1^2 + b2^2)*R1 # = 0
(b1^2 + b2^2)*S^2 - 2*(R2*(b1 + b2)*S - b1*b2*R1 - R2^2) - (b1^2 + b2^2)*R1 # = 0
### Eq:
(b1^2 + b2^2)*S^2 - 2*R2*(b1 + b2)*S - (b1^2 + b2^2 - 2*b1*b2)*R1 + R2^3 # = 0


### Solver:
solve.AsymSimple.P2 = function(R, b, debug=TRUE) {
	bs = b[1] + b[2]; bsq = b[1]^2 + b[2]^2;
	coeff = c(bsq, - 2*R[2]*bs, - (bsq - 2*b[1]*b[2])*R[1] + R[2]^3)
	S = roots(coeff);
	if(debug) print(S);
	#
	bd = b[2] - b[1];
	x = (b[2]*S - R[2]) / bd;
	y = - (b[1]*S - R[2]) / bd;
	cbind(x=as.vector(x), y=as.vector(y))
}

### Examples:

R = c(-1, 2)
b = c(-1, 3)
sol = solve.AsymSimple.P2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2 + y^2 # - R[1]
b[1]*x + b[2]*y # - R[2]


###############
### Order 3 ###
### Simple  ###

# x^3 + y^3 = R1
# b1*x + b2*y = R2

### Solution:

### Eq 2 =>
b1*b2*(S^2 - 2*x*y) + (b1^2 + b2^2)*x*y - R2*(b1 + b2)*S + R2^2 # = 0
b1*b2*S^2 + (b1^2 + b2^2 - 2*b1*b2)*x*y - R2*(b1 + b2)*S + R2^2 # = 0
# (b1^2 + b2^2 - 2*b1*b2)*x*y = - b1*b2*S^2 + R2*(b1 + b2)*S - R2^2

### Eq 1 =>
S^3 - 3*x*y*S - R1 # = 0
(b1 - b2)^2*S^3 - 3*(b1 - b2)^2*x*y*S - (b1 - b2)^2*R1 # = 0
(b1 - b2)^2*S^3 + 3*(b1*b2*S^2 - R2*(b1 + b2)*S + R2^2)*S - (b1 - b2)^2*R1 # = 0
(b1^2 + b2^2 + b1*b2)*S^3 - 3*R2*(b1 + b2)*S^2 + 3*R2^2*S - (b1 - b2)^2*R1 # = 0


### Solver:
solve.AsymSimple.P3 = function(R, b, debug=TRUE) {
	bs = b[1] + b[2]; bsq = b[1]^2 + b[2]^2;
	coeff = c((bsq + b[1]*b[2]), - 3*R[2]*bs, 3*R[2]^2, - (b[1] - b[2])^2*R[1])
	S = roots(coeff);
	if(debug) print(S);
	#
	bd = b[2] - b[1];
	x = (b[2]*S - R[2]) / bd;
	y = - (b[1]*S - R[2]) / bd;
	cbind(x=as.vector(x), y=as.vector(y))
}

### Examples:

R = c(-1, 2)
b = c(-1, 3)
sol = solve.AsymSimple.P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3 + y^3 # - R[1]
b[1]*x + b[2]*y # - R[2]

poly.calc(x) * 7


###############
### Order 4 ###
### Simple  ###

# x^4 + y^4 = R1
# b1*x + b2*y = R2

### Solution:

### Eq 2 =>
b1*b2*S^2 + (b1^2 + b2^2 - 2*b1*b2)*x*y - R2*(b1 + b2)*S + R2^2 # = 0
# (b1 - b2)^2*x*y = - b1*b2*S^2 + R2*(b1 + b2)*S - R2^2

### Eq 1 =>
S^4 - 4*x*y*S^2 + 2*(x*y)^2 - R1 # = 0
(b1 - b2)^4*S^4 - 4*(b1 - b2)^4*x*y*S^2 + 2*(b1-b2)^4*(x*y)^2 - (b1-b2)^4*R1 # = 0
(b1 - b2)^4*S^4 + 4*(b1 - b2)^2*(b1*b2*S^2 - R2*(b1 + b2)*S + R2^2)*S^2 +
	+ 2*(b1*b2*S^2 - R2*(b1 + b2)*S + R2^2)^2 - (b1 - b2)^4*R1 # = 0
(b1-b2)^2*(b1+b2)^2*S^4 - 4*(b1 - b2)^2*(b1 + b2)*R2*S^3 + 4*(b1 - b2)^2*R2^2*S^2 +
	+ 2*(b1*b2*S^2 - R2*(b1 + b2)*S + R2^2)^2 - (b1 - b2)^4*R1 # = 0
(b1^4 + b2^4)*S^4 - 4*(b1^2 + b2^2 - b1*b2)*(b1+b2)*R2*S^3 +
	+ 6*(b1^2 + b2^2)*R2^2*S^2 - 4*(b1 + b2)*R2^3*S + 2*R2^4 - (b1 - b2)^4*R1 # = 0


### Solver:
solve.AsymSimple.P4 = function(R, b, debug=TRUE) {
	bs = b[1] + b[2]; bp = b[1]*b[2];
	bsq = bs^2 - 2*bp;
	coeff = c(bsq^2 - 2*bp^2, - 4*(bsq - bp)*bs*R[2],
		6*bsq*R[2]^2, - 4*bs*R[2]^3, 2*R[2]^4 - (b[1] - b[2])^4*R[1])
	S = roots(coeff);
	if(debug) print(S);
	#
	bd = b[2] - b[1];
	x = (b[2]*S - R[2]) / bd;
	y = - (b[1]*S - R[2]) / bd;
	cbind(x=as.vector(x), y=as.vector(y))
}

### Examples:

R = c(-1, 2)
b = c(-1, 3)
sol = solve.AsymSimple.P4(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^4 + y^4 # - R[1]
b[1]*x + b[2]*y # - R[2]

poly.calc(x) * 82


##################
##################

##################
### Asymmetry: ###
### Order 2    ###
##################

### b1*x^2 + b2*y^2 = R2

### =>
(b1*x^2 + b2*y^2)*(b2*x^2 + b1*y^2) - R2*(b2*x^2 + b1*y^2) # = 0
b1*b2*(x^4 + y^4) + (b1^2 + b2^2)*x^2*y^2 - R2*(b2*x^2 + b1*y^2) # = 0
### + R2*Eq2 =>
b1*b2*(x^4 + y^4) + (b1^2 + b2^2)*x^2*y^2 - R2*(b1 + b2)*(x^2 + y^2) + R2^2 # = 0
b1*b2*(S^4 - 4*x*y*S^2 + 2*(x*y)^2) + (b1^2 + b2^2)*x^2*y^2 - R2*(b1 + b2)*(S^2 - 2*x*y) + R2^2 # = 0
b1*b2*S^4 + (b1 + b2)^2*x^2*y^2 - 4*b1*b2*x*y*S^2 + 2*R2*(b1 + b2)*x*y - R2*(b1 + b2)*S^2 + R2^2 # = 0

### Auxiliary Eq:
### *(b2*x + b1*y) =>
b1*b2*(x^3 + y^3) + x*y*(b1^2*x + b2^2*y) - R2*(b2*x + b1*y) # = 0
b1*b2*(x^3 + y^3) + (x*y*b1^2 - b2*R2)*x + (x*y*b2^2 - b1*R2)*y # = 0


###############
### Order 3 ###
###############

### x^3 + y^3 = R1
### b1*x^2 + b2*y^2 = R2

### Solution:

### Eq 1:
S^3 - 3*x*y*S - R1 # = 0
# 3*x*y*S = S^3 - R1

### Eq 2 =>
9*b1*b2*S^6 + 9*(b1 + b2)^2*x^2*y^2*S^2 - 4*9*b1*b2*x*y*S^4 +
	+ 2*9*R2*(b1 + b2)*x*y*S^2 - 9*R2*(b1 + b2)*S^4 + 9*R2^2*S^2 # = 0
9*b1*b2*S^6 + (b1 + b2)^2*(S^3 - R1)^2 - 12*b1*b2*(S^3 - R1)*S^3 +
	+ 6*R2*(b1 + b2)*(S^3 - R1)*S - 9*R2*(b1 + b2)*S^4 + 9*R2^2*S^2 # = 0
	
### Eq:
((b1 + b2)^2 - 3*b1*b2)*S^6 - 3*R2*(b1 + b2)*S^4 - 2*((b1 + b2)^2 - 6*b1*b2)*R1*S^3 +
	+ 9*R2^2*S^2 - 6*(b1 + b2)*R1*R2*S + (b1 + b2)^2*R1^2 # = 0


### Solver:
solve.AsymSimple.P3A2 = function(R, b, b.ext=c(0,0), debug=TRUE) {
	bs = b[1] + b[2]; bp = b[1]*b[2];
	bsq = bs^2;
	coeff = c((bsq - 3*bp), 0, - 3*bs*R[2], - 2*(bsq - 6*bp)*R[1],
		9*R[2]^2, - 6*bs*R[1]*R[2], bsq*R[1]^2)
	if(length(b.ext) < 2) b.ext = c(b.ext, 0)
	if(any(b.ext != 0)) {
		coeff = coeff + c(0, 3*b.ext[2]*bs, 2*(bsq - 6*bp)*b.ext[1] + 9*b.ext[2]^2,
			-18*R[2]*b.ext[2] - 6*bs*b.ext[1]*b.ext[2], 6*bs*(R[1]*b.ext[2] + R[2]*b.ext[1]) + bsq*b.ext[1]^2,
			-2*bsq*R[1]*b.ext[1], 0)
	}
	S = roots(coeff);
	if(debug) print(S);
	S = S[ ! (S == 0)]
	#
	R1 = R[1] - b.ext[1]*S; R2 = R[2] - b.ext[2]*S;
	xy = (S^3 - R1) / (3*S);
	# b1*b2*(x^3 + y^3) + (x*y*b1^2 - b2*R2)*x + (x*y*b2^2 - b1*R2)*y # = 0
	T0 = - b[1]*b[2]*R1; bd = b[1] - b[2];
	diff = - xy*bs*bd - R2*bd;
	x = ((xy*b[2]^2 - b[1]*R2)*S - T0) / diff;
	y = - ((xy*b[1]^2 - b[2]*R2)*S - T0) / diff;
	cbind(x=as.vector(x), y=as.vector(y))
}

### Examples:

R = c(-1, 2)
b = c(-1, 3)
sol = solve.AsymSimple.P3A2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3 + y^3 # - R[1]
b[1]*x^2 + b[2]*y^2 # - R[2]


round0.p(poly.calc(x)) * 13
round0.p(poly.calc(y)) * 13
round0.p(poly.calc(x+y)) * 13


### Extensions:

R = c(-1, 2)
b = c(-1, 3)
b.ext = c(2, -3)
sol = solve.AsymSimple.P3A2(R, b, b.ext=b.ext)
x = sol[,1]; y = sol[,2];

### Test
x^3 + y^3 + b.ext[1]*(x+y) # - R[1]
b[1]*x^2 + b[2]*y^2 + b.ext[2]*(x+y) # - R[2]


round0.p(poly.calc(x)) * 13
round0.p(poly.calc(y)) * 13
round0.p(poly.calc(x + y)) * 13


##################
##################

##################
### Transforms ###

### Simple Transforms:
# - simple linear transforms:
#   e.g. combinations of initial roots;
# - other simple transforms;

### Initial:
# x^2 + y^2 = R1
# x*y = R2

### Transformed:
# x^2 + b*y^2 = Rf1
# y^2 + x*y = Rf2

### Transform:
# x => x + k*y
# y => y

### Eq 1:
(x + k*y)^2 + y^2 - R1 # = 0
x^2 + 2*k*x*y + (k^2 + 1)*y^2 - R1

### Eq 2:
(x + k*y)*y - R2 # = 0
x*y + k*y^2 - R2

### Eq 1 - 2*k*Eq2 =>
x^2 + (1 - k^2)*y^2 - R1 + 2*k*R2 # = 0
k*y^2 + x*y - R2 # = 0

### y => y/k
x^2 + (1 - k^2)/k^2*y^2 - R1 + 2*k*R2 # = 0
y^2 + x*y - k*R2 # = 0

### Solver:
solve.AsymSimple.P2 = function(R, b.y) {
	k = sqrt(1/(b.y + 1) + 0i)
	R2 = R[2] / k;
	R1 = R[1] + 2*R[2];
	S = sqrt(R1 + 2*R2 + 0i)
	S = c(S, -S)
	xy.d = sqrt(R1 - 2*R2 + 0i)
	x0 = (S + xy.d) / 2;
	y0 = (S - xy.d) / 2;
	x = c(x0, y0)
	y = c(y0, x0)
	y = y * k;
	x = x - y;
	cbind(x=as.vector(x), y=as.vector(y))
}
test.AsymSimple.P2 = function(R, b.y) {
	err1 = x^2 + b.y*y^2 # - R[1]
	err2 = y^2 + x*y # - R[2]
	err = rbind(err1, err2)
	round0(err)
}

### Examples:
R = c(1, 3)
b.y = 2
sol = solve.AsymSimple.P2(R, b.y)
x = sol[,1]; y = sol[,2];

### Test
test.AsymSimple.P2(R, b.y)


### Ex 2:
R = c(0, 3)
b.y = 4
sol = solve.AsymSimple.P2(R, b.y)
x = sol[,1]; y = sol[,2];

### Test
test.AsymSimple.P2(R, b.y)


### Ex 3:
R = c(1, -1)
b.y = 1
sol = solve.AsymSimple.P2(R, b.y)
x = sol[,1]; y = sol[,2];

### Test
test.AsymSimple.P2(R, b.y)


########################
########################

##############
### Hetero ###
##############

###############
### Order 2 ###
### Partial ###

# x^2 + b*y = R1
# y^2 + b*x = R2

### Sum =>
x^2 + y^2 + b*(x+y) - R1 - R2 # = 0
S^2 - 2*x*y + b*S - R1 - R2
# 2*x*y = S^2 + b*S - R1 - R2

### Sum(x*...) =>
x^3 + y^3 + 2*b*x*y - R1*x - R2*y # = 0
x^3 + y^3 + b*(S^2 + b*S - R1 - R2) - R1*x - R2*y # = 0
S^3 - 3*x*y*S + b*(S^2 + b*S - R1 - R2) - R1*x - R2*y # = 0
2*S^3 - 3*(S^2 + b*S - R1 - R2)*S + 2*b*(S^2 + b*S - R1 - R2) - 2*R1*x - 2*R2*y # = 0
-S^3 - b*S^2 + (3*R1 + 3*R2 + 2*b^2)*S - 2*b*(R1 + R2) - 2*R1*x - 2*R2*y # = 0
# 2*R1*x + 2*R2*y = -S^3 - b*S^2 + (3*R1 + 3*R2 + 2*b^2)*S - 2*b*(R1 + R2)

### Sum(y*...) =>
x^2*y + y^2*x + b*(x^2+y^2) - R1*y - R2*x # = 0
x*y*(x+y) + b*(S^2 - 2*x*y) - R1*y - R2*x # = 0
x*y*S + b*(S^2 - (S^2 + b*S - R1 - R2)) - R1*y - R2*x # = 0
x*y*S - b*(b*S - R1 - R2) - R1*y - R2*x # = 0
(S^2 + b*S - R1 - R2)*S - 2*b*(b*S - R1 - R2) - 2*R1*y - 2*R2*x # = 0
S^3 + b*S^2 - (R1 + R2 + 2*b^2)*S + 2*b*(R1 + R2) - 2*R1*y - 2*R2*x # = 0
# 2*R2*x + 2*R1*y = S^3 + b*S^2 - (R1 + R2 + 2*b^2)*S + 2*b*(R1 + R2)

### Eq 2 +/- 3:
# 2*(R1^2 - R2^2)*x =
R1 * (-S^3 - b*S^2 + (3*R1 + 3*R2 + 2*b^2)*S - 2*b*(R1 + R2)) +
	- R2 * (S^3 + b*S^2 - (R1 + R2 + 2*b^2)*S + 2*b*(R1 + R2)) # =
-(R1+R2)*S^3 - b*(R1+R2)*S^2 + (3*R1^2 + R2^2 + 4*R1*R2 + 2*b^2*R1+ 2*b^2*R2)*S - 2*b*(R1 + R2)^2
# 2*(R1^2 - R2^2)*y =
R1*(S^3 + b*S^2 - (R1 + R2 + 2*b^2)*S + 2*b*(R1 + R2)) +
	- R2*(-S^3 - b*S^2 + (3*R1 + 3*R2 + 2*b^2)*S - 2*b*(R1 + R2)) # =
(R1+R2)*S^3 + b*(R1+R2)*S^2 - (R1^2 + 3*R2^2 + 4*R1*R2 + 2*b^2*R1 + 2*b^2*R2)*S + 2*b*(R1 + R2)^2

### Diff =>
(x-y)*(x+y) - b*(x-y) # = R1 - R2
(x-y)*S - b*(x-y) - R1 + R2 # = 0
(y - x)*(S - b) + R1 - R2 # = 0
2*(R1^2 - R2^2)*(y - x)*(S - b) + 2*(R1^3 - R1^2*R2 - R1*R2^2 + R2^3) # = 0
# =>
((R1+R2)*S^3 + b*(R1+R2)*S^2 - 2*(R1^2 + R2^2 + 2*R1*R2 + b^2*R1 + b^2*R2)*S + 2*b*(R1 + R2)^2)*(S-b) +
	+ (R1^3 - R1^2*R2 - R1*R2^2 + R2^3) # = 0
(R1+R2)*S^4 + b*(R1+R2)*S^3 - 2*(R1^2 + R2^2 + 2*R1*R2 + b^2*R1 + b^2*R2)*S^2 + 2*b*(R1 + R2)^2*S +
	- b*((R1+R2)*S^3 + b*(R1+R2)*S^2 - 2*(R1^2 + R2^2 + 2*R1*R2 + b^2*R1 + b^2*R2)*S + 2*b*(R1 + R2)^2) +
	+ (R1^3 - R1^2*R2 - R1*R2^2 + R2^3) # = 0
(R1+R2)*S^4 - (R1+R2)*(2*R1 + 2*R2 + 3*b^2)*S^2 +
	+ 2*b*(R1+R2)*(2*R1 + 2*R2 + b^2)*S +
	+ (R1+R2)*(R1^2+R2^2 - 2*R1*R2) - 2*b^2*(R1+R2)^2 # = 0
S^4 - (2*R1 + 2*R2 + 3*b^2)*S^2 + 2*b*(2*R1 + 2*R2 + b^2)*S +
	+ (R1^2+R2^2 - 2*R1*R2) - 2*b^2*(R1+R2) # = 0

### Diff(x*...) =>
x^3 - y^3 - R1*x + R2*y # = 0
(x-y)*(S^2 - x*y) - R1*x + R2*y # = 0
(S-b)*(x-y)*(S^2 - x*y) - (R1*x - R2*y)*(S - b) # = 0
(R2 - R1)*(S^2 - b*S + R1 + R2) + 2*(R1*x - R2*y)*(S - b) # = 0


### TODO:
# - solve efficiently?
# - evaluate higher powers;


### Solution:
solve.asym.S2P2 = function(R, b) {
	R12 = R[1] + R[2]
	coeff = c(1, 0, - (2*R[1] + 2*R[2] + 3*b[1]^2), 2*b[1]*(2*R[1] + 2*R[2] + b[1]^2),
		(R[1]^2+R[2]^2 - 2*R[1]*R[2]) - 2*b[1]^2*(R[1]+R[2]))
	S = roots(coeff)
	len = length(S)
	#
	div = 2*(R[1] - R[2])*R12
	x = -R12*S^3 - b*R12*S^2 + (2*R12^2 + R12*(R[1] - R[2]) + 2*b^2*R12)*S - 2*b*R12^2
	x = x / div;
	y = R12*S^3 + b*R12*S^2 - (2*R12^2 - R12*(R[1] - R[2]) + 2*b^2*R12)*S + 2*b*R12^2
	y = y / div;
	return(cbind(x=x, y=y))
}

### Examples:
R = c(1,2)
b = 1
sol = solve.asym.S2P2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2 + b[1]*y
y^2 + b[1]*x


#########
### Ex 2:
# special Case:
# S^4 - 8*S^2 + 32;
R = c(1, -3)
b = 2
sol = solve.asym.S2P2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2 + b[1]*y
y^2 + b[1]*x

### Classic Polynomial:
round0.p(poly.calc(x))
round0.p(poly.calc(x + y))

err = 13 + 8*x - 2*x^2 + x^4
round0(err)


#########
### Ex 3:
# special Case: b0 == 0
# S^4 - 16*S^2 + 32*S;
R = c(3, -1)
b = 2
sol = solve.asym.S2P2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2 + b[1]*y
y^2 + b[1]*x

### Classic Polynomial:
round0.p(poly.calc(x))
round0.p(poly.calc(x + y))

err = 13 + 8*x - 6*x^2 + x^4
round0(err)


#########
### Ex 4:
# special Case: b0 == 0
# S^4 - 31*S^2 + 78*S;
R = c(4, -2)
b = 3
sol = solve.asym.S2P2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2 + b[1]*y
y^2 + b[1]*x

### Classic Polynomial:
round0.p(poly.calc(x))
round0.p(poly.calc(x + y))

err = 34 + 27*x - 8*x^2 + x^4
round0(err)


###########

### Debug
R = c(2,3)
b = 3
x = 0.9563705005
y = 0.3617851553
S = x+y; R1 = R[1]; R2 = R[2];

### Test
x^2 + b[1]*y
y^2 + b[1]*x


######################

###############
### Order 3 ###

# x^3 + b*y = R1
# y^3 + b*x = R2

### Sum =>
x^3 + y^3 + b*(x+y) - R1 - R2 # = 0
S^3 - 3*x*y*S + b*S - R1 - R2
# 3*x*y*S = S^3 + b*S - R1 - R2

### Sum(x*...) =>
x^4 + y^4 + 2*b*x*y - R1*x - R2*y # = 0
S^4 - 4*x*y*S^2 + 2*(x*y)^2 + 2*b*x*y - R1*x - R2*y # = 0
# R1*x + R2*y = ...;
S^6 - 4*x*y*S^4 + 2*(x*y)^2*S^2 + 2*b*x*y*S^2 - R1*x*S^2 - R2*y*S^2 # = 0
-3*S^6 - 6*b*S^4 + 6*(2*R1 + 2*R2)*S^3 + 2*(S^3 + b*S - R1 - R2)^2 +
	+ 6*b^2*S^2 - 6*b*(R1 + R2)*S - 9*R1*x*S^2 - 9*R2*y*S^2 # = 0
-S^6 - 2*b*S^4 + (8*R1 + 8*R2)*S^3 + 8*b^2*S^2 - 10*b*(R1 + R2)*S +
	- 9*R1*x*S^2 - 9*R2*y*S^2 + 2*(R1+R2)^2 # = 0
# 9*(R1*x + R2*y)*S^2 =
	-S^6 - 2*b*S^4 + 8*(R1+R2)*S^3 + 8*b^2*S^2 - 10*b*(R1+R2)*S + 2*(R1+R2)^2

### Sum(y*...) =>
x*y*(x^2+y^2) + b*(x^2+y^2) - R2*x - R1*y # = 0
x*y*S^2 - 2*(x*y)^2 + b*S^2 - 2*b*x*y - R2*x - R1*y # = 0
# R2*x + R1*y = ...;

### Liniar (x, y) =>
### (R1^2 - R2^2) * x =
R1*(S^4 - 4*x*y*S^2 + 2*(x*y)^2 + 2*b*x*y) +
	- R2*(x*y*S^2 - 2*(x*y)^2 + b*S^2 - 2*b*x*y)
R1*S^4 - (4*R1+R2)*x*y*S^2 - b*R2*S^2 + 2*(R1+R2)*(x*y)^2 + 2*b*(R1+R2)*x*y
### -(R1^2 - R2^2) * y =
R2*(S^4 - 4*x*y*S^2 + 2*(x*y)^2 + 2*b*x*y) +
	- R1*(x*y*S^2 - 2*(x*y)^2 + b*S^2 - 2*b*x*y)
R2*S^4 - (4*R2+R1)*x*y*S^2 - b*R1*S^2 + 2*(R1+R2)*(x*y)^2 + 2*b*(R1+R2)*x*y

### =>
# (R1^2 - R2^2)*(x+y) =
(R1-R2)*S^4 - 3*(R1-R2)*x*y*S^2 + b*(R1-R2)*S^2
# =>
(R1-R2)*S^4 - 3*(R1-R2)*x*y*S^2 + b*(R1-R2)*S^2 - (R1^2 - R2^2)*S # = 0
S^3 - 3*x*y*S + b*S - (R1 + R2) # = 0 # cyclic redundancy!
# (R1^2 - R2^2)*(x-y) =
(R1+R2)*S^4 - 5*(R1+R2)*x*y*S^2 - b*(R1+R2)*S^2 + 4*(R1+R2)*(x*y)^2 + 4*b*(R1+R2)*x*y
# (R1 - R2)*(x-y) =
S^4 - 5*x*y*S^2 - b*S^2 + 4*(x*y)^2 + 4*b*x*y


### Diff(x*...) =>
x^3 - y^3 - b*(x - y) - R1 + R2 # = 0
(x-y)*(x^2 + y^2 + x*y) - b*(x - y) - R1 + R2 # = 0
d * (S^2 - x*y - b) - R1 + R2
d * (3*S^3 - (S^3 + b*S - R1 - R2) - 3*b*S) - (3*R1 - 3*R2)*S
d * (2*S^3 - 4*b*S + R1 + R2) - (3*R1 - 3*R2)*S


### Diff(x*...) =>
x^4 - y^4 - R1*x + R2*y # = 0
(x-y)*S*(x^2 + y^2) - R1*x + R2*y
d*S*(S^2 - 2*x*y) - R1*x + R2*y
d*S*(S^3 - 2*b*S + 2*R1 + 2*R2) - 3*R1*x*S + 3*R2*y*S
# =>
d*S*(2*S^3 - 4*b*S + 4*R1 + 4*R2) - 6*R1*x*S + 6*R2*y*S
(R1 - R2)*S^2 + (R1 + R2)*d*S - (2*R1*x - 2*R2*y)*S
# (2*R1*x - 2*R2*y)*S =
	(R1 - R2)*S^2 + (R1 + R2)*d*S
### x =>
# 36*R1*S^2*x # =
-2*S^6 - 4*b*S^4 + (25*R1 + 7*R2)*S^3 + 9*(R1 + R2)*d*S^2 +
	+ 16*b^2*S^2 - 20*b*(R1 + R2)*S + 4*(R1+R2)^2
-2*S^6 - 4*b*S^4 + (25*R1 + 7*R2)*S^3 + 27*(R1 + R2)*(R1 - R2)*S^3 / (2*S^3 - 4*b*S + R1 + R2) +
	+ 16*b^2*S^2 - 20*b*(R1 + R2)*S + 4*(R1+R2)^2
### y =>
# 36*R2*S^2*y # =
-2*S^6 - 4*b*S^4 + (7*R1 + 25*R2)*S^3 - 9*(R1 + R2)*d*S^2 +
	+ 16*b^2*S^2 - 20*b*(R1 + R2)*S + 4*(R1+R2)^2
-2*S^6 - 4*b*S^4 + (7*R1 + 25*R2)*S^3 - 27*(R1 + R2)*(R1 - R2)*S^3 / (2*S^3 - 4*b*S + R1 + R2) +
	+ 16*b^2*S^2 - 20*b*(R1 + R2)*S + 4*(R1+R2)^2

### Sum x + y =>
# 36*R1*R2*S^2*(x+y)*(2*S^3 - 4*b*S + R1 + R2) # =
(2*S^3 - 4*b*S + R1 + R2)*(R1 + R2)*(-2*S^6 - 4*b*S^4 + 16*b^2*S^2 - 20*b*(R1 + R2)*S + 4*(R1+R2)^2) +
	+ (2*S^3 - 4*b*S + R1 + R2)*(50*R1*R2 + 7*R1^2 + 7*R2^2)*S^3 +
	- 27*(R1 + R2)*(R1 - R2)^2*S^3
(2*S^3 - 4*b*S + R1 + R2)*(R1 + R2)*(-2*S^6 - 4*b*S^4 + 16*b^2*S^2 - 20*b*(R1 + R2)*S + 4*(R1+R2)^2) +
	+ (2*S^3 - 4*b*S)*(50*R1*R2 + 7*R1^2 + 7*R2^2)*S^3 +
	- (R1 + R2)*(20*R1^2 + 20*R2^2 - 104*R1*R2)*S^3
# =>
(2*S^3 - 4*b*S + R1 + R2)*(R1 + R2)*(-2*S^6 - 4*b*S^4 + 16*b^2*S^2 - 20*b*(R1 + R2)*S + 4*(R1+R2)^2) +
	+ (2*S^3 - 4*b*S)*(50*R1*R2 + 7*R1^2 + 7*R2^2)*S^3 - 36*R1*R2*(2*S^3 - 4*b*S + R1 + R2)*S^3+
	- (R1 + R2)*(20*R1^2 + 20*R2^2 - 104*R1*R2)*S^3 # = 0
(2*S^3 - 4*b*S + R1 + R2)*(-S^6 - 2*b*S^4 + 8*b^2*S^2 - 10*b*(R1 + R2)*S + 2*(R1+R2)^2) +
	+ 7*(S^3 - 2*b*S)*(R1 + R2)*S^3 +
	- (10*R1^2 + 10*R2^2 - 34*R1*R2)*S^3 # = 0

### Prod =>
(x*y)^3 + b*(x^4 + y^4) + b^2*x*y - R1*R2 # = 0
(x*y)^3 + b*S^4 - 4*b*x*y*S^2 + 2*b*(x*y)^2 + b^2*x*y - R1*R2 # = 0
# 3*x*y*S = S^3 + b*S - R1 - R2 =>
# S^4 = 3*x*y*S^2 - b*S^2 + (R1 + R2)*S
(x*y)^3 - b*x*y*S^2 - b^2*S^2 + b*(R1 + R2)*S + 2*b*(x*y)^2 + b^2*x*y - R1*R2 # = 0

### Dual System: (S, x*y)
S^3 - 3*x*y*S + b*S - R1 - R2 # = 0
b*S^4 + (x*y)^3 - 4*b*x*y*S^2 + 2*b*(x*y)^2 + b^2*x*y - R1*R2 # = 0
# Q: How is the Dual System helpful?

### Eq:
S^9 - 3*(R1 + R2)*S^6 - 12*b^2*S^5 + 18*b*(R1 + R2)*S^4 +
	+ (3*(R1 + R2)^2 - 27*R1*R2 + 16*b^3)*S^3 - 24*b^2*(R1 + R2)*S^2 +
	+ 9*b*(R1 + R2)^2*S - (R1 + R2)^3

#############

### Solution:
solve.asym.S2P3 = function(R, b) {
	R12 = R[1] + R[2]
	coeff = c(1, 0, 0, - 3*R12, - 12*b[1]^2, 18*b[1]*R12,
		(3*R12^2 - 27*R[1]*R[2] + 16*b[1]^3), - 24*b[1]^2*R12,
		9*b[1]*R12^2, - R12^3)
	S = roots(coeff)
	len = length(S)
	#
	div = (2*S^3 - 4*b[1]*S + R12);
	### TODO: div == 0
	# d = 3*(R[1] - R[2])*S / div;
	if(round0(R12) == 0) {
		isZero = round0(S == 0)
		S = S[ ! isZero]
		d = 3*(R[1] - R[2])*S / div[ ! isZero]
		d0 = roots(c(1,0,-b[1], -R[1]))
		d = c(d, 2*d0)
		S = c(S, rep(0, length(d0)))
	} else {
		d = ifelse(round0(div) == 0,
			NA, # TODO
			3*(R[1] - R[2])*S / div)
	}
	x = (S + d)/2;
	y = (S - d)/2;
	return(cbind(x=x, y=y))
}

### Examples:
R = c(1,2)
b = 1
sol = solve.asym.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3 + b[1]*y
y^3 + b[1]*x


#########
### Ex 2:
# special case: R1 + R2 == 0
# Note: can be decomposed more efficiently using a different approach;
R = c(1, -1)
b = 1
sol = solve.asym.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3 + b[1]*y
y^3 + b[1]*x


#########
### Ex 3:
# special case: b3 == 0 & R1 + R2 == 0
R = c(4, -4)
b = -3
sol = solve.asym.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3 + b[1]*y
y^3 + b[1]*x


#########
### Ex 3:
R = c(1, 2)
b = -3
sol = solve.asym.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3 + b[1]*y
y^3 + b[1]*x

round0.p(poly.calc(x))
x = sol[,1] # Note: x^1!
err = -55 - 81*x + 3*x^3 - 3*x^6 + x^9
round0(err)
#
round0.p(poly.calc(x + y))
x = sol[,1] + sol[,2]
err = -27 - 243*x - 648*x^2 - 459*x^3 - 162*x^4 - 108*x^5 - 9*x^6 + x^9
round0(err)


#########
### Ex 4:
R = c(-1, 2)
b = -3
sol = solve.asym.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3 + b[1]*y
y^3 + b[1]*x

round0.p(poly.calc(x))
x = sol[,1] # Note: x^1!
err = -53 - 81*x + 3*x^3 + 3*x^6 + x^9
round0(err)
#
round0.p(poly.calc(y))
x = sol[,2] # Note: x^1!
err = 19 - 81*x + 12*x^3 - 6*x^6 + x^9
round0(err)
# special relationship between polynomials:
round0.p(poly.calc(x + y))
x = sol[,1] + sol[,2]
err = -1 - 27*x - 216*x^2 - 375*x^3 - 54*x^4 - 108*x^5 - 3*x^6 + x^9
round0(err)


###########

### Debug
R = c(2,3)
b = 3
x = 0.9860634751
y = 0.3470765347
S = x+y; d = x - y; R1 = R[1]; R2 = R[2];

### Test
x^3 + b[1]*y
y^3 + b[1]*x

### alternative Eq:
(- R1*R2*b^6 - R1^3*R2^3 + R1^4*b^3 + R2^4*b^3) +
(- R1^2*R2^2*b^2 + b^8)*E2 +
(- R1*R2*b^4)*E2^2 +
(3*R1^2*R2^2 - 4*b^6)*E2^3 +
(5*R1*R2*b^2)*E2^4 +
(6*b^4)*E2^5 +
(- 3*R1*R2)*E2^6 +
(- 4*b^2)*E2^7 + E2^9 


#######################

################
### Extended ###
################

# x^n + b*y + bm1*x*y = R1
# y^n + b*x + bm2*x*y = R2

### Solution:

### Prod =>
# x^n + b*y = R1 - bm1*x*y
# y^n + b*x = R2 - bm2*x*y
### Prod =>
(x*y)^n + b*(x^(n+1) + y^(n+1)) + (b^2 - bm1*bm2)*(x*y)^2 + (bm1*R2+bm2*R1)*(x*y) - R1*R2 # = 0

### Sum =>
x^n + y^n + b*(x+y) + (bm1+bm2)*x*y - R1 - R2 # = 0

###############
### Order 3 ###
###############

# x^3 + b*y + bm1*x*y = R1
# y^3 + b*x + bm2*x*y = R2

### Solution:

### Sum =>
x^3 + y^3 + b*S + (bm1+bm2)*x*y - R1 - R2 # = 0
S^3 - 3*x*y*S + b*S + (bm1+bm2)*x*y - R1 - R2 # = 0
# (3*S - (bm1+bm2))*x*y =
S^3 + b*S - R1 - R2;

### Prod =>
(x*y)^3 + b*(x^4 + y^4) + b^2*x*y - bm1*bm2*(x*y)^2 + (bm1*R2+bm2*R1)*(x*y) - R1*R2 # = 0
(x*y)^3 + b*S^4 - 4*b*x*y*S + (2*b - bm1*bm2)*(x*y)^2 + (b^2 + bm1*R2+bm2*R1)*(x*y) - R1*R2

### TODO


#######################
#######################


#######################
#######################

### Existence of Symmetrization

### x^n + b1*x*y^(n-1) = R
### y^n + b2*x^(n-1)*y = R

###############
### Order 3 ###
###############

# x^3 + b1*x*y^2 = R
# y^3 + b2*x^2*y = R

### Solution:

### Sum(x*...) =>
x^4 + y^4 + (b1+b2)*x^2*y^2 - R*(x+y) # = 0
S^4 - 4*x*y*S^2 + (b1+b2+2)*(x*y)^2 - R*S # = 0

### Prod:
# x^3 - R = -b1*x*y^2 # Prod =>
(x*y)^3 - R*(x^3 + y^3) + R^2 - b1*b2*(x*y)^3 # = 0
(x*y)^3 - R*(S^3 - 3*x*y*S) + R^2 - b1*b2*(x*y)^3 # = 0
(x*y)^3 - b1*b2*(x*y)^3 - R*S^3 + 3*R*x*y*S + R^2 # = 0

### TODO: solve;


### Case: b1*b2 = 1
(x^3 + y^3) - R # = 0
S^3 - 3*x*y*S - R # = 0
# 3*x*y*S = S^3 - R
### =>
9*S^6 - 36*x*y*S^4 + 9*(b1+b2+2)*(x*y)^2*S^2 - 9*R*S^3 # = 0
9*S^6 - 12*(S^3 - R)*S^3 + (b1+b2+2)*(S^3 - R)^2 - 9*R*S^3 # = 0
(b1+b2-1)*S^6 - (2*b1+2*b2+1)*R*S^3 + (b1+b2+2)*R^2 # = 0
### alternative:
(S^3 - R + R)*S - 4*x*y*S^2 + (b1+b2+2)*(x*y)^2 - R*S
S^2 - (b1+b2+2)*(x*y) # = 0
3*S^3 - (b1+b2+2)*(S^3 - R) # = 0
(b1+b2-1)*S^3 - (b1+b2+2)*R # = 0

### Solver:
solve.special.S2P3Asym = function(R, b, b.ext=0, debug=TRUE) {
	if(length(b) == 1) {
		b[2] = 1 / b[1];
	} else if(round0(b[1]*b[2] - 1) != 0) stop("Only special case implemented!");
	coeff = c((b[1]+b[2]-1), 0, (b[1]+b[2]+2)*b.ext[1], - (b[1]+b[2]+2)*R)
	S = roots(coeff);
	if(debug) print(S);
	R1 = R[1] - b.ext[1]*S;
	xy = (S^3 - R1) / (3*S);
	isZero = round0(xy) == 0;
	if(any(isZero)) print("Invalid cases: x*y == 0!")
	xy = xy[ ! isZero]; S = S[ ! isZero]; R1 = R1[ ! isZero];
	xy.sb = (2*R1 - (S^3 - 3*xy*S)) / xy;
	x = (b[1]*S - xy.sb) / (b[1] - b[2]);
	y = S - x;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	return(sol);
}

### Examples:

### Special Case:
R = 2
b = 2
#
sol = solve.special.S2P3Asym(R, b);
x = sol[,1]; y = sol[,2];

### Test
x^3 + b[1]*x*y^2 # - R
y^3 + 1/b[1]*x^2*y # - R

###########
### Special: Ex 2
R = 2
b = 2
b.ext = 3
#
sol = solve.special.S2P3Asym(R, b, b.ext=b.ext);
x = sol[,1]; y = sol[,2];

### Test
x^3 + b[1]*x*y^2 + b.ext[1]*(x+y) # - R
y^3 + 1/b[1]*x^2*y + b.ext[1]*(x+y) # - R

#######################

### Variants:

# x^2*y + b1*x*y^2 = R
# y^2*x + b2*x^2*y = R


### Solution:
# - method using linear algebra: NOT generalizable;

### Prod:
# x^2*y - R = -b1*x*y^2
(x*y)^3 - R*x*y*S + R^2 - b1*b2*(x*y)^3 # = 0
### Case: b1*b2 = 1
# - has NO solutions!
x*y*S - R # = 0
# x*y*S = R;

### Diff =>
(1-b2)*x + (b1-1)*y # = 0
# =>
(b2-1)*x^2*S - (b1-1)*R # = 0
(b2-1)*(b1+b2-2)*x^3 - (b1-1)^2*R # = 0

### Sum(x*...) =>
x*y*(x^2 + y^2) + (b1+b2)*x^2*y^2 - R*(x+y) # = 0
x*y*(S^2 - 2*x*y) + (b1+b2)*x^2*y^2 - R*S # = 0
x*y*S^2 + (b1+b2-2)*x^2*y^2 - R*S # = 0
x*y*S*S^3 + (b1+b2-2)*x^2*y^2*S^2 - R*S^3 # = 0
(b1+b2-2)*R^2 # = 0
# NO solution: ?? except for extensions ??

### Solver:
solve.special.S2P3Asym = function(R, b, b.ext=1, debug=TRUE) {
	if(length(b) == 1) {
		b[2] = 1 / b[1];
	} else if(round0(b[1]*b[2] - 1) != 0) stop("Only special case implemented!");
	if(all(b.ext == 0)) stop("NO solutions!");
	S = (b[1]-1)*R[1] / (b.ext[1]*(b[1] - 1));
	x = (b[1] - 1) * S / (b[1]+b[2]-2);
	y = S - x;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	return(sol);
}

### Examples:

### Special Case:
R = 2
b = 2
b.ext = 1;
#
sol = solve.special.S2P3Asym(R, b, b.ext = b.ext);
x = sol[,1]; y = sol[,2];

### Test
x^2*y + b[1]*x*y^2 + b.ext[1]*(x+y) # - R
x*y^2 + 1/b[1]*x^2*y + b.ext[1]*(x+y) # - R


#######################
#######################

######################
### Cross-Products ###
######################

# x^2*y + b1*x = R
# y^2*x + b2*y = R

### Variant 1:
# x^2*y + b3*x*y + b1*x = R
# y^2*x + b3*x*y + b2*y = R

### Variant 2:
# x^2*y + b4*(x*y)^2 + b3*x*y + b1*x = R
# y^2*x + b4*(x*y)^2 + b3*x*y + b2*y = R

### Sum(y*...) =>
2*(x*y)^2 + (b1+b2)*(x*y) - R*S # = 0
# for Variants:
2*(x*y)^2 + (b1+b2)*(x*y) - (R - b3*x*y)*S # = 0 # V1
2*(x*y)^2 + (b1+b2)*(x*y) - (R - b4*(x*y)^2 - b3*x*y)*S # = 0 # V2

### Prod =>
(x*y)^3 + (b1+b2)*(x*y)^2 + b1*b2*(x*y) - R^2 # = 0
# for Variants:
# x^2*y + b1*x = R - b3*x*y - b4*(x*y)^2
(x*y)^3 + (b1+b2-b3^2)*(x*y)^2 + (b1*b2+2*b3*R)*(x*y) - R^2 # = 0 # V1
b4^2*(x*y)^4 + (2*b3*b4 - 1)*(x*y)^3 - (b1+b2-b3^2+2*b4*R)*(x*y)^2 +
	- (b1*b2+2*b3*R)*(x*y) + R^2 # = 0 # V2

### Solver:
solve.Pr.S2P21 = function(R, b, debug=TRUE) {
	b.s = b[1] + b[2];
	len = length(b)
	if(len < 4) b = c(b, rep(0, 4 - len));
	if(len < 4) {
		coeff = c(1, b.s - b[3]^2, b[1]*b[2] + 2*b[3]*R, - R[1]^2);
	} else {
		coeff = c(b[4]^2, (2*b[3]*b[4] - 1), - (b.s - b[3]^2 + 2*b[4]*R[1]),
			- (b[1]*b[2] + 2*b[3]*R), R[1]^2);
	}
	xy = roots(coeff);
	# TODO: Case R[1] - b[3]*xy == 0;
	bxy = b[3]*xy + b[4]*xy^2;
	S = (2*(xy)^2 + b.s*xy) / (R[1] - bxy);
	if(debug) print(xy);
	if(debug) print(S);
	# sum =>
	xy.sb = 2*R[1] - xy*S - 2*bxy;
	x = (b[2]*S - xy.sb) / (b[2] - b[1]);
	y = S - x;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	return(sol);
}

### Examples:
R = -1
b = c(2,3)
#
sol = solve.Pr.S2P21(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2*y + b[1]*x # - R
y^2*x + b[2]*y # - R


#########
### Ex 2:
R = -1
b = c(2,3, -1)
#
sol = solve.Pr.S2P21(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2*y + b[3]*x*y + b[1]*x # - R
y^2*x + b[3]*x*y + b[2]*y # - R


#########
### Ex 3:
R = -1
b = c(2,3, -1, -1)
#
sol = solve.Pr.S2P21(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2*y + b[4]*(x*y)^2 + b[3]*x*y + b[1]*x # - R
y^2*x + b[4]*(x*y)^2 + b[3]*x*y + b[2]*y # - R


#######################

### Cross-Products

################
### Order 3+ ###
################

# x^3*y + b3*(x*y)^2 + b1*x^2 = R
# y^3*x + b3*(x*y)^2 + b2*y^2 = R

### Solution:

### Prod =>
(x*y)^4 - b3^2*(x*y)^4 + (b1+b2)*(x*y)^3 + (b1*b2 + 2*b3*R)*(x*y)^2 - R^2 # = 0

### Sum(y^2*...) =>
2*(x*y)^3 + b3*(x*y)^2*(x^2+y^2) + (b1+b2)*(x*y)^2 - R*(x^2+y^2) # = 0
2*(x*y)^3 + b3*(x*y)^2*(S^2 - 2*x*y) + (b1+b2)*(x*y)^2 - R*(S^2 - 2*x*y) # = 0
(b3*(x*y)^2 - R)*S^2 + 2*(x*y)^3 - 2*b3*(x*y)^3 + (b1+b2)*(x*y)^2 + 2*R*x*y # = 0


### Solver:
solve.Pr.S2P31 = function(R, b, debug=TRUE) {
	b.s = b[1] + b[2];
	len = length(b)
	coeff = c(1 - b[3]^2, b.s, (b[1]*b[2] + 2*b[3]*R[1]), 0, - R[1]^2);
	xy = roots(coeff);
	S = sapply(xy, function(xy) {
			coeff = c((b[3]*(xy)^2 - R[1]), 0, 2*(1 - b[3])*(xy)^3 + b.s*(xy)^2 + 2*R[1]*xy);
			roots(coeff);
		});
	if(debug) print(xy);
	if(debug) print(S);
	len = 2; # length(S) / length(xy);
	xy = rep(xy, each=len);
	# robust:
	xb3 = b[3]*(xy)^2;
	x2 = (R[1] - xb3) / (xy + b[1]);
	y2 = S^2 - x2 - 2*xy;
	xy.diff = (x2 - y2) / S;
	x = (S + xy.diff)/2;
	y = S - x;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	return(sol);
}

### Examples:
R = -1
b = c(2,3, 2)
#
sol = solve.Pr.S2P31(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3*y + b[3]*(x*y)^2 + b[1]*x^2 # - R
y^3*x + b[3]*(x*y)^2 + b[2]*y^2 # - R

### degenerate P[8]
round0.p(poly.calc(x)) * 27
round0.p(poly.calc(y)) * 19*3


#######################
#######################

### Generalized

# a11*x^3*y + a12*x*y^3 + b12*(x*y)^2 + b11*(x*y) = R1
# a21*x^3*y + a22*x*y^3 + b22*(x*y)^2 + b21*(x*y) = R2

### Solution:

### Prod (x^3*y)*(x*y^3)
(x*y)^4 - (r10 + ra12*(x*y)^2 + ra11*(x*y))*(r20 + ra22*(x*y)^2 + ra21*(x*y))
(x*y)^4 - ra12*ra22*(x*y)^4 - (ra11*ra22 + ra12*ra21)*(x*y)^3 +
	- (ra11*ra21 + r10*ra22 + r20*ra12)*(x*y)^2 +
	- (r10*ra21 + r20*ra11)*(x*y) - r10*r20;

### Sum =>
x^3*y + x*y^3 - (r10 + r12*(x*y)^2 + r11*(x*y)) - (r20 + r22*(x*y)^2 + r21*(x*y)) # = 0
x*y*(x^2+y^2) - (r10+r20 + (r12+r22)*(x*y)^2 + (r11+r21)*(x*y))

### Solver:
solve.Gen.S2P31 = function(R, b, a, debug=TRUE) {
	coeff.m = solve(a, cbind(R, -b));
	r1 = coeff.m[1,]; r10 = r1[1]; r11 = r1[2]; r12 = r1[3];
	r2 = coeff.m[2,]; r20 = r2[1]; r21 = r2[2]; r22 = r2[3];
	coeff = c(1 - r12*r22, - (r11*r22 + r12*r21),
		- (r11*r21 + r10*r22 + r20*r12),
		- (r10*r21 + r20*r11), - r10*r20);
	xy = roots(coeff);
	if(debug) print(xy);
	# TODO: x*y = 0;
	S2 = (r10+r20 + (r12+r22 + 2)*(xy)^2 + (r11+r21)*(xy)) / xy;
	S = sqrt(S2 + 0i); S = c(S, -S);
	xy = c(xy, xy);
	# TODO: robust;
	xy.diff = sqrt(S^2 - 4*xy)
	x = (S + xy.diff) / 2
	y = S - x;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	return(sol);
}

### Examples:
R = c(-1, 2)
a = matrix(c(2,3, 1, -2), nrow=2, byrow=T)
b = matrix(c(1,2, -2,3), nrow=2, byrow=T)
#
sol = solve.Gen.S2P31(R, b=b, a=a)
x = sol[,1]; y = sol[,2];

### Test
a[1,1]*x^3*y + a[1,2]*x*y^3 + b[1,2]*(x*y)^2 + b[1,1]*x*y # - R[1]
a[2,1]*x^3*y + a[2,2]*x*y^3 + b[2,2]*(x*y)^2 + b[2,1]*x*y # - R[2]


#######################

### Order 3 Mixt/Mixt-Term

# x^3*y^3 + b1*x^2*y + b2*x*y^2 = R1
# x^2*y + b3*x*y^2 = R2

### Solution

# alternative solution:
# X = x^2*y; Y = x*y^2;
# solve for (X, Y);

### {1,b3}*Eq 1 - {b1, b2}*Eq 2 =>
x^3*y^3 + (b2-b1*b3)*x*y^2 = R1 - b1*R2
b3*x^3*y^3 - (b2-b1*b3)*x^2*y = b3*R1 - b2*R2

### Re-order & Prod =>
b3*(x*y)^6 - (2*b3*R1 - (b1*b3+b2)*R2)*(x*y)^3 + (b2-b1*b3)^2*(x*y)^3  +
	+ (R1 - b1*R2)*(b3*R1 - b2*R2) # = 0

### Solver:
solve.Asym.S2P33 = function(R, b, debug=TRUE) {
	coeff = c(b[3], 0, 0, - (2*b[3]*R[1] - (b[1]*b[3]+b[2])*R[2]) + (b[2]-b[1]*b[3])^2,
		0, 0, (R[1] - b[1]*R[2])*(b[3]*R[1] - b[2]*R[2]))
	xy = roots(coeff);
	if(debug) print(xy);
	xy.sb2 = (R[1] - xy^3) / xy;
	bdiv = 1;
	if(b[1] != 0) {
		xy.sb2 = xy.sb2 / b[1];
		xy.sb3 = R[2] / xy;
		xy.sb2 = xy.sb2 - xy.sb3;
		bdiv = b[2]/b[1] - b[3];
	}
	y = xy.sb2 / bdiv;
	x = xy / y;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	return(sol);
}

### Examples:

R = c(-1, 3)
b = c(2,-1,2)
sol = solve.Asym.S2P33(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3*y^3 + b[1]*x^2*y + b[2]*x*y^2 # - R[1]
x^2*y + b[3]*x*y^2 # - R[2]

# trivial P[6] Poly
round0.p(poly.calc(x)) * 7

##################

### Generalization:
x^3*y^3 + b13*x^3 + b12*x^2*y + b11*x*y^2 = R1
x^3*y^3 + b23*y^3 + b22*x^2*y + b21*x*y^2 = R2

### Solution:

# - if (x, y) is a solution, then:
#   (x*m, y*m), (x*m^2, y*m^2) are also solutions;

### Y * ... =>
# - where: X = x^2*y, Y = x*y^2;
# - reduction of order;
X*Y^2 + b13*X^2 + b12*X*Y + b11*Y^2 - R1*Y # = 0
X^2*Y + b23*Y^2 + b22*X^2 + b21*X*Y - R2*X # = 0
# false solution: X = 0, Y = 0, but mat help to simplify system;


##########################
### Cross-Product Variant:

x^2*y*(x^3 + b1*x*y^2 + b10) = R1*(y^3 + b2*x^2*y + b20)
x*y^2*(y^3 + b2*x^2*y + b20) = R2*(x^3 + b1*x*y^2 + b10)

### System Decomposition
### Sub-Sys 1:
x^3 + b1*x*y^2 + b10 # = 0
y^3 + b2*x^2*y + b20 # = 0

### Sub-Sys 2:
(x*y)^3 - R1*R2 # = 0
# let: X = x^2*y, Y = x*y^2 =>
X*Y - R1*R2 # = 0
X^2*(X^2 + b1*Y^2 + b10*Y) - R1*Y*(Y^2 + b2*X^2 + b20*X) # = 0

### TODO


#######################
#######################

######################
### Cross-Products ###
######################

### Simple Order 3
# - decomposable

# x^2*(x + b) = R1*(y + b)
# y^2*(y + b) = R2*(x + b)

### Solution:

### Case 1:
# x = -b; y = -b;

### Case 2:
(x*y)^2 - R1*R2 # = 0

### *x =>
x^4 + b*x^3 - b*R1*x - R1*(x*y) # = 0


### Solver:
solve.Asym.S2P3 = function(R, b, debug=TRUE) {
	xy = sqrt(R[1]*R[2] + 0i);
	xy = c(xy, - xy);
	x = sapply(xy, function(xy) roots(c(1, b[1], 0, - b[1]*R[1], - R[1]*(xy))));
	xy = rep(xy, each=4);
	y = xy / x;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	return(sol);
}

### Examples:
R = c(2,-1)
b = -1
sol = solve.Asym.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test:
x^2*(x + b) - R[1]*y # = R[1]*b
y^2*(y + b) - R[2]*x # = R[2]*b


#######################

### Cross-Products

###############
### Order 4 ###

x^4 + b13*x^3 - R1/b23*y = R1
y^4 + b23*y^3 - R2/b13*x = R2

### Solution:

### Prod:
# b13*x^3 - R1/b23*y = R1 - x^4
b13*b23*(x*y)^3 + R1*R2/(b13*b23)*(x*y) - R1*R2 - (x*y)^4 # = 0
(x*y)^4 - b13*b23*(x*y)^3 - R1*R2/(b13*b23)*(x*y) + R1*R2 # = 0
# decomposable:
# Case 1: (x*y)^3 = R1*R2/(b1*b2);
# Case 2: x = -b1; y = -b2;

### Solver:
solver.Asym.S2P4 = function(R, b, debug=TRUE) {
	coeffs = c(b[1]*b[2], 0, 0, -R[1]*R[2]);
	xy = roots(coeffs);
	if(debug) print(xy);
	# x^5 + b1*x^4 - R1*x - R1/b2*(xy)
	x = sapply(xy, function(xy) roots(c(1, b[1], 0, 0, -R[1], -R[1]/b[2]*xy)));
	xy = rep(xy, each=5);
	y = xy / x;
	sol = cbind(x= c(-b[1], as.vector(x)), y= c(-b[2], as.vector(y)));
	return(sol)
}

### Examples:
R = c(-1, 3)
b = c(2, 3)
sol = solver.Asym.S2P4(R, b)
x = sol[,1]; y = sol[,2];

### Test:
x^4 + b[1]*x^3 - R[1]/b[2]*y # - R[1]
y^4 + b[2]*y^3 - R[2]/b[1]*x # - R[2]

#######################

### Full Cross-Products

x^2*(x^2 + b11*y + b10) = R1*(y^2 + b21*x + b20)
y^2*(y^2 + b21*x + b20) = R2*(x^2 + b11*y + b10)

### Solution:

### Case 1:
# (x*y)^2 = R1*R2;

### Case 2:
# x^2 + b11*y + b10 = 0
# y^2 + b21*x + b20 = 0
# - special case: b = symmetric, only R different;


#######################

### Cross-Product

x*(x^3 + y^3 + b1) = R1*(x*y + b2)
y*(x*y + b2) = R2*(x^3 + y^3 + b1)

### Solution:

### Decomposition =>
### Sys 1:
# x^3 + y^3 + b1 = 0
# x*y + b2 = 0

### Sys 2:
# x*y = R1*R2;
# =>
x^6 + b1*x^3 - R1*(x*y + b2)*x^2 + (R1*R2)^3 # = 0

### Solver:
solve.AsymDecomp.S2PS3 = function(R, b, debug=TRUE) {
	# Note: does NOT solve sub-system 1!
	xy = R[1]*R[2];
	coeff = c(1, 0, 0, b[1], - R[1]*(xy + b[2]), 0, (R[1]*R[2])^3);
	x = roots(coeff);
	y = xy / x;
	sol = cbind(x=x, y=y);
	return(sol);
}

### Examples:

R = c(-1, 3)
b = c(2, -1)
sol = solve.AsymDecomp.S2PS3(R, b)
x = sol[,1]; y = sol[,2];

### Test
round0( x*(x^3 + y^3 + b[1]) - R[1]*(x*y + b[2]) )
round0( y*(x*y + b[2]) - R[2]*(x^3 + y^3 + b[1]) )


#######################
#######################

##############
### Simple ###
##############

### x^n + b1*y = R
### y^n + b2*x = R

###############
### Order 2 ###
###############

# x^2 + b1*y = R
# y^2 + b2*x = R

### Sum(x*...) =>
x^3 + y^3 + (b1+b2)*x*y - R*(x+y) # = 0
S^3 - 3*x*y*S + (b1+b2)*x*y - R*S # = 0
# x*y*(3*S - (b1+b2)) = S^3 - R*S

### Prod:
# x^2 - R = -b1*y # Prod =>
(x*y)^2 - R*(x^2 + y^2) + R^2 - b1*b2*x*y # = 0
(x*y)^2 - R*(S^2 - 2*x*y) + R^2 - b1*b2*x*y # = 0
R*S^2 - (x*y)^2 - (2*R - b1*b2)*x*y - R^2 # = 0

### Dual System:
S^3 - 3*x*y*S + (b1+b2)*x*y - R*S # = 0
R*S^2 - (x*y)^2 - (2*R - b1*b2)*x*y - R^2 # = 0

### Auxilliary Eq:
### Sum =>
x^2 + y^2 + b2*x + b1*y - 2*R # = 0
S^2 + 2*x*y + b2*x + b1*y - 2*R # = 0
### Diff =>
x*(S - b[2]) - y*(S - b[1]) # = 0

### Eq:
S^4 - (4*R + 3*b1*b2)*S^2 + (b1 + b2)*(4*R + b1*b2)*S - R*(b1 + b2)^2

### Solution:
solve.asymCoeff.S2P2 = function(R, b) {
	coeff = c(1, 0, - (4*R[1] + 3*b[1]*b[2]),
		(b[1] + b[2])*(4*R + b[1]*b[2]), - R[1]*(b[1] + b[2])^2)
	S = roots(coeff)
	#
	# div = 3*S - (b[1]+b[2]);
	# isZero = round0(div) == 0;
	# if(isZero) print("Division by 0!")
	# xy = (S^3 - R[1]*S) / div
	x = S*(S - b[1]) / (2*S - b[1] - b[2])
	y = S*(S - b[2]) / (2*S - b[1] - b[2])
	return(cbind(x=x, y=y))
}

### Examples:
R = 1
b = c(1, 2)
#
sol = solve.asymCoeff.S2P2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2 + b[1]*y # - R
y^2 + b[2]*x # - R


#########
### Ex 2:
# special Case: b1[S] == 0
R = c(-3)
b = c(3, 4)
#
sol = solve.asymCoeff.S2P2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2 + b[1]*y # - R
y^2 + b[2]*x # - R

### Classic Polynomial:
round0.p(poly.calc(x))
# S-Polynomial:
round0.p(poly.calc(x+y))


#########
### Ex 3:
# special Case: b1[S] == 0
R = c(2)
b = c(-1, 8)
#
sol = solve.asymCoeff.S2P2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2 + b[1]*y # - R
y^2 + b[2]*x # - R

### Classic Polynomial:
round0.p(poly.calc(x))
err = 2 + 8*x - 4*x^2 + x^4
round0(err)
# S-Polynomial:
round0.p(poly.calc(x+y))


########################

###############
### Order 3 ###
###############

# x^3 + b1*y = R
# y^3 + b2*x = R

### Sum(x*...) =>
x^4 + y^4 + (b1+b2)*x*y - R*(x+y) # = 0
S^4 - 4*x*y*S^2 + 2*(x*y)^2 + (b1+b2)*x*y - R*S # = 0
# 2*(x*y)^2 - x*y*(4*S^2 - (b1+b2)) = - S^4 + R*S

### Prod:
# x^3 - R = -b1*y # Prod =>
(x*y)^3 - R*(x^3 + y^3) + R^2 - b1*b2*x*y # = 0
(x*y)^3 - R*(S^3 - 3*x*y*S) + R^2 - b1*b2*x*y # = 0
R*S^3 - (x*y)^3 - (3*R*S - b1*b2)*x*y - R^2 # = 0

### Dual System:
S^4 - 4*x*y*S^2 + 2*(x*y)^2 + (b1+b2)*x*y - R*S # = 0
R*S^3 - (x*y)^3 - (3*R*S - b1*b2)*x*y - R^2 # = 0

### Eq:
S^9 - 6*R*S^6 - 12*b1*b2*S^5 + 18*R*(b1 + b2)*S^4 - (15*R^2 - 8*b1*b2*(b1+b2))*S^3 +
	- R*(9*(b1 + b2)^2 + 12*b1*b2)*S^2 +
	+ (18*R^2*(b1 + b2) - b1*b2*(b1^2 - 2*b1*b2 + b2^2))*S +
	- R*(4*b1*b2*(b1 + b2) - (b1 + b2)^3 + 8*R^2)


### Auxilliary Eq:
### x*y =>
# xy = - R*S*b1 - R*S*b2 + 4*R^2 + S^4*b1 + S^4*b2 - 4*S^6
# div = 14*R*S - 8*S^2*b1 - 8*S^2*b2 + 14*S^4 - 2*b1*b2 + b1^2 + b2^2
# xy = - xy / div;
### Diff =>
x^3 - y^3 - b2*x + b1*y # = 0
(x - y)*(S^2 - x*y) - b2*x + b1*y
x*(S^2 - x*y - b2) - y*(S^2 - x*y - b1) # = 0
# x = S*(S^2 - xy - b1) / (2*S^2 - 2*xy - b1 - b2)
# y = S*(S^2 - xy - b2) / (2*S^2 - 2*xy - b1 - b2)


### Solution:
solve.asymCoeff.S2P3 = function(R, b) {
	b12 = b[1] + b[2]; bp = b[1]*b[2];
	coeff = c(1, 0, 0, - 6*R[1], - 12*bp, 18*R[1]*b12,
		(- 15*R[1]^2 + 8*bp*b12), - R[1]*(9*b12^2 + 12*bp),
		(18*R^2*b12 - bp*(b[1] - b[2])^2),
		- R*(4*bp*b12 - b12^3 + 8*R^2))
	S = roots(coeff)
	# x*y
	xy = - R[1]*S*b12 + 4*R^2 + S^4*b12 - 4*S^6
	div = 14*R*S - 8*S^2*b12 + 14*S^4 - 4*bp + b12^2
	xy = - xy / div
	# x*(S^2 - x*y - b2) - y*(S^2 - x*y - b1) # = 0
	x = S*(S^2 - xy - b[1]) / (2*S^2 - 2*xy - b12)
	y = S*(S^2 - xy - b[2]) / (2*S^2 - 2*xy - b12)
	return(cbind(x=x, y=y))
}

### Examples:

R = 1
b = c(2, 3)
sol = solve.asymCoeff.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3 + b[1]*y # - R
y^3 + b[2]*x # - R


#########
### Ex 2:
R = 1
b = c(-1, 1)
sol = solve.asymCoeff.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3 + b[1]*y # - R
y^3 + b[2]*x # - R


##################
### Special Cases:
### b1 + b2 == 0
b1 = b[1]
S^9 - 6*R*S^6 + 12*b1^2*S^5 - 15*R^2*S^3 + 12*b1^2*R*S^2 + 4*b1^4*S - 8*R^3
x^9 - 3*R*x^6 + 3*R^2*x^3 + b1^4*x + b1^3*R - R^3

### alternative Solution (special Case):
### Diff =>
x^3 - y^3 + b1*x + b1*y # = 0
(x - y)*(S^2 - x*y) + b1*S # = 0

### Sum =>
x^3 + y^3 - b1*(x - y) - 2*R # = 0
S^3 - 3*x*y*S - b1*(x - y) - 2*R # = 0
# =>
S^3*(S^2 - x*y) - 3*x*y*S*(S^2 - x*y) - b1*(x - y)*(S^2 - x*y) - 2*R*(S^2 - x*y) # = 0
S^5 - 4*x*y*S^3 - 2*R*S^2 + 3*(x*y)^2*S + b1^2*S + 2*R*x*y # = 0
### Dual System:
S^4 - 4*x*y*S^2 + 2*(x*y)^2 - R*S # = 0
S^5 - 4*x*y*S^3 + 3*(x*y)^2*S - 2*R*S^2 + b1^2*S + 2*R*x*y # = 0
### =>
S^9 - 6*R*S^6 + 12*b1^2*S^5 - 15*R^2*S^3 + 12*b1^2*R*S^2 + 4*b1^4*S - 8*R^3

### Sum(x*...) =>
x^4 + y^4 - R*(x+y) # = 0
S^4 - 4*x*y*S^2 + 2*(x*y)^2 - R*S # = 0
### Diff(x*...) =>
x^4 - y^4 + 2*b1*x*y - R*(x - y) # = 0
(x - y)*(S^3 - 2*x*y*S - R) + 2*b1*x*y # = 0
# - (x - y)*(S^3 - 3*x*y*S - 2*R - b1*(x - y)) =>
(x - y)*(x*y*S + R) + 2*b1*x*y + b1*(x - y)^2 # = 0
(x - y)*(x*y*S + R) + b1*(S^2 - 2*x*y) # = 0 => [redundant]


#########

### Debug
R = 1; b = c(2, 3);
b1 = b[1]; b2 = b[2];
x = 0.2947874543
y = 0.4871915377
S = x + y;


######################
######################

### Complex Transforms

# x^3 - 3*x*y^2 + b*x = R
# y^3 - 3*x^2*y + b*y = 0

### Base System:
# x^3 + b*y = R
# y^3 + b*x = R

### Solver:
solve.htC.S2P3 = function(R, b, debug=TRUE) {
	coeff = c(1, 0, -2*b[1], R[1])
	S = roots(coeff)
	if(debug) print(S); # Debug
	xy = S^2 - b[1];
	#
	x.diff = sqrt(S^2 - 4*xy + 0i);
	x.diff = c(x.diff, - x.diff);
	S = c(S, S);
	x0 = (S + x.diff)/2;
	y0 = S - x0 # TODO: include also x = y cases
	# modified roots
	x = S / 2;
	y = - (x0 - y0) * 0.5i;
	sol = cbind(x=as.vector(x), y=as.vector(y))
	sol
}

### Examples
R = 2
b = 3
#
sol = solve.htC.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3 - 3*x*y^2 + b*x # - R[1]
y^3 - 3*x^2*y + b*y # == 0


#########
### Ex 2:
R = -1
b = 4
#
sol = solve.htC.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^3 - 3*x*y^2 + b*x # - R[1]
y^3 - 3*x^2*y + b*y # == 0


######################

######################
### Asymmetric Ext ###

# x^3 + b*y + bc1*(x+y) = R
# y^3 + b*x + bc2*(x+y) = R

### Sum =>
x^3 + y^3 + (b + bc1 + bc2)*S - 2*R # = 0
S^3 - 3*x*y*S + (b + bc1 + bc2)*S - 2*R # = 0
# 3*x*y*S = S^3 + (b + bc1 + bc2)*S - 2*R;

### =>
# x^3 + b*y = R - bc1*S
### Prod =>
(x*y)^3 + b*(x^4 + y^4) + b^2*(x*y) - bc1*bc2*S^2 + R*(bc1+bc2)*S - R^2 # = 0
(x*y)^3 + b*(S^4 - 4*x*y*S + 2*(x*y)^2) + b^2*(x*y) - bc1*bc2*S^2 + R*(bc1+bc2)*S - R^2
b*S^4 - bc1*bc2*S^2 + R*(bc1+bc2)*S + (x*y)^3 + 2*b*(x*y)^2 - 4*b*x*y*S + b^2*(x*y) - R^2

### TODO



###########################
###########################

### Section C

###########################
### Binomial Expansions ###
###########################

# x^3 + 3*x*y^2 = R1
# y^3 + 3*x^2*y = R2

### Extensions:
# x^3 + 3*x*y^2 + b2*y^2 + b2*x*y + b1*y = R1
# y^3 + 3*x^2*y + b2*x^2 + b2*x*y + b1*x = R2

### Solution:

### Sum =>
# S^3 = R1 + R2

### Subst =>
# y = S - x
x^3 + 3*x*(S - x)^2 - R1 # = 0
4*x^3 - 6*S*x^2 + 3*S^2*x - R1 # = 0


### Solver:
solve.Bin.S2P3 = function(R, b=c(0,0), debug=TRUE) {
	S = roots(c(1, b[2], b[1], - (R[1] + R[2])));
	if(debug) print(S);
	x = sapply(S, function(S) roots(c(4, - 6*S, 3*S^2 - b[2]*S - b[1],
		b[2]*S^2 + b[1]*S - R[1])));
	y = rep(S, each=3) - x;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	return(sol);
}

### Examples:
R = c(-1,2)
#
sol = solve.Bin.S2P3(R);
x = sol[,1]; y = sol[,2];

### Test
x^3 + 3*x*y^2 # - R[1]
y^3 + 3*x^2*y # - R[2]

### degenerate Polynomial:
round0.p(poly.calc(x)) * 2^8

###############
### Extensions:

### Examples:
R = c(-1,2)
b = c(1,2)
#
sol = solve.Bin.S2P3(R, b=b);
x = sol[,1]; y = sol[,2];

### Test
x^3 + 3*x*y^2 + b[2]*y^2 + b[2]*x*y + b[1]*y # - R[1]
y^3 + 3*x^2*y + b[2]*x^2 + b[2]*x*y + b[1]*x # - R[2]

### Classic Polynomial:
round0.p(poly.calc(x)) * 2^6


############################

### Binomial Expansions:
### Derivatives

# x^3 + y^3 + 3*x^2*y + 3*x*y^2 - b2*y^2 + b1*y = R1
# b2*y^2 + b1*x = R2

### Solution:

### Sum =>
S^3 + b1*S - R1 - R2 # = 0

### Step 2:
b2*(S-x)^2 + b1*x - R2 # = 0
b2*x^2 + (b1 - 2*b2*S)*x + b2*S^2 - R2 # = 0


### Solver:
solve.Bin.S2P3 = function(R, b=c(1,1), debug=TRUE) {
	S = roots(c(1, 0, b[1], - (R[1] + R[2])));
	if(debug) print(S);
	x = sapply(S, function(S) roots(c(b[2], (b[1] - 2*b[2]*S), b[2]*S^2 - R[2])));
	y = rep(S, each=2) - x; # assumes b[2] != 0;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	return(sol);
}

### Examples:
R = c(-1/2, -1/2)
b = c(-1, 1)
#
sol = solve.Bin.S2P3(R, b=b);
x = sol[,1]; y = sol[,2];

### Test
x^3 + y^3 + 3*x^2*y + 3*x*y^2 - b[2]*y^2 + b[1]*y # - R[1]
b[2]*y^2 + b[1]*x # - R[2]


### P[6]
round0.p(poly.calc(x))

######################

### Derivatives

# x^3 + y^3 + b1*y = R1
# x*y*(x+y) + b1/3 * x = R2

### Solution:

### Sum =>
S^3 + b1*S - R1 - 3*R2 # = 0

### Step 2:
x*y*S + b1/3 * x - R2 # = 0
x*(S-x)*S + b1/3 * x - R2 # = 0
S*x^2 - x*(S^2 + b1/3) + R2 # = 0

### Solver:
solve.Bin.S2P3 = function(R, b=-1, debug=TRUE) {
	S = roots(c(1, 0, b[1], - (R[1] + 3*R[2])));
	if(debug) print(S);
	x = sapply(S, function(S) roots(c(S, - (S^2 + b[1]/3), R[2])));
	y = rep(S, each=2) - x; # assumes b[2] != 0;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	return(sol);
}

### Examples:
R = c(-1, -2)
b = c(-1)
#
sol = solve.Bin.S2P3(R, b=b);
x = sol[,1]; y = sol[,2];

### Test
x^3 + y^3 + b[1]*y # - R[1]
x*y*(x+y) + b[1]/3 * x # - R[2]

### P[6]
round0.p(poly.calc(x)) * 21

######################

### Derivatives 2: 

# x^3 + y^3 + 3*x*y^2 - 3*b2*y^2 + b1*y = R1
# x^2*y + b2*y^2 + b1/3 * x = R2

### Solution:

### Sum =>
S^3 + b1*S - R1 - 3*R2 # = 0

### Step 2:
x^2*y + b2*y^2 + b1/3 * x - R2 # = 0
x^2*(S-x) + b2*(S-x)^2 + b1/3 * x - R2 # = 0
x^3 - (S+b2)*x^2 - (b1/3 - 2*b2*S) * x - b2*S^2 + R2 # = 0

### Solver:
solve.Bin.S2P3 = function(R, b=c(-1, 0), debug=TRUE) {
	S = roots(c(1, 0, b[1], - (R[1] + 3*R[2])));
	if(debug) print(S);
	if(length(b) < 2) b = c(b, 0);
	x = sapply(S, function(S) roots(c(1, - (S+b[2]), - (b[1]/3 - 2*b[2]*S), R[2] - b[2]*S^2)));
	y = rep(S, each=3) - x; # assumes b[2] != 0;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	return(sol);
}

### Examples:
R = c(-1, -2)
b = c(-1)
#
sol = solve.Bin.S2P3(R, b=b);
x = sol[,1]; y = sol[,2];

### Test
x^3 + y^3 + 3*x*y^2 + b[1]*y # - R[1]
x^2*y + b[1]/3 * x # - R[2]

### P[9]
round0.p(poly.calc(x)) * 9


#########
### Ex 2:
R = c(-1, -2)
b = c(-1, 1)
#
sol = solve.Bin.S2P3(R, b=b);
x = sol[,1]; y = sol[,2];

### Test
x^3 + y^3 + 3*x*y^2 - 3*b[2]*y^2 + b[1]*y # - R[1]
x^2*y + b[2]*y^2 + b[1]/3 * x # - R[2]

### P[9]
round0.p(poly.calc(x)) * 9


#########
### Ex 3:
R = c(1, 0)
b = c(-1, -1)
#
sol = solve.Bin.S2P3(R, b=b);
x = sol[,1]; y = sol[,2];

### Test
x^3 + y^3 + 3*x*y^2 - 3*b[2]*y^2 + b[1]*y # - R[1]
x^2*y + b[2]*y^2 + b[1]/3 * x # - R[2]

### P[9] = simple P[3] o P[3]
round0.p(poly.calc(x)) * 9


###########################
###########################
###########################

### Section D

########################
### Asymmetric Shift ###
########################


### S2: d given
# (x + d)^n + y^n = R1
# x^n + (y + d)^n = R2


### S3: d unknown
# (x + d)^n + y^n = R1
# x^n + (y + d)^n = R2
# x^n + y^n = R3


###############
### Order 2 ###
###############


### S2: d given
# (x + d)^2 + y^2 = R1
# x^2 + (y + d)^2 = R2

### Sum Eq 2 + Eq 3:
(x + d)^2 + (y + d)^2 + x^2 + y^2 - R1 - R2 # = 0
2*(x^2 + y^2) + 2*d*(x + y) + 2*d^2 - R1 - R2 # = 0
2*S^2 + 2*d*S - 4*x*y + 2*d^2 - R1 - R2
# 4*x*y = 2*S^2 + 2*d*S + 2*d^2 - R1 - R2

### Diff =>
2*(x - y)*(x + y + d) - R1 + R2 # = 0
2*(x - y)*(S + d) - R1 + R2 # = 0

### Case: R1 == R2 && x != y
S + d # = 0
# S = -d

### Case: R1 != R2
# TODO


### Solver:
solve.AsymShift.S2P2 = function(R, d) {
	if(R[1] != R[2]) stop("Not yet impleemnted!")
	S = -d;
	xy = (2*S^2 + 2*d*S + 2*d^2 - R[1] - R[2]) / 4;
	#
	xy.d = sqrt(S^2 - 4*xy + 0i)
	x = (S + xy.d)/2;
	y = (S - xy.d)/2;
	cbind(x, y)
}

### Examples:
R = c(2,2)
d = 2
sol = solve.AsymShift.S2P2(R, d)
x = sol[1]; y = sol[2];

### Test
(x + d)^2 + y^2 # - R[1]
x^2 + (y + d)^2 # - R[2]


#################

### S3: d unknown
# (x + d)^2 + y^2 = R1
# x^2 + (y + d)^2 = R2
# x^2 + y^2 = R3

### Sum Eq 1 + Eq 2:
(x + d)^2 + (y + d)^2 + x^2 + y^2 - R1 - R2 # = 0
2*(x^2 + y^2) + 2*d*(x + y) + 2*d^2 - R1 - R2 # = 0
2*d*S + 2*d^2 + 2*R3 - R1 - R2
# 2*d*S = - 2*d^2 - 2*R3 + R1 + R2


### Case 1: R1 == R2 && x != y
### Diff: Eq 1 - Eq 2
2*(x - y)*(x + y + d) - R1 + R2 # = 0
2*(x - y)*(S + d) - R1 + R2 # = 0
S + d # = 0
# S = -d
# NO solution: if 2*R3 != R1 + R2;


### Case 2: R1 == R3
### Diff: Eq 1 - Eq 3
d*(2*x + d) - R1 + R3 # = 0
# d = 0 (requires R1 == R2) OR
# d = -2*x;


### Case 3: R1 != R2 != R3
# TODO


### Eq 3 =>
S^2 - 2*x*y - R3 # = 0
4*d^2*S^2 - 8*d^2*x*y - 4*d^2*R3


### Solver:
solve.AsymShift.S2P2 = function(R) {
	# TODO:
	d = 0
	#
	xy.d = sqrt(S^2 - 4*xy + 0i)
	x = (S + xy.d)/2;
	y = (S - xy.d)/2;
	cbind(x, y, d)
}

### Examples:
R = c(1,2,3)
#
sol = solve.AsymShift.S2P2(R)
x = sol[1]; y = sol[2];

### Test
(x + d)^2 + y^2 # - R[1]
x^2 + (y + d)^2 # - R[2]
x^2 + y^2 # - R[3]


########################
########################
########################

########################
### Composed Systems ###
########################

### x*y*(x+y)

# x^4*y^2 + b2*x^3*y^3 + b1*x*y^2 = R1
# x^2*y^4 + (2-b2)*x^3*y^3 + b1*x^2*y = R2

### Solution:

# 2 approaches:
# x1 = x^2*y; y1 = x*y^2;
# - solve for x1 & y1:
#   x1^2 + b2*x1*y1 + b1*y1 = R1;
# - decouples b2 from (2-b2);

### Approach 2:
### Sum =>
# V = x*y*(x+y)
V^2 + b1*V - R1 - R2 # = 0

### Step 2:
# - uglier;
x^4*y^2 + b2*x^3*y^3 + b1*x*y^2 - R1 # = 0
x^3*x*y*(S-x) + b2*x^3*y^3 + b1*x*y*(S-x) - R1 # = 0
x^3*x*y*S - x^5*y + b2*x^3*y^3 + b1*x*y*S - b1*x^2*y - R1 # = 0
V*x^3 - x^5*y + b2*x^2*y*x*y*(S - x) - b1*x^2*y + b1*V - R1 # = 0
V*x^3 - x^5*y + b2*V*x^2*y - b2*x^4*y^2 - b1*x^2*y + b1*V - R1 # = 0

### alternative Step 2:
### Diff =>
x^2*y^2*S*(x-y) + 2*(b2-1)*x^3*y^3 - b1*x*y*(x-y) - R1 + R2 # = 0
V*x*y*(2*x-S) + 2*(b2-1)*x^3*y^3 - b1*x*y*(2*x-S) - R1 + R2 # = 0
2*V*x^2*y + 2*(b2-1)*x^3*y^3 - 2*b1*x^2*y - V^2 + b1*V - R1 + R2
2*V*x^2*y + 2*(b2-1)*x^2*y*x*y*(S - x) - 2*b1*x^2*y - V^2 + b1*V - R1 + R2
2*V*x^2*y + 2*(b2-1)*V*x^2*y - 2*(b2-1)*x^4*y^2 - 2*b1*x^2*y - V^2 + b1*V - R1 + R2
2*(b2*V - b1)*x^2*y - 2*(b2-1)*x^4*y^2 - V^2 + b1*V - R1 + R2
### solve for x^2*y


### Solver:
solve.AsM.S2P42 = function(R, b, debug=TRUE) {
	V = roots(c(1, b[1], - (R[1] + R[2])));
	if(debug) print(V);
	x2y = sapply(V, function(V) roots(c(
		- 2*(b[2]-1), 2*(b[2]*V - b[1]), - V^2 + b[1]*V - R[1] + R[2])));
	len = length(x2y) / length(V)
	V = rep(V, each=len);
	# x, y
	ydx = (V - x2y) / x2y;
	x3 = x2y / ydx;
	x = rootn(x3, 3);
	m = unity(3, all=TRUE);
	x = sapply(x, function(x) x*m);
	y = x * rep(ydx, each=3);
	sol = cbind(x=as.vector(x), y=as.vector(y));
	return(sol);
}

### Examples:
R = c(-1, 2)
b = c(3, -3)
#
sol = solve.AsM.S2P42(R, b);
x = sol[,1]; y = sol[,2];

### Test
x^4*y^2 + b[2]*x^3*y^3 + b[1]*x*y^2 # - R1
x^2*y^4 + (2-b[2])*x^3*y^3 + b[1]*x^2*y # - R2


#########
### Ex 2:
R = c(-1, 2)
b = c(3, 2)
#
sol = solve.AsM.S2P42(R, b);
x = sol[,1]; y = sol[,2];

### Test
x^4*y^2 + b[2]*x^3*y^3 + b[1]*x*y^2 # - R1
x^2*y^4 + (2-b[2])*x^3*y^3 + b[1]*x^2*y # - R2

# degenerate P[12]
round0.p(poly.calc(x)) * 13


### Debug
R = c(-1, 2)
b = c(3, -3)
x =  0.5829666585 - 0.7755174525i;
y = -0.6425343366 + 0.3544215388i;
b1 = b[1]; b2 = b[2]; R1 = R[1]; R2 = R[2];
S = x+y; V = x*y*S;


#########################
#########################

#####################
### Special Cases ###
#####################

#######################
### Root Symmetries ###
#######################

### (x, y) & (-x, -y)

##############
### Order: 3+1

x^3*y + b*y^2 = R1
y^3*x + b*x^2 = R2

### Solution:

# - if (x, y) is a root, then:
#   (-x, -y) is also a root;

### Sum =>
# let: s = x^2 + y^2, p = x*y;
p*s + b*s - R1 - R2 # = 0
# s = (R1 + R2) / (p + b)

### Prod =>
p^4 + b*p*(s^2 - 2*p^2) + b^2*p^2 - R1*R2 # = 0
p^4*(p + b)^2 + b*p*((R1+R2)^2 - 2*p^2*(p + b)^2) + b^2*p^2*(p + b)^2 - R1*R2*(p + b)^2 # = 0
p^6 - 2*b^2*p^4 + (b^4 - R1*R2)*p^2 +
	+ b*p*(R1^2+R2^2) - R1*R2*b^2 # = 0


### Solver:
solve.AsymMinus.S2P31 = function(R, b, debug=TRUE) {
	if(R[1] == R[2]) print("Special case R1 == R2: NOT implemented!");
	coeff = c(1, 0, - 2*b[1]^2, 0, (b[1]^4 - R[1]*R[2]),
		b[1]*(R[1]^2+R[2]^2), - R[1]*R[2]*b^2);
	xy = roots(coeff);
	s2 = (R[1] + R[2]) / (xy + b[1]);
	if(debug) print(xy);
	# robust: (x*y)*x^2 + b*y^2 = R1;
	x2 = (R[1] - b[1]*s2) / (xy - b[1]);
	x = sqrt(x2); x = c(x, -x);
	y = rep(xy, 2) / x;
	sol = cbind(x = as.vector(x), y = as.vector(y));
	# sol = rbind(sol, - sol);
	sol = sol[order(sol[,1]),];
	return(sol);
}

### Examples:
R = c(-1, 3)
b = 2
sol = solve.AsymMinus.S2P31(R, b);
x = sol[,1]; y = sol[,2];

### Test
x^3*y + b[1]*y^2 # - R[1]
y^3*x + b[1]*x^2 # - R[2]


### degenerate P[12]
round0.p(poly.calc(x))


##############
### Extension:

x^3*y + b1*y^2 = R1
y^3*x + b2*x^2 = R2

### Symmetrization:
# let: y => k*y;
x^3*y + k*b1*y^2 - R1/k
y^3*x + b2/k^3*x^2 - R2/k^3
# k^4 = b2 / b1;


##################
##################

##################
###  Equality  ###
##################

### Equality Order 1
### x == y

### System Order 3

# x^3 + b11*(x^k1 - y^k1) = R
# y^3 + b21*(x^k2 - y^k2) = R

### Solution:

### Case 1: x == y
# - trivial;
x^3 - R # = 0

### Case 2: x != y
# for k1 = k2 = 1

### Diff =>
(x-y)*(S^2 - x*y + b11 - b21) # = 0
S^2 - x*y + b11 - b21 # = 0

### TODO:
# - properly solve for the remaining roots;


####################
####################

### Equality Order 3
### x^3 == y^3

### System Order 4
# x^4 + x^3*y - x*y^3 + b2*y^3 + b1*y = R
# y^4 + b2*x^3 + b1*y = R

### Solution:

### Case 1: x^3 == y^3
# - trivial;
y^4 + b2*y^3 + b1*y - R # = 0

### Case 2: x^3 != y^3

### Diff =>
(x^3 - y^3)*(x + y - b2) # = 0
# S = b2;

### Eq y:
y^4 - b2*y^3 + 3*b2^2*y^2 + b1*y - 3*b2^3*y - R + b2^4 # = 0


### Solver:
solve.S2Equal = function(R, b, sort=TRUE) {
	# [solved] removed false roots & duplicates;
	sol = solve.Case1.S2Equal(R, b, sort=sort)
	# set with 4 distinct (correct) roots;
	sol2 = solve.Case2.S2Equal(R, b);
	sol = rbind(sol2, sol);
	return(sol);
}
solve.Case1.S2Equal = function(R, b, sort=TRUE, add.duplicates=FALSE) {
	y = solve.Case1Eq.S2Equal(R, b, permute=FALSE)[,2];
	x = sapply(y, function(y) roots(c(b[2], 0, 0, y^4 + b[1]*y - R)));
	y = rep(y, each=3);
	sol = cbind(x=as.vector(x), y=as.vector(y));
	if(add.duplicates) {
		# probably duplicates: TODO: check!
		sol2 = solve.Case1Eq.S2Equal(R, b, permute=FALSE);
		sol = rbind(sol2, sol);
	}
	if(sort) {
		id = order.cm(sol, 5);
		sol = sol[id,];
	}
	return(sol);
}
solve.Case1Eq.S2Equal = function(R, b, permute=FALSE) {
	y = roots(c(1, b[2], 0, b[1], -R[1]));
	if(permute) {
		sol = roots.y.S2Equal(y, R=R, b=b, n=3);
	} else sol = cbind(x=y, y=y);
	return(sol);
}
solve.Case2.S2Equal = function(R, b) {
	y = roots(c(1, -b[2], 3*b[2]^2, b[1] - 3*b[2]^3, b[2]^4 - R[1]));
	x = b[2] - y;
	return(cbind(x=x, y=y));
}
### Case 2:
# - does NOT compute "all" roots;
# - does NOT work for Case 1;
roots.x.Case2.S2Equal = function(y, R, b, n=3) {
	# b1 = b[1]; b2 = b[2];
	x = b[2] - y;
	return(cbind(x=x, y=y));
}
### [old] Simple
# - NOT robust!
roots.y.S2Equal = function(x, R, b, n=3, sort=TRUE) {
	sol = cbind(x=x, y=x); # base-set;
	m = unity(3, all=FALSE);
	sol2 = cbind(x=rep(x, 2), y=c(x*m, x*m^2));
	sol2 = rbind(sol2, sol2[,2:1]);
	sol = rbind(sol, sol2);
	if(sort) {
		id = order.cm(sol, 5);
		sol = sol[id,];
	}
	return(sol);
}

### Examples:
R = 2
b = c(2,-1)
# [solved] removed false roots & duplicates;
sol = solve.S2Equal(R, b)
x = sol[,1]; y = sol[,2];


### Test:
b1 = b[1]; b2 = b[2];
x^4 + x^3*y - x*y^3 + b2*y^3 + b1*y # = R
y^4 + b2*x^3 + b1*y # = R


### Derivation:
### Case: x^3 != y^3
### x = b2 - y =>
y^4 - b2*y^3 + 3*b2^2*y^2 + b1*y - 3*b2^3*y - R + b2^4

### Eq y: Py-Case2[4] * PyEq[4]^3

b1 = b[1]; b2 = b[2];
# actually PyEq^2
coeff = c(1, 2*b2, b2^2, 2*b1, - 2*R + 2*b1*b2, - 2*R*b2, b1^2, - 2*R*b1, R^2);
y = roots(coeff)
x = sapply(y, function(y) roots(c(b[2], 0, 0, y^4 + b[1]*y - R)));
y = rep(y, each=3)

