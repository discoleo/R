########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Asymmetric S2:
### Base Types
###
### draft v.0.2f


### Asymmetric Polynomial Systems: 2 Variables
### Base Types

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


### draft v.0.2f:
# - generalized approach to Order 2 Asymmetric:
#   b1*x^2 + b2*y^2 = R2;
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


##########################

##########################
### Polynomial Systems ###
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
solve.AsymSimple.P3A2 = function(R, b, debug=TRUE) {
	bs = b[1] + b[2]; bp = b[1]*b[2];
	bsq = bs^2;
	coeff = c((bsq - 3*bp), 0, - 3*R[2]*bs, - 2*(bsq - 6*bp)*R[1],
		9*R[2]^2, - 6*bs*R[1]*R[2], bsq*R[1]^2)
	S = roots(coeff);
	if(debug) print(S);
	#
	R1 = R[1]; R2 = R[2];
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


##################
##################

##################
### Transforms ###

### Initial:
# x^2 + y^2 = R1
# x*y = R2

### Final:
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


###############

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

### Eq 2 / 3:
# 2*(R1^2 - R2^2)*x =
R1 * (-S^3 - b*S^2 + (3*R1 + 3*R2 + 2*b^2)*S - 2*b*(R1 + R2)) +
	- R2 * (S^3 + b*S^2 - (R1 + R2 + 2*b^2)*S + 2*b*(R1 + R2))
-(R1+R2)*S^3 - b*(R1+R2)*S^2 + (3*R1^2 + R2^2 + 4*R1*R2 + 2*b^2*R1+ 2*b^2*R2)*S - 2*b*(R1 + R2)^2
# 2*(R1^2 - R2^2)*y =
R1*(S^3 + b*S^2 - (R1 + R2 + 2*b^2)*S + 2*b*(R1 + R2)) +
	- R2*(-S^3 - b*S^2 + (3*R1 + 3*R2 + 2*b^2)*S - 2*b*(R1 + R2))
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
# - solve;
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
S^6 - 4*x*y*S^4 + 2*(x*y)^2*S^2 + 2*b*x*y*S^2 - R1*x*S^2 - R2*y*S^2 # = 0
-3*S^6 - 6*b*S^4 + 6*(2*R1 + 2*R2)*S^3 + 2*(S^3 + b*S - R1 - R2)^2 +
	+ 6*b^2*S^2 - 6*b*(R1 + R2)*S - 9*R1*x*S^2 - 9*R2*y*S^2 # = 0
-S^6 - 2*b*S^4 + (8*R1 + 8*R2)*S^3 + 8*b^2*S^2 - 10*b*(R1 + R2)*S +
	- 9*R1*x*S^2 - 9*R2*y*S^2 + 2*(R1+R2)^2 # = 0
# (9*R1*x + 9*R2*y)*S^2 =
	-S^6 - 2*b*S^4 + (8*R1 + 8*R2)*S^3 + 8*b^2*S^2 - 10*b*(R1 + R2)*S + 2*(R1+R2)^2


### Diff =>
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
#######################

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


###########################
###########################

########################
### Asymmetric Shift ###

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

