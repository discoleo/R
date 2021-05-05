########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogeneous Symmetric S3:
### Mixed Type: Composite
###
### draft v.0.1a-var


### Heterogeneous Symmetric
### Polynomial Systems: 3 Variables
### Mixed/Composite: Hetero + Symmetric

### Example:
# x^p + b*y = Ru # where Ru = unknown
# y^p + b*z = Ru
# z^p + b*x = Ru
# x*y + x*z + y*z = R2

### Variants:
# x^n + y^n + z^n = R2


######################
######################

###############
### History ###
###############


### draft v.0.1a - v.0.1a-var:
# - simple linear system:
#   x^2 + b*y = Ru;
# - variant: x^2 + y^2 + z^2 = R;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R


######################
######################

###################
### Liniar Type ###
###################

###############
### Order 2 ###
###############

x^2 + b*y = Ru # where Ru = unknown
y^2 + b*z = Ru
z^2 + b*x = Ru
x*y + x*z + y*z = R

### Solution

### Sum =>
x^2 + y^2 + z^2 + b1*(x+y+z) - 3*Ru # = 0
S^2 - 2*E2 + b1*S - 3*Ru # = 0
# 3*Ru = S^2 - 2*E2 + b1*S;

### Sum(x[i]*...) =>
6*E3 - S^3 - 2*b1*S^2 + 7*Ru*S + b1^2*S - 3*b1*Ru # = 0

### Sum(x[i+1]^2*...) =>
### Eq:
(S^2 + 3*b1*S - 9*Ru)*(S^2 - b1*S - Ru + 2*b1^2) # = 0
(S^2 - 3*E2)*(S^2 - 2*b1*S + E2 + 3*b1^2) # = 0


### Solver:
solve.CompLin.S3P2 = function(R, b, debug=TRUE) {
	coeff = c(1, - 2*b[1], R[1] + 3*b[1]^2);
	S  = roots(coeff);
	S2 = roots(c(1, 0, -3*R[1]));
	S = c(S, S2);
	if(debug) print(S);
	len = length(S);
	E2 = rep(R[1], len);
	Ru = (S^2 - 2*E2 + b[1]*S) / 3;
	E3 = (S^3 + 2*b[1]*S^2 - 7*Ru*S - b[1]^2*S + 3*b[1]*Ru) / 6;
	x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	#
	Ru = rep(Ru, each=3); S = rep(S, each=3);
	y = (Ru - x^2) / b[1];
	z = S - x - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z));
	return(sol);
}

### Examples:

R = -1
b = 1
sol = solve.CompLin.S3P2(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

cbind(
x^2 + b*y,
y^2 + b*z,
z^2 + b*x)
x*y + x*z + y*z # - R

###
round0.p(poly.calc(x[1:6]))
x = x[1:6]
err = 1 - x^2 + 2*x^3 - 2*x^5 + x^6
round0(err)


###############

###############
### Variant ###

x^2 + b*y = Ru # where Ru = unknown
y^2 + b*z = Ru
z^2 + b*x = Ru
x^2 + y^2 + z^2 = R

### Solution

### Sum =>
x^2 + y^2 + z^2 + b1*(x+y+z) - 3*Ru # = 0
# 3*Ru = b1*S + R;
# 2*E2 = S^2 - R;

### Eq:
(S^2 + 3*b1*S - 9*Ru)*(S^2 - b1*S - Ru + 2*b1^2) # = 0
(S^2 - 3*R)*(3*S^2 - 4*b1*S - R + 6*b1^2) # = 0


### Solver:
solve.CompLin.S3P2 = function(R, b, debug=TRUE) {
	coeff = c(3, - 4*b[1], -R[1] + 6*b[1]^2);
	S  = roots(coeff);
	S2 = roots(c(1, 0, -3*R[1]));
	S = c(S, S2);
	if(debug) print(S);
	len = length(S);
	E2 = (S^2 - R) / 2;
	Ru = (S^2 + b[1]*S - 2*E2) / 3;
	E3 = (S^3 + 2*b[1]*S^2 - 7*Ru*S - b[1]^2*S + 3*b[1]*Ru) / 6;
	x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	#
	Ru = rep(Ru, each=3); S = rep(S, each=3);
	y = (Ru - x^2) / b[1];
	z = S - x - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z));
	return(sol);
}

### Examples:

R = -1
b = 1
sol = solve.CompLin.S3P2(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

cbind(
x^2 + b*y,
y^2 + b*z,
z^2 + b*x)
x^2 + y^2 + z^2 # - R

###
round0.p(poly.calc(x[1:6]) * 27)
x = x[1:6]
err = 7 - 20*x + 32*x^2 - 38*x^3 + 51*x^4 - 36*x^5 + 27*x^6
round0(err)


#########
### Ex 2:
R = 3
b = -3
sol = solve.CompLin.S3P2(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

cbind(
x^2 + b*y,
y^2 + b*z,
z^2 + b*x)
x^2 + y^2 + z^2 # - R

###
round0.p(poly.calc(x[1:6]))
x = x[1:6]
err = 173 + 28*x - 8*x^2 + 6*x^3 + 5*x^4 + 4*x^5 + x^6
round0(err)
