########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogeneous Symmetric S3:
### Mixed Type: Composite
###
### draft v.0.1d


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
# or
# x*y*z = R2


######################
######################

###############
### History ###
###############


### draft v.0.1d:
# - variant: x^2*y + y^2*z + z^2*x = R;
#   [started work]
### draft v.0.1c:
# - variant: x*y*z = R;
### draft v.0.1b:
# - A2 extension:
#   x*y + x*z + y*z + be*S = R;
# - P[6] example:
#   11 + 2*x + 5*x^2 - 2*x^3 + x^6 = 0;
### draft v.0.1a - v.0.1a-eq:
# - simple linear system:
#   x^2 + b*y = Ru;
# - classic P[6] example:
#   1 - x^2 + 2*x^3 - 2*x^5 + x^6 = 0;
# - variant: x^2 + y^2 + z^2 = R; [v.0.1a-var]
# - explicit handling of cases: x = y = z; [v.0.1a-eq]


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

### Extensions:
### A1: has NO impact;
x^2 + b*y + be*S = Ru
### A2:
x*y + x*z + y*z + be*S = R

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
solve.CompLin.S3P2 = function(R, b, be=0, debug=TRUE) {
	coeff = c(1, - 2*b[1] - be, R[1] + 3*b[1]^2);
	S  = roots(coeff);
	if(debug) print(S);
	len = length(S);
	E2 = rep(R[1], len) - be[1]*S;
	Ru = (S^2 - 2*E2 + b[1]*S) / 3;
	E3  = (S^3 + 2*b[1]*S^2 - 7*Ru*S - b[1]^2*S + 3*b[1]*Ru) / 6;
	x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	#
	Ru = rep(Ru, each=3); S = rep(S, each=3);
	y = (Ru - x^2) / b[1];
	z = S - x - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z));
	sol = sol[order(abs(sol[,1])), ];
	### Case: x == y == z
	x = roots(c(1, 3*be[1], -3*R[1])) / 3;
	solEq = cbind(x=x, y=x, z=x);
	sol = rbind(sol, solEq);
	return(sol);
}

### Examples:

R = -1
b = 1
sol = solve.CompLin.S3P2(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
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


#########
### Ex 2:
# A2 Extension
R = -1
b = 1
be = -2
sol = solve.CompLin.S3P2(R, b, be)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
S = (x+y+z);
cbind(
x^2 + b*y,
y^2 + b*z,
z^2 + b*x)
x*y + x*z + y*z + be[1]*S # - R

###
round0.p(poly.calc(x[1:6]))
x = x[1:6]
err = 11 + 2*x + 5*x^2 - 2*x^3 + x^6
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
	sol = sol[order(abs(sol[,1])), ];
	### Case: x == y == z
	x = roots(c(1, 0, -3*R[1])) / 3;
	solEq = cbind(x=x, y=x, z=x);
	sol = rbind(sol, solEq);
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


###############

###############
### Variant ###

x^2 + b*y = Ru # where Ru = unknown
y^2 + b*z = Ru
z^2 + b*x = Ru
x*y*z = R

### Solution

### Sum =>
x^2 + y^2 + z^2 + b1*(x+y+z) - 3*Ru # = 0
# 3*Ru = S^2 - 2*E2 + b1*S;

### Sum(x[i]*...) =>
6*E3 - S^3 - 2*b1*S^2 + 7*Ru*S + b1^2*S - 3*b1*Ru # = 0
# (7*S - 3*b1)*Ru = (S^3 + 2*b1*S^2 - b1^2*S - 6*E3);

### Eq:
(S^2 + 3*b1*S - 9*Ru)*(S^2 - b1*S - Ru + 2*b1^2) # = 0
(7*S - 3*b1)*(S^2 - b1*S + 2*b1^2) - (S^3 + 2*b1*S^2 - b1^2*S - 6*E3)
S^3 - 2*b1*S^2 + 3*b1^2*S - b1^3 + E3

### Solver:
solve.CompLin.S3P2 = function(R, b, be=0, debug=TRUE) {
	coeff = c(1, - 2*b[1], 3*b[1]^2 - be[1], R[1] - b[1]^3);
	S  = roots(coeff);
	if(debug) print(S);
	len = length(S);
	E3 = R[1] - be[1]*S;
	Ru = (S^3 + 2*b[1]*S^2 - b[1]^2*S - 6*E3) / (7*S - 3*b[1]);
	E2 = (S^2 + b[1]*S - 3*Ru) / 2;
	x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	#
	Ru = rep(Ru, each=3); S = rep(S, each=3);
	y = (Ru - x^2) / b[1];
	z = S - x - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z));
	sol = sol[order(
		round(abs(sol[,1]), 10),
		round(abs(Re(sol[,1])), 10)), ];
	### Case: x == y == z
	x = roots(c(1, 0, 3*be[1], -R[1]));
	solEq = cbind(x=x, y=x, z=x);
	sol = rbind(sol, solEq);
	return(sol);
}

### Examples:

R = -1
b = 1
sol = solve.CompLin.S3P2(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
cbind(
x^2 + b*y,
y^2 + b*z,
z^2 + b*x)
x*y*z # - R

###
round0.p(poly.calc(x[1:9]))
x = x[1:9]
err = 1 - 3*x + 2*x^2 - x^3 + x^4 - x^5 + x^6 - 2*x^8 + x^9
round0(err)


#########
### Ex 2:
R = -1
b = 2
be = -1
sol = solve.CompLin.S3P2(R, b, be)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
cbind(
x^2 + b*y,
y^2 + b*z,
z^2 + b*x)
x*y*z + be[1]*(x+y+z) # - R

###
round0.p(poly.calc(x[1:9]))
x = x[1:9]
err = 1 - 77*x + 119*x^2 - 200*x^3 + 94*x^4 - 53*x^5 + 7*x^6 + 3*x^7 - 4*x^8 + x^9
round0(err)


###############

###############
### Variant ###

x^2 + b*y = Ru # where Ru = unknown
y^2 + b*z = Ru
z^2 + b*x = Ru
x^2*y + y^2*z + z^2*x = R

### Solution

### Sum =>
x^2 + y^2 + z^2 + b1*(x+y+z) - 3*Ru # = 0
# 3*Ru = S^2 - 2*E2 + b1*S;
3*(S^2 - b1*S + 2*b1^2) - (S^2 - 2*E2 + b1*S) # = 0
S^2 - 2*b1*S + 3*b1^2 + E2 # = 0
# E2 = - (S^2 - 2*b1*S + 3*b1^2);

### Eq 2:
(S^2 + 3*b1*S - 9*Ru)*(S^2 - b1*S - Ru + 2*b1^2) # = 0
# Ru = S^2 - b1*S + 2*b1^2

### Sum(x[i]*...) =>
6*E3 - S^3 - 2*b1*S^2 + 7*Ru*S + b1^2*S - 3*b1*Ru # = 0
# (7*S - 3*b1)*Ru = (S^3 + 2*b1*S^2 - b1^2*S - 6*E3);
(7*S - 3*b1)*(S^2 - 2*E2 + b1*S) - 3*(S^3 + 2*b1*S^2 - b1^2*S - 6*E3) # = 0
4*S^3 - 14*E2*S - 2*b1*S^2 + 6*b1*E2 + 18*E3 # = 0
4*S^3 + 14*(S^2 - 2*b1*S + 3*b1^2)*S - 2*b1*S^2 +
	- 6*b1*(S^2 - 2*b1*S + 3*b1^2) + 18*E3 # = 0
S^3 - 2*b1*S^2 + 3*b1^2*S - b1^3 + E3 # = 0
# E3 = - (S^3 - 2*b1*S^2 + 3*b1^2*S - b1^3);


### R*(x^2*y + y^2*z + z^2*x) + R*(x^2*z + y^2*x + z^2*y) =
R^2 + (x^2*y + y^2*z + z^2*x)*(x^2*z + y^2*x + z^2*y) # =
((x*y)^3 + (x*z)^3 + (y*z)^3) + E3*(x^3 + y^3 + z^3) + 3*E3^2 + R^2
E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2 + R^2
# =>
E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2 + R^2 - R*(E2*S - 3*E3) # = 0
E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2 + 3*R*E3 - R*E2*S + R^2 # = 0
### TODO: seems to be S^6;

