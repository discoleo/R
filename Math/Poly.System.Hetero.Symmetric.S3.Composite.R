########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogeneous Symmetric S3:
### Mixed Type: Composite
###
### draft v.0.1e


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


### draft v.0.1e:
# - P[6] classic poly (for the case: E2 = R);
### draft v.0.1d - v.0.1d-true:
# - variant: x^2*y + y^2*z + z^2*x = R;
#   [solved & factorized: true solutions]
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
solve.CompLin.S3P2 = function(R, b, be=0, all.roots=TRUE, debug=TRUE) {
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
	if(all.roots) {
		x = roots(c(1, 3*be[1], -3*R[1])) / 3;
		solEq = cbind(x=x, y=x, z=x);
		sol = rbind(sol, solEq);
	}
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


#########
### Ex 3:
b = 1
be = -3
R = - (3*b + be)*(b - be) / 3;
sol = solve.CompLin.S3P2(R, b, be)
x = sol[,1]; y = sol[,2]; z = sol[,3];

###
round0.p(poly.calc(x[1:6]))
x = x[1:6]
err = 67 + 30*x + 17*x^2 - 5*x^3 + x^5 + x^6
round0(err)


### P[6]
b1 = b[1]; be = be;
x^6 - (2*b1 + be)*x^5 +
	+ (3*R + 3*b1^2 - 2*b1*be - be^2)*x^4 + # (3*b1 + be)*(b1 - be)
	- (2*(2*b1 + be)*R + 2*b1^3 - 4*b1^2*be - 4*b1*be^2 - be^3)*x^3 +
	+ (3*R^2 + 6*b1^2*R + 2*b1^4 - 4*b1*be*R - 5*b1^3*be - be^2*R)*x^2 +
	- (2*b1*R^2 + 2*b1^3*R + be*R^2 - 4*b1^2*be*R - 2*b1^4*be - 2*b1*be^2*R +
		+ 5*b1^3*be^2 + 3*b1^2*be^3)*x +
	+ R^3 + 144 - 108*b1 + 3*b1^2*R^2 - 88*b1^2 + 39*b1^3 + 2*b1^4*R + 17*b1^4 - 3*b1^5 +
		- 2*b1*be*R^2 - 3*b1^3*be*R + 2*b1^5*be + 3*b1^2*be^2*R + 5*b1^4*be^2 - b1^3*be^3


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
S^6 - 4*b1*S^5 + 6*b1^2*S^4 - 3*b1^3*S^3 - 2*R*S^3 - 12*b1^4*S^2 +
	4*b1*R*S^2 + 18*b1^5*S - 6*b1^2*R*S - 18*b1^6 + 3*b1^3*R + R^2
### P[6] = (S^3 + 3*b1^3 - R) * P[3];
(S^3 - 4*b1*S^2 + 6*b1^2*S - 6*b1^3 - R) # true roots


### Solver:
solve.CompLin.S3P2 = function(R, b, be=0, sort=TRUE, debug=TRUE) {
	coeff = c(1, - 4*b[1], 6*b[1]^2, - 6*b[1]^3 - R[1]);
	S = roots(coeff);
	if(debug) print(S);
	len = length(S);
	E3 = - (S^3 - 2*b[1]*S^2 + 3*b[1]^2*S - b[1]^3);
	Ru = S^2 - b[1]*S + 2*b[1]^2;
	E2 = - (S^2 - 2*b[1]*S + 3*b[1]^2);
	x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	#
	Ru = rep(Ru, each=3); S = rep(S, each=3);
	y = (Ru - x^2) / b[1];
	z = S - x - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z));
	if(sort) sol = sol[order(
		round(abs(sol[,1]), 10),
		round(abs(Re(sol[,1])), 10)), ];
	### Case: x == y == z
	x = roots(c(3, 0, 3*be[1], -R[1]));
	solEq = cbind(x=x, y=x, z=x);
	sol = rbind(sol, solEq);
	return(sol);
}

### Examples:

R = 2
b = -1
sol = solve.CompLin.S3P2(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
cbind(
x^2 + b*y,
y^2 + b*z,
z^2 + b*x)
x^2*y + y^2*z + z^2*x # - R

###
round0.p(poly.calc(x[1:9]))
x = x[1:9]
err = -5 - 3*x + 7*x^2 + 14*x^3 + 4*x^4 - 15*x^5 - 13*x^6 + x^7 + 4*x^8 + x^9
round0(err)


### Ex 2:
R = 5
b = -1
sol = solve.CompLin.S3P2(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
cbind(
x^2 + b*y,
y^2 + b*z,
z^2 + b*x)
x^2*y + y^2*z + z^2*x # - R

###
round0.p(poly.calc(x[1:9]))
x = x[1:9]
err = 19 - 90*x - 56*x^2 + 44*x^3 + 43*x^4 - 9*x^5 - 16*x^6 + x^7 + 4*x^8 + x^9
round0(err)


### Debug
R = 2; b = -1; b1 = b[1];
x =  0.7952618850 + 0.1596169599i;
y = -0.3930361082 - 0.7461254313i;
z = -1.4022257769 - 0.4134915286i;
S = x+y+z; E2 = x*(y+z)+y*z; E3 = x*y*z;
Ru = x^2 + b[1]*y;


### Derivation
p3 = list(
	S = 3:0,
	b1 = 0:3,
	coeff = c(-1, 2, -3, 1)
)
# - (S^2 - 2*b1*S + 3*b1^2)
p2 = list(
	S = 2:0,
	b1 = 0:2,
	coeff = c(-1, 2, -3)
)
pS = list(S=1, coeff=1)
pR = list(R=1, coeff=1)
#
pRez = add.pm(mult.pm(p3, pow.pm(pS, 3)), pow.pm(p2, 3))
pRez = add.pm(pRez, mult.all.pm(list(p3, -6, p2, pS)))
pRez = add.pm(pRez, mult.sc.pm(pow.pm(p3, 2), 9))
pRez = add.pm(pRez, mult.all.pm(list(p3, 3, pR)))
pRez = add.pm(pRez, mult.all.pm(list(p2, -1, pR, pS)))
pRez

print.poly(pR)

