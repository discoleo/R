########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Hetero-Symmetric S4: Mixed
### Basic Types
###
### draft v.0.2a


##############
### Theory ###
##############

### System:
# x1^n + x2^n + x3^n + x4^n = R1
# x1*x2 + x2*x3 + x3*x4 + x4*x1 = R2
# x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 = R3
# x1*x2*x3*x4 = R4

### Symmetries:
# - if (x1, x2, x3, x4) is a solution,
#   so is every cyclic permutation, e.g. (x2, x3, x4, x1);
# Note:
# - Eq. 2 is NOT fully symmetric: it is not the full E2!

### Solution:

### Step 1:
# - solve a P[2*n] in S = x1 + x2 + x3 + x4;

### Step 2:
# - solve for x1: either using a P[2] o {P[2], P[2]} approach,
#   or using a P[4] (the lazy approach);
# - solve the remaining variables;


### Basic Derivations:
# - moved to file:
#   Poly.System.S4.HtMixed.Basic.Derivation.R;


### Q:
# Naming of Functions:
# Use: E2 vs E2a, E22 vs E22a, etc ?


####################
####################

### Helper Functions

source("Poly.System.S4.HtMixed.Basic.Helper.R")

### Elementary Polynomials:
# source("Polynomials.Helper.EP.R")

# - Helper functions & Base-Solvers:
#   moved to file:
#   Poly.System.S4.HtMixed.Basic.Helper.R;


###############
###############

### Equations:

### Other:
E3^2 - E121a*(E2 - E2a) + E4*S^2 - 4*E2*E4 # = 0
# Derivation:
E222 - E121a*(E2 - E2a) + E4*(S^2 - 2*E2) # = 0
# Epoly.gen(2, v=4, e=3)
E222 - E3^2 + 2*E2*E4 # = 0

#####################
#####################

###############
### Order 1 ###
###############

x1 + x2 + x3 + x4 - R1 # = 0
x1*x2 + x2*x3 + x3*x4 + x4*x1 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0


### Solution:

### P[2] for E2:
S = R[1]; E2a = R[2]; E3 = R[3]; E4 = R[4];
E2a*E2^2 - (S*E3 + 2*E2a^2)*E2 +
	+ S^2*E4 + S*E2a*E3 + E2a^3 - 4*E2a*E4 + E3^2 # = 0


### Solver:
coeff.S4Ht.E2P1 = function(R) {
	# - coefficients for E2, based on E2a;
	S = R[1]; E2a = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(E2a, - (S*E3 + 2*E2a^2),
		S^2*E4 + S*E2a*E3 + E2a^3 - 4*E2a*E4 + E3^2);
	return(coeff);
}
solve.S4HtM.E2P1 = function(R, sort=TRUE, debug=TRUE) {
	S = R[1]; E2a = R[2]; E3 = R[3]; E4 = R[4];
	coeff = coeff.S4Ht.E2P1(R);
	E2 = roots(coeff);
	# Note:
	# - classic approach needs only P[2] o P[2];
	x1 = sapply(E2, function(e2) roots(c(1, -S, e2, -E3, E4)));
	x1 = as.vector(x1); E2 = rep(E2, each=4);
	# robust:
	div = (x1^2*(2*x1-S) + (E2-E2a)*x1);
	x3 = (E4 - x1^4 + x1^3*S - E2a*x1^2) / div;
	# x2, x4:
	xs = S - x1 - x3; x24 = E4 / (x1*x3);
	xd = sqrt(xs^2 - 4*x24);
	x2 = (xs + xd)/2; x4 = (xs - xd)/2;
	sol = cbind(x1=x1, x2=as.vector(x2), x3=as.vector(x3), x4=as.vector(x4))
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}

### Examples:

R = c(-1,-2,3,1)
sol = solve.S4HtM.E2P1(R)

test.S4HtMixed(sol, n=1)


### Ex 2:
R = c(-1,-5,3,2)
sol = solve.S4HtM.E2P1(R)

test.S4HtMixed(sol, n=1)

### Classic Poly:
round0(poly.calc(sol[,1]) * abs(R[2]))


### Derivation:
# - moved to file:
#   Poly.System.S4.HtMixed.Basic.Derivation.R;


########################
########################

###############
### Order 2 ###
###############

x1^2 + x2^2 + x3^2 + x4^2 - R1 # = 0
x1*x2 + x2*x3 + x3*x4 + x4*x1 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

### Abbrev:
E2a = x1*x2 + x2*x3 + x3*x4 + x4*x1;
E2b = x1*x3 + x2*x4;
# x1*x3*E2b - (x1*x3)^2 - R4 = 0;

### P[4]
R2*S^4 - 2*E3*S^3 + 2*(2*E4 - 2*R2^2 - R1*R2)*S^2 +
	+ 2*E3*(R1+2*R2)*S - 16*R2*E4 + (4*R2^3 + 4*R1*R2^2 + R1^2*R2) + 4*E3^2 # = 0


### Solver:
coeff.S4Ht.E2P2 = function(R) {
	# - technically are coeffs for S^4, not for E2;
	# - but are computed using E2a:
	#   R2 = E2a;
	R1 = R[1]; R2 = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(R2, - 2*E3, 2*(2*E4 - 2*R2^2 - R1*R2),
		2*E3*(R1+2*R2), - 16*R2*E4 + (4*R2^3 + 4*R1*R2^2 + R1^2*R2) + 4*E3^2);
	return(coeff);
}
solve.S4HtM.E2P2 = function(R, sort=TRUE, debug=TRUE) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	coeff = coeff.S4Ht.E2P2(R);
	S = roots(coeff);
	if(debug) print(S);
	#
	E2 = (S^2 - R1) / 2;
	len = length(S);
	x1 = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -R3, R4)));
	S = rep(S, each=4); E2 = rep(E2, each=4);
	x1 = as.vector(x1);
	# robust:
	div = (x1^2*(2*x1-S) + (E2-R2)*x1);
	x3 = (R4 - x1^4 + x1^3*S - R2*x1^2) / div;
	# x1 = rep(x1, each=2); S = rep(S, each = 2);
	xs = S - x1 - x3; x24 = R4 / (x1*x3);
	xd = sqrt(xs^2 - 4*x24);
	x2 = (xs + xd)/2; x4 = (xs - xd)/2;
	sol = cbind(x1=x1, x2=as.vector(x2), x3=as.vector(x3), x4=as.vector(x4))
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}

### Examples:

### Ex 1:
R = c(1,-1,2,3)
sol = solve.S4HtM.E2P2(R)

test.S4HtMixed(sol)


### Ex 2:
R = c(0,-2,2,1)
sol = solve.S4HtM.E2P2(R)

test.S4HtMixed(sol)


### Ex 3:
R = c(-1,-2,0,3)
sol = solve.S4HtM.E2P2(R)

test.S4HtMixed(sol)


### Ex 4:
# R2 = 0;
R = c(-1,0,5,3)
sol = solve.S4HtM.E2P2(R)

test.S4HtMixed(sol)


########################
########################

###############
### Order 3 ###
###############

n = 3
x1^n + x2^n + x3^n + x4^n - R1 # = 0
x1*x2 + x2*x3 + x3*x4 + x4*x1 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

### Eq 1 =>
S^3 - 3*E2*S + 3*E3 - R1 # = 0

### P[6]
R2*S^6 - 3*E3*S^5 + (9*E4 - 6*R2^2)*S^4 + R2*(15*E3 - 2*R1)*S^3 +
	- (36*E4*R2 - 3*R1*E3 - 9*R2^3)*S^2 - 6*R2^2*(3*E3 - R1)*S + 9*E3^2*R2 - 6*R1*E3*R2 + R1^2*R2


### Solver:
coeff.S4Ht.E2P3 = function(R) {
	# R2 = E2a;
	R1 = R[1]; R2 = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(R2, - 3*E3, 9*E4 - 6*R2^2, 15*E3*R2 - 2*R1*R2,
		- 36*E4*R2 + 3*R1*E3 + 9*R2^3,
		- 18*E3*R2^2 + 6*R1*R2^2,
		9*E3^2*R2 - 6*R1*E3*R2 + R1^2*R2);
	return(coeff);
}
solve.S4HtM.E2P3 = function(R, sort=TRUE, debug=TRUE) {
	coeff = coeff.S4Ht.E2P3(R);
	S = roots(coeff);
	if(debug) print(S);
	hasZero = (round0(S) == 0);
	if(any(hasZero)) S = S[ ! hasZero];
	# TODO: S = 0;
	len = length(S);
	#
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	E2 = (S^3 + 3*R3 - R1) / (3*S);
	#
	x1 = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -R3, R4)));
	x1 = as.vector(x1);
	S  = rep(S, each=4); E2 = rep(E2, each=4);
	# robust:
	div = (x1^2*(2*x1-S) + (E2-R2)*x1);
	x3 = (R4 - x1^4 + x1^3*S - R2*x1^2) / div;
	xs = S - x1 - x3; x24 = R4 / (x1*x3);
	xd = sqrt(xs^2 - 4*x24);
	x2 = (xs + xd)/2; x4 = (xs - xd)/2;
	sol = cbind(x1=x1, x2=as.vector(x2), x3=as.vector(x3), x4=as.vector(x4))
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}

### Examples:

R = c(-2,-1,2,3)
sol = solve.S4HtM.E2P3(R);

test.S4HtMixed(sol, n=3)


### Ex 2:
# b0 = 0
R = c(-2,0,3,1)
sol = solve.S4HtM.E2P3(R);

test.S4HtMixed(sol, n=3)


### Derivation:

p1 = toPoly.pm("S^3 - 3*E2*S + 3*E3 - R1")
p2 = toPoly.pm("R2*E2^2 - (S*E3 + 2*R2^2)*E2 +
	+ S^2*E4 + S*R2*E3 + R2^3 - 4*R2*E4 + E3^2")
#
pR = solve.pm(p1, p2, "E2")
str(pR)
pR$Rez = sort.pm(pR$Rez, "S", xn2=c("E4", "E3", "R1"))
print.pm(pR$Rez, lead="S")
print.coeff(pR$Rez, "S")


########################
########################

###############
### Order 4 ###
###############

n = 4
x1^n + x2^n + x3^n + x4^n - R1 # = 0
x1*x2 + x2*x3 + x3*x4 + x4*x1 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

### Eq 1 =>
S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 - 4*E4 - R1 # = 0

### P[8]
R2^2*S^8 - 4*E3*R2*S^7 + (12*E4*R2 + 2*E3^2 - 8*R2^3)*S^6 - E3*(8*E4 - 28*R2^2)*S^5 +
	+ 4*E4^2*S^4 - 72*E4*R2^2*S^4 - 12*E3^2*R2*S^4 - 2*R2^2*R1*S^4 + 20*R2^4*S^4 +
	+ 40*E3*E4*R2*S^3 + 4*E3*R2*R1*S^3 - 56*E3*R2^3*S^3 +
	- 16*E4^2*R2*S^2 + 4*E4*R2*R1*S^2 + 104*E4*R2^3*S^2 - 2*E3^2*R1*S^2 + 20*E3^2*R2^2*S^2 + 8*R2^3*R1*S^2 - 16*R2^5*S^2 +
	- 16*E3*E4*R2^2*S - 8*E3^3*R2*S - 12*E3*R2^2*R1*S + 24*E3*R2^4*S +
	+ 16*E4^2*R2^2 - 16*E3^2*E4*R2 - 8*E4*R2^2*R1 - 48*E4*R2^4 + 4*E3^4 + 4*E3^2*R2*R1 + 8*E3^2*R2^3 + R2^2*R1^2 - 4*R2^4*R1 + 4*R2^6


### Solver:
coeff.S4Ht.E2P4 = function(R) {
	# R2 = E2a;
	R1 = R[1]; R2 = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(R2^2, - 4*E3*R2, 12*E4*R2 + 2*E3^2 - 8*R2^3,
		- 8*E3*E4 + 28*E3*R2^2,
		4*E4^2 - 72*E4*R2^2 - 12*E3^2*R2 - 2*R2^2*R1 + 20*R2^4,
		40*E3*E4*R2 + 4*E3*R2*R1 - 56*E3*R2^3,
		- 16*E4^2*R2 + 4*E4*R2*R1 + 104*E4*R2^3 - 2*E3^2*R1 + 20*E3^2*R2^2 + 8*R2^3*R1 - 16*R2^5,
		- 16*E3*E4*R2^2 - 8*E3^3*R2 - 12*E3*R2^2*R1 + 24*E3*R2^4,
		16*E4^2*R2^2 - 16*E3^2*E4*R2 - 8*E4*R2^2*R1 - 48*E4*R2^4 + 4*E3^4 + 4*E3^2*R2*R1 + 8*E3^2*R2^3 +
			+ R2^2*R1^2 - 4*R2^4*R1 + 4*R2^6);
	return(coeff);
}
solve.S4HtM.E2P4 = function(R, sort=TRUE, debug=TRUE) {
	coeff = coeff.S4Ht.E2P4(R);
	S = roots(coeff);
	if(debug) print(S);
	hasZero = (round0(S) == 0);
	if(any(hasZero)) S = S[ ! hasZero];
	# TODO: S = 0;
	len = length(S);
	#
	R1 = R[1]; R2 = R[2]; E3 = R[3]; E4 = R[4]; R3 = R[3]; R4 = R[4];
	E2 = (2*E3^2 - 4*E4*R2 + 2*R2^3 - 2*E3*R2*S + 2*E4*S^2 - R2*S^4 + R2*R1);
	E2 = E2 / (4*R2^2 + 2*E3*S - 4*R2*S^2);
	#
	x1 = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -R3, R4)));
	x1 = as.vector(x1);
	S  = rep(S, each=4); E2 = rep(E2, each=4);
	# robust:
	div = (x1^2*(2*x1-S) + (E2-R2)*x1);
	x3 = (R4 - x1^4 + x1^3*S - R2*x1^2) / div;
	xs = S - x1 - x3; x24 = R4 / (x1*x3);
	xd = sqrt(xs^2 - 4*x24);
	x2 = (xs + xd)/2; x4 = (xs - xd)/2;
	sol = cbind(x1=x1, x2=as.vector(x2), x3=as.vector(x3), x4=as.vector(x4))
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}

### Examples:

R = c(-2,-1,2,3)
sol = solve.S4HtM.E2P4(R);

test.S4HtMixed(sol, n=4)


### Ex 2:
R = c(-1,-1,0,3)
sol = solve.S4HtM.E2P4(R);

test.S4HtMixed(sol, n=4)


### Derivation:

### EP:
source("Polynomials.Helper.EP.R")
p = Epoly.base(n=4, v=4)[[4]]
p = sort.pm(p, "S")
print.pm(p, lead="S")

### S4:
p1 = toPoly.pm("S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 - 4*E4 - R1")
p2 = toPoly.pm("R2*E2^2 - (S*E3 + 2*R2^2)*E2 +
	+ S^2*E4 + S*R2*E3 + R2^3 - 4*R2*E4 + E3^2")
#
pR = solve.pm(p1, p2, "E2")
str(pR)
pR$Rez = sort.pm(pR$Rez, "S", xn2=c("E4", "E3", "R1"))
print.pm(pR$Rez, lead="S")
print.coeff(pR$Rez, "S")


########################
########################

####################
### E22a:        ###
### E2a: Order 2 ###
####################

###############
### Order 1 ###
###############

x1 + x2 + x3 + x4 - R1 # = 0
(x1*x2)^2 + (x2*x3)^2 + (x3*x4)^2 + (x4*x1)^2 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

E2a  = (x1*x2) + (x2*x3) + (x3*x4) + (x4*x1);
E22a = (x1*x2)^2 + (x2*x3)^2 + (x3*x4)^2 + (x4*x1)^2;

### E2: P[4]
(4*E4 - E22a)*E2^4 - 2*(E3^2 + S^2*E4)*E2^3 +
	+ (4*(S*E3 - 2*E4)*E22a - 8*S*E3*E4 + S^2*E3^2 + 2*E22a^2)*E2^2 +
	+ (4*S^3*E3*E4 + 4*S*E3^3 + 2*(S^2*E4 + E3^2)*E22a)*E2 +
	- S^4*E4^2 - 2*S^3*E3^3 - E22a^3 - E3^4 + 4*E22a^2*E4 + 2*S^2*E3^2*E4 +
		+ 8*S*E22a*E3*E4 - 4*S*E22a^2*E3 - 5*S^2*E22a*E3^2; # = 0


### Solver:
solve.S4HtM.E22P1 = function(R, sort=TRUE, all.sol=TRUE, debug=TRUE) {
	S = R[1]; E22a = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c((4*E4 - E22a),   - 2*(E3^2 + S^2*E4),
		(4*(S*E3 - 2*E4)*E22a - 8*S*E3*E4 + S^2*E3^2 + 2*E22a^2),
		(4*S^3*E3*E4 + 4*S*E3^3 + 2*(S^2*E4 + E3^2)*E22a),
		- S^4*E4^2 - 2*S^3*E3^3 - E22a^3 - E3^4 + 4*E22a^2*E4 + 2*S^2*E3^2*E4 +
			+ 8*S*E22a*E3*E4 - 4*S*E22a^2*E3 - 5*S^2*E22a*E3^2);
	E2 = roots(coeff);
	if(debug) print(E2);
	return(solve.S4HtM.E22Base(R, E2, sort=sort, all.sol=all.sol));
}

### Examples:

R = c(1,-1,2,1)
sol = solve.S4HtM.E22P1(R)

test.S4HtMixed(sol, n=1, nE2=2)


### Ex 2:
R = c(-2,-3,2,-1)
sol = solve.S4HtM.E22P1(R)

test.S4HtMixed(sol, n=1, nE2=2)


########################
########################

################
### E22a:    ###
### Order 2  ###
################

x1^2 + x2^2 + x3^2 + x4^2 - R1 # = 0
(x1*x2)^2 + (x2*x3)^2 + (x3*x4)^2 + (x4*x1)^2 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

E2a  = (x1*x2) + (x2*x3) + (x3*x4) + (x4*x1);
E22a = (x1*x2)^2 + (x2*x3)^2 + (x3*x4)^2 + (x4*x1)^2;

### E2: P[8]
# - see coeff();


### Solver
coeff.S4Ht.E22P2 = function(R) {
	R1 = R[1]; E22a = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(E22a, 0, 4*R1*E4 - 4*R1*E22a, - 16*E3*E22a,
		16*E4^2 - 12*R1^2*E4 + 16*E22a*E4 - 4*R1*E3^2 + 6*R1^2*E22a - 8*E22a^2,
		- 32*R1*E3*E4 + 32*R1*E3*E22a,
		- 32*E3^2*E4 + 12*R1^3*E4 - 48*R1*E22a*E4 + 8*R1^2*E3^2 + 64*E3^2*E22a - 4*R1^3*E22a + 16*R1*E22a^2,
		32*R1^2*E3*E4 - 128*E3*E22a*E4 + 32*R1*E3^3 - 16*R1^2*E3*E22a + 64*E3*E22a^2,
		- 4*R1^4*E4 + 32*R1^2*E22a*E4 - 64*E22a^2*E4 + 16*E3^4 - 4*R1^3*E3^2 + 16*R1*E3^2*E22a +
			+ R1^4*E22a - 8*R1^2*E22a^2 + 16*E22a^3);
	return(coeff);
}
solve.S4HtM.E22P2 = function(R, sort=TRUE, all.sol=FALSE, debug=TRUE) {
	coeff = coeff.S4Ht.E22P2(R);
	S = roots(coeff);
	if(debug) print(S);
	len = length(S);
	E2 = (S^2 - R[1])/2;
	#
	sol = lapply(seq(len), function(id) {
		RS = R;
		RS[1] = S[id];
		solve.S4HtM.E22Base(RS, E2[id], sort=sort, all.sol=all.sol)
	})
	sol = do.call(rbind, sol);
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}

### Examples:

### Ex 1:
R = c(3,-1,2,1)
sol = solve.S4HtM.E22P2(R)

test.S4HtMixed(sol, n=2, nE2=2)


### Ex 2:
R = c(-3,2,2,-1)
sol = solve.S4HtM.E22P2(R)

test.S4HtMixed(sol, n=2, nE2=2)


### Ex 3:
# E22a = 0: only 6*4 = 24 solutions;
# (*2 with all=T);
R = c(-3,0,1,2)
sol = solve.S4HtM.E22P2(R)

test.S4HtMixed(sol, n=2, nE2=2)


### Derivation:

pE2 = polyE2Ord2();
p1  = toPoly.pm("S^2 - 2*E2 - R1")

pR = solve.pm(p1, pE2, "E2")
pR$Rez$coeff = - pR$Rez$coeff;
pR$Rez = sort.pm(pR$Rez, xn="S", xn2=c("E4", "E3"))
print.pm(pR$Rez, lead="S")
print.coeff(pR$Rez, "S")


########################
########################

################
### E22a:    ###
### Order 3  ###
################

x1^3 + x2^3 + x3^3 + x4^3 - R1 # = 0
(x1*x2)^2 + (x2*x3)^2 + (x3*x4)^2 + (x4*x1)^2 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

E2a  = (x1*x2) + (x2*x3) + (x3*x4) + (x4*x1);
E22a = (x1*x2)^2 + (x2*x3)^2 + (x3*x4)^2 + (x4*x1)^2;

### E2: P[12]
# - see coeff();


### Solver
coeff.S4Ht.E22P3 = function(R) {
	R1 = R[1]; E22a = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(2*E4 + E22a, 0, - 3*E3^2,  - 30*E3*E4 - 2*R1*E4 - 24*E3*E22a - 4*R1*E22a,
		81*E4^2 + 18*E22a*E4 - 18*E22a^2,
		54*E3^3,  - 108*E3^2*E4 - 6*R1^2*E4 + 189*E3^2*E22a + 36*R1*E3*E22a + 6*R1^2*E22a,
		- 378*E3*E22a*E4 - 90*R1*E22a*E4 + 216*E3*E22a^2 + 36*R1*E22a^2,
		- 324*E22a^2*E4 - 162*E3^4 + 54*R1*E3^3 + 9*R1^2*E3^2 + 81*E22a^3,
		378*E3^3*E4 - 162*R1*E3^2*E4 - 18*R1^2*E3*E4 + 10*R1^3*E4 - 378*E3^3*E22a + 162*R1*E3^2*E22a - 4*R1^3*E22a,
		648*E3^2*E22a*E4 - 432*R1*E3*E22a*E4 + 72*R1^2*E22a*E4 - 162*E3^2*E22a^2 + 108*R1*E3*E22a^2 - 18*R1^2*E22a^2,
		162*E3^5 - 162*R1*E3^4 + 54*R1^2*E3^3 - 6*R1^3*E3^2,
		- 324*E3^4*E4 + 432*R1*E3^3*E4 - 216*R1^2*E3^2*E4 + 48*R1^3*E3*E4 - 4*R1^4*E4 + 81*E3^4*E22a +
			- 108*R1*E3^3*E22a + 54*R1^2*E3^2*E22a - 12*R1^3*E3*E22a + R1^4*E22a);
	return(coeff);
}
solve.S4HtM.E22P3 = function(R, sort=TRUE, all.sol=FALSE, debug=TRUE) {
	coeff = coeff.S4Ht.E22P3(R);
	S = roots(coeff);
	if(debug) print(S);
	len = length(S);
	E2 = (S^3 + 3*R[3] - R[1]) / (3*S);
	#
	sol = lapply(seq(len), function(id) {
		RS = R;
		RS[1] = S[id];
		solve.S4HtM.E22Base(RS, E2[id], sort=sort, all.sol=all.sol)
	})
	sol = do.call(rbind, sol);
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}

### Examples:

### Ex 1:
R = c(3,-1,2,1)
sol = solve.S4HtM.E22P3(R)

test.S4HtMixed(sol, n=3, nE2=2)


### Ex 2:
# 2*E4 + E22a = 0: only 10*4 = 40 solutions;
# (*2 with all=T);
R = c(-3,2,2,-1)
sol = solve.S4HtM.E22P3(R)

test.S4HtMixed(sol, n=3, nE2=2)


### Ex 3:
# 2*E4 + E22a = 0: only 10*4 = 40 solutions;
# (*2 with all=T);
R = c(-3,2,1,-1)
sol = solve.S4HtM.E22P3(R)

test.S4HtMixed(sol, n=3, nE2=2)


### Derivation:

pE2 = polyE2Ord2();
p1  = toPoly.pm("S^3 - 3*E2*S + 3*E3 - R1")

pR = solve.pm(p1, pE2, "E2")
pR$Rez$coeff = - pR$Rez$coeff;
pR$Rez = sort.pm(pR$Rez, xn="S", xn2=c("E4", "E3"))
print.pm(pR$Rez, lead="S")
print.coeff(pR$Rez, "S")


########################
########################
########################

####################
### Type: E121a  ###
### Order 1      ###
####################

x1 + x2 + x3 + x4 - R1 # = 0
x1*x2^2*x3 + x2*x3^2*x4 + x3*x4^2*x1 + x4*x1^2*x2 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

### Abbreviations:
E121a = x1*x2^2*x3 + x2*x3^2*x4 + x3*x4^2*x1 + x4*x1^2*x2;

### E2 vs E121a:
4*(4*E4 + E121a)*E4*E2^2 - (E121a + 8*E4)*(S^2*E4 + E3^2)*E2 +
	+ S^4*E4^2 + S*E121a^2*E3 - E121a^3 + E3^4 - 4*E121a^2*E4 + 2*S^2*E3^2*E4;


### Solver
coeff.S4Ht.E121P1 = function(R) {
	S = R[1]; E121a = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(4*(4*E4 + E121a)*E4,  - (E121a + 8*E4)*(S^2*E4 + E3^2),
		S^4*E4^2 + S*E121a^2*E3 - E121a^3 + E3^4 - 4*E121a^2*E4 + 2*S^2*E3^2*E4);
	return(coeff);
}
solve.S4HtM.E121P1 = function(R, sort=TRUE, all.sol=FALSE, debug=TRUE) {
	coeff = coeff.S4Ht.E121P1(R);
	E2 = roots(coeff);
	if(debug) print(E2);
	sol = solve.S4HtM.E121Base(R, E2=E2, sort=sort, all.sol=all.sol);
	return(sol)
}

### Examples:

### Ex 1:
R = c(3,-1,2,1)
sol = solve.S4HtM.E121P1(R)

test.S4HtMixed.En3(sol, n=1, nE=c(1,2,1))


### Ex 2:
R = c(-1,-3,2,2)
sol = solve.S4HtM.E121P1(R)

test.S4HtMixed.En3(sol, n=1, nE=c(1,2,1))


### Ex 3:
R = c(1,-1,3,-1)
sol = solve.S4HtM.E121P1(R)

test.S4HtMixed.En3(sol, n=1, nE=c(1,2,1))


### Ex 4:
# 4*E4 + E121a = 0;
R = c(1,4,3,-1)
sol = solve.S4HtM.E121P1(R)

test.S4HtMixed.En3(sol, n=1, nE=c(1,2,1))


########################
########################

####################
### Type: E121a  ###
### Order 2      ###
####################

x1^2 + x2^2 + x3^2 + x4^2 - R1 # = 0
x1*x2^2*x3 + x2*x3^2*x4 + x3*x4^2*x1 + x4*x1^2*x2 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:


### Eq S:
E4*(E121a + 2*E4)*S^4 - (E3^2*E121a + 3*R1*E4*E121a + 4*E3^2*E4 + 8*R1*E4^2)*S^2 +
	+ 2*E3*E121a^2*S + R1*E3^2*E121a + 2*R1^2*E4*E121a + 2*E3^4 - 2*E121a^3 + 8*R1^2*E4^2 + 8*R1*E3^2*E4 - 8*E4*E121a^2 # = 0


### Solver:
coeff.S4HtM.E121P2 = function(R) {
	R1 = R[1]; E121a = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(E4*(E121a + 2*E4),  0,  - (E3^2*E121a + 3*R1*E4*E121a + 4*E3^2*E4 + 8*R1*E4^2),
		2*E3*E121a^2,
		R1*E3^2*E121a + 2*R1^2*E4*E121a + 2*E3^4 - 2*E121a^3 +
			+ 8*R1^2*E4^2 + 8*R1*E3^2*E4 - 8*E4*E121a^2);
	return(coeff);
}
solve.S4HtM.E121P2 = function(R, sort=TRUE, all.sol=FALSE, debug=TRUE) {
	coeff = coeff.S4HtM.E121P2(R);
	S = roots(coeff);
	if(debug) print(S);
	#
	len = length(S);
	E2  = (S^2 - R[1]) / 2;
	sol = lapply(seq(len), function(id) {
			RS = R; RS[1] = S[id];
			solve.S4HtM.E121Base(RS, E2[id], sort=sort, all.sol=all.sol);
		});
	sol = do.call(rbind, sol);
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}

### Examples:

### Ex 1:
R = c(3,-1,2,1)
sol = solve.S4HtM.E121P2(R)

test.S4HtMixed.En3(sol, n=2, nE=c(1,2,1))


### Ex 2:
R = c(-1,-3,2,2)
sol = solve.S4HtM.E121P2(R)

test.S4HtMixed.En3(sol, n=2, nE=c(1,2,1))


round0(poly.calc(sol[,1]))


### Derivation:

p1 = toPoly.pm("S^2 - 2*E2 - R1");
p2 = polyE2_E121P1();

pR = solve.pm(p1, p2, "E2");
pR = sort.pm(pR$Rez, "S", xn2=c("E4", "E3", "R1"), sort.coeff=c(5, 6:8))
print.pm(pR, lead="S")


########################
########################

####################
### Type: E121a  ###
### Order 3      ###
####################

n = 3
x1^n + x2^n + x3^n + x4^n - R1 # = 0
x1*x2^2*x3 + x2*x3^2*x4 + x3*x4^2*x1 + x4*x1^2*x2 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:


### Eq S:
E4*(E4 + E121a)*S^6 - 3*E3^2*(E121a + 2*E4)*S^4 +
	- (5*R1*E4*E121a + 8*R1*E4^2 - 9*E3*E121a^2 - 15*E3*E4*E121a - 24*E3*E4^2)*S^3 +
	+ 9*(E3^4 - E121a^3 - 4*E4*E121a^2)*S^2 +
	+ 3*E3^2*(R1 - 3*E3)*(E121a + 8*E4)*S +
	+ 4*R1^2*E4*E121a + 16*R1^2*E4^2 - 24*R1*E3*E4*E121a + 36*E3^2*E4*E121a - 96*R1*E3*E4^2 + 144*E3^2*E4^2 # = 0


### Solver:
coeff.S4HtM.E121P3 = function(R) {
	R1 = R[1]; E121a = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(E4*(E4 + E121a), 0, - 3*E3^2*(E121a + 2*E4),
		- (5*R1*E4*E121a + 8*R1*E4^2 - 9*E3*E121a^2 - 15*E3*E4*E121a - 24*E3*E4^2),
		9*(E3^4 - E121a^3 - 4*E4*E121a^2), 3*E3^2*(R1 - 3*E3)*(E121a + 8*E4),
		4*R1^2*E4*E121a + 16*R1^2*E4^2 - 24*R1*E3*E4*E121a + 36*E3^2*E4*E121a +
			- 96*R1*E3*E4^2 + 144*E3^2*E4^2);
	return(coeff);
}
solve.S4HtM.E121P3 = function(R, sort=TRUE, all.sol=FALSE, debug=TRUE) {
	coeff = coeff.S4HtM.E121P3(R);
	S = roots(coeff);
	if(debug) print(S);
	#
	len = length(S);
	E3 = R[3];
	E2  = (S^3 + 3*E3 - R[1]) / (3*S);
	sol = lapply(seq(len), function(id) {
			RS = R; RS[1] = S[id];
			solve.S4HtM.E121Base(RS, E2[id], sort=sort, all.sol=all.sol);
		});
	sol = do.call(rbind, sol);
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}

### Examples:

### Ex 1:
# R121a + E4 = 0
R = c(3,-1,2,1)
sol = solve.S4HtM.E121P3(R)

test.S4HtMixed.En3(sol, n=3, nE=c(1,2,1))


### Ex 2:
R = c(-1,-3,-2,2)
sol = solve.S4HtM.E121P3(R)

test.S4HtMixed.En3(sol, n=3, nE=c(1,2,1))


round0(poly.calc(sol[,1]))


### Derivation:

p1 = toPoly.pm("S^3 - 3*E2*S + 3*E3 - R1");
p2 = polyE2_E121P1();

pR = solve.pm(p1, p2, "E2");
pR = sort.pm(pR$Rez, "S", xn2=c("E4", "E3", "R1"), sort.coeff=c(5, 6:8))
print.pm(pR, lead="S")


########################
########################
########################

####################
### Type: E212a  ###
####################

###############
### Order 1 ###
###############

x1 + x2 + x3 + x4 - R1 # = 0
x1^2*x2*x3^2 + x2^2*x3*x4^2 + x3^2*x4*x1^2 + x4^2*x1*x2^2 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

### Eq: for E2
E3*((E4*S + E212a)^2 - 4*E4*E3^2)*E2 +
	- S^3*E4^3 + E3^5 - S*E3^3*E212a + 4*S*E4^2*E3^2 - 3*S^2*E4^2*E212a +
		+ 4*E4*E3^2*E212a - 3*S*E4*E212a^2 - E212a^3 # = 0

### Solver:
coeff.S4HtM.E212P1 = function(R) {
	# coefficients for E2;
	S = R[1]; E212a = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(E3*((E4*S + E212a)^2 - 4*E4*E3^2),
		- S^3*E4^3 + E3^5 - S*E3^3*E212a + 4*S*E4^2*E3^2 - 3*S^2*E4^2*E212a +
			+ 4*E4*E3^2*E212a - 3*S*E4*E212a^2 - E212a^3);
	return(coeff);
}
solve.S4HtM.E212P1 = function(R, sort=TRUE, all.sol=FALSE, debug=TRUE) {
	coeff = coeff.S4HtM.E212P1(R);
	if(coeff[1] == 0) stop("No solution!");
	E2 = - coeff[2] / coeff[1];
	if(debug) print(E2);
	#
	S = R[1]; E212a = R[2]; E3 = R[3]; E4 = R[4];
	E2b = (E212a + E4*S) / E3;
	E2a = E2 - E2b;
	# robust based on (x1 + x3):
	xs  = roots(c(1, -S, E2a));
	x13 = (E3 - E2b*xs) / (S - 2*xs);
	xd = sqrt(xs^2 - 4*x13 + 0i);
	x1 = (xs + xd)/2; x3 = (xs - xd)/2;
	x24 = E2b - x13; xs = S - xs;
	xd = sqrt(xs^2 - 4*x24 + 0i);
	x2 = (xs + xd)/2; x4 = (xs - xd)/2;
	sol = cbind(x1=x1, x2=x2, x3=x3, x4=x4);
	# TODO: all.sol + permutations;
	return(sol)
}

### Examples:

### Ex 1:
R = c(3,-1,2,1)
sol = solve.S4HtM.E212P1(R)

test.S4HtMixed.En3(sol, n=1, nE=c(2,1,2))


### Ex 2:
R = c(5,2,3,-1)
sol = solve.S4HtM.E212P1(R)

test.S4HtMixed.En3(sol, n=1, nE=c(2,1,2))


###############
### Derivation:
E212a - E3*E2b + E4*S # = 0
# =>
E212a - E3*E2 + E3*E2a + E4*S # = 0

p1 = polyE2a()
p2 = toPoly.pm("E212a - E3*E2 + E3*E2a + E4*S")
pR = solve.pm(p2, p1, "E2a")
str(pR)
pR = pR$Rez; pR$coeff = - pR$coeff;
pR = sort.pm(pR, "E2")
print.pm(pR, lead="E2")


########################
########################
########################

###################
### Type: E21a  ###
###################

###############
### Order 1 ###
###############

x1 + x2 + x3 + x4 - R1 # = 0
x1^2*x2 + x2^2*x3 + x3^2*x4 + x4^2*x1 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

# TODO: complicated;


### Derivation:

E21b - (E2 - E2a)*S + E3 # = 0
# ???
# Epoly.distinct(c(2,1), v=4)
E21 - S*E2 + 3*E3 = 0
 
