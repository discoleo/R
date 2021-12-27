########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Hetero-Symmetric S4: Mixed
### Basic Types
###
### draft v.0.1k-test


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


####################
####################

### Helper Functions

source("Polynomials.Helper.R")


### Other

test.S4HtMixed = function(sol, n=2, nE2 = 1, R = NULL) {
	x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];
	err1 = x1^n + x2^n + x3^n + x4^n;
	# Ht
	ht = cbind(x1*x2, x2*x3, x3*x4, x4*x1);
	ht = if(nE2 == 1) ht else ht^nE2;
	err2 = apply(ht, 1, sum);
	err3 = x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4;
	err4 = x1*x2*x3*x4;
	err = rbind(err1, err2, err3, err4);
	if( ! is.null(R)) {
		for(id in 1:4) err[id,] = err[id,] - R[id];
	}
	err = round0(err);
	return(err);
}

test.S4HtMixed.En3 = function(sol, R=NULL, n=2, nE=c(1,2,1)) {
	# Ht
	ht.f = function(x) {
		sum(x^nE[1] * (x^nE[2])[c(2,3,4,1)] * (x^nE[3])[c(3,4,1,2)]);
	}
	err2 = apply(sol, 1, ht.f);
	#
	x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];
	err1 = x1^n + x2^n + x3^n + x4^n;
	err3 = x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4;
	err4 = x1*x2*x3*x4;
	err = rbind(err1, err2, err3, err4);
	if( ! is.null(R)) {
		for(id in 1:4) err[id,] = err[id,] - R[id];
	}
	err = round0(err);
	return(err);
}

e2.f = function(x) {
	e2.f0 = function(x) x[1]*sum(x, -x[1]) + x[2]*(x[3]+x[4]) + x[3]*x[4];
	e2 = if(is.matrix(x)) apply(x, 1, e2.f0) else e2.f0(x);
	sort.sol(matrix(e2, ncol=1), useRe=TRUE);
}
e3.f = function(x) {
	e3.f0 = function(x) {
		x34 = x[3]*x[4];
		x[1]*(x[2]*(x[3] + x[4]) + x34) + x[2]*x34;
	}
	e3 = if(is.matrix(x)) apply(x, 1, e3.f0) else e3.f0(x);
	sort.sol(matrix(e3, ncol=1), useRe=TRUE);
}
e2a.f = function(x) {
	e2.f0 = function(x) sum(x * x[c(2,3,4,1)]);
	e2 = if(is.matrix(x)) apply(x, 1, e2.f0) else e2.f0(x);
	sort.sol(matrix(e2, ncol=1), useRe=TRUE);
}
e3a.f = function(x, pow=c(1,2,1)) {
	e2.f0 = function(x) sum(x^pow[1] * x[c(2,3,4,1)]^pow[2] * x[c(3,4,1,2)]^pow[3]);
	e2 = if(is.matrix(x)) apply(x, 1, e2.f0) else e2.f0(x);
	sort.sol(matrix(e2, ncol=1), useRe=TRUE);
}

### Solve Coefficients:
# - hack the formulas;
which.sq = function(x, sq=2, iter=1000, digits=6, pow=2) {
	if(is.na(x)) return(NA);
	if(round(x) == round(x, digits)) return(0);
	i = seq(iter);
	sq = if(pow == 2) sqrt(sq) else rootn(sq, n=pow);
	d = round(i * sq - x, digits);
	id = which(d == round(d));
	if(length(id) > 0) return(id);
	d = round(i * sq + x, digits);
	id = which(d == round(d));
	if(length(id) == 0) return(NA);
	return(- id);
}

### Formulas

polyE2Ord1 = function() {
	p = toPoly.pm("E2a*E2^2 - (S*E3 + 2*E2a^2)*E2 +
		+ S^2*E4 + S*E2a*E3 + E2a^3 - 4*E2a*E4 + E3^2");
	return(p);
}
polyE2Ord2 = function() {
	p = toPoly.pm("(4*E4 - E22a)*E2^4 - 2*(E3^2 + S^2*E4)*E2^3 +
	+ (4*(S*E3 - 2*E4)*E22a - 8*S*E3*E4 + S^2*E3^2 + 2*E22a^2)*E2^2 +
	+ (4*S^3*E3*E4 + 4*S*E3^3 + 2*(S^2*E4 + E3^2)*E22a)*E2 +
	- S^4*E4^2 - 2*S^3*E3^3 - E22a^3 - E3^4 + 4*E22a^2*E4 + 2*S^2*E3^2*E4 +
		+ 8*S*E22a*E3*E4 - 4*S*E22a^2*E3 - 5*S^2*E22a*E3^2");
	return(p);
}

###############
###############

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
solve.S4Ht.P1 = function(R, sort=TRUE, debug=TRUE) {
	S = R[1]; E2a = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(E2a, - (S*E3 + 2*E2a^2),
		S^2*E4 + S*E2a*E3 + E2a^3 - 4*E2a*E4 + E3^2);
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
sol = solve.S4Ht.P1(R)

test.S4HtMixed(sol, n=1)

### Classic Poly:
round0(poly.calc(sol[,1]) * abs(R[2]))


### Derivation:

# - classic approach: P[2] o P[2];

### E2a:
(x1+x3)*(S - x1 - x3) - R2 # = 0
x13s^2 - S*x13s + R2 # = 0

### E3 =>
x1*x3*(S - x1 - x3) + x2*x4*(x1+x3) - R3 # = 0
(x1*x3)^2*(S - x1 - x3) - R3*x1*x3 + R4*(x1+x3) # = 0


### Solution: based on "classic" approach
solve.S4Ht.P1old = function(R, debug=FALSE) {
	xs  = roots(c(1, -R[1], R[2]));
	x13 = sapply(seq(length(xs)), function(id) roots(c(R[1] - xs[id], -R[3], R[4]*xs[id])));
	xs = rep(xs, each=2); x13 = as.vector(x13);
	xd = sqrt(xs^2 - 4*x13 + 0i);
	x1 = (xs + xd)/2; x3 = (xs - xd)/2;
	# x2, x4:
	xs = R[1] - xs; x24 = R[4] / x13;
	xd = sqrt(xs^2 - 4*x24 + 0i);
	x2 = (xs + xd)/2; x4 = (xs - xd)/2;
	sol = cbind(x1, x2, x3, x4)
	return(sol)
}
e2.f = function(x) x[1]*sum(x, -x[1]) + x[2]*(x[3]+x[4]) + x[3]*x[4]

###
R = c(1,-1,2,3)
sol = solve.S4Ht.P1old(R)
round0(poly.calc(apply(sol, 1, e2.f)[1:2]) * R[2])
R[1]^2*R[4] + R[1]*R[2]*R[3] + R[2]^3 - 4*R[2]*R[4] + R[3]^2


R[2]*x^2 - (R[1]*R[3] + 2*R[2]^2)*x +
	+ R[1]^2*R[4] + R[1]*R[2]*R[3] + R[2]^3 - 4*R[2]*R[4] + R[3]^2

test.S4HtMixed(sol, n=1)


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

### E3 =>
x1*x3*(x2 + x4) + x2*x4*(x1 + x3) - R3 # = 0
# * (x1 + x3) =>
R2*x1*x3 + x2*x4*(x1+x3)^2 - R3*(x1+x3) # = 0
# * x1*x3 =>
R2*(x1*x3)^2 + R4*(x1+x3)^2 - R3*x1*x3*(x1+x3) # = 0

### S2 + 2*E2b =>
(x1+x3)^2 + (x2+x4)^2 - R1 - 2*E2b # = 0
x1*x3*(x1+x3)^2 + x1*x3*(x2+x4)^2 - R1*x1*x3 - 2*((x1*x3)^2 + R4) # = 0
x1*x3*(x1+x3)^4 + x1*x3*R2^2 - R1*x1*x3*(x1+x3)^2 - 2*((x1*x3)^2 + R4)*(x1+x3)^2 # = 0

# TODO: derive properly the solution based on S^4;

### P[4]
R2*S^4 - 2*R3*S^3 + 2*(2*R4 - 2*R2^2 - R1*R2)*S^2 +
	+ 2*R3*(R1+2*R2)*S - 16*R2*R4 + (4*R2^3 + 4*R1*R2^2 + R1^2*R2) + 4*R3^2 # = 0


### Solver:
solve.S4Ht.P2 = function(R, sort=TRUE, debug=TRUE) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	# S^4:
	coeff = c(R2, - 2*R3, 2*(2*R4 - 2*R2^2 - R1*R2),
		2*R3*(R1+2*R2), - 16*R2*R4 + (4*R2^3 + 4*R1*R2^2 + R1^2*R2) + 4*R3^2);
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

R = c(1,-1,2,3)
sol = solve.S4Ht.P2(R)

test.S4HtMixed(sol)


### Ex 2:
R = c(0,-2,2,1)
sol = solve.S4Ht.P2(R)

test.S4HtMixed(sol)


###############

### Derivation:
R2*y^2 + R4*x^2 - R3*x*y # = 0
y*x^4 + y*R2^2 - R1*y*x^2 - 2*(y^2 + R4)*x^2 # = 0

p1 = toPoly.pm("R2*y^2 + R4*x^2 - R3*x*y")
p2 = toPoly.pm("y*x^4 + y*R2^2 - R1*y*x^2 - 2*(y^2 + R4)*x^2")
pR = solve.pm(p1, p2, "y")
pR$Rez = sort.pm(pR$Rez, "x")
print.pm(pR$Rez, lead="x")

### P[8]
# - is NOT the minimal polynomial!
R2*x^8 - 2*R3*x^7 - 2*R2*R1*x^6 + 4*R4*x^6 - 2*R3*R2*x^5 + 2*R3*R1*x^5 +
	+ 2*R2^3*x^4 + R2*R1^2*x^4 + 4*R3^2*x^4 - 8*R4*R2*x^4 - 2*R3*R2^2*x^3 + 2*R3*R2*R1*x^3 +
	- 2*R2^3*R1*x^2 + 4*R4*R2^2*x^2 - 2*R3*R2^3*x + R2^5

coeffs.S4 = function(R) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	coeff = c(R2, - 2*R3, - 2*R2*R1 + 4*R4,
		- 2*R3*R2 + 2*R3*R1,
		2*R2^3 + R2*R1^2 + 4*R3^2 - 8*R4*R2,
		- 2*R3*R2^2 + 2*R3*R2*R1,
		- 2*R2^3*R1 + 4*R4*R2^2,
		- 2*R3*R2^3,
		R2^5);
	return(coeff);
}
x13.S4 = function(x, R) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	x0 = 2*x^4*R4 - 2*x^2*R4*R2;
	div = - x^4*R2 + x^2*R1*R2 - R2^3 + 2*x^3*R3;
	return( x0 / div);
}


### Solution: based on P[8]
solve.S4Ht.P2old = function(R, debug=FALSE) {
	coeff = coeffs.S4(R)
	xs = roots(coeff);
	x13 = x13.S4(xs, R);
	xd = sqrt(xs^2 - 4*x13 + 0i);
	x1 = (xs + xd)/2; x3 = (xs - xd)/2;
	xs = R[2] / xs; x24 = R[4] / x13;
	xd = sqrt(xs^2 - 4*x24 + 0i);
	x2 = (xs + xd)/2; x4 = (xs - xd)/2;
	sol = cbind(x1, x2, x3, x4)
	return(sol)
}
### Explore Coefficient of x^2
coeffx2 = function(R) 2*(2*R[4] - 2*R[2]^2 - R[1]*R[2]);

###
R = c(1,-1,2,3)
sol = solve.S4Ht(R)
s = apply(sol, 1, sum) # TODO: proper solution;
s = sort.sol(matrix(s, ncol=1), mod.first=FALSE)
s
round0(poly.calc(s[c(1,3,5,7)])) * R[2]
coeffx2(R)

sol
test.S4HtMixed(sol)

-3  + 6*x - 2*x^2 - 2*x^3 + x^4 # R = c(1,1,1,1)
-12       + 4*x^2 - 2*x^3 + x^4 # R = c(-2,1,1,1)
-11 - 2*x + 6*x^2 - 2*x^3 + x^4 # R = c(-3,1,1,1)
-9 +  3*x + 4*x^2 +   x^3 + x^4 # R = c(1,-2,1,1)
 9 - 12*x - 2*x^2 + 4*x^3 + x^4 # R = c(1,1,-2,1)
45 + 6*x - 14*x^2 - 2*x^3 + x^4 # R = c(1,1,1,-2)
-63 + 4*x - 10*x^2 + 4*x^3 + x^4 # R = c(1,-1,2,3)

R2*x^4 - 2*R3*x^3 + 2*(2*R4 - 2*R2^2 - R1*R2)*x^2 +
	+ 2*R[3]*(R[1]+2*R[2])*x - 16*R[2]*R[4] + (4*R[2]^3 + 4*R[1]*R[2]^2 + R[1]^2*R[2]) + 4*R[3]^2;

R = c(1,1,1,1)
coeffs(1, R)

coeffs = function(r.id=1, R, c0=c(-3,-2,-1,1,2,3,4,5)) {
	sapply(c0, function(Rv) {
	Rm = R; Rm[r.id] = Rv;
	sol = solve.S4Ht(Rm);
	s = apply(sol, 1, sum)
	s = sort.sol(matrix(s, ncol=1), mod.first=FALSE)
	(round0(poly.calc(s[c(1,3,5,7)])) * Rm[2])[1] # ID
})
}


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
coeff.S4Ht.P3 = function(R) {
	R1 = R[1]; R2 = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(R2, - 3*E3, 9*E4 - 6*R2^2, 15*E3*R2 - 2*R1*R2,
		- 36*E4*R2 + 3*R1*E3 + 9*R2^3,
		- 18*E3*R2^2 + 6*R1*R2^2,
		9*E3^2*R2 - 6*R1*E3*R2 + R1^2*R2);
	return(coeff);
}
solve.S4Ht.P3 = function(R, sort=TRUE, debug=TRUE) {
	coeff = coeff.S4Ht.P3(R);
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
sol = solve.S4Ht.P3(R);

test.S4HtMixed(sol, n=3)


### Ex 2:
R = c(-2,0,3,1)
sol = solve.S4Ht.P3(R);

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
coeff.S4Ht.P4 = function(R) {
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
solve.S4Ht.P4 = function(R, sort=TRUE, debug=TRUE) {
	coeff = coeff.S4Ht.P4(R);
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
sol = solve.S4Ht.P4(R);

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
solve.S4HtM.Ord2.P1 = function(R, sort=TRUE, all.sol=TRUE, debug=TRUE) {
	S = R[1]; E22a = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c((4*E4 - E22a),   - 2*(E3^2 + S^2*E4),
		(4*(S*E3 - 2*E4)*E22a - 8*S*E3*E4 + S^2*E3^2 + 2*E22a^2),
		(4*S^3*E3*E4 + 4*S*E3^3 + 2*(S^2*E4 + E3^2)*E22a),
		- S^4*E4^2 - 2*S^3*E3^3 - E22a^3 - E3^4 + 4*E22a^2*E4 + 2*S^2*E3^2*E4 +
			+ 8*S*E22a*E3*E4 - 4*S*E22a^2*E3 - 5*S^2*E22a*E3^2);
	E2 = roots(coeff);
	if(debug) print(E2);
	return(solve.S4HtM.Ord2Base(R, E2, sort=sort, all.sol=all.sol));
}

solve.S4HtM.Ord2Base = function(R, E2, sort=TRUE, all.sol=TRUE) {
	S = R[1]; E22a = R[2]; E3 = R[3]; E4 = R[4];
	#
	len = length(E2);
	x1 = sapply(seq(len), function(id) roots(c(1, -S, E2[id], -E3, E4)));
	x1 = as.vector(x1);
	E2 = rep(E2, each=4);
	len = length(x1);
	# fully robust:
	x13T1 = x1^2*(E2^2 - E22a - 2*S*E3 + 2*E4);
	x3sq = x1^6*(x1^2 + 2*E2 - S^2) + x1^4*E22a - E4^2;
	div =  x1^4*(2*x1^2 + 2*E2 - S^2) + x13T1;
	x3sq = - x3sq / div;
	x3 = (x3sq*(S - x1) + E4/x1) / (x3sq + E2 - x1*(S - x1));
	x3 = as.vector(x3);
	# x2, x4:
	xs = S - x1 - x3; x24 = E4 / (x1*x3);
	# Note: both sqrt() values are valid;
	xd = sqrt(xs^2 - 4*x24 + 0i);
	x2 = (xs + xd)/2; x4 = (xs - xd)/2;
	sol = cbind(x1=x1, x2=x2, x3=x3, x4=x4);
	if(all.sol) sol = rbind(sol, sol[, c(1,4,3,2)]); # all roots
	#
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}

### Examples:

R = c(1,-1,2,1)
sol = solve.S4HtM.Ord2.P1(R)

test.S4HtMixed(sol, n=1, nE2=2)


### Ex 2:
R = c(-2,-3,2,-1)
sol = solve.S4HtM.Ord2.P1(R)

test.S4HtMixed(sol, n=1, nE2=2)


###############
### Derivation:

# - classic approach: P[2] o P[2];

### E2a:
# let: xs = x1 + x3 => x2 + x4 = S - xs;
(x1^2 + x3^2)*((S - xs)^2 - 2*x2*x4) - R2 # = 0
(xs^2 - 2*x1*x3)*(x1*x3*(S - xs)^2 - 2*R4) - R2*x1*x3 # = 0


### E3 =>
x1*x3*(S - x1 - x3) + x2*x4*(x1+x3) - R3 # = 0
(x1*x3)^2*(S - xs) - R3*x1*x3 + R4*xs # = 0


### Solution: based on "classic" approach
coeff.S4HtM.Ord2.P1 = function(R) {
	S = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	coeff = c(1, - 4*S, 6*S^2, - 4*S^3, S^4 - 2*S*R3 - 8*R4 - 2*R2,
		4*S^2*R3 + 16*R4*S + 4*R2*S,
		- 2*S^3*R3 - 12*R4*S^2 - 2*R2*S^2 - 4*R3^2,
		4*R4*S^3 + 4*S*R3^2,
		- 8*R4*S*R3 + 2*R2*S*R3 + 16*R4^2 - 8*R2*R4 + R2^2);
	return(coeff);
}
solve.S4HtM.Ord2.P1old = function(R, debug=TRUE) {
	coeff = coeff.S4HtM.Ord2.P1(R);
	xs  = roots(coeff);
	if(debug) print(xs);
	#
	S = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	x13 = 4*R4*xs^3 - 6*R4*xs^2*S + 2*R4*xs*S^2;
	div = - R2*xs + 4*R4*xs + xs^5 + R2*S - 4*R4*S - 3*xs^4*S + 3*xs^3*S^2 - xs^2*S^3 +
		+ 2*xs^2*R3 - 4*xs*S*R3 + 2*S^2*R3;
	x13 = x13 / div;
	#
	xd = sqrt(xs^2 - 4*x13 + 0i);
	x1 = (xs + xd)/2; x3 = (xs - xd)/2;
	# x2, x4:
	xs = R[1] - xs; x24 = R[4] / x13;
	xd = sqrt(xs^2 - 4*x24 + 0i);
	x2 = (xs + xd)/2; x4 = (xs - xd)/2;
	sol = cbind(x1, x2, x3, x4)
	return(sol)
}

###
R = c(1,-1,2,3)
sol = solve.S4HtM.Ord2.P1old(R)

test.S4HtMixed(sol, n=1, nE2=2)


round0(poly.calc(e2.f(sol)[c(1,3,5,7)])) * (4*R[4] - R[2])


R = c(1,  1,  1,  2,   7)
R = c(1,  1,  1,  3,  11)
R = c(1,  1,  1,  4,  15)
R = c(1,  2,  1,  2,   6)
R = c(1,  2,  1,  3,  10)
R = c(1,  2,  1,  4,  14)
R = c(1, -7,  1,  9,  43)
R = c(2, -3,  3, -5,  17)
R = c(2, -3,  3, sqrt(2)) # - 24
R = c(2, -3,  3, 2^(1/3)) # - 24
R = c(2, -3,  3, sqrt(3)) # - 24
R = c(2, -3,  3, 3^(1/3)) # - 24
R = c(2, -3,  3, sqrt(5)) # - 24
#
R = c(1, sqrt(2), 1, sqrt(2)) # -4
R = c(1, sqrt(2), 1, -sqrt(2)) # 12
R = c(1, -sqrt(2), 1, sqrt(2)) # -12
R = c(1, -sqrt(2), 3, sqrt(2)) # -36
R = c(1, -sqrt(3), 3, sqrt(3)) # -36
R = c(1, -sqrt(3), 5, sqrt(3)) # -60
#
R = c(1, 2*sqrt(3), 5, 1) # 24
which.coeff(R, sq=3)
#
which.coeff(R <- c(1, 2*sqrt(3), 5, 1), sq=3)
which.coeff(R <- c(1, 2*2^(1/3), 5, 1), sq=2^2, pow=3, DIFF=DIFF)

sol = solve.S4HtM.Ord2.P1old(R)
round0(poly.calc(e2.f(sol)[c(1,3,5,7)])) * (4*R[4] - R[2])


(4*R[4] - R[2])*x^4 - 2*(R[3]^2 + R[1]^2*R[4])*x^3 +
	+ (4*(R[1]*R[3] - 2*R[4])*R[2] - 8*R[1]*R[3]*R[4] + R[1]^2*R[3]^2 + 2*R[2]^2)*x^2 +
	+ (4*R[1]^3*R[3]*R[4] + 4*R[1]*R[3]^3 + 2*(R[1]^2*R[4] + R[3]^2)*R[2])*x +
	- R[1]^4*R[4]^2 - 2*R[1]^3*R[3]^3 - R[2]^3 - R[3]^4 + 4*R[2]^2*R[4] + 2*R[1]^2*R[3]^2*R[4] +
		+ 8*R[1]*R[2]*R[3]*R[4] - 4*R[1]*R[2]^2*R[3] - 5*R[1]^2*R[2]*R[3]^2;

whichHasPower(4, id=2, type=2)
whichHasPower(R <- c(1,1,-2,-3), id=2, type=2)
polyR(R)
which.sq(DIFF(R), sq=2)


### Solve Coefficient
whichHasPower = function(R, id=2, type=1, FUN=NULL, print=FALSE, digits=5, iter=1000) {
	if(length(R) == 1) {len = R; R0 = rep(1, len); }
	else {len = length(R); R0 = R; }
	vals = c(2,3,5); vsqrt = sqrt(vals);
	if(type == 1) {
		vsqrt = c(vsqrt, - vsqrt); vals = c(vals, vals);
	} else if(type == 2) {
		vsqrt = c(vsqrt, 2*vsqrt); vals = c(vals, vals);
	} else if(type == 3) {
		vsqrt = c(vsqrt, - vsqrt) + 1; vals = c(vals, vals);
	} else if(type == 4) {
		vsqrt = c(vsqrt, 2*vsqrt) - 1; vals = c(vals, vals);
	}
	VLEN = length(vals);
	m = array(NA, c(len, VLEN));
	f0 = if( ! is.null(FUN)) {
		function(vid, nr) {
			R = R0;
			R[nr] = vsqrt[vid];
			which.sq(FUN(R), sq=vals[vid], pow=2, digits=digits, iter=iter)
		}
	} else function(vid, nr) {
		R = R0;
		R[nr] = vsqrt[vid];
		which.coeff(R, sq=vals[vid], id=id, pow=2, digits=digits, iter=iter, print=print)
	}
	for(nr in seq(len)) {
		tmp = sapply(seq(VLEN), f0, nr);
		m[nr, ] = tmp;
	}
	return(m);
}
polyR = function(R) {
	sol = solve.S4HtM.Ord2.P1old(R, debug=FALSE)
	p = round0(poly.calc(e2.f(sol)[c(1,3,5,7)])) * (4*R[4] - R[2]);
	return(p);
}
which.coeff = function(R, sq=2, id=3, pow=2, DIFF=NULL, print=TRUE, digits=6, iter=1000) {
	sol = solve.S4HtM.Ord2.P1old(R, debug=FALSE)
	p = round0(poly.calc(e2.f(sol)[c(1,3,5,7)])) * (4*R[4] - R[2]);
	if(print) print(p);
	x = p[id];
	if( ! is.null(DIFF)) {
		x = x - DIFF(R);
	}
	return(which.sq(x, sq=sq, pow=pow, digits=digits, iter=iter))
}

# S^2:
DIFF = function(R) 4*(R[1]*R[3] - 2*R[4])*R[2] - 8*R[1]*R[3]*R[4];
# S^1:
DIFF = function(R) 4*R[1]^3*R[3]*R[4] + 4*R[1]*R[3]^3 + 2*(R[1]^2*R[4] + R[3]^2)*R[2];
# S^0:
DIFF = function(R) - R[2]^3 + 4*R[2]^2*R[4] - 2*R[1]^3*R[3]^3 +
	+ 2*R[1]^2*R[3]^2*R[4] - R[1]^4*R[4]^2 - R[3]^4 + 8*R[1]*R[2]*R[3]*R[4] +
	- 4*R[1]*R[2]^2*R[3] - 5*R[1]^2*R[2]*R[3]^2;

###
p1 = toPoly.pm("(xs^2 - 2*x13)*(x13*(S - xs)^2 - 2*R4) - R2*x13")
p2 = toPoly.pm("x13^2*(S - xs) - R3*x13 + R4*xs")
#
pR = solve.pm(p2, p1, "x13")
str(pR)
pR = div.pm(pR$Rez, toPoly.pm("xs^2 - 2*S*xs + S^2"), "xs")
pR$Rez = sort.pm(pR$Rez, "xs", xn2=c("S", "R4", "R3"))
print.pm(pR$Rez, lead="xs")
print.coeff(pR$Rez, "xs")


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
coeff.S4Ht.Ord2P2 = function(R) {
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
solve.S4HtM.Ord2P2 = function(R, sort=TRUE, all.sol=FALSE, debug=TRUE) {
	coeff = coeff.S4Ht.Ord2P2(R);
	S = roots(coeff);
	if(debug) print(S);
	len = length(S);
	E2 = (S^2 - R[1])/2;
	#
	sol = lapply(seq(len), function(id) {
		RS = R;
		RS[1] = S[id];
		solve.S4HtM.Ord2Base(RS, E2[id], sort=sort, all.sol=all.sol)
	})
	sol = do.call(rbind, sol);
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}

### Examples:

### Ex 1:
R = c(3,-1,2,1)
sol = solve.S4HtM.Ord2P2(R)

test.S4HtMixed(sol, n=2, nE2=2)


### Ex 2:
R = c(-3,2,2,-1)
sol = solve.S4HtM.Ord2P2(R)

test.S4HtMixed(sol, n=2, nE2=2)


### Ex 3:
# E22a = 0: only 6*4 = 24 solutions;
# (*2 with all=T);
R = c(-3,0,1,2)
sol = solve.S4HtM.Ord2P2(R)

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
coeff.S4Ht.Ord2P3 = function(R) {
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
solve.S4HtM.Ord2P3 = function(R, sort=TRUE, all.sol=FALSE, debug=TRUE) {
	coeff = coeff.S4Ht.Ord2P3(R);
	S = roots(coeff);
	if(debug) print(S);
	len = length(S);
	E2 = (S^3 + 3*R[3] - R[1]) / (3*S);
	#
	sol = lapply(seq(len), function(id) {
		RS = R;
		RS[1] = S[id];
		solve.S4HtM.Ord2Base(RS, E2[id], sort=sort, all.sol=all.sol)
	})
	sol = do.call(rbind, sol);
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}

### Examples:

### Ex 1:
R = c(3,-1,2,1)
sol = solve.S4HtM.Ord2P3(R)

test.S4HtMixed(sol, n=3, nE2=2)


### Ex 2:
R = c(-3,2,2,-1)
sol = solve.S4HtM.Ord2P3(R)

test.S4HtMixed(sol, n=3, nE2=2)


### Ex 3:
# 2*E4 + E22a = 0: only 10*4 = 40 solutions;
# (*2 with all=T);
R = c(-3,2,1,-1)
sol = solve.S4HtM.Ord2P3(R)

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

################
### E121a:   ###
### Order 1  ###
################

x1 + x2 + x3 + x4 - R1 # = 0
x1*x2^2*x3 + x2*x3^2*x4 + x3*x4^2*x1 + x4*x1^2*x2 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

### Abbreviations:
E121a = x1*x2^2*x3 + x2*x3^2*x4 + x3*x4^2*x1 + x4*x1^2*x2;


### Solver
coeff.S4Ht.E121aP1 = function(R) {
	p = replace.pm(pR$Rez, c(S = R[1], E121a = R[2], E3 = R[3], E4 = R[4]));
	coeff = coef.pm(p, "E2");
	return(coeff);
}
E22a.E121aP1 = function(R, E2) {
	len = length(E2);
	# c(E3, S, E4, E2, E121a)
	e22a = sapply(seq(len),
		function(id) eval.pm(pR$x0, list(E3=R[3], S=R[1], E4=R[4], E2=E2[id], E121a=R[2])));
	eDiv = sapply(seq(len),
		function(id) eval.pm(pR$div, list(E3=R[3], S=R[1], E4=R[4], E2=E2[id], E121a=R[2])));
	return( e22a / eDiv);
}
solve.S4HtM.E121aP1 = function(R, sort=TRUE, all.sol=FALSE, debug=TRUE) {
	coeff = coeff.S4Ht.E121aP1(R);
	E2 = roots(coeff);
	if(debug) print(E2);
	E22a = E22a.E121aP1(R, E2);
	len  = length(E2);
	#
	sol = lapply(seq(len), function(id) {
		RS = R; RS[2] = E22a[id];
		solve.S4HtM.Ord2Base(RS, E2[id], sort=sort, all.sol=all.sol)
	})
	sol = do.call(rbind, sol);
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}

### Examples:

### Ex 1:
R = c(3,-1,2,1)
sol = solve.S4HtM.E121aP1(R)

test.S4HtMixed.En3(sol, n=1, nE=c(1,2,1))

apply(sol[c(19:22, 34:37), ], 1, e2.f) * 9*4


### Ex 2:
R = c(-1,-3,2,2)
sol = solve.S4HtM.E121aP1(R)

test.S4HtMixed.En3(sol, n=1, nE=c(1,2,1))

apply(sol[c(8,9, 16,17, 28,29, 42,43), ], 1, e2.f) * 20


### Derivation:

pE121a = toPoly.pm("E2a^2 - E22a - 2*E121a - 4*E4")
pE21 = polyE2Ord1(); # E2, E2a
pE22 = polyE2Ord2(); # E2, E22a

pR1 = solve.pm(pE121a, pE21, "E2a");
pR = solve.pm(pR1$Rez, pE22, "E22a");

# pR$Rez$coeff = - pR$Rez$coeff;
pR$Rez = sort.pm(pR$Rez, xn="E2", xn2=c("E4", "E3", "S"))

# TODO:
# Note:
# - 1200 Monomials;
# - E2^12, S^14;
# - correct roots: should be only E2^{1 or 2}?
# print.pm(pR$Rez, lead="E2")
# print.coeff(pR$Rez, "E2")


