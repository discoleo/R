########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Hetero-Symmetric S4: Mixed
### Basic Types
###
### draft v.0.2e


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

solve.vandermond = function(A, R) {
	if(is.matrix(A)) {
		m = A;
	} else {
		len = length(A);
		m = matrix(1, nrow=len, ncol=len);
		for(nc in seq(2, len)) {
			m[,nc] = A^(nc - 1);
		}
	}
	solve(m, R);
}
applyChoose = function(x, seq.lst, n, FUN=prod, print.id=NULL) {
	r = sapply(seq(length(n)), function(id) {
		len = length(seq.lst[[id]]);
		m = combn(seq(len), n[[id]]);
		m = seq.lst[[id]][m];
		m = matrix(m, nrow=n[[id]]);
		# print relevant values: but not trivial;
		if(! is.null(print.id) && length(print.id) > 1) print(x[m[ , print.id[[id]] ]]);
		apply(m, 2, function(m) FUN(x[m]));
	})
	r = expand.grid(r);
	if(! is.null(print.id) && length(print.id) == 1) print(r[print.id, ]);
	r = apply(r, 1, FUN);
	return(r);
}
# "Extract" a simplified polynomial:
simple = function(p, vals, sort.by=NULL) {
	p = replace.pm(p, vals);
	p = drop.pm(p);
	xgcd = gcd.vpm(p);
	if(xgcd > 1) {
		p$coeff = p$coeff / xgcd;
		print(paste0("GCD = ", xgcd));
	}
	print(paste0("Monomials = ", nrow(p)));
	if( ! is.null(sort.by)) {
		p = sort.pm(p, xn=sort.by);
	}
	return(p)
}
cmp.pm = function(p1, p2, by=NULL, all=TRUE) {
	#
	if(is.null(by)) {
		by = intersect(names(p1), names(p2));
		idc = match("coeff", by);
		if(is.na(idc)) stop("Not a polynomial!");
		by = by[ - idc];
	}
	pR = merge(p1, p2, by=by, all=all);
	nms = names(pR);
	names(pR)[grepl("^coeff\\.x", nms)] = "coeff";
	xn2 = if(length(by) < 2) NULL else by[2];
	pR = sort.pm(pR, by[1], xn2=xn2);
	return(pR);
}
cmp = function(p1, p2, by=c("E313a", "E3", "E4")) {
	cmp.pm(p1, p2, by=by);
}

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
	+ E3^5 - E3^3*E212a*S + 4*E4*E3^2*(E4*S + E212a) +
	- (E4*S + E212a)^3 # = 0

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
	S = R[1];
	return(solve.S4HtM.E212Base(R, S, E2, sort=sort, all.sol=all.sol))
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


###############
###############

###############
### Order 2 ###
###############

x1^2 + x2^2 + x3^2 + x4^2 - R1 # = 0
x1^2*x2*x3^2 + x2^2*x3*x4^2 + x3^2*x4*x1^2 + x4^2*x1*x2^2 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

### Eq S:
E4^2*E3*S^4 + 2*E4*(E212a*E3 - E4^2)*S^3 - (4*E4*E3^3 + R1*E4^2*E3 - E212a^2*E3 + 6*E212a*E4^2)*S^2 +
	- 2*(E212a*E3^3 - 4*E4^2*E3^2 + R1*E212a*E4*E3 + 3*E212a^2*E4)*S +
	+ 2*E3^5 + 4*R1*E4*E3^3 + 8*E212a*E4*E3^2 - R1*E212a^2*E3 - 2*E212a^3 # = 0


### Solver:
coeff.S4HtM.E212P2 = function(R) {
	# coefficients for S;
	R1 = R[1]; E212a = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(E4^2*E3, 2*E4*(E212a*E3 - E4^2),
		- (4*E4*E3^3 + R1*E4^2*E3 - E212a^2*E3 + 6*E212a*E4^2),
		- 2*(E212a*E3^3 - 4*E4^2*E3^2 + R1*E212a*E4*E3 + 3*E212a^2*E4),
		2*E3^5 + 4*R1*E4*E3^3 + 8*E212a*E4*E3^2 - R1*E212a^2*E3 - 2*E212a^3);
	return(coeff);
}
solve.S4HtM.E212P2 = function(R, sort=TRUE, all.sol=FALSE, debug=TRUE) {
	coeff = coeff.S4HtM.E212P2(R);
	S = roots(coeff);
	if(debug) print(S);
	#
	len = length(S);
	E2 = (S^2 - R[1]) / 2;
	sol = lapply(seq(len), function(id) {
		solve.S4HtM.E212Base(R, S[id], E2[id], sort=sort, all.sol=all.sol)
	})
	sol = do.call(rbind, sol);
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol)
}

### Examples:

### Ex 1:
R = c(3,-1,2,1)
sol = solve.S4HtM.E212P2(R)

test.S4HtMixed.En3(sol, n=2, nE=c(2,1,2))


### Ex 2:
R = c(5,2,3,-1)
sol = solve.S4HtM.E212P2(R)

test.S4HtMixed.En3(sol, n=2, nE=c(2,1,2))

### [filter: all 32 solutions]
round0(poly.calc(sol[seq(1, 32, by=2), 1]) * 9)


### Derivation:

p1 = toPoly.pm("S^2 - 2*E2 - R1");
p2 = polyE2_E212P1();
pR = solve.pm(p1, p2, "E2")
pR = pR$Rez;
pR = sort.pm(pR, c("S", "E3", "E4"), sort.coeff=10:12)
print.pm(pR, lead="S")


###############
###############

###############
### Order 3 ###
###############

n = 3
x1^n + x2^n + x3^n + x4^n - R1 # = 0
x1^2*x2*x3^2 + x2^2*x3*x4^2 + x3^2*x4*x1^2 + x4^2*x1*x2^2 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

### Eq S: P[5]
# - see coeffs;


### Solver:
coeff.S4HtM.E212P3 = function(R) {
	# coefficients for S;
	R1 = R[1]; E212a = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(E3*E4^2,  2*E3*E212a*E4 - 3*E4^3,
		- 4*E3^3*E4 + E3*E212a^2 - 9*E212a*E4^2,
		- 3*E3^3*E212a + 15*E3^2*E4^2 - R1*E3*E4^2 - 9*E212a^2*E4,
		3*E3^5 + 18*E3^2*E212a*E4 - 2*R1*E3*E212a*E4 - 3*E212a^3,
		- 12*E3^4*E4 + 4*R1*E3^3*E4 + 3*E3^2*E212a^2 - R1*E3*E212a^2);
	return(coeff);
}
solve.S4HtM.E212P3 = function(R, sort=TRUE, all.sol=FALSE, debug=TRUE, tol=1E-8) {
	coeff = coeff.S4HtM.E212P3(R);
	S = roots(coeff);
	if(debug) print(S);
	#
	E3 = R[3];
	isZero = (round0(S, tol=tol) == 0);
	isAnyZero = any(isZero);
	if(isAnyZero) {
		warning("Some values S = 0");
		S = S[ ! isZero];
	}
	len = length(S);
	E2 = (S^3 + 3*E3 - R[1]) / (3*S);
	sol = lapply(seq(len), function(id) {
		solve.S4HtM.E212Base(R, S[id], E2[id], sort=sort, all.sol=all.sol)
	})
	sol = do.call(rbind, sol);
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol)
}

### Examples:

### Ex 1:
R = c(3,-1,2,1)
sol = solve.S4HtM.E212P3(R)

test.S4HtMixed.En3(sol, n=3, nE=c(2,1,2))


### Ex 2:
R = c(5,2,3,-1)
sol = solve.S4HtM.E212P3(R)

test.S4HtMixed.En3(sol, n=3, nE=c(2,1,2))


### Ex 3:
# E3 = 0 => S^4
# TODO:
R = c(2,3,0,-4)
sol = solve.S4HtM.E212P3(R)

test.S4HtMixed.En3(sol, n=3, nE=c(2,1,2))


### Ex 4:
# S == 0
# E4 = E212a^2 / (4*E3^2);
R = c(3,4,-1, NA)
R[4] = R[2]^2 / (4*R[3]^2);
# R[4] = R[2]^2 / (4*R[3]^2) + 1E-5;
sol = solve.S4HtM.E212P3(R)

test.S4HtMixed.En3(sol, n=3, nE=c(2,1,2))
# isZ = apply(sol, 1, function(x) any(abs(x) > 10))


### Derivation:

p1 = toPoly.pm("S^3 - 3*E2*S + 3*E3 - R1");
p2 = polyE2_E212P1();
pR = solve.pm(p1, p2, "E2")
pR = pR$Rez;
pR = sort.pm(pR, c("S", "E3", "E4"), sort.coeff=10:12)
print.pm(pR, lead="S")
print.coeff(pR, "S")


### Case: S = 0
x13 = roots(c(1, - E212a/E3, E4));
x24 = E4 / x13;
# xs = E3 / (x24 - x13);
# however: (x24 - x13) == 0!
# => xs = Inf! => NO solution;


########################
########################
########################

####################
### Type: E313a  ###
####################

###############
### Order 1 ###
###############

x1 + x2 + x3 + x4 - R1 # = 0
x1^3*x2*x3^3 + x2^3*x3*x4^3 + x3^3*x4*x1^3 + x4^3*x1*x2^3 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:


### Solver:
# TODO!
# - requires also the polynomial: pX3;
solve.S4HtM.E313P1 = function(R, sort=FALSE, all.sol=FALSE, debug=TRUE) {
	coeff = coeff.S4HtM.E313P1(R);
	E2 = roots(coeff);
	if(debug) print(E2);
	#
	len = length(E2);
	S = R[1];
	x1 = sapply(seq(len), function(id) {
		roots(c(1, -S, E2[id], -R[3], R[4]));
	});
	E2 = rep(E2, each=4);
	x1 = as.vector(x1); len = length(x1);
	x3 = sapply(seq(len), function(id) {
			lst = list(x1=x1[id], E2=E2[id], S=S, E313a=R[2], E4=R[4]);
			eval.pm(pX3$x0, lst) / eval.pm(pX3$div, lst)
		});
	xs = S - x1 - x3; x24 = R[4] / (x1*x3);
	xd = sqrt(xs^2 - 4*x24 + 0i);
	x2 = (xs + xd)/2; x4 = (xs - xd)/2;
	#
	sol = cbind(x1=x1, x2=x2, x3=x3, x4=x4);
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol)
}
# Note: only some special cases;
coeff.S4HtM.E313P1 = function(R) {
	S = R[1]; E313a = R[2]; E3 = R[3]; E4 = R[4];
	b2 = 4*S^2*E4^3 - E313a^2 - 9*E3^2*E4^2 + 6*E313a*E3*E4;
	coeffLead = c(b2 * E3);
	# Special cases:
	coeff = c(coeffLead, - 8941, 14400); # R = c(3,-1,2,3)
	# coeff = c(coeffLead, 282, -1512); # R = c(-2,3,1,2)
	# coeff = c(coeffLead, -21942, -73392); # R = c(-1,3,5,-1)
	return(coeff);
}
# Note: requires the computed polynomial: pR;
coeff.S4HtM.E313P1 = function(R) {
	# coefficients for E2;
	S = R[1]; E313a = R[2]; E3 = R[3]; E4 = R[4];
	vals = list(S=S, E313a=E313a, E3=E3, E4=E4);
	p = replace.pm(pR, vals);
	coeff = coef.pm(p, "E2");
	return(coeff);
}
coeffLead.S4HtM.E313 = function(R) {
	S = R[1]; E313a = R[2]; E3 = R[3]; E4 = R[4];
	b2 = 4*S^2*E4^3 - E313a^2 - 9*E3^2*E4^2 + 6*E313a*E3*E4;
	return(b2);
}

### Examples:

### Ex 1:
R = c(3,-1,2,3)
sol = solve.S4HtM.E313P1(R)

sol = sol[c(15, 40), ];
test.S4HtMixed.En3(sol, n=1, nE=c(3,1,3))
E2  = Re(apply(sol, 1, e2.f))
poly.calc(E2) * 611 * 2 # pD2
# Errors due to inaccuracies;
# 14400 - 8941*E2 + 1222*E2^2


### Ex 2:
R = c(-2,3,1,2)
sol = solve.S4HtM.E313P1(R)

sol = sol[c(14, 20), ];
test.S4HtMixed.En3(sol, n=1, nE=c(3,1,3))
E2  = apply(sol, 1, e2.f)
poly.calc(E2) * 119 # pD2
# -1512 + 282*E2 + 119*E2^2


### Ex 3:
R = c(-1,3,5,-1)
sol = solve.S4HtM.E313P1(R)
test.S4HtMixed.En3(sol, n=1, nE=c(3,1,3))

sol2 = sol[c(1,6), ];
E2  = apply(sol2, 1, e2.f)
poly.calc(E2) * coeffLead.S4HtM.E313(R) * R[3];
# -73392 - 21942*E2 - 1640*E2^2
E3 = 5; # =>
(E3^6 + 2*E3^5 + 4*E3^4 - 2*E3^3 + 8*E3^2 + 24*E3 + 19)*(E3 - 2)
E3^7 - 10*E3^4 + 12*E3^3 + 8*E3^2 - 29*E3 - 38

#
R = c(-1,3, E3, -2)
E3^7 - 31*E3^4 + 78*E3^3 - 26*E3^2 + 22*E3 - 107
#
R = c(-1,3, E3, 2)
E3^7 - 31*E3^4 - 78*E3^3 - 10*E3^2 + 202*E3 + 85
#
R = c(-1,3, E3, 3)
E3^7 - 66*E3^4 - 252*E3^3 - 36*E3^2 + 567*E3 + 378
#
R = c(-1,3, E3, -3)
E3^7 - 66*E3^4 + 252*E3^3 - 90*E3^2 + 297*E3 - 270
#
E3^7 - (7*E4^2 + 3)*E3^4 - (9*E4^3 + 3*E4)*E3^3 +
	+ (E4^3 - 3*E313a*E4^2 + 18)*E3^2 +
	+ (4*E4^4 + 12*E4^2 + 5*E4*E313a^2)*E3 +
	+ (E4^4 + 12*E4^3 - E313a^3)
solve.vandermond(c(-1,-2,3,-3), c(12,78,-252,252))
solve.vandermond(c(-1,-2,3,-3), c(8,-26,-36,-90) + 3*(3)*c(-1,-2,3,-3)^2)
solve.vandermond(c(-1,-2,2,3,-3), c(-29,22,202,567,297) - 5*(3^2)*c(-1,-2,2,3,-3))
solve.vandermond(c(-1,-2,2,3,-3), c(-38,-107,85,378,-270) + (3^3))
# TODO: real powers of S!
E3^7 + (7*E4^2*S^1 - E313a*S^2)*E3^4 - (9*E4^3 - E313a*E4*S^1)*E3^3 +
	- (E4^3*S^1 + 3*E313a*E4^2 + 2*E313a^2*S^1)*E3^2 +
	+ (4*E4^4*S^2 - 4*E4^2*E313a*S^1 + 5*E4*E313a^2)*E3 +
	- (E4^4*S^1 - 4*E4^3*E313a*S^2 + E313a^3)


### Derivation:

p1 = polyE2a();
p2 = polyE2_E212P1();
pE313 = toPoly.pm("E212a*(E2 - E2a) - E4*E3 - E313a");
# pR = solve.pm(p2, pE313, "E212a")
# pR = solve.pm(pR$Rez, p1, "E2a") # 2.800 Monomials;
pR = solve.pm(p1, pE313, "E2a")
pR = solve.pm(pR$Rez, p2, "E212a") # 2.400 Monomials;
pR = pR$Rez;
pR = sort.pm(pR, c("S", "E3", "E4"), sort.coeff=10:12)
str(pR)
top.pm(pR, "E2")
# > 2000 Monomials!
# correct solutions: probably only E2^2;
# print.pm(pR, lead="S")
# print.coeff(pR, "S")
pLead = simplify.spm(top.pm(pR, "E2"), TRUE);
rownames(pLead) = NULL; pLead = drop.pm(pLead);
pLead = div.pm(pLead, toPoly.pm("S^2*E4 - 4*E3^2"), "S")$Rez;
lst = list(S=R[1], E313a=R[2], E3=R[3], E4=R[4]);
eval.pm(pLead, lst)
# roots(coef.pm(replace.pm(pLead, list(S=1, E3=2, E4=-1)), "E313a"))
# toPoly.pm("E4^4*(4*S^2*E4 - 9*E3^2)*(4*S^2*E4 - 25*E3^2)*(S^2*E4 - 4*E3^2)")
div.pm(pLead, toPoly.pm("S^2*E4 - 4*E3^2"), "S")
pD2 = toPoly.pm("4*S^2*E4^3 - E313a^2 - 9*E3^2*E4^2 + 6*E313a*E3*E4")
pD3 = toPoly.pm("4*S^2*E4^3 - E313a^2 - 25*E3^2*E4^2 - 10*E313a*E3*E4")



###
p1 = toPoly.pm("(x1*x3)^6*(S - x1 - x3) + E4^3*(x1 + x3) - E313a*(x1*x3)^3")
p2 = toPoly.pm("x1^2*x3*(S - x1) + x1*x3^2*(S - x1 - x3) + E4 - E2*x1*x3")
p3 = toPoly.pm("(x1*x3)^2*(S - x1 - x3) + E4*(x1 + x3) - E3*x1*x3") # redundant;
pX3 = solve.pm(p2, p1, "x3")
str(pX3)

### Factorize B0 Coeff
# library(gmp)
pT = B0.pm(pR, "E2")
pT$coeff = as.bigz(pT$coeff)
eval.pm(pT, list(S=R[1], E2=1, E313a=R[2], E3=R[3], E4=R[4]))


pT2 = simple(B0.pm(pR, "E2"), c(S=0), sort.by="E313a")


pT3.f = function() {
	pT3 = toPoly.pm("25*E3^6*E4^6 - 10*E3^10*E4^3 + E3^14 + E313a^6 + 14*E3*E4*E313a^5 + 71*E3^2*E4^2*E313a^4 +
			+ 164*E3^3*E4^3*E313a^3 + 191*E3^4*E4^4*E313a^2 + 110*E3^5*E4^5*E313a - 2*E3^7*E313a^3 - 14*E3^8*E4*E313a^2 +
			- 22*E3^9*E4^2*E313a") *
		toPoly.pm("E4^2*(9*E3^5*E4^5 - E3^9*E4^2 + E313a^5 - 3*E3*E4*E313a^4 - 6*E3^2*E4^2*E313a^3 + 10*E3^3*E4^3*E313a^2 +
			+ 21*E3^4*E4^4*E313a - E3^7*E313a^2 - 2*E3^8*E4*E313a)");
}
pT3 = pT3.f(); cmp(pT2, pT3)

pT3 = simple(B0.pm(pR, "E2"), c(S=0, E3=1, E4=-1), sort.by="E313a")
r = roots(pT3$coeff); r;
poly.calc(r[-c(3,4,9, 7,8)])

toPoly.pm("E4^2*(E313a + E3*E4)^2*(E313a^3 - 5*E3*E4*E313a^2 + 3*E3^2*E4^2*E313a - E3^7 + 9*E3^3*E4^3)")

###
pT2 = simple(B0.pm(pR, "E2"), c(S=1, E3=1, E4=1), sort.by="E313a")
div.pm(pT2, toPoly.pm("(E313a - 1)^7"), "E313a")


########################
########################
########################

