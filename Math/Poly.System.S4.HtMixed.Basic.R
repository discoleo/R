########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Hetero-Symmetric S4: Mixed
### Basic Types
###
### draft v.0.1c



####################

### Helper Functions

source("Polynomials.Helper.R")


### Other

test.S4HtMixed = function(sol, n=2, R = NULL) {
	x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];
	err1 = x1^n + x2^n + x3^n + x4^n;
	err2 = x1*x2 + x2*x3 + x3*x4 + x4*x1; # Ht!
	err3 = x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4;
	err4 = x1*x2*x3*x4;
	err = rbind(err1, err2, err3, err4);
	if( ! is.null(R)) {
		for(id in 1:4) err[id,] = err[id,] - R[id];
	}
	err = round0(err);
	return(err);
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
### Order 2 ###
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
R2*S^6 - 3*E3*S^5 + (9*E4 - 6*R2^2)*S^4 + (15*E3*R2 - 2*R1*R2)*S^3 +
	- 36*E4*R2*S^2 + 3*R1*E3*S^2 + 9*R2^3*S^2 - 18*E3*R2^2*S + 6*R1*R2^2*S + 9*E3^2*R2 - 6*R1*E3*R2 + R1^2*R2


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

