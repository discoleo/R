########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Hetero-Symmetric S4: Mixed
### Basic Types



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

# TODO: proper solution (based on S^4);

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
solve.S4Ht = function(R, debug=FALSE) {
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

R2*x^4 - 2*R3*x^3 + 2*(2*R4 - 2*R2^2 - R1*R2)*x^2


