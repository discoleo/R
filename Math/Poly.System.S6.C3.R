########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S6: C3-Hetero-Symmetric
### Basic Types
###
### draft v.0.1a


####################

### Helper Functions


source("Polynomials.Helper.R")


#####################
#####################

##############
### Theory ###
##############

### Basic Systems:

### Eqs 1 & 2:
# x1^n1 + x2^n1 + x3^n1 = R1
# y1^n2 + y2^n2 + y3^n2 = R2

### A.) Type 2-Ht:
# x1*y1 + x2*y2 + x3*y3 = R3
# x1*y2 + x2*y3 + x3*y1 = R4
# x1*x2*x3 = R5
# y1*y2*y3 = R6

### Symmetries:
# - if {(x1, x2, x3), (y1, y2, y3)} is a solution,
#   so are also all concordant C3-cyclic permutations:
#   e.g. {(x2, x1, x3), (y2, y1, y3)};


### B.) Type 1-Ht & E2:
# - has additional symmetries;
# x1*y1 + x2*y2 + x3*y3 = R3
# E2x + E2y = R4
# x1*x2*x3 = R5
# y1*y2*y3 = R6
### Variant Eq 4:
# E2x * E2y = R4


### C.) Type 1-Ht & 2 x E2:
# - has additional symmetries;
# x1*y1 + x2*y2 + x3*y3 = R3
# E2x = R4
# E2y = R5
# x1*x2*x3*y1*y2*y3 = R6


###################
###################

#################
### Type 2-Ht ###
#################

### Basic System: Eqs 3 & 4
# x1*y1 + x2*y2 + x3*y3 = R3
# x1*y2 + x2*y3 + x3*y1 = R4


### Solution:

### Sum:
A1 + B1 + C1 - sx*sy # = 0
# C1 = sx*sy - A1 - B1;

### E2:
E2ABC + 3*E2x*E2y - E2y*sx^2 - E2x*sy^2 # = 0

### E3:
A1*B1*C1 - py*(sx^3 - 3*E2x*sx) - px*(sy^3 - 3*E2y*sy) - 18*px*py +
	- 2*A2a*B2a + (E2x*sx - 3*px)*B2a + (E2y*sy - 3*py)*A2a +
	- E2x*sx*E2y*sy + 3*E2x*py*sx + 3*E2y*px*sy # = 0
# with:
A2a^2 - (E2x*sx - 3*px)*A2a - (6*E2x*px*sx - E2x^3 - 9*px^2 - px*sx^3) # = 0
B2a^2 - (E2y*sy - 3*py)*B2a - (6*E2y*py*sy - E2y^3 - 9*py^2 - py*sy^3) # = 0


### Solver:

solver.S6C3.P1Ht2 = function(R, debug=TRUE, all=FALSE) {
	coeff = coeff.S6C3.P1Ht2(R);
	E2x = roots(coeff);
	if(debug) print(E2x);
	#
	sx = R[1]; sy = R[2]; px = R[5]; py = R[6];
	A1 = R[3]; B1 = R[4];
	C1 = sx*sy - A1 - B1;
	E2ABC = (A1 + B1)*C1 + A1*B1;
	E2y = (sy^2*E2x - E2ABC) / (3*E2x - sx^2);
	#
	len = length(E2x);
	x1 = sapply(seq(len), function(id) {
		roots(c(1, -sx, E2x[id], -px));
	})
	print(x1)
	# robust:
	E2x = rep(E2x, each=3); E2y = rep(E2y, each=3);
	x23 = sx - x1; px23 = E2x - x1*x23;
	# TODO
}
coeff.S6C3.P1Ht2 = function(R) {
	sx = R[1]; sy = R[2]; px = R[5]; py = R[6];
	A1 = R[3]; B1 = R[4];
	C1 = sx*sy - A1 - B1;
	E2ABC = (A1 + B1)*C1 + A1*B1;
	coeff = c(sy^6 - 54*sy^3*py + 729*py^2,
		- 6*sy^4*E2ABC + 54*sy^3*sx^2*py + 162*sy*E2ABC*py - 1458*sx^2*py^2,
		sy^4*sx^2*E2ABC + 9*sy^2*E2ABC^2 - 18*sy^3*sx^4*py - 189*sy*sx^2*E2ABC*py + 1215*sx^4*py^2 +
			+ 9*sy^3*sx*C1*B1*A1 - 243*sx*py*C1*B1*A1,
		- 2*sy^2*sx^2*E2ABC^2 - 4*E2ABC^3 - 2*sy^6*sx^3*px + 9*sy^4*sx*E2ABC*px + 2*sy^3*sx^6*py +
			+ 81*sy*sx^4*E2ABC*py + 54*sy^3*sx^3*px*py - 243*sy*sx*E2ABC*px*py - 540*sx^6*py^2 +
			- 6*sy^3*sx^3*C1*B1*A1 - 9*sy*sx*E2ABC*C1*B1*A1 - 27*sy^3*px*C1*B1*A1 + 297*sx^3*py*C1*B1*A1 +
			+ 729*px*py*C1*B1*A1 - 27*C1^2*B1^2*A1^2,
		sx^2*E2ABC^3 + 6*sy^4*sx^3*E2ABC*px - 27*sy^2*sx*E2ABC^2*px - 15*sy*sx^6*E2ABC*py - 54*sy^3*sx^5*px*py +
			+ 243*sy*sx^3*E2ABC*px*py + 135*sx^8*py^2 + sy^3*sx^5*C1*B1*A1 + 6*sy*sx^3*E2ABC*C1*B1*A1 +
			+ 81*sy*E2ABC*px*C1*B1*A1 - 135*sx^5*py*C1*B1*A1 - 729*sx^2*px*py*C1*B1*A1 + 27*sx^2*C1^2*B1^2*A1^2,
		- sy^4*sx^5*E2ABC*px + 18*sx*E2ABC^3*px + sy*sx^8*E2ABC*py + 18*sy^3*sx^7*px*py - 81*sy*sx^5*E2ABC*px*py +
			- 18*sx^10*py^2 - sy*sx^5*E2ABC*C1*B1*A1 + 9*sy^3*sx^4*px*C1*B1*A1 - 54*sy*sx^2*E2ABC*px*C1*B1*A1 +
			+ 27*sx^7*py*C1*B1*A1 + 243*sx^4*px*py*C1*B1*A1 - 9*sx^4*C1^2*B1^2*A1^2,
		sy^2*sx^5*E2ABC^2*px - 4*sx^3*E2ABC^3*px + sy^6*sx^6*px^2 - 9*sy^4*sx^4*E2ABC*px^2 +
			+ 27*sy^2*sx^2*E2ABC^2*px^2 - 27*E2ABC^3*px^2 - 2*sy^3*sx^9*px*py + 9*sy*sx^7*E2ABC*px*py +
			+ sx^12*py^2 - 2*sy^3*sx^6*px*C1*B1*A1 + 9*sy*sx^4*E2ABC*px*C1*B1*A1 - 2*sx^9*py*C1*B1*A1 +
			- 27*sx^6*px*py*C1*B1*A1 + sx^6*C1^2*B1^2*A1^2
	);
	return(coeff);
}

### Examples

R = c(2,3,-1,4,5,6)
sol = solver.S6C3.P1Ht2(R)


### Debug

R = c(2,3,-1,4,5,6)
x1 = -0.0751764371 - 1.2207329994i;
x2 = -0.3636407258 + 1.6144160902i;
x3 =  2.4388171629 - 0.3936830909i;
y1 =  0.1828061717 + 1.2261038046i;
y2 =  2.8373554636 + 0.4578248858i;
y3 = -0.0201616352 - 1.6839286904i;


sx = x1 + x2 + x3; sy = y1 + y2 + y3;
px = x1 * x2 * x3; py = y1 * y2 * y3;
E2x = (x1 + x2)*x3 + x1*x2;
E2y = (y1 + y2)*y3 + y1*y2;
#
A1 = x1*y1 + x2*y2 + x3*y3;
B1 = x1*y2 + x2*y3 + x3*y1;
C1 = x1*y3 + x2*y1 + x3*y2;
E2ABC = (A1 + B1)*C1 + A1*B1;
# [intermediary]
A2a = x1^2*x2 + x2^2*x3 + x3^2*x1;
A2b = x1*x2^2 + x2*x3^2 + x3*x1^2;
B2a = y1^2*y2 + y2^2*y3 + y3^2*y1;
B2b = y1*y2^2 + y2*y3^2 + y3*y1^2;


### Test:

x1 + x2 + x3 # = R1
y1 + y2 + y3 # = R2
x1*y1 + x2*y2 + x3*y3 # = R3
x1*y2 + x2*y3 + x3*y1 # = R4
x1*x2*x3 # = R5
y1*y2*y3 # = R6


### Derivation:
p2 = toPoly.pm(...); # see formula for E3 above;
p1 = toPoly.pm("E2ABC + 3*E2x*E2y - E2y*sx^2 - E2x*sy^2");
pA2a = toPoly.pm("A2a^2 - (E2x*sx - 3*px)*A2a - (6*E2x*px*sx - E2x^3 - 9*px^2 - px*sx^3)");
pB2a = toPoly.pm("B2a^2 - (E2y*sy - 3*py)*B2a - (6*E2y*py*sy - E2y^3 - 9*py^2 - py*sy^3)");

# Note: B2a cancels out:
pR = solve.lpm(p2, pA2a, pB2a, xn=c("A2a", "B2a"))
pR = drop.pm(pR[[2]]$Rez)

pR2 = solve.pm(p1, pR, "E2y")
str(pR2)
toCoeff(pR2$Rez, "E2x")

