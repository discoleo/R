########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Hetero-Symmetric S4: Mixed
### Basic Types: Derivation
###
### draft v.0.1a


### Derivation of the Formulas


########################

### Helper Functions

# ...

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


### x = E2:
(4*R[4] - R[2])*x^4 - 2*(R[3]^2 + R[1]^2*R[4])*x^3 +
	+ (4*(R[1]*R[3] - 2*R[4])*R[2] - 8*R[1]*R[3]*R[4] + R[1]^2*R[3]^2 + 2*R[2]^2)*x^2 +
	+ (4*R[1]^3*R[3]*R[4] + 4*R[1]*R[3]^3 + 2*(R[1]^2*R[4] + R[3]^2)*R[2])*x +
	- R[1]^4*R[4]^2 - 2*R[1]^3*R[3]^3 - R[2]^3 - R[3]^4 + 4*R[2]^2*R[4] + 2*R[1]^2*R[3]^2*R[4] +
		+ 8*R[1]*R[2]*R[3]*R[4] - 4*R[1]*R[2]^2*R[3] - 5*R[1]^2*R[2]*R[3]^2;

whichHasPower(4, id=2, type=2)
whichHasPower(R <- c(1,1,-2,-3), id=2, type=2)
polyR(R)
which.sq(DIFF(R), sq=2)


### Hack Coefficients
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

