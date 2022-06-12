########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S2
### Asymmetric Derived from Symmetric
### Transform: Based on Roots of Unity
###
### draft v.0.1e


#######################

### Helper Functions

source("Polynomials.Helper.R")


#######################

### Base System:
# x^3 + y^3 + b*(x+y) = R1
# x*y = R2

### Extension:
# x*y + b2*(x+y) = R2

### Transform:
# x => x + m*y
# y => x + m^2*y
# where m^3 = 1

### Derived System
2*x^3 + 2*y^3 - 3*x*y*(x + y) + b1*(2*x - y) - R[1] # = 0
x^2 + y^2 - x*y - R[2] # = 0
# Extension:
x^2 + y^2 - x*y + b2*(2*x - y) - R[2] # = 0


### Solver
# Symmetric System:
solve.S2Symm.P3 = function(R, b, all=TRUE, debug=TRUE) {
	# Step 1:
	# S^3 - (3*R2 - b)*S - R1 = 0
	if(length(b) == 1) {
		S = roots(c(1,0, - 3*R[2] + b[1], -R[1]));
		R2 = R[2];
	} else {
		S = roots(c(1, 3*b[2]*R[2], - 3*R[2] + b[1], -R[1]));
		R2 = R[2] - b[2]*S;
	}
	if(debug) print(S);
	# Step 2:
	d = rootn(S^2 - 4*R2, 2);
	x = (S + d)/2;
	y = S - x;
	# Solution:
	sol = cbind(x=x, y=y);
	if(all) sol = rbind(sol, sol[, c(2,1)]);
	return(sol);
}
# Asymmetric System:
solve.S2As.P3 = function(R, b, debug=TRUE) {
	# Step 1 & 2:
	sol = solve.S2Symm.P3(R, b=b, debug=debug, all=TRUE);
	d = sol[,1] - sol[,2];
	S = sol[,1] + sol[,2];
	# Step 3:
	m = unity(3, all=FALSE);
	y = d / (m - m^2);
	x = (S + y) / 2;
	cbind(x=x, y=y);
}
test.S2As.P3 = function(sol, R=NULL, b) {
	if(length(b) < 2) b = c(b, 0);
	err1 = 2*x^3 + 2*y^3 - 3*x*y*(x + y) + b[1]*(2*x - y);
	err2 = x^2 + y^2 - x*y + b[2]*(2*x - y);
	err = rbind(err1, err2);
	if( ! is.null(R)) {
		err = err - R;
	}
	err = round0(err);
	return(err);
}

### Examples:
R = c(1,2)
b = 3
sol = solve.S2As.P3(R, b);
x = sol[,1]; y = sol[,2];

### Test
test.S2As.P3(sol, R, b)

### Classic Poly:
round0.p(poly.calc(x) * 27)


### Ex 2:
R = c(1,-2)
b = 6*R[2]
sol = solve.S2As.P3(R, b);
x = sol[,1]; y = sol[,2];

### Test
test.S2As.P3(sol, R, b)

### Classic Poly:
round0.p(poly.calc(x) * 27)


### Ex 3:
R = c(5, 1)
b = 6*R[2]
sol = solve.S2As.P3(R, b);
x = sol[,1]; y = sol[,2];

### Test
test.S2As.P3(sol, R, b)

### Classic Poly:
round0.p(poly.calc(x) * 3)
1 - 10*x + 5*x^2 + 3*x^6


### Ex 4:
R = c(5, 1)
b = c(6*R[2], -1)
sol = solve.S2As.P3(R, b);
x = sol[,1]; y = sol[,2];

### Test
test.S2As.P3(sol, R, b)

### Classic Poly:
round0.p(poly.calc(x) * 27)


### Derivation:

### Eq 1:
(x+m*y)^3 + (x+m^2*y)^3 + b*(2*x + (m+m^2)*y) - R1 # = 0
2*x^3 + 2*y^3 + 3*x*y*(m*x + m^2*y + m^2*x + m*y) + b*(2*x - y) - R1 # = 0
2*x^3 + 2*y^3 - 3*x*y*(x + y) + b*(2*x - y) - R1 # = 0

### Eq 2:
(x + m*y)*(x + m^2*y) - R2 # = 0
x^2 + y^2 - x*y - R2 # = 0

### Classic Poly
R1 = R[1]; R2 = R[2];
27*x^6 - 9*(6*R2 - b)*x^4 + 3*(9*R2^2 - 5*b*R2 + b^2)*x^2 - 3*R1*b*x +
	+ R1^2 - 4*R2^3 + 4*b*R2^2 - b^2*R2


##############

### Variant 1:

### Derived:
# x => x + i*y
# y => x - i*y

### Derived System
x^3 - 3*x*y^2 + b*x - R1/2 # = 0
x^2 + y^2 - R2 # = 0
# [trivial]

### Derivation:

### Eq 1:
(x + 1i*y)^3 + (x - 1i*y)^3 + 2*b*x - R1 # = 0
2*x^3 - 6*x*y^2 + 2*b*x - R1 # = 0

### Eq 2:
x^2 + y^2 - R2 # = 0

### Poly:
4*x^3 + (b - 3*R2)*x - R1/2 # = 0


#######################
#######################

### Base System: P[5]
# x^5 + y^5 + b*(x+y) = R1
# x*y = R2

### Derived:
# x => x + m*y
# y => x + m^2*y
# where m^3 = 1

### Derived System
2*x^5 - y^5 - 5*x*y*(x^3 + y^3) - 10*x^2*y^2*(x - 2*y) + b*(2*x - y) - R1 # = 0
x^2 + y^2 - x*y - R2 # = 0
# =>
2*x^5 - y^5 - 10*x^2*y^2*(x - 2*y) - 5*R2*x*y*(x+y) + b*(2*x - y) - R1 # = 0
x^2 + y^2 - x*y - R2 # = 0

# Note:
# - Base-System is easy to solve for b = 0;


### Solver
solve.S2As.P5 = function(R, b, debug=TRUE) {
	m = unity(3, all=FALSE);
	# Step 1:
	# S^5 - 5*R2*S^3 + (5*R2^2 + b)*S - R1 = 0
	S = roots(c(1,0, - 5*R[2], 0, 5*R[2]^2 + b, -R[1]));
	if(debug) print(S);
	# Step 2:
	d = rootn(S^2 - 4*R[2], 2);
	d = c(d, -d); S = c(S, S);
	# Step 3:
	y = d / (m - m^2);
	x = (S + y) / 2;
	cbind(x=x, y=y);
}
test.S2As.P5 = function(sol, R=NULL, b) {
	err1 = 2*x^5 - y^5 - 10*x^2*y^2*(x - 2*y) +
		- 5*R[2]*x*y*(x+y) + b*(2*x - y);
	err2 = x^2 + y^2 - x*y;
	err = rbind(err1, err2);
	if( ! is.null(R)) {
		err = err - R;
	}
	err = round0(err);
	return(err);
}

### Examples:
R = c(1,2)
b = 3
sol = solve.S2As.P5(R, b);
x = sol[,1]; y = sol[,2];

### Test
test.S2As.P5(sol, R, b)

### Classic Poly:
round0.p(poly.calc(x) * 3^5)


### Derivation:

### Eq 1:
(x+m*y)^5 + (x+m^2*y)^5 + b*(2*x + (m+m^2)*y) - R1 # = 0
2*x^5 + (m + m^2)*y^5 + 5*x*y*(m*x^3 + m*y^3 + m^2*x^3 + m^2*y^3) +
	+ 10*x^2*y^2*(m^2*x + y + m*x + y) + b*(2*x - y) - R1 # = 0
2*x^5 - y^5 - 5*x*y*(x^3 + y^3) - 10*x^2*y^2*(x - 2*y) + b*(2*x - y) - R1 # = 0

### Eq 2:
(x + m*y)*(x + m^2*y) - R2 # = 0
x^2 + y^2 - x*y - R2 # = 0


########################
########################

########################
### Hetero-Symmetric ###
########################

### Base System: P[3]
# x^3 + b*y = R
# y^3 + b*x = R

### Derived:
# x => x + m*y
# y => x + m^2*y
# where m^3 = 1

### Derived System
2*x^3 + 2*y^3 - 3*x*y*(x+y) + b*(2*x - y) - 2*R # = 0
3*x*(x - y) - b # = 0


### Solver

# - for exact solution of Ht-System, see:
#   Poly.System.Hetero.Symmetric.R;
solve.Ht.P3 = function(R, b, debug=TRUE) {
	if(length(b) == 1) {
		coeff = c(1, 0, -2*b[1], R);
	} else {
		coeff = c(1, - b[2], -2*b[1], R + b1[1]*b[2]);
	}
	S = roots(coeff);
	if(debug) print(S);
	xy = S^2 - b[1];
	d  = sqrt(S^2 - 4*xy + 0i);
	x = (S + d)/2;
	y = (S - d)/2;
	sol = cbind(x, y);
	sol = rbind(sol, sol[,2:1]);
	return(sol);
}
solve.DerHt.P3 = function(R, b, debug=TRUE) {
	sol = solve.Ht.P3(R, b, debug=debug);
	m = unity(3, all=FALSE);
	y = (sol[,1] - sol[,2]) / (m - m^2);
	x = (sol[,1] + sol[,2] + y) / 2;
	return(cbind(x=x, y=y));
}
solve.DerHtTrig.P3 = function(R, b, debug=TRUE) {
	sol = solve.Ht.P3(R, b, debug=debug);
	r = 2*cos(1:2 * 2*pi/5);
	y = (sol[,1] - sol[,2]) / (r[1] - r[2]);
	x = sol[,1] - r[1]*y;
	return(cbind(x=x, y=y));
}
test.DerHt.P3 = function(sol, R = NULL, b) {
	x = sol[,1]; y = sol[,2];
	if(length(b) > 1) stop("Not yet implemented!");
	err1 = 2*x^3 + 2*y^3 - 3*x*y*(x+y) + b*(2*x - y);
	err2 = 3*x*(x - y) - b[1];
	err = rbind(err1, err2);
	if( ! is.null(R)) {
		err[1,] = err[1,] - 2*R;
	}
	err = round0(err);
	return(err);
}
test.DerHtTrig.P3 = function(sol, R = NULL, b) {
	x = sol[,1]; y = sol[,2];
	b1 = b[1];
	b2 = if(length(b) < 2) 0 else b[2];
	xy = x*y;
	err1 = 2*x^3 - 4*y^3 - 3*xy*(x - 3*y) + 2*b2*(x^2 - xy - y^2) + 2*b1*x - b1*y;
	err2 = 3*x^2 - 3*xy + 2*y^2 - b1;
	err = rbind(err1, err2);
	if( ! is.null(R)) {
		err[1,] = err[1,] - 2*R;
	}
	err = round0(err);
	return(err);
}

### Examples:
R = 2
b = -1
sol = solve.DerHt.P3(R, b);

test.DerHt.P3(sol, R, b)


### Derivation:

### Sum =>
(x + m*y)^3 + (x + m^2*y)^3 + b*(2*x + m*y + m^2*y) - 2*R # = 0
2*x^3 + 2*y^3 - 3*x*y*(x+y) + b*(2*x - y) - 2*R # = 0

### Diff =>
3*x*y*(m*x + m^2*y - m^2*x - m*y) + b*(m^2-m)*y # = 0
# y != 0 =>
3*x*(m - m^2)*(x - y) + b*(m^2-m) # = 0
3*x*(x - y) - b # = 0

###
p1 = toPoly.pm("2*x^3 + 2*y^3 - 3*x*y*(x+y) + b*(2*x - y) - 2*R")
p2 = toPoly.pm("3*x*(x - y) - b")
pR = solve.pm(p1, p2, "y")
pR = pR$Rez;
pR$coeff = - pR$coeff;
pR = sort.pm(pR, "x")
print.pm(pR, sort=FALSE, lead="x")

27*x^6 - 27*b*x^4 + 27*R*x^3 - 9*b^2*x^2 + b^3
# only "resembles" (3*x^2 +/- b)^3 + 27*R*x^3
# b => 3*b =>
x^6 - 3*b*x^4 + R*x^3 - 3*b^2*x^2 + b^3
(x^2 + b)*(x^4 - 4*b*x^2 + b^2) + R*x^3

### y: [trivial]
27*y^6 - 36*b^2*y^2 + 16*b^3 - 27*R^2
# b => 3*b =>
y^6 - 12*b^2*y^2 + 16*b^3 - R^2


###################

### Extension:
### Base System: P[3]
# x^3 + b2*x*y + b1*y = R
# y^3 + b2*x*y + b1*x = R

### Derived:
# x => x + m*y
# y => x + m^2*y
# where m^3 = 1

### Derived System
2*x^3 + 2*y^3 - 3*x*y*(x+y) + 2*b2*(x^2 + y^2 - x*y) + b1*(2*x - y) - 2*R # = 0
3*x*(x - y) - b1 # = 0

# TODO: check;

###
p1 = toPoly.pm("2*x^3 + 2*y^3 - 3*x*y*(x+y) + 2*b2*(x^2 + y^2 - x*y) + b1*(2*x - y) - 2*R")
p2 = toPoly.pm("3*x*(x - y) - b1")
pR = solve.pm(p1, p2, "y")
pR = pR$Rez;
pR$coeff = - pR$coeff;
pR = sort.pm(pR, "x")
print.pm(pR, sort=FALSE, lead="x")

27*x^6 - 27*b2*x^5 - 27*b1*x^4 + 9*b1*b2*x^3 + 27*R*x^3 - 9*b1^2*x^2 - 3*b1^2*b2*x + b1^3

### y:
27*y^6 - 27*b2^2*y^4 - 18*b1*b2^2*y^2 - 36*b1^2*y^2 + 54*R*b2*y^2 - 3*b1^2*b2^2 + 16*b1^3 + 18*R*b1*b2 - 27*R^2


#############

### Variants:
### V.1.) Trigonometric

### Derived:
# x => x + r1*y
# y => x + r2*y
# where r[j] = 2*cos(2*j*pi/5);

### System
2*x^3 - 3*x^2*y + 9*x*y^2 - 4*y^3 + 2*b2*x^2 - 2*b2*x*y - 2*b2*y^2 + 2*b1*x - b1*y - 2*R # = 0
3*x^2 - 3*x*y + 2*y^2 - b1 # = 0

### Solver
# see above (base-system);

### Examples

R = 5
b = c(-10, 15)
sol = solve.DerHtTrig.P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
test.DerHtTrig.P3(sol, R=R, b=b)

round0.p(poly.calc(x)) * 5


### Derivation:

r = 2*cos(1:2 * 2*pi/5)
# sum(r) = -1
# prod(r) = -1

### Sum =>
S^3 - 3*xy*S + 2*b2*xy + b1*S - 2*R # = 0
(2*x - y)^3 - 3*(x^2 - y^2 - x*y)*(2*x - y) + 2*b2*(x^2 - y^2 - x*y) + b1*(2*x - y) - 2*R # = 0
#
pS = toPoly.pm("2*x - y")
pXY = toPoly.pm("x^2 - y^2 - x*y")
p1 = toPoly.pm("pS()^3 - 3*pXY()*pS() + 2*b2*pXY() + b1*pS() - 2*R")
p1 = sort.pm(p1, c("x", "y"))
print.pm(p1, sort=F, lead=c("x", "y"))

### Diff =>
S^2 - xy - b1 # = 0
#
p2 = toPoly.pm("pS()^2 - pXY() - b1")
p2 = sort.pm(p2, c("x", "y"))
print.pm(p2, sort=F, lead=c("x", "y"))


### Classic Poly
pR = solve.pm(p1, p2, "y")
pR = pR$Rez
pR = sort.pm(pR, "x", xn2=c("R"))
print.pm(pR, sort=F, lead=c("x"))

125*x^6 - 125*b2*x^5 + 25*(2*b2^2 - 5*b1)*x^4 - 25*(R - 5*b1*b2)*x^3 +
	- 5*(2*R*b2 + 6*b1*b2^2 - 13*b1^2)*x^2 +
	- 5*b1*(2*R + 7*b1*b2)*x +
	+ 8*R^2 + 8*R*b1*b2 + 2*b1^2*b2^2 - 9*b1^3



