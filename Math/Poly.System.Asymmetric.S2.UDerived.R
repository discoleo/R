########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S2
### Asymmetric Derived from Symmetric
### Transform: Based on Roots of Unity
###
### draft v.0.2c


#######################

### Helper Functions

source("Polynomials.Helper.R")


#######################

### Theory

### Base-System:
# - is the starting system;
# - various transformations are applied to the base system,
#   e.g. specific substitutions;


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
# Generalized Asymmetric:
solve.S2AsGen.P3 = function(R, b, k, debug=TRUE) {
	# Step 1 & 2:
	sol = solve.S2Symm.P3(R, b=b, debug=debug, all=TRUE);
	# Step 3a:
	r = roots(c(1, -k[1], k[2]));
	# Step 3b:
	y = (sol[,1] - sol[,2]) / (r[1] - r[2]);
	x = sol[,1] - r[1]*y;
	sol = cbind(x=x, y=y);
	return(sol);
}
test.S2As.P3 = function(sol, b, R=NULL) {
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
test.S2AsGen.P3 = function(sol, b, k, R=NULL) {
	if(length(b) < 2) b = c(b, 0);
	k1 = k[1]; k2 = k[2];
	err1 = 2*x^3 + k1*(k1^2 - 3*k2)*y^3 + 3*k1*x^2*y + 3*(k1^2 - 2*k2)*x*y^2 +
		+ b[1]*(2*x + k1*y);
	err2 = x^2 + k2*y^2 + k1*x*y + b[2]*(2*x + k1*y);
	err  = rbind(err1, err2);
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
test.S2As.P3(sol, b=b, R=R)

### Classic Poly:
round0.p(poly.calc(x) * 27)


### Ex 2:
R = c(1,-2)
b = 6*R[2]
sol = solve.S2As.P3(R, b);
x = sol[,1]; y = sol[,2];

### Test
test.S2As.P3(sol, b=b, R=R)

### Classic Poly:
round0.p(poly.calc(x) * 27)


### Ex 3:
R = c(5, 1)
b = 6*R[2]
sol = solve.S2As.P3(R, b);
x = sol[,1]; y = sol[,2];

### Test
test.S2As.P3(sol, b=b, R=R)

### Classic Poly:
round0.p(poly.calc(x) * 3)
1 - 10*x + 5*x^2 + 3*x^6


### Ex 4:
R = c(5, 1)
b = c(6*R[2], -1)
sol = solve.S2As.P3(R, b);
x = sol[,1]; y = sol[,2];

### Test
test.S2As.P3(sol, b=b, R=R)

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

### Transform:
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


##############

### Variant 2:
### Generalized

### Transform:
# x => x + r1*y
# y => x + r2*y
# where: r1 + r2 = k1, r1*r2 = k2;

### Derived System
k1 = k[1]; k2 = k[2]; b1 = b[1]; b2 = b[2];
2*x^3 + k1*(k1^2 - 3*k2)*y^3 + 3*k1*x^2*y + 3*(k1^2 - 2*k2)*x*y^2 +
	+ b1*(2*x + k1*y) - R[1] # = 0
x^2 + k2*y^2 + k1*x*y + b2*(2*x + k1*y) - R[2] # = 0

### Examples:

###
k = c(-2, 2)
# =>
x^3 + 2*y^3 - 3*x^2*y + b[1]*(x - y) - R[1]/2 # = 0
x^2 + 2*y^2 - 2*x*y + 2*b[2]*(x - y) - R[2] # = 0

### Ex 1:
R = c(-2, 1)
b = c(-1, -1)
k = c(-2, 2)
#
sol = solve.S2AsGen.P3(R, b, k=k)
x = sol[,1]; y = sol[,2];

# Test
test.S2AsGen.P3(sol, b=b, k=k, R=R)

### Classic Poly
round0.p(poly.calc(x))
round0.p(poly.calc(y)) # trivial

R1h = R[1] / 2; R2 = R[2];
b1 = b[1]; b2d = 2*b[2];
k1 = -2; k2 = 2; # fixed values / only special Case!
8*x^6 + 12*b2d*x^5 - 3*(8*R2 - b2d^2)*x^4 + 2*(4*R1h - 9*R2*b2d - b1*b2d)*x^3 +
	+ (12*b2d*R1h + 18*R2^2 - 3*b2d^2*R2 - 4*b1*R2 - b1*b2d^2 + 2*b1^2)*x^2 +
	- 2*(6*R1h*R2 + 2*b1*R1h - 3*b2d*R2^2 - b1*b2d*R2)*x +
	+ 4*R1h^2 - 6*b2d*R1h*R2 - b2d^3*R1h + 2*b1*b2d*R1h +
	- 2*R2^3 + 4*b1*R2^2 + b1*b2d^2*R2 - 2*b1^2*R2


### Ex 2:
R = c(-2, 1)
b = c(-1, -1)
k = c(1, 3)
#
sol = solve.S2AsGen.P3(R, b, k=k)
x = sol[,1]; y = sol[,2];

# Test
test.S2AsGen.P3(sol, b=b, k=k, R=R)

round0.p(poly.calc(x) * 11^3)


### Ex 3:
R = c(-2, 1)
b = c(-1, 0)
k = c(1, 3)
#
sol = solve.S2AsGen.P3(R, b, k=k)
x = sol[,1]; y = sol[,2];

# Test
test.S2AsGen.P3(sol, b=b, k=k, R=R)

round0.p(poly.calc(x) * 11^3)


### only for Ex 1:
# k = c(-2, 2)!
p1 = toPoly.pm("x^3 + 2*y^3 - 3*x^2*y + b1*(x - y) - R1h")
p2 = toPoly.pm("x^2 + 2*y^2 - 2*x*y + b2d*(x - y) - R2")
pR = solve.pm(p1, p2, "y")
pR = pR$Rez;
pR = sort.pm(pR, "x", c("R1h", "R2"))
print.pm(pR, sort=F, lead="x")

pR = toPoly.pm("p2()*(x+y) - p1()")
pR = sort.pm(pR, c("x", "y"), c("R1h", "R2"))
print.pm(pR, sort=F, lead=NA)
# alternative Eq:
2*x^2*y + b2d*x^2 - b2d*y^2 - (R2 + b1)*x - (R2 - b1)*y + R1h # = 0


##############

### Variant 3:
### Generalized / with radicals

### Transform:
# x => (s10 + k)*x + (s20 - 2*k)*y
# y => (s10 - k)*x + (s20 + 2*k)*y
# where:
#   k^2 = K,
#   s = c(s10, s20) given parameters,
#   while (s11, s21) = c(1, -2) (parameters omitted);

# TODO


### Derivation:
# x + y => 2*(s10*x + s20*y);
# x*y => (s10^2 - K)*x^2 + (s20^2 - 4*K)*y^2 + 2*(s10*s20 + 2*K)*x*y

### Eq 1:
8*(s10*x + s20*y)^3 +
	- 6*((s10^2 - K)*x^2 + (s20^2 - 4*K)*y^2 + 2*(s10*s20 + 2*K)*x*y)*(s10*x + s20*y) +
	+ 2*b1*(s10*x + s20*y) - R1
# =>
s10*(s10^2 + 3*K)*x^3 + s20*(s20^2 + 12*K)*y^3 +
	+ 3*(s10^2*s20 + s20*K - 4*s10*K)*x^2*y +
	+ 3*(s10*s20^2 - 4*s20*K + 4*s10*K)*x*y^2 + s10*b1*x + s20*b1*y - R1/2 # = 0

### Eq 2:
(s10^2 - K)*x^2 + (s20^2 - 4*K)*y^2 + 2*(s10*s20 + 2*K)*x*y + 2*b2*(s10*x + s20*y) - R2 # = 0


#######################
#######################

### Base System: P[5]
# x^5 + y^5 + b*(x+y) = R1
# x*y = R2

### Transform:
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

########################
### Hetero-Symmetric ###
########################

### Base System: P[3]
# x^3 + b*y = R
# y^3 + b*x = R

### Transform:
# x => x + m*y
# y => x + m^2*y
# where m^3 = 1

### Derived System:
# Sum & Diff:
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


####################

### Variant: Product

# Sum =>
2*x^3 + 2*y^3 - 3*x*y*(x+y) + b*(2*x - y) - 2*R # = 0

# Original =>
(x + m*y)^3 - R + b*(x + m^2*y) # = 0
x^3 + y^3 + b*x - R + 3*m*x*y*(x + m*y) + b*m^2*y # = 0
# Prod
(x^3 + y^3 + b*x - R)^2 - m^2*(3*x*y*(x + m*y) + b*m*y)*(3*m*x*y*(x + m^2*y) + b*y) # = 0
x^6 - 9*x^4*y^2 + 11*x^3*y^3 - 9*x^2*y^4 + y^6 + 2*b*x^4 + 3*b*x^2*y^2 - 4*b*x*y^3 +
	- 2*R*x^3 - 2*R*y^3 + b^2*x^2 - b^2*y^2 - 2*b*R*x + R^2 # = 0
# Reduction =>
12*x^5*y - 33*x^4*y^2 + 18*x^3*y^3 - 33*x^2*y^4 + 12*x*y^5 + 16*b*x^3*y + 18*b*x^2*y^2 - 30*b*x*y^3 + 4*b*y^4 +
	- 12*R*x^2*y - 12*R*x*y^2 + 4*b^2*x*y - 5*b^2*y^2 - 4*b*R*y # = 0
27*x^4*y^2 - 54*x^3*y^3 + 27*x^2*y^4 - 4*b*x^3*y - 12*b*x^2*y^2 + 24*b*x*y^3 - 4*b*y^4 +
	- 4*b^2*x*y + 5*b^2*y^2 + 4*b*R*y # = 0

### Tr. System:
2*x^3 + 2*y^3 - 3*x*y*(x+y) + b*(2*x - y) - 2*R # = 0
x^6 - 9*x^4*y^2 + 11*x^3*y^3 - 9*x^2*y^4 + y^6 + 2*b*x^4 + 3*b*x^2*y^2 - 4*b*x*y^3 +
	- 2*R*x^3 - 2*R*y^3 + b^2*x^2 - b^2*y^2 - 2*b*R*x + R^2 # = 0


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



