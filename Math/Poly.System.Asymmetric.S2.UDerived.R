########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S2
### Asymmetric Derived from Symmetric
###   Based on Roots of Unity
###
### draft v.0.1c


#######################

### Helper Functions

source("Polynomials.Helper.R")


#######################

### Base System:
# x^3 + y^3 + b*(x+y) = R1
# x*y = R2

### Derived:
# x => x + m*y
# y => x + m^2*y
# where m^3 = 1

### Derived System
2*x^3 + 2*y^3 - 3*x*y*(x + y) + b*(2*x - y) - R[1] # = 0
x^2 + y^2 - x*y - R[2] # = 0


### Solver
solve.S2As.P3 = function(R, b, debug=TRUE) {
	m = unity(3, all=FALSE);
	# Step 1:
	# S^3 - (3*R2 - b)*S - R1 = 0
	S = roots(c(1,0, - 3*R[2] + b, -R[1]));
	if(debug) print(S);
	# Step 2:
	d = rootn(S^2 - 4*R[2], 2);
	d = c(d, -d); S = c(S, S);
	# Step 3:
	y = d / (m - m^2);
	x = (S + y) / 2;
	cbind(x=x, y=y);
}
test.S2As.P3 = function(sol, R=NULL, b) {
	err1 = 2*x^3 + 2*y^3 - 3*x*y*(x + y) + b*(2*x - y);
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


#######################
#######################

### Hetero-Symmetric

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
	S = roots(c(1, 0, -2*b[1], R));
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
test.DerHt.P3 = function(sol, R = NULL, b) {
	x = sol[,1]; y = sol[,2];
	err1 = 2*x^3 + 2*y^3 - 3*x*y*(x+y) + b*(2*x - y);
	err2 = 3*x*(x - y) - b;
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

