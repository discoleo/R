########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogeneous Symmetric S3:
### Mixed Type: Non-Oriented
###
### draft v.0.1d-sol


### Heterogeneous Symmetric
### Polynomial Systems: 3 Variables
### Mixed: Hetero + Symmetric
### Non-Oriented (or Undirected)


### Non-Oriented Ht:
# => all permutations are valid solutions!

### Example: Simple
# (x^n*y^p + y^n*z^p + z^n*x^p) - (x^p*y^n + y^p*z^n + z^p*x^n) = 0
# x*y + x*z + y*z = R2
# x*y*z = R3

### Example: 2 Ht
# (x^n*y^p + y^n*z^p + z^n*x^p) - (x^p*y^n + y^p*z^n + z^p*x^n) = 0
# x^n2*y + y^n2*z + z^n2*x = R2
# x*y*z = R3


####################
####################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")

### Numerical solver
source("Polynomials.Helper.Solvers.Num.R")

# library(polynom)
# library(pracma)

### Other
E2n.f = function(sol, n) {
	if(length(n) == 1) {
		n1 = n; n2 = 1;
	} else { n1 = n[1]; n2 = n[2]; }
	if(is.null(dim(sol))) sol = matrix(sol, nrow=1);
	apply(sol, 1, function(sol) {
		x = sol[1]; y = sol[2]; z = sol[3];
		x^n1*y^n2 + y^n1*z^n2 + z^n1*x^n2;
	})
}

### Tests:
test.S3HtM.Simple = function(sol, b.ext=0, R=NULL, n, a=1, tol=1E-8) {
	x = sol[,1]; y = sol[,2]; z = sol[,3];
	n1 = n[1]; n2 = n[2];
	err1 = (x^n1*y^n2 + y^n1*z^n2 + z^n1*x^n2) - a*(x^n2*y^n1 + y^n2*z^n1 + z^n2*x^n1);
	err2 = x*y + x*z + y*z;
	err3 = x*y*z;
	err = rbind(err1, err2, err3);
	if( ! is.null(R)) {
		err = err - R;
	}
	err = round0(err, tol=tol);
	return(err);
}
test.S3HtM.Dual = function(sol, b.ext=0, R=NULL, n, a=1, tol=1E-8) {
	x = sol[,1]; y = sol[,2]; z = sol[,3];
	n1 = n[1]; n2 = n[2]; n3 = n[3]; n4 = n[4];
	err1 = (x^n1*y^n2 + y^n1*z^n2 + z^n1*x^n2) - a*(x^n2*y^n1 + y^n2*z^n1 + z^n2*x^n1);
	err2 = (x^n3*y^n4 + y^n3*z^n4 + z^n3*x^n4);
	err3 = x*y*z;
	err = rbind(err1, err2, err3);
	if( ! is.null(R)) {
		err = err - R;
	}
	err = round0(err, tol=tol);
	return(err);
}

######################

####################
### Order E[2,1] ###
####################

# (x^2*y + y^2*z + z^2*x) - (x*y^2 + y*z^2 + z*x^2) = 0

### Solution:

### Note:
# - all solutions are of the form:
#   x == y or x == z or y == z;

### Sum =>
2*E21a - (E2*S - 3*E3) # = 0

### Prod =>
E21a^2 - (x^3*y^3+x^3*z^3+y^3*z^3) - 3*E3^2 - E3*(x^3 + y^3 + z^3) # = 0
E21a^2 - E2^3 - 9*E3^2 - E3*S^3 + 6*E2*E3*S # = 0
(E2*S - 3*E3)^2 - 4*E2^3 - 36*E3^2 - 4*E3*S^3 + 24*E2*E3*S # = 0
E2^2*S^2 - 4*E2^3 - 27*E3^2 - 4*E3*S^3 + 18*E2*E3*S # = 0

### Eq S:
4*E3*S^3 - E2^2*S^2 - 18*E2*E3*S + 4*E2^3 + 27*E3^2 # = 0


### Solver:

# Note:
# - includes automatically all permutations;
solve.S3HtM.P21NonD = function(R, b.ext=0, debug=TRUE, all=TRUE) {
	# assumes R[1] == 0;
	E2 = R[2]; E3 = R[3];
	coeff = c(4*E3, - E2^2, - 18*E2*E3, 4*E2^3 + 27*E3^2);
	S = roots(coeff);
	if(debug) print(S);
	x = sapply(S, function(S) roots(c(1, -S, E2, -E3)));
	x = as.vector(x);
	# quasi-Robust: all triplets are valid roots;
	# - but a lot of numerical instabilities!
	S = rep(S, each=3);
	s = S - x; e2 = E2 - x*s;
	len = length(x);
	y12 = sapply(seq(len), function(id) roots(c(1, -s[id], e2[id])));
	y = y12[1,]; z = y12[2,];
	#
	sol = cbind(x, y, z);
	return(sol);
}
test.S3HtM.P21 = function(sol, b.ext=0, R=NULL, tol=1E-8) {
	test.S3HtM.Simple(sol, b.ext=b.ext, R=R, n=c(2,1), a=1, tol=tol);
}

### Examples:

R = c(0, 2, 3)
sol = solve.S3HtM.P21NonD(R)

# quite some numeric instability!
test.S3HtM.P21(sol, tol=1E-6)


### Ex 2:
R = c(0, -1, 4)
sol = solve.S3HtM.P21NonD(R)

# quite some numeric instability!
test.S3HtM.P21(sol, tol=1E-6)


#########
### Debug
R = c(0, 2, 3)

x =  1.089990536315 - 1.250695049316i;
y = -0.148968812161 + 1.079762780452i;
z =  1.089990536315 - 1.250695049316i;
#
S = x + y + z; E2 = R[2]; E3 = R[3];
E21a = x^2*y + y^2*z + z^2*x;
#
sol = cbind(x, y, z)
# every permutation is valid;
sol = rbind(sol, sol[c(1,3,2)])
test.S3HtM.P21(sol)


### Test
(x*y^2 + y*z^2 + z*x^2) - (x*z^2 + y*x^2 + z*y^2) # = 0
x*y + x*z + y*z # - R[2] # = 0
x*y*z # - R[3] # = 0


### Numerical solver

solve.S3HtM.Num = function(x, R) {
	x = matrix(x, ncol=3);
	xc = x[2,]; x = x[1,] + 1i*xc;
	x = matrix(x, nrow=1);
	y = test.S3HtM.P21(matrix(x, nrow=1), R=R, tol=1E-15);
	y = rbind(Re(y), Im(y));
	y = as.vector(y);
	return(y);
}

x0 = c(1.09-1.2507i, -0.149+1.0798i, 1.09-1.2507i);
x = solve.all(solve.S3HtM.Num, x0, R=R)

### Derivation:
s = S - x; e2 = E2 - s*x;
# Diff =>
x^2*(y-z) + y*z*(y-z) + x*z^2 - x*y^2 # = 0
# Redundant!
# x^2 - s*x + e2 == 0 (6 cases) vs y == z (3 cases);
(x^2 - s*x + e2)*(y-z) # = 0

# E21a:
2*(x^2*y + y^2*z + z^2*x) - E2*S + 3*E3 # = 0
2*(x^2*y + (s*y - e2)*z + (s*z - e2)*x) - E2*S + 3*E3 # = 0
2*(x^2*y + (s*y - e2)*(s - y) + (s*(s - y) - e2)*x) - E2*S + 3*E3 # = 0
# Redundant!
2*(x^2 - s*x + e2)*y + 3*E3 - S*E2 - 2*e2*x + 2*s^2*x # = 0

# Alternative Solution:
# x = y =>
# System:
x^2 + 2*x*z - R2 # = 0
x^2*z - R3 # = 0
# =>
x^3 - R2*x + 2*R3 # = 0
# z = R3 / x^2;
# Note: formula for z works only for this method (x == y);


####################
####################

####################
### Order E[3,1] ###
####################

# (x^3*y + y^3*z + z^3*x) - (x^3*z + y^3*x + z^3*y) = 0

### Solution:

### Note:
# Case 1: x == y or x == z or y == z;
# Case 2: S = 0;

### Case 1: x == y;
# - same as E[2,1], see previous section;

### Case 2: S = 0;
# - trivial P[3] in x;

### Solver:

solve.S3HtM.P31 = function(R, b.ext=NULL, debug=TRUE, all=FALSE) {
	# Note: does NOT include solutions for E[2,1]!
	S = 0; E2 = R[2]; E3 = R[3];
	x = roots(c(1, -S, E2, -E3));
	#
	s = S - x; e2 = E2 - s*x;
	len = length(x);
	y12 = sapply(seq(len), function(id) roots(c(1, -s[id], e2[id])) );
	y = y12[1,]; z = y12[2,];
	sol = cbind(x, y, z);
	if(all) {
		sol = rbind(sol, sol[, c(1,3,2)]);
	}
	return(sol);
}
test.S3HtM.P31 = function(sol, b.ext=0, R=NULL, tol=1E-8) {
	test.S3HtM.Simple(sol, b.ext=b.ext, R=R, n=c(3,1), a=1, tol=tol);
}

### Examples:

R = c(0,2,3)
sol = solve.S3HtM.P31(R);

test.S3HtM.P31(sol)


### Ex 2:
R = c(0,2,-1)
sol = solve.S3HtM.P31(R);

test.S3HtM.P31(sol)


##########
### Debug:
R = c(0,2,3)
x =  1.089990536315 - 1.250695049316i;
y = -0.148968812161 + 1.079762780452i
z =  1.089990536315 - 1.250695049316i;
sol = cbind(x, y, z); # x == z;

###
x = -0.5 + 1.658312395178i;
y = -0.5 - 1.658312395178i;
z =  1;
sol = rbind(sol, c(x, y, z));

test.S3HtM.P31(sol)


### Numerical solver

solve.S3HtM.Num = function(x, R, isSZero=FALSE) {
	x = matrix(x, ncol=3);
	xc = x[2,]; x = x[1,] + 1i*xc;
	x = matrix(x, nrow=1);
	y = test.S3HtM.P31(matrix(x, nrow=1), R=R, tol=1E-15);
	if(isSZero) y[1] = sum(x);
	y = rbind(Re(y), Im(y));
	y = as.vector(y);
	return(y);
}

x0 = c(1.09-1.2507i, -0.149+1.0798i, 1.09-1.2507i); # x == z;
x0 = c(-1.09-1.2507i, 0.149+0.0798i, 1.09-1.2507i); # S = 0;
x = solve.all(solve.S3HtM.Num, x0, R=R, isSZero=TRUE)


##########################
##########################

##########################
### Dual: Non-Directed ###
###  & Simple          ###
##########################

### Order E[2,1] & E[3,1]
# E21a - E21b = 0
# E31 = R2
# x*y*z = R3


### Solution:

###
(E21a - E21b)*S # = 0
E31a - E31b # = 0
# =>
E31a + E31b - (E2*S^2 - E3*S - 2*E2^2) # = 0
E2*S^2 - E3*S - 2*E2^2 - 2*R2 # = 0

### Prod =>
E31a*E31b - (2*E3^2*S^2 - 4*E2^2*E3*S + E2^4 + 4*E2*E3^2) +
	- E3^2*E2 - E3*(x^5 + y^5 + z^5) # = 0
2*E3^2*S^2 - 4*E2^2*E3*S + E2^4 + 5*E2*E3^2 +
	+ E3*(x^5 + y^5 + z^5) - R2^2 # = 0
E3*S^5 - 5*E2*E3*S^3 + 7*E3^2*S^2 + E2^2*E3*S + E2^4 - R2^2 # = 0

### Eq S:
R2 = R[2];
E3*S^9 - 26*E3^2*S^6 - 26*R2*E3*S^5 - R2^2*S^4 + 119*E3^3*S^3 + 372*R2*E3^2*S^2 + 168*R2^2*E3*S +
	+ 16*R2^3 + 729*E3^4 # = 0


### Solver:

solve.S3HtM.nE21E31 = function(R, debug=TRUE, all=TRUE) {
	# assumes R1 == 0;
	R2 = R[2]; E3 = R[3];
	coeff = c(E3, 0, 0, - 26*E3^2, - 26*R2*E3, - R2^2, 119*E3^3, 372*R2*E3^2, 168*R2^2*E3,
		16*R2^3 + 729*E3^4);
	S = roots(coeff);
	if(debug) print(S);
	E2 = - (7*E3*S^3 - 2*R2*S^2 + 54*E3^2) / (S^4 - 40*E3*S - 8*R2);
	#
	len = length(S);
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3)));
	x = as.vector(x);
	S = rep(S, each=3); E2 = rep(E2, each=3);
	#
	s = S - x; e2 = E2 - s*x;
	len = length(x);
	y12 = sapply(seq(len), function(id) roots(c(1, -s[id], e2[id])));
	y = y12[1,]; z = y12[2,];
	#
	sol = cbind(x, y, z);
	return(sol);
}
# quite some numerical instability;
test.S3HtM.nE21E31 = function(sol, R=NULL, tol=1E-6) {
	test.S3HtM.Dual(sol, R=R, n=c(2,1,3,1), tol=tol);
}

### Examples:

R = c(0, 2, 3)
sol = solve.S3HtM.nE21E31(R)

test.S3HtM.nE21E31(sol)

### Ex 2:
R = c(0, -1, 4)
sol = solve.S3HtM.nE21E31(R)

test.S3HtM.nE21E31(sol)


#########
### Debug
R = c(0, 2, 3)

x =  0.725162671785 - 1.419122695717i;
y = -0.692069115990 + 0.957233337703i;
z =  0.725162671785 - 1.419122695717i;
#
sol = cbind(x, y, z)
S = x + y + z; E3 = R[3];
E2 = x*y + x*z + y*z;
E21a = E2n.f(sol, n=c(2,1));
E21b = E2n.f(sol, n=c(1,2));
E31a = E2n.f(sol, n=c(3,1));
E31b = E2n.f(sol, n=c(1,3));

# every permutation is valid;
sol = rbind(sol, sol[c(1,3,2)])
test.S3HtM.nE21E31(sol)



### Numerical solver

solve.S3HtM.nE21E31Num = function(x, R, a=1) {
	x = matrix(x, ncol=3);
	xc = x[2,]; x = x[1,] + 1i*xc;
	x = matrix(x, nrow=1);
	y = test.S3HtM.nE21E31(x, R=R, tol=1E-15);
	y = rbind(Re(y), Im(y));
	y = as.vector(y);
	return(y);
}

x0 = c(0.72-1.4i, -0.69+0.96i, 0.72-1.4i);
x = solve.all(solve.S3HtM.nE21E31Num, x0, R=R)


### Derivation

p1 = toPoly.pm("E2*S^2 - E3*S - 2*E2^2 - 2*R2")
p2 = toPoly.pm("E3*S^5 - 5*E2*E3*S^3 + 7*E3^2*S^2 + E2^2*E3*S + E2^4 - R2^2")

pR = solve.pm(p1, p2, "E2")
pR$Rez$coeff = - pR$Rez$coeff;
print.pm(pR$Rez, "S")


### Classic Poly
# P[9] * P[9]

# P[9] = under the assumption: x == y;
# - Note: fails for the 9 triples where x is the distinct value;
x = sol[,1]; R2 = R[2]; R3 = R[3];
err = x^9 + R3*x^6 - R2*x^5 + R3^3 # = 0
round0(err, tol=1E-4)

