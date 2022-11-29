########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogeneous Symmetric S3:
### Mixed Type: Non-Oriented
###
### draft v.0.1b


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

solve.S3HtM.P21NonD = function(R, b.ext=0, debug=TRUE, all=FALSE) {
	# assumes R[1] == 0;
	E2 = R[2]; E3 = R[3];
	coeff = c(4*E3, - E2^2, - 18*E2*E3, 4*E2^3 + 27*E3^2);
	S = roots(coeff);
	if(debug) print(S);
	x = sapply(S, function(S) roots(c(1, -S, E2, -E3)));
	x = as.vector(x);
	# Robust:
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
test.S3HtM.P21 = function(sol, a=1, b.ext=0, R=NULL, tol=1E-8) {
	x = sol[,1]; y = sol[,2]; z = sol[,3];
	err1 = (x*y^2 + y*z^2 + z*x^2) - a*(x*z^2 + y*x^2 + z*y^2);
	err2 = x*y + x*z + y*z;
	err3 = x*y*z;
	err = rbind(err1, err2, err3);
	if( ! is.null(R)) {
		err = err - R;
	}
	err = round0(err, tol=tol);
	return(err);
}

### Examples:

R = c(0, 2, 3)
sol = solve.S3HtM.P21NonD(R)

test.S3HtM.P21(sol, a=a)


### Debug
R = c(0, 2, 3)
a = 1

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
test.S3HtM.P21(sol, a=a)


### Test
(x*y^2 + y*z^2 + z*x^2) - (x*z^2 + y*x^2 + z*y^2) # = 0
x*y + x*z + y*z # - R[2] # = 0
x*y*z # - R[3] # = 0


### Numerical solver

solve.S3HtM.Num = function(x, R, a=1) {
	x = matrix(x, ncol=3);
	xc = x[2,]; x = x[1,] + 1i*xc;
	x = matrix(x, nrow=1);
	y = test.S3HtM.P21(matrix(x, nrow=1), R=R, a=a, tol=1E-15);
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
# redundant!
(x^2 - s*x + e2)*(y-z) # = 0

# E21a:
2*(x^2*y + y^2*z + z^2*x) - E2*S + 3*E3 # = 0
2*(x^2*y + (s*y - e2)*z + (s*z - e2)*x) - E2*S + 3*E3 # = 0
2*(x^2*y + (s*y - e2)*(s - y) + (s*(s - y) - e2)*x) - E2*S + 3*E3 # = 0
# redundant!
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

