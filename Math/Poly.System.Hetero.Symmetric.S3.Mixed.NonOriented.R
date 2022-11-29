########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogeneous Symmetric S3:
### Mixed Type: Non-Oriented
###
### draft v.0.1a


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


### Solver:

solve.S3HtM.P21NonD = function(R, b.ext=0, debug=TRUE) {
	# TODO
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
sol = solve.S3HtM.DualSum21(R, a=a)

test.S3HtM.P21(sol, a=a)

### Debug
R = c(0, 2, 3)
a = 1

x =  1.089990536315 - 1.250695049316i;
y = -0.148968812161 + 1.079762780452i;
z =  1.089990536315 - 1.250695049316i;
sol = cbind(x, y, z)
# every permutation is valid;
sol = rbind(sol, sol[c(1,3,2)])
test.S3HtM.P21(sol, a=a)


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


### Test
(x*y^2 + y*z^2 + z*x^2) - (x*z^2 + y*x^2 + z*y^2) # = 0
x*y + x*z + y*z # - R[2] # = 0
x*y*z # - R[3] # = 0

