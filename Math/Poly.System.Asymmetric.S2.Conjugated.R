########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Asymmetric S2:
### Conjugated Types
###
### draft v.0.1b-fix


### Asymmetric Polynomial Systems: 2 Variables
### Conjugated Types


###############
### History ###
###############


### draft v.0.1a - v.0.1b:
# - [solved] Simple Order z=2;
# - [solved] Simple Order z=3;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R;
# e.g. round0(), round0.p;
source("Polynomials.Helper.R")

### Other
# [...]


##########################

##########################
### Polynomial Systems ###
##########################

###################
### Order z = 2 ###
###################

# based on:
# [Michael Penn] Something is hiding in this problem...
# https://www.youtube.com/watch?v=zBWgNf4Jl_A

### System:
x*(x^2 + y^2) - R1*(x^2 + y^2) + b*x # = 0
y*(x^2 + y^2) - R2*(x^2 + y^2) - b*y # = 0

### Solution:

### Trivial solution:
x = y = 0

### Non-Trivial solutions:

### let:
z = x + 1i*y; zbar = x - 1i*y;
# =>
x*z*zbar - R1*z*zbar + b*x # = 0
y*z*zbar - R2*z*zbar - b*y # = 0
### Sum(Eq 1 + 1i * Eq 2) =>
z^2*zbar - (R1+R2*1i)*z*zbar + b*zbar # = 0
# zbar != 0 =>
z^2 - (R1+R2*1i)*z + b # = 0

### x = z - 1i*y # =>
y*z*(z - 2i*y) - R2*z*(z - 2i*y) - b*y # = 0
- 2i*z*y^2 + (z^2 + 2i*R2*z - b)*y - R2*z^2 # = 0

### Solver:
solve.S2Conj2 = function(R, b) {
	zb1 = (R[1] + 1i*R[2]);
	det = sqrt(zb1^2 - 4*b[1]);
	z = c((zb1 + det) / 2, (zb1 - det)/2);
	# Simple roots:
	# x = Re(z); y = Im(z);
	# sol1 = cbind(x=x, y=y);
	# All solutions:
	y = sapply(z, function(z) {
		coeff = c(- 2i*z, (z^2 + 2i*R2*z - b), - R2*z^2);
		roots(coeff);
	})
	z = rep(z, each=2);
	x = z - 1i*y;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	# sol = rbind(sol, sol1);
	return(sol);
}

### Examples:
R = c(1,2)
b = 3;
#
sol = solve.S2Conj2(R, b);
x = sol[,1]; y = sol[,2];

### Test:
R1 = R[1]; R2 = R[2];
err1 = x*(x^2 + y^2) - R1*(x^2 + y^2) + b*x
err2 = y*(x^2 + y^2) - R2*(x^2 + y^2) - b*y
round0(rbind(err1, err2))


###################

###################
### Order z = 3 ###
###################

### System:
(x^2 - y^2)*(x^2 + y^2) + b2*x*(x^2 + y^2) - R1*(x^2 + y^2) + b1*x # = 0
2*x*y*(x^2 + y^2) + b2*y*(x^2 + y^2) - R2*(x^2 + y^2) - b1*y # = 0

### Solution:

### Trivial solution:
x = y = 0
# TODO: explore more;

### Non-Trivial solutions:

### let:
z = x + 1i*y; zbar = x - 1i*y;
# =>
(x^2 - y^2)*z*zbar + b2*x*z*zbar - R1*z*zbar + b1*x # = 0
2*x*y*z*zbar + b2*y*z*zbar - R2*z*zbar - b1*y # = 0
### Sum(Eq 1 + 1i * Eq 2) =>
z^3*zbar + b2*z^2*zbar - (R1+R2*1i)*z*zbar + b1*zbar # = 0
# zbar != 0 =>
z^3 + b2*z^2 - (R1+R2*1i)*z + b1 # = 0

### x = z - 1i*y # =>
2*(z - 1i*y)*y*z*(z - 2i*y) + b2*y*z*(z - 2i*y) - R2*z*(z - 2i*y) - b1*y # = 0
2*z*y*(2*y^2 + 3i*z*y - z^2) + b2*z*y*(2i*y - z) - R2*z*(2i*y - z) + b1*y # = 0
4*z*y^3 + 2i*(3*z^2 + b2*z)*y^2 - (2*z^3 + b2*z^2 + 2i*R2*z - b1)*y + R2*z^2 # = 0

### Solver:
solve.S2Conj3 = function(R, b) {
	zb1 = (R[1] + 1i*R[2]);
	coeff = c(1, b[2], -zb1, b[1]);
	z = roots(coeff);
	# Simple roots:
	# x = Re(z); y = Im(z);
	# sol1 = cbind(x=x, y=y);
	# All solutions:
	b1 = b[1]; b2 = b[2]; R1 = R[1]; R2 = R[2];
	y = sapply(z, function(z) {
		coeff = c(4*z, 2i*(3*z^2 + b2*z), - (2*z^3 + b2*z^2 + 2i*R2*z - b1), R2*z^2);
		roots(coeff);
	})
	z = rep(z, each=3);
	x = z - 1i*y;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	# sol = rbind(sol, sol1);
	return(sol);
}
test.S2Conj3 = function(sol, R, b) {
	x = sol[,1]; y = sol[,2];
	R1 = R[1]; R2 = R[2]; b1 = b[1]; b2 = b[2];
	err1 = (x^2 - y^2)*(x^2 + y^2) + b2*x*(x^2 + y^2) - R1*(x^2 + y^2) + b1*x # = 0
	err2 = 2*x*y*(x^2 + y^2) + b2*y*(x^2 + y^2) - R2*(x^2 + y^2) - b1*y # = 0
	err = round0(rbind(err1, err2));
	return(err);
}

### Examples:
R = c(1,2)
b = c(-2,3);
#
sol = solve.S2Conj3(R, b);
x = sol[,1]; y = sol[,2];

### Test:
test.S2Conj3(sol, R, b)

### Debug:
R1 = R[1]; R2 = R[2]; b1 = b[1]; b2 = b[2];


