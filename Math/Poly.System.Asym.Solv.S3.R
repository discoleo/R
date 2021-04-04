########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Asymmetric: solvable
###
### draft v.0.1d



####################

###############
### History ###
###############

### draft v.0.1d:
# - extensions:
#   A * x*y*z + B %*% c(x^n, y^n, z^n) = R;
### draft v.0.1c:
# - system:
#   x*y*z + B %*% c(x*y, x*z, y*z) = R;
### draft v.0.1b:
# - system:
#   x*y*z + B %*% c(x^2, y^2, z^2) = R;
### draft v.0.1a:
# - system:
#   x*y*z + B %*% c(x,y,z) = R;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R


################################
################################

### x*y*z + B %*% c(x,y,z) = R

x*y*z + b11*x + b12*y + b13*z = R1
# ...

### Solution:

### Step 1: solve for: x, y, z;

### Step 2:
(r1 + E3)*(r2 + E3)*(r3 + E3) - a*E3 # = 0
E3^3 + (r1+r2+r3)*E3^2 + (r1*r2+r1*r3+r2*r3 - a)*E3 + r1*r2*r3 # = 0

### Solver
solve.E3.S3P1 = function(R, B, a=c(1,1,1), debug=TRUE) {
	m.coeff = solve(B, cbind(R, -1 * a))
	# print(m.coeff)
	r = m.coeff[,1] / m.coeff[,2];
	a = 1 / prod(m.coeff[,2]);
	E3 = roots(c(1, sum(r), (r[1]*r[2]+r[1]*r[3]+r[2]*r[3] - a), prod(r)));
	if(debug) print(E3);
	sol = sapply(E3, function(e3) m.coeff[,1] + m.coeff[,2]*e3)
	return(sol)
}

### Examples:
R = c(1,2,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.E3.S3P1(R, B);

### Test
rep(apply(sol, 2, prod), each=3) + B %*% sol


### Ex 2:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.E3.S3P1(R, B);

### Test
round0(rep(apply(sol, 2, prod), each=3) + B %*% sol)


### Ex 3:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)
a = c(1,1,0)

sol = solve.E3.S3P1(R, B, a=a);

### Test
round0(a * rep(apply(sol, 2, prod), each=3) + B %*% sol)



############################

### x*y*z + B %*% c(x^2, y^2, z^2) = R

x*y*z + b11*x^2 + b12*y^2 + b13*z^2 = R1
# ...

### Solution:

### Step 1: solve for: x^2, y^2, z^2;

### Step 2:
(r1 + E3)*(r2 + E3)*(r3 + E3) - a*E3^2 # = 0
E3^3 + (r1+r2+r3 - a)*E3^2 + (r1*r2+r1*r3+r2*r3)*E3 + r1*r2*r3 # = 0

### Solver
solve.E3.S3P1 = function(R, B, a=c(1,1,1), debug=TRUE) {
	m.coeff = solve(B, cbind(R, -1 * a))
	r = m.coeff[,1] / m.coeff[,2];
	a = 1 / prod(m.coeff[,2]);
	E3 = roots(c(1, sum(r) - a, (r[1]*r[2]+r[1]*r[3]+r[2]*r[3]), prod(r)));
	if(debug) print(E3);
	sol = sapply(E3, function(e3) m.coeff[,1] + m.coeff[,2]*e3)
	return(-sqrt(sol + 0i)) # TODO: robust!
}

### Examples:
R = c(1,2,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.E3.S3P1(R, B);

### Test
rep(apply(sol, 2, prod), each=3) + B %*% (sol^2)


### Ex 2:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.E3.S3P1(R, B);

### Test
round0(rep(apply(sol, 2, prod), each=3) + B %*% sol^2)


### Ex 3:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)
a = c(1,1,0)

sol = solve.E3.S3P1(R, B, a=a);

### Test: TODO: robust!
round0(a * rep(apply(sol, 2, prod), each=3) + B %*% sol^2)


############################
############################

### x*y*z + B %*% c(x*y, x*z, y*z) = R

x*y*z + b11*x*y + b12*x*z + b13*y*z = R1
# ...

### Solution:

### Step 1: solve for: x*y, x*z, y*z;

### Step 2:
(r1 + E3)*(r2 + E3)*(r3 + E3) - a*E3^2 # = 0
E3^3 + (r1+r2+r3 - a)*E3^2 + (r1*r2+r1*r3+r2*r3)*E3 + r1*r2*r3 # = 0

### Step 3:
# x = E3 / (y*z);

### Solver
solve.E3.S3P1 = function(R, B, debug=TRUE) {
	m.coeff = solve(B, cbind(R, -1 * c(1,1,1)))
	r = m.coeff[,1] / m.coeff[,2];
	a = 1 / prod(m.coeff[,2]);
	E3 = roots(c(1, sum(r) - a, (r[1]*r[2]+r[1]*r[3]+r[2]*r[3]), prod(r)));
	if(debug) print(E3);
	sol = sapply(E3, function(e3) m.coeff[,1] + m.coeff[,2]*e3)
	E3 = rep(E3, each=3);
	sol = E3 / sol[c(3,2,1),];
	rownames(sol) = c("x","y","z")
	return(sol)
}
E2.f = function(sol) {
	matrix(c(sol[1,]*sol[2,], sol[1,]*sol[3,], sol[2,]*sol[3,]), ncol=3, byrow=TRUE)
}

### Examples:
R = c(1,2,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.E3.S3P1(R, B);

### Test
rep(apply(sol, 2, prod), each=3) + B %*% E2.f(sol)


### Ex 2:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.E3.S3P1(R, B);

### Test
round0(rep(apply(sol, 2, prod), each=3) + B %*% E2.f(sol))

