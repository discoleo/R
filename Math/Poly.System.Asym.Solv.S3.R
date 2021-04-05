########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Asymmetric: solvable
###
### draft v.0.2e



####################

###############
### History ###
###############

### draft v.0.2e:
# - solved mixed-leading term system:
#   A * x^4*y^4*z + B %*% c(x^2*y, x*y^2, x*y*z) = R;
### draft v.0.2d-pre - v.0.2d-fix:
# - started work on system of Order 4:
#   A * (x^4+y^4+z^4) + B %*% c(x*y, x*z, y*z) = R;
# - fixed bugs & tested;
### draft v.0.2c:
# - system Order 2:
#   A * (x^2+y^2+z^2) + B %*% c(x*y, x*z, y*z) = R;
### draft v.0.2b:
# - system Order 3:
#   A * (x^3+y^3+z^3) + B %*% c(x,y,z) = R;
### draft v.0.2a:
# - system Order 2:
#   A * (x^2+y^2+z^2) + B %*% c(x,y,z) = R;
### draft v.0.1d - v.0.1d-ext:
# - extensions:
#   A * x*y*z + B %*% c(x^n, y^n, z^n) = R;
#   A * x*y*z + B %*% c(x*y, x*z, y*z) = R;
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
	sol = -sqrt(sol + 0i) # TODO: robust!
	rownames(sol) = c("x","y","z");
	return(sol)
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
solve.E3.S3P1 = function(R, B, a=c(1,1,1), debug=TRUE) {
	m.coeff = solve(B, cbind(R, -1 * a))
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
	matrix(c(sol[1,]*sol[2,], sol[1,]*sol[3,], sol[2,]*sol[3,]), nrow=3, byrow=TRUE)
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


#########
### Ex 3:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)
a = c(1,1,0)

sol = solve.E3.S3P1(R, B, a=a);

### Test
round0(a * rep(apply(sol, 2, prod), each=3) + B %*% E2.f(sol))


################################
################################

### A * (x^n+y^n+z^n) + B %*% c(x,y,z) = R

### Order 2:
### A * (x^2+y^2+z^2) + B %*% c(x,y,z) = R

a1*(x^2+y^2+z^2) + b11*x + b12*y + b13*z = R1
# ...

### Solution:

### Step 1:
# - solve for: x, y, z as f(S2);

### Step 2:
S2 - x^2 - y^2 - z^2 # = 0
S2 - (ra1*S2 + r1)^2 - (ra1*S2 + r1)^2 - (ra1*S2 + r1)^2 # = 0
(ra1^2+ra2^2+ra3^2)*S2^2 + (2*(ra1*r1+ra2*r2+ra3*r3) - 1)*S2 + (r1^2 + r2^2 + r3^2) # = 0


### Solver:
solve.S2.S3P2 = function(R, B, a=c(1,1,1), debug=TRUE) {
	m.coeff = solve(B, cbind(R, -1 * a))
	# solve S2
	coeff = c(sum(m.coeff[,2]^2), (2*sum(m.coeff[,1]*m.coeff[,2]) - 1), sum(m.coeff[,1]^2))
	S2 = roots(coeff);
	if(debug) print(S2);
	# solve (x,y,z)
	sol = sapply(S2, function(s2) m.coeff[,1] + m.coeff[,2]*s2);
	rownames(sol) = c("x","y","z")
	return(sol)
}

### Examples:
R = c(1,2,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.S2.S3P2(R, B);

### Test
rep(apply(sol^2, 2, sum), each=3) + B %*% sol


### Ex 2:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.S2.S3P2(R, B);

### Test
round0(rep(apply(sol^2, 2, sum), each=3) + B %*% sol)


#########
### Ex 3:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)
a = c(1,1,0)

sol = solve.S2.S3P2(R, B, a=a);

### Test
round0(a * rep(apply(sol^2, 2, sum), each=3) + B %*% sol)


############################

### Order 3:
### A * (x^3+y^3+z^3) + B %*% c(x,y,z) = R

a1*(x^3+y^3+z^3) + b11*x + b12*y + b13*z = R1
# ...

### Solution:

### Step 1:
# - solve for: x, y, z as f(S3);

### Step 2:
S3 - x^3 - y^3 - z^3 # = 0
S3 - (ra1*S3 + r1)^3 - (ra1*S3 + r1)^3 - (ra1*S3 + r1)^3 # = 0
(ra1^3+ra2^3+ra3^3)*S3^3 + 3*(ra1^2*r1+ra2^2*r2+ra3^2*r3)*S3^2 +
	+ (3*(ra1*r1^2+ra2*r2^2+ra3*r3^2) - 1)*S3 +
	+ (r1^3 + r2^3 + r3^3) # = 0


### Solver:
solve.S3.S3P3 = function(R, B, a=c(1,1,1), debug=TRUE) {
	m.coeff = solve(B, cbind(R, -1 * a))
	# solve S2
	coeff = c(sum(m.coeff[,2]^3), 3*sum(m.coeff[,1]*m.coeff[,2]^2),
		3*sum(m.coeff[,1]^2*m.coeff[,2]) - 1, sum(m.coeff[,1]^3))
	S3 = roots(coeff);
	if(debug) print(S3);
	# solve (x,y,z)
	sol = sapply(S3, function(s3) m.coeff[,1] + m.coeff[,2]*s3);
	rownames(sol) = c("x","y","z")
	return(sol)
}

### Examples:
R = c(1,2,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.S3.S3P3(R, B);

### Test
rep(apply(sol^3, 2, sum), each=3) + B %*% sol


### Ex 2:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.S3.S3P3(R, B);

### Test
round0(rep(apply(sol^3, 2, sum), each=3) + B %*% sol)


#########
### Ex 3:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)
a = c(1,1,0)

sol = solve.S3.S3P3(R, B, a=a);

### Test
round0(a * rep(apply(sol^3, 2, sum), each=3) + B %*% sol)




################################
################################

### A * (x^n+y^n+z^n) + B %*% c(x*y, x*z, y*z) = R

### Order 2:
### A * (x^2+y^2+z^2) + B %*% c(x*y, x*z, y*z) = R

a1*(x^2+y^2+z^2) + b11*x*y + b12*x*z + b13*y*z = R1
# ...

### Solution:

### Step 1:
# - solve for: x*y, x*z, y*z as f(S2);

### Step 2:
S2 - (x*y)*(x*z)/(y*z) - (x*y)*(y*z)/(x*z) - (x*z)*(y*z)/(x*y) # = 0
(x*y)*(x*z)*(y*z)*S2 - ((x*y)*(x*z))^2 - ((x*y)*(y*z))^2 - ((x*z)*(y*z))^2 # = 0
(ra1*S2+r1)*(ra2*S2+r2)*(ra3*S2+r3)*S2 +
	- ((ra1*S2+r1)*(ra2*S2+r2))^2 - ((ra1*S2+r1)*(ra3*S2+r3))^2 - ((ra2*S2+r2)*(ra3*S2+r3))^2 # = 0
ra1*ra2*ra3*S2^4 + (ra1*ra2*r3+ra1*ra3*r2+ra2*ra3*r1)*S2^3 +
	+ (ra1*r2*r3+ra2*r1*r3+ra3*r1*r2)*S2^2 + r1*r2*r3*S2 +
	- ((ra1*ra2)^2+(ra1*ra3)^2+(ra2*ra3)^2)*S2^4 +
	- 2*(ra1*ra2*(ra1*r2+ra2*r1) + ra1*ra3*(ra1*r3+ra3*r1) + ra2*ra3*(ra2*r3+ra3*r2))*S2^3 +
	- ((ra1*r2+ra2*r1)^2 + (ra1*r3+ra3*r1)^2 + (ra2*r3+ra3*r2)^2 +
		+ 2*ra1*ra2*r1*r2 + 2*ra1*ra3*r1*r3 + 2*ra2*ra3*r2*r3)*S2^2 +
	- 2*(r1*r2*(ra1*r2+ra2*r1) + r1*r3*(ra1*r3+ra3*r1) + r2*r3*(ra2*r3+ra3*r2))*S2 +
	- ((r1*r2)^2 + (r1*r3)^2 + (r2*r3)^2) # = 0
p = prod(m.coeff[,1]); pa = prod(m.coeff[,2]);
# TODO:
(pa - pa^2 * sum(1 / m.coeff[,2]^2))*S2^4 +
	+ pa * sum(m.coeff[,1] / m.coeff[,2])*S2^3 +
	- 2*(sum(m.coeff[,2]^2)*sum(m.coeff[,1]*m.coeff[,2]) - sum(m.coeff[,1]*m.coeff[,2]^3))*S2^3 +
	+ p * sum(m.coeff[,2] / m.coeff[,1])*S2^2 +
	- (sum(m.coeff[,1]^2)*sum(m.coeff[,2]^2) - sum(m.coeff[,1]^2*m.coeff[,2]^2) +
		+ 4*p*pa * sum(1/(m.coeff[,1]*m.coeff[,2])))*S2^2 +
	- 2*(sum(m.coeff[,1]^2)*sum(m.coeff[,1]*m.coeff[,2]) - sum(m.coeff[,1]^3*m.coeff[,2]))*S2 +
	+ p*S2 +
	- p^2 * sum(1 / m.coeff[,1]^2) # = 0


### Solver:
solve.S2E2.S3P2 = function(R, B, a=c(1,1,1), debug=TRUE) {
	m.coeff = solve(B, cbind(R, -1 * a))
	# solve S2
	# TODO: case any() == 0;
	p = prod(m.coeff[,1]); pa = prod(m.coeff[,2]);
	coeff = c(
		(pa - pa^2 * sum(1 / m.coeff[,2]^2)),
		pa * sum(m.coeff[,1] / m.coeff[,2]) +
		- 2*(sum(m.coeff[,2]^2)*sum(m.coeff[,1]*m.coeff[,2]) - sum(m.coeff[,1]*m.coeff[,2]^3)),
		p * sum(m.coeff[,2] / m.coeff[,1]) +
			- (sum(m.coeff[,1]^2)*sum(m.coeff[,2]^2) - sum(m.coeff[,1]^2*m.coeff[,2]^2) +
			+ 4*p*pa * sum(1/(m.coeff[,1]*m.coeff[,2]))),
		- 2*(sum(m.coeff[,1]^2)*sum(m.coeff[,1]*m.coeff[,2]) - sum(m.coeff[,1]^3*m.coeff[,2])) +
			+ p,
		- p^2 * sum(1 / m.coeff[,1]^2)
	)
	S2 = roots(coeff);
	if(debug) print(S2);
	# solve (xy, xz, yz)
	sol = sapply(S2, function(s2) m.coeff[,1] + m.coeff[,2]*s2);
	# x, y, z; TODO: robust
	E3 = apply(sol, 2, prod);
	E3 = sqrt(E3 + 0i); E3 = c(E3, -E3);
	E3 = rep(E3, each=3); sol = cbind(sol, sol);
	sol = E3 / sol[c(3,2,1),]
	rownames(sol) = c("x","y","z")
	return(sol)
}
E2.f = function(sol) {
	matrix(c(sol[1,]*sol[2,], sol[1,]*sol[3,], sol[2,]*sol[3,]), nrow=3, byrow=TRUE)
}

### Examples:
R = c(1,2,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.S2E2.S3P2(R, B);

### Test
rep(apply(sol^2, 2, sum), each=3) + B %*% E2.f(sol)


### Ex 2:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.S2E2.S3P2(R, B);

### Test
round0(rep(apply(sol^2, 2, sum), each=3) + B %*% E2.f(sol))


#########
### Ex 3:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)
a = c(1,1,0)

sol = solve.S2E2.S3P2(R, B, a=a);

### Test
round0(a * rep(apply(sol^2, 2, sum), each=3) + B %*% E2.f(sol))



################################
################################

### Order 4:
### A * (x^4+y^4+z^4) + B %*% c(x*y, x*z, y*z) = R

a1*(x^4+y^4+z^4) + b11*x*y + b12*x*z + b13*y*z = R1
# ...

### Solution:

### Step 1:
# - solve for: x*y, x*z, y*z as f(S4);

### Step 2:
S4 - ((x*y)*(x*z)/(y*z))^2 - ((x*y)*(y*z)/(x*z))^2 - ((x*z)*(y*z)/(x*y))^2 # = 0
((x*y)*(x*z)*(y*z))^2*S4 - ((x*y)*(x*z))^4 - ((x*y)*(y*z))^4 - ((x*z)*(y*z))^4 # = 0
((ra1*S4+r1)*(ra2*S4+r2)*(ra3*S4+r3))^2*S4 +
	- ((ra1*S4+r1)*(ra2*S4+r2))^4 - ((ra1*S4+r1)*(ra3*S4+r3))^4 - ((ra2*S4+r2)*(ra3*S4+r3))^4 # = 0
(ra1*ra2*ra3)^2*S4^7 +
	+ 2*(ra1*ra2*ra3)*(ra1*ra2*r3+ra1*ra3*r2+ra2*ra3*r1)*S4^6 +
	+ (ra1*ra2*r3+ra1*ra3*r2+ra2*ra3*r1)^2*S4^5 +
	+ 2*(ra1*ra2*ra3)*(ra1*r2*r3+ra2*r1*r3+ra3*r1*r2)*S4^5 +
	+ 2*(ra1*ra2*ra3)*(r1*r2*r3)*S4^4 +
	+ 2*(ra1*ra2*r3+ra1*ra3*r2+ra2*ra3*r1)*(ra1*r2*r3+ra2*r1*r3+ra3*r1*r2)*S4^4 +
	+ (ra1*r2*r3+ra2*r1*r3+ra3*r1*r2)^2*S4^3 +
	+ 2*(ra1*ra2*r3+ra1*ra3*r2+ra2*ra3*r1)*(r1*r2*r3)*S4^3 +
	+ 2*(ra1*r2*r3+ra2*r1*r3+ra3*r1*r2)*(r1*r2*r3)*S4^2 +
	+ (r1*r2*r3)^2*S4 +
	# (ra1^2*ra2^2*S4^2 + (ra1*r2+ra2*r1)*S4 + r1*r2)^4
	- ((ra1*ra2)^4+(ra1*ra3)^4+(ra2*ra3)^4)*S4^8 +
	- 4*((ra1*ra2)^3*(ra1*r2+ra2*r1) + (ra1*ra3)^3*(ra1*r3+ra3*r1) +
		+ (ra2*ra3)^3*(ra2*r3+ra3*r2))*S4^7 +
	- 4*((ra1*ra2)^3*r1*r2 + (ra1*ra3)^3*r1*r3 + (ra2*ra3)^3*r2*r3)*S4^6 +
	- 6*((ra1*ra2)^2*(ra1*r2+ra2*r1)^2 + (ra1*ra3)^2*(ra1*r3+ra3*r1)^2 +
		+ (ra2*ra3)^2*(ra2*r3+ra3*r2)^2)*S4^6 +
	- 12*((ra1*ra2)^2*(ra1*r2+ra2*r1)*r1*r2 + (ra1*ra3)^2*(ra1*r3+ra3*r1)*r1*r3 +
		+ (ra2*ra3)^2*(ra2*r3+ra3*r2)*r2*r3)*S4^5 +
	- 4*((ra1*ra2)*(ra1*r2+ra2*r1)^3 + (ra1*ra3)*(ra1*r3+ra3*r1)^3 +
		+ (ra2*ra3)*(ra2*r3+ra3*r2)^3)*S4^5 +
	- 6*((ra1*ra2)^2*(r1*r2)^2 + (ra1*ra3)^2*(r1*r3)^2 + (ra2*ra3)^2*(r2*r3)^2)*S4^4 +
	- 12*((ra1*ra2)*(ra1*r2+ra2*r1)^2*r1*r2 + (ra1*ra3)*(ra1*r3+ra3*r1)^2*r1*r3 +
		+ (ra2*ra3)*(ra2*r3+ra3*r2)^2*r2*r3)*S4^4 +
	- ((ra1*r2+ra2*r1)^4 + (ra1*r3+ra3*r1)^4 + (ra2*r3+ra3*r2)^4)*S4^4 +
	- 12*((ra1*ra2)*(ra1*r2+ra2*r1)*(r1*r2)^2 + (ra1*ra3)*(ra1*r3+ra3*r1)*(r1*r3)^2 +
		+ (ra2*ra3)*(ra2*r3+ra3*r2)*(r2*r3)^2)*S4^3 +
	- 4*((ra1*r2+ra2*r1)^3*r1*r2 + (ra1*r3+ra3*r1)^3*r1*r3 + (ra2*r3+ra3*r2)^3*r2*r3)*S4^3 +
	- 4*((ra1*ra2)*(r1*r2)^3 + (ra1*ra3)*(r1*r3)^3 + (ra2*ra3)*(r2*r3)^3)*S4^2 +
	- 6*((ra1*r2+ra2*r1)^2*(r1*r2)^2 + (ra1*r3+ra3*r1)^2*(r1*r3)^2 +
		+ (ra2*r3+ra3*r2)^2*(r2*r3)^2)*S4^2 +
	- 4*((ra1*r2+ra2*r1)*(r1*r2)^3 + (ra1*r3+ra3*r1)*(r1*r3)^3 +
		+ (ra2*r3+ra3*r2)*(r2*r3)^3)*S4 +
	- ((r1*r2)^4 + (r1*r3)^4 + (r2*r3)^4) # = 0
p = prod(m.coeff[,1]); pa = prod(m.coeff[,2]);
# TODO: compact coefficients;


### Solver:
solve.S4E2.S3P4 = function(R, B, a=c(1,1,1), debug=TRUE) {
	m.coeff = solve(B, cbind(R, -1 * a))
	# solve S4
	# TODO: case any() == 0;
	p = prod(m.coeff[,1]); pa = prod(m.coeff[,2]);
	# TODO: coeffs;
	coeff = coeff.S4E2.S3P4(m.coeff);
	if(debug) print(m.coeff);
	if(debug) print(coeff);
	S4 = roots(coeff);
	if(debug) print(S4);
	# solve (xy, xz, yz)
	sol = sapply(S4, function(s) m.coeff[,1] + m.coeff[,2]*s);
	# x, y, z; TODO: robust
	E3 = apply(sol, 2, prod);
	E3 = sqrt(E3 + 0i); E3 = c(E3, -E3);
	E3 = rep(E3, each=3); sol = cbind(sol, sol);
	sol = E3 / sol[c(3,2,1),]
	rownames(sol) = c("x","y","z")
	return(sol)
}
E2.f = function(sol) {
	matrix(c(sol[1,]*sol[2,], sol[1,]*sol[3,], sol[2,]*sol[3,]), nrow=3, byrow=TRUE)
}

### Examples:
R = c(1,2,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.S4E2.S3P4(R, B);

### Test
rep(apply(sol^4, 2, sum), each=3) + B %*% E2.f(sol)


### Ex 2:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.S4E2.S3P4(R, B);

### Test
round0(rep(apply(sol^4, 2, sum), each=3) + B %*% E2.f(sol))


#########
### Ex 3:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)
a = c(1,1,0)

sol = solve.S4E2.S3P4(R, B, a=a);

### Test
round0(a * rep(apply(sol^4, 2, sum), each=3) + B %*% E2.f(sol))


################################
################################

### Special/Artificial Systems

### Order 4: Mixt Term
### A * x^4*y^4*z + B %*% c(x^2*y, x*y^2, x*y*z) = R

a1*x^4*y^4*z + b11*x^2*y + b12*x*y^2 + b13*x*y*z = R1
# ...

### Solution:

### Step 1:
# - solve for: x^2*y, x*y^2, x*y*z as f(Sa),
#   where Sa = x^4*y^4*z;

### Step 2:
Sa - (x^2*y) * (x*y^2) * (x*y*z) # = 0
Sa - (ra1*Sa + r1)*(ra2*Sa + r2)*(ra2*Sa + r2) # = 0
# rescale r[i] by 1/ra[i]; a = 1 / (ra1*ra2*ra3);
a*Sa - (Sa + r1)*(Sa + r2)*(Sa + r2) # = 0
Sa^3 + (r1+r2+r3)*Sa^2 + (r1*r2+r1*r3+r2*r3 - a)*Sa + r1*r2*r3 # = 0

### Step 3:
# - solve:
# x^2*y = ra1*Sa + r1;
# x*y^2 = ra2*Sa + r2;
# x*y*z = ra3*Sa + r3;


### Solver
solve.E3V.S3P441 = function(R, B, a=c(1,1,1), debug=TRUE) {
	m.coeff = solve(B, cbind(R, -1 * a))
	# print(m.coeff)
	r = m.coeff[,1] / m.coeff[,2];
	a = 1 / prod(m.coeff[,2]);
	E3 = roots(c(1, sum(r), (r[1]*r[2]+r[1]*r[3]+r[2]*r[3] - a), prod(r)));
	if(debug) print(E3);
	sol = sapply(E3, function(e3) m.coeff[,1] + m.coeff[,2]*e3)
	m = unity(3, all=TRUE);
	p = sol[1,] * sol[2,];
	xy = rootn(p, 3); xy = sapply(xy, function(xy) xy*m);
	sol = sol[, rep(seq(ncol(sol)), each=3)];
	x = as.vector(sol[1,] / xy);
	y = as.vector(sol[2,] / xy);
	z = as.vector(sol[3,] / xy);
	sol = rbind(x=x, y=y, z=z);
	return(sol)
}
calc = function(v) {
	c(v[1]^2*v[2], v[1]*v[2]^2, prod(v))
}

### Examples:
R = c(1,2,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.E3V.S3P441(R, B);

### Test
rep(apply(sol^c(4,4,1), 2, prod), each=3) + B %*% apply(sol, 2, calc)


### Ex 2:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)

sol = solve.E3V.S3P441(R, B);

### Test
round0(rep(apply(sol^c(4,4,1), 2, prod), each=3) + B %*% apply(sol, 2, calc))


#########
### Ex 3:
R = c(0,-1,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)
a = c(1,1,0)

sol = solve.E3V.S3P441(R, B, a=a);

### Test
round0(a * rep(apply(sol^c(4,4,1), 2, prod), each=3) + B %*% apply(sol, 2, calc))


