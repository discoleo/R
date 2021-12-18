########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Symmetric S4
###
### draft v.0.1b



####################

### Helper Functions

source("Polynomials.Helper.R")


### Other

test.S4Symm = function(sol, n=2, R = NULL) {
	x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];
	err1 = x1^n + x2^n + x3^n + x4^n;
	err2 = x1*x2 + x1*x3 + x1*x4 + x2*x3 + x2*x4 + x3*x4;
	err3 = x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4;
	err4 = x1*x2*x3*x4;
	err = rbind(err1, err2, err3, err4);
	if( ! is.null(R)) {
		for(id in 1:4) err[id,] = err[id,] - R[id];
	}
	err = round0(err);
	return(err);
}

###############
###############

###############
### Order 2 ###
###############

x1^2 + x2^2 + x3^2 + x4^2 - R1 # = 0
x1*x2 + x1*x3 + x1*x4 + x2*x3 + x2*x4 + x3*x4 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

### Solver:
solve.S4Symm.P2 = function(R) {
	S2 = R[1] + 2*R[2];
	S = rootn(S2, 2); S = c(S, -S);
	x = sapply(S, function(S) roots(c(1, -S, R[2], -R[3], R[4])));
	# All permutations possible:
	# - only 1 variant generated!
	perm.f = function(id) {
		id = c(seq(id+1, 4), seq(1, id));
		as.vector(x[id, ]);
	}
	sol = cbind(x1=as.vector(x), x2 = perm.f(1), x3 = perm.f(2), x4 = perm.f(3));
	return(sol);
}

### Examples:

R = c(1,-1,2,3)
sol = solve.S4Symm.P2(R);
test.S4Symm(sol, n=2);

### Classic Poly:
round0(poly.calc(sol[,1]))


### Ex 2:
R = c(0,-1,1,2)
sol = solve.S4Symm.P2(R);
test.S4Symm(sol, n=2);

### Classic Poly:
round0(poly.calc(sol[,1]))


####################
####################

###############
### Order 3 ###
###############

x1^3 + x2^3 + x3^3 + x4^3 - R1 # = 0
x1*x2 + x1*x3 + x1*x4 + x2*x3 + x2*x4 + x3*x4 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

### Solver:
solve.S4Symm.P3 = function(R, sort=TRUE) {
	S = roots(c(1, 0, -3*R[2], - R[1] + 3*R[3]));
	x = sapply(S, function(S) roots(c(1, -S, R[2], -R[3], R[4])));
	# All permutations possible:
	# - only 1 variant generated!
	perm.f = function(id) {
		id = c(seq(id+1, 4), seq(1, id));
		as.vector(x[id, ]);
	}
	sol = cbind(x1=as.vector(x), x2 = perm.f(1), x3 = perm.f(2), x4 = perm.f(3));
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}

### Examples:

R = c(1,-1,2,3)
sol = solve.S4Symm.P3(R);
test.S4Symm(sol, n=3);

### Classic Poly:
round0(poly.calc(sol[,1]))

