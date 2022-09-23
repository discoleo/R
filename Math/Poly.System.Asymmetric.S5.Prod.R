########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Asymmetric
### Product Type
###
### draft v.0.1a


### Product Type
# - can be solved using the product of the eqs;

# and now for the many variants
# and roots of the asymmetric types!

####################
####################

### Helper Functions

source("Polynomials.Helper.R")


####################

### History

### v.0.1a
# - moved specific code from file:
#   Poly.System.Hetero.Symmetric.S4.R


#############################
#############################

### Lead: x1^3*x2

# x1^3*x2 + b1*x1*x2*x3*x4 = R1
# x2^3*x3 + b2*x1*x2*x3*x4 = R2
# x3^3*x4 + b3*x1*x2*x3*x4 = R3
# x4^3*x1 + b4*x1*x2*x3*x4 = R4

### Solution:

### =>
# x1^3*x2 = R1 - b1*x1*x2*x3*x4
### Prod =>
(x1*x2*x3*x4)^4 - b1*b2*b3*b4*(x1*x2*x3*x4)^4 +
	+ b1*b2*b3*b4*(R1/b1 + R2/b2 + R3/b3 + R4/b4)*(x1*x2*x3*x4)^3 +
	- (b1*b2*R3*R4 + b1*b3*R2*R4 + b1*b4*R2*R3 + b2*b3*R1*R4 + b2*b4*R1*R3 + b3*b4*R1*R2)*
		(x1*x2*x3*x4)^2 +
	+ R1*R2*R3*R4*(b1/R1 + b2/R2 + b3/R3 + b4/R4)*(x1*x2*x3*x4) +
	- R1*R2*R3*R4 # = 0

### Special Case: b1*b2*b3*b4 = 1
b1*b2*b3*b4*(R1/b1 + R2/b2 + R3/b3 + R4/b4)*(x1*x2*x3*x4)^3 +
	- (b1*b2*R3*R4 + b1*b3*R2*R4 + b1*b4*R2*R3 + b2*b3*R1*R4 + b2*b4*R1*R3 + b3*b4*R1*R2)*
		(x1*x2*x3*x4)^2 +
	+ R1*R2*R3*R4*(b1/R1 + b2/R2 + b3/R3 + b4/R4)*(x1*x2*x3*x4) +
	- R1*R2*R3*R4 # = 0

### Solver
solve.Pr2.S4P31 = function(R, b, debug=TRUE) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	b1 = b[1]; b2 = b[2]; b3 = b[3]; b4 = b[4];
	coeff = c(1- b1*b2*b3*b4,
		b1*b2*b3*b4*(R1/b1 + R2/b2 + R3/b3 + R4/b4),
		- (b1*b2*R3*R4 + b1*b3*R2*R4 + b1*b4*R2*R3 + b2*b3*R1*R4 + b2*b4*R1*R3 + b3*b4*R1*R2),
		R1*R2*R3*R4*(b1/R1 + b2/R2 + b3/R3 + b4/R4),
		- R1*R2*R3*R4)
	p = roots(coeff);
	if(debug) print(p);
	len = length(p)
	Xij = sapply(p, function(p) R - b*p);
	X13sq = Xij[1,]*Xij[3,] / p;
	# X24sq = p^2 / X13sq;
	X13 = sqrt(X13sq + 0i); X13 = c(X13, -X13);
	p = c(p, p); len = length(p); Xij = cbind(Xij, Xij);
	X24 = p / X13;
	# TODO: 10 roots per p;
	x1 = rootn(Xij[1,]^3 / Xij[2,] * X13, 10);
	# NOT robust!
	# sapply(seq(len),
		# function(id) rootn(Xij[1,id]^27 * Xij[3,id]^3 / Xij[2,id]^9 / Xij[4,id], 80));
	x2 = Xij[1,] / x1^3; x3 = Xij[2,] / x2^3; x4 = Xij[3,] / x3^3;
	sol = cbind(x1=as.vector(x1), x2=as.vector(x2), x3=as.vector(x3), x4=as.vector(x4))
	invisible(sol);
}

### Examples:
R = c(1,2,3,4)
b = c(1,2,-2, -1/4)
#
sol = solve.Pr2.S4P31(R, b)
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];

### Test
# TODO: debug vs robust ???
x1^3*x2 + b[1]*x1*x2*x3*x4 # - R1
x2^3*x3 + b[2]*x1*x2*x3*x4 # - R2
x3^3*x4 + b[3]*x1*x2*x3*x4 # - R3
x4^3*x1 + b[4]*x1*x2*x3*x4 # - R4


###############################

### Variant

# x1^2*x2*x3 + b1*x1*x2*x3*x4 = R1
# x2^2*x3*x4 + b2*x1*x2*x3*x4 = R2
# x3^2*x4*x1 + b3*x1*x2*x3*x4 = R3
# x4^2*x1*x2 + b4*x1*x2*x3*x4 = R4

### Solution:
# Step 1: same as above;


### Solver
solve.Pr3.S4P211 = function(R, b, debug=TRUE) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	b1 = b[1]; b2 = b[2]; b3 = b[3]; b4 = b[4];
	coeff = c(1- b1*b2*b3*b4,
		b1*b2*b3*b4*(R1/b1 + R2/b2 + R3/b3 + R4/b4),
		- (b1*b2*R3*R4 + b1*b3*R2*R4 + b1*b4*R2*R3 + b2*b3*R1*R4 + b2*b4*R1*R3 + b3*b4*R1*R2),
		R1*R2*R3*R4*(b1/R1 + b2/R2 + b3/R3 + b4/R4),
		- R1*R2*R3*R4)
	p = roots(coeff);
	if(debug) print(p);
	len = length(p)
	Xij = sapply(p, function(p) R - b*p);
	X13sq = Xij[1,]*Xij[3,] / p;
	# X24sq = p^2 / X13sq;
	X13 = sqrt(X13sq + 0i); X13 = c(X13, -X13);
	p = c(p, p); len = length(p); Xij = cbind(Xij, Xij);
	X24 = p / X13; X23 = Xij[2,] / X24; X12 = Xij[1,] / X13;
	x1 = sqrt(Xij[1,] / X23);
	x1 = c(x1, -x1);
	x2 = c(X12, X12) / x1; x3 = c(X13, X13) / x1; x4 = c(X24, X24) / x2;
	sol = cbind(x1=as.vector(x1), x2=as.vector(x2), x3=as.vector(x3), x4=as.vector(x4))
	invisible(sol);
}

### Examples:
R = c(1,2,3,4)
b = c(1,2,-2, -1/4)
#
sol = solve.Pr3.S4P211(R, b)
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];

### Test
x1^2*x2*x3 + b[1]*x1*x2*x3*x4 # - R1
x2^2*x3*x4 + b[2]*x1*x2*x3*x4 # - R2
x3^2*x4*x1 + b[3]*x1*x2*x3*x4 # - R3
x4^2*x1*x2 + b[4]*x1*x2*x3*x4 # - R4


### Ex 2:
R = c(-1,2,3,4)
b = c(1,2,-2, -3)
#
sol = solve.Pr3.S4P211(R, b)
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];

### Test
x1^2*x2*x3 + b[1]*x1*x2*x3*x4 # - R1
x2^2*x3*x4 + b[2]*x1*x2*x3*x4 # - R2
x3^2*x4*x1 + b[3]*x1*x2*x3*x4 # - R3
x4^2*x1*x2 + b[4]*x1*x2*x3*x4 # - R4

# degenerate Polynomial
round0.p(poly.calc(x1)) * prod(b) * 11

