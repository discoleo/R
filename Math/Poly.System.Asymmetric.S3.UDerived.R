########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Asymmetric Derived from Symmetric
###   Based on Roots of Unity
###
### draft v.0.1a


#######################

### Helper Functions

source("Polynomials.Helper.R")


#######################

### Base System:
# x^2 + y^2 + z^2 = R1
# x*y + x*z + y*z = R2
# x*y*z = R3

### Derived:
# x => x + y + z
# y => x + m*y + m^2*z
# z => x + m^2*y + m*z
# where m^3 = 1

### Derived System
# TODO


### Solver:
solve.S3.AsDer.P2 = function(R, b=0, debug=TRUE) {
	# Step 1:
	S = sqrt(R[1] + 2*R[2] + 0i);
	S = c(S, -S);
	if(debug) print(S);
	# Step 2:
	x3 = sapply(seq(2), function(id) {
		roots(c(1, -S[id], R[2], -R[3]));
	})
	S = rep(S, each=3);
	x3 = as.vector(x3);
	yz.s = S - x3;
	yz = R[2] - x3*yz.s;
	d3 = sqrt(yz.s^2 - 4*yz);
	y3 = (yz.s + d3)/2;
	z3 = (yz.s - d3)/2;
	# Step 3:
	m = unity(3, all=FALSE);
	mm = matrix(c(1,1,1, 1,m,m^2, 1,m^2,m), ncol=3, nrow=3, byrow=TRUE)
	mms = rbind(x3, y3, z3);
	sol = solve(mm, mms);
	return(sol);
}

### Examples:
R = c(-1,2,3)
sol = solve.S3.AsDer.P2(R)
m = unity(3, all=FALSE)
x = sol[1,]; y = sol[2,]; z = sol[3,];

### Test
(x + y + z)^2 + (x + m*y + m^2*z)^2 + (x + m^2*y + m*z)^2
(x + y + z)*(x + m*y + m^2*z) + (x + y + z)*(x + m^2*y + m*z) + (x + m*y + m^2*z)*(x + m^2*y + m*z)
(x + y + z)*(x + m*y + m^2*z)*(x + m^2*y + m*z)

### Classic Poly
round0.p(poly.calc(x) * 27)
