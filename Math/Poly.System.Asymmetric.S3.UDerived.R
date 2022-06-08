########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Asymmetric Derived from Symmetric
###   Based on Roots of Unity
###
### draft v.0.1a-ext2


#######################

### Helper Functions

source("Polynomials.Helper.R")


#######################

### Base System:
# x^2 + y^2 + z^2 = R1
# x*y + x*z + y*z = R2
# x*y*z = R3

### Extension:
# x^2 + y^2 + z^2 + b1*(x+y+z) = R1
# x*y + x*z + y*z + b2*(x+y+z) = R2


### Derived:
# x => x + y + z
# y => x + m*y + m^2*z
# z => x + m^2*y + m*z
# where m^3 = 1

### Derived System
x^2 + 2*y*z + b1*x # = R1 / 3
x^2 - y*z + b2*x # = R2 / 3
x^3 + y^3 + z^3 - 3*x*y*z # = R3


### Solver:
solve.S3.AsDer.P2 = function(R, b=0, debug=TRUE) {
	# Step 1:
	if(all(b == 0)) {
		S = sqrt(R[1] + 2*R[2] + 0i);
		S = c(S, -S);
		E2 = c(R[2], R[2]);
	} else {
		S = roots(c(1, b[1] + 2*b[2], -2*R[2] - R[1]))
		E2 = R[2] - b[2]*S;
	}
	if(debug) print(S);
	# Step 2:
	x3 = sapply(seq(2), function(id) {
		roots(c(1, -S[id], E2[id], -R[3]));
	})
	S = rep(S, each=3); E2 = rep(E2, each=3);
	x3 = as.vector(x3);
	yz.s = S - x3;
	yz = E2 - x3*yz.s;
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
R = c(-3,2,3)
b = c(1, -2)
sol = solve.S3.AsDer.P2(R, b=b)
m = unity(3, all=FALSE)
x = sol[1,]; y = sol[2,]; z = sol[3,];

### Test
3*(x^2 + 2*y*z + b[1]*x) # = R[1]
3*(x^2 - y*z + b[2]*x) # = R[2]
x^3 + y^3 + z^3 - 3*x*y*z # = R[3]

### Classic Poly
round0.p(poly.calc(x) * 27 * 27)


### Derivation

### Eq 1:
p1 = toPoly.pm("(x + y + z)^2 + (x + m*y + m^2*z)^2 + (x + m^2*y + m*z)^2")
p1 = replace.pm(p1, 1, "m", pow=3)
p1 = replace.pm(p1, toPoly.pm("-m-1"), "m", pow=2)

3*(x^2 + 2*y*z)

### Eq 2:
p2 = toPoly.pm(paste0(
	"(x + y + z)*(x + m*y + m^2*z) + (x + y + z)*(x + m^2*y + m*z) +",
	"(x + m*y + m^2*z)*(x + m^2*y + m*z)"))
p2 = replace.pm(p2, 1, "m", pow=3)
p2 = replace.pm(p2, toPoly.pm("-m-1"), "m", pow=2)

3*(x^2 - y*z)

### Eq 3:
p3 = toPoly.pm("(x + y + z)*(x + m*y + m^2*z)*(x + m^2*y + m*z)")
p3 = replace.pm(p3, 1, "m", pow=3)
p3 = replace.pm(p3, toPoly.pm("-m-1"), "m", pow=2)

x^3 + y^3 + z^3 - 3*x*y*z

