########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S2
### Asymmetric Derived from Symmetric
###   Based on Roots of Unity
###
### draft v.0.1a


#######################

### Helper Functions

source("Polynomials.Helper.R")


#######################

### Base System:
# x^3 + y^3 + b*(x+y) = R1
# x*y = R2

### Derived:
# x => x + m*y
# y => x + m^2*y
# where m^3 = 1

### Derived System
2*x^3 + 2*y^3 - 3*x*y*(x + y) + b*(2*x - y) - R[1] # = 0
x^2 + y^2 - x*y - R[2] # = 0


### Solver
solve.S2As.P3 = function(R, b, debug=TRUE) {
	m = unity(3, all=FALSE);
	# Step 1:
	# S^3 - (3*R2 - b)*S - R1 = 0
	S = roots(c(1,0, - 3*R[2] + b, -R[1]));
	if(debug) print(S);
	# Step 2:
	d = rootn(S^2 - 4*R[2], 2);
	d = c(d, -d); S = c(S, S);
	# Step 3:
	y = d / (m - m^2);
	x = (S + y) / 2;
	cbind(x=x, y=y);
}
test.S2As.P3 = function(sol, R=NULL, b) {
	err1 = 2*x^3 + 2*y^3 - 3*x*y*(x + y) + b*(2*x - y);
	err2 = x^2 + y^2 - x*y;
	err = rbind(err1, err2);
	if( ! is.null(R)) {
		err = err - R;
	}
	err = round0(err);
	return(err);
}

### Examples:
R = c(1,2)
b = 3
sol = solve.S2As.P3(R, b);
x = sol[,1]; y = sol[,2];

### Test
test.S2As.P3(sol, R, b)

### Classic Poly:
round0.p(poly.calc(x) * 27)


### Derivation:

### Eq 1:
(x+m*y)^3 + (x+m^2*y)^3 + b*(2*x + (m+m^2)*y) - R1 # = 0
2*x^3 + 2*y^3 + 3*x*y*(m*x + m^2*y + m^2*x + m*y) + b*(2*x - y) - R1 # = 0
2*x^3 + 2*y^3 - 3*x*y*(x + y) + b*(2*x - y) - R1 # = 0

### Eq 2:
(x + m*y)*(x + m^2*y) - R2 # = 0
x^2 + y^2 - x*y - R2 # = 0

