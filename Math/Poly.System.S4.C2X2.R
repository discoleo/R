########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S4: Hetero-Symmetric
### Type: C2 X2
###
### draft v.0.1a


### S4 Systems
### Type: C2 + X2


### Symmetries

# - if (x1,x2,y1,y2) is a solution, then
#   the following permutations are also solutions:
#   C2 permutation: {(x2,x1}, (y2,y1)};
#   X2 permutation: {(y1,y2), (x2,x1)};
### Note:
# - in the X2 permutation: the x's are swapped (x2,x1);
# - NO solutions of type: x1 = y1, x2 = y2;
#   (only in special cases)


### Examples:

### A. Simple: s&s
# x1^p1 + x2^p1 = R1
# y1^p1 + y2^p1 = R1
# x1^n1*y1^n2 + x2^n1*y2^n2 = R2
# y1^n1*x2^n2 + y2^n1*x1^n2 = R2

### B. SP-Type
# x1^p1 + x2^p1 + y1^p1 + y2^p1 = R1
# x1^n1*y1^n2 + x2^n1*y2^n2 = R2
# y1^n1*x2^n2 + y2^n1*x1^n2 = R2
# x1*x2*y1*y2 = R3

### C. Non-Oriented
# x1^p1 + x2^p1 = R1
# y1^p1 + y2^p1 = R1
# x1^n1*y1^n2 + x2^n1*y2^n2 - (y1^n1*x2^n2 + y2^n1*x1^n2) = 0
# x1*x2*y1*y2 = R3


####################

### Helper Functions

### Solver Tools
source("Polynomials.Helper.Solvers.Num.R")

# - is loaded automatically in "Solvers.Num.R";
# source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


### Other

test.S4C2X2.Simple = function(x, R=NULL, n=c(2,1)) {
	if(is.null(dim(x))) {
		x1 = x[1]; x2 = x[2]; y1 = x[3]; y2 = x[4];
	} else {
		x1 = x[,1]; x2 = x[,2]; y1 = x[,3]; y2 = x[,4];
	}
	n1 = n[1]; n2 = n[2];
	err1 = x1 + x2;
	err2 = y1 + y2;
	err3 = x1^n1*y1 + x2^n1*y2;
	err4 = y1^n1*x2 + y2^n1*x1;
	err = rbind(err1, err2, err3, err4);
	if( ! is.null(R)) err = err - rep(R, each=2);
	return(err)
}


######################
######################

###################
### Simple Type ###
###################

### Order 2+1
# x1 + x2 = R1
# y1 + y2 = R1
# x1^2*y1 + x2^2*y2 = R2
# y1^2*x2 + y2^2*x1 = R2

### Solution:

###
E11a = x1*y1 + x2*y2;
E11b = x1*y2 + x2*y1;

### Eq 3 =>
E11a*s1 - p1*s2 - R2 # = 0
### Eq 4 =>
E11b*s2 - p2*s1 - R2 # = 0

### Sum:
E11a + E11b - s1*s2 # = 0
E11a + E11b - ps # = 0

### Prod:
E11a * E11b - p1*(s2^2 - 2*p2) - p2*(s1^2 - 2*p1) # = 0
E11a * E11b - (p1*s2^2 + p2*s1^2) + 4*E4 # = 0
E11a * E11b - (p1*s2 + p2*s1)*S + ps*sp + 4*E4 # = 0

### Transformed System:
E11a*s1 - p1*s2 - R2 # = 0
E11b*s2 - p2*s1 - R2 # = 0
E11a + E11b - ps # = 0
E11a * E11b - (p1*s2 + p2*s1)*S + ps*sp + 4*E4 # = 0
# with: sp = p1 + p2 and E4 still unknown;
# ps = R1^2 and S = 2*R1 are known;

### Sum Eqs 1b & 2b:
E11a*s1 + E11b*s2 - (p1*s2 + p2*s1) - 2*R2 # = 0
### Sum: s1*Eq1b + s2*Eq2b =>
E11a*s1^2 + E11b*s2^2 - ps*sp - R2*S # = 0
### Prod Eqs 1b & 2b:
(E11a*s1 - R2)*(E11b*s2 - R2) - ps*E4 # = 0
ps*E11a*E11b - R2*(E11a*s1 + E11b*s2) + R2^2 - ps*E4 # = 0

# Subst =>
E11a * E11b - (E11a*s1 + E11b*s2 - 2*R2)*S + ps*sp + 4*E4 # = 0
E11a * E11b - (E11a*s1 + E11b*s2)*S + E11a*s1^2 + E11b*s2^2 + R2*S + 4*E4 # = 0
E11a * E11b - ps^2 + R2*S + 4*E4 # = 0
5*ps*E11a * E11b - 4*R2*(E11a*s1 + E11b*s2) - ps^3 + R2*ps*S + 4*R2^2 # = 0


### Eq E11a:
5*ps*E11a^2 - (5*ps^2 + 4*R2*s2 - 4*R2*s1)*E11a + ps^3 - 4*R2^2 - R2*ps*S + 4*R2*ps*s2 # = 0
#
5*R1^2*E11a^2 - 5*R1^4*E11a + R1^6 + 2*R1^3*R2 - 4*R2^2 # = 0


### Solver:

solve.S4C2X2.E21 = function(R, debug=TRUE, all=TRUE) {
	R1 = R[1]; R2 = R[2];
	s1 = s2 = R1; S = s1 + s2; ps = s1 * s2;
	if(round0(s1) == 0) {
		if(round0(R2) != 0) {
			warning("No solution!");
			return(array(numeric(0), c(0,4)));
		}
		warning("Infinitely many solutions!");
		x = 3; y = -4; # just an example;
		x1 = c(x, NA); x2 = c(-x, NA);
		y1 = c(y, NA); y2 = c(-y, NA);
		return(cbind(x1=x1, x2=x2, y1=y1, y2=y2));
	}
	coeff = c(5*ps, - 5*ps^2 - 4*R2*s2 + 4*R2*s1,
		ps^3 - 4*R2^2 - ps*R2*S + 4*ps*R2*s2);
	E11a = roots(coeff);
	E11b = ps - E11a;
	if(debug) print(E11a);
	# Case: s1 = s2 = 0 (is already covered);
	p1 = (E11a*s1 - R2) / s2;
	p2 = (E11b*s2 - R2) / s1;
	#
	len = length(p1);
	x12 = sapply(seq(len), function(id) roots(c(1, -s1, p1[id])));
	x1 = x12[1,]; x2 = x12[2,];
	y1 = (s2*x1 - s1*s2 + E11a) / (2*x1 - s1);
	x2 = s1 - x1; y2 = s2 - y1;
	#
	sol = cbind(x1, x2, y1, y2);
	if(all) {
		sol = rbind(sol, sol[, c(2,1,4,3)]);
	}
	return(sol);
}
test.S4C2X2.E21 = function(sol, R=NULL) {
	test.S4C2X2.Simple(sol, n=c(2,1));
}

### Examples:

### Ex 1:
R = c(2, -3)
sol = solve.S4C2X2.E21(R)

test.S4C2X2.E21(sol)


### Ex 2:
R = c(-4, 1)
sol = solve.S4C2X2.E21(R)

test.S4C2X2.E21(sol)


### Debug

R = c(5, -3)
x1 = 2.5 + 0.850396226268i;
x2 = 2.5 - 0.850396226268i;
y1 = 2.5 + 3.60233622228i;
y2 = 2.5 - 3.60233622228i;
R1 = R[1]; R2 = R[2];
s1 = x1 + x2; s2 = y1 + y2;
p1 = x1 * x2; p2 = y1 * y2;
S  = s1 + s2; E4 = p1 * p2;
ps = s1 * s2; sp = p1 + p2;


### Numeric:
R = c(5, -3)
x0 = c(2.5+0.8504i, 2.5-0.8504i, 2.5+3.6023i, 2.5-3.6023i);
sol = solve.all(wrap(test.S4C2X2.Simple), x0=x0, R=R, atol=1E-12, rtol=1E-10)

### Derivation:

# robust: but NOT necessary;
E11a = x1*y1 + x2*y2;
E11b = x1*y2 + x2*y1;
x1*y1 + (s1 - x1)*(s2 - y1) - E11a # = 0
2*x1*y1 - s2*x1 - s1*y1 + s1*s2 - E11a # = 0
# =>
(s2^2 - 4*p2)*x1^2 - s1*(s2^2 - 4*p2)*x1 - E11a^2 + E11a*s2*s1 - s1^2*p2 # = 0
# =>
(s2^2 - 4*p2)*(s1*x1 - p1) - s1*(s2^2 - 4*p2)*x1 - E11a^2 + E11a*s2*s1 - s1^2*p2 # = 0
# redundant: x1 cancels out;
- p1*(s2^2 - 4*p2) - E11a^2 + E11a*s2*s1 - s1^2*p2 # = 0


###
p1y = toPoly.pm("2*x1*y1 - s2*x1 - s1*y1 + s1*s2 - E11a")
p2y = toPoly.pm("y1^2 - s2*y1 + p2")

# redundant:
x1*(s2 - y1) + (s1 - x1)*y1 - E11b # = 0
s2*x1 + y1*s1 - 2*x1*y1 - E11b # = 0

