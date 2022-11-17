########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Diff Hetero-Symmetric
### Derivations
###
### draft v.0.1a


### Derivations
# - Basic derivations;
# - Numerical approaches: basic intuition;


####################

### Helper Functions

### Solver Tools
source("Polynomials.Helper.Solvers.Num.R")

# - is loaded automatically in "Solvers.Num.R";
# source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")

### Other
test.S4HtDiff.P1 = function(x, R=NULL) {
	if(is.null(dim(x))) {
		x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4];
	} else {
		x1 = x[,1]; x2 = x[,2]; x3 = x[,3]; x4 = x[,4];
	}
	s1 = x1 + x3; s2 = x2 + x4; S  = s1 + s2;
	p1 = x1 * x3; p2 = x2 * x4; E4 = p1 * p2;
	ps = s1 * s2; sp = p1 + p2;
	E3 = p1*s2 + p2*s1;

	### E2:
	E11d = x1*x2 - x2*x3 + x3*x4 - x4*x1;
	#
	err = rbind(S, E11d, E3, E4);
	if( ! is.null(R)) err = err - R;
	return(err);
}

####################
####################

### Base System
x1 + x2 + x3 + x4 # = R1
x1*x2 - x2*x3 + x3*x4 - x4*x1 # = R2
x1*x3*(x2 + x4) + x2*x4*(x1 + x3) # = R3
x1*x2*x3*x4 # = R4

### Solution:

### C2-Transform
# => Transformed System
S - R1 # = 0
(x1 - x3)*(x2 - x4) - R2 # = 0
p1*s2 + p2*s1 - R3 # = 0
E4 - R4 # = 0

### Eq 2:
(s1^2 - 4*p1)*(s2^2 - 4*p2) - R2^2 # = 0
ps^2 - 4*(p1*s2^2 + p2*s1^2) + 16*E4 - R2^2 # = 0

### E3*S
E3*S - (p1*s2^2 + p2*s1^2) - ps*sp # = 0
# =>
ps^2 - 4*(E3*S - ps*sp) + 16*E4 - R2^2 # = 0

### Eq E3:
# - see file: Poly.System.S4.C2.Formulas.R;
E3^2 - sp*S*E3 + ps*sp^2 + E4*S^2 - 4*ps*E4 # = 0


### Transformed System
ps^2 + 4*ps*sp - 4*E3*S + 16*E4 - E11d^2 # = 0
E3^2 - sp*S*E3 + ps*sp^2 + E4*S^2 - 4*ps*E4 # = 0


### Eq ps:
ps^4 - 2*E11d^2*ps^2 - 32*E4*ps^2 - 4*S*E3*ps^2 + 16*E4*S^2*ps + 16*E3^2*ps + E11d^4 - 32*E11d^2*E4 +
	+ 256*E4^2 + 4*E11d^2*S*E3 - 64*E4*S*E3 # = 0


### Derivation:

### Numerical Solution:

solve.S4HtDiff.Num = function(x, R) {
	x = matrix(x, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
	x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4];
	s1 = x1 + x3; s2 = x2 + x4; S  = s1 + s2;
	p1 = x1 * x3; p2 = x2 * x4; E4 = p1 * p2;
	ps = s1 * s2; sp = p1 + p2;
	E3 = p1*s2 + p2*s1;

	### E2:
	E11d = x1*x2 - x2*x3 + x3*x4 - x4*x1;
	#
	y = c(S, E11d, E3, E4) - R;
	y = rbind(Re(y), Im(y));
	return(y);
}


### Debug
R = c(5,3,-1,2)
x0 = c(4.5-1.4i, 0.65+0.95i, -0.25-0.4i, 0.06+0.8i)
x  = solve.all(solve.S4HtDiff.Num, x0, R=R, debug=T)
x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4];
#
R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
s1 = x1 + x3; s2 = x2 + x4; S  = s1 + s2;
p1 = x1 * x3; p2 = x2 * x4; E4 = p1 * p2;
ps = s1 * s2; sp = p1 + p2;
E3 = p1*s2 + p2*s1; E11d = (x1 - x3)*(x2 - x4);


test.S4HtDiff.P1(x)


x1 =  4.53993054-1.40670948i;
x2 =  0.65407124+0.95114529i;
x3 = -0.25001533-0.36559021i;
x4 =  0.05601355+0.82115439i;


### Derivation:
p1 = toPoly.pm("ps^2 + 4*ps*sp - 4*E3*S + 16*E4 - E11d^2");
p2 = toPoly.pm("E3^2 - sp*S*E3 + ps*sp^2 + E4*S^2 - 4*ps*E4");

pR = solve.pm(p1, p2, "sp")
str(pR)

print.pm(pR$Rez, "ps")

