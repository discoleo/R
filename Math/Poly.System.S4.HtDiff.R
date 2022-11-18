########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Diff Hetero-Symmetric
###
### draft v.0.1b-special


### Hetero-Symmetric: Diff-Ht
# x1*x2 - x2*x3 + x3*x4 - x4*x1 = R2


####################

### Helper Functions

source("Polynomials.Helper.R")


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

# - based on C2-Transform;
# - see derivation in file:
#   Poly.System.S4.HtDiff.Derivation.R;

### Transformed System
ps^2 + 4*ps*sp - 4*E3*S + 16*E4 - E11d^2 # = 0
E3^2 - sp*S*E3 + ps*sp^2 + E4*S^2 - 4*ps*E4 # = 0


### Eq ps:
ps^4 - 2*E11d^2*ps^2 - 32*E4*ps^2 - 4*S*E3*ps^2 + 16*E4*S^2*ps + 16*E3^2*ps + E11d^4 - 32*E11d^2*E4 +
	+ 256*E4^2 + 4*E11d^2*S*E3 - 64*E4*S*E3 # = 0

### Eq for Special case:
64*ps^3 + 16*S^2*ps^2 - 4*S^4*ps - 4096*E4*ps + 1024*E3^2 - S^6 # = 0


### Solver:
solve.S4HtDiff.P1 = function(R, debug=TRUE, all=FALSE) {
	coeff = coeff.S4HtDiff.P1(R, debug=debug);
	# Special Cases:
	isSpecial = attr(coeff, "Special");
	if(is.null(isSpecial)) { isSpecial = FALSE; }
	else coeff = as.vector(coeff);
	ps = roots(coeff);
	if(debug) print(ps);
	S = R[1]; E11d = R[2]; E3 = R[3]; E4 = R[4];
	sp = - (ps^2 - E11d^2 + 16*E4 - 4*E3*S) / (4*ps);
	# Step 2:
	len = length(ps);
	s1 = sapply(seq(len), function(id) roots(c(1, -S, ps[id])));
	s1 = as.vector(s1);
	s2 = S - s1;
	sp = rep(sp, each=2);
	# Step 3: robust
	p1 = (sp*s1 - E3) / (s1 - s2);
	p2 = sp - p1;
	# Step 4:
	len = length(p1);
	x13 = sapply(seq(len), function(id) roots(c(1, - s1[id], p1[id])));
	x1  = x13[1,]; x3 = x13[2,];
	d24 = E11d / (x1 - x3);
	x2 = (s2 + d24)/2;
	x4 = (s2 - d24)/2;
	#
	sol = cbind(x1, x2, x3, x4);
	if(isSpecial) {
		sol2 = solve.S4HtDiff.P1.Special(R, debug=debug);
		sol  = rbind(sol, sol2);
	}
	if(all) {
		sol = rbind(sol, sol[, c(3,4,1,2)]);
	}
	return(sol);
}
solve.S4HtDiff.P1.Special = function(R, debug=TRUE) {
	S = R[1]; E11d = R[2]; E3 = R[3]; E4 = R[4];
	s1 = s2 = S/2;
	sp = E3 / s1; # TODO: S == 0;
	p1 = roots(c(1, -sp, E4));
	# p1 = p1[1]; # p2 = p1[2]; p1 = p1[1];
	# Step 4:
	len = length(p1)
	x13 = sapply(seq(len), function(id) roots(c(1, - s1, p1[id])));
	x1  = x13[1,]; x3 = x13[2,];
	d24 = E11d / (x1 - x3);
	x2 = (s2 + d24)/2;
	x4 = (s2 - d24)/2;
	#
	sol = cbind(x1, x2, x3, x4);
	return(sol);
}
coeff.S4HtDiff.P1 = function(R, debug=TRUE) {
	S = R[1]; E11d = R[2]; E3 = R[3]; E4 = R[4];
	# Special Cases:
	z = S^4 - 32*E3*S + 256*E4 - 16*E11d^2;
	isSpecial = (round0(z) == 0);
	if(isSpecial) {
		if(debug) warning("Special Case!");
		coeff = c(64, 16*S^2, - 4*S^4 - 4096*E4, 1024*E3^2 - S^6);
		attr(coeff, "Special") = TRUE;
		return(coeff);
	}
	coeff = c(1, 0, - 2*E11d^2 - 32*E4 - 4*S*E3, 16*E4*S^2 + 16*E3^2,
		E11d^4 - 32*E11d^2*E4 + 256*E4^2 + 4*E11d^2*S*E3 - 64*E4*S*E3);
	return(coeff);
}

### Examples:
R = c(5,3,-1,2)
x = solve.S4HtDiff.P1(R)

test.S4HtDiff.P1(x)

# Leading: (... + 16*E4);
print(poly.calc(x[,c(1,3)]) * 23, 16)


### Ex 2:
R = c(5,3,-2,2)
x = solve.S4HtDiff.P1(R)

test.S4HtDiff.P1(x)

print(poly.calc(x[,c(1,3)]) * 23, 16)


### Ex 3:
R = c(4,2,3,-1)
x = solve.S4HtDiff.P1(R)

test.S4HtDiff.P1(x)

print(poly.calc(x[,c(1,3)]) * 17, 16)


### Ex 4:
R = c(4,2,1,-3)
x = solve.S4HtDiff.P1(R)

test.S4HtDiff.P1(x)

print(poly.calc(x[,c(1,3)]) * 13, 16)


###  Special Cases

### Ex 5:
R = c(2,-1,4,1)
x = solve.S4HtDiff.P1(R)

test.S4HtDiff.P1(x)

### Ex 6:
R = c(2,3,2,1)
# TODO: 2 pairs of double roots OR still missing?
x = solve.S4HtDiff.P1(R)

test.S4HtDiff.P1(x)

