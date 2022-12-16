########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S5: Hetero-Symmetric
### Useful Formulas
###
### draft v.0.1c


### Formulas:

# - Formulas & derivations;
# - Useful for S5 Ht Systems;

# - Applicable for systems described in:
#   TODO


### Sections

### Basic:
# A.) E11a
# B.) Higher Powers


####################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


### Debug
x = sqrt(c(2,3,5,7,11));
x[1] = - x[1]; x[5] = - x[5];
x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4]; x5 = x[5];

### Notation:
S  = sum(x);
E5 = prod(x);

s1 = x1 + x3; s2 = x2 + x4;
p1 = x1 * x3; p2 = x2 * x4;
ps = s1 * s2; sp = p1 + p2;
# E2 = x1*(S - x1) + x2*(x3 + x4 + x5) + x3*(x4 + x5) + x4*x5;
E2 = sp + ps + x5*(S - x5);
# E3 = x1*x2*(x3 + x4 + x5) + x3*x4*(x1 + x2 + x5) + (x1 + x2)*(x3 + x4)*x5;
E3 = p1*s2 + p2*s1 + x5*(sp + ps);
# E4 = x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x2*x4*x5 + x1*x3*x4*x5 + x2*x3*x4*x5;
E4 = p1*p2 + x5*(p1*s2 + p2*s1);

### E2:
E11a = x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1;
E11b = x1*x3 + x2*x4 + x3*x5 + x4*x1 + x5*x2;

E11a + E11b - E2 # = 0

### E3:
E111a = x1*x2*x3 + x2*x3*x4 + x3*x4*x5 + x4*x5*x1 + x5*x1*x2;
E111b = x1*x2*x4 + x2*x3*x5 + x3*x4*x1 + x4*x5*x2 + x5*x1*x3;

E111a + E111b - E3 # = 0

### Note:
# - the remaining cyclic permutations (E2) equal E11b & E11a;
#   Perm(S5, by = 3) = rev(E11b) = E11b;
#   Perm(S5, by = 4) = rev(E11a) = E11a;


### Mixed C2 & C3-Decomposition

s1 = x1 + x2; p1 = x1 * x2;
s2 = S - s1; e2 = (x3 + x4)*x5 + x3*x4; e3 = x3*x4*x5;

### Transformed P[5] System:
s1 + s2 - S # = 0
s1*s2 + p1 + e2 - E2 # = 0
s1*e2 + s2*p1 + e3 - E3 # = 0
s1*e3 + p1*e2 - E4 # = 0
p1*e3 - E5 # = 0

### Note:
# Order: 5! / (2*6) = 120/12 = 10;
# - much better than 5!, but still not enough;


#######################

A1 = E11a*E111a + E11b*E111b;
B1 = E11b*E111a + E11a*E111b;

### Sum:
A1 + B1 - E2*E3 # = 0

### Mult:
A1 * B1 - E11a*E11b*(E111a^2 + E111b^2) - E111a*E111b*(E11a^2 + E11b^2) # = 0
A1 * B1 - E11a*E11b*(E3^2 - 2*E111a*E111b) - E111a*E111b*(E2^2 - 2*E11a*E11b) # = 0
A1 * B1 - (E11a*E11b*E3^2 + E111a*E111b*E2^2) + 4*E11a*E11b*E111a*E111b # = 0

# TODO:
# - unfortunately is not yet symmetric;
#   E4321 = 60 (not 5! = 120);
#   E43111 = 5*4; OK
#   E42211 = 20 w coeff=3; 5 * choose(4,2) = 30;
#   E42211 = 10 w coeff=2; 5 * choose(4,2) = 30;
E11a*E11b*E111a*E111b


pE2a = toPoly.pm("x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1")
pE2b = toPoly.pm("x1*x3 + x2*x4 + x3*x5 + x4*x1 + x5*x2")
pE3a = toPoly.pm("x1*x2*x3 + x2*x3*x4 + x3*x4*x5 + x4*x5*x1 + x5*x1*x2")
pE3b = toPoly.pm("x1*x2*x4 + x2*x3*x5 + x3*x4*x1 + x4*x5*x2 + x5*x1*x3")
pR = prod.pm(pE2a, pE2b, pE3a, pE3b)

pc = countMonoms(pR)
pc = sort.pm(pc, "V1")
pc

# TODO
A2 = E11a*E11b*E3^2 + E111a*E111b*E2^2;
B2 = E11a*E11b*E2^2 + E111a*E111b*E3^2;

### Sum =>
A2 + B2 - (E11a*E11b + E111a*E111b)*(E2^2 + E3^2) # = 0


#######################

### Disaster:
pE1 = toPoly.pm("x1 + x2 + x3 + x4 + x5 - 1"); # S = 1; !!!
pE2a = toPoly.pm("x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1"); # E11a = 0; !!!
pE2b = toPoly.pm("x1*x3 + x2*x4 + x3*x5 + x4*x1 + x5*x2"); # E11b = 0; !!!
# E3 = 0; E4 = 0; !
pE3 = toPoly.pm("x1*x2*(x3 + x4 + x5) + x3*x4*(x1 + x2 + x5) + (x1 + x2)*x4*x5 + x3*x5*(x1 + x2)");
pE4 = toPoly.pm("x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x2*x4*x5 + x1*x3*x4*x5 + x2*x3*x4*x5");

pR = solve.lpm(pE1, pE2a, pE2b, pE3, pE4, xn=c("x5", "x4", "x3", "x2"))


### Overflows!
pE1 = toPoly.pm("s1 + s2 + x5 - 1"); # S = 1; !!!
pE2 = toPoly.pm("p1 + p2 + s1*s2 + x5*(1 - x5)"); # E2 = 0; S = 1;
pE3 = toPoly.pm("p1*s2 + p2*s1 + x5*(s1*s2 + p1 + p2)");
pE4 = toPoly.pm("p1*p2 + x5*(p1*s2 + p2*s1)");
pE5 = toPoly.pm("p1*p2*x5 - 1"); # E5 = 1;

pR = solve.lpm(pE1, pE5, pE2, pE3, pE4, xn=c("x5", "s2", "s1", "p2"))
max(pR[[4]]$Rez$p1)


####################

### Standard P[5]

### quasi-redundant
# spp = sp + ps; pp4 = p1*p2;
# - Solvable, but just ordinary P[5];
pE2 = toPoly.pm("spp + x5*(S - x5) - E2");
pE3 = toPoly.pm("A + x5*spp - E3");
pE4 = toPoly.pm("pp4 + x5*A - E4");
pE5 = toPoly.pm("pp4*x5 - E5");

pE4 = toPoly.pm("x5^2*(E3 - x5*spp) + E5 - E4*x5")

pR = solve.pm(pE2, pE4, xn=c("spp"))
pR = pR$Rez;
print.pm(pR, lead="x5") # trivial: ordinary P[5];


#####################

### Test

# any permutation:
x = 2*cos(2*pi* c(1,5,3,4,2) /11)
# *OR*
m = unity(5, all=T)
k = rootn(2, 5); x = sapply(m, function(m) (k*m)^3 - (k*m)^2 + 3*k*m)
#
x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4]; x5 = x[5];
S = sum(x);
s1 = x1 + x2; p1 = x1 * x2;
s2 = S - s1; e2 = (x3 + x4)*x5 + x3*x4; e3 = x3*x4*x5;

x12 = roots(c(1, -s1, p1));
x35 = roots(c(1, -s2, e2, -e3));

# the same:
poly.calc(c(x12, x35))
poly.calc(x)


##########################
##########################

### E11a & E11b:
# - a great discovery in mathematics;
# - newest formula is in file: Poly.System.S5.Ht.Formulas.Derivation.R;
# TODO: work out E3, E4;

27*(E11a^7 + E11b^7)*E5^2 + 81*E11a*E11b*(E11a^5 + E11b^5)*E5^2 + 27*(E11a*E11b)^2*(E11a^3 + E11b^3)*E5^2 +
	- 10*(E11a*E11b)^3*E2*E5^2 - (E11a*E11b)^6 +
	# - 54*(E11a*E11b)^5*(E11a^1 + E11b^1)*S^2 + 19*5^3*(E11a^4 + E11b^4)*E5^3*S +
	# TODO: + ... +
	- 5^5*E5^4*(E11b^2 + 3*E11a*E11b + E11a^2) # = 0


### Solver:

# Coeffs for x2 & x3;
source("Poly.System.S5.Ht.Formulas.CoeffX.R")
# Coeffs for Characterstic polynomial:
source("Poly.System.S5.Ht.Formulas.Derivation.Coeffs.R")

# Note:
# - some R-values are still fixed!
# - the functions: f6(), ..., f0() are in file:
#   Poly.System.S5.Ht.Formulas.Derivation.Coeffs.R;
solve.S5HtMixed = function(R, debug=TRUE, all=FALSE) {
	coeff = coeff.S5HtMixed(R);
	E11b = roots(coeff);
	if(debug) print(E11b);
	S = R[1]; E11a = R[2]; E3 = R[3]; E5 = R[5];
	E2 = R[2] + E11b;
	x1 = sapply(E2, function(E2) roots(c(1, -S, E2, -E3, R[4], -E5)));
	x1 = as.vector(x1);
	return(solve.S5HtMixed.x2(x1, E11b, R, all=all));
}
coeff.S5HtMixed = function(R) {
	S = R[1]; E11a = R[2]; E3 = R[3]; E4 = R[4]; E5 = R[5];
	# Note: still assumes E4 == 0 || E11a == 0;
	if(E4 != 0 && E11a != 0) {
		warning("The results are NOT exact!");
		# but may be still useful to indicate where the true roots are;
	}
	# Note: f[i]() are divided by E5^2;
	coeff = c(27, f6(R), f5(R), f4(R), f3(R), f2(R), f1(R), f0(R)
	);
	return(coeff);
}
# Classic Solver:
solve.S5HtMixed.Classic = function() {
	coeff = c(27, 0, 109, 0, 114, -189, -110, -654, -355, -570, 3332, 440, -1609, 1065, -1984, -6475, -660,
		1064, -1065, -1140, 3710, 440, 1635, 355, 570, -567, -110, -654, 0, -114, 189, 0, 109, 0, 0, -27);
	x1 = roots(coeff);
	# TODO:
	return(x1)
}
###

solve.S5HtMixed.x2 = function(x1, E11b, R, all=FALSE) {
	S = R[1]; E11a = R[2]; E3 = R[3]; E4 = R[4]; E5 = R[5];
	E11b = rep(E11b, 5);
	E2 = E11a + E11b;
	s = S - x1; e2 = E2 - s*x1; e3 = E3 - e2*x1; e4 = E4 - e3*x1;
	# Robust:
	# TODO
	return(x1); # temporary;
	### x2:
	# function is defined in file:
	# Poly.System.S5.Ht.Formulas.CoeffX.R;
	xCoeff = coeff.S5HtMixed.x2(x1, cbind(s, e2, e3, e4), E11a, E11b);
	v = xCoeff$v; dq = xCoeff$dq;
	v0r = v[,1]; v1r = v[,2]; v2r = v[,3]; v3r = v[,4];
	dq0r = dq[,1]; dq1r = dq[,2]; dq2r = dq[,3]; dq3r = dq[,4];
	# TODO: x2;
	# ...
	### x3:
	x3 = - (v3r*x2^3 + v2r*x2^2 + v1r*x2 + v0r) / (dq3r*x2^3 + dq2r*x2^2 + dq1r*x2 + dq0r);
	### x4:
	x4 = x2^2 - s*x2 + 2*x3*x2 - s*x3 + x3^2 + E11b - x3*x1;
	x4 = x4 / (x1 - x3);
	x5 = s - x2 - x3 - x4;
	sol = cbind(x1, x2, x3, x4, x5);
	return(sol);
}

### Examples:



