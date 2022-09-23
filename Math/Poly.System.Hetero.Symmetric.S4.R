########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
###
### draft v.0.4a-ref3



# "The many moods of an Irish setter"
#  and the many variants of S4!



### V1: x1^n + b*x2 = R
### V2a: x1^n + b*x1*x2 = R
### V2b: x1^n + b*x2*x3 = R
### V3a: x1^n + b*x1*x2*x3 = R
### V3b: x1^n + b*x2*x3*x4 = R
### V4: x1^n + b*x1*x2*x3*x4 = R
### ...

### TODO:
# - some proper classification;


### Note: C2-Decomposition
# - it is possible to use the C2-Decomposition
#   (see Poly.System.S4.C2.Formulas.R),
#   but it tends to inflate the formulas;
# - E3 = has a quadratic formula in C2,
#   and all derived formulas include this additional
#   quadratic step;

### Note: Refactor
# - Types 2a & 2b moved to separate file:
#   Poly.System.Hetero.Symmetric.S4.L1NLm.R;
# - Types 3a & 3b moved to separate file:
#   Poly.System.Hetero.Symmetric.S4.L1V3.R;
# - Type Asymmetric (Product-type) moved to:
#   Poly.System.Asymmetric.S5.Prod.R;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")

# library(polynom)
# library(pracma)

# the functions are in the file:
# Polynomials.Helper.R
# - e.g. round0(), round0.p(),
#   solve.EnAll(), solveEn();

test.S4.Simple = function(sol, R=NULL, b, n=2) {
	x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];
	err1 = x1^n + b*x2 # - R
	err2 = x2^n + b*x3 # - R
	err3 = x3^n + b*x4 # - R
	err4 = x4^n + b*x1 # - R
	err = rbind(err1, err2, err3, err4);
	if( ! is.null(R)) err = err - R;
	err = round0(err);
	return(err);
}


###########################
###########################
###########################

##############
### Simple ###
##############

### x[i]^n + b*x[i+1] = R

# x1^n + b*x2 = R
# x2^n + b*x3 = R
# x3^n + b*x4 = R
# x4^n + b*x1 = R

### Solution:

### Case 1: x1=x2=x3=x4
# P[n]: n solutions;

### Case 2: x1=x3, x2=x4, x1 != x3
# P[2] o P[(n^2 - n)/2];

### Case 3: all distinct;
# P[4] o P[(n^4 - n^2)/4];

# - System is decomposable into 3 subsystems:
#   P[n] * (P[2] o P[(n^2 - n)/2]) * (P[4] o P[(n^4 - n^2)/4]);


###############
### Order 2 ###
###############

### x[i]^2 + b*x[i+1] = R

### Solution:

### Case 1: x1=x2=x3=x4
# - 2 solutions;

### Case 2: x1=x3, x2=x4
# - 4 solutions: 2 overlap Case 1;

### Case 3: all distinct;
# Classic Poly: P[16 - 4] = P[12];
# S: P[3];

### Derivation:
# - see file: Poly.System.Hetero.Symmetric.S4.Derivation.R;

### Sum =>
S^2 - 2*E2 + b*S - 4*R # = 0
### Diff =>
(-E3^2 + E3*E2*S - E4*S^2) + b^6 # = 0
### Eq 3:
E4^2 - b^4*E4 + b*R^3*S - b^2*R^2*E2 + b^3*R*E3 - R^4 # = 0
### Eq 4:
E4^2 - b^4*E4 + 2*R^2*E4 + 2*R*E2*E4 - R*E3^2 - 2*R^2*E3*S +
	+ 2*R^3*E2 + R^2*E2^2 - R^3*S^2 + R^4 # = 0

### Eq:
S^3 + (3*b^2 - 4*R)*S - 4*b^3 # = 0

### Solver:
solve.Simple.S4P2 = function(R, b, debug=TRUE) {
	coeff = c(1, 0, 3*b[1]^2 - 4*R, - 4*b[1]^3);
	S = roots(coeff);
	if(debug) print(S);
	E2 = (S^2 + b*S - 4*R) / 2;
	E3 = E3.helper(S, R, b);
	E4 = (-E3^2 + E3*E2*S + b^6) / S^2;
	#
	len = length(S);
	x1 = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id], E4[id])));
	x1 = as.vector(x1);
	x2 = xi.f(x1, R, b, n=2);
	x3 = xi.f(x2, R, b, n=2);
	x4 = xi.f(x3, R, b, n=2);
	sol = cbind(x1=x1, x2=x2, x3=x3, x4=x4);
	id = order(abs(x1), abs(Re(x1)));
	sol = sol[id,];
	return(sol);
}
xi.f = function(x, R, b, n=2) {
	(R - x^n) / b[1];
}
E3.helper = function(S, R, b) {
	pE3 = ((87*b^16 + 57*R*b^14 + 2987*R^2*b^12 - 14958*R^3*b^10 + 16796*R^4*b^8 + 1052*R^5*b^6 +
			- 11312*R^6*b^4 + 6464*R^7*b^2 - 1280*R^8)*S^2 +
		(727*R*b^15 - 9711*R^2*b^13 + 11120*R^3*b^11 + 30776*R^4*b^9 - 68788*R^5*b^7 + 49168*R^6*b^5 +
			- 13760*R^7*b^3 + 768*R^8*b + 111*b^17)*S +
		(- 948*R*b^16 + 6132*R^2*b^14 + 6224*R^3*b^12 - 35360*R^4*b^10 + 34736*R^5*b^8 - 11904*R^6*b^6 +
			+ 768*R^7*b^4 - 196*b^18));
	pDiv = ((- 28*b^13 - 262*R*b^11 + 2450*R^2*b^9 - 4616*R^3*b^7 + 2596*R^4*b^5 + 32*R^5*b^3 - 448*R^6*b)*S^2 +
		(114*b^14 + 116*R*b^12 - 12*R^2*b^10 - 7044*R^3*b^8 + 16236*R^4*b^6 - 14768*R^5*b^4 + 6208*R^6*b^2 +
			- 1280*R^7)*S +
		(80*R*b^13 - 2256*R^2*b^11 + 8432*R^3*b^9 - 10384*R^4*b^7 + 5248*R^5*b^5 - 1280*R^6*b^3 - 88*b^15));
	return(- pE3/pDiv);
}

### Examples:

R = -1
b = 3
sol = solve.Simple.S4P2(R, b);
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];

test.S4.Simple(sol, b=b, n=2)


### Test
x1^2 + b*x2 # - R
x2^2 + b*x3 # - R
x3^2 + b*x4 # - R
x4^2 + b*x1 # - R


######################
######################

###############
### Order 3 ###
###############

### x[i]^3 + b*x[i+1] = R

### Solution:

### Case 1: x1=x2=x3=x4
# - 3 solutions;

### Case 2: x1=x3, x2=x4
# - 9 solutions: 3 overlap Case 1;

### Case 3: all distinct;
# Classic Poly: P[81 - 9] = P[72];
# S: P[18];

### Eq 1:
### Sum =>
S^3 - 3*E2*S + 3*E3 + b*S - 4*R # = 0

### Eq 2:
### Sum(x2*...) =>
E31a + b*(S^2 - 2*E2) - R*S # = 0
### Sum(x1^3*...) =>
S6 + b*E31a - R*S3 # = 0
### Diff: Eq2b - b*Eq2a =>
S6 - R*S3 - b^2*(S^2 - 2*E2) + b*R*S # = 0
# =>
S^6 - 6*E2*S^4 + 6*E3*S^3 - R*S^3 + 9*E2^2*S^2 - 6*E4*S^2 - b^2*S^2 +
	- 12*E2*E3*S + 3*R*E2*S - 2*E2^3 + 3*E3^2 + 6*E2*E4 +
	- 3*R*E3 + 2*b^2*E2 + b*R*S # = 0
# alternative: div by Eq 1 =>
3*E4*S^2 + 3*R*b*S - 3*E2*E3*S - 6*R^2 - b^2*E2 + E2^3 + 3*E3^2 - 3*E2*E4 # = 0

### Eq 3:
### Sum(x1*...) =>
S4 + b*E11a - R*S # = 0
### Sum(x2^3*...) =>
E33a + b*S4 - R*S3 # = 0
E33 - E33b + b*S4 - R*(S^3 - 3*E2*S + 3*E3) # = 0
E2^3 + 3*E3^2 - 3*E2*E3*S + 3*E4*S^2 - 3*E2*E4 - E11b^3 + 3*E4*E11b +
	+ b*S4 - R*(S^3 - 3*E2*S + 3*E3) # = 0
b^3*E2^3 + 3*b^3*E3^2 - 3*b^3*E2*E3*S - (b*E2 + S4 - R*S)^3 + 3*b^2*E4*(S4 - R*S) +
	+ 3*b^3*E4*S^2 + b^4*(S^4 - 4*S^2*E2 + 2*E2^2 + 4*S*E3 - 4*E4) +
	- b^3*R*(S^3 - 3*E2*S + 3*E3) # = 0
# 92 monomials

# [with some reductions]
54*E3*E2*S^7 + 27*E2^3*S^6 + 24*E2^2*b*S^6 + 10*E2*b^2*S^6 + b^3*S^6 - 33*E3^2*S^6 +
	- 9*E3*R*S^6 + 17*R^2*S^6 - 84*E2*E4*S^6 - 159*E3*E2^2*S^5 - 21*E2^2*R*S^5 +
	+ 21*E3*E2*b*S^5 - 60*E2*R*b*S^5 - 3*E3*b^2*S^5 + 3*R*b^2*S^5 + 96*E3*E4*S^5 +
	- 24*R*E4*S^5 - 108*E2^4*S^4 - 60*E2^3*b*S^4 - 3*E2^2*b^2*S^4 + b^4*S^4 +
	+ 192*E3^2*E2*S^4 - 123*E3*E2*R*S^4 + 48*E2*R^2*S^4 - 18*E3^2*b*S^4 +
	+ 45*E3*R*b*S^4 - 28*R^2*b*S^4 + 240*E2^2*E4*S^4 + 24*E2*b*E4*S^4 + 3*b^2*E4*S^4 +
	+ 192*E3*E2^3*S^3 - 48*E2^3*R*S^3 + 96*E3*E2^2*b*S^3 - 24*E2^2*R*b*S^3 - R*b^3*S^3 +
	- 64*E3^3*S^3 + 48*E3^2*R*S^3 - 12*E3*R^2*S^3 + R^3*S^3 - 432*E3*E2*E4*S^3 +
	+ 96*E2*R*E4*S^3 + 48*R*b*E4*S^3 + 48*E2^5*S^2 + 48*E2^4*b*S^2 - 4*E2*b^4*S^2 +
	+ 12*E2^3*b^2*S^2 - 176*E2^3*E4*S^2 + 3*b^3*E4*S^2 - 96*E3^2*E2^2*S^2 +
	+ 48*E3*E2^2*R*S^2 - 6*E2^2*R^2*S^2 - 48*E3^2*E2*b*S^2 + 24*E3*E2*R*b*S^2 +
	- 3*E2*R^2*b*S^2 - 96*E2^2*b*E4*S^2 - 28*E2*b^2*E4*S^2 + 240*E3^2*E4*S^2 +
	- 96*E3*R*E4*S^2 - 84*R^2*E4*S^2 + 144*E2*E4^2*S^2 - 48*E3*E2^4*S + 12*E2^4*R*S +
	+ 4*E3*b^4*S - 48*E3*E2^3*b*S + 12*E2^3*R*b*S - 3*E3*E2*b^3*S + 3*E2*R*b^3*S +
	- 12*E3*E2^2*b^2*S + 3*E2^2*R*b^2*S + 192*E3*E2^2*E4*S - 48*E2^2*R*E4*S +
	+ 12*E3*b^2*E4*S - 3*R*b^2*E4*S + 96*E3*E2*b*E4*S - 24*E2*R*b*E4*S - 192*E3*E4^2*S +
	+ 48*R*E4^2*S - 8*E2^6 - 12*E2^5*b - 6*E2^4*b^2 + 2*E2^2*b^4 + 48*E2^4*E4 - 4*b^4*E4 +
	+ 3*E3^2*b^3 - 3*E3*R*b^3 + 48*E2^3*b*E4 + 18*E2^2*b^2*E4 - 96*E2^2*E4^2 +
	- 48*E2*b*E4^2 - 12*b^2*E4^2 + 64*E4^3 # = 0

### Eq 4:
### Sum(x1*...) =>
S4 + b*E11a - R*S # = 0
### Sum(x3^3*...) =>
2*E33b + b*E13a - R*S3 # = 0
### Sum(x4*...) =>
E13a + 2*b*E11b - R*S # = 0

# Eq 4d: Eq 4b - b*Eq 4c =>
2*E33b - R*(S^3 - 3*E2*S + 3*E3) - 2*b^2*E11b + b*R*S # = 0
2*((E2 - E11a)^3 - 3*E4*(E2 - E11a)) - R*(S^3 - 3*E2*S + 3*E3) - 2*b^2*(E2 - E11a) + b*R*S # = 0
2*(b*E2 + S4 - R*S)^3 - 6*b^2*E4*(b*E2 + S4 - R*S) +
	+ 2*b^4*(R*S - S4) - b^3*R*(S^3 - 3*E2*S + 3*E3) - 2*b^5*E2 + b^4*R*S # = 0

# TODO

### Alternative Eq:
# - works only with the Simple system;
#  (with no additional terms)
### Diff: (Eq 1 - Eq 3) & (Eq 2 - Eq 4) =>
(x1 - x3)*(s1^2 - p1) + b*(x2 - x4) # = 0
(x2 - x4)*(s2^2 - p2) + b*(x3 - x1) # = 0
# Case: x[i] distinct
# Prod =>
(s1^2 - p1)*(s2^2 - p2) + b^2 # = 0
ps^2 - (p1*s2^2 + p2*s1^2) + E4 + b^2 # = 0
ps^2 - (S*(p1*s2 + p2*s1) - ps*sp) + E4 + b^2 # = 0
ps^2 + ps*sp - E3*S + E4 + b^2 # = 0
E11a*E2 - E3*S + E4 + b^2 # = 0
(R*S - S4)*E2 - b*E3*S + b*E4 + b^3 # = 0
# =>
E2*S^4 - 4*E2^2*S^2 + 2*E2^3 + 4*E2*E3*S - 4*E2*E4 - R*E2*S + b*E3*S - b*E4 - b^3 # = 0
# Reduction =>
E2^2*S^2 + b*E2*S^2 - 2*E2^3 - E2*E3*S + 4*E2*E4 - 3*R*E2*S - b*E3*S + b*E4 + b^3 # = 0

### S4:
# S4 = - E2*S^2 - b*S^2 + 2*E2^2 + E3*S - 4*E4 + 4*R*S;
# partly "reduced" version; but has 6 monomials!

### Eq S:
S^18 - 15*R*S^15 + 48*b^2*S^14 - 126*R*b*S^13 + (222*R^2 - 256*b^3)*S^12 + 609*R*b^2*S^11 +
	- (1764*R^2*b + 540*b^4)*S^10 + (2158*R^3 + 5061*R*b^3)*S^9 +
	- (8433*R^2*b^2 + 960*b^5)*S^8 + (6048*R^3*b - 1170*R*b^4)*S^7 +
	- (7671*R^4 - 3435*R^2*b^3 - 5800*b^6)*S^6 + (6099*R^3*b^2 - 18840*R*b^5)*S^5 +
	- (6300*R^4*b - 16632*R^2*b^4 + 3600*b^7)*S^4  + (8049*R^5 - 23297*R^3*b^3 + 10080*R*b^6)*S^3 +
	+ (1677*R^4*b^2 + 3672*R^2*b^5 - 10125*b^8)*S^2 +
	+ (2142*R^5*b - 7470*R^3*b^4 + 8100*R*b^7)*S +
	- 2744*R^6 + 9225*R^4*b^3 - 8910*R^2*b^6;


###########
### Solver:

solve.S4P3.Simple = function(R, b, debug=TRUE) {
	coeff = coeff.S4P3.Simple(R, b=b);
	S = roots(coeff);
	if(debug) print(S);
	#
	Ex = Ex.S4P3.Simple(S, R, b);
	E2 = Ex$E2; E3 = Ex$E3; E4 = Ex$E4;
	#
	len = length(S);
	x1 = sapply(seq(len), function(id) {
		roots(c(1, -S[id], E2[id], -E3[id], E4[id]));
	})
	x1 = as.vector(x1);
	#
	x2 = (R - x1^3) / b;
	x3 = (R - x2^3) / b;
	x4 = (R - x3^3) / b;
	#
	sol = cbind(x1, x2, x3, x4);
	return(sol);
}
coeff.S4P3.Simple = function(R, b) {
	coeffs = c(
		1, 0, 0, - 15*R, 48*b^2, - 126*R*b, (222*R^2 - 256*b^3), 609*R*b^2, # S^11
		- (1764*R^2*b + 540*b^4), (2158*R^3 + 5061*R*b^3), # S^9
		- (8433*R^2*b^2 + 960*b^5), (6048*R^3*b - 1170*R*b^4), # S^7
		- (7671*R^4 - 3435*R^2*b^3 - 5800*b^6), (6099*R^3*b^2 - 18840*R*b^5), # S^5
		- (6300*R^4*b - 16632*R^2*b^4 + 3600*b^7), (8049*R^5 - 23297*R^3*b^3 + 10080*R*b^6), # S^3
		(1677*R^4*b^2 + 3672*R^2*b^5 - 10125*b^8),
		(2142*R^5*b - 7470*R^3*b^4 + 8100*R*b^7),
		- 2744*R^6 + 9225*R^4*b^3 - 8910*R^2*b^6
	);
	return(coeffs);
}
Ex.S4P3.Simple = function(S, R, b) {
	### Coefficients
	c1 = S^3 + b*S - 4*R;
	c20 = 3*R*b*S - 6*R^2;
	c2 = c1^2 + 3*c20;
	c3 = - 3*b^3 - S*b*c1;
	c4 = (R*c1^3 + 9*R^2*c1^2 + 27*R^3*c1 - 27*S^3*R^3 + 27*R^4) / 9;
	c5 = - b*c2 - 3*S^2*c3;
	c33 = 18*S^2 + 3*b;
	c32 = 27*S*R + 12*b^2 + 9*S*c1;
	c31 = - 9*b^3 - 27*S^3*R - 4*c2 + 3*b^3 + 3*S^3*c1;
	#
	c410 = 3*(3*R^2*S^2 - b^4);
	c411 = 3*(c1*R + 3*R^2);
	c412 = 9*R*S;
	c403 = R^2 - R*S^3;
	c401 = 3*R^2*S*c1 + R*S*c1^2;
	#
	c57 = - 324*R*S - 108*S*c1 + 288*c412;
	c56 = - 324*R*S*b + 108*c3 + 1728*c403 - 288*c411 + 144*b*c412;
	c55 = - 2430*R^2*S^2 + 1836*R*S^2*c1 + 18*S^2*c1^2 + 1296*b*c403 + 288*c410 - 144*b*c411 +
		+ 432*R*S*c412 + 18*b^2*c412 - 48*S*c1*c412;
	c54 = - 972*R^2*S^2*b + 1404*R*S^2*b*c1 - 108*R*S*c3 - 36*S*c1*c3 - 576*c401 +
		+ 324*b^2*c403 + 144*b*c410 - 432*R*S*c411 - 18*b^2*c411 + 48*S*c1*c411 +
		+ 216*R*S*b*c412 - 24*S*b*c1*c412 + 48*c3*c412;
	c53 = - 2187*R^3*S^3 + 405*R^2*S^3*c1 + 324*R*S^2*b^2*c1 - 9*R*S^3*c1^2 - S^3*c1^3 +
		- 108*R*S*b*c3 + 18*c3^2 + 576*c4 - 432*b*c401 + 27*b^3*c403 + 432*R*S*c410 +
		+ 18*b^2*c410 - 48*S*c1*c410 - 216*R*S*b*c411 + 24*S*b*c1*c411 - 48*c3*c411 +
		+ 27*R*S*b^2*c412 - 3*S*b^2*c1*c412 + 24*b*c3*c412;
	c52 = - 729*R^3*S^3*b + 162*R^2*S^3*b*c1 + 27*R*S^2*b^3*c1 - 9*R*S^3*b*c1^2 - 405*R^2*S^2*c3 +
		+ 18*R*S^2*c1*c3 + 3*S^2*c1^2*c3 + 432*b*c4 - 108*b^2*c401 + 216*R*S*b*c410 +
		- 24*S*b*c1*c410 + 48*c3*c410 - 27*R*S*b^2*c411 + 3*S*b^2*c1*c411 - 24*b*c3*c411 + 3*b^2*c3*c412;
	c51 = - 162*R^2*S^2*b*c3 + 18*R*S^2*b*c1*c3 - 9*R*S*c3^2 - 3*S*c1*c3^2 + 108*b^2*c4 +
		- 9*b^3*c401 + 27*R*S*b^2*c410 - 3*S*b^2*c1*c410 + 24*b*c3*c410 - 3*b^2*c3*c411;
	c50 = - 9*R*S*b*c3^2 + c3^3 + 9*b^3*c4 + 3*b^2*c3*c410;
	#
	c57 = c57 - 36*c32 + 6*c33^2;
	c56 = c56 - 36*c31 - 6*c33*c32 + c57*c33/6;
	c55 = c55 - 36*c5 - 6*c33*c31 - c57*c32/6 + c56*c33/6;
	c54 = c54 - 6*c33*c5 - c57*c31/6 - c56*c32/6 + c55*c33/6;
	c53 = c53 - c57*c5/6 - c56*c31/6 - c55*c32/6 + c54*c33/6;
	c52 = c52 - c56*c5/6 - c55*c31/6 - c54*c32/6;
	c51 = c51 - c55*c5/6 - c54*c31/6;
	c50 = c50 - c54*c5/6;
	### E2:
	E2div = 36*c51^3 - 72*c50*c51*c52 + 36*c50^2*c53 - 6*c52^2*c53*c5 + 6*c51*c53^2*c5 +
		- 6*c52^3*c31 + 18*c51*c52*c53*c31 - 12*c50*c53^2*c31 + c53^3*c31^2 + 6*c51*c52^2*c32 +
		- 12*c51^2*c53*c32 - c53^3*c5*c32 - c52*c53^2*c31*c32 + c51*c53^2*c32^2 +
		+ 6*c51^2*c52*c33 - 6*c50*c52^2*c33 - 6*c50*c51*c53*c33 - c52*c53^2*c5*c33 +
		- c52^2*c53*c31*c33 + 2*c51*c53^2*c31*c33 + c51*c52*c53*c32*c33 - c50*c53^2*c32*c33 +
		+ c51^2*c53*c33^2 - c50*c52*c53*c33^2;
	#
	E2x0 = 36*c50*c51^2 - 36*c50^2*c52 - 6*c52^3*c5 + 12*c51*c52*c53*c5 - 6*c50*c53^2*c5 +
		+ 6*c50*c52*c53*c31 + c53^3*c5*c31 + 6*c50*c52^2*c32 - 12*c50*c51*c53*c32 +
		- c52*c53^2*c5*c32 + c50*c53^2*c32^2 + 6*c50*c51*c52*c33 - 6*c50^2*c53*c33 +
		- c52^2*c53*c5*c33 + c51*c53^2*c5*c33 + c50*c53^2*c31*c33 + c50*c52*c53*c32*c33 +
		+ c50*c51*c53*c33^2;
	E2 = -E2x0 / E2div;
	### E4
	E4 = (- 3*E2^3 + 3*(c1*S + b^2)*E2 - c2) / (9*(S^2 - E2));
	E3 = E2*S - c1/3;
	#
	return(list(E2=E2, E3=E3, E4=E4));
}
test.S4P3.Simple = function(sol, b, R=NULL) {
	test.S4.Simple(sol, R=R, b=b, n=3);
}

# TODO:
# - debug wrong/FALSE roots: 24 pairs?
# - Diff eqs: some of the roots may be wrong:
#   Case: x1 == x3 or x2 == x4; may need a separate solver;

### Examples:

### Ex 1:
R = -2; b = 3;
sol = solve.S4P3.Simple(R, b)

test.S4P3.Simple(sol, b=b)

### Ex 2:
R = 3; b = 1;
sol = solve.S4P3.Simple(R, b)

test.S4P3.Simple(sol, b=b)


### Debug:
b = -2; R = 3;
x1 = -0.4111558165 - 1.4630834100i;
x2 = -0.2145644927 + 1.1949483185i;
x3 = -1.0453736160 - 0.7706148327i;
x4 = -1.1400069541 - 1.0343850276i;
x = c(x1, x2, x3, x4);
#
s1 = x1 + x3; s2 = x2 + x4;
p1 = x1 * x3; p2 = x2 * x4;
ps = s1 * s2; sp = p1 + p2;
E11a = ps; E11b = sp;
#
S  = s1 + s2; E4 = p1 * p2;
E2 = sp + ps;
E3 = p1*s2 + p2*s1;
#
S4 = sum(x^4);


#####################
#####################
#####################

