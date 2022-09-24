########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
###  == Derivation ==
###  Type: L1 V3
###
### draft v.0.1c


####################
####################

### Helper Functions

source("Polynomials.Helper.R")


####################
####################

### History

# draft v.0.1a:
# - moved specific code from file:
#   Poly.System.Hetero.Symmetric.S4.Derivation.R; [draft v.0.3b]

####################
####################

################
### Type V3b ###
################

###############
### Order 2 ###
###############

### V3b:
### x1^2 + b*x2*x3*x4 = R
x1^2 + b*x2*x3*x4 # - R
x2^2 + b*x1*x3*x4 # - R
x3^2 + b*x1*x2*x4 # - R
x4^2 + b*x1*x2*x3 # - R

### Solution:

### Diff Eq[i] - Eq[i+1] =>
(x1 - x2)*(x1 + x2 - b*x3*x4) # = 0
(x2 - x3)*(x2 + x3 - b*x1*x4) # = 0
(x3 - x4)*(x3 + x4 - b*x1*x2) # = 0
(x4 - x1)*(x1 + x4 - b*x2*x3) # = 0
### Diff Eq[i] - Eq[i+2], etc. =>
(x1 - x3)*(x1 + x3 - b*x2*x4) # = 0
(x2 - x4)*(x2 + x4 - b*x1*x3) # = 0
# Case: x[i] != x[j]: Sum =>
3*S - b*E2 # = 0

### Sum =>
S^2 - 2*E2 + b*E3 - 4*R # = 0

### Sum(x1*...) =>
(x1^3 + x2^3 + x3^3 + x4^3) + 4*b*E4 - R*S # = 0
S^3 - 3*E2*S + 3*E3 + 4*b*E4 - R*S # = 0

### Diff(x[i]*Eq[i] - x[i+1]*Eq[i+1]) =>
(x1 - x2)*(x1^2 + x2^2 + x1*x2 - R) # = 0
(x2 - x3)*(x2^2 + x3^2 + x2*x3 - R) # = 0
# ...

### Sum(Diff 1)^2 vs Sum(Diff 2):
x1^2 + x2^2 + 2*x1*x2 - b^2*(x3*x4)^2 # = 0
3*S^2 - 6*E2 + 2*E2 - b^2*(E2^2 - 2*S*E3 + 2*E4) # = 0
### Sum(Diff 2):
3*S^2 - 5*E2 - 6*R # = 0
# =>
E2 - b^2*(E2^2 - 2*S*E3 + 2*E4) + 6*R = 0
# TODO:
# - verify!
# - analyze existence of solutions;
# - possibility to simplify solution using this Diff;


### =>
#   b*E2 = 3*S
# - b^2*E3 = b*S^2 - 6*S - 4*b*R
# - 4*b^3*E4 = b^2*S^3 - 12*b*S^2 - b^2*R*S + 18*S + 12*b*R

### Sum(x1^2*...) =>
(x1^4 + x2^4 + x3^4 + x4^4) + b*E4*S - R*(S^2 - 2*E2) # = 0
S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 - 4*E4 + b*E4*S - R*(S^2 - 2*E2) # = 0
S^4 - R*S^2 + 2*E2^2 - 4*E2*S^2 + 2*R*E2 + 4*E3*S + b*E4*S - 4*E4 # = 0
#
3*b^3*S^4 - 48*b^2*S^3 - 3*b^3*R*S^2 + 102*b*S^2 + 72*b^2*R*S + 72*S + 48*b*R # = 0


#########
### [old] [probably cyclic redundancy]

### Sum(x2*x3*x4*...) =>
# - possibly/probably redundant!
E4*S + b*Sum((x2*x3*x4)^2 ) - R*E3 # = 0
E4*S + b*(E3^2 - 2*E4*E2) - R*E3 # = 0
E4*S - 2*b*E2*E4 + b*E3^2 - R*E3 # = 0

### SEq1:
(10*R*S^2 - 14*R*S^3*b + 40*R^2*S*b - S^4 + S^5*b) +
(- 22*R*S*b^2 + 32*R*b + 6*S - 9*S^2*b + 4*S^3*b^2)*E3^1 +
(3*S*b^3 - 14*b^2)*E3^2

### SEq2:
(- 8*R*S^3 + 12*R*S^4*b - 8*R^2*S - 24*R^2*S^2*b - 32*R^3*b + S^5 - S^6*b) +
(8*R - 28*R*S*b + 10*R*S^2*b^2 + 32*R^2*b^2 - 8*S^2 + 10*S^3*b - 3*S^4*b^2)*E3^1 +
(- 10*R*b^3 + 9*S*b^2 - S^2*b^3 - 8*b)*E3^2 +
(b^4)*E3^3


### E3 =>
Subst = 560*R*S^2 - 780*R*S^3*b - 224*R*S^4*b^2 + 213*R*S^5*b^3 - 9*R*S^6*b^4 + 3024*R^2*S*b - 564*R^2*S^2*b^2 - 768*R^2*S^3*b^3 + 24*R^2*S^4*b^4 + 816*R^3*S*b^3 - 16*R^3*S^2*b^4 + 3136*R^3*b^2 - 56*S^4 + 36*S^5*b + 37*S^6*b^2 - 18*S^7*b^3 + S^8*b^4;
Subst = - Subst;

E3Div = - 324*R*S*b^2 - 256*R*S^2*b^3 + 249*R*S^3*b^4 - 5*R*S^4*b^5 + 1008*R*b - 252*R^2*S*b^4 + 4*R^2*S^2*b^5 - 1408*R^2*b^3 + 336*S - 188*S^2*b - 240*S^3*b^2 + 199*S^4*b^3 - 48*S^5*b^4 + S^6*b^5;

### =>
(10*R*S^2 - 14*R*S^3*b + 40*R^2*S*b - S^4 + S^5*b) +
(- 22*R*S*b^2 + 32*R*b + 6*S - 9*S^2*b + 4*S^3*b^2) * (Subst/E3Div) +
(3*S*b^3 - 14*b^2) * (Subst/E3Div)^2

### Eq:
b^5*S^10 +
	+ 2*b^4*S^9 - (2*b^3 + 6*R*b^5)*S^8  + (19*b^2 - 110*R*b^4)*S^7 +
	+ b*(-34 + 16*R*b^2 + 9*R^2*b^4)*S^6 + (-56 - 432*R*b^2 + 523*R^2*b^4)*S^5 +
	+ b*(568*R + 2026*R^2*b^2 - 4*R^3*b^4)*S^4 + (1120*R + 2472*R^2*b^2 - 596*R^3*b^4)*S^3 +
	- (2624*b*R^2 + 6608*b^3*R^3)*S^2 + (-3584*R^2 - 13184*b^2*R^3 + 256*b^4*R^4)*S +
	+ -7168*b*R^3 + 256*b^3*R^4
# (b*S^3 + 4*S^2 - 64*R) * (b*S + 1) * P[6]
### b*S + 1: Solution to "distinct" system (but x3 == x4);
### P[6]: Solution to degenerate System
# P[6] = (b*S^2 - 2*S - 4*b*R) * P[4]
(- 28*R + R^2*b^2) +
(- 24*R*b)*S^1 +
(7 - 2*R*b^2)*S^2 +
(- b)*S^3 +
(b^2)*S^4

##################

### Debug:
R = 2;
b = 3;
#
x1 =  1.6128981492115229;
x2 = -0.5852711593612232;
x3 = x4 = x2;


### Test
x1^2 + b*x2*x3*x4 # - R
x2^2 + b*x1*x3*x4 # - R
x3^2 + b*x1*x2*x4 # - R
x4^2 + b*x1*x2*x3 # - R


#######################
#######################

###############
### Order 3 ###
###############

### V3:
### x1^3 + b*x2*x3*x4 = R
x1^3 + b*x2*x3*x4 # - R
x2^3 + b*x1*x3*x4 # - R
x3^3 + b*x1*x2*x4 # - R
x4^3 + b*x1*x2*x3 # - R

### Solution:

# TODO:
# - general solver fails;
# - check all eqs: especially the Diff-Eqs;

# Note:
# - if (x1, x2, x3, x4) is a solution,
#   then m & m^2 * (x1, x2, x3, x4) are also solutions;

### Eq 1: Sum =>
S^3 - 3*E2*S + 3*E3 + b*E3 - 4*R # = 0

### Eq 2: Sum(x1*...) =>
(x1^4 + x2^4 + x3^4 + x4^4) + 4*b*E4 - R*S # = 0
S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 - 4*E4 + 4*b*E4 - R*S # = 0
# Reduction =>
E2*S^2 - 2*E2^2 + (b-1)*E3*S - 4*(b-1)*E4 - 3*R*S # = 0

### Eq 3: Diff Eq[i] - Eq[i+1] =>
(x1 - x2)*(x1^2 + x2^2 + x1*x2 - b*x3*x4) # = 0
# ...
# Case: x[i] != x[j]: Sum all (6) =>
3*(S^2 - 2*E2) + E2 - b*E2 # = 0
3*S^2 - 5*E2 - b*E2 # = 0

### Eq 4: Diff(x[i]*Eq[i] - x[i+1]*Eq[i+1]) =>
(x1 - x2)*(x1^3 + x2^3 + x1^2*x2 + x1*x2^2 - R) # = 0
(x2 - x3)*(x2^3 + x3^3 + x2^2*x3 + x2*x3^2 - R) # = 0
# ...
# sum(Diff) =>
3*(x1^3 + x2^3 + x3^3 + x4^3) + (E2*S - 3*E3) - 6*R # = 0
3*(S^3 - 3*E2*S + 3*E3) + E2*S - 3*E3 - 6*R # = 0
3*S^3 - 9*E2*S + 6*E3 + E2*S - 6*R # = 0
# Reduction =>
3*(b + 1)*E3 - E2*S - 6*R # = 0

### Alternatives:

### Diff(Eq 1 - Eq 3) =>
x1^3 - x3^3 - b*x2*x4*(x1 - x3) # = 0
x2^3 - x4^3 - b*x1*x3*(x2 - x4) # = 0
# Case: x1 != x3; x2 ! = x4;
s1^2 - p1 - b*p2 # = 0
s2^2 - p2 - b*p1 # = 0
# Sum =>
S^2 - 2*ps - (b+1)*sp # = 0
S^2 - 2*E2 - (b-1)*sp # = 0
S^2 - 2*E2 - bd*sp # = 0
# =>
S^6 - (bd + 6)*E2*S^4 + bd^2*E3*S^3 + 4*(bd + 3)*E2^2*S^2 - bd^2*(bd + 4)*E4*S^2 +
	- 2*bd^2*E2*E3*S - 4*(bd + 2)*E2^3 + 4*(2*bd^2 + bd^3)*E2*E4 - bd^3*E3^2 # = 0
# Reduction =>
bd^2*(bd + 4)*E4*S^2 +
	- bd*(bd + 1)*E2*E3*S + (bd - 9)*R*E2*S +
	+ 2*(bd + 1)*E2^3 - 4*bd*(bd^2 + 3*bd + 3)*E2*E4 +
	+ (2*bd^3 + 3*bd^2 - 8*bd - 16)*E3^2 - 4*(bd^2 - 2*bd - 8)*R*E3 - 16*R^2 # = 0

### Sum(x[i]^2*Eq[i]) =>
# TODO: is it necessary?


### Relations:
# (b+5)*E2 = 3*S^2
# 3*(b + 1)*E3 = E2*S + 6*R;
# Alternative: independent of b
# 6*E3 = -(3*S^3 - 9*E2*S + E2*S - 6*R)
# (b+3)*E3 = -(S^3 - 3*E2*S - 4*R) # check for cyclic redundancy?
# 4*(b-1)*E4 = E2*S^2 - 2*E2^2 + (b-1)*E3*S - 3*R*S;


###
pP1 = toPoly.pm("S^3 - 3*E2*S + 3*E3 + b*E3 - 4*R")
pP2 = toPoly.pm("E2*S^2 - 2*E2^2 + (b-1)*E3*S - 4*(b-1)*E4 - 3*R*S")
pP3 = toPoly.pm("3*S^2 - 5*E2 - b*E2")
pP4 = toPoly.pm("3*(b + 1)*E3 - E2*S - 6*R")


pR = solve.lpm(pP3, pP4, pP1, xn=c("E2", "E3"))
pR = pR[[2]]
str(pR)

print.pm(pR$Rez, lead="S")
(5 + 11*b - 3*b^2 - b^3)*S^3 - 50*R + 30*b*R + 18*b^2*R + 2*b^3*R # = 0

###
pP1 = toPoly.pm("x2^4 - R*x2 + b*E4")
pP2 = toPoly.pm("x2^3 - s3*x2^2 + e2*x2 - p3")
pR = solve.pm(pP1, pP2, "x2", stop.at=1)

toCoeff(pR[[2]], "x2")

# tends to 0/0 for the Special cases;
x2_div = p3^2 + e2^3 - 2*p3*e2*s3 + e2*E4*b - s3^2*E4*b - 2*p3*R + 3*e2*s3*R - s3^3*R + R^2;
x2_0 = - p3*e2^2 + p3^2*s3 + p3*E4*b - 2*e2*s3*E4*b + s3^3*E4*b - p3*s3*R - E4*b*R;
- x2_0 / x2_div

### tends to 0/0 for Special cases;
pP1 = toPoly.pm("x3^4 - R*x3 + b*E4")
pP2 = toPoly.pm("x3^2 - s2*x3 + p2")
pR = solve.pm(pP1, pP2, "x3", stop.at=1)

toCoeff(pR[[2]], "x3")


solve.S4Ht.L1V3aP3 = function(R, b, debug=TRUE) {
	S3 = (50 - 30*b - 18*b^2 + 2*b^3)*R;
	S3 = - S3 / ((b + 5)*(b^2 - 2*b - 1));
	S = rootn(S3, 3);
	m = unity(3, all=TRUE);
	S = S * m;
	if(debug) print(S);
	#
	E2 = 3*S^2 / (b + 5);
	E3 = (E2*S + 6*R) / (3*(b + 1));
	E4 = (E2*S^2 - 2*E2^2 + (b-1)*E3*S - 3*R*S) / (4*(b-1));
	#
	len = length(S);
	x1 = sapply(seq(len), function(id) {
		roots(c(1, -S[id], E2[id], -E3[id], E4[id]));
	})
	x1 = as.vector(x1);
	S  = rep(S, each=4);
	E3 = rep(E3, each=4);
	E4 = rep(E4, each=4);
	#
	s3 = S - x1;
	p3 = (R - x1^3) / b;
	e2 = (E3 - p3) / x1;
	#
	div = p3^2 + e2^3 - 2*p3*e2*s3 + e2*E4*b - s3^2*E4*b - 2*p3*R + 3*e2*s3*R - s3^3*R + R^2;
	x2 = - p3*e2^2 + p3^2*s3 + p3*E4*b - 2*e2*s3*E4*b + s3^3*E4*b - p3*s3*R - E4*b*R;
	x2 = - x2 / div;
	#
	s2 = s3 - x2;
	p2 = p3 / x2;
	x3 = (E4*b + p2^2 - p2*s2^2) / (R + 2*p2*s2 - s2^3);
	x4 = s2 - x3;
	#
	sol = cbind(x1, x2, x3, x4);
	return(sol);
}
# Special Case: x1 == x3;
solve.S4Ht.L1V3aP3.Case13 = function(R, b, debug=TRUE) {
	coeff = c((b^6 + b^4 - b^2 - 1), - (b^4 - 2*b^2 - 3)*R, - (b^3 + b^2 + 3)*R^2, R^3);
	x1_3 = roots(coeff);
	x1 = rootn(x1_3, 3);
	m = unity(3, all=TRUE);
	x1 = sapply(x1, function(x) x*m);
	if(debug) print(x1);
	x1 = as.vector(x1);
	#
	p2 = (R - x1^3) / (b*x1);
	s2 = b*R*x1 / ((b^2 + 1)*x1^3 - R);
	#
	len = length(s2);
	x24 = sapply(seq(len), function(id) {
		roots(c(1, -s2[id], p2[id]));
	})
	x24 = t(x24);
	x2 = x24[,1]; x4 = x24[,2];
	#
	sol = cbind(x1, x2, x3=x1, x4);
	return(sol);
}

###
R = 3
b = -7
sol = solve.S4Ht.L1V3aP3.Case13(R, b)
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];

###
x1^3 + b*x2*x3*x4 # - R
x2^3 + b*x1*x3*x4 # - R
x3^3 + b*x1*x2*x4 # - R
x4^3 + b*x1*x2*x3 # - R


### Debug:
R = 4; b = -3; bd = b - 1;
# Note: x2 == x4, which breaks Eq 3 above (temporarily);
x1 =  0.2924017738 + 0.5064547285i;
x2 =  1.5329574364 - 0.2097804173i;
x3 =  0.2924017738 + 0.5064547285i;
x4 = -0.9481538887 + 1.2226898741i;
x = c(x1,x2,x3,x4)
s1 = x1 + x3; s2 = x2 + x4;
p1 = x1 * x3; p2 = x2 * x4;
sp = p1 + p2; ps = s1 * s2;
S = s1 + s2; E4 = p1 * p2;
E2 = sp + ps;
E3 = p1*s2 + p2*s1;

###
s3 = S - x1;
p3 = (R - x1^3) / b;
e2 = (E3 - p3) / x1;

### Case: x1 == x3;
pP1 = toPoly.pm("s2^3 - 2*p2*s2 - R");
pP2 = toPoly.pm("s2^2 - p2 - b*x1^2");
pP3 = toPoly.pm("x1^3 + b*x1*p2 - R");
pR = solve.lpm(pP1, pP2, pP3, xn=c("p2", "s2"))

str(pR)

(b^6 + b^4 - b^2 - 1)*x1^9 - (b^4 - 2*b^2 - 3)*R*x1^6 +
	- (b^3 + b^2 + 3)*R^2*x1^3 + R^3 # = 0

