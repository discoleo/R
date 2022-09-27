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
### draft v.0.1e


####################
####################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


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

### Eq 1: Sum =>
S^2 - 2*E2 + b*E3 - 4*R # = 0

### Eq 2: Sum(x1*...) =>
(x1^3 + x2^3 + x3^3 + x4^3) + 4*b*E4 - R*S # = 0
S^3 - 3*E2*S + 3*E3 + 4*b*E4 - R*S # = 0
# Reduction =>
E2*S + b*E3*S - 4*b*E4 - 3*E3 - 3*R*S # = 0

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
# Reduction =>
b*E2*E3 - 2*R*E2 + E3*S - 3*b*E4*S - 4*E4 # = 0

### Sum(x1^3*...) =>
(x1^5 + x2^5 + x3^5 + x4^5) + b*E4*(S^2 - 2*E2) - R*(S^3 - 3*E2*S + 3*E3) # = 0
S^5 - 5*S^3*E2 + 5*S*E2^2 + 5*S^2*E3 - 5*E2*E3 - 5*S*E4 +
	+ b*E4*(S^2 - 2*E2) - R*(S^3 - 3*E2*S + 3*E3) # = 0
5*E3*S^2 + b*E4*S^2 - E2^2*S - 3*R*E2*S + 12*R^2*S - 5*E4*S + b*E3*E2*S - 7*b*R*E3*S +
	+ b^2*E3^2*S - 5*E3*E2 - 3*R*E3 - 2*b*E2*E4 # = 0


### Stable Eqs:
S^2 - 2*E2 + b*E3 - 4*R # = 0
E2*S + b*E3*S - 4*b*E4 - 3*E3 - 3*R*S # = 0
b*E2*E3 - 2*R*E2 + E3*S - 3*b*E4*S - 4*E4 # = 0
5*E3*S^2 + b*E4*S^2 - E2^2*S - 3*R*E2*S + 12*R^2*S - 5*E4*S + b*E3*E2*S - 7*b*R*E3*S +
	+ b^2*E3^2*S - 5*E3*E2 - 3*R*E3 - 2*b*E2*E4 # = 0

### Eq S:
(b*S^3 + 4*S^2 - 64*R) * (b*S^2 - 2*S - 4*b*R) * (b*S + 1) *
(b^2*S^4 - b*S^3 + 7*S^2 - 2*b^2*R*S^2 - 24*b*R*S + b^2*R^2 - 28*R) # = 0


###
pP1 = toPoly.pm("S^2 - 2*E2 + b*E3 - 4*R");
pP2 = toPoly.pm("E2*S + b*E3*S - 4*b*E4 - 3*E3 - 3*R*S");
pP3 = toPoly.pm("b*E2*E3 - 2*R*E2 + E3*S - 3*b*E4*S - 4*E4");
pP4 = toPoly.pm("5*E3*S^2 + b*E4*S^2 - E2^2*S - 3*R*E2*S + 12*R^2*S - 5*E4*S + b*E3*E2*S - 7*b*R*E3*S +
	+ b^2*E3^2*S - 5*E3*E2 - 3*R*E3 - 2*b*E2*E4");

pR = solve.lpm(pP1, pP2, pP3, pP4, xn=c("E3", "E4", "E2"))
pR = pR[[3]];
pR$Rez$coeff = - pR$Rez$coeff;
pR$Rez = sort.pm(pR$Rez, "S")
pR$Rez = orderVars.pm(pR$Rez, c("b", "R"), last=FALSE)

pR$Rez = div.pm(pR$Rez, toPoly.pm("b*S^3 + 4*S^2 - 64*R"), "S")$Rez;
pR$Rez = div.pm(pR$Rez, toPoly.pm("b*S^2 - 2*S - 4*b*R"), "S")$Rez;

# print.pm(pR$Rez, lead="S")

str(pR)
toCoeff(pR$Rez, "S")

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

### [old]
### full Polynomial:
# - overinflated P[16];
# - P[10] is the correct polynomial;
# - P[16]: [checked] coefficients are correct,
#   but has been deleted as it is redundant (and overinflated);
# - P[16] can be factored into: P[10] * P[6]
# - P[6]: 9*b^4*S^6 + ...; # is NOT part of solution!
# - P[10]: can be factored itself;
#   P[10] = P[3]*P[1]*P[2]*P[4];

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


### Case: x2 = x3 = x4
x1^2 + b*x2^3 - R # = 0
x2^2 + b*x1*x2^2 - R # = 0

### Classic Poly:
b^2*x1^2*x2^4 + b^3*x2^7 - b^2*R*x2^4 # = 0
(R - x2^2)^2 + b^3*x2^7 - b^2*R*x2^4 # = 0
# =>
b^3*x2^7 - (b^2*R - 1)*x2^4 - 2*R*x2^2 + R^2 # = 0
(x2^2 + b*x2^3 - R) * (b^2*x2^4 - b*x2^3 + x2^2 - R)

### x1:
(x1^2 + b*x1^3 - R) * (b^2*x1^4 + 2*b*x1^3 + x1^2 - 2*R*b^2*x1^2 - 3*R*b*x1 + R^2*b^2 - R) # = 0


############
### Workout:

### [old]
# - not used anymore;

debug.old = function(R, b, x, E, tol) {
	len = length(E$S)
	S1 = E$S; # 1 copy;
	S = matrix(E$S, ncol=len, nrow=4, byrow=T)
	E2 = matrix(E$E2, ncol=len, nrow=4, byrow=T)
	E3 = matrix(E$E3, ncol=len, nrow=4, byrow=T)
	E4 = matrix(E$E4, ncol=len, nrow=4, byrow=T)
	isZero = round0(E4/x - (R - x^2)/b[1], tol=tol) == 0
	isZ = apply(isZero, 2, all)
	return(list(sol=cbind(x=as.vector(x)), S=S1, isZ=isZ, isZero=isZero))
}
solve.old = function(x, E) {
	len = length(E$S);
	S = matrix(E$S, ncol=len, nrow=4, byrow=T)
	E2 = matrix(E$E2, ncol=len, nrow=4, byrow=T)
	E3 = matrix(E$E3, ncol=len, nrow=4, byrow=T)
	E4 = matrix(E$E4, ncol=len, nrow=4, byrow=T)
	SS3  = S - x; E2S3 = E2 - x*SS3; E3S3 = E4 / x;
	# x2
	x2 = sapply(seq_along(x), function(id) roots(c(1, -SS3[id], E2S3[id], -E3S3[id])))
	#
	x2 = as.vector(x2);
	x = rep(as.vector(x), each=3);
	SS2 = rep(as.vector(SS3), each=3) - x2;
	E2S2 = rep(as.vector(E3S3), each=3) / x2;
	# TODO: root[2]
	x3 = sapply(seq_along(x2), function(id) roots(c(1, -SS2[id], E2S2[id]))[1])
	x4 = rep(SS2, each=1) - x3; sol=cbind(x1=x, x2=x2, x3=x3, x4=x4);
}

b = -4:4
b = b[b != 0]
R = 2;

sapply(b, function(b) {
	sol = solve.S4(R=R, b, tol=5E-2); # tol=5E-2
	table(sol$isZ)[1]
	} )

R = 2
b = 3
sol = solve.S4(R=R, b=b, tol=5E-2)
poly.calc(sol$S[ ! sol$isZ]) * 9 *b^4

S = sol$S; # ...


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


### Eq 3: Sum(x1^2*...) =>
(x1^5 + x2^5 + x3^5 + x4^5) + b*E4*S - R*(S^2 - 2*E2) # = 0
S^5 - 5*E2*S^3 + 5*E3*S^2 - R*S^2 + 5*E2^2*S + b*E4*S - 5*E4*S - 5*E2*E3 + 2*R*E2 # = 0
# Reduction =>
2*E2*S^3 - 3*R*S^2 - 2*E3*S^2 + b*E3*S^2 - 5*E2^2*S + 5*E4*S - E4*b*S - 2*E2*R + 5*E2*E3 # = 0
(bd + 1)*E3*S^2 - 3*R*S^2 + E2^2*S - (7*bd + 4)*E4*S + 2*R*E2 - 5*E2*E3 # = 0


### Alternatives: Diff

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

### Alternatives: more Diff

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


pP1 = toPoly.pm("S^3 - 3*E2*S + 3*E3 + b*E3 - 4*R")
pP2 = toPoly.pm("E2*S^2 - 2*E2^2 + (b-1)*E3*S - 4*(b-1)*E4 - 3*R*S")
pP3 = toPoly.pm("b*E3*S^2 - 3*R*S^2 + E2^2*S - (7*b - 3)*E4*S + 2*R*E2 - 5*E2*E3")
pP4 = toPoly.pm("3*b*E4*S^2 - (b + 1)*E2*E3*S + 3*E3^2 +
	+ (2*b + 2)*E2*E4 + 2*R*E2*S - 3*R*E3")

pE3 = toPoly.pm("S^3 - 3*E2*S - 4*R")
pE3div = toPoly.pm("- b - 3")

pP2 = replace.fr.pm(pP2, pE3, pE3div, "E3")
pP3 = replace.fr.pm(pP3, pE3, pE3div, "E3")
pP4 = replace.fr.pm(pP4, pE3, pE3div, "E3")

###
pE4 = toPoly.pm("(b - 1)*S^4 - 4*b*E2*S^2 - b*R*S + 13*R*S + 2*b*E2^2 + 6*E2^2")
pE4div = toPoly.pm("- (4*b^2 + 8*b - 12)")

pP3 = replace.fr.pm(pP3, pE4, pE4div, "E4")
pP4 = replace.fr.pm(pP4, pE4, pE4div, "E4")

c0 = 8*(b + 3)*(b - 1)*(b - 7)*R; # (168 - 136*b - 40*b^2 + 8*b^3)*R;
c1 = (90 - 18*b + 38*b^2 + 18*b^3);
c2 = - (9 - 159*b - 45*b^2 + 3*b^3)*R;
c3 = - (60 - 40*b + 28*b^2 + 16*b^3);
c4 = 9 - 15*b + 3*b^2 + 3*b^3;
#
d0 = 144*R^2 - 240*b*R^2 + 48*b^2*R^2 + 48*b^3*R^2;
d1 = 108 + 216*b + 144*b^2 + 40*b^3 + 4*b^4;
d2 = (846 - 192*b + 28*b^2 + 80*b^3 + 6*b^4)*R;
d3 = 216 - 198*b + 30*b^2 + 70*b^3 + 10*b^4;
d4 = (- 180 + 507*b + 243*b^2 + 9*b^3 - 3*b^4)*R;
d5 = - 198 + 156*b - 52*b^2 - 84*b^3 - 14*b^4;
d6 = 36 - 51*b - 3*b^2 + 15*b^3 + 3*b^4;

c1*E2^2*S + c3*E2*S^3 + c0*E2 + c4*S^5 + c2*S^2 # = 0
d1*E2^3 + d6*S^6 + d5*E2*S^4 + d4*S^3 + d3*E2^2*S^2 + d2*E2*S + d0 # = 0

### Reduction:
d6 = c1*d6;
d5 = c1*d5 - c4*d1;
d4 = c1*d4;
d3 = c1*d3 - c3*d1;
d2 = c1*d2 - c2*d1;
d0 = c1*d0;
d7 = - c0*d1;
d1 = c1*d1; # Last!
#
d7 = c1*d7;
d5 = c1*d5 - d3*c3;
d2 = c1*d2 - d3*c0;
d6 = c1*d6 - d3*c4;
d4 = c1*d4 - d3*c2;
d0 = c1*d0;

d7*E2^2 + d5*E2*S^5 + d2*E2*S^2 + d6*S^7 + d4*S^4 + d0*S # = 0

###
pP3 = toPoly.pm("c4*S^5 + c3*E2*S^3 + c2*S^2 + c1*E2^2*S + c0*E2")
pP4 = toPoly.pm("d7*E2^2 + d5*E2*S^5 + d2*E2*S^2 + d6*S^7 + d4*S^4 + d0*S")

pR = solve.pm(pP3, pP4, "E2")
str(pR)

# P[15] = P[3] * P[12];
# P[12] = ((b+1)*S^3 - 64*R) * ((b-1)*S^3 - 8*(b+1)*R) * P[6];
# P[3]: NO solutions;
# P[6]: probably Case: x2 = x3 = x4;
# TODO: Method to factor the compact P[15] (?);
# - done on the expanded P[15] (not that big after reductions);

coeff.S4Ht.L1V3bP3.old = function(R, b) {
	# [old] NOT needed anymore;
	cc.all = coeff0.S4Ht.L1V3bP3(R, b=b);
	with(cc.all, {
		coeff = c(c1^2*d6^2 - c1*c3*d6*d5 + c1*c4*d5^2,
			2*c1^2*d4*d6 - c1*c3*d6*d2 - c1*c3*d4*d5 - c0*c1*d6*d5 + 2*c1*c4*d2*d5 +
				+ c1*c2*d5^2 + c3^2*d6*d7 - 2*c1*c4*d6*d7 - c3*c4*d5*d7,
			c1^2*d4^2 + 2*c1^2*d0*d6 - c1*c3*d4*d2 - c0*c1*d6*d2 + c1*c4*d2^2 - c1*c3*d0*d5 +
				- c0*c1*d4*d5 + 2*c1*c2*d2*d5 + c3^2*d4*d7 - 2*c1*c4*d4*d7 - 2*c1*c2*d6*d7 +
				+ 2*c0*c3*d6*d7 - c3*c4*d2*d7 - c2*c3*d5*d7 - c0*c4*d5*d7 + c4^2*d7^2,
			2*c1^2*d0*d4 - c1*c3*d0*d2 - c0*c1*d4*d2 + c1*c2*d2^2 - c0*c1*d0*d5 + c3^2*d0*d7 +
				- 2*c1*c4*d0*d7 - 2*c1*c2*d4*d7 + 2*c0*c3*d4*d7 + c0^2*d6*d7 - c2*c3*d2*d7 +
				- c0*c4*d2*d7 - c0*c2*d5*d7 + 2*c2*c4*d7^2,
			c1^2*d0^2 - c0*c1*d0*d2 - 2*c1*c2*d0*d7 + 2*c0*c3*d0*d7 + c0^2*d4*d7 +
				- c0*c2*d2*d7 + c2^2*d7^2,
			c0^2*d0*d7);
		return(coeff);
	})
}
coeff0.S4Ht.L1V3bP3 = function(R, b) {
	c0 = 8*(b + 3)*(b - 1)*(b - 7)*R; # (168 - 136*b - 40*b^2 + 8*b^3)*R;
	c1 = (90 - 18*b + 38*b^2 + 18*b^3);
	c2 = - (9 - 159*b - 45*b^2 + 3*b^3)*R;
	c3 = - (60 - 40*b + 28*b^2 + 16*b^3);
	c4 = 9 - 15*b + 3*b^2 + 3*b^3;
	#
	d0 = 144*R^2 - 240*b*R^2 + 48*b^2*R^2 + 48*b^3*R^2;
	d1 = 108 + 216*b + 144*b^2 + 40*b^3 + 4*b^4;
	d2 = (846 - 192*b + 28*b^2 + 80*b^3 + 6*b^4)*R;
	d3 = 216 - 198*b + 30*b^2 + 70*b^3 + 10*b^4;
	d4 = (- 180 + 507*b + 243*b^2 + 9*b^3 - 3*b^4)*R;
	d5 = - 198 + 156*b - 52*b^2 - 84*b^3 - 14*b^4;
	d6 = 36 - 51*b - 3*b^2 + 15*b^3 + 3*b^4;
	### Reduction:
	d6 = c1*d6; d5 = c1*d5 - c4*d1;
	d4 = c1*d4; d3 = c1*d3 - c3*d1;
	d2 = c1*d2 - c2*d1;
	d0 = c1*d0; d7 = - c0*d1;
	d1 = c1*d1; # Last!
	#
	d7 = c1*d7; d0 = c1*d0;
	d5 = c1*d5 - d3*c3; d2 = c1*d2 - d3*c0;
	d6 = c1*d6 - d3*c4; d4 = c1*d4 - d3*c2;
	c.l = list(c0=c0, c1=c1, c2=c2, c3=c3, c4=c4);
	d.l = list(d0=d0, d2=d2, d4=d4, d5=d5, d6=d6, d7=d7);
	return(c(c.l, d.l));
}

##################

### [NOT working!]
pP1 = toPoly.pm("S^3 - 3*E2*S + 3*E3 + b*E3 - 4*R")
pP2 = toPoly.pm("E2*S^2 - 2*E2^2 + (b-1)*E3*S - 4*(b-1)*E4 - 3*R*S")
pP3 = toPoly.pm("3*S^2 - 5*E2 - b*E2")
pP4 = toPoly.pm("3*(b + 1)*E3 - E2*S - 6*R")

### [NOT working] Redundancy (?)
pP1 = toPoly.pm("S^3 - 3*E2*S + 3*E3 + b*E3 - 4*R")
pP2 = toPoly.pm("E2*S^2 - 2*E2^2 + (b-1)*E3*S - 4*(b-1)*E4 - 3*R*S")
pP3 = toPoly.pm("b*E3*S^2 - 3*R*S^2 + E2^2*S - (7*b - 3)*E4*S + 2*R*E2 - 5*E2*E3")
pP4 = toPoly.pm("bd^2*(bd + 4)*E4*S^2 +
	- bd*(bd + 1)*E2*E3*S + (bd - 9)*R*E2*S +
	+ 2*(bd + 1)*E2^3 - 4*bd*(bd^2 + 3*bd + 3)*E2*E4 +
	+ (2*bd^3 + 3*bd^2 - 8*bd - 16)*E3^2 - 4*(bd^2 - 2*bd - 8)*R*E3 - 16*R^2")

pE3 = toPoly.pm("S^3 - 3*E2*S - 4*R")
pE3div = toPoly.pm("- b - 3")

pP2 = replace.fr.pm(pP2, pE3, pE3div, "E3")
pP3 = replace.fr.pm(pP3, pE3, pE3div, "E3")
pP4 = replace.fr.pm(pP4, pE3, pE3div, "E3")


bd = b - 1;
dd = 4*b^2 + 8*b - 12;
c0 = 32*bd^3*R^2 - 16*b*bd^2*R^2 - 16*b^2*R^2 + 32*b*bd*R^2 + 32*b*R^2 - 32*bd*R^2 - 16*R^2;
c1 = - 4*b^2*bd^3 - 24*b*bd^3 - 12*b^2*bd^2 - 36*bd^3 - 12*b^2*bd - 72*b*bd^2 - 108*bd^2 +
	- 72*b*bd - 108*bd;
c2 = 2*b^2*bd + 2*b^2 + 12*b*bd + 12*b + 18*bd + 18;
c3 = (48*bd^3 + b^2*bd - 16*b*bd^2 - 9*b^2 + 24*bd^2 + 26*b*bd + 42*b - 123*bd - 177)*R;
c4 = b^2*bd^3 + 6*b*bd^3 + 4*b^2*bd^2 + 9*bd^3 + 24*b*bd^2 + 36*bd^2;
c5 = 18*bd^3 - 3*b*bd^2 - 3*b*bd + 18*bd^2 - 81*bd - 144;
c6 = (4*b*bd^2 - 16*bd^3 - 8*b*bd - 12*bd^2 - 32*b + 40*bd + 32)*R;
c7 = b*bd^2 - 12*bd^3 + b*bd - 15*bd^2 + 51*bd + 96;
c8 = 2*bd^3 + 3*bd^2 - 8*bd - 16;
#
c7 = bd*c7 + 4*b*c8;
c6 = bd*c6 + (b - 13)*R*c8;
c5 = bd*c5 - (2*b + 6)*c8;
c4 = bd*c4 - (4*b^2 + 8*b - 12)*c8;
c3 = bd*c3;
c2 = bd*c2;
c1 = bd*c1;
c0 = bd*c0;

c7*E2*S^4 + c6*S^3 + c5*E2^2*S^2 +
	+ c4*E4*S^2 + c3*E2*S + c2*E2^3 + c1*E4*E2 + c0 # = 0

pP2 = toPoly.pm("(b - 1)*S^4 - 4*b*E2*S^2 - b*R*S + 13*R*S + (4*b^2 + 8*b - 12)*E4 +
	+ (2*b + 6)*E2^2")
pP3 = toPoly.pm("b*S^5 - (3*b + 5)*E2*S^3 - (b - 9)*R*S^2 - (b - 12)*E2^2*S +
	+ (7*b^2 + 18*b - 9)*E4*S - (2*b - 14)*R*E2")
pP4 = toPoly.pm("c7*E2*S^4 + c6*S^3 + c5*E2^2*S^2 +
	+ c4*E4*S^2 + c3*E2*S + c2*E2^3 + c1*E4*E2 + c0")


pE4 = toPoly.pm("bd*S^4 - 4*b*E2*S^2 - b*R*S + 13*R*S + (2*b + 6)*E2^2")
pE4div = toPoly.pm("- (4*b^2 + 8*b - 12)")
pE4div = toPoly.pm("- dd")

pP3 = replace.fr.pm(pP3, pE4, pE4div, "E4")
pP4 = replace.fr.pm(pP4, pE4, pE4div, "E4")

### Eq 3: Redundancy (?)
# - E4 does NOT appear in pP4 anymore;
d0 = (168 - 136*b - 40*b^2 + 8*b^3)*R;
d1 = 90 - 18*b + 38*b^2 + 18*b^3;
d2 = (- 9 + 159*b + 45*b^2 - 3*b^3)*R;
d3 = - 60 + 40*b - 28*b^2 - 16*b^3;
d4 = 9 - 15*b + 3*b^2 + 3*b^3;
#
e0 = - dd*c0;
e1 = 6*c1*E2^3 - dd*c2*E2^3 + 2*b*c1;
e2 = - dd*c3 + 13*c1*R - b*c1*R;
e3 = 6*c4 - dd*c5 - 4*b*c1 + 2*b*c4;
e4 = - dd*c6 + 13*c4*R - b*c4*R;
e5 = - dd*c7 - 4*b*c4 + c1*bd;
e6 = c4*bd;

d4*S^5 + d3*E2*S^3 + d2*R*S^2 + d1*E2^2*S + d0*E2 # = 0
e6*S^6 + e5*E2*S^4 + e4*S^3 + e3*E2^2*S^2 + e2*E2*S + e1*E2^3 + e0 # = 0

# d1*pP4*S - e1*pP3*E2
e6 = d1*e6;
e5 = d1*e5 - e1*d4;
e4 = d1*e4;
e3 = d1*e3 - e1*d3;
e2 = d1*e2;
e0 = d1*e0;

d1*E2^2*S + d3*E2*S^3 + d0*E2 + d4*S^5 + d2*R*S^2 # = 0
e3*E2^2*S^3 - e1*d0*E2^2 + e2*E2*S^2 + e6*S^7 + e5*E2*S^5 - e1*d2*R*E2*S^2 + e4*S^4 + e0*S # = 0


pP3 = toPoly.pm("d4*S^5 + d3*E2*S^3 + d2*R*S^2 + d1*E2^2*S + d0*E2")
pP4 = toPoly.pm("e6*S^6 + e5*E2*S^4 + e4*S^3 + e3*E2^2*S^2 + e2*E2*S + e1*E2^3 + e0")
	
pR = solve.lpm(pP3, pP4, xn=c("E2"))
pR = pR[[1]]
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


solve.S4Ht.L1V3bP3 = function(R, b, debug=TRUE) {
	S3 = (50 - 30*b - 18*b^2 - 2*b^3)*R;
	S3 = - S3 / ((b + 5)*(b^2 - 2*b - 1));
	S = rootn(S3, 3);
	m = unity(3, all=TRUE);
	S = S * m;
	if(debug) print(S);
	#
	E2 = 3*S^2 / (b + 5);
	# E3 = (E2*S + 6*R) / (3*(b + 1));
	E3 = - (S^3 - 3*E2*S - 4*R) / (b + 3); print(E3)
	E4 = (E2*S^2 - 2*E2^2 + (b-1)*E3*S - 3*R*S) / (4*(b-1)); print(E4)
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
solve.S4Ht.L1V3bP3.Case13 = function(R, b, debug=TRUE, all=FALSE) {
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
	if(all) {
		sol = rbind(sol, sol[ , c(2,4,1,3)]);
		sol = rbind(sol, sol[ , c(3,1,4,2)]);
	}
	return(sol);
}

###
R = 3
b = -7
sol = solve.S4Ht.L1V3bP3.Case13(R, b)
sol = solve.S4Ht.L1V3bP3(R, b)
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

###########

source("Polynomials.Helper.BigNumbers.R")

pP1 = toPoly.pm("S^3 - 3*E2*S + 3*E3 + b*E3 - 4*R")
pP2 = toPoly.pm("E2*S^2 - 2*E2^2 + (b-1)*E3*S - 4*(b-1)*E4 - 3*R*S")
pP3 = toPoly.pm("b*E3*S^2 - 3*R*S^2 + E2^2*S - (7*b - 3)*E4*S + 2*R*E2 - 5*E2*E3")
pP4 = toPoly.pm("3*b*E4*S^2 - (b + 1)*E2*E3*S + 3*E3^2 +
	+ (2*b + 2)*E2*E4 + 2*R*E2*S - 3*R*E3")

pP1$coeff = as.bigz(pP1$coeff);
pP2$coeff = as.bigz(pP2$coeff);
pP3$coeff = as.bigz(pP3$coeff);
pP4$coeff = as.bigz(pP4$coeff);

pR = solve.lpm(pP1, pP2, pP3, pP4, xn=c("E3", "E4", "E2"), asBigNum=TRUE)

pR = pR[[3]]$Rez

tmp = div.pm(pR, toPoly.pm("((b+1)*S^3 - 64*R) * ((b-1)*S^3 - 8*(b+1)*R)"), c("S", "b"))
tmp = tmp$Rez;

tmp2 = tmp[tmp$S == 9, c("b", "coeff")]
pDiv = div.pm(tmp2, toPoly.pm("(b^2+1)^2*(b^2-1)*(b-1)^2"), "b")$Rez
tmp = div.pm(tmp, pDiv, "b")
tmp = tmp$Rez;
toCoeff(tmp, "S")

