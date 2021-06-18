########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
###  == Derivation ==
###
### draft v.0.3a-Eq3


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R
# - e.g. round0(), round0.p(),
#   solve.EnAll(), solveEn();

### other

### Simple:
xi.f = function(x, R, b, n=2) {
	(R - x^n) / b[1];
}
xip.f = function(x, R, b, n=2, p=1) {
	(R - x^n) / b[1] / x^p;
}
test.S4.Simple = function(sol, R, b, n=2) {
	x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];
	err1 = x1^n + b*x2 # - R
	err2 = x2^n + b*x3 # - R
	err3 = x3^n + b*x4 # - R
	err4 = x4^n + b*x1 # - R
	err = rbind(err1, err2, err3, err4);
	if( ! missing(R)) err = err - R;
	err = round0(err);
	return(err);
}
debug.E = function(x) {
	x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4];
	S = sum(x)
	E4 = prod(x)
	E3 = prod(x)*sum(1/x)
	m = perm2(4)
	E2 = sum(sapply(seq(nrow(m)), function(id) prod(x[which(m[id,] != 0)])))
	data.frame(S=S, E2=E2, E3=E3, E4=E4);
}

########################

###############
### Order 2 ###
###############

### V3:
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

### Diff Eq[i] - Eq[i+1] =>
(x1 - x2)*(x1^2 + x2^2 + x1*x2 - b*x3*x4) # = 0
# ...
# Case: x[i] != x[j]: Sum =>
3*S^2 - 5*E2 - b*E2 # = 0

### Sum =>
S^3 - 3*E2*S + 3*E3 + b*E3 - 4*R # = 0

### Sum(x1*...) =>
(x1^4 + x2^4 + x3^4 + x4^4) + 4*b*E4 - R*S # = 0
S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 - 4*E4 + 4*b*E4 - R*S # = 0

### Diff(x[i]*Eq[i] - x[i+1]*Eq[i+1]) =>
(x1 - x2)*(x1^3 + x2^3 + x1^2*x2 + x1*x2^2 - R) # = 0
(x2 - x3)*(x2^3 + x3^3 + x2^2*x3 + x2*x3^2 - R) # = 0
# ...
# sum(Diff) =>
3*(x1^3 + x2^3 + x3^3 + x4^3) + (E2*S - 3*E3) - 6*R # = 0
3*(S^3 - 3*E2*S + 3*E3) + E2*S - 3*E3 - 6*R # = 0
3*S^3 - 9*E2*S + 6*E3 + E2*S - 6*R # = 0

### Sum(x[i]^2*Eq[i]) =>
# TODO: is it necessary?


### Relations:
# (b+5)*E2 = 3*S^2
# 6*E3 = -(3*S^3 - 9*E2*S + E2*S - 6*R)
# (b+3)*E3 = -(S^3 - 3*E2*S - 4*R) # check for cyclic redundancy?
# 4*(b-1)*E4 = S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 - R*S


###########################
###########################
###########################

##############
### Simple ###
##############

###############
### Order 2 ###
###############

### x[i]^2 + b*x[i+1] = R

### Solution:

### Case: all x[i] different;
# - root structure: P[4] o P[3];
### Special Sub-Case:
# x1 = Conj(x3); x2 = Conj(x4);
# - 4 roots;
# - properties: Re(x1) * Re(x2) = - b^2 / 4;
# Remaining roots: 8;

### Sum =>
S^2 - 2*E2 + b*S - 4*R # = 0

### Diff =>
(x1+x2)*(x1+x3)*(x1+x4)*(x2+x3)*(x2+x4)*(x3+x4) + b^6 # = 0
(-E3^2 + E3*E2*S - E4*S^2) + b^6 # = 0

### Eq 3:
# x1^2 = R - b*x2 => Prod =>
(x1*x2*x3*x4)^2 - R^4 + b*S*R^3 - b^2*E2*R^2 + b^3*E3*R - b^4*E4 # = 0
E4^2 - b^4*E4 + b*R^3*S - b^2*R^2*E2 + b^3*R*E3 - R^4 # = 0

### Eq 4:
# b*x2 = R - x1^2 => Prod =>
b^4*E4 - R^4 + R^3*(x1^2+x2^2+x3^2+x4^2) - R^2*E2_2 + R*E3_2 - E4^2 # = 0
b^4*E4 - R^4 + R^3*(S^2 - 2*E2) - R^2*(E2^2 - 2*E3*S + 2*E4) + R*(E3^2 - 2*E4*E2) - E4^2 # = 0
E4^2 - b^4*E4 + 2*R^2*E4 + 2*R*E2*E4 - R*E3^2 - 2*R^2*E3*S +
	+ 2*R^3*E2 + R^2*E2^2 - R^3*S^2 + R^4 # = 0
# Diff(Eq 4 - Eq 3) =>
2*R^2*E4 + 2*R*E2*E4 - R*E3^2 - 2*R^2*E3*S - b^3*R*E3 +
	+ R^2*E2^2 + 2*R^3*E2 + b^2*R^2*E2 - R^3*S^2 - b*R^3*S + 2*R^4 # = 0

### Eq:
S^3 + (3*b^2 - 4*R)*S - 4*b^3 # = 0

### Auxiliary Eqs:
# E2 = (S^2 + b*S - 4*R) / 2;
# E3 = ...; [see file Poly.System.Hetero.Symmetric.S4.R]
# E4 = (-E3^2 + E3*E2*S + b^6) / S^2;

### [old]
# - simplified version available in the primary file;
E3.helper.old = function(S, E2, R, b) {
pE3 = - 4*E2*R*S^6*b^12 + 4*E2*R*S^8*b^10 + 8*E2*R^2*S^6*b^10 - 12*E2*R^2*S^7*b^9 + E2*R^2*S^8*b^8 +
	- E2*R^2*S^10*b^6 + 4*E2*R^3*S^6*b^8 - 4*E2*R^3*S^7*b^7 - 22*E2*R^3*S^8*b^6 + 4*E2*R^3*S^10*b^4 +
	+ 2*E2*R^3*S^11*b^3 + E2*R^3*S^12*b^2 + 12*E2*R^4*S^6*b^6 + 8*E2*R^4*S^7*b^5 - 8*E2*R^4*S^8*b^4 +
	- 10*E2*R^4*S^9*b^3 + 4*E2*R^4*S^10*b^2 - 2*E2*R^4*S^11*b + 4*E2*R^4*S^12 - 24*E2*R^5*S^6*b^4 +
	+ 4*E2*R^5*S^7*b^3 + 10*E2*R^5*S^8*b^2 - 24*E2*R^5*S^9*b - 8*E2*R^5*S^10 - 8*E2*R^6*S^7*b +
	+ 32*E2*R^6*S^8 + 4*E2^2*R*S^6*b^10 - 4*E2^2*R*S^7*b^9 + 4*E2^2*R^2*S^6*b^8 +
	- 13*E2^2*R^2*S^8*b^6 - 2*E2^2*R^2*S^9*b^5 - E2^2*R^2*S^10*b^4 + 4*E2^2*R^3*S^6*b^6 +
	+ 4*E2^2*R^3*S^7*b^5 - 9*E2^2*R^3*S^8*b^4 - 8*E2^2*R^3*S^9*b^3 + 2*E2^2*R^3*S^10*b^2 +
	+ E2^2*R^3*S^11*b + E2^2*R^3*S^12 - 30*E2^2*R^4*S^6*b^4 + 4*E2^2*R^4*S^7*b^3 +
	+ 24*E2^2*R^4*S^8*b^2 - 10*E2^2*R^4*S^9*b - 12*E2^2*R^4*S^10 + 8*E2^2*R^5*S^6*b^2 +
	- 12*E2^2*R^5*S^7*b + 16*E2^2*R^5*S^8 - 8*E2^2*R^6*S^6 - 2*E2^3*R*S^8*b^6 - 4*E2^3*R^2*S^8*b^4 +
	- 2*E2^3*R^2*S^9*b^3 - E2^3*R^2*S^10*b^2 - 18*E2^3*R^3*S^6*b^4 + 10*E2^3*R^3*S^8*b^2 +
	- 6*E2^3*R^3*S^10 + 12*E2^3*R^4*S^6*b^2 - 4*E2^3*R^4*S^7*b - 16*E2^3*R^5*S^6 +
	- 4*E2^4*R^2*S^6*b^4 - E2^4*R^2*S^10 + 4*E2^4*R^3*S^6*b^2 - E2^4*R^3*S^8 - 10*E2^4*R^4*S^6 +
	- 2*E2^5*R^3*S^6 - R*S^8*b^12 + R*S^10*b^10 - 4*R^2*S^6*b^12 + 4*R^2*S^8*b^10 + 4*R^3*S^6*b^10 +
	- 8*R^3*S^7*b^9 - R^3*S^9*b^7 - R^3*S^10*b^6 + R^3*S^11*b^5 + R^3*S^12*b^4 - 4*R^4*S^7*b^7 +
	- 10*R^4*S^8*b^6 + 4*R^4*S^9*b^5 + 6*R^4*S^10*b^4 + 4*R^4*S^11*b^3 - R^4*S^13*b + 8*R^5*S^6*b^6 +
	+ 4*R^5*S^7*b^5 - 4*R^5*S^8*b^4 - 8*R^5*S^9*b^3 - R^5*S^10*b^2 - 4*R^5*S^11*b + 4*R^5*S^12 +
	- 8*R^6*S^6*b^4 - 2*R^6*S^8*b^2 - 12*R^6*S^9*b + 16*R^7*S^8;
pE3div = - 4*E2*R*S^6*b^9 + 4*E2*R*S^8*b^7 + E2*R*S^11*b^4 + 8*E2*R^2*S^6*b^7 - 14*E2*R^2*S^7*b^6 +
	+ 2*E2*R^2*S^8*b^5 + 12*E2*R^2*S^9*b^4 + 2*E2*R^2*S^10*b^3 + 4*E2*R^3*S^6*b^5 +
	+ 16*E2*R^3*S^7*b^4 - 32*E2*R^3*S^8*b^3 + 4*E2*R^3*S^9*b^2 - 2*E2*R^3*S^10*b + 6*E2*R^3*S^11 +
	- 8*E2*R^4*S^6*b^3 + 8*E2*R^4*S^7*b^2 - 12*E2*R^4*S^8*b - 8*E2*R^4*S^9 + 40*E2*R^5*S^7 +
	+ 4*E2^2*R*S^6*b^7 - 6*E2^2*R*S^7*b^6 + 4*E2^2*R*S^9*b^4 + E2^2*R*S^10*b^3 + 4*E2^2*R^2*S^6*b^5 +
	+ 16*E2^2*R^2*S^7*b^4 - 22*E2^2*R^2*S^8*b^3 + 2*E2^2*R^2*S^9*b^2 + 2*E2^2*R^2*S^11 +
	- 12*E2^2*R^3*S^6*b^3 + 12*E2^2*R^3*S^7*b^2 - 4*E2^2*R^3*S^8*b - 12*E2^2*R^3*S^9 +
	+ 40*E2^2*R^4*S^7 + 4*E2^3*R*S^7*b^4 - 4*E2^3*R*S^8*b^3 - 4*E2^3*R^2*S^6*b^3 +
	+ 4*E2^3*R^2*S^7*b^2 - 8*E2^3*R^2*S^9 + 20*E2^3*R^3*S^7 - 2*E2^4*R*S^9 + 4*E2^4*R^2*S^7 +
	- R*S^8*b^9 + R*S^10*b^7 - 4*R^2*S^6*b^9 + 4*R^2*S^8*b^7 + 2*R^2*S^9*b^6 + 2*R^2*S^11*b^4 +
	- R^2*S^12*b^3 + 4*R^3*S^6*b^7 - 8*R^3*S^7*b^6 + 6*R^3*S^9*b^4 + 4*R^3*S^10*b^3 +
	+ 4*R^4*S^7*b^4 - 12*R^4*S^8*b^3 - 4*R^4*S^10*b + 4*R^4*S^11 - 8*R^5*S^8*b + 16*R^6*S^7;
- pE3 / pE3div;
}


### Solver:
R = -1
b = 3
sol = solve.Simple.S4P2(R, b);
E = do.call(rbind, lapply(seq(nrow(sol)), function(id) debug.E(sol[id,])));
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];
S = E$S; E2 = E$E2; E3 = E$E3; E4 = E$E4;

test.S4P2.Simple(sol, b=b)


### Helper Derivation:

### Elementary Polynomials:
n = 4
p = prod.perm.poly(n)
p = sort.pm(p);
p

(x1^3*x2^2*x3 + ...) +
	+ 2*(x1^2*x2^2*x3^2 + ...) +
	+ 2*x1*x2*x3*x4*(x1^2+x2^2+x3^2+x4^2) +
	+ 4*x1*x2*x3*x4*E2;
(-3*E3^2 + E3*E2*S - 3*E4*S^2 + 6*E4*E2 - 2*E4*E2) +
	+ 2*(E3^2 - 2*E4*E2) +
	+ 2*E4*(S^2 - 2*E2) +
	+ 4*E4*E2;

p = perm.poly(4)
# sort.pm(pow.pm(p, 3))
eval.pm(p, x)^3

pow = c(3,3)
p3 = perm.poly(4, p=pow)
eval.pm(p3, x)

p321 = perm3(4, p=c(3,2,1))
p321 = as.data.frame(p321)
names(p321) = paste0("x", seq(4));
p321$coeff = 1;
eval.pm(p321, x)
(-3*E3^2 + E3*E2*S - 3*E4*S^2 + 6*E4*E2 - 2*E4*E2)

### Debug:
x = sqrt(2:5)
x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4];
S = sum(x)
E4 = prod(x)
E3 = prod(x)*sum(1/x)
E2 = eval.pm(perm.poly(4, c(1,1)), x)


### Test
x1^2 + b*x2 # - R
x2^2 + b*x3 # - R
x3^2 + b*x4 # - R
x4^2 + b*x1 # - R


### Classic Polynomial:
b^3*x3 = b^2*R - (x1^2 - R)^2
b^7*x4 = b^6*R - ((x1^2 - R)^2 - b^2*R)^2
(((x1^2 - R)^2 - b^2*R)^2 - b^6*R)^2 + b^15*x1 - b^14*R # = 0

p1 = list(
	x = c(2,0),
	b = c(0,0),
	R = c(0,1),
	coeff = c(1,-1)
)
bR.gen = function(pb, pR=1) list(b = pb, R = pR, coeff = 1)
bx.gen = function(pb, px=1) list(b = pb, x = px, coeff = 1)
p1 = pow.pm(p1, 2)
p1 = diff.pm(p1, bR.gen(2, pR=1))
p1 = pow.pm(p1, 2)
p1 = diff.pm(p1, bR.gen(6, pR=1))
p1 = pow.pm(p1, 2)
p1 = diff.pm(p1, bR.gen(14, pR=1))
p1 = add.pm(p1, bx.gen(15, px=1))
p1 = sort.pm(p1, sort.coeff=c(4,2,3,1), xn="x")
p1

print.p(p1)

x^16 - 8*R*x^14 - 4*b^2*R*x^12 + 28*R^2*x^12 + 24*b^2*R^2*x^10 - 56*R^3*x^10 +
	- 2*b^6*R*x^8 + 6*b^4*R^2*x^8 - 60*b^2*R^3*x^8 + 70*R^4*x^8 +
	+ 8*b^6*R^2*x^6 - 24*b^4*R^3*x^6 + 80*b^2*R^4*x^6 - 56*R^5*x^6 +
	+ 4*b^8*R^2*x^4 - 16*b^6*R^3*x^4 + 36*b^4*R^4*x^4 - 60*b^2*R^5*x^4 + 28*R^6*x^4 +
	- 8*b^8*R^3*x^2 + 16*b^6*R^4*x^2 - 24*b^4*R^5*x^2 + 24*b^2*R^6*x^2 - 8*R^7*x^2 +
	+ b^15*x - b^14*R + b^12*R^2 - 2*b^10*R^3 + 5*b^8*R^4 - 6*b^6*R^5 + 6*b^4*R^6 - 4*b^2*R^7 + R^8

(x^4 - 2*R*x^2 + b^3*x - b^2*R + R^2) # * P[12]
(x^12 - 6*R*x^10 - b^3*x^9 + 3*(5*R^2 - b^2*R)*x^8 + 4*b^3*R*x^7 +
	+ (b^6 + 12*b^2*R^2 - 20*R^3)*x^6 + b^3*(2*b^2*R - 6*R^2)*x^5 +
	+ (15*R^4 - 18*b^2*R^3 + 3*b^4*R^2 - 4*b^6*R)*x^4 +
	+ b^3*(4*R^3 - 4*b^2*R^2 - b^6)*x^3 +
	- (6*R^5 - 12*b^2*R^4 + 6*b^4*R^3 - 5*b^6*R^2 + b^8*R)*x^2 +
	- b^3*(R^4 - 2*b^2*R^3 + b^4*R^2 - 2*b^6*R)*x +
	+ b^12 + 2*b^8*R^2 - 3*b^6*R^3 + 3*b^4*R^4 - 3*b^2*R^5 + R^6)


p4 = list(
	x = c(4,2,1,0,0),
	b = c(0,0,3,2,0),
	R = c(0,1,0,1,2),
	coeff = c(1,-2,1,-1,1)
)
p12 = list(
	x = c(12,10, 9, 8, 8, 7, 6, 6, 6, 5, 5, 4, 4, 4, 4,
		3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
	b = c( 0, 0, 3, 2, 0, 3, 6, 2, 0, 5, 3, 6, 4, 2, 0,
		9, 5, 3, 8, 6, 4, 2, 0, 9, 7, 5, 3,12, 8, 6, 4, 2, 0),
	R = c( 0, 1, 0, 1, 2, 1, 0, 2, 3, 1, 2, 1, 2, 3, 4,
		0, 2, 3, 1, 2, 3, 4, 5, 1, 2, 3, 4, 0, 2, 3, 4, 5, 6),
	coeff = c(1,-6,-1,-3,15, 4, 1,12,-20, 2,-6,-4, 3,-18,15,
		-1,-4, 4,-1, 5,-6,12,-6, 2, -1, 2, -1, 1, 2,-3, 3,-3, 1)
)
pR = diff.pm(p1, mult.pm(p12, p4))
pR = sort.pm(pR, sort.coeff=c(4,2,3,1), xn="x")
pR
print.p(as.data.frame(p12))

xi.f = function(x, R, b, n=2) {
	(R - x^n) / b[1];
}
solve.Simple.Classic.S4P2 = function(R, b, debug=TRUE) {
	coeff = c(1, 0, - 6*R, - b^3, 3*(5*R^2 - b^2*R), 4*b^3*R,
		(b^6 + 12*b^2*R^2 - 20*R^3), b^3*(2*b^2*R - 6*R^2),
		(15*R^4 - 18*b^2*R^3 + 3*b^4*R^2 - 4*b^6*R), # x^4
		b^3*(4*R^3 - 4*b^2*R^2 - b^6), # x^3
		- (6*R^5 - 12*b^2*R^4 + 6*b^4*R^3 - 5*b^6*R^2 + b^8*R), # x^2
		- b^3*(R^4 - 2*b^2*R^3 + b^4*R^2 - 2*b^6*R),
		b^12 + 2*b^8*R^2 - 3*b^6*R^3 + 3*b^4*R^4 - 3*b^2*R^5 + R^6
	)
	x1 = roots(coeff);
	x2 = xi.f(x1, R, b, n=2);
	x3 = xi.f(x2, R, b, n=2);
	x4 = xi.f(x3, R, b, n=2);
	sol = cbind(x1=x1, x2=x2, x3=x3, x4=x4);
	return(sol);
}

###
R = -1
b = 3
### distinct:
sol = solve.Simple.Classic.S4P2(R, b);
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];
S = x1 + x2 + x3 + x4;

### equal:
x = roots(c(1, 0, - 2*R, b^3, R^2 - b^2*R))
x1 = x;
x2 = xi.f(x1, R, b); x3 = xi.f(x2, R, b); x4 = xi.f(x3, R, b);

### E2_2
sapply(seq(12), function(id) sum(apply(perm2(4), 1, function(id2) prod(sol[id, id2 !=0]^2))))

### E3: Simplification (but NO progress on further simplification)
87*b^16 + 57*R*b^14 + 2987*R^2*b^12 - 14958*R^3*b^10 + 16796*R^4*b^8 + 1052*R^5*b^6 - 11312*R^6*b^4 +
	+ 6464*R^7*b^2 - 1280*R^8
pE3S2 = data.frame(
	b = seq(16, 0, by=-2), R = 0:8, coeff = c(87, 57, 2987, - 14958, 16796, 1052, - 11312, 6464, - 1280) )
pE3DS2 = data.frame(
	b = seq(13, 1, by=-2), R = 0:6, coeff = c(28, 262, -2450, 4616, -2596, -32, 448) ) # *-1

eval.pm(pE3S2, c(b, R))
eval.pm(pE3DS2, c(b, R))


#########################
#########################

###############
### Order 3 ###
###############

### x[i]^3 + b*x[i+1] = R

### Solution:

### Case: all x[i] different;
# - root structure: P[4] o P[18];
### Special Sub-Case:
# x1 = Conj(x3); x2 = Conj(x4);
# - 8 roots;
# - properties:
#  (x1^2+x3^2+x1*x3)*(x2^2+x4^2+x2*x4) = - b^2; [valid for all roots]
#  (3*Re(x1)^2 - Im(x1)^2)*(3*Re(x2)^2 - Im(x2)^2) = - b^2;
### Remaining roots: 64;

### Sum =>
S^3 - 3*E2*S + 3*E3 + b*S - 4*R # = 0

### Square => Sum =>
(x1^6 + x2^6 + x3^6 + x4^6) - b^2*(S^2 - 2*E2) + 2*b*R*S - 4*R^2 # = 0
S^6 - 6*E2*S^4 + 6*E3*S^3 + 9*E2^2*S^2 - 6*E4*S^2 - 12*E2*E3*S - 2*E2^3 + 3*E3^2 + 6*E2*E4 +
	- b^2*(S^2 - 2*E2) + 2*b*R*S - 4*R^2 # = 0


### Alternatives:
# - unfortunately NOT simpler;
### other Square => Sum =>
(x1^6 + x2^6 + x3^6 + x4^6) - 2*R*(x1^3 + x2^3 + x3^3 + x4^3) - b^2*(S^2 - 2*E2) + 4*R^2 # = 0
### Diff(Sq1 - Sq2) =>
2*R*(x1^3 + x2^3 + x3^3 + x4^3) + b^2*(S^2 - 2*E2) - b^2*(S^2 - 2*E2) + 2*b*R*S - 8*R^2 # = 0
2*R*(S^3 - 3*E2*S + 3*E3) + 2*b*R*S - 8*R^2 # = 0 # redundancy!
### Alternative:
### Diff =>
# PROD((x1^2+x2^2+x1*x2)) + b^6 # = 0
(x1^6*x2^4*x3*x4 + ...) + (x1^6*x2^4*x3^2 + ...) + (x1^6*x2^3*x3^3 + ...) +
	+ 2*(x1^6*x2^3*x3^2*x4 + ...) + 3*(x1^6*x2^2*x3^2*x4^2 + ...) +
	+ (x1^5*x2^5*x3*x4 + ...) + (x1^5*x2^5*x3^2 + ...) +
	+ 2*(x1^5*x2^4*x3^3 + ...) + 4*(x1^5*x2^4*x3^2*x4 + ...) +
	+ 5*(x1^5*x2^3*x3^3*x4 + ...) + 7*(x1^5*x2^3*x3^2*x4^2 + ...) +
	+ 3*(x1^4*x2^4*x3^4 + ...) + 7*(x1^4*x2^4*x3^3*x4 + ...) +
	+ 10*(x1^4*x2^4*x3^2*x4^2 + ...) + 12*(x1^4*x2^3*x3^3*x4^2 + ...) +
	+ 15*E4^3 + b^6 # = 0
# TODO;

### Eq 3:
# x1^3 = R - b*x2 => Prod =>
(x1*x2*x3*x4)^3 - R^4 + b*S*R^3 - b^2*E2*R^2 + b^3*E3*R - b^4*E4 # = 0
E4^3 - b^4*E4 + b*R^3*S - b^2*R^2*E2 + b^3*R*E3 - R^4 # = 0

### Eq 4:
# b*x2 = R - x1^3 => Prod =>
b^4*E4 - R^4 + R^3*(x1^3+x2^3+x3^3+x4^3) - R^2*E2_3 + R*E3_3 - E4^3 # = 0
b^4*E4 + R^3*(S^3 - 3*E2*S + 3*E3) - R^2*(E2^3 + 3*E3^2 - 3*E3*E2*S + 3*E4*S^2 - 3*E2*E4) +
	+ R*(E3^3 - 3*E4*E3*E2 + 3*E4^2*S) - E4^3 - R^4 # = 0

### TODO:
# - solve;
# - solve special Sub-Case;

### Special Sub-Case:
# x1 = Conj(x3); x2 = Conj(x4);
r1 = Re(x1); r2 = Re(x2); z1 = Im(x1); z2 = Im(x2);
### Sum: Eq[i] + Eq[i+2]
r1^3 - 3*r1*z1^2 + b*r2 - R # = 0
r2^3 - 3*r2*z2^2 + b*r1 - R # = 0
### Diff: Eq[i] - Eq[i+2] # Anti-Symmetric!
z1*(3*r1^2 - z1^2) + b*z2 # = 0
z2*(3*r2^2 - z2^2) - b*z1 # = 0
# if (z1, z2) is a solution => so is: (-z1, -z2);

### Diff(r1*Sum) =>
(r1-r2)*(r1^3 + r2^3 + r1*r2*(r1+r2) - R) - 3*(r1*z1 - r2*z2)*(r1*z1 + r2*z2) # = 0
### Diff(r2*Sum) =>
(r1-r2)*(r1*r2*(r1+r2) - b*(r1+r2) + R) - 3*r1*r2*(z1-z2)*(z1+z2) # = 0

########

### Test
x1^3 + b*x2 # - R
x2^3 + b*x3 # - R
x3^3 + b*x4 # - R
x4^3 + b*x1 # - R

### Classic Solver:
solve.Simple.Classic.S4P3 = function(R, b, debug=TRUE) {
	# coeff: generate using print.coeff(p1);
	# from Classic Polynomial! [see below]
	coeff = coeff.S4P3.Classic(R, b);
	x1 = roots(coeff);
	n = 3;
	x2 = xi.f(x1, R, b, n=n);
	x3 = xi.f(x2, R, b, n=n);
	x4 = xi.f(x3, R, b, n=n);
	sol = cbind(x1=x1, x2=x2, x3=x3, x4=x4);
	return(sol);
}

### Test
R = -1
b = 2
sol = solve.Simple.Classic.S4P3(R, b)
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];

E = do.call(rbind, lapply(seq(nrow(sol)), function(id) debug.E(sol[id,])));
E$S = round0(E$S);
S = E$S; E2 = E$E2; E3 = E$E3; E4 = E$E4;

test.S4.Simple(sol, R, b, n=3)

# S = sort(apply(sol, 1, sum))[ ! duplicated(round(sort(apply(sol, 1, sum)), 5))]
round0.p(poly.calc(S))


### Classic Polynomial:
n = 3
b^(n+1)*x3 = b^n*R - (R - x1^n)^n
b^(n^2+n+1)*x4 = b^(n^2+n)*R - (b^n*R - (R - x1^n)^n)^n
(b^(n^2+n)*R - (b^n*R - (R - x1^n)^n)^n)^n +
	+ b^((n^4-1)/(n-1))*x1 - b^(n^3+n^2+n)*R # = 0

### generate Classic Polynomial:
bR.gen = function(pb, pR=1) list(b = pb, R = pR, coeff = 1)
bx.gen = function(pb, px=1) list(b = pb, x = px, coeff = 1)
S4P3.gen.Classic.P81 = function() {
	n = 3
	p1 = list(
		x = c(n,0),
		b = c(0,0),
		R = c(0,1),
		coeff = c(-1, 1)
	)
	p1 = pow.pm(p1, n)
	p1 = diff.pm(bR.gen(n, pR=1), p1)
	p1 = pow.pm(p1, n)
	p1 = diff.pm(bR.gen(n^2+n, pR=1), p1)
	p1 = pow.pm(p1, n)
	p1 = diff.pm(p1, bR.gen((n^4-1)/(n-1) - 1, pR=1))
	p1 = add.pm(p1, bx.gen((n^4-1)/(n-1), px=1))
	p1 = sort.pm(p1, sort.coeff=c(4,2,3,1), xn="x")
	rownames(p1) = seq(nrow(p1))
	p1 = mult.sc.pm(p1, -1)
	return(p1);
}
p1 = S4P3.gen.Classic.P81();

# print with line breaks:
pprint.m = matrix(c(
	 1, 12,	13, 23,
	24, 28,	29, 39,
	40, 51,	52, 58,
	59, 65,	66, 72,
	73, 81,	82, 91,
	92, 99,	100, 109,
	110, 119,	120, 129,
	130, 140,	141, 151,
	152, 162,	163, 176
), nrow=2)

apply(pprint.m, 2, function(rw.id) print.p(p1[rw.id[1]:rw.id[2], ], leading="x"))

### Coefficients
print.coeff(p1)

### P[3] * P[6] * P[54]
# P[9]:
# (x^3 + b*x - R) * (x^6 - b*x^4 - 2*R*x^3 + b^2*x^2 + b*R*x + R^2 - b^3)
# x^9 - 3*R*x^6 + 3*R^2*x^3 - b^4*x - R^3 + b^3*R
pDiv = list(
	x = c(9, 6, 3, 1, 0, 0),
	b = c(0, 0, 0, 4, 0, 3),
	R = c(0, 1, 2, 0, 3, 1),
	coeff = c(1,-3, 3,-1, -1, 1)
)
pR = div.pm(p1, pDiv)
pR = sort.pm(pR$Rez, sort.coeff=c(4,2,3,1), xn="x")
print.coeff(pR)

# significant numerical error!
eval.pm(pR$Rez, c(b, R, sol[10,1]))


### Deriving Diff:
# PROD((x1^2+x2^2+x1*x2))
pp1 = data.frame(
	x1 = c(2, 0, 1),
	x2 = c(0, 2, 1),
	coeff = c(1, 1, 1)
)
pAll = perm2.pm(pp1, paste0("x", 1:4))
pR = mult.all.pm(pAll);
pR = sort.pm(pR, sort.coeff=c(2,3,1,5,4), xn="x1")
pR

### E polynomials
### E2_3
E2^3 - 3*(-3*E3^2 + E3*E2*S - 3*E4*S^2 + 4*E4*E2) +
	- 6*E4*(S^2 - 2*E2) - 15*E2*E4 - 6*(E3^2 - 2*E2*E4)
E2^3 + 3*E3^2 - 3*E3*E2*S + 3*E4*S^2 - 3*E2*E4

p33  = perm.poly(4, p=c(3,3))
p321 = perm.poly(4, p=c(3,2,1))
diff.lpm(pow.pm(perm.poly(4), 3), list(p33, mult.sc.pm(p321, 3)))

### Test
E = do.call(rbind, lapply(seq(nrow(sol)), function(id) debug.E(sol[id,])));
E$S = round0(E$S);
S = E$S; E2 = E$E2; E3 = E$E3; E4 = E$E4;

eval.pm(p33, x) # E2_3: seems correct

### E3_3
E3^3 - 3*E4*E3*E2 + 3*E4^2*S
diff.lpm(pow.pm(perm.poly(4, c(1,1,1)), 3), # E3^3
	list(perm.poly(4, c(3,3,3)),
		mult.pm(Eprod.pm(4, 2), Esum.pm(4), sc=-3), # + 3*E4^2*S
		mult.all.pm( # - 3*E4*E3*E2
			list(Eprod.pm(4, 1), perm.poly(4, c(1,1,1)), perm.poly(4, c(1,1)), sc=3))))


#############################
#############################

#############################
### Mixed Terms: Type V2a ###
#############################

###############
### Order 2 ###
###############

### x1^2 + b*x1*x2 = R

### Solution:

### Case: all x[i] different;

### Eq 1: Sum(Prod(x2*x3*x4*...)) =>
E4*S + b*E4*S - R*E3 # = 0
(b+1)*E4*S - R*E3 # = 0

### Eq 2:
### Sum =>
S^2 - 2*E2 + b*(x1*x2+x2*x3+x3*x4+x4*x1) - 4*R # = 0
# P2a = (x1*x2+x2*x3+x3*x4+x4*x1); P2a = round0(P2a);
# b*P2a = 2*E2 - S^2 + 4*R;

### x1^2 = R - b*x1*x2 => Prod =>
E4^2 - R^4 + b*R^3*P2a - b^2*R^2*(P3a + 2*E4) + b^3*R*E4*P2a - b^4*E4^2 # = 0
# Eq 2:
(b^4 - 1)*(b^2 + 2)*E4^2 - 2*b^2*(b^2 + 4)*R^2*E4 +
	- 2*b^2*(b^2 + 2)*R*E4*E2 + b^2*(b^2 + 2)*R*E4*S^2 +
	+ 4*R^2*E2^2 + 2*(b^2 + 6)*R^3*E2 - 4*R^2*E2*S^2 +
	+ (b^2+10)*R^4 - (b^2+6)*R^3*S^2 + R^2*S^4 # = 0

### Eq 3:
### b*x1*x2 = R - x1^2 => Prod =>
# [see Simple version: E4 => E4^2]
E4^2 - b^4*E4^2 + 2*R^2*E4 + 2*R*E2*E4 - R*E3^2 - 2*R^2*E3*S +
	+ 2*R^3*E2 + R^2*E2^2 - R^3*S^2 + R^4 # = 0

### Eq 4:
### Alternative 1:
### Sum(x3^2*x4^2*...) =>
(E3^2 - 2*E4*E2) + b*E4*P2a - R*(P2a^2 - 2*P3a - 4*E4) # = 0
b^2*E3^2 - b^2*E4*S^2 + 8*b^2*R*E4 - R*(2*E2 - S^2 + 4*R)^2 + 2*b^2*R*P3a # = 0
b^2*E3^2 - b^2*E4*S^2 + 8*b^2*R*E4 - 4*R*E2^2 + 4*R*E2*S^2 - 16*R^2*E2 +
	- R*(S^4 - 8*R*S^2 + 16*R^2) + 2*b^2*R*P3a # = 0
(b^2 + 2)*E3^2 - (b^2 + 2)*E4*S^2 + 8*(b^2 + 1)*R*E4 +
	- 4*R*E2^2 + 4*R*E2*S^2 - 8*R^2*E2 +
	- R*S^4 + 4*R^2*S^2 - 8*R^3 # = 0


#######
### Eq:
S^3 - R*(b^4 - 2*b^3 + 4*b^2 - 4*b + 4)*S

### Auxiliary eqs:
### Special Case: S = 0
E3 = 0; E2 = -2*R; E4 = R^2 / (b^2 + 1);
### Case: S != 0
E2 = -b*R*(b^2 - b + 2);
E3 = - (b+1) * R * S;
E4 = R*E3 / ((b+1)*S); # - R^2;


### [Eq 4] Alternative 2:
### Sum(x4*...) =>
(x1^2*x4 + x2^2*x1 + x3^2*x2 + x4^2*x3) + b*E3 - R*S # = 0
### Sum(x3*x4*...) =>
P3b = x1^2*x3*x4 + x2^2*x1*x4 + x3^2*x1*x2 + x4^2*x2*x3;
P3b + 4*b*E4 - R*P2a # = 0
# [...]

############
### Workout:
### Half-Elementary Polynomials
P3a = x1*x2^2*x3 + x1^2*x2*x4 + x2*x3^2*x4 + x1*x3*x4^2
p2a = perm.poly(4, c(1,1))[c(1,3,4,6),]
pR = data.frame(x1=0, coeff=1);
for(nr in seq(nrow(p2a))) {
	p1 = p2a[nr,];
	p1 = rbind(p1, 0); p1$coeff[2] = -1;
	p1$R = c(0, 1);
	pR = mult.pm(pR, p1);
}
pR = sort.pm(pR, c(4,2,3,1), xn="R")
pR

### P3a:
### x1^2 = R - b*x1*x2 => Sum(Prod(2 eqs)) =>
b^2*P3a - 2*b*R*P2a + 4*R^2 - (P2a^2 - 2*P3a - 4*E4) # = 0
(b^2 + 2)*P3a - 2*R*(2*E2 - S^2 + 4*R) - P2a^2 + 4*E4 + 4*R^2
b^2*(b^2 + 2)*P3a - 2*b^2*R*(2*E2 - S^2 + 4*R) - (2*E2 - S^2 + 4*R)^2 + 4*b^2*E4 + 4*b^2*R^2
b^2*(b^2 + 2)*P3a - 4*E2^2 + 4*E2*S^2 - 4*(b^2+4)*R*E2 +
	+ 4*b^2*E4 - S^4 + 2*(b^2+4)*R*S^2 - 4*(b^2+4)*R^2


### Test
x1^2 + b*x1*x2 # - R
x2^2 + b*x2*x3 # - R
x3^2 + b*x3*x4 # - R
x4^2 + b*x4*x1 # - R

### Debug:
solve.S4P2V2.classic = function(R, b, debug=FALSE) {
	x1 = roots(coeff.S4P2V2(R, b)); # see below for coeff.S4P2V2();
	x2 = xip.f(x1, R, b, p=1);
	x3 = xip.f(x2, R, b, p=1);
	x4 = xip.f(x3, R, b, p=1);
	sol = cbind(x1, x2, x3, x4);
	return(sol);
}
#
R = -1; b = 2;
sol = solve.S4P2V2.classic(R, b);
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];
E = do.call(rbind, lapply(seq(nrow(sol)), function(id) debug.E(sol[id,])));
E$S = round0(E$S);
S = E$S; E2 = E$E2; E3 = E$E3; E4 = E$E4;
E2 = -b*R*(b^2 - b + 2);

### Classic Polynomial:
b^n*x1^(n-1)*(R - x1^n)*x3 = b^n*x1^n*R - (R - x1^n)^n;
b^3*x1*(R - x1^n)*(b^2*x1^2*R - (R - x1^n)^2)*x4 =
	b^4*x1^2*(R - x1^n)^2*R - (b^2*x1^2*R - (R - x1^n)^2)^2
(b^4*R*x1^2*(R - x1^n)^2 - (b^2*R*x1^2 - (R - x1^n)^2)^2)^2 +
	+ b^4*x1^2*(R - x1^n)*(b^2*R*x1^2 - (R - x1^n)^2)*
		(b^4*R*x1^2*(R - x1^n)^2 - (b^2*R*x1^2 - (R - x1^n)^2)^2) +
	- b^6*R*x1^2*(R - x1^n)^2*(b^2*R*x1^2 - (R - x1^n)^2)^2
# pRx = (R - x1^n); pRx2 = pRx^2;
# pBd = (b^2*R*x1^2 - pRx2); pBd2 = pBd^2;
# pBB = (b^4*R*x1^2*pRx2 - pBd2);
pBB^2 + b^4*x1^2*pRx*pBd*pBB - b^6*R*x1^2*pRx2*pBd2

n = 2
pRx = data.frame(
	x = c(n,0),
	b = c(0,0),
	R = c(0,1),
	coeff = c(-1, 1)
)
# (b^2-1)*x^4 + 2*R*x^2 - R^2
pDiv = data.frame(
	x = c(4, 4, 2, 0),
	b = c(2, 0, 0, 0),
	R = c(0, 0, 1, 2),
	coeff = c(1,-1, 2,-1)
)
bRx.gen = function(pb, px=n, pR=1) list(x = px, b = pb, R = pR, coeff = 1)
pRx2 = pow.pm(pRx, 2);
pBd  = diff.pm(bRx.gen(2, px=2, pR=1), pRx2);
pBd2 = pow.pm(pBd, 2);
pBB  = diff.pm(mult.pm(bRx.gen(4, px=2, pR=1), pRx2), pBd2);
#
p1 = mult.all.pm(list(pRx, pBd, pBB, bRx.gen(4, px=2, pR=0)))
p2 = mult.all.pm(list(pRx2, pBd2, bRx.gen(6, px=2, pR=1)))
p1 = diff.pm(p1, p2);
p1 = add.pm(p1, pow.pm(pBB, 2));
p1 = mult.sc.pm(p1, -1);
p1 = sort.pm(p1, sort.coeff=c(4,2,3,1), xn="x")
p1

print.p(p1)
p12 = div.pm(p1, pDiv, by="b")
print.coeff(p12$Rez)

### P[16] = P[2] * P[2] * P[12]
(b^4-1)*x^16 + (- b^8*R - 2*b^6*R - 5*b^4*R + 4*b^2*R + 8*R)*x^14 +
	+ (b^10*R^2 + 5*b^8*R^2 + 5*b^6*R^2 + 3*b^4*R^2 - 24*b^2*R^2 - 28*R^2)*x^12 +
	+ (- 3*b^10*R^3 - 5*b^8*R^3 + 5*b^6*R^3 + 19*b^4*R^3 + 60*b^2*R^3 + 56*R^3)*x^10 +
	+ (2*b^10*R^4 - 4*b^8*R^4 - 22*b^6*R^4 - 41*b^4*R^4 - 80*b^2*R^4 - 70*R^4)*x^8 +
	+ (6*b^8*R^5 + 20*b^6*R^5 + 33*b^4*R^5 + 60*b^2*R^5 + 56*R^5)*x^6 +
	+ (- 2*b^8*R^6 - 7*b^6*R^6 - 11*b^4*R^6 - 24*b^2*R^6 - 28*R^6)*x^4 +
	+ (b^6*R^7 + b^4*R^7 + 4*b^2*R^7 + 8*R^7)*x^2 - R^8
# P[2] * P[2] * P[12]
((b^2-1)*x^4 + 2*R*x^2 - R^2) * P[12]
((b+1)*x^2 - R) * ((b-1)*x^2 + R) * P[12]

### Classic solver:
coeff.S4P2V2 = function(R, b) {
	# P[12]
	coeff = c(b^2 + 1, 0,
		- R*(b^2+1)*(b^4 + 2*b^2 + 6), 0,
		15*R^2 + 22*b^2*R^2 + 13*b^4*R^2 + 6*b^6*R^2 + b^8*R^2, 0,
		- 20*R^3 - 28*b^2*R^3 - 18*b^4*R^3 - 10*b^6*R^3 - 3*b^8*R^3, 0,
		15*R^4 + 17*b^2*R^4 + 9*b^4*R^4 + 5*b^6*R^4 + 2*b^8*R^4, 0,
		- 6*R^5 - 4*b^2*R^5 - b^4*R^5 - b^6*R^5, 0,
		R^6);
	return(coeff);
}
coeff.S4P2V2_P16 = function(R, b) {
	# P[16]
	coeff = c(b^4 - 1, 0,
		- b^8*R - 2*b^6*R - 5*b^4*R + 4*b^2*R + 8*R, 0,
		b^10*R^2 + 5*b^8*R^2 + 5*b^6*R^2 + 3*b^4*R^2 - 24*b^2*R^2 - 28*R^2, 0,
		- 3*b^10*R^3 - 5*b^8*R^3 + 5*b^6*R^3 + 19*b^4*R^3 + 60*b^2*R^3 + 56*R^3, 0,
		2*b^10*R^4 - 4*b^8*R^4 - 22*b^6*R^4 - 41*b^4*R^4 - 80*b^2*R^4 - 70*R^4, 0,
		6*b^8*R^5 + 20*b^6*R^5 + 33*b^4*R^5 + 60*b^2*R^5 + 56*R^5, 0,
		- 2*b^8*R^6 - 7*b^6*R^6 - 11*b^4*R^6 - 24*b^2*R^6 - 28*R^6, 0,
		b^6*R^7 + b^4*R^7 + 4*b^2*R^7 + 8*R^7, 0, - R^8)
	return(coeff);
}


### Auxiliary Eqs:
### Special Case: S = 0; E3 = 0;
- 2*(b^4 + 3*b^2 - 2)*R*E4 +
	- 2*(b^4 + b^2 - 2)*E4*E2 + (b^2+6)*R*E2^2 + 4*(b^2 + 4)*R^2*E2 + 2*(b^2+6)*R^3 # = 0
2*(b^2 + 1)*R*E4 - R*E2^2 - 2*R^2*E2 - 2*R^3 # = 0
# =>
(b^2-1)*E4*E2 - 4*R*E4 - R^2*E2
# => E2 = -2*R; E4 = R^2 / (b^2 + 1);


pEEq1 = data.frame(
	E4 = c(1, 1, 1,   1, 1, 1,  0, 0,   0, 0,   0, 0),
	E2 = c(0, 0, 0,   1, 1, 1,  2, 2,   1, 1,   0, 0),
	b  = c(4, 2, 0,   4, 2, 0,  2, 0,   2, 0,   2, 0),
	R  = c(2, 2, 2,   1, 1, 1,  2, 2,   3, 3,   4, 4) - 1,
	coeff = c(-2,-6,4, -2,-2,4,  1, 6,  4,16,   2,12)
)
pEEq2 = data.frame(
	E4 = c(1, 1,   0, 0, 0),
	E2 = c(0, 0,   2, 1, 0),
	b  = c(2, 0,   0, 0, 0),
	R  = c(2, 2,   2, 3, 4) - 2,
	coeff = c(2, 2,  -1,-2,-2)
)
pb.gen = function(pB, pCoeff, pR=0) {
	data.frame(
		b = pB, R = pR, coeff = pCoeff
	)
}
round0(sapply(seq_along(E4), function(id) eval.pm(pEEq1, c(E4[id], E2[id], b, R))))
round0(sapply(seq_along(E4), function(id) eval.pm(pEEq2, c(E4[id], E2[id], b, R))))

pEEq3 = add.pm(pEEq1, mult.pm(pEEq2, pb.gen(c(2,0), c(1,6), pR=1)))


### Full Case:
### Eq 1:
(b+1)*E4*S - R*E3 # = 0
### Eq 3:
E4^2 - b^4*E4^2 + 2*R^2*E4 + 2*R*E2*E4 - R*E3^2 - 2*R^2*E3*S +
	+ 2*R^3*E2 + R^2*E2^2 - R^3*S^2 + R^4 # = 0
### Eq S:
S^3 - R*(b^4 - 2*b^3 + 4*b^2 - 4*b + 4)*S


pEq3 = data.frame(
	E4 = c(2, 2, 1, 1, 0, 0, 0, 0, 0, 0),
	E3 = c(0, 0, 0, 0, 2, 1, 0, 0, 0, 0),
	E2 = c(0, 0, 0, 1, 0, 0, 1, 2, 0, 0),
	S  = c(0, 0, 0, 0, 0, 1, 0, 0, 2, 0),
	b  = c(0, 4, 0, 0, 0, 0, 0, 0, 0, 0),
	R  = c(0, 0, 2, 1, 1, 2, 3, 2, 3, 4),
	coeff = c(1,-1, 2, 2, -1, -2, 2, 1, -1, 1)
)
pEq1   = data.frame(E3=1, R=1, coeff=1); # E4 vs E3
pEq1fr = data.frame(S=c(1,1), b=c(1,0), coeff=c(1,1));
pDiv = data.frame(
	S = c(2, 0, 0, 0, 0, 0),
	b = c(0, 4, 3, 2, 1, 0),
	R = c(0, 1, 1, 1, 1, 1),
	coeff = c(1,-1, 2, -4, 4, -4)
)
pEq3Coeff = data.frame(b=5:0, coeff=c(1,-1,3,-1,1,3));
#
pEq3r = replace.fr.pm(pEq3, pEq1, pEq1fr, "E4", pow=1)
lP3 = div.pm(pEq3r, pDiv, by="S")
# lP3$Rem
lP3 = div.pm(lP3$Rem, data.frame(b=c(1,0), coeff=c(1,1)), by="b")
lP3 = lP3$Rez; lP3$R = lP3$R - 2; lP3$coeff = - lP3$coeff;
id = order( - lP3$E3, - lP3$E2, - lP3$b); lP3 = lP3[id,];
print.p(lP3, "E3")
lP3

# Debug:
eval.pm(pEq3, c(E4[10], E3[10], E2[10], S[10], b, R))

(3 + b - b^2 + 3*b^3 - b^4 + b^5)*E3^2 +
	- 2*(E2 - 3*R - 2*R*b^3 + R*b^4 - R*b^5)*S*E3 +
	- R*(b^5 - b^4 + 2*b^3 + 4)*(E2^2 + 2*R*E2) +
	+ (b^9 - 3*b^8 + 8*b^7 - 12*b^6 + 15*b^5 - 7*b^4 - 2*b^3 + 16*b^2 - 16*b + 12)*R^3


### Eq 4:
(b^2 + 2)*E3^2 - (b^2 + 2)*E4*S^2 + 8*(b^2 + 1)*R*E4 +
	- 4*R*E2^2 + 4*R*E2*S^2 - 8*R^2*E2 +
	- R*S^4 + 4*R^2*S^2 - 8*R^3 # = 0

pEq4 = data.frame(
	E4 = c(0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
	E3 = c(2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
	E2 = c(0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0),
	S  = c(0, 0, 2, 2, 0, 0, 0, 2, 0, 4, 2, 0),
	b  = c(2, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0),
	R  = c(0, 0, 0, 0, 1, 1, 1, 1, 2, 1, 2, 3),
	coeff = c(1, 2,-1,-2, 8, 8,-4, 4,-8,-1, 4,-8)
)
pEq4Coeff = data.frame(S=1, b=c(3,2,1,0), coeff=c(1,1,2,2));
#
pEq4r = replace.fr.pm(pEq4, pEq1, pEq1fr, "E4", pow=1)
lP4 = div.pm(pEq4r, pDiv, by="S")
lP4 = lP4$Rem;
id = order( - lP4$E3, - lP4$E2, -lP4$b); lP4 = lP4[id,];
print.p(lP4, "E3")
# lP4

pR = diff.pm(mult.pm(lP3, pEq4Coeff), mult.pm(lP4, pEq3Coeff))
lP2r = div.pm(pR, pDiv, by="S")
lP2r = lP2r$Rem;
id = order( - lP2r$E3, - lP2r$E2, -lP2r$b); lP2r = lP2r[id,];
print.p(lP2r, "E3")


# Debug:
eval.pm(pEq4, c(E4[10], E3[10], E2[10], S[10], b, R))

(b+1)*(b^2 + 2)*E3^2*S +
	+ R^2 * (- b^6 + 2*b^5 - 6*b^4 + 8*b^3 - 4*b^2 + 8*b)*E3 +
	- 4*R*(b+1)*E2^2*S +
	+ 4*R^2 * (b+1)*(b^4 - 2*b^3 + 4*b^2 - 4*b + 2)*E2*S +
	- R^3 * (b+1)*(b^8 - 4*b^7 + 12*b^6 - 24*b^5 + 36*b^4 - 40*b^3 + 32*b^2 - 16*b + 8)*S
	# (b^9 - 3*b^8 + 8*b^7 - 12*b^6 + 12*b^5 - 4*b^4 - 8*b^3 + 16*b^2 - 8*b + 8)*S

((- 2*b^7 + 2*b^6 - 8*b^5 + 4*b^4 - 8*b^3 - 8*b^2 - 16)*E2 +
	+ R*(2*b^12 - 3*b^11 + 11*b^10 - 5*b^9 + 7*b^8 + 31*b^7 - 39*b^6 + 84*b^5 - 46*b^4 +
		+ 44*b^3 + 28*b^2 - 24*b + 48))*E3 +
	- (b^8 - b^6 + 2*b^5 - 6*b^4 + 4*b^2 - 8*b - 4)*E2^2*S +
	- (4*b^10 - 8*b^9 + 26*b^8 - 24*b^7 + 30*b^6 + 20*b^5 - 32*b^4 + 72*b^3 - 8*b^2 + 40)*R*E2*S +
	+ (b^14 - 4*b^13 + 15*b^12 - 32*b^11 + 59*b^10 - 68*b^9 + 56*b^8 + 12*b^7 - 71*b^6 + 134*b^5 - 86*b^4 +
		+ 48*b^3 + 44*b^2 - 24*b + 48)*R^2*S

### Eq 2:
(b^4 - 1)*(b^2 + 2)*E4^2 - 2*b^2*(b^2 + 4)*R^2*E4 +
	- 2*b^2*(b^2 + 2)*R*E4*E2 + b^2*(b^2 + 2)*R*E4*S^2 +
	+ 4*R^2*E2^2 + 2*(b^2 + 6)*R^3*E2 - 4*R^2*E2*S^2 +
	+ (b^2+10)*R^4 - (b^2+6)*R^3*S^2 + R^2*S^4 # = 0

pEq2 = data.frame(
	E4 = c(2, 2, 2, 2, 1, 1,   1, 1, 1, 1,   0, 0, 0, 0,   0, 0, 0, 0, 0),
	E2 = c(0, 0, 0, 0, 0, 0,   1, 1, 0 ,0,   2, 1, 1, 1,   0, 0, 0, 0, 0),
	S  = c(0 ,0, 0, 0, 0, 0,   0, 0, 2, 2,   0, 0, 0, 2,   0, 0, 2, 2, 4),
	b  = c(6, 4, 2, 0, 4, 2,   4, 2, 4, 2,   0, 2, 0, 0,   2, 0, 2, 0, 0),
	R  = c(0, 0, 0, 0, 2, 2,   1, 1, 1, 1,   2, 3, 3, 2,   4, 4, 3, 3, 2),
	coeff = c(1, 2,-1,-2,-2,-8,  -2,-4,1, 2,  4, 2,12,-4,  1,10,-1,-6, 1)
)
bDiv4 = data.frame(b=0:9, coeff=c(8,16,12,12,10,4,5,3,1,1))
bDiv3 = data.frame(b=0:8, coeff=c(8,8,4,8,2,2,3,0,1))

pEq2r = replace.fr.pm(pEq2, pEq1, pEq1fr, "E4", pow=1)
lP2 = div.pm(pEq2r, pDiv, by="S")
id = order( - lP2$Rem$E3, - lP2$Rem$E2, -lP2$Rem$b); lP2$Rem = lP2$Rem[id,];
print.p(lP2$Rem, "E3")
# lP$Rem

# lP2r: E3 vs E2
ncE3 = match("E3", names(lP2r))
isDivE3 = (lP2r[, ncE3] > 0);
E3fr = lP2r[isDivE3, -ncE3];
E3p  = lP2r[ ! isDivE3, -ncE3];
E3p$coeff = - E3p$coeff;
#
lP2 = replace.fr.pm(lP2$Rem, E3p, E3fr, "E3", pow=1)
lP2 = div.pm(lP2, pDiv, by="S")
lP2 = lP2$Rem;
lP2$S = lP2$S - min(lP2$S);
lP2$R = lP2$R - min(lP2$R);
lP2 = div.pm(lP2, bDiv3, by="b");
lP2 = lP2$Rez;
id = order( - lP2$E2, -lP2$b); lP2 = lP2[id,];
rownames(lP2) = seq(nrow(lP2));
lP2

# lP2 & lP4b are similar: redundancy or simplification method?
pEq4r = replace.fr.pm(lP4, E3p, E3fr, "E3", pow=1)
lP4b = div.pm(pEq4r, pDiv, by="S")
lP4b = lP4b$Rem;
lP4b$S = lP4b$S - min(lP4b$S);
lP4b$R = lP4b$R - min(lP4b$R);
id = order( - lP4b$E2, -lP4b$b); lP4b = lP4b[id,];
lP4b = div.pm(lP4b, bDiv3, by="b");
lP4b = lP4b$Rez;
id = order( - lP4b$E2, -lP4b$b); lP4b = lP4b[id,];
rownames(lP4b) = seq(nrow(lP4b));
lP4b
# print.p(lP4, "E2")

### simplification of E3:
# -b*(b^2 - b + 2)
pRepl = data.frame(b=c(3,2,1), coeff=c(-1,1,-2));
pDiv = add.pm(mult.pm(E3fr[E3fr$E2 == 1, c("b", "coeff")], pRepl),
	E3fr[E3fr$E2 == 0, c("b", "coeff")])
print.p(pDiv, "b")

pE3p = add.lpm(
	list(mult.pm(E3p[E3p$E2 == 2, c("b", "coeff")], pow.pm(pRepl, 2)),
	mult.pm(E3p[E3p$E2 == 1, c("b", "coeff")], pRepl),
	E3p[E3p$E2 == 0, c("b", "coeff")]))
pE3p$coeff = - pE3p$coeff; # use: - E3 !!!
print.p(pE3p, "b")

gcd.pm(pE3p, pDiv, by="b");
# (b+1) !!!

# Debug:
eval.pm(pEq2, c(E4[10],E2[10],S[10],b,R))
eval.pm(E3p, c(R,E2[10],S[10],b))

gcd.pm(lP2[lP2$E2 == 4, c("b", "coeff")], lP4b[lP4b$E2 == 4, c("b", "coeff")], by="b", div.sc=5974.576)
div.pm(lP2[lP2$E2 == 4, c("b", "coeff")], bDiv4, "b")
div.pm(lP4b[lP4b$E2 == 4, c("b", "coeff")], bDiv4, "b")

diff.pm(mult.pm(lP2, lP4b[lP4b$E2 == 4, c("b", "coeff")]), mult.pm(lP4b, lP2[lP2$E2 == 4, c("b", "coeff")]))

### [old]
getE3.old = function(S, E2, R, b) {
	# only for Case: S != 0;
	# E3 = 0 for S == 0; [this formula actually works as well]
	pDiv = ((2*b^7 - 2*b^6 + 8*b^5 - 4*b^4 + 8*b^3 + 8*b^2 + 16)*E2 +
		- R*(2*b^12 - 3*b^11 + 11*b^10 - 5*b^9 + 7*b^8 + 31*b^7 - 39*b^6 + 84*b^5 - 46*b^4 +
			+ 44*b^3 + 28*b^2 - 24*b + 48));
	#
	pE3 = - (b^8 - b^6 + 2*b^5 - 6*b^4 + 4*b^2 - 8*b - 4)*E2^2 +
		- (4*b^10 - 8*b^9 + 26*b^8 - 24*b^7 + 30*b^6 + 20*b^5 - 32*b^4 + 72*b^3 - 8*b^2 + 40)*R*E2 +
		+ (b^14 - 4*b^13 + 15*b^12 - 32*b^11 + 59*b^10 - 68*b^9 + 56*b^8 + 12*b^7 - 71*b^6 + 134*b^5 - 86*b^4 +
			+ 48*b^3 + 44*b^2 - 24*b + 48)*R^2;
	return(pE3 * S / pDiv);
}


################################
################################

################################
### Mixed Leading Terms: L2  ###
### Simple Chain: Type V1    ###
################################

###############
### Order 2 ###
###############

### x1^2*x2^2 + b*x3 = R

# x1^2*x2^2 + b*x3 = R
# x2^2*x3^2 + b*x4 = R
# x3^2*x4^2 + b*x1 = R
# x4^2*x1^2 + b*x2 = R

### Solution:

### Case: all x[i] different;

### Eq 1:
# x1^2*x2^2 = R - b*x3 => Prod =>
### Eq 1:
E4^4 - b^4*E4 + b^3*R*E3 - b^2*R^2*E2 + b*R^3*S - R^4 # = 0

### Eq 2:
### x1^2*x2^2 = R - b*x3 => Prod(Eq1, Eq 3) =>
E4^2 - R^2 + b*R*(x1+x3) - b^2*x1*x3 # = 0
E4^2 - R^2 + b*R*(x2+x4) - b^2*x2*x4 # = 0
### Sum =>
2*E4^2 - b^2*(x1*x3 + x2*x4) + b*R*S - 2*R^2 # = 0
# b^2*(x1*x3 + x2*x4) = 2*E4^2 + b*R*S - 2*R^2

### Sum =>
((x1*x2)^2+(x2*x3)^2+(x3*x4)^2+(x4*x1)^2) + b*S - 4*R # = 0
E2^2 - 2*E3*S + 2*E4 + b*S - 4*R - (2*E4^2 + b*R*S - 2*R^2)^2 / b^4 + 2*E4 # = 0
b^4*E2^2 - 2*b^4*E3*S + 4*b^4*E4 + b^5*S - 4*b^4*R - (2*E4^2 + b*R*S - 2*R^2)^2 # = 0
b^4*E2^2 - 2*b^4*E3*S + 4*b^4*E4 + b^5*S - 4*b^4*R +
	- (4*E4^4 + b^2*R^2*S^2 + 4*R^4 + 4*b*R*E4^2*S - 8*R^2*E4^2 - 4*b*R^3*S) # = 0
### Eq 2:
4*E4^4 + 4*b*R*E4^2*S - 8*R^2*E4^2 - 4*b^4*E4 - b^4*E2^2 + 2*b^4*E3*S +
	+ b^2*R^2*S^2 - b^5*S - 4*b*R^3*S + 4*b^4*R + 4*R^4 # = 0
### Eq 2 (simplified)
# Diff: Eq 2 - 4*Eq 1 =>
4*b*R*E4^2*S - 8*R^2*E4^2 - 4*b^3*R*E3 - b^4*E2^2 + 4*b^2*R^2*E2 + 2*b^4*E3*S +
	- 4*b*R^3*S - b^5*S + b^2*R^2*S^2 - 4*b*R^3*S + 4*b^4*R + 8*R^4 # = 0

### Eq 3:
### Sum(x3^2*...) =>
((x1*x2*x3)^2+(x2*x3*x4)^2+(x3*x4*x1)^2+(x4*x1*x2)^2) + b*(S^3 - 3*E2*S + 3*E3) - R*(S^2 - 2*E2) # = 0
(E3^2 - 2*E2*E4) + b*(S^3 - 3*E2*S + 3*E3) - R*(S^2 - 2*E2) # = 0
### Eq 3:
E3^2 - 2*E2*E4 + 3*b*E3 - 3*b*E2*S + 2*R*E2 + b*S^3 - R*S^2 # = 0


### Eqs:
### Eq 1:
E4^4 - b^4*E4 + b^3*R*E3 - b^2*R^2*E2 + b*R^3*S - R^4 # = 0
### Eq 2:
4*b*R*E4^2*S - 8*R^2*E4^2 - 4*b^3*R*E3 - b^4*E2^2 + 4*b^2*R^2*E2 + 2*b^4*E3*S +
	- 4*b*R^3*S - b^5*S + b^2*R^2*S^2 - 4*b*R^3*S + 4*b^4*R + 8*R^4 # = 0
### Eq 3:
E3^2 - 2*E2*E4 + 3*b*E3 - 3*b*E2*S + 2*R*E2 + b*S^3 - R*S^2 # = 0


### Test:
x1^2*x2^2 + b*x3 # - R
x2^2*x3^2 + b*x4 # - R
x3^2*x4^2 + b*x1 # - R
x4^2*x1^2 + b*x2 # - R


R = -2
b = 3
# Special Sub-Case type:
x1 = 1.3715291620 + 0.2900893695i;
x2 = 0.4884236632 + 1.7041597366i;
x3 = 1.3715291620 - 0.2900893695i;
x4 = 0.4884236632 - 1.7041597366i;
x = c(x1, x2, x3, x4);
E = debug.E(x)
S = E$S; E2 = E$E2; E3 = E$E3; E4 = E$E4;


