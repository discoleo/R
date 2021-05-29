########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
###  == Derivation ==
###
### draft v.0.1d-P16


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R
# - e.g. round0(), round0.p(),
#   solve.EnAll(), solveEn();


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

##############
### Simple ###
##############

###############
### Order 2 ###
###############

### x[i]^2 + b*x[i+1] = R

### Solution:

### Case: all x[i] different;

### Sum =>
S^2 - 2*E2 + b*S - 4*R # = 0

### Diff =>
(x1+x2)*(x1+x3)*(x1+x4)*(x2+x3)*(x2+x4)*(x3+x4) + b^6 # = 0
(-E3^2 + E3*E2*S - E4*S^2) + b^6 # = 0

### Eq 3:
### TODO!


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
m = perm2(4)
E2 = sum(sapply(seq(nrow(m)), function(id) prod(x[which(m[id,] != 0)])))


### Test
x1^2 + b*x2 # - R
x2^2 + b*x3 # - R
x3^2 + b*x4 # - R
x4^2 + b*x1 # - R

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


###
R = -1
b = 3
xi.f = function(x, R, b, n=2) {
	(R - x^n) / b[1];
}
x = roots(c(1, 0, - 2*R, b^3, R^2 - b^2*R))
x1 = x;
x2 = xi.f(x1, R, b); x3 = xi.f(x2, R, b); x4 = xi.f(x3, R, b);

