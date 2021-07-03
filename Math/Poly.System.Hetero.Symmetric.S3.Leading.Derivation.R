########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Heterogeneous Symmetric
###  Mixed Leading Term
###  == Derivation ==
###
### draft v.0.2b-clPoly


###############
### History ###
###############


### draft v.0.2b - v.0.2b-clPoly:
# - [started work] x^4*y^4 + b*z = R;
# - classic Poly: P[28] for Case x == y or y == z;
### draft v.0.2a:
# - solved: x^3*y^3 + b*z = R;
### draft v.0.1a:
# - cleanup;
# - moved derivation from file:
#   Poly.System.Hetero.Symmetric.S3.Leading.R
#   to this file;

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

library(gmp)

solve.3pm = function(p, v, bigz=FALSE, xn=NULL, val=1, stop.at=NULL) {
	if( ! is.null(xn)) {
		for(i in seq_along(xn)) {
			p[[1]] = replace.pm(p[[1]], val, x=xn[i]);
			p[[2]] = replace.pm(p[[2]], val, x=xn[i]);
			p[[3]] = replace.pm(p[[3]], val, x=xn[i]);
		}
	}
	#
	pE3R = solve.pm(p[[1]], p[[2]], v[[1]])
	p3sub = replace.fr.pm(p[[3]], pE3R$x0, pE3R$div, v[[1]]);
	p3sub = sort.pm(p3sub, c(4,3), v[[2]]);
	if(bigz) {
		pE3R$Rez$coeff = as.bigz(pE3R$Rez$coeff);
		p3sub$coeff = as.bigz(p3sub$coeff);
	}
	pS = solve.pm(pE3R$Rez, p3sub, v[[2]], asBigNum=bigz);
	return(pS);
}

### Classic Polynomial:
### Type: x^n*y^n
classic.P3Lnn = function(n, type="x", div=TRUE) {
	id = match(type, c("x", "z"));
	if(is.na(id)) stop("Invalid variable type!")
	p0 = data.frame(x=c(2*n,1,0), b=c(0,1,0), R=c(0,0,1), coeff=c(1,1,-1))
	p1 = data.frame(x=c(2*n,0,0), z=c(0,1,0), b=c(0,1,0), R=c(0,0,1), coeff=c(1,1,-1))
	p2 = data.frame(x=c(n,1,0),   z=c(n,0,0), b=c(0,1,0), R=c(0,0,1), coeff=c(1,1,-1))
	var.name = type
	p = solve.pm(p1, p2, x=var.name)
	if(type == "x") {
		var.other = "z"
		names(p$Rez)[names(p$Rez) == var.other] = "x";
	}
	p$Rez = sort.pm(p$Rez, c(4,3), xn="x")
	#
	if(div) {
		pR = div.pm(p$Rez, p0, "x")$Rez;
		pR = sort.pm(pR, c(4,3), xn="x");
	}
	return(list(pL=p$Rez, p=pR))
}


########################
########################

########################
### Leading: (X*Y)^n ###
########################

###############
### Order 1 ###
###############

### x*y + b*y = R

# x*y + b1*y = R
# y*z + b1*z = R
# z*x + b1*x = R

### Solution:

### Sum =>
E2 + b1*S - 3*R # = 0

### Sum(z*...) =>
3*E3 + b1*E2 - R*S # = 0
# b1*E2 = R*S - 3*E3

### Diff =>
# y*(x-z) = -b1*(y-z)
# z*(x-y) = -b1*(x-z)
# x*(y-z) =  b1*(x-y)
### Prod =>
# E3 = b1^3

### Eq:
b1*E2 + b1^2*S - 3*b1*R # = 0
R*S - 3*E3 + b1^2*S - 3*b1*R # = 0
(R + b1^2)*S - 3*b1*(R + b1^2) # = 0
### Eq:
(R + b1^2)*(S - 3*b1) # = 0
# Note:
# S = 3*b1 is a FALSE solution;


########################
########################

########################
### Mixed-Order: 2+2 ###
########################

### x[i]^2*x[j]^2 + b*x[k]

# x^2*y^2 + b*z = R
# y^2*z^2 + b*x = R
# z^2*x^2 + b*y = R

### Solution:

### Case 1: (x, y, z) distinct;
# - NO solutions;
### Diff =>
(x-z)*(y^2*(x+z) - b) # = 0
y^2*(x+z) - b # = 0
### Sum(...) =>
x^2*(y+z) + y^2*(x+z) + z^2*(x+y) - 3*b # = 0
E2*S - 3*E3 - 3*b # = 0


### Case 2:
# x = y, but z distinct;

### Sum =>
E2^2 - 2*E3*S + b*S - 3*R # = 0

### Sum(z*...) =>
E2*E3 + b*(S^2 - 2*E2) - R*S # = 0

### Sum(z^2*...) =>
3*E3^2 + b*(x^3+y^3+z^3) - R*(x^2+y^2+z^2) # = 0
3*E3^2 + b*(S^3 - 3*E2*S + 3*E3) - R*(S^2 - 2*E2) # = 0
3*E3^2 + 3*b*E3 - 3*b*E2*S + 2*R*E2 + b*S^3 - R*S^2 # = 0


### Eq:
# S = 0;
b^2*S^6 - 2*b*R*S^5 + R^2*S^4 - 9*b^3*S^3 + 9*b^2*R*S^2 - 3*b*R^2*S + 27*b^4 - R^3


### Derivation:

###
pE3x0 = data.frame(E2=c(2,0,0), S=c(0,1,0), b=c(0,1,0), R=c(0,0,1), coeff=c(1,1,-3))
pE3Div = data.frame(S=1, coeff=2)
pz2 = data.frame(
	E3 = c(2, 1, 0, 0, 0, 0),
	E2 = c(0, 0, 1, 1, 0, 0),
	S  = c(0, 0, 1, 0, 3, 2),
	b  = c(0, 1, 1, 0, 1, 0),
	R  = c(0, 0, 0, 1, 0, 1),
	coeff = c(3, 3, -3, 2, 1, -1)
)
p3 = data.frame(
	E3 = c(1, 0, 0, 0),
	E2 = c(1, 0, 1, 0),
	S  = c(0, 2, 0, 1),
	b  = c(0, 1, 1, 0),
	R  = c(0, 0, 0, 1),
	coeff = c(1, 1, -2, -1)
)
pTriv = data.frame(S=c(4, 1, 0), b=c(0,1,0), R=c(0,0,1), coeff=c(1,27,-81))
# 7*9*b^3*S^3 - 139*b^2*R*S^2 + 97*b*R^2*S - 21*R^3
pDiv = data.frame(S=3:0, b=3:0, R=0:3, coeff=c(7*9,-139,97,-21))
#
p2r = replace.fr.pm(pz2, pE3x0, pE3Div, "E3", 1)
p3r = replace.fr.pm(p3, pE3x0, pE3Div, "E3", 1)
#
pE3r = solve.pm(p2r, p3r, "E2");
pS = pE3r[[1]]
pS$coeff = pS$coeff / 384;
pS$S = pS$S - min(pS$S); # S^6!
pS = div.pm(pS, pTriv, "S")$Rez;
pS = div.pm(pS, pDiv, "S")$Rez;
pS = sort.pm(pS, c(4,2), xn="S");
pS
print.p(pS, "S")
### Aux:
pE3r$x0$coeff = pE3r$x0$coeff / 8;
pE3r$div$coeff = pE3r$div$coeff / 8;
S.min = min(pE3r$x0$S, pE3r$div$S);
pE3r$x0$S = pE3r$x0$S - S.min;
pE3r$div$S = pE3r$div$S - S.min;
pE3r$x0 = diff.pm(pE3r$x0, mult.sc.pm(pS, 9))
print.p(pE3r$x0, "S")
print.p(pE3r$div, "S")
cat(paste(toCoeff(pS, "S"), collapse=",\n"))


# FALSE "solutions"
9*(E3 + b)^2 - 2*E3*S^3 + b*S^3 - 3*R*S^2 # = 0
9*E3^2 + 18*b*E3 - 2*E3*S^3 + b*S^3 - 3*R*S^2 + 9*b^2 # = 0
#
p2f = data.frame(
	E3 = c(2, 1, 1,   0, 0, 0),
	S  = c(0, 0, 3,   3, 2, 0),
	b  = c(0, 1, 0,   1, 0, 2),
	R  = c(0, 0, 0,   0, 1, 0),
	coeff = c(9, 18, -2, 1, -3, 9)
)
pE3x0f = data.frame(E3=c(1,0), b=c(0,1), coeff=c(3,3))
pE3Divf = data.frame(S=1, coeff=3)
#
p3rf = replace.fr.pm(p3, pE3x0f, pE3Divf, "E2", 1)
p3rf$coeff = p3rf$coeff / 3;
p3rf
pE3rf = solve.pm(p3rf, p2f, "E3");
pSf = pE3rf[[1]]
pSf$S = pSf$S - min(pSf$S); # S^2!
pSf = sort.pm(pSf, c(4,2), xn="S")
pSf
print.p(pSf, "S")
print.p(pE3rf$x0, "S")


### Test
x^2*y^2 + b*z # - R
y^2*z^2 + b*x # - R
z^2*x^2 + b*y # - R

### Debug:
R = 2; b = -1;
x = -0.6359778;
y = -1.8364059538;
z = -0.6359778;
S = x+y+z; E2 = x*y+x*z+y*z; E3 = x*y*z;


### Case: all distinct
# but NO solutions;
# S = 0;
4*b*S^7 - 4*R*S^6 - 36*b^2*S^4 + 36*R^2*S^2 + 243*b^3*S - 243*b^2*R

### Case: S == 0
# - still NO solutions;
E2^2 - 3*R # = 0
# =>
# E3 = b;


### Classic Polynomial:
R = -1;
b = 3
x0 = roots(c(1, 0, 0, - b, - R, 0, b^2));
x = x0; y = (R - x0^4) / b; S = 2*x + y;
err = x^6 - b*x^3 - R*x^2 + b^2;
round0(err)

x = y;
err = b^2*x^6 - 2*b*R*x^5 + R^2*x^4 - b^3*x^3 + 3*b^2*R*x^2 - b*R^2*x + b^4 - R^3;
round0(err)
# S = 2*x0 + (R - x0^4)/b;
err = b^2*S^6 - 2*b*R*S^5 + R^2*S^4 - 9*b^3*S^3 + 9*b^2*R*S^2 - 3*b*R^2*S + 27*b^4 - R^3;
round0(err)


########################
########################

########################
### Mixed-Order: 3+3 ###
########################

### x[i]^3*x[j]^3 + b*x[k]

# x^3*y^3 + b*z = R
# y^3*z^3 + b*x = R
# z^3*x^3 + b*y = R

### Solution:

### Sum =>
(x^3*y^3 + y^3*z^3 + z^3*x^3) + b*S - 3*R # = 0
E2^3 - 3*E3*E2*S + 3*E3^2 + b*S - 3*R

### Sum(z*...) =>
E3*(x^2*y^2 + y^2*z^2 + z^2*x^2) + b*(S^2 - 2*E2) - R*S # = 0
E3*(E2^2 - 2*E3*S) + b*(S^2 - 2*E2) - R*S # = 0
2*E3^2*S - E2^2*E3 + 2*b*E2 - b*S^2 + R*S # = 0

### Diff =>
y^3*(x^3-z^3) - b*(x-z) # = 0
### Case: (x, y, z) distinct =>
y^3*(x^2 + z^2 + x*z) - b # = 0
z^3*(x^2 + y^2 + x*y) - b # = 0
x^3*(y^2 + z^2 + y*z) - b # = 0
### Diff =>
x^2*(y^3-z^3) + y^2*z^2*(y-z) + E3*(y^2-z^2) # = 0
x^2*(y^2 + z^2 + y*z) + y^2*z^2 + E3*(y+z) # = 0
### Sum =>
3*(x^2*y^2 + y^2*z^2 + z^2*x^2) + 3*E3*S # = 0
E2^2 - E3*S # = 0
# E3*S = E2^2

### Auxiliary Eqs:
E2x0 = (32*b^2*S^7 - 56*b*R*S^6 + 24*R^2*S^5 + 144*b^3*S^2 - 216*b^2*R*S);
E2Div = (104*b^2*S^5 - 168*b*R*S^4 + 72*R^2*S^3 + 216*b^3);
E2 = E2x0 / E2Div;

### Eq S:
b^3*S^5 - 3*b^2*R*S^4 + 3*b*R^2*S^3 - R^3*S^2 + b^4
### P[5] * (S^6 - 27*b*S + 81*R)
b^3*S^11 - 3*b^2*R*S^10 + 3*b*R^2*S^9 - R^3*S^8 - 26*b^4*S^6 + 162*b^3*R*S^5 - 324*b^2*R^2*S^4 +
	+ 270*b*R^3*S^3 - 81*R^4*S^2 - 27*b^5*S + 81*b^4*R


### Solve:
### Eq 1:
3*E2^4 - 2*E2^3*S^2 + b*S^3 - 3*R*S^2 # = 0
### Eq 2:
E2^4 + 2*b*E2*S - b*S^3 + R*S^2 # = 0

p1 = data.frame(
	E2 = c(4, 3, 0, 0), S = c(0, 2, 3, 2),
	b  = c(0, 0, 1, 0), R = c(0, 0, 0, 1), coeff = c(3,-2,1,-3)
)
p2 = data.frame(
	E2 = c(4, 1, 0, 0), S = c(0, 1, 3, 2),
	b  = c(0, 1, 1, 0), R = c(0, 0, 0, 1), coeff = c(1, 2,-1,1)
)
pSr = solve.pm(p1, p2, xn="E2")
pSr$Rez$coeff = pSr$Rez$coeff / 1152;
pSr$Rez$S = pSr$Rez$S - min(pSr$Rez$S);
pSr$Rez$b = pSr$Rez$b - min(pSr$Rez$b);
pS = sort.pm(pSr$Rez, c(4,3), "S")
print.p(pS, "S")
pS = div.pm(pS, data.frame(S=c(6,1,0), b=c(0,1,0), R=c(0,0,1), coeff=c(1,-3^3, 3^4)), "S")
pS = sort.pm(pS$Rez, c(4,3), "S")
print.p(pS, "S")
#
pE2x0 = pSr$x0; pE2Div = pSr$div;
Scmm = min(pE2x0$S, pE2Div$S);
pE2x0$S = pE2x0$S - Scmm; pE2Div$S = pE2Div$S - Scmm;
pE2x0 = sort.pm(pE2x0, c(4,3), "S");
pE2Div = sort.pm(pE2Div, c(4,3), "S");
pE2x0; pE2Div;
print.p(pE2x0, "S");
print.p(pE2Div, "S");


### Classic Polynomial:

### Case: x == y
# x^6 + b*z = R
# x^3*z^3 + b*x = R
# (but one has to know this)

n = 3
p = classic.P3Lnn(n)
print.p(p$pL, "x")
print.p(p$p, "x")

### Case: x != y != z
E2^2 - E3*S # = 0
#
p3 = data.frame(
	x=c(2,2,0,2,1,1), y=c(2,0,2,1,2,1), z=c(0,2,2,1,1,2),
	coeff = c(1,1,1,1,1,1)
)
p1 = data.frame(x=c(n,0,0), y=c(n,0,0), z=c(0,1,0), b = c(0,1,0), R=c(0,0,1), coeff=c(1,1,-1))
p2 = data.frame(x=c(n,0,0), y=c(0,1,0), z=c(n,0,0), b = c(0,1,0), R=c(0,0,1), coeff=c(1,1,-1))
#
# pR = solve.pm(p1, p2, "y")
# p3 = replace.fr.pm(p3, pR$x0, pR$div, pR$xn);
# solve.pm(pR$Rez, p3, "z")
pR = solve.pm(p3, p2, "y")
p1 = replace.fr.pm(p1, pR$x0, pR$div, pR$xn);
pR = solve.pm(pR$Rez, p1, "z")
pR = sort.pm(pR$Rez, c(4,3), xn="x")
print.p(pR, "x")


# Case x == y;
x^21 - 3*R*x^15 + 3*R^2*x^9 - R^3*x^3 - b^4*x + b^3*R
# Case y == z;
b^3*x^21 - 3*b^2*R*x^20 + 3*b*R^2*x^19 - R^3*x^18 + 3*b^2*R^2*x^14 - 6*b*R^3*x^13 + 3*R^4*x^12 +
	- 2*b^5*x^11 + 4*b^4*R*x^10 - 2*b^3*R^2*x^9 + 3*b*R^4*x^7 - 3*R^5*x^6 + 6*b^4*R^2*x^4 +
	- 6*b^3*R^3*x^3 + b^7*x - b^6*R + R^6
# Case: all distinct
b^3*x^15 - 3*b^2*R*x^14 + 3*b*R^2*x^13 - R^3*x^12 + b^4*x^10 - b^3*R*x^9 +
	- b*R^3*x^7 + R^4*x^6 + b^5*x^5 + b^4*R*x^4 - 2*b^3*R^2*x^3 + b^6
# All Cases:
(x^6 + b*x - R) * (x^15 - b*x^10 - 2*R*x^9 + b^2*x^5 + b*R*x^4 + R^2*x^3 - b^3) *
(b^3*x^15 - 3*b^2*R*x^14 + 3*b*R^2*x^13 - R^3*x^12 - b^4*x^10 + 4*b^3*R*x^9 - 3*b^2*R^2*x^8 +
	- 2*b*R^3*x^7 + 2*R^4*x^6 - b^5*x^5 - b^4*R*x^4 + 5*b^3*R^2*x^3 - b^2*R^3*x^2 - b*R^4*x + b^6 - R^5) *
(b^3*x^15 - 3*b^2*R*x^14 + 3*b*R^2*x^13 - R^3*x^12 + b^4*x^10 - b^3*R*x^9 +
	- b*R^3*x^7 + R^4*x^6 + b^5*x^5 + b^4*R*x^4 - 2*b^3*R^2*x^3 + b^6)


### Test
x^3*y^3 + b*z # - R
y^3*z^3 + b*x # - R
z^3*x^3 + b*y # - R

### Debug
R = -1; b = 3;
x = roots(c(1,0,0,0,0, - b, - 2*R, 0,0,0, b^2, b*R, R^2, 0,0, - b^3))
x = x[1]; y = x; z = (R - x^6) / b;
S = x+y+z; E2 = x*y+x*z+y*z; E3 = x*y*z;


########################
########################

########################
### Mixed-Order: 4+4 ###
########################

### x[i]^4*x[j]^4 + b*x[k]

# x^4*y^4 + b*z = R
# y^4*z^4 + b*x = R
# z^4*x^4 + b*y = R

### Solution:

### Sum =>
(x^4*y^4 + y^4*z^4 + z^4*x^4) + b*S - 3*R # = 0
4*E2*E3^2 + 2*E3^2*S^2 - 4*E2^2*E3*S + E2^4 + b*S - 3*R

### Sum(z*...) =>
E3*(x^3*y^3 + y^3*z^3 + z^3*x^3) + b*(S^2 - 2*E2) - R*S # = 0
E3*(3*E3^2 - 3*E3*E2*S + E2^3) + b*(S^2 - 2*E2) - R*S # = 0
3*E3^3 - 3*E3^2*E2*S + E2^3*E3 - 2*b*E2 + b*S^2 - R*S # = 0

### Diff =>
y^4*(x^4-z^4) - b*(x-z) # = 0
### Case: (x, y, z) distinct =>
y^4*(x^3 + z^3 + x*z*(x+z)) - b # = 0
z^4*(x^3 + y^3 + x*y*(x+y)) - b # = 0
x^4*(y^3 + z^3 + y*z*(y+z)) - b # = 0
### Diff =>
x^3*(y^4-z^4) + y^3*z^3*(y-z) + E3*(x*(y^2+z^2+y*z)*(y-z) + y*z*(y+z)*(y-z)) # = 0
(x^3*y^3 + x^3*z^3 + y^3*z^3) + E3*(x^2*y + x^2*z + x*y^2 + x*z^2 + y^2*z + y*z^2 + E3) # = 0
(3*E3^2 - 3*E3*E2*S + E2^3) + E3*(E2*S - 2*E3) # = 0
E3^2 - 2*E3*E2*S + E2^3 # = 0


### Test
x^4*y^4 + b*z # - R
y^4*z^4 + b*x # - R
z^4*x^4 + b*y # - R

### Debug
R = -2; b = 3;
x =  0.5665027820 + 0.8399445356i;
y = -0.7530613076 + 0.8767452519i;
z = -1.2267275591 - 0.2810649213i;
S = x+y+z; E2 = x*y+x*z+y*z; E3 = x*y*z;


### Eq S Derivation:
4*E2*E3^2 + 2*E3^2*S^2 - 4*E2^2*E3*S + E2^4 + b*S - 3*R # = 0
3*E3^3 - 3*E2*E3^2*S + E2^3*E3 - 2*b*E2 + b*S^2 - R*S # = 0
E3^2 - 2*E2*E3*S + E2^3 # = 0
# =>
2*E3^2*S^2 + 4*E2^2*E3*S - 3*E2^4 + b*S - 3*R # - 2*S^2 =>
3*E2*E3^2*S - 2*E2^3*E3 - 2*b*E2 + b*S^2 - R*S # - 3*E2*S =>
# =>
4*E2^2*E3*S + 4*E2*E3*S^3 - 3*E2^4 - 2*E2^3*S^2 + b*S - 3*R
3*E2^4*S + 2*E2^3*E3 - 6*E2^2*E3*S^2 + 2*b*E2 - b*S^2 + R*S
# =>
(2*E2^3*S^2 + 4*E2^2*S^4 - 3*E2^4 + b*S - 3*R)^2 - (4*E2^2*S + 4*E2*S^3)^2*(E2^2*S^2 - E2^3)
(5*E2^4*S - 6*E2^3*S^3 + 2*b*E2 - b*S^2 + R*S)^2 - (2*E2^3 - 6*E2^2*S^2)^2*(E2^2*S^2 - E2^3)
#
p1 = data.frame(
	E3 = c(2, 2, 1, 0, 0, 0),
	E2 = c(1, 0, 2, 4, 0, 0),
	S  = c(0, 2, 1, 0, 1, 0),
	b  = c(0, 0, 0, 0, 1, 0),
	R  = c(0, 0, 0, 0, 0, 1),
	coeff = c(4, 2, -4, 1, 1, -3)
)
p2 = data.frame(
	E3 = c(3, 2, 1, 0, 0, 0),
	E2 = c(0, 1, 3, 1, 0, 0),
	S  = c(0, 1, 0, 0, 2, 1),
	b  = c(0, 0, 0, 1, 1, 0),
	R  = c(0, 0, 0, 0, 0, 1),
	coeff = c(3,-3, 1,-2, 1,-1)
)
p3 = data.frame(
	E3 = c(2, 1, 0),
	E2 = c(0, 1, 3),
	S  = c(0, 1, 0),
	coeff = c(1,-2,1)
)
p11 = diff.pm(p1, mult.pm(p3, data.frame(E2=c(1,0), S=c(0,2), coeff=c(4,2))))
p21 = diff.pm(p2, mult.pm(p3, data.frame(E3=c(1,0), E2=c(0,1), S=c(0,1), coeff=c(3,3))))
#

solve.P3L44 = function(pS, debug=TRUE) {
	if(is.data.frame(pS)) {
		S = roots(rev(pS$coeff));
	} else {
		S = roots(pS);
	}
	if(debug) print(length(S));
	len = nrow(pE2x0);
	print("Starting E2");
	E2 = sapply(S, function(S) {
			Spow = S^seq(len);
			sum(pE2x0$coeff * Spow) / sum(pE2div$coeff * c(0, head(Spow, -2)))
		});
	R = 1; b = 1;
	E3 = (3*E2^4 + 2*E2^3*S^2 - b*S + 3*R) / (4*E2^2*S + 4*E2*S^3);
	len = length(S);
	print("Starting x");
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	S = rep(S, each=3); E3 = rep(E3, each=3);
	yz = E3 / x; yz.s = S - x;
	len = length(x);
	if(debug) print(len);
	print("Starting y");
	y = sapply(seq(len), function(id) {
		coeff = c(x[id]^4*yz[id], b, R - 2*b*yz.s[id], b*yz.s[id]^2 - R*yz.s[id]);
		roots(coeff);
	});
	x = as.vector(x); x = rep(x, each=3);
	yz.s = as.vector(yz.s); yz.s = rep(yz.s, each=3);
	y = as.vector(y); z = yz.s - y;
	sol = cbind(x=x, y=y, z=z);
	return(sol);
}
R = 1; b = 1;
sol = solve.P3L44(pS);
x = sol[,1]; y = sol[,2]; z = sol[,3];
sol = cbind(sol, x+y+z);
# true roots: only 18 pairs!
sol2 = sol[abs(round0(x^4*y^4 + b*z - R)) < 1E-4,]


### from csv:
library(gmp);
pE2x0 = read.csv("S3L44.E2x0.csv", colClasses=c("numeric", "character"))
pE2x0$coeff = as.bigz(pE2x0$coeff)
pE2div = read.csv("S3L44.E2div.csv", colClasses=c("numeric", "character"))
pE2div$coeff = as.bigz(pE2div$coeff)
r = toDouble.lpm(list(pE2x0, pE2div));
pE2x0 = r[[1]]; pE2div = r[[2]];
### S
pS = read.csv("S3L44.S.csv", colClasses=c("numeric", "character"))
pS$coeff = as.bigz(pS$coeff);
pS = toDouble.pm(pS, scale=1E+50);


### Variable elimination:
xn = c("R", "b")
# slow
# pR = solve.3pm(list(p2, p3, p1), c("E2", "E3"), bigz=TRUE, xn=xn)

# the actual method used [~1 hour]
pR = solve.3pm(list(p11, p21, p3), c("E3", "E2"), bigz=TRUE, xn=xn)

# [failed as well]
pI1 = data.frame(x=c(4,0,0), y=c(4,0,0), z=c(0,1,0), b=c(0,1,0), R=c(0,0,1), coeff=c(1,1,-1))
pI2 = pI3 = pI1; nm = c("x", "y", "z");
names(pI2)[1:3] = nm[c(2,3,1)];
names(pI3)[1:3] = nm[c(3,1,2)];
#
pR = solve.3pm(list(pI1, pI3, pI2), c("z", "y"), bigz=TRUE)
pR2 = div.pm(pR$Rez, classic.P3Lnn(n, type="x")$pL, "x")
pR2 = div.pm(pR2$Rez, classic.P3Lnn(n, type="z")$p, "x")


### Classic Polynomial:

### Case: x == y
# Note: z = (R - x^8) / b; # is the solution for the Case y == z;
x^36 - 4*R*x^28 + 6*R^2*x^20 - 4*R^3*x^12 + R^4*x^4 + b^5*x - b^4*R
# P[8] * P[28]
x^28 - b*x^21 - 3*R*x^20 + b^2*x^14 + 2*b*R*x^13 + 3*R^2*x^12 - b^3*x^7 - b^2*R*x^6 - b*R^2*x^5 - R^3*x^4 + b^4

### Case: y == z
# x = z; # test polynomial
b^4*x^36 - 4*b^3*R*x^35 + 6*b^2*R^2*x^34 - 4*b*R^3*x^33 + R^4*x^32 + 4*b^3*R^2*x^27 - 12*b^2*R^3*x^26 +
	+ 12*b*R^4*x^25 - 4*R^5*x^24 + 6*b^2*R^4*x^18 - 12*b*R^5*x^17 + 6*R^6*x^16 - 8*b^6*R*x^14 +
	+ 16*b^5*R^2*x^13 - 8*b^4*R^3*x^12 + 4*b*R^6*x^9 - 4*R^7*x^8 + 8*b^5*R^3*x^5 - 8*b^4*R^4*x^4 + b^9*x +
	- b^8*R + R^8
# P[8] * P[28]
b^4*x^28 - 4*b^3*R*x^27 + 6*b^2*R^2*x^26 - 4*b*R^3*x^25 + R^4*x^24 - b^5*x^21 + 5*b^4*R*x^20 - 6*b^3*R^2*x^19 +
	- 2*b^2*R^3*x^18 + 7*b*R^4*x^17 - 3*R^5*x^16 + b^6*x^14 - 6*b^5*R*x^13 + 11*b^4*R^2*x^12 +
	- 4*b^3*R^3*x^11 - 3*b^2*R^4*x^10 - 2*b*R^5*x^9 + 3*R^6*x^8 - b^7*x^7 - b^6*R*x^6 - b^5*R^2*x^5 +
	+ 7*b^4*R^3*x^4 - b^3*R^4*x^3 - b^2*R^5*x^2 - b*R^6*x + b^8 - R^7

### Case: (x,y,z) all distinct
### TODO

### Derivation:

### Case: x == y
# x^8 + b*z = R
# x^4*z^4 + b*x = R
# (but one has to know this)

n = 4
p = classic.P3Lnn(n, type="z")
print.p(p$pL, "x")
print.p(p$p, "x")

R = -1
b = 3
coeff = rev(eval(parse(text=paste0("c(", paste(toCoeff(p$p, "x"), collapse=", "), ")"))));
x = roots(coeff);
y = x; z = (R - x^8)/b;


########################
########################
########################

########################
### Mixed-Order: 2+1 ###
########################

### x[i]^2*x[j] + b*x[j]

# x^2*y + b*y = R
# y^2*z + b*z = R
# z^2*x + b*x = R

### Solution:

# - shortcut: using eq. from the HtSymmetric Mixed system;

### Sum =>
# shortcut: Eq1[Mixed](R1 = 3*R - b*S) =>
E3*S^3 - ((3*R - b*S) + 6*E3)*E2*S + (3*R - b*S)^2 + E2^3 + 9*E3^2 + 3*(3*R - b*S)*E3 # = 0

# long:
(x^2*y + y^2*z + z^2*x) + b*S - 3*R # = 0 # Eq 1-bis
(x^2*y + y^2*z + z^2*x)*(x*y^2 + y*z^2 + z*x^2) +
	+ (b*S - 3*R)*(x*y^2 + y*z^2 + z*x^2) # = 0
E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2 +
	+ (b*S - 3*R)*(x*y^2 + y*z^2 + z*x^2) # = 0 # Eq 2-bis
### Sum Eq 1-bis + Eq 2-bis =>
(b*S - 3*R)*(x*y^2 + y*z^2 + z*x^2 + x^2*y + y^2*z + z^2*x) +
	+ (b*S - 3*R)^2 + E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2 # = 0
(b*S - 3*R)*(E2*S - 3*E3) +
	+ b^2*S^2 - 6*b*R*S + 9*R^2 + E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2 # = 0
E3*S^3 - 3*b*E3*S + 9*E3^2 + 9*R*E3 - 6*E3*E2*S +
	+ E2^3 + b*E2*S^2 - 3*R*E2*S + b^2*S^2 - 6*b*R*S + 9*R^2 # = 0

### Sum(z*...) =>
x*y*z*(x+y+z) + b*E2 - R*S # = 0
E3*S + b*E2 - R*S # = 0

### Sum(y*...) =>
(x^2*y^2 + y^2*z^2 + z^2*x^2) + b*(x^2 + y^2 + z^2) - R*S # = 0
E2^2 - 2*E3*S + b*(S^2 - 2*E2) - R*S # = 0
E2^2 - 2*E3*S + b*S^2 - 2*b*E2 - R*S # = 0

### Auxiliary:
E3Subst = - 27*R*S*b^5 - 12*R*S^3*b^4 + 6*R^2*S^4*b^2 + 9*S^2*b^6 + 5*S^4*b^5;
E3Div = - 27*R*S^2*b^3 - 6*R*S^4*b^2 + 3*S^3*b^4 - S^5*b^3;
# E3 = - E3Subst / E3Div;


### Eq S:
((R^2 + b^3)*S^2 - b^2*R*S + b^4) * S^4 * (S^3 + 9*b*S - 27*R) * P[9]

### P[9]: false solution;
(- 6561*R^2*b^3) +
(2916*R*b^4)*S^1 +
(- 6561*R^2*b^2 - 243*b^5)*S^2 +
(2916*R*b^3 - 729*R^3)*S^3 +
(- 972*R^2*b - 189*b^4)*S^4 +
(621*R*b^2)*S^5 +
(- 9*b^3)*S^6 +
(54*R*b)*S^7 + b^2*S^8 + R*S^9


### Debug:
R = 2; b = -1;
x =  1.5907409211 - 0.9236008907i;
y =  0.1489943605 + 0.6462891180i;
z = -1.4064019508 - 0.1940927468i;
S = x+y+z; E2 = x*y+x*z+y*z; E3 = x*y*z;

