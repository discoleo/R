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
### draft v.0.2d


###############
### History ###
###############


### draft v.0.2d:
# - cleanup of S3Ht L22;
### draft v.0.2c:
# - solved S3L33M:
#   x^3*y^3 + b*x*y = R;
### draft v.0.2b - v.0.2b-clPoly-full:
# - [started work] x^4*y^4 + b*z = R;
# - classic Poly:
#   P[28] for Case x == y or y == z;
#   P[42] for Case (x,y,z) all distinct; [v.0.2b-clPoly-full]
# - Eq S: for cases R = 1, b = +/- 1; [v.0.2b-S-cases]
# - Eq S: full: any parameter b & R; [v.0.2b-Eq-S-b & v.0.2b-Eq-S-full]
### draft v.0.2a:
# - solved: x^3*y^3 + b*z = R;
### draft v.0.1a:
# - cleanup;
# - moved derivation from file:
#   Poly.System.Hetero.Symmetric.S3.Leading.R
#   to this file;

####################
####################

### Helper Functions


source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


# library(polynom)
# library(pracma)

# the functions are in the file:
# Polynomials.Helper.R
# - e.g. round0(), round0.p(),
#   solve.EnAll(), solveEn();

### other

if(FALSE) {
	library(gmp)
	source("Polynomials.Helper.BigNumbers.R")
}

### Dynamic Computation of En:
# - performs elimination of the variable from 2 polynomials;
# - cc, dd = coefficients of the 2 polynomials (univariate);
# - Result: Val + X = 0 => X = - Val of fraction;
ExDynamic = function(cc, dd, debug=TRUE) {
	cc = round0(cc); dd = round0(dd);
	len = length(cc);
	strip = function(x) {
		pos = 0;
		for(i in seq(length(x), 1, by=-1)) {
			if(x[i] == 0) { pos = i; }
			else break;
		}
		if(pos > 0) x = x[seq(pos - 1)];
		x = x / x[length(x)];
		return(x);
	}
	for(i in seq(len - 1)) {
		# Debug
		if(debug) print(paste0("Step: ", i));
		cc = strip(cc);
		dd = strip(dd);
		len1 = length(cc); len2 = length(dd);
		if(len2 <= 2) return(dd);
		while(len1 >= len2) {
			cc = cc - c(rep(0, len1 - len2), dd);
			cc = cc[ - len1];
			cc = strip(cc); len1 = length(cc);
		}
		if(len1 <= 2) return(cc);
		while(len2 >= len1) {
			dd = dd - c(rep(0, len2 - len1), cc);
			dd = dd[ - len2];
			dd = strip(dd); len2 = length(dd);
		}
		if(len2 <= 2) return(dd);
	}
	warning("Some error!");
	return(NA);
}

# TODO: use solve.lpm();
solve.3pm = function(p, v, bigz=FALSE, xn=NULL, val=1, stop.at=NULL) {
	if( ! is.null(xn)) {
		if(length(xn) > length(val)) val = rep(val, length(xn));
		for(i in seq_along(xn)) {
			p[[1]] = replace.pm(p[[1]], val[i], x=xn[i]);
			p[[2]] = replace.pm(p[[2]], val[i], x=xn[i]);
			p[[3]] = replace.pm(p[[3]], val[i], x=xn[i]);
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
classic.P3Lnn = function(n, type="x") {
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
	# Div:
	pR = div.pm(p$Rez, p0, "x")$Rez;
	pR = sort.pm(pR, c(4,3), xn="x");
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
pz2 = toPoly.pm("3*E3^2 + 3*b*E3 - 3*b*E2*S + 2*R*E2 + b*S^3 - R*S^2");
p3 = toPoly.pm("E2*E3 + b*(S^2 - 2*E2) - R*S");
pTriv = data.frame(S=c(4, 1, 0), b=c(0,1,0), R=c(0,0,1), coeff=c(1,27,-81))
# 7*9*b^3*S^3 - 139*b^2*R*S^2 + 97*b*R^2*S - 21*R^3
pDiv = data.frame(S=3:0, b=3:0, R=0:3, coeff=c(7*9,-139,97,-21))
#
p2r = replace.fr.pm(pz2, pE3x0, pE3Div, "E3", 1)
p3r = replace.fr.pm(p3, pE3x0, pE3Div, "E3", 1)
#
pE3r = solve.pm(p2r, p3r, "E2");
pS = pE3r[[1]]
str(pS);
# pS$coeff = pS$coeff / 384;
pS$S = pS$S - min(pS$S); # was S^6!
pS = div.pm(pS, pTriv, "S")$Rez;
pS = div.pm(pS, pDiv, "S")$Rez;
pS = sort.pm(pS, xn="S");
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


### [old] Solver:
# TODO: clean;
# - robust version is in file:
#   Poly.System.Hetero.Symmetric.S3.Leading.R;
solve.S3Ht.L22.Simple.old = function(R, b, be=0, debug=TRUE) {
	coeff = coeff.S3L22.Simple(R, b, be=be);
	S = roots(coeff);
	if(debug) print(S);
	R1 = R - be[1]*S; # Extension
	E2x0 = 2*b*R1*S^5 - 2*R1^2*S^4 - 9*b^3*S^3 + 27*b^2*R1*S^2 + 36*b*R1^2*S - 243*b^4 - 18*R1^3;
	E2Div = 30*b^2*S^4 - 48*b*R1*S^3 + 20*R1^2*S^2 - 189*b^3*S + 81*b^2*R1;
	E2 = E2x0 / E2Div;
	# E2*E3 + b*(S^2 - 2*E2) - R*S
	E3 = - (b*S^2 - 2*b*E2 - R1*S) / E2;
	#
	len = length(S);
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	S = rep(S, each=3); E3 = rep(E3, each=3);
	yz.s = S - x;
	yz = E3 / x;
	# yz.d = sqrt(yz.s^2 - 4*yz)
	# robust
	isEq = round0(x^2*yz.s - b, tol=1E-3) != 0;
	y = z = as.vector(yz.s[isEq]/2);
	sol = cbind(x=as.vector(x[isEq]), y=y, z=z);
	# TODO: each sol2 is duplicated
	# y = x[ ! isEq]; z = yz.s[ ! isEq] - y;
	# sol2 = cbind(x=as.vector(x[ ! isEq]), y=as.vector(y), z=as.vector(z))
	sol2 = sol[, c(2,1,3)];
	sol = rbind(sol, sol2, sol2[, c(1,3,2)]);
	### Case: S == 0
	# does NOT seem to be a valid solution!
	# S = 0; E3 = 2*b; E2 = - (27*b^4 + 2*R^3) / (9*b^2*R);
	# x = roots(c(1, 0, E2, -E3));
	# sol = rbind(sol, x[c(1,2,3)], x[c(1,3,2)]);
	return(sol);
}
coeff.S3Ht.L22.Simple = function(R, b, be=0) {
	coeff = c(b^2, - 2*b*R, R^2, - 9*b^3, 9*b^2*R, - 3*b*R^2, 27*b^4 - R^3);
	if(any(be != 0)) {
		coeff = coeff +
			c(2*b*be[1] + be[1]^2, -2*R*be[1], 0, be[1]^3 - 3*b*be[1]^2 - 9*b^2*be[1],
				-3*R*be[1]^2 + 6*b*R*be[1], 3*R^2*be[1], 0);
	}
	return(coeff);
}

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
E2x0 = (4*b^2*S^7 - 7*b*R*S^6 + 3*R^2*S^5 + 18*b^3*S^2 - 27*b^2*R*S);
E2Div = (13*b^2*S^5 - 21*b*R*S^4 + 9*R^2*S^3 + 27*b^3);
E2 = E2x0 / E2Div;
# TODO: S^6
E2x0 = (5*b^2*R*S^6 - 9*b*R^2*S^5 + 4*R^3*S^4 + 14*b^4*S^2 - 27*b^3*R*S)
E2Div = (18*b^2*R*S^4 - 30*b*R^2*S^3 + 13*R^3*S^2 + 14*b^4);
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
# Reduction (Diff-based) =>
4*E2^2*E3*S + 4*E2*E3*S^3 - 2*E2^3*S^2 - 3*E2^4 + b*S - 3*R # = 0


### Sum(z*...) =>
E3*(x^3*y^3 + y^3*z^3 + z^3*x^3) + b*(S^2 - 2*E2) - R*S # = 0
E3*(3*E3^2 - 3*E3*E2*S + E2^3) + b*(S^2 - 2*E2) - R*S # = 0
3*E3^3 - 3*E3^2*E2*S + E2^3*E3 - 2*b*E2 + b*S^2 - R*S # = 0
# Reduction =>
E2^3*E3 - E2^3*S^3 - E2^2*E3*S^2 + 2*E2*E3*S^4 + b*E2 - R*S # = 0
E3^2*S^3 - E2^2*E3*S^2 + E2^3*E3 + b*E2 - R*S # = 0
E2^3*E3 - E2^3*S^3 - E2^2*E3*S^2 + 2*E2*E3*S^4 + b*E2 - R*S # = 0
4*E2^3*E3 - 3*E2^4*S - 6*E2^3*S^3 + 12*E2*E3*S^4 + 4*b*E2 + b*S^2 - 7*R*S # = 0


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


### Alternatives:
# Diff(z^3*...) =>
E3^3*y - b*(x^3 + z^3 + x*z*(x+z)) + R*(x^2 + z^2 + x*z) # = 0
# Sum(...) =>
E3^3*S - b*(2*(S^3 - 3*E2*S + 3*E3) + E21) + R*(2*S^2 - 3*E2) # = 0
E3^3*S - 2*b*S^3 + 5*b*E2*S - 3*b*E3 + 2*R*S^2 - 3*R*E2 # = 0
# Reduction =>
2*E2*E3^2*S^2 - E2^3*E3*S - 2*b*S^3 + 5*b*E2*S - 3*b*E3 + 2*R*S^2 - 3*R*E2 # = 0
E2^3*E3*S + 19*b*E2*S - 9*b*E3 - 9*R*E2 - 8*b*S^3 + 8*R*S^2 # = 0
E2^3*S^4 + E2^2*E3*S^3- 2*E2*E3*S^5 + 18*b*E2*S - 9*b*E3 - 9*R*E2 - 8*b*S^3 + 9*R*S^2 # = 0


### Eq S:
# P[14] (true roots):
-R^7 - b^8 + 9*b*R^6*S - 39*b^2*R^5*S^2 + 103*b^3*R^4*S^3 - 175*b^4*R^3*S^4 + 187*b^5*R^2*S^5 +
	- 113*b^6*R*S^6 + 29*b^7*S^7 + R^6*S^8 - 6*b*R^5*S^9 + 15*b^2*R^4*S^10 - 20*b^3*R^3*S^11 +
	+ 15*b^4*R^2*S^12 - 6*b^5*R*S^13 + b^6*S^14


### Derivation:
### Case: R = 1; b = 1;
-2 + 9*S - 39*S^2 + 103*S^3 - 175*S^4 + 187*S^5 - 113*S^6 + 29*S^7 + S^8 - 6*S^9 + 15*S^10 +
	- 20*S^11 + 15*S^12 - 6*S^13 + S^14
### Case: R = 1; b = -1;
-2 - 9*S - 39*S^2 - 103*S^3 - 175*S^4 - 187*S^5 - 113*S^6 - 29*S^7 + S^8 + 6*S^9 + 15*S^10 +  
	+ 20*S^11 + 15*S^12 + 6*S^13 + S^14

# R = 1; b = 1; P[30] = P[14] * P[16]
5*S^30 - 30*S^29 + 75*S^28 - 100*S^27 + 75*S^26 - 30*S^25 + 5*S^24 - 1961*S^23 + 18389*S^22 - 68563*S^21 +
	+ 136015*S^20 - 157435*S^19 + 107211*S^18 - 39969*S^17 - 56953*S^16 + 447444*S^15 - 1238976*S^14 +
	+ 1908684*S^13 - 1913058*S^12 + 1336500*S^11 - 641520*S^10 + 128871*S^9 + 595350*S^8 - 2462562*S^7 +
	+ 5060718*S^6 - 6202332*S^5 + 4881384*S^4 - 2558790*S^3 + 890109*S^2 - 203391*S + 39366

### Test
x^4*y^4 + b*z # - R
y^4*z^4 + b*x # - R
z^4*x^4 + b*y # - R

### Debug
R = -2; b = 3;
x =  0.566502800447 + 0.839944522063i;
y = -0.753061286566 + 0.876745266246i;
z = -1.226727559580 - 0.281064956040i;
S = x+y+z; E2 = x*y+x*z+y*z; E3 = x*y*z;
E21 = x*y*(x+y) + x*z*(x+z) + y*z*(y+z);

R = 1; b = 2;
x =  0.3802393865 - 1.0711999097i;
y =  0.5185862609 + 1.0269618578i;
z = -0.7789741406 + 0.7088699732i;
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
p1 = toPoly.pm("4*E2*E3^2 + 2*E3^2*S^2 - 4*E2^2*E3*S + E2^4 + b*S - 3*R");
p2 = toPoly.pm("3*E3^3 - 3*E2*E3^2*S + E2^3*E3 - 2*b*E2 + b*S^2 - R*S");
p3 = toPoly.pm("E3^2 - 2*E2*E3*S + E2^3");
# Alternatives:
p1 = toPoly.pm("4*E2^2*E3*S + 4*E2*E3*S^3 - 2*E2^3*S^2 - 3*E2^4 + b*S - 3*R");
p2 = toPoly.pm("E2^3*E3 - E2^3*S^3 - E2^2*E3*S^2 + 2*E2*E3*S^4 + b*E2 - R*S");
p4 = toPoly.pm("E2^3*S^4 + E2^2*E3*S^3- 2*E2*E3*S^5 + 18*b*E2*S - 9*b*E3 - 9*R*E2 - 8*b*S^3 + 9*R*S^2");

p1$coeff = as.bigz(p1$coeff)
p2$coeff = as.bigz(p2$coeff)
p4$coeff = as.bigz(p4$coeff)
# relatively fast (seconds), but still ugly!
pR = solve.lpm(p1, p2, p4, xn=c("E3", "E2"), stop.at=1, asBigNum=TRUE)
pR = pR[[2]]
str(pR)
table(pR[[2]]$E2)
# Fraction: 264 monomials / 274 monomials with ugly coefficients;

# [old]
p11 = diff.pm(p1, mult.pm(p3, data.frame(E2=c(1,0), S=c(0,2), coeff=c(4,2))))
p21 = diff.pm(p2, mult.pm(p3, data.frame(E3=c(1,0), E2=c(0,1), S=c(0,1), coeff=c(3,3))))
#

solve.S3L44 = function(pS, R=1, b=1, debug=TRUE) {
	if(missing(pS)) {
		coeff = c(b^6, -6*b^5*R, 15*b^4*R^2, -20*b^3*R^3, 15*b^2*R^4, -6*b*R^5, R^6,
				29*b^7, -113*b^6*R, 187*b^5*R^2, -175*b^4*R^3, 103*b^3*R^4, -39*b^2*R^5,
				9*b*R^6, - R^7 - b^8);
		S = roots(coeff);
	} else if(is.data.frame(pS)) {
		S = roots(rev(pS$coeff));
	} else {
		S = roots(pS);
	}
	if(debug) print(length(S));
	print("Starting E2");
	E2 = E2.S3L44(S, R=R, b=b);
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
		yz = yz[id]; yz.s = yz.s[id];
		# coeff = c(x[id]^4*yz, b, R - 2*b*yz.s, b*yz.s^2 - R*yz.s);
		# robust
		b2 = 3*x[id]^4*yz*yz.s;
		coeff = c(b2, -b2*yz.s, b2*yz.s^2 / 3 + b*(yz.s^2 - 2*yz) - R*yz.s);
		roots(coeff);
	});
	x = as.vector(x); x = rep(x, each=2);
	yz.s = as.vector(yz.s); yz.s = rep(yz.s, each=2);
	y = as.vector(y); z = yz.s - y;
	sol = cbind(x=x, y=y, z=z);
	return(sol);
}
E2.S3L44 = function(S, R=1, b=1, type=E2.type) {
	len = nrow(pE2x0);
	E2 = sapply(S, function(S) {
		if(type == 199) {
			Spow = S^seq(len);
			sum(pE2x0$coeff * Spow) / sum(pE2div$coeff * c(1, head(Spow, -2)));
		} else {
			# Spow = S^seq(0, len-1);
			Spow = S^seq(len-1, 0);
			sum(pE2x0$coeff * Spow) / sum(pE2div$coeff * Spow);
		}
		}); print(E2);
	return(E2);
}
### [OLD]
init.E2.S3L44 = function(b=1, type=13, toDouble=TRUE) {
	strS = if(type == 13) { strS = "S13"; } else strS = c("S199", "S197");
	fn = c("E2x0.", "E2div.");
	if(b == 1) {
		E2.files = paste0("S3L44.", fn, strS, ".csv");
	} else if(b == -1) {
		E2.files = paste0("S3L44.", fn, strS, ".b-1.csv");
	} else if(b == -2) {
		E2.files = paste0("S3L44.", fn, strS, ".b-2.csv");
	} else if(b == -3) {
		E2.files = paste0("S3L44.", fn, strS, ".b-3.csv");
	}
	pE2x0 = read.csv(E2.files[1], colClasses=c("numeric", "character"))
	pE2x0$coeff = as.bigz(pE2x0$coeff)
	pE2div = read.csv(E2.files[2], colClasses=c("numeric", "character"))
	pE2div$coeff = as.bigz(pE2div$coeff)
	# *2 only for b = -2;
	if(b == -2 && type == 13) pE2div$coeff = 2 * pE2div$coeff;
	# E2 factors: c(1, 1)
	if(toDouble) r = toDouble.lpm(list(pE2x0, pE2div))
	else r = list(pE2x0 = pE2x0, pE2div = pE2div);
	return(r)
}


library(gmp);

# values are NOT fixed anymore;
# but still problems with E2!
R = 1;
b = -3;
#
E2.type = 199; # 13; # 199;
r = init.E2.S3L44(b=b, type=E2.type)
pE2x0 = r[[1]]; pE2div = r[[2]];
#
solAll = solve.S3L44(b=b);
# solAll = solve.S3L44(p1, b=b);
x = solAll[,1]; y = solAll[,2]; z = solAll[,3];
solAll = cbind(solAll, x+y+z);
# true roots: only 14*6 pairs;
sol = solAll[abs(round0(x^4*y^4 + b*z - R)) < 1E-3,]
x = sol[,1]; y = sol[,2]; z = sol[,3];
nrow(sol);


### S: not needed anymore
pS = read.csv("S3L44.S30.csv")

# pS = read.csv("S3L44.S30.csv", colClasses=c("numeric", "character"))
# pS$coeff = as.bigz(pS$coeff);
### [initial]
# pS = read.csv("S3L44.S.csv", colClasses=c("numeric", "character"))
# pS$coeff = as.bigz(pS$coeff);
# pS = toDouble.pm(pS, scale=1E+50);


# for R = 1; b = 1;
# for R = 1; b = -1;
pS = read.csv("S3L44.S474.b-1.csv", colClasses=c("numeric", "character"))
pS$coeff = as.bigz(pS$coeff);
# pS = pR$Rez;
pS = factorize.p(pS, xn="S")
# write.csv(pS[[1]]$p1, file="S3L44.S30.b-3.csv", row.names=FALSE)
# write.csv(pS[[1]]$p1, file="S3L44.S30.b-1.csv", row.names=FALSE)
# write.csv(pS[[1]]$p1, file="S3L44.S30.csv", row.names=FALSE)

# pGCD = read.csv("_R.Temp.GCD.1.csv", colClasses=c("numeric", "character"))
# pGCD$coeff = as.bigz(pGCD$coeff);
# pAll = read.csv("_R.Temp.ALL.1.csv", colClasses=c("numeric", "character"))
# pAll$coeff = as.bigz(pAll$coeff);


### Reduce E2:
b = -2;
pS = data.frame(S=14:0, coeff = c(b^6, -6*b^5, 15*b^4, -20*b^3, 15*b^2, -6*b, 1,
	29*b^7, -113*b^6, 187*b^5, -175*b^4, 103*b^3, -39*b^2, 9*b, -1 - b^8));
pS$coeff = as.bigz(pS$coeff);
E2.files = c("S3L44.E2x0.S199.b-2.csv", "S3L44.E2div.S197.b-2.csv")
pE2x0 = read.csv(E2.files[1], colClasses=c("numeric", "character"))
pE2x0$coeff = as.bigz(pE2x0$coeff)
pE2div = read.csv(E2.files[2], colClasses=c("numeric", "character"))
pE2div$coeff = as.bigz(pE2div$coeff)
#
pE2x0Red = divByZero.pm(pE2x0, pS, "S")
pE2divRed = divByZero.pm(pE2div, pS, "S")
pE2x0Red$p = pE2x0Red$p[order(-pE2x0Red$p$S), ]
pE2divRed$p = pE2divRed$p[order(-pE2divRed$p$S), ]
# write.csv(pE2x0Red$p, file="S3L44.E2x0.S13.b-2.csv", row.names=FALSE)
# write.csv(pE2divRed$p, file="S3L44.E2div.S13.b-2.csv", row.names=FALSE)
# 46768052394588893382517914646921056628989841375232
# 23384026197294446691258957323460528314494920687616


### Variable elimination:
xn = c("R", "b")
# slow
# pR = solve.3pm(list(p2, p3, p1), c("E2", "E3"), bigz=TRUE, xn=xn)


xn = c("R", "b")
# the actual method used [~1 hour]
# pR = solve.3pm(list(p11, p21, p3), c("E3", "E2"), bigz=TRUE, xn=xn)
# b = -3;
pR = solve.3pm(list(p11, p21, p3), c("E3", "E2"), bigz=TRUE, xn=xn, val=c(1,-3))
# write.csv(pR$Rez, file="S3L44.S474.b-3.csv", row.names=FALSE)
# write.csv(pR$x0, file="S3L44.E2x0.S199.b-3.csv", row.names=FALSE)
# write.csv(pR$div, file="S3L44.E2div.S197.b-3.csv", row.names=FALSE)


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
b^12 - 3*b^8*R^3*x^4 + b^9*R^2*x^5 + b^10*R*x^6 + b^11*x^7 + 3*b^4*R^6*x^8 - 2*b^5*R^5*x^9 - b^6*R^4*x^10 +
	- (4*b^8*R^2 + R^9)*x^12 + b*R*(4*b^8 + R^7)*x^13 - b^4*R^5*x^16 + 5*b^5*R^4*x^17 - 2*b^6*R^3*x^18 +
	- 6*b^7*R^2*x^19 + (4*b^8*R+3*R^8)*x^20 - 8*b*R^7*x^21 + 6*b^2*R^6*x^22 - 5*b^4*R^4*x^24 +
	+ 12*b^5*R^3*x^25 - 12*b^6*R^2*x^26 + 4*b^7*R*x^27 - 3*R^7*x^28 + 13*b*R^6*x^29 - 21*b^2*R^5*x^30 +
	+ 15*b^3*R^4*x^31 - 5*b^4*R^3*x^32 + 3*b^5*R^2*x^33 - 3*b^6*R*x^34 + b^7*x^35 + R^6*x^36 - 6*b*R^5*x^37 +
	+ 15*b^2*R^4*x^38 - 20*b^3*R^3*x^39 + 15*b^4*R^2*x^40 - 6*b^5*R*x^41 + b^6*x^42



### Derivation:

### Case: (x,y,z) all distinct
round(poly.calc(x[seq(1, 84, by=2)]) * b^6, 3)

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


# robust solution:
x^4*(y^3 + z^3 + y*z*(y+z)) - b;
x^4*(y*z*(y^2+z^2) + (y*z)^2) - b*(y+z) + R
x^4*(y*z)^2*(y+z) - b*(y+z)^2 + R*(y+z) + b*y*z
x^8*(y*z)^3*(y+z) + R*x^4*(y*z)*(y+z) - b^2*(y+z) + b*R
(y+z)*(x^8*(y*z)^3 + R*x^4*(y*z) - b^2) + b*R
# =>
- b*R*x^4*(y*z)^2*(x^8*(y*z)^3 + R*x^4*(y*z) - b^2) +
	- b^3*R^2 - b*R^2*(x^8*(y*z)^3 + R*x^4*(y*z) - b^2) +
	+ b*y*z*(x^8*(y*z)^3 + R*x^4*(y*z) - b^2)^2
- b*R*x^4*(R*x^4*(y*z)^3 - b^2*(y*z)^2 + x^8*(R-b*x)*(y*z)) +
	- b*R^2*(x^8*(y*z)^3 + R*x^4*(y*z) - b^2) +
	+ b*y*z*(x^8*(y*z)^3 + R*x^4*(y*z) - b^2)^2 - b^3*R^2
- 2*b*R^2*x^8*(y*z)^3 + b^3*R*x^4*(y*z)^2 - (b*R*x^12*(R-b*x) + b*R^3*x^4)*(y*z) +
	+ b*y*z*(x^8*(y*z)^3 + R*x^4*(y*z) - b^2)^2
(b*x^16*(R-b*x) - b*R^2*x^8)*(y*z)^3 - b^3*R*x^4*(y*z)^2 + (b*R*x^12*(R-b*x) - b*R^3*x^4 + b^5)*(y*z) +
	- 2*b^3*x^8*(R - b*x)
# (y*z)^4 = R - b*x =>
pyz4 = data.frame(yz=c(4,0,0), x=c(0,1,0), b=c(0,1,0), R=c(0,0,1), coeff=c(1,1,-1))
pyz3 = data.frame(
	yz = c( 3, 3,   2,   1, 1, 1, 0),
	x  = c(16, 8,   4,  12, 4, 0, 8),
	b  = c( 1, 1,   3,   1, 1, 5, 3),
	R  = c( 0, 2,   1,   1, 3, 0, 0),
	Rbx= c( 1, 0,   0,   1, 0, 0, 1),
	coeff = c(1,-1, -1,  1,-1, 1,-2)
)
solve.pm(pyz4, pyz3, "yz")


### Derivation of E2:
# - functional;
p1 = toPoly.pm("E3^2 - 2*E2*E3*S + E2^3")
p2 = toPoly.pm("3*E2*E3^2 + 2*E3^2*S^2 - 2*E2^2*E3*S + b*S - 3*R")
p3 = toPoly.pm("3*E2*E3^2*S - 2*E2^3*E3 - 2*b*E2 + b*S^2 - R*S")

pR = solve.pm(p1, p2, "E3")
p2 = replace.fr.pm(p2, pR$x0, pR$div, "E3")
p3 = replace.fr.pm(p3, pR$x0, pR$div, "E3")

toCoeff(p2, "E2")
toCoeff(p3, "E2")

c8 = 30*S^2 / 27;
c7 = - 4*S^4 / 27;
c6 = - 8*S^6 / 27;
c5 = (54*R - 18*S*b) / 27;
c4 = 0;
c3 = (- 96*R*S^4 + 32*S^5*b) / 27;
c2 = (- 48*R*S^6 + 16*S^7*b) / 27;
c1 = (27*R^2 - 18*R*S*b + 3*S^2*b^2) / 27;
c0 = (18*R^2*S^2 - 12*R*S^3*b + 2*S^4*b^2) / 27;

# [excluded] Case: S = 0
d7 = - 4*S^2 / 3;
d6 = - 4*S^4 / 3;
d5 = 0;
d4 = (30*R - 42*b*S) / 3;
d3 = - (4*R*S^2 + 52*b*S^3) / 3;
d2 = - 32*R*S^4 / 3;
d1 = 16*(S*b - R)*S^6 / 3;
d0 = (27*R^2 - 18*b*R*S + 3*b^2*S^2) / 3;

#
c8 = c8 - d7;
c7 = c7 - d6 - c8*d7;
c6 = c6 - d5 - c8*d6;
c5 = c5 - d4 - c8*d5;
c4 = c4 - d3 - c8*d4;
c3 = c3 - d2 - c8*d3;
c2 = c2 - d1 - c8*d2;
c1 = c1 - d0 - c8*d1;
c0 = c0 - c8*d0;

#
d7 = c7*d7 - c6;
d6 = c7*d6 - c5;
d5 = c7*d5 - c4;
d4 = c7*d4 - c3;
d3 = c7*d3 - c2;
d2 = c7*d2 - c1;
d1 = c7*d1 - c0;
d0 = c7*d0;

p2 = toPoly.pm("c7*E2^7 + c6*E2^6 + c5*E2^5 + c4*E2^4 +
	+ c3*E2^3 + c2*E2^2 + c1*E2 + c0");
p3 = toPoly.pm("d7*E2^7 + d6*E2^6 + d5*E2^5 + d4*E2^4 +
	+ d3*E2^3 + d2*E2^2 + d1*E2 + d0");
p2$coeff = as.bigz(p2$coeff);
p3$coeff = as.bigz(p3$coeff);

# still NO chance;
pR = solve.pm(p2, p3, "E2", stop.at=1)

# alternative:
Ex.S3Ht.L44ChZ = function(S, R, b) {
	cc = c((18*R^2*S^2 - 12*R*S^3*b + 2*S^4*b^2) / 27,
		(27*R^2 - 18*R*S*b + 3*S^2*b^2) / 27,
		(- 48*R + 16*b*S)*S^6 / 27,
		(- 96*R + 32*b*S)*S^4 / 27, 0,
		(54*R - 18*S*b) / 27,
		- 8*S^6 / 27, - 4*S^4 / 27, 30*S^2 / 27, 1);
	# [excluded] Case: S = 0
	dd = c((27*R^2 - 18*b*R*S + 3*b^2*S^2) / 3,
		16*(S*b - R)*S^6 / 3,
		- 32*R*S^4 / 3, - (4*R*S^2 + 52*b*S^3) / 3,
		(30*R - 42*b*S) / 3, 0, - 4*S^4 / 3, - 4*S^2 / 3, 1);
	cc = cc - c(0, dd);
	return(ExDynamic(cc, dd));
}

Ex.S3Ht.L44ChZ(S, R, b)


### Classic Solver:
solver.S3L44.classic = function(R, b) {
	coeff = coeff.S3L44.classic(R, b);
	x = roots(coeff);
	# robust:
	yz = solve.yz.S3L44.classic(x, R, b);
	# x = rep(x, each=4);
	yz.s = - b*R / (x^8*(yz)^3 + R*x^4*(yz) - b^2);
	len = length(x);
	y = sapply(seq(len), function(id) {
		yz = yz[id]; yz.s = yz.s[id];
		# robust
		b2 = 3*x[id]^4*yz*yz.s;
		coeff = c(b2, -b2*yz.s, b2*yz.s^2 / 3 + b*(yz.s^2 - 2*yz) - R*yz.s);
		roots(coeff);
	});
	x = as.vector(x); x = rep(x, each=2);
	yz.s = as.vector(yz.s); yz.s = rep(yz.s, each=2);
	y = as.vector(y); z = yz.s - y;
	sol = cbind(x=x, y=y, z=z);
	return(sol);
}
coeff.S3L44.classic = function(R, b) {
	coeff = c(b^6, - 6*b^5*R, 15*b^4*R^2, - 20*b^3*R^3, 15*b^2*R^4, - 6*b*R^5, R^6, b^7,
		- 3*b^6*R, 3*b^5*R^2, - 5*b^4*R^3, 15*b^3*R^4, - 21*b^2*R^5, 13*b*R^6, - 3*R^7,
		4*b^7*R, - 12*b^6*R^2, 12*b^5*R^3, -5*b^4*R^4, 0, 6*b^2*R^6, - 8*b*R^7,
		(4*b^8*R + 3*R^8), - 6*b^7*R^2, - 2*b^6*R^3, 5*b^5*R^4, -b^4*R^5, 0, 0,
		b*R*(4*b^8 + R^7), - (4*b^8*R^2 + R^9), 0, - b^6*R^4, - 2*b^5*R^5, 3*b^4*R^6, b^11,
		b^10*R, b^9*R^2, - 3*b^8*R^3, 0, 0, 0, b^12);
	return(coeff);
}
solve.yz.S3L44.classic = function(x, R, b) {
	Rbx = R - b*x;
	yz0 = 2*b^3*Rbx^5*x^57 - 2*b^2*R*Rbx^5*x^56 - 10*b^3*R^2*Rbx^4*x^49 + 10*b^2*R^3*Rbx^4*x^48 +
		- 2*b^2*R^2*Rbx^5*x^48 + 20*b^3*R^4*Rbx^3*x^41 - 20*b^2*R^5*Rbx^3*x^40 + 8*b^2*R^4*Rbx^4*x^40 +
		- 2*b^7*R*Rbx^3*x^37 + 2*b^6*R^2*Rbx^3*x^36 - 20*b^3*R^6*Rbx^2*x^33 + 20*b^2*R^7*Rbx^2*x^32 +
		- 12*b^2*R^6*Rbx^3*x^32 + 7*b^7*R^3*Rbx^2*x^29 - 7*b^6*R^4*Rbx^2*x^28 + 4*b^6*R^3*Rbx^3*x^28 +
		+ 10*b^3*R^8*Rbx*x^25 - 10*b^2*R^9*Rbx*x^24 + 8*b^2*R^8*Rbx^2*x^24 - 2*b^10*Rbx^3*x^24 +
		- 8*b^7*R^5*Rbx*x^21 + 8*b^6*R^6*Rbx*x^20 - 8*b^6*R^5*Rbx^2*x^20 - 2*b^3*R^10*x^17 +
		+ 2*b^2*R^11*x^16 - 2*b^2*R^10*Rbx*x^16 + 4*b^10*R^2*Rbx^2*x^16 + 3*b^7*R^7*x^13 +
		- 3*b^6*R^8*x^12 + 4*b^6*R^7*Rbx*x^12 - 2*b^10*R^4*Rbx*x^8;
	yzdiv = - b*R*Rbx^5*x^61 + R^2*Rbx^5*x^60 + 5*b*R^3*Rbx^4*x^53 - 5*R^4*Rbx^4*x^52 - R^3*Rbx^5*x^52 +
		- b^5*Rbx^4*x^49 + b^4*R*Rbx^4*x^48 - 4*b^4*Rbx^5*x^48 - 10*b*R^5*Rbx^3*x^45 + 10*R^6*Rbx^3*x^44 +
		+ 5*R^5*Rbx^4*x^44 + 5*b^5*R^2*Rbx^3*x^41 - 5*b^4*R^3*Rbx^3*x^40 + 13*b^4*R^2*Rbx^4*x^40 +
		+ 10*b*R^7*Rbx^2*x^37 - 10*R^8*Rbx^2*x^36 - 10*R^7*Rbx^3*x^36 - 9*b^5*R^4*Rbx^2*x^33 +
		+ 9*b^4*R^5*Rbx^2*x^32 - 12*b^4*R^4*Rbx^3*x^32 - 5*b*R^9*Rbx*x^29 + 5*R^10*Rbx*x^28 +
		+ 10*R^9*Rbx^2*x^28 + b^8*R*Rbx^3*x^28 + 7*b^5*R^6*Rbx*x^25 - 7*b^4*R^7*Rbx*x^24 +
		- 2*b^4*R^6*Rbx^2*x^24 + b*R^11*x^21 - R^12*x^20 - 5*R^11*Rbx*x^20 + b^8*R^3*Rbx^2*x^20 +
		- 2*b^5*R^8*x^17 + 2*b^4*R^9*x^16 + 8*b^4*R^8*Rbx*x^16 - b^12*Rbx^2*x^16 + R^13*x^12 +
		- 5*b^8*R^5*Rbx*x^12 - 3*b^4*R^10*x^8 + 2*b^12*R^2*Rbx*x^8 + 3*b^8*R^7*x^4 - b^12*R^4;
	return(yz0 / yzdiv);
}

R = 2;
b = -5
sol = solver.S3L44.classic(R, b);
x = sol[,1]; y = sol[,2]; z = sol[,3];
S = x+y+z; E2 = x*y+x*z+y*z; E3 = x*y*z;


#########################
#########################

#########################
### Variant:          ###
### (x*y)^n + b*(x*y) ###
#########################

###############
### Order 3 ###
###############

# x^3*y^3 + b*x*y = R
# y^3*z^3 + b*y*z = R
# z^3*x^3 + b*z*x = R

### Solution:

### Case: (x,y,z) all distinct

### Diff =>
y^3*(x^2 + z^2 + x*z) + b*y # = 0
# for R != 0 => (x,y,z) != 0 =>
y^2*(x^2 + z^2 + x*z) + b # = 0
x^2*(y^2 + z^2 + y*z) + b # = 0
z^2*(x^2 + y^2 + y*z) + b # = 0
### Diff =>
x^2*(y+z) + E3 # = 0
y^2*(x+z) + E3 # = 0
z^2*(x+y) + E3 # = 0
### Diff =>
E2 # = 0

### Sum =>
((x*y)^3 + (x*z)^3 + (y*z)^3) + b*E2 - 3*R # = 0
(3*E3^2 - 3*E3*E2*S + E2^3) + b*E2 - 3*R # = 0
### E2 = 0 =>
E3^2 - R # = 0

### Sum(z^3*...) =>
3*E3^3 + b*E3*(x^2+y^2+z^2) - R*(x^3+y^3+z^3) # = 0
3*E3^3 + b*E3*(S^2 - 2*E2) - R*(S^3 - 3*E2*S + 3*E3) # = 0
3*E3^3 - 2*b*E2*E3 + b*E3*S^2 - 3*R*E3 + 3*R*E2*S - R*S^3 # = 0
### E2 = 0 =>
3*E3^3 + b*E3*S^2 - 3*R*E3 - R*S^3 # = 0
### E3^2 = R =>
b*E3*S^2 - R*S^3 # = 0
b^2*R*S^4 - R^2*S^6 # = 0
R*S^2 - b^2 # = 0


### Eq S:
R*S^2 - b^2 # = 0

### Auxiliary Eqs:
# E2 = 0
# b*E3 = R*S


#########################

### Structural Extension:
### + b2*x*y*z

# x^3*y^3 + b2*x*y*z + b1*x*y = R
# y^3*z^3 + b2*x*y*z + b1*y*z = R
# z^3*x^3 + b2*x*y*z + b1*z*x = R

### Solution:

### Sum(z^3*...) with R => R - b2*E3 =>
3*E3^3 + b1*E3*S^2 - 3*(R - b2*E3)*E3 - (R - b2*E3)*S^3 # = 0
3*E3^3 + b1*E3*S^2 + 3*b2*E3^2 - 3*R*E3 - R*S^3 + b2*E3*S^3 # = 0
### E3^2 = R - b2*E3 =>
-3*b2*E3^2 + 3*b2*E3^2 + b2*E3*S^3 + b1*E3*S^2 - R*S^3 # = 0
b2*E3*S^3 + b1*E3*S^2 - R*S^3 # = 0

### Auxiliary Eqs:
# (b2*S + b1)*E3 = R*S;

### Eq S:
# R => R - b2*E3
(R - b2*E3)*S^2 - b1^2 # = 0
# (b2*S + b1)*E3 = R*S;
(R*(b2*S + b1) - b2*R*S)*S^2 - b1^2*(b2*S + b1) # = 0
b1*R*S^2 - b1^2*b2*S - b1^3 # = 0
R*S^2 - b1*b2*S - b1^2 # = 0

### Debug:
R = -3; b = c(-1, -2);
b1 = b[1]; b2 = b[2];
x = -1.0192703451 + 0.9885857104i;
y =  1.0137384347 + 0.5076021802i;
z = -0.3278013984 - 1.0247833637i;
S = x+y+z; E2 = x*y+x*z+y*z; E3 = x*y*z;


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

