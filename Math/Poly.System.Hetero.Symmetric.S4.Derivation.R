########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
###  == Derivation ==
###
### draft v.0.3d-Refactor3


### History

# draft v.0.3d:
# - moved L1 V2 code from this file to:
#   Poly.System.Hetero.Symmetric.S4.Derivation.L2V2.R;
# - moved L1 V3 code from this file to:
#   Poly.System.Hetero.Symmetric.S4.Derivation.L1V3.R;
# - moved L2 V1 code from this file to:
#   Poly.System.Hetero.Symmetric.S4.Derivation.L2V1.R;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")

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
isConj = function(x, y, tol=1E-3) {
	isConj = (abs(Re(x) - Re(y)) < tol) & (abs(Im(x) + Im(y)) < tol);
	return(isConj);
}
### generate Classic Polynomial:
bR.gen = function(pb, pR=1) data.frame(b = pb, R = pR, coeff = 1)
bx.gen = function(pb, px=1) data.frame(b = pb, x = px, coeff = 1)
classic.BaseSimple.gen = function(n) data.frame(x=c(n, 1, 0), b=c(0,1,0), R=c(0,0,1), coeff=c(1,1,-1));
classic.S2Simple.gen = function(n=3) {
	p1 = data.frame(
		x = c(n,0), b = c(0,0), R = c(0,1),
		coeff = c(-1, 1)
	)
	p1 = pow.pm(p1, n);
	p1 = diff.pm(p1, bR.gen(n, pR=1));
	p1 = add.pm(p1, bx.gen(n+1, px=1));
	if(n %% 2 == 1) p1$coeff = - p1$coeff;
	p1 = sort.pm(p1, sort.coeff=c(4,2,3,1), xn="x")
	rownames(p1) = seq(nrow(p1))
	return(p1);
}
classic.S4Simple.gen = function(n = 3) {
	p1 = data.frame(
		x = c(n,0), b = c(0,0), R = c(0,1),
		coeff = c(-1, 1)
	)
	p1 = pow.pm(p1, n)
	p1 = diff.pm(bR.gen(n, pR=1), p1)
	p1 = pow.pm(p1, n)
	p1 = diff.pm(bR.gen(n^2+n, pR=1), p1)
	p1 = pow.pm(p1, n)
	p1 = diff.pm(p1, bR.gen((n^4-1)/(n-1) - 1, pR=1))
	p1 = add.pm(p1, bx.gen((n^4-1)/(n-1), px=1))
	if(n %% 2 == 1) p1$coeff = - p1$coeff;
	p1 = sort.pm(p1, sort.coeff=c(4,2,3,1), xn="x")
	rownames(p1) = seq(nrow(p1))
	return(p1);
}

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


### Eq S:
S^18 - 15*R*S^15 + 48*b^2*S^14 - 126*R*b*S^13 + (222*R^2 - 256*b^3)*S^12 + 609*R*b^2*S^11 +
	- (1764*R^2*b + 540*b^4)*S^10 + (2158*R^3 + 5061*R*b^3)*S^9 +
	- (8433*R^2*b^2 + 960*b^5)*S^8 + (6048*R^3*b - 1170*R*b^4)*S^7 +
	- (7671*R^4 - 3435*R^2*b^3 - 5800*b^6)*S^6 + (6099*R^3*b^2 - 18840*R*b^5)*S^5 +
	- (6300*R^4*b - 16632*R^2*b^4 + 3600*b^7)*S^4  + (8049*R^5 - 23297*R^3*b^3 + 10080*R*b^6)*S^3 +
	+ (1677*R^4*b^2 + 3672*R^2*b^5 - 10125*b^8)*S^2 +
	+ (2142*R^5*b - 7470*R^3*b^4 + 8100*R*b^7)*S +
	- 2744*R^6 + 9225*R^4*b^3 - 8910*R^2*b^6;


### TODO:
# - check thoroughly Eq for S:
#  -- substantial numerical instability for coeffs of S;
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
classic.S4Simple.gen = function(n = 3) {
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
### P[81]
p1 = classic.S4Simple.gen();

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

###
pP1 = toPoly.pm("S^3 - 3*E2*S + 3*E3 + b*S - 4*R")
pP2 = toPoly.pm("3*E4*S^2 + 3*R*b*S - 3*E2*E3*S - 6*R^2 - b^2*E2 + E2^3 + 3*E3^2 - 3*E2*E4")
pP3 = toPoly.pm("E2^2*S^2 + b*E2*S^2 - 2*E2^3 - E2*E3*S + 4*E2*E4 - 3*R*E2*S - b*E3*S + b*E4 + b^3")
pP4 = toPoly.pm("b^4*E4 + R^3*(S^3 - 3*E2*S + 3*E3) +
	- R^2*(E2^3 + 3*E3^2 - 3*E3*E2*S + 3*E4*S^2 - 3*E2*E4) +
	+ R*(E3^3 - 3*E4*E3*E2 + 3*E4^2*S) - E4^3 - R^4") # the Mult (from Derivation) variant;

pR = solve.lpm(pP1, pP3, pP2, pP4, xn=c("E3", "E4", "E2"), stop.at=1, asBigNum = TRUE)
#
pP1 = toPoly.pm("c1 - 3*E2*S + 3*E3")
pP2 = toPoly.pm("c20 + 3*E4*S^2 - 3*E2*E3*S - b^2*E2 + E2^3 + 3*E3^2 - 3*E2*E4")
#
pP2 = toPoly.pm("c2 - 3*S*E2*c1 + 3*E2^3 + 9*S^2*E4 - 9*E2*E4 - 3*E2*b^2")
pP3 = toPoly.pm("c3 - S*E2*c1 + 6*E2^3 - 12*E2*E4 - 3*b*E4 + 9*S*E2*R")
pP4 = toPoly.pm("c4 - 9*S*E2*R*c1^2 + 27*S^2*E2^2*R*c1 - 27*S*E2*R^2*c1 +
	- 27*E2*R*E4*c1 - 27*S^3*E2^3*R + 27*E2^3*R^2 +
	+ 81*S*E2^2*R*E4 + 81*S^2*R^2*E4 - 81*E2*R^2*E4 - 81*S*R*E4^2 + 27*E4^3 - 27*E4*b^4")
pP4 = toPoly.pm("3*E4^3 - 9*R*S*E4^2 + (c412*E2^2 - c411*E2 + c410)*E4 +
	+ 3*c403*E2^3 + 3*c1*R*S^2*E2^2 - c401*E2 + c4")
#
pP2 = toPoly.pm("6*E2^4 - c33*E2^3 + c32*E2^2 + c31*E2 + c5")
pP4 = toPoly.pm("216*E2^9 + c57*E2^7 + c56*E2^6 + c55*E2^5 + c54*E2^4 +
	+ c53*E2^3 + c52*E2^2 + c51*E2 + c50")
#
# * (36*E2^5 + 6*c33*E2^4 + c57/6*E2^3 + c56/6*E2^2 + c55/6*E2 + c54/6)
pP2 = toPoly.pm("6*E2^4 - c33*E2^3 + c32*E2^2 + c31*E2 + c5")
pP4 = toPoly.pm("c53*E2^3 + c52*E2^2 + c51*E2 + c50")

# Overflows!
library(gmp)
source("Polynomials.Helper.BigNumbers.R")
pP1$coeff = as.bigz(pP1$coeff); pP2$coeff = as.bigz(pP2$coeff);
pP3$coeff = as.bigz(pP3$coeff); pP4$coeff = as.bigz(pP4$coeff);
pR = solve.lpm(pP2, pP4, xn=c("E2"), stop.at=1, asBigNum = TRUE)
str(pR)
toCoeff(pR[[1]][[2]], "E2")


c1 = S^3 + b*S - 4*R;
c20 = 3*R*b*S - 6*R^2;
c2 = c1^2 + 3*c20;
c3 = - 3*b^3 - S*b*c1;
c4 = (R*c1^3 + 9*R^2*c1^2 + 27*R^3*c1 - 27*S^3*R^3 + 27*R^4) / 9;
c5 = - b*c2 - 3*S^2*c3;
c33 = 18*S^2 + 3*b;
c32 = 27*S*R + 12*b^2 + 9*S*c1;
# c31 = 3*c3 - 27*S^3*R - 4*c2 + 3*b^3 + 3*S^3*c1 + 3*S*b*c1;
c31 = - 9*b^3 - 27*S^3*R - 4*c2 + 3*b^3 + 3*S^3*c1;
#
c410 = 27*(3*R^2*S^2 - b^4)/9;
c411 = 27*(c1*R + 3*R^2)/9;
c412 = 81*R*S/9;
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
	- 108*R*S*b*c3 + 18*c3^2 + 576*c4 - 432*b*c401 + 27*b^3*c403 + 432*R*S*c410 + 18*b^2*c410 +
	- 48*S*c1*c410 - 216*R*S*b*c411 + 24*S*b*c1*c411 - 48*c3*c411 + 27*R*S*b^2*c412 +
	- 3*S*b^2*c1*c412 + 24*b*c3*c412;
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

E2div = 36*c51^3 - 72*c50*c51*c52 + 36*c50^2*c53 - 6*c52^2*c53*c5 + 6*c51*c53^2*c5 +
	- 6*c52^3*c31 + 18*c51*c52*c53*c31 - 12*c50*c53^2*c31 + c53^3*c31^2 + 6*c51*c52^2*c32 +
	- 12*c51^2*c53*c32 - c53^3*c5*c32 - c52*c53^2*c31*c32 + c51*c53^2*c32^2 +
	+ 6*c51^2*c52*c33 - 6*c50*c52^2*c33 - 6*c50*c51*c53*c33 - c52*c53^2*c5*c33 +
	- c52^2*c53*c31*c33 + 2*c51*c53^2*c31*c33 + c51*c52*c53*c32*c33 - c50*c53^2*c32*c33 +
	+ c51^2*c53*c33^2 - c50*c52*c53*c33^2;

E2x0 = 36*c50*c51^2 - 36*c50^2*c52 - 6*c52^3*c5 + 12*c51*c52*c53*c5 - 6*c50*c53^2*c5 +
	+ 6*c50*c52*c53*c31 + c53^3*c5*c31 + 6*c50*c52^2*c32 - 12*c50*c51*c53*c32 +
	- c52*c53^2*c5*c32 + c50*c53^2*c32^2 + 6*c50*c51*c52*c33 - 6*c50^2*c53*c33 +
	- c52^2*c53*c5*c33 + c51*c53^2*c5*c33 + c50*c53^2*c31*c33 + c50*c52*c53*c32*c33 +
	+ c50*c51*c53*c33^2;


#############################
#############################

