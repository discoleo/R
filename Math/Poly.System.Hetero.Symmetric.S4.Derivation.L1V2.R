########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
###  == Derivation ==
###  Type: L1 V2
###
### draft v.0.1f-clean


### Types:

# V2a: x1^n + b*x1*x2 = R
# V2b: x1^n + b*x2*x3 = R
# V2c: x1^n + b*x1*x3 = R

### Note:
# - V2c: (x1, x3) and (x2, x4) are actually independent systems;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")

### Solver:
xi.f = function(x, R, b, n=2) {
	(R - x^n) / b[1];
}
xip.f = function(x, R, b, n=2, p=1) {
	(R - x^n) / b[1] / x^p;
}

### En:
E2af = function(x, n=2) {
	if(is.null(dim(x))) {
		x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4];
	} else {
		stop("TODO");
	}
	if(length(n) == 1) {
		E2a = x1^n*x2 + x2^n*x3 + x3^n*x4 + x4^n*x1;
	} else {
		p = n[2];
		n = n[1];
		E2a = x1^n*x2^p + x2^n*x3^p + x3^n*x4^p + x4^n*x1^p;
	}
	return(E2a);
}
E3af = function(x, n=2) {
	if(is.null(dim(x))) {
		x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4];
	} else {
		stop("TODO");
	}
	if(length(n) == 1) {
		E3a = x1^n*x2*x3 + x2^n*x3*x4 + x3^n*x4*x1 + x4^n*x1*x2;
	} else {
		p1 = n[2]; p2 = n[3];
		n  = n[1];
		E3a = x1^n*x2^p1*x3^p2 + x2^n*x3^p1*x4^p2 + x3^n*x4^p1*x1^p2 + x4^n*x1^p1*x2^p2;
	}
	return(E3a);
}
### Debug:
debug.E = function(x) {
	if(is.matrix(x)) {
		x1 = x[,1]; x2 = x[,2]; x3 = x[,3]; x4 = x[,4];
	} else {
		x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4];
	}
	s1 = x1 + x3; s2 = x2 + x4;
	p1 = x1 * x3; p2 = x2 * x4;
	sp = p1 + p2; ps = s1 * s2;
	S = s1 + s2; E4 = p1 * p2;
	E2 = sp + ps;
	E3 = p1*s2 + p2*s1;
	data.frame(S=S, E2=E2, E3=E3, E4=E4, E11a = ps);
}

### Classic Polynomial
polyGen.S4Ht.V2a = function(n, factorize=TRUE) {
	pR = toPoly.pm("x4^n + b*x4*x1 - R");
	pSub = data.frame(x3=c(n, 0), R=c(0,1), coeff=c(1,-1));
	pDiv = data.frame(x3=1, b=1, coeff = -1);
	for(i in 4:2) {
		xn = paste0("x", i);
		pR = replace.fr.pm(pR, pSub, pDiv, xn);
		if(i > 1) {
			xn = paste0("x", i - 2);
			names(pSub)[1] = xn;
			names(pDiv)[1] = xn;
		}
	}
	if(factorize) {
		pDiv = toPoly.pm("(R - x1^n)^n + b^n*x1^n*(R - x1^n) - b^n*R*x1^n");
		pR = div.pm(pR, pDiv, c("x1", "b"));
		pR = pR$Rez;
	}
	return(pR);
}

### Helper Equations:

### E21a + E12a
E21a + E12a - s2*(x1^2 + x3^2) - s1*(x2^2 + x4^2) # = 0
E21a + E12a - s2*(s1^2 - 2*p1) - s1*(s2^2 - 2*p2) # = 0
E21a + E12a - s1*s2*S + 2*(p1*s2 + p2*s1) # = 0
E21a + E12a - E11a*S + 2*E3 # = 0

### E31a + E13a
E31a + E13a - s2*(x1^3 + x3^3) - s1*(x2^3 + x4^3) # = 0
E31a + E13a - s2*(s1^3 - 3*p1*s1) - s1*(s2^3 - 3*p2*s2) # = 0
E31a + E13a - s1*s2*(s1^2 + s2^2) + 3*s1*s2*(p1 + p2) # = 0
E31a + E13a - E11a*(S^2 - 2*E11a) + 3*E11a*(E2 - E11a) # = 0
E31a + E13a - E11a^2 + 3*E11a*E2 - E11a*S^2 # = 0


####################
####################

################
### Type V2a ###
################

###############
### Order 2 ###
###############

### V2a:
### x1^2 + b*x1*x2 = R
x1^2 + b*x1*x2 # - R
x2^2 + b*x2*x3 # - R
x3^2 + b*x3*x4 # - R
x4^2 + b*x4*x1 # - R

### Solution:

### Eq 1: Sum(Prod(x2*x3*x4*...)) =>
E4*S + b*E4*S - R*E3 # = 0
(b+1)*E4*S - R*E3 # = 0

### Eq 3:
E4^2 - b^4*E4^2 + 2*R^2*E4 + 2*R*E2*E4 - R*E3^2 - 2*R^2*E3*S +
	+ 2*R^3*E2 + R^2*E2^2 - R^3*S^2 + R^4 # = 0


### Helper Eqs:

### Eq H1: E11a
# Sum =>
S^2 - 2*E2 + b*E11a - 4*R # = 0
# b*E11a = - S^2 + 2*E2 + 4*R;

### Eq H2: E21a
# Sum(x1*...) =>
S^3 - 3*E2*S + 3*E3 + b*E21a - R*S # = 0
# b*E21a = - S^3 + 3*E2*S - 3*E3 + R*S;

### Eq H3: E12a
# Sum(x4*...) =>
E12a + b*E3 - R*S # = 0
# E12a = - b*E3 + R*S;

### Eq H4: E21a + b*E12a
# Sum(x2*...) =>
E21a + b*E12a - R*S # = 0

### Eq H5: E31a
# Sum(x1^2*...) =>
S4 + b*E31a - R*(S^2 - 2*E2) # = 0
b*E31a + 2*E2^2 - 4*E2*S^2 + 2*R*E2 + 4*E3*S - 4*E4 + S^4 - R*S^2 # = 0

### Eq H6: E31a (also)
# Sum(x1*x2*...) =>
E31a + b*E22a - R*E11a # = 0
b*E31a + b^2*(E22 - (E2 - E11a)^2 + 2*E4) + R*(S^2 - 2*E2 - 4*R) # = 0
b*E31a - (S^2 - 2*E2 - 4*R)^2 +
	+ 4*b^2*E4 - 2*b^2*E3*S+ 4*b*E2^2 - 2*b*E2*S^2 + (8*b - 2)*R*E2 + R*S^2 - 4*R^2 # = 0
b*E31a + 4*b^2*E4 - 2*b^2*E3*S +
	+ 4*(b - 1)*E2^2 - (2*b - 4)*E2*S^2 + (8*b - 18)*R*E2 +
	- S^4 + 9*R*S^2 - 20*R^2 # = 0

### Eq H7: E13a
# Sum(x2^2*...) =>
E22a + b*E13a - R*(S^2 - 2*E2) # = 0
b*E22a + b^2*E13a - b*R*S^2 + 2*b*R*E2 # = 0
# Diff: Eq H6 =>
b^2*E13a - E31a + R*E11a - b*R*S^2 + 2*b*R*E2 # = 0
b^3*E13a - b*E31a + 2*(b^2 + 1)*R*E2 - (b^2 + 1)*R*S^2 + 4*R^2 # = 0
# Subst: H 5 =>
b^3*E13a + 2*E2^2 - 4*E2*S^2 + 2*R*E2 + 4*E3*S - 4*E4 + S^4 - R*S^2 +
	+ 2*(b^2 + 1)*R*E2 - (b^2 + 1)*R*S^2 + 4*R^2 # = 0
b^3*E13a - 4*E4 + 4*E3*S + 2*E2^2 - 4*E2*S^2 + 2*(b^2 + 2)*R*E2 +
	+ S^4 - (b^2 + 2)*R*S^2 + 4*R^2 # = 0

### Eq H8: E211a
# Sum(x2^2*x3*x4*...) =>
E4*E11a + b*E4*(S^2 - 2*E2) - R*E211a # = 0

### Eq H9: E112a
# Sum(x3*x4*...) =>
E112a + 4*b*E4 - R*E11a # = 0


#######
### Eq:
S^3 - R*(b^4 - 2*b^3 + 4*b^2 - 4*b + 4)*S


### Debug:
R = 5; b = -2;
x1 = -5.0222909065 + 0.0000000000i;
x2 = -2.0133646503 + 0.0000000000i;
x3 =   0.2350202147 + 0.0000000000i;
x4 = -10.5198727342 + 0.0000000000i;

x = c(x1,x2,x3,x4)
s1 = x1 + x3; s2 = x2 + x4;
p1 = x1 * x3; p2 = x2 * x4;
sp = p1 + p2; ps = s1 * s2;
S = s1 + s2; E4 = p1 * p2;
E2 = sp + ps;
E3 = p1*s2 + p2*s1;
#
E11a = ps;
E21  = S*E2 - 3*E3;
E21a = E2af(x, n=2);
E12a = E2af(x, n=c(1,2));
E31a = E2af(x, n=3);
E13a = E2af(x, n=c(1,3));
#
E22  = E2^2 - 2*S*E3 + 2*E4;
E22a = E2af(x, n=c(2,2));
S4 = sum(x^4);
E121a = E3af(x, n=c(1,2,1));
E112a = E3af(x, n=c(1,1,2));
E211a = E3af(x, n=c(2,1,1));


### Classic Poly:
pR = polyGen.S4Ht.V2a(2)
toCoeff(pR, "x1")

(b^2 + 1)*x1^12 - (6 + 8*b^2 + 3*b^4 + b^6)*R*x1^10 +
	+ (15 + 22*b^2 + 13*b^4 + 6*b^6 + b^8)*R^2*x1^8 +
	- (20 + 28*b^2 + 18*b^4 + 10*b^6 + 3*b^8)*R^3*x1^6 +
	+ (15 + 17*b^2 + 9*b^4 + 5*b^6 + 2*b^8)*R^4*x1^4 +
	- (6 + 4*b^2 + b^4 + b^6)*R^5*x1^2 + R^6 # = 0

### Case: b = 1i
x1^8 - 3*R*x1^6 + 4*R^2*x1^4 - 2*R^3*x1^2 + R^4 # = 0


### [old]
# - moved from file:
#   Poly.System.Hetero.Symmetric.S4.Derivation.R;

### Eq 1: Sum(Prod(x2*x3*x4*...)) =>
E4*S + b*E4*S - R*E3 # = 0
(b+1)*E4*S - R*E3 # = 0

### Eq 2:
### Sum =>
S^2 - 2*E2 + b*E11a - 4*R # = 0

### x1^2 = R - b*x1*x2 => Prod =>
E4^2 - R^4 + b*R^3*E11a - b^2*R^2*(E121a + 2*E4) + b^3*R*E4*E11a - b^4*E4^2 # = 0
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
(E3^2 - 2*E4*E2) + b*E4*E11a - R*(E11a^2 - 2*E121a - 4*E4) # = 0
b^2*E3^2 - b^2*E4*S^2 + 8*b^2*R*E4 - R*(2*E2 - S^2 + 4*R)^2 + 2*b^2*R*E121a # = 0
b^2*E3^2 - b^2*E4*S^2 + 8*b^2*R*E4 - 4*R*E2^2 + 4*R*E2*S^2 - 16*R^2*E2 +
	- R*(S^4 - 8*R*S^2 + 16*R^2) + 2*b^2*R*E121a # = 0
(b^2 + 2)*E3^2 - (b^2 + 2)*E4*S^2 + 8*(b^2 + 1)*R*E4 +
	- 4*R*E2^2 + 4*R*E2*S^2 - 8*R^2*E2 +
	- R*S^4 + 4*R^2*S^2 - 8*R^3 # = 0


#######
### Eq:
S^3 - R*(b^4 - 2*b^3 + 4*b^2 - 4*b + 4)*S


### Auxiliary Eqs:
### Special Case: S = 0
E3 = 0; E2 = -2*R; E4 = R^2 / (b^2 + 1);
### Case: S != 0
E2 = -b*R*(b^2 - b + 2);
E3 = - (b+1) * R * S;
E4 = R*E3 / ((b+1)*S); # = - R^2;


### [old approach]
### Workout:
### Half-Elementary Polynomials

### E121a:
### x1^2 = R - b*x1*x2 => Sum(Prod(2 eqs)) =>
b^2*E121a - 2*b*R*E11a + 4*R^2 - (E11a^2 - 2*E121a - 4*E4) # = 0
(b^2 + 2)*E121a - 2*R*(2*E2 - S^2 + 4*R) - E11a^2 + 4*E4 + 4*R^2
b^2*(b^2 + 2)*E121a - 2*b^2*R*(2*E2 - S^2 + 4*R) - (2*E2 - S^2 + 4*R)^2 + 4*b^2*E4 + 4*b^2*R^2
b^2*(b^2 + 2)*E121a - 4*E2^2 + 4*E2*S^2 - 4*(b^2+4)*R*E2 +
	+ 4*b^2*E4 - S^4 + 2*(b^2+4)*R*S^2 - 4*(b^2+4)*R^2 # = 0


### Debug:

### Classic solver:
solve.S4P2V2a.classic = function(R, b, debug=FALSE) {
	x1 = roots(coeff.S4P2V2a(R, b)); # see below for coeff.S4P2V2();
	x2 = xip.f(x1, R, b, p=1);
	x3 = xip.f(x2, R, b, p=1);
	x4 = xip.f(x3, R, b, p=1);
	sol = cbind(x1, x2, x3, x4);
	return(sol);
}
coeff.S4P2V2a = function(R, b) {
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
coeff.S4P2V2a_P16 = function(R, b) {
	# P[16]: complete version;
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
#
R = -1; b = 2;
sol = solve.S4P2V2a.classic(R, b);
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];

E = debug.E(sol);
E = sapply(E, round0);
S = E$S; E2 = E$E2; E3 = E$E3; E4 = E$E4;
E2 = -b*R*(b^2 - b + 2);

### [old]
# TODO: clean;

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


####################
####################

################
### Type V2b ###
################

###############
### Order 2 ###
###############

### V2b:
### x1^2 + b*x2*x3 = R
x1^2 + b*x2*x3 # - R
x2^2 + b*x3*x4 # - R
x3^2 + b*x4*x1 # - R
x4^2 + b*x1*x2 # - R

### Solution:

### Eq 1: Sum(x1*...) =>
S^3 - 3*E2*S + (b + 3)*E3 - R*S # = 0

### Eq 2: Eq H2 + Eq H3 =>
E21 - (E21a + E12a) - (sp*S - E3) # = 0
S*E2 - 3*E3 - (E21a + E12a) - ((E2 - E11a)*S - E3) # = 0
b*(E21a + E12a) + S^3 - 2*E2*S - 4*R*S + 2*b*E3 # = 0
b*E21a + S^3 - 2*E2*S + (b - 4)*R*S - (b^2 - 2*b)*E3 # = 0
(b+1)*S^3 - 2*(b+1)*E2*S +
	+ (b^2 - 2*b - 4)*R*S - (b+1)*(b^2 - 2*b)*E3 # = 0
# Reduction =>
(b + 1)*E2*S - (b + 1)*(b^2 - b + 3)*E3 + (b^2 - b - 3)*R*S # = 0

### Eq 3: Eq H4 + Eq H5
E211a + E112a - ps*sp # = 0
E211a + E112a - E11a*(E2 - E11a) # = 0
b^3*E211a + b^3*E112a + b*(S^2 - 2*E2 - 4*R)*(S^2 + (b - 2)*E2 - 4*R) # = 0
b^2*(S^4 - R*S^2 + 2*E2^2 - 4*E2*S^2 + 2*R*E2 + 4*E3*S - 4*E4) +
	+ (4*(b - 1)*E2^2 - (2*b - 4)*E2*S^2 + (2*b^2 + 8*b - 16)*R*E2 - 2*b^2*E3*S +
		- S^4 - (b^2 - 8)*R*S^2 + 4*b^2*E4 - 16*R^2) +
	- b*(S^2 - 2*E2 - 4*R)*(S^2 + (b - 2)*E2 - 4*R) # = 0
4*(b^2 - 1)*E2^2 - (5*b^2 - 2*b - 4)*E2*S^2 + (8*b^2 - 8*b - 16)*R*E2 +
	+ 2*b^2*E3*S + (b^2 - b - 1)*S^4 - (2*b^2 - 8*b - 8)*R*S^2 - 16*(b + 1)*R^2 # = 0
# Reduction =>
4*(b^2 - 1)*E2^2 - (2*b^2 + b - 1)*E2*S^2 + (8*b^2 - 8*b - 16)*R*E2 +
	- (b^3 - 4*b - 3)*E3*S - (b^2 - 7*b - 7)*R*S^2 - 16*(b + 1)*R^2 # = 0

### Eq 4: Eq H6 + Eq H7
b^3*(E31a + E13a - ps*(s1^2 + s2^2 - 3*sp)) # = 0
b^3*(E31a + E13a - ps*(S^2 - 2*ps - 3*sp)) # = 0
b^3*(E31a + E13a - E11a*(S^2 - 2*E2 - (E2 - E11a))) # = 0
b^3*(E31a + E13a) - b^2*E11a*(b*S^2 - 3*b*E2 + b*E11a) # = 0
b^3*(E31a + E13a) + b*(S^2 - 2*E2 - 4*R)*((b - 1)*S^2 - (3*b - 2)*E2 + 4*R) # = 0
b^3*E31a + b^3*E13a +
	+ (6*b^2 - 4*b)*E2^2 - (5*b^2 - 4*b)*E2*S^2 + (12*b^2 - 16*b)*R*E2 +
	+ (b^2 - b)*S^4 - (4*b^2 - 8*b)*R*S^2 - 16*b*R^2 # = 0
b^3*E31a - 4*b^4*E4 +
	+ (6*b^2 - 4*b)*E2^2 - (5*b^2 - 4*b)*E2*S^2 + (14*b^2 - 16*b)*R*E2 +
	+ (b^2 - b)*S^4 - (5*b^2 - 8*b)*R*S^2 + (4*b^2 - 16*b)*R^2 # = 0
2*b^2*E3*S - 4*(b^4 + b^2)*E4 +
	+ (6*b^2 - 8*b + 4)*E2^2 - (5*b^2 - 6*b + 4)*E2*S^2 + (12*b^2 - 24*b + 16)*R*E2 +
	+ (b^2 - b + 1)*S^4 - (4*b^2 - 8*b + 8)*R*S^2 + (4*b^2 - 16*b + 16)*R^2 # = 0
# Reduction =>
4*(b^4 + b^2)*E4 + (b^3 - 2*b + 3)*E3*S +
	- (6*b^2 - 8*b + 4)*E2^2 + (2*b^2 - 3*b + 1)*E2*S^2 - (12*b^2 - 24*b + 16)*R*E2 +
	+ (3*b^2 - 7*b + 7)*R*S^2 - (4*b^2 - 16*b + 16)*R^2 # = 0


### Helper Eqs:

### Eq H1: E11a
# Sum =>
S^2 - 2*E2 + b*E11a - 4*R # = 0
# b*E11a = - S^2 + 2*E2 + 4*R;

### Eq H2: E21a
# Sum(x2*...) =>
(b+1)*E21a - R*S # = 0
# (b+1)*E21a = R*S;

### Eq H3: E12a
# Sum(x4*...) =>
E12a + b*E3 - R*S # = 0
# E12a = - b*E3 + R*S;

### Eq H4: E211a
# Sum(x1^2*...) =>
S4 + b*E211a - R*(S^2 - 2*E2) # = 0
b*E211a + S^4 - R*S^2 + 2*E2^2 - 4*E2*S^2 + 2*R*E2 + 4*E3*S - 4*E4 # = 0

### Eq H5: E112a
# Sum(x4^2*...) =>
E22a + b*E112a - R*(S^2 - 2*E2) # = 0
E22 - ((E2 - E11a)^2 - 2*E4) + b*E112a - R*(S^2 - 2*E2) # = 0
2*E2*E11a - 2*E3*S - E11a^2 + b*E112a - R*S^2 + 2*R*E2 + 4*E4 # = 0
b^3*E112a + 4*b*E2^2 - 2*b*E2*S^2 + (2*b^2 + 8*b)*R*E2 - 2*b^2*E3*S +
	- (S^2 - 2*E2 - 4*R)^2 - b^2*R*S^2 + 4*b^2*E4 # = 0
b^3*E112a + 4*(b - 1)*E2^2 - (2*b - 4)*E2*S^2 + (2*b^2 + 8*b - 16)*R*E2 - 2*b^2*E3*S +
	- S^4 - (b^2 - 8)*R*S^2 + 4*b^2*E4 - 16*R^2 # = 0


### Eq H6: E31a
# Sum(x2^2*...) =>
E22a + b*E31a - R*(S^2 - 2*E2) # = 0
E22 - E22b + b*E31a - R*(S^2 - 2*E2) # = 0
E22 - E22b + b*E31a - R*(S^2 - 2*E2) # = 0
E2^2 - 2*S*E3 + 4*E4 - (E2 - E11a)^2 + b*E31a - R*S^2 + 2*R*E2 # = 0
b*E31a - 2*S*E3 + 4*E4 - E11a^2 + 2*E2*E11a - R*S^2 + 2*R*E2 # = 0
b^3*E31a + 4*b*E2^2 - 2*b*E2*S^2 + (2*b^2 + 8*b)*R*E2 - 2*b^2*S*E3 + 4*b^2*E4 +
	- (S^2 - 2*E2 - 4*R)^2 - b^2*R*S^2 # = 0
b^3*E31a + 4*(b - 1)*E2^2 - (2*b - 4)*E2*S^2 + (2*b^2 + 8*b - 16)*R*E2 +
	- 2*b^2*S*E3 + 4*b^2*E4 +
	- S^4 - (b^2 - 8)*R*S^2 - 16*R^2 # = 0


### Eq H7: E13a
# Sum(x1*x4*...) =>
E13a + 4*b*E4 - R*E11a # = 0
b*E13a + 4*b^2*E4 - 2*R*E2 + R*S^2 - 4*R^2 # = 0


### Alternatives:

# - Diffs usually do NOT work with low powers;


### Debug:
R = 3; b = -2;
x1 =  0.6180339888 - 0.5257311121i;
x2 = -1.6180339887 + 0.8506508083i;
x3 =  0.6180339888 + 0.5257311121i;
x4 = -1.6180339887 - 0.8506508083i;

R = 3; b = -5;
x1 = -1.5115402800 + 0.8512706373i;
x2 =  0.2902350198 + 0.1770534310i;
x3 = -1.5115402799 - 0.8512706374i;
x4 =  0.2902350198 - 0.1770534310i;

x = c(x1,x2,x3,x4)
s1 = x1 + x3; s2 = x2 + x4;
p1 = x1 * x3; p2 = x2 * x4;
sp = p1 + p2; ps = s1 * s2;
S = s1 + s2; E4 = p1 * p2;
E2 = sp + ps;
E3 = p1*s2 + p2*s1;
#
E11a = ps;
E21  = S*E2 - 3*E3;
E21a = E2af(x, n=2);
E12a = E2af(x, n=c(1,2));
E31a = E2af(x, n=3);
E13a = E2af(x, n=c(1,3));
#
E22  = E2^2 - 2*S*E3 + 2*E4;
E22a = E2af(x, n=c(2,2));
E22b = (x1*x3)^2 + (x2*x4)^2;
E33a = E2af(x, n=c(3,3));
#
S4 = sum(x^4);
E211a = E3af(x, n=2);
E112a = E3af(x, n=c(1,1,2));

###
pP1 = toPoly.pm("S^3 - 3*E2*S + (b + 3)*E3 - R*S")
pP2 = toPoly.pm("(b + 1)*E2*S - (b + 1)*(b^2 - b + 3)*E3 + (b^2 - b - 3)*R*S")
pP3 = toPoly.pm("4*(b^2 - 1)*E2^2 - (2*b^2 + b - 1)*E2*S^2 + (8*b^2 - 8*b - 16)*R*E2 +
	- (b^3 - 4*b - 3)*E3*S - (b^2 - 7*b - 7)*R*S^2 - 16*(b + 1)*R^2")
# pP4 = toPoly.pm("") # not needed;

pR = solve.lpm(pP1, pP2, pP3, xn=c("E3", "E2"));
pR[[2]]$Rez = orderVars.pm(pR[[2]]$Rez, c("b", "R", "S"));
# TODO: debug div.pm;
# pR[[2]]$Rez = div.pm(pR[[2]]$Rez, toPoly.pm("(b + 3)*(b + 1)"), "b")$Rez

(b+1)^2*(b-1)*(b^2+1)*S^4 +
	- (b+1)*(3*b^4 + 22*b^3 - 12*b^2 + 12*b - 20)*R*S^2 +
	+ 16*(b^2 + 2*b + 2)*(3*b^2 - 2)*R^2 # = 0
# P[2] * P[2] =>
((b+1)*S^2 - 16*R) * ((b-1)*(b+1)*(b^2+1)*S^2 - (b^2 + 2*b + 2)*(3*b^2 - 2)*R)

