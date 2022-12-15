########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S5: Hetero-Symmetric
### Useful Formulas
###
### draft v.0.1j


### Derivations
# - Basic derivations;
# - Numerical approaches: basic intuition;


####################

### Helper Functions

### Solver Tools
source("Polynomials.Helper.Solvers.Num.R")

# - is loaded automatically in "Solvers.Num.R";
# source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


### Other

# - check if 2 solutions are permutations of each other;
is.perm.S5 = function(s1, s2, tol=1E-6) {
	if(length(s1) != 5 || length(s2) != 5) stop("Not solutions of S5!");
	# E11b is different:
	id  = c(3,4,5,1,2);
	e21 = sum(s1 * s1[id]);
	e22 = sum(s2 * s2[id]);
	d = round0(e21 - e22, tol=tol);
	return(d == 0);
}
which.perm.S5 = function(s, tol=1E-6, verbose=TRUE) {
	nr = nrow(s);
	if(nr <= 1) return(array(0, c(0, 2)));
	# ID of permuted solutions;
	id = numeric(0);
	for(i1 in seq(nr, 2, by=-1)) {
		for(i2 in seq(i1 - 1)) {
			if(is.perm.S5(s[i1,], s[i2,])) {
				id = c(id, i2, i1); break;
			}
		}
	}
	id = matrix(id, nrow=2);
	if(ncol(id) == 0) {
		if(verbose) cat("No duplicates!\n");
		return(invisible(id));
	}
	return(id);
}

as.conj.S5 = function(x, nrow, rm.rows=0) {
	id = seq(7);
	if(any(rm.rows > 0)) {
		id = id[-rm.rows];
	}
	r = x[id,];
	r = rbind(r, as.conj(x[nrow,]));
	return(r);
}

neg = function(x, id=1) {
	x[id] = - x[id];
	return(x);
}

# Compare specific coefficients of E11b;
cmpE11b = function(p1, p2, n=NULL, xn="E11b") {
	if( ! is.null(n)) {
		id = match(xn, names(p1));
		if( ! is.na(id)) {
			p1 = p1[p1[,id] == n, - id, drop=FALSE];
		}
		id = match(xn, names(p2));
		if( ! is.na(id)) {
			p2 = p2[p2[,id] == n, - id, drop=FALSE];
		}
	}
	pR = diff.pm(p1, p2);
	return(pR);
}


#######################
#######################

### Base-System

### Numerical Solution:

# Helper Functions: for Baseline system
source("Poly.System.S5.Ht.Formulas.Derivation.HS0.R")

# Functions with the coefficients
source("Poly.System.S5.Ht.Formulas.Derivation.Coeffs.R")

### List of Initial values
# - is loaded automatically in file HS0.R;
# source("Poly.System.S5.Ht.Formulas.Derivation.x0.R")


###################
###################

# Note:
# - the double permutation (x1, x3) & (x4, x5) like c(x3, x2, x1, x5, x4),
#   and equivalent permutations (e.g. (x1, x2) & (x3, x5)),
#   are also a solution;

# - Basic solution moved to file:
#   Poly.System.S5.Ht.Formulas.Derivation.SolS0.R;

###################

### Case 3:
# S = 1; E11a = 0;

### Cases 4 - 6:
# S =  1; E11a =  1;
# S =  1; E11a = -1;
# S = -1; E11a =  1;

### List of Initial values
# - moved to file:
#   Poly.System.S5.Ht.Formulas.Derivation.x0.R;


### Examples

### Using given Path:
# - deriving solution for: R2 = c(-2,-1,0,0,2.5)
#   using a path from R2 = c(-5/3,-1,0,0,2.5);
R2 = c(-2,-1,0,0,2.5)
path = lapply(c(-5/3, -1.8, -2), function(R1) c(R1, -1,0,0,2.5));
x0 = x0All$Vn1fn12f
x.all = solve.path(solve.S5HtMixed.Num, x0, path=path, debug=T)
x.all = x.all[-2,]

### Example plot.path:
tmp = plot.path.S5(c(5,3,1,0,2), c(1,3,1,0,2), "E3V131", steps=101)
tmp = plot.path.S5(c(-5,3,1,0,2), c(-1,3,1,0,2), "E3Vn131", steps=101)


### Ex 2:
R2 = c(-1,1,0,0,2)
x0 = x0All$Vn11
x.all = solve.all(solve.S5HtMixed.Num, x0, R=R2, debug=F)
poly.calc(apply(x.all, 1, function(x) sum(x * x[c(3,4,5,1,2)]))) * 27

### E3:
R2 = c(0,0,1,0,2)
path = lapply(c(1/3, 2/3, 1), function(x) c(0, 1 - x, x,0,2));
x0 = x0All$V01
x.all = solve.path(solve.S5HtMixed.Num, x0, path=path, debug=T)
which.perm.S5(x.all) # all roots OK;

x.all = polyS(R2, "E3V001")
x.all = polyS(R2, "E3V011")


R2 = c(0,0,0,1,2)
x.all = polyS(R2, "E4V0001")

max.pow.S(c("E4V1011", "E4Vn1011", "E4V2011", "E4Vn2011"), pow=0, FUN=f0, R=c(1,0,1,1,2), R2=c(3,0,1,1,2), skip.path=T)
max.pow.S(c("E4V1011", "E4V101n1", "E4V1011", "E4V101n1"), pow=0, FUN=f0, npos=4, R=c(1,0,1,1,2), R2=c(1,0,1,1.3,2), skip.path=T)


fn = function(R) {
	Rn = R; Rn[1] = - Rn[1];
	R4 = R; R4[4] = - R4[4]; R4n = R4; R4n[1] = - R4n[1];
	(f0(R) - f0(Rn) - (f0(R4) - f0(R4n)))/4;
}

cc = list(c(1,0,1,1,2), c(2,0,1.1,1.3,2), c(1.3,0,0.9,1,2.1),
	c(2.1,0,0.9,1.1,2.3), c(2.2,0,1,0.8,1.9), c(3,0,1.2,0.9,2.2), c(3.3,0,4/5,1.2,2.2));
tmp = sapply(cc, function(R) polyS(R, "E4V1011"));
c0 = c() / 2;
c1 = c() / 2;
c0 = (c0 - c1)/2;


fc = function(R) { S = R[1]; E3 = R[3]; E4 = R[4]; E5 = R[5];
	c() / E5^2; }
m = sapply(cc, fc)
solve(t(m), c0 - sapply(cc, fn))


27*(E11a^7 + E11b^7)*E5^2 +
	# x^6
	- (E11a^6 + E11b^6)*(27*E5^2*S^2 - 4*E4^3 + 18*E3*E4*E5) - (E11a*E11b)^6 +
	+ 81*E11a*E11b*(E11a^5 + E11b^5)*E5^2 + 9*(E11a*E11b)^3*(E11a^3 + E11b^3)*E5*S +
	# x^5:
	+ (E11a^5 + E11b^5)*(9*E5^2*S^4 - E4^3*S^2 + 36*E3*E4*E5*S^2 - 90*E3*E5^2*S - 3*E4^2*E5*S +
		+ 4*E3^3*E5 - 225*E4*E5^2 - E3^2*E4^2) +
	+ (E11a*E11b)*(E11a^4 + E11b^4)*(9*E3^2*E5*S) +
	+ (E11a*E11b)^2*(E11a^3 + E11b^3)*(27*E5^2 - 18*E3*E5*S^2) +
	- (E11a*E11b)^3*(E11a^2 + E11b^2)*(2*E5*S^3 + 15*E3*E5) +
	- (E11a*E11b)^4*(E11a + E11b)*(3*E5*S + 2*E3^2) +
	+ 4*(E11a*E11b)^5*E3*S +
	# x^4: Note: 1 term from x^5 contributes as well;
	+ (E11a^4 + E11b^4)*( - E5^2*S^6 - 16*E3*E4*E5*S^4 + 48*E3*E5^2*S^3 - 31*E4^2*E5*S^3 +
		- 10*E3^3*E5*S^2 + 345*E4*E5^2*S^2 - 12*E3^2*E4^2*S^2 - 13*E3*E4^3*S + 5*E3^2*E4*E5*S +
		+ 150*E3^2*E5^2 - 33*E4^4 + 155*E3*E4^2*E5) +
	- (E11a*E11b)*(E11a^3 + E11b^3)*(21*E5^2*S^4 + 11*E3^3*E5 - 7*E3^2*E5*S^3 + 270*E3*E5^2*S) +
	+ (E11a*E11b)^2*(E11a^2 + E11b^2)*(18*E5^2*S^2 - E3^4 + 4*E3*E5*S^4 + 36*E3^2*E5*S) +
	- (E11a*E11b)^3*(E11a + E11b)*(10*E5^2 - 6*E3^3*S) +
	+ (E11a*E11b)^4*(4*E5*S^3 + 20*E3*E5 - 6*E3^2*S^2) +
	# x^3:
	+ (E11a^3 + E11b^3)*(2*E3*E4*E5*S^6 - 6*E3*E5^2*S^5 + 16*E4^2*E5*S^5 +
		+ 2*E3^3*E5*S^4 - 122*E4*E5^2*S^4 + 7*E3^2*E4^2*S^4 +
		- 50*E5^3*S^3 + 27*E3*E4^3*S^3 - E3^2*E4*E5*S^3 +
		+ 5*E3^2*E5^2*S^2 + 29*E4^4*S^2 - 189*E3*E4^2*E5*S^2 + 7*E3^4*E4*S^2 +
		+ 15*E4^3*E5*S + 6*E3^4*E5*S + 27*E3^3*E4^2*S +
		+ 375*E4^2*E5^2 + 29*E3^2*E4^3 - 105*E3^3*E4*E5) +
	+ (E11a*E11b)*(E11a^2 + E11b^2)*(4*E5^2*S^6 - 2*E3^2*E5*S^5 + 68*E3*E5^2*S^3 +
		- 27*E3^3*E5*S^2 + 2*E3^5*S - 375*E5^3*S + 275*E3^2*E5^2) +
	+ (E11a*E11b)^2*(E11a + E11b)*(12*E5^2*S^4 - 3*E3^3*E5 +
		+ 140*E3*E5^2*S - 6*E3^4*S^2 + 11*E3^2*E5*S^3) +
	- (E11a*E11b)^3*(4*E3^4 + 68*E5^2*S^2 + 8*E3*E5*S^4 - 4*E3^3*S^3 + 26*E3^2*E5*S) +
	# x^2:
	+ (E11a^2 + E11b^2)*(- 2*E4^2*E5*S^7 - E3^2*E4^2*S^6 + 12*E4*E5^2*S^6 +
		+ 14*E5^3*S^5 + 2*E3^2*E4*E5*S^5 - 14*E3*E4^3*S^5 +
		- 21*E3^2*E5^2*S^4 - 17*E4^4*S^4 - 2*E3^4*E4*S^4 + 102*E3*E4^2*E5*S^4 +
		+ 8*E3^4*E5*S^3 + 59*E4^3*E5*S^3 - 29*E3^3*E4^2*S^3 - 280*E3*E4*E5^2*S^3 +
		- E3^6*S^2 + 500*E3*E5^3*S^2 - 500*E4^2*E5^2*S^2 - 29*E3^2*E4^3*S^2 + 87*E3^3*E4*E5*S^2 +
		- 200*E3^3*E5^2*S + 625*E4*E5^3*S - 14*E3^5*E4*S + 115*E3^2*E4^2*E5*S - 11*E3*E4^4*S +
		+ 12*E3^5*E5 - 5^5*E5^4 + 54*E4^5 - 17*E3^4*E4^2 - 250*E3^2*E4*E5^2 - 275*E3*E4^3*E5) +
	+ (E11a*E11b)*(E11a + E11b)*(6*E3*E5^2*S^5 - 8*E3^3*E5*S^4 + 200*E5^3*S^3 + 2*E3^5*S^3 +
		- 20*E3^2*E5^2*S^2 + 30*E3^4*E5*S + 5^4*E3*E5^3 - 2*E3^6) +
	- (E11a*E11b)^2*(6*E5^2*S^6 - 4*E3^2*E5*S^5 + E3^4*S^4 + 60*E3*E5^2*S^3 + 24*E3^3*E5*S^2 +
		+ 750*E5^3*S - 8*E3^5*S + 375*E3^2*E5^2) +
	# x^1:
	+ (E11a + E11b)*(2*E3*E4^3*S^7 + 7*E4^4*S^6 - 18*E3*E4^2*E5*S^6 +
		- 62*E4^3*E5*S^5 + 6*E3^3*E4^2*S^5 + 86*E3*E4*E5^2*S^5 +
		- 150*E3*E5^3*S^4 + 420*E4^2*E5^2*S^4 + 23*E3^2*E4^3*S^4 - 44*E3^3*E4*E5*S^4 +
		+ 110*E3^3*E5^2*S^3 - 1250*E4*E5^3*S^3 - 3*E3*E4^4*S^3 - 108*E3^2*E4^2*E5*S^3 + 6*E3^5*E4*S^3 +
		+ 5^5*E5^4*S^2 - 26*E3^5*E5*S^2 - 35*E4^5*S^2 + 23*E3^4*E4^2*S^2 +
			+ 475*E3^2*E4*E5^2*S^2 + 185*E3*E4^3*E5*S^2 +
		+ 2*E3^7*S - 1250*E3^2*E5^3*S + 100*E4^4*E5*S - 375*E3*E4^2*E5^2*S +
			- 3*E3^3*E4^3*S - 40*E3^4*E4*E5*S +
		+ 125*E3^4*E5^2 - 625*E4^3*E5^2 +
		+ 7*E3^6*E4 + 175*E3^3*E4^2*E5 - 35*E3^2*E4^4 + 5^5*E3*E4*E5^3) +
	- E11a*E11b*(28*E5^3*S^5 + 44*E3^2*E5^2*S^4 - 28*E3^4*E5*S^3 +
		+ 4*E3^6*S^2 - 250*E3*E5^3*S^2 + 34*E3^5*E5 + 3*5^5*E5^4) +
	# B0:
	- E4^4*S^8 + 12*E4^3*E5*S^7 - 86*E4^2*E5^2*S^6 - 4*E3^2*E4^3*S^6 +
		+ 300*E4*E5^3*S^5 - 2*E3*E4^4*S^5 + 44*E3^2*E4^2*E5*S^5 +
		+ 10*E4^5*S^4 - 5^4*E5^4*S^4 - 6*E3^4*E4^2*S^4 - 220*E3^2*E4*E5^2*S^4 - 38*E3*E4^3*E5*S^4 +
		+ 500*E3^2*E5^3*S^3 - 60*E4^4*E5*S^3 + 250*E3*E4^2*E5^2*S^3 - 4*E3^3*E4^3*S^3 + 52*E3^4*E4*E5*S^3 +
		- 150*E3^4*E5^2*S^2 + 250*E4^3*E5^2*S^2 + 19*E3^2*E4^4*S^2 - 80*E3^3*E4^2*E5*S^2 +
			- 4*E3^6*E4*S^2 - 1250*E3*E4*E5^3*S^2 +
		+ 20*E3^6*E5*S - 2*E3^5*E4^2*S + 10*E3*E4^5*S - 150*E3^2*E4^3*E5*S + 500*E3^3*E4*E5^2*S +
	- E3^8 - 25*E4^6 - 50*E3^5*E4*E5 + 10*E3^4*E4^3 - 625*E3^2*E4^2*E5^2 + 250*E3*E4^4*E5 # = 0



max.pow.S(c("E3V101", "E3Vn101", "E3V501", "E3Vn501"), pow=0, FUN=f0, skip.path=T, R=c(1.1,0,3/4,0,2), R2=c(5,0,3/4,0,2))

solve.coeff(c(1.1,0,3/4,0,2), c(5,0,3/4,0,2), c(-2967.058 + 4468.348, -1493365 + 1634008) / 2,
"c(E3^6*S/E5, E3^2*E5*S^3)", function(R) { Rn = R; Rn[1] = - Rn[1]; (f0(R) - f0(Rn))/2; })

# before f3 was updated;
solve.coeff(c(1,0,1,0,2.2), c(5,0,1,0,2.2), c(- 107.3636 - 119.1818, - 31793.18 - 33179.55) / 2,
	"c(E3^4*S/E5, E3*S^5)", function(R) { Rn = R; Rn[1] = - Rn[1]; (f3(R) - f3(Rn))/2; })
solve.coeff(c(1,0,1,0,2), c(5,0,1,0,2), c(- 97 - 109, - 30485 - 31985) / 2,
	"c(E3^4*S/E5, E3*S^5)", function(R) { Rn = R; Rn[1] = - Rn[1]; (f3(R) - f3(Rn))/2; })

### Simple Examples:

### S = 0
R2 = c(0,1,0,0,2)
x0 = x0All$V01;

x.all = solve.all(solve.S5HtMixed.Num, x0, R=R2, debug=F)
poly.calc(apply(x.all, 1, function(x) sum(x * x[c(3,4,5,1,2)]))) * 27
# round0(poly.calc(x.all)) * 27

-12473 - 37419*x - 12473*x^2 - 10*x^3 - 10*x^4 + 27*x^5 + 80.75*x^6 + 27*x^7

### E11a = 0
R2 = c(1,0,0,0,2)
-2500 + 12500*x - 12472*x^2 - 100*x^3 - 1*x^4 + 9*x^5 - 27*x^6 + 27*x^7

R2 = c(1,0,0,0,3/2)
-1406.25 + 7031.25*x - 7010.25*x^2 - 75*x^3 - 1*x^4 + 9*x^5 - 27*x^6 + 27*x^7

### Examples: S & E11a
R2 = c(1,1/3,0,0,2)
277.185 - 2*x - 12505.22*x^2 - 349.9733*x^3 - 6.351853*x^4 + 11.94444*x^5 + 0.1663215*x^6 + 27*x^7
R2 = c(1,-1/3,0,0,2)
-8048.84 + 25090.59*x - 12773*x^2 + 152.4053*x^3 + 8.401235*x^4 + 12.01852*x^5 - 54.16701*x^6 + 27*x^7
R2 = c(-1,1/3,0,0,2)
278.3705 + 2*x - 12494.56*x^2 + 350.0226*x^3 - 6.388888*x^4 + 12.05556*x^5 - 0.1670072*x^6 + 27*x^7

R2 = c(1,1,0,0,2)
-2564 - 25342*x - 13521*x^2 - 908.5*x^3 - 13.5*x^4 + 33.5*x^5 + 58.25*x^6 + 27*x^7
R2 = c(1,-1,0,0,2)
-27436 + 51262*x - 14399*x^2 + 721.5*x^3 + 51.5*x^4 + 35.5*x^5 - 112.75*x^6 + 27*x^7
R2 = c(-1,1,0,0,2)
-2420 - 24530*x - 11377*x^2 + 784.5*x^3 - 14.5*x^4 + 38.5*x^5 + 49.25*x^6 + 27*x^7
R2 = c(-1,-1,0,0,2)
-27692 + 48850*x - 10655*x^2 - 589.5*x^3 + 44.5*x^4 + 36.5*x^5 - 103.75*x^6 + 27*x^7
R2 = c(1,1,0,0,3)
-5725 - 56795*x - 29682*x^2 - 1334.667*x^3 - 13.66667*x^4 + 34.33333*x^5 + 56.88889*x^6 + 27*x^7
R2 = c(-1,1,0,0,3)
-5509 - 55577*x - 26466*x^2 + 1210.667*x^3 - 14.33333*x^4 + 37.66667*x^5 + 50.88889*x^6 + 27*x^7

### E3:
R2 = c(0,0,1,0,2)
-0.25 + 125*x - 12494*x^2 + 150*x^4 + 2*x^5 + 27*x^7

### E11a & E3:
R2 = c(0,1,1,0,2)
-12190.25 - 35792*x - 11594.25*x^2 + 255*x^3 + 143.75*x^4 + 21*x^5 + 80.75*x^6 + 27*x^7

### Other:
R2 = c(5,1/3,0,0,2)
x0 = x0All$V50f
x.all = solve.all(solve.S5HtMixed.Num, x0, R=R2, debug=F)
poly.calc(apply(x.all, 1, function(x) sum(x * x[c(3,4,5,1,2)]))) * 27

###
R2 = c(5/4, -1/3,0,0,2)
x0 = x0All$V10
x.all = solve.all(solve.S5HtMixed.Num, x0, R=R2, debug=F)
poly.calc(apply(x.all, 1, function(x) sum(x * x[c(3,4,5,1,2)]))) * 27

######################
######################

### Robust Derivation:


x = x.all[1,]; E3 = E4 = 0;
x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4]; x5 = x[5];
s1 = x1 + x2; p1 = x1 * x2;
s2 = x3 + x4 + x5; e2 = (x3 + x4)*x5 + x3*x4; e3 = x3*x4*x5;
S = s1 + s2; E5 = p1*e3;
E11a = x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1;
E11b = x1*x3 + x2*x4 + x3*x5 + x4*x1 + x5*x2;
E2 = E11a + E11b;


###
### Transformed P[5] System:
s1 + s2 - S # = 0
#
s1*S - s1^2 + p1 + e2 - E2 # = 0
s1*e2 - s1*p1 + S*p1 + e3 - E3 # = 0
s1*e3 + p1*e2 - E4 # = 0
p1*e3 - E5 # = 0
x5^3 + s1*x5^2 + e2*x5 - e3 - S*x5^2 # = 0
#
s1*S - s1^2 + p1 + e2 - E2 # = 0
s1*p1*e2 - s1*p1^2 + S*p1^2 - p1*E3 + E5 # = 0
s1*E5 + p1^2*e2 - p1*E4 # = 0
p1*x5^3 + p1*s1*x5^2 + p1*e2*x5 - p1*S*x5^2 - E5 # = 0


p1 = toPoly.pm("s1*S - s1^2 + p1 + e2 - E2")
p2 = toPoly.pm("s1*p1*e2 - s1*p1^2 + S*p1^2 - p1*E3 + E5")
p3 = toPoly.pm("s1*E5 + p1^2*e2 - p1*E4")
p4 = toPoly.pm("p1*x5^3 + p1*s1*x5^2 + p1*e2*x5 - p1*S*x5^2 - E5")

pR1 = solve.lpm(p1, p4, p2, xn=c("p1", "e2"))
pR2 = solve.lpm(p1, p4, p3, xn=c("p1", "e2"))
pR1 = pR1[[2]]$Rez; pR1$coeff = - pR1$coeff;
pR2 = pR2[[2]]$Rez; pR2$coeff = - pR2$coeff;
table(pR2$s1)

tmp = gcd.pm(pR1, pR2, by="s1")
pR2 = diff.pm(pR2, mult.pm(pR1, toPoly.pm("x5^3")))

# Note: coeff a == 0!
x5^2*(S - x5)*(x5^5 - S*x5^4 + E2*x5^3 - E3*x5^2 + E4*x5 - E5)*s1^2 +
	- x5^2*(S - x5)^2*(x5^5 - S*x5^4 + E2*x5^3 - E3*x5^2 + E4*x5 - E5)*s1 +
	- E5^2 + 2*E4*E5*x5 - E4^2*x5^2 - E5*E2*S*x5^2 + E5*E2*x5^3 + E4*E2*S*x5^3 + E5*S^2*x5^3 - E4*E2*x5^4 +
	- 2*E5*S*x5^4 - E4*S^2*x5^4 + E5*x5^5 + 2*E4*S*x5^5 + E2^2*S*x5^5 - E4*x5^6 - 2*E2*S^2*x5^6 +
	+ 2*E2*S*x5^7 + S^3*x5^7 - 2*S^2*x5^8 + S*x5^9 - E2*S*x5^4*E3 - E2*x5^5*E3 + S^2*x5^5*E3 +
	- x5^7*E3 + x5^4*E3^2 # = 0

