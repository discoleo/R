########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Asymmetric S2:
### Binomial Expansions
###
### draft v.0.2g-auto


### Asymmetric Polynomial Systems: 2 Variables
### Binomial Expansions


### Base Polynomials:
# - Class 1 Polynomials
# - Class 2 Polynomials
# - Class 3 Polynomials


###############
### History ###
###############


### draft v.0.2g - v.0.2g-auto:
# - another variant of Ht-system: Dual/Double-variant;
# - automatic generation of System & all roots; [v.0.2g-auto]
### draft v.0.2f:
# - reordered sections: non-Correlated Variants;
### draft v.0.2e - v.0.2e-var2:
# - Ht system: automatic Generator for Class 1 polynomials;
# - Generator: for Ht-SumDiff variant;
### draft v.0.2c - v.0.2d-varP-ex:
# - Ht system with Class 3 polynomials;
# - Ht system: automatic Generator for Class 3 polynomials & base-roots; [v.0.2d-sol]
# - Ht system: Class 3 variant based on Powers; [v.0.2d-varP]
#   & example; [v.0.2d-varP-ex]
### draft v.0.2b-ht - v.0.2b-sol:
# - Ht-variant for Class 1 Order 3;
# - some concrete & special cases; (v.0.2b-sp)
# - Entangled variants: Ht & Diff-type; (v.0.2b-vars)
# - more work on the Ht-solver; (v.0.2b-sol)
### draft v.0.2a:
# - System derived from Class 2 polynomials;
### draft v.0.1f:
# - solver for the Order 5 / Simple system;
### draft v.0.1e:
# - Multiplicative entanglement;
### draft v.0.1d - v.0.1d-Eq2:
# - Base Order 3 variants:
#   1st equation & 2nd Eq; (v.0.1d - v.0.1d-Eq2)
# - analysis of the variants; (v.0.1d-var)
### draft v.0.1c:
# - exact solution to the Order 3 system;
### draft v.0.1b:
# - Simple Order 5: Cardano-type;
### draft v.0.1a:
# - Order 3;


####################
####################

### helper functions

### fast load:
# source("Polynomials.Helper.R")
# source("Polynomials.Helper.Generators.R")

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R;
# Polynomials.Helper.Generators.R;
# e.g. round0(), round0.p;

solve.Cardano = function(c, d, n=3) {
	m = unity(n=n, all=TRUE);
	xdet = rootn(d^2 - c^n, 2);
	p = d + xdet; p = rootn(p, n);
	q = d - xdet; q = rootn(q, n);
	if(n %% 2 == 0 && c < 0) q = -q;
	sol = p*m + q/m;
	return(sol);
}
roots.Cl1 = function(K, s, n=3, debug=TRUE) {
	m = unity(n, all=TRUE);
	k = rootn(K, n=n);
	k = k * m;
	len = length(s);
	if(debug) {
		if(len > n) warning("Length of shifts > n!");
	}
	r = sapply(k, function(k) sum(s * k^seq(0, len-1)) );
	return(r);
}
roots.Cl3 = function(s, n=3) {
	div = 2*n + 1;
	cs  = 2*cos(seq(n) * (2*pi/div));
	len = length(s) %% (n+1);
	if(len != 0) s = c(s, rep(0, (n+1 - len)));
	s0 = s[seq(1, length(s), by=n+1)]; s = s[ - seq(1, length(s), by=n+1)];
	r = sapply(seq(n), function(id) sum(s0, s * shift(cs, by=id)) );
	return(r);
}
roots.Cl3P = function(s, n=3) {
	# based on Powers;
	div = 2*n + 1;
	cs  = 2*cos(seq(n) * (2*pi/div));
	len = length(s) %% (n+1);
	if(len != 0) s = c(s, rep(0, (n+1 - len)));
	s0 = s[seq(1, length(s), by=n+1)]; s = s[ - seq(1, length(s), by=n+1)];
	r = sapply(seq(n), function(id) sum(s0, s * cs[id]^seq(n)) );
	return(r);
}

### Generators

### Base: Class 1
system.S2Cl1Ht = function(K, s1, s2, n=3, type="Ht", tol=1E-10, withBase=FALSE, debug=TRUE) {
	type = pmatch(type, c("Ht", "Sum", "HtSumDiff", "HtDual"));
	if(is.na(type)) stop("Unsupported type!");
	if(length(K) > 1) stop("Parameter K must have only 1 value!")
	FUN = toPoly.Class1S.pm;
	p.gen = function(s, xn="x") {
		p = FUN(b=s, n=n, kn="K", xn=xn);
		p = replace.pm(p, K, "K");
		return(p);
	}
	s = round0(s1 + s2, tol=tol);
	if(debug) print(s);
	pS = p.gen(s);
	pS = replace.pm(pS, data.frame(x=1:0, y=0:1, coeff=1), "x");
	#
	px = p.gen(s1);
	py = p.gen(s2, xn="y");
	#
	diff.p = function() {
		s.d = round0(s1 - s2);
		pD = p.gen(s.d);
		pD = replace.pm(pD, data.frame(x=1:0, y=0:1, coeff=c(1,-1)), "x");
	}
	if(type == 1) {
		p1 = diff.pm(pS, py);
		p2 = diff.pm(pS, px);
	} else if(type == 2) {
		p1 = pS;
		p2 = diff.pm(sum.pm(px, py), pS);
	} else if(type == 3) {
		pD = diff.p();
		p1 = diff.pm(pS, py);
		p2 = diff.pm(px, pD);
	} else if(type == 4) {
		pD = diff.p();
		p1 = diff.pm(sum.pm(pS, pD), px);
		p2 = diff.pm(diff.pm(pS, pD), py);
	}
	rez = list(p1=p1, p2=p2);
	if(withBase) rez = c(rez, list(px=px, py=py));
	return(rez);
}

### Base: Class 3
system.S2Cl3Ht = function(s1, s2, n=3, type="Simple", tol=1E-10, debug=TRUE) {
	type = pmatch(type, c("Simple", "Powers"));
	if(is.na(type)) stop("Unsupported type!");
	FUN = if(type == 1) toPoly.Class3.pm else toPoly.Class3P.pm;
	p.gen = function(s, xn="x") {
		s.id = which(s != 0);
		s = s[s.id];
		p = FUN(n, s.id = (s.id - 1), xn=xn, sn = "s");
		p = replace.pm(p, s, paste0("s", (s.id - 1)) );
		return(p);
	}
	s = round0(s1 + s2, tol=tol);
	if(debug) print(s);
	pS = p.gen(s);
	pS = replace.pm(pS, data.frame(x=1:0, y=0:1, coeff=1), "x");
	#
	px = p.gen(s1);
	py = p.gen(s2, xn="y");
	#
	p1 = diff.pm(pS, py);
	p2 = diff.pm(pS, px);
	rez = list(p1=p1, p2=p2);
	return(rez);
}

### Classic Polynomial

clPoly.S2Cl1 = function(K, s1, s2, n=3, div=NULL, type="Ht", tol=1E-10) {
	pL = system.S2Cl1Ht(K=K, s1=s1, s2=s2, n=n, type=type, tol=tol, withBase=TRUE);
	pR = solve.pm(pL[[1]], pL[[2]], "y");
	xgcd = gcd.vpm(pR$Rez);
	pR$Rez$coeff = pR$Rez$coeff / xgcd;
	pR2 = divOK.pm(pR$Rez, pL$px, "x", warn=TRUE);
	if(pR2$isDiv) {
		if( ! is.null(div)) {
			pR2 = divOK.pm(pR2$Rez, div, "x", warn=TRUE);
		}
		pR$Rez = pR2$Rez;
	}
	return(pR);
}
sysAll.S2Cl1 = function(K, s1, s2, n=3, div=toPoly.pm("x^2-2*x+1"), type="Ht", allRoots=TRUE) {
	# Px
	p = system.S2Cl1Ht(K, s1, s2, n=n, type=type);
	# roots
	x = roots.Cl1(K, s1, n=n);
	y = roots.Cl1(K, s2, n=n);
	sol = cbind(x=x, y=y);
	# remaining roots:
	if(allRoots) {
		px = clPoly.S2Cl1(K, s1, s2, n=n, div=div, type=type);
		x = roots(evalCoeff(px$Rez, "x", c(), c()));
		y = sapply(x, function(x) eval.pm(px$x0, x, "x"));
		y = y / sapply(x, function(x) eval.pm(px$div, x, "x"));
		sol = rbind(sol, cbind(x=x, y=y));
	}
	# px = the classic polynomial;
	sol = c(p, list(px=px$Rez, sol=sol));
	return(sol)
}


##########################

##########################
### Polynomial Systems ###
##########################

#####################
### Base: Class 1 ###
#####################

# Binomial expansions of Class 1 polynomials;

###############
### Order 3 ###
###############

### Base:
# x = k^2 + b11*k
# y = k^2 + b21*k
# where k^3 = K;

### Simple System:
x^3 + y^3 - 3*b11*K*x - 3*b21*K*y - 2*K^2 - (b11^3 + b21^3)*K # = 0
x*y*(x+y) - (b11 + 2*b21)*K*x - (2*b11 + b21)*K*y - 2*K^2 - b11*b21*(b11+b21)*K # = 0
# - linked to the Binomial expansion of S^3, S = (x+y);

### Entangled Variants:

### 1.) Ht-Variant:
x^3 + 3*x*y*(x+y) - 6*(b11+b21)*K*(x+y) + 3*b21*K*y - 7*K^2 + b21^3*K - (b11+b21)^3*K # = 0
y^3 + 3*x*y*(x+y) - 6*(b11+b21)*K*(x+y) + 3*b11*K*x - 7*K^2 + b11^3*K - (b11+b21)^3*K # = 0

### 2.) Diff-Variant
x^3 - y^3 - 3*b11*K*x + 3*b21*K*y - (b11^3 - b21^3)*K # = 0
x*y*(x+y) - (b11 + 2*b21)*K*x - (2*b11 + b21)*K*y - 2*K^2 - b11*b21*(b11+b21)*K # = 0
# variant: y => -y;

### 3.) Other Variants:
# - see section on Multiplicative variants;


### Derivation:

### P(x^3) + P(y^3)
x^3 + y^3 - 3*b11*K*x - 3*b21*K*y - 2*K^2 - (b11^3 + b21^3)*K # = 0

### (x+y)^3 - (x^3+y^3)
3*x*y*(x+y) - 6*(b11+b21)*K*(x+y) + 3*b11*K*x + 3*b21*K*y - 6*K^2 - 3*b11*b21*(b11+b21)*K # = 0
x*y*(x+y) - 2*(b11+b21)*K*(x+y) + b11*K*x + b21*K*y - 2*K^2 - b11*b21*(b11+b21)*K # = 0
x*y*(x+y) - (b11 + 2*b21)*K*x - (2*b11 + b21)*K*y - 2*K^2 - b11*b21*(b11+b21)*K # = 0

### Note:
### Eq 1 + 3*Eq 2 => S = (x+y) =>
#   S^3 - 6*(b11+b21)*K*S + ... = 0;
# - enables solving for S = (x+y);
### Step 2:
# x + y = S;
# x*y*S - b11*K*x - b21*K*y = ... - S^3 / 3;


### Solver:
solve.DP3 = function(K, b, type="Simple", all=TRUE) {
	# parameter "all" = not used anymore;
	d.f = function(b) 4*K^2 + b^3*K / 2;
	# can also use the direct formulas;
	bs = b[1] + b[2];
	S = solve.Cardano(2*bs*K, d.f(bs), n=3);
	type = pmatch(type, c("Simple", "Ht", "Diff"));
	if(is.na(type)) stop("Type NOT supported!");
	if(type == 1) {
		x = sapply(S, function(S) {
			roots(c(-S, S^2 + (b[1] - b[2])*K, -2*K^2 - b[1]*b[2]*bs*K - (bs + b[1])*K*S))
		})
		S = rep(S, each=2);
		y = S - x;
	} else if(type == 2) {
		# TODO: correct solution?
		if(FALSE) {
			x = sapply(S, function(S) {
				coeff = c(-3*S, 3*S^2 + 3*(b[1]-b[2])*K, S^3 + 3*b[2]*K*S - 12*bs*K*S - 14*K^2 - (bs^3 + 3*b[1]*b[2]*bs)*K);
				roots(coeff);
			})
			S = rep(S, each=2);
			y = S - x;
		} else {
			# only 3 of the roots!
			m = unity(3, all=TRUE); k = rootn(K, 3);
			x = (k*m)^2 + b[1]*k*m;
			y = (k*m)^2 + b[2]*k*m; # direct, but only 3 roots;
			# y = solve.y.HtP3(x, K, b); # all roots;
		}
	}
	sol = cbind(x=as.vector(x), y=as.vector(y));
	return(sol);
}
solve.y.HtP3 = function(x, K, b) {
	bs = b[1] + b[2]; b11 = b[1]; b21 = b[2];
	y0 = - 6*x^5 - 3*K*(7*bs - 8*b11)*x^3 - 3*K*(7*K + bs^3 - 3*b11^3 + 2*b21^3)*x^2 +
		+ 18*bs*K^2*(bs + b11)*x +
		+ 21*(bs + b11)*K^3 + 3*(2*bs^4 - bs^3*b21 - 2*bs*b21^3 + b21^4)*K^2;
	div = - 6*x^4 + 9*K*(2*bs + b21)*x^2 - 3*K*(7*K + bs^3 - b21^3)*x - 9*K^2*(bs + b11)^2;
	return(y0 / div);
}
### [old]
# simple variant of solver (non-robust);
solve.DP3.old = function(K, b, all=FALSE) {
	c.f = function(b) b*K;
	d.f = function(b) (K^2 + b^3*K) / 2;
	# can also use the direct formulas;
	x = solve.Cardano(c.f(b[1]), d.f(b[1]), n=3);
	y = solve.Cardano(c.f(b[2]), d.f(b[2]), n=3);
	if(all) {
		# is actually NOT correct
		sol = expand.grid(x, y);
	} else sol = cbind(x=x, y=y);
	return(sol);
}

### Examples:

### Ex 1:
b = c(1,-2)
K = 3
#
sol = solve.DP3(K, b);
x = sol[,1]; y = sol[,2];
b11 = b[1]; b21 = b[2];
# concrete Example
round0(x^3 + y^3 - 3*K*x + 6*K*y - 2*K^2 + 7*K) # = 0
round0(x*y*(x+y) + 3*K*x - 2*K^2 - 2*K) # = 0

### Ht-Variant / concrete:
sol = solve.DP3(K, b, type="Ht");
x = sol[,1]; y = sol[,2];
#
x^3 + 3*x*y*(x+y) + 6*K*(x+y) - 6*K*y - 7*K^2 - 7*K # = 0
y^3 + 3*x*y*(x+y) + 6*K*(x+y) + 3*K*x - 7*K^2 + 2*K # = 0
### Classic Poly:
x = roots(c(28, 0, - 378, 558, 5103, - 11826, 2943, 0, - 381024, - 592704));
y = solve.y.HtP3(x, K, b);
round0.p(poly.calc(x[c(1:4, 6,7)]) * 28);


p2  = toPoly.pm("y^3 + 3*x^2*y + 3*x*y^2 - 6*bs*K*x - 6*bs*K*y + 3*b11*K*x - 7*K^2 + b11^3*K - bs^3*K")
p2s = p2
p2s = replace.pm(p2s, c(b[1], b[1]+b[2], K), c("b11", "bs", "K"))
p2s = shift.pm(p2s, c(x[5], y[5]), c("x", "y"))
print.p(round.pm(p2s),  c("x", "y"))


### Ex 2:
b = c(-2,3)
K = 3
#
sol = solve.DP3(K, b, all=T);
x = sol[,1]; y = sol[,2];
b11 = b[1]; b21 = b[2];
### Test
# [with non-robust solver]
# - only 3 solutions are correct, but not necessarily the base-set;
x^3 + y^3 - 3*b11*K*x - 3*b21*K*y # = ...
2*K^2 + (b11^3 + b21^3)*K
err = x*y*(x+y) - (b11 + 2*b21)*K*x - (2*b11 + b21)*K*y - 2*K^2 - b11*b21*(b11+b21)*K
round0(err)


### Ex 3:
n = 3
K = 3
s1 = c(1,-1,2); s2 = c(0,2,-1);
x = roots.Cl1(K, s1, n=n);
y = roots.Cl1(K, s2, n=n);
p = system.S2Cl1Ht(K, s1, s2, n=n);
print.p(p[[1]], c("x","y"))
print.p(p[[2]], c("x","y"))
x^3 + 3*x^2*y + 3*x*y^2 - 3*x^2 - 6*x*y - 3*y^2 - 6*x - 24*y + 11 # = 0
y^3 + 3*x^2*y + 3*x*y^2 - 6*x*y - 3*y^2 - 27*x - 6*y + 84 # = 0


### Ex 4:
n = 3
K = 3
s1 = c(1,-1,2); s2 = c(0,2,-1);
x = roots.Cl1(K, s1, n=n);
y = roots.Cl1(K, s2, n=n);
p = system.S2Cl1Ht(K, s1, s2, n=n, type="HtSumDiff");
print.p(p[[1]], c("x","y"))
print.p(p[[2]], c("x","y"))
x^3 + 3*x^2*y + 3*x*y^2 - 3*x^2 - 6*x*y - 3*y^2 - 6*x - 24*y + 11 # = 0
y^3 + 3*x^2*y - 3*x*y^2 - 6*x*y + 3*y^2 - 63*x + 84*y + 156 # = 0


### Ex 5:
n = 3
K = 3
s1 = c(1,-1,2); s2 = c(0,2,-1);
x = roots.Cl1(K, s1, n=n);
y = roots.Cl1(K, s2, n=n);
p = system.S2Cl1Ht(K, s1, s2, n=n, type="HtDual");
print.p(p[[1]], c("x","y"))
print.p(p[[2]], c("x","y"))
x^3 + 6*x*y^2 - 3*x^2 - 6*y^2 + 57*x - 90*y - 160 # = 0
y^3 + 6*x^2*y - 12*x*y - 90*x + 60*y + 255 # = 0

# remaining roots:
px = clPoly.S2Cl1(K, s1, s2, n=n, div=toPoly.pm("x^2-2*x+1"), type="HtDual");
x = roots(evalCoeff(px$Rez, "x", c(), c()));
y = sapply(x, function(x) eval.pm(px$x0, x, "x"));
y = y / sapply(x, function(x) eval.pm(px$div, x, "x"));


### Ex 6:
n = 3
K = 3
s1 = c(1,-1, 2)
s2 = c(1,-1, 4)
p = sysAll.S2Cl1(K, s1, s2, n=n, type="Ht")
sol = p$sol; x = sol[,1]; y = sol[,2];
#
print.p(p[[1]], c("x","y"))
print.p(p[[2]], c("x","y"))


#############
### Variants:

### Multiplicative Variants

### Base:
# x = k^2 + b11*k
# y = k^2 + b21*k
# where k^3 = K;

### System:
x^3 + y^3 - 3*b11*K*x - 3*b21*K*y - (2*K^2 + (b11^3 + b21^3)*K) # = 0
(b11+b21)*(x*y)^2 - K*((b11+b21)^2 + 2*b11*b21)*x*y +
	- b11*K*(K + b21^3)*x - b21*K*(K + b11^3)*y # = 0

### Special Case:
# b21 = - b11
x^3 + y^3 - 3*b11*K*x + 3*b11*K*y - 2*K^2 # = 0
2*b11*x*y - (K - b11^3)*x + (K + b11^3)*y # = 0

### Solver:
test.SP3M = function(sol, K, b) {
	b11 = b[1]; b21 = b[2];
	if(is.matrix(sol)) {x = sol[,1]; y = sol[,2];}
	else {x = sol[1]; y = sol[2];}
	err1 = x^3 + y^3 - 3*b11*K*x - 3*b21*K*y - (2*K^2 + (b11^3 + b21^3)*K);
	err2 = (b11+b21)*(x*y)^2 - K*((b11+b21)^2 + 2*b11*b21)*x*y +
		- b11*K*(K + b21^3)*x - b21*K*(K + b11^3)*y
	err = rbind(err1, err2); err = round0(err);
	return(err);
}

### Derivation:

### Prod =>
# (x*y) = b11*b21*k^2 + K*k + (b11+b21)*K;
# =>
(x*y)^3 - 3*(b11+b21)*K*(x*y)^2 + 3*(b11+b21)^2*K^2*x*y - 3*b11*b21*K^2*x*y +
	- K^4 + 3*b11*b21*(b11+b21)*K^3 - (b11+b21)^3*K^3 - (b11*b21)^3*K^2 # = 0
(3*b11*K*x + K^2 + b11^3*K)*(3*b21*K*y + K^2 + b21^3*K) - 3*(b11+b21)*K*(x*y)^2 + 3*(b11+b21)^2*K^2*x*y - 3*b11*b21*K^2*x*y +
	- K^4 + 3*b11*b21*(b11+b21)*K^3 - (b11+b21)^3*K^3 - (b11*b21)^3*K^2
K*(3*b11*x + K + b11^3)*(3*b21*y + K + b21^3) - 3*(b11+b21)*(x*y)^2 + 3*(b11+b21)^2*K*x*y - 3*b11*b21*K*x*y +
	- K^3 + 3*b11*b21*(b11+b21)*K^2 - (b11+b21)^3*K^2 - (b11*b21)^3*K
K*(9*b11*b21*x*y + 3*b11*K*x + 3*b11*b21^3*x + 3*b21*K*y + 3*b21*b11^3*y + K^2 + (b11^3+b21^3)*K + b11^3*b21^3) +
	- 3*(b11+b21)*(x*y)^2 + 3*((b11+b21)^2 - b11*b21)*K*x*y +
	- K^3 + 3*b11*b21*(b11+b21)*K^2 - (b11+b21)^3*K^2 - (b11*b21)^3*K
(b11+b21)*(x*y)^2 - ((b11+b21)^2 + 2*b11*b21)*K*x*y +
	- b11*K*(K + b21^3)*x - b21*K*(K + b11^3)*y # = 0


### Examples:

### Ex 1:
b = c(1,-2)
K = 3
#
b11 = b[1]; b21 = b[2];
k = rootn(K, 3);
x = k^2 + b[1]*k;
y = k^2 + b[2]*k;
test.SP3M(c(x, y), K, b)


### Ex 2: Special case
b = c(2,-2)
K = 3
#
b11 = b[1]; b21 = b[2];
k = rootn(K, 3);
x = k^2 + b[1]*k;
y = k^2 + b[2]*k;
test.SP3M(c(x, y), K, b)


#############
### Variants:

### non-Correlated Variants

### Base
# A^3 - 3*c1*A - 2*d1 = 0
# B^3 - 3*c2*B - 2*d2 = 0
# x = A + B
# y = A - B

# Note:
# - Simple variants seem to be reducible;
#   (at least for base-order 3)

### System
x^4 - y^4 - x*y*(x^2-y^2) - 3*(c1+3*c2)*(x^2-y^2) - 4*(d1+3*d2)*x + 4*(d1-3*d2)*y # = 0
x^4 - y^4 + x*y*(x^2-y^2) - 3*(3*c1+c2)*(x^2-y^2) - 4*(3*d1+d2)*x + 4*(3*d1-d2)*y # = 0

### Transformed
x^4 - y^4 - 6*(c1+c2)*(x^2-y^2) - 8*(d1+d2)*x + 8*(d1-d2)*y # = 0
x*y*(x^2-y^2) - 3*(c1-c2)*(x^2-y^2) - 4*(d1-d2)*x + 4*(d1+d2)*y # = 0
# Case: (x,y) != 0 => (x != y) =>
x^3 + y^3 + 3*x*y*(x+y) - 12*c1*(x+y) - 16*d1 # = 0

### Transformed: Variant
x^4 - y^4 - 6*(c1+c2)*(x^2-y^2) - 8*(d1+d2)*x + 8*(d1-d2)*y # = 0
x^3 + y^3 + 3*x*y*(x+y) + x*y*(x^2-y^2) - 3*(c1-c2)*(x^2-y^2) - 4*(d1-d2+3*c1)*x + 4*(d1+d2-3*c1)*y - 16*d1 # = 0
# Sum(Eq 1 + 2*Eq 2) =>
(x - y + 2)*(x^3 + y^3 + 3*x*y*(x+y) - 12*c1*(x+y) - 16*d1) # = 0


### Solution:
### Trivial solution
# x = y = 0;
### Non-trivial solution:
# see above & Derivation;

### Derivation:
### Eq 1: x^3 + Y^3:
x^3 + y^3 - 2*A^3 - 6*A*B^2 # = 0
x^3 + y^3 - 2*(3*c1*A + 2*d1) - 6*A*B^2 # = 0
# [simple/redundant]
4*x^3 + 4*y^3 - 12*c1*(x+y) - 16*d1 - 3*(x+y)*(x-y)^2 # = 0
4*x^3 + 4*y^3 - 12*c1*(x+y) - 16*d1 - 3*(x+y)*(x^2 + y^2 - 2*x*y) # = 0
x^3 + y^3 + 3*x*y*(x + y) - 12*c1*(x+y) - 16*d1 # = 0
# variant:
B*(x^3 + y^3) - 2*B*(3*c1*A + 2*d1) - 6*A*B^3 # = 0
(x-y)*(x^3 + y^3) - (x-y)*(3*c1*(x+y) + 4*d1) - 3*(x+y)*(3*c2*(x-y) + 4*d2) # = 0
x^4 - y^4 - x*y*(x^2-y^2) - 3*c1*(x^2-y^2) - 4*d1*(x-y) - 9*c2*(x^2-y^2) - 12*d2*(x+y) # = 0
x^4 - y^4 - x*y*(x^2-y^2) - 3*(c1+3*c2)*(x^2-y^2) - 4*(d1+3*d2)*x + 4*(d1-3*d2)*y # = 0

### Eq 2: x^3 - y^3
x^3 - y^3 - 2*B^3 - 6*A^2*B # = 0
x^3 - y^3 - 2*(3*c2*B + 2*d2) - 6*A^2*B # = 0
A*(x^3 - y^3) - 2*A*(3*c2*B + 2*d2) - 6*A^3*B # = 0
(x+y)*(x^3 - y^3) - 2*(x+y)*(3*c2*B + 2*d2) - 6*(3*c1*A + 2*d1)*(x-y) # = 0
(x+y)*(x^3 - y^3) - 2*(x+y)*(3*c2*B + 2*d2) - 3*(3*c1*(x+y) + 4*d1)*(x-y) # = 0
x^4 - y^4 + x*y*(x^2-y^2) - 9*c1*(x^2-y^2) - (x+y)*(3*c2*(x-y) + 4*d2) - 12*d1*(x-y) # = 0
x^4 - y^4 + x*y*(x^2-y^2) - 3*(3*c1+c2)*(x^2-y^2) - 4*(3*d1+d2)*x + 4*(3*d1-d2)*y # = 0


### Examples
c = c(1, 2)
d  = c(2, 1)
#
c1 = c[1]; c2 = c[2]; d1 = d[1]; d2 = d[2];
A = solve.Cardano(c[1], d[1], n=3)
B = solve.Cardano(c[2], d[2], n=3)
x = A + B; y = A - B;
### "all" roots
sol = expand.grid(A, B);
x = sol[,1] + sol[,2]; y = sol[,1] - sol[,2];


### Test


#####################
#####################

###############
### Order 5 ###
###############

### Simple case:
### Base:
# x = k^4 + b11*k
# y = k^4 + b21*k
# where k^5 = K;

### System:
# x^5 + y^5 - 5*b11*K*x^3 - 5*b21*K*y^3 + 5*b11^2*K^2*x + 5*b21^2*K^2*y = 2*K^4 + (b11^5 + b21^5)*K
# Eq 2: see derivation;


### Derivation:

### P(x^3) + P(y^3)
x^5 + y^5 - 5*b11*K*x^3 - 5*b21*K*y^3 + 5*b11^2*K^2*x + 5*b21^2*K^2*y - 2*K^4 - (b11^5 + b21^5)*K # = 0

### (x+y)^5 - (x^5+y^5)
5*x*y*(x^3 + y^3) + 10*x^2*y^2*(x+y) - 10*(b11+b21)*K*(x+y)^3 + 20*(b11+b21)^2*K^2*(x+y) +
	+  5*b11*K*x^3 + 5*b21*K*y^3 - 5*b11^2*K^2*x - 5*b21^2*K^2*y - 30*K^4 - (b11+b21)^5*K + (b11^5 + b21^5)*K # = 0
5*x*y*(x^3 + y^3) + 10*x^2*y^2*(x+y) - 10*(b11+b21)*K*(x^3+y^3) - 30*(b11+b21)*K*x*y*(x+y) + 20*(b11+b21)^2*K^2*(x+y) +
	+  5*b11*K*x^3 + 5*b21*K*y^3 - 5*b11^2*K^2*x - 5*b21^2*K^2*y - 30*K^4 - (b11+b21)^5*K + (b11^5 + b21^5)*K # = 0
x*y*(x^3 + y^3) + 2*x^2*y^2*(x+y) - (b11+2*b21)*K*x^3 - (2*b11+b21)*K*y^3 - 6*(b11+b21)*K*x*y*(x+y) +
	+ (3*b11^2 + 4*b21^2 + 8*b11*b21)*K^2*x + (4*b11^2 + 3*b21^2 + 8*b11*b21)*K^2*y +
	- 6*K^4 - b11*b21*((b11+b21)^2 - b11*b21)*(b11+b21)*K # = 0

### Solution:
x*y*(S^3 - 3*x*y*S) + 2*x^2*y^2*S - (b11+2*b21)*K*x^3 - (2*b11+b21)*K*y^3 - 6*(b11+b21)*K*x*y*S +
	+ (3*b11^2 + 4*b21^2 + 8*b11*b21)*K^2*x + (4*b11^2 + 3*b21^2 + 8*b11*b21)*K^2*y +
	- 6*K^4 - b11*b21*((b11+b21)^2 - b11*b21)*(b11+b21)*K # = 0
S*x^4 - (2*S^2 + (b11-b21)*K)*x^3 + (2*S^3 - 3*b21*K*S)*x^2 +
	- S^4*x + (b11^2 - b21^2)*K^2*x + 3*b21*K*S^2*x +
	- (4*b11^2 + 3*b21^2 + 8*b11*b21)*K^2*S +
	+ 6*K^4 + b11*b21*((b11+b21)^2 - b11*b21)*(b11+b21)*K + (2*b11+b21)*K*S^3 # = 0


### Solver:
solve.DP5 = function(K, b) {
	bs = b[1] + b[2]; bd = b[1] - b[2]; bp = b[1]*b[2];
	S = solve.Cardano(c=2*bs*K, d=(16*K^4 + bs^5*K/2), n=5);
	x = sapply(S, function(S) {
		coeff = c(S, - (2*S^2 + bd*K), (2*S^3 - 3*b[2]*K*S),
			- S^4 + bs*bd*K^2 + 3*b[2]*K*S^2,
			- (b[1]^2 + 3*bs^2 + 2*bp)*K^2*S +
				+ 6*K^4 + bp*(bs^2 - bp)*bs*K + (b[1]+bs)*K*S^3);
			return(roots(coeff));
	})
	S = rep(S, each=4);
	y = S - x;
	sol = cbind(x = as.vector(x), y = as.vector(y));
	return(sol);
}
test.DP5 = function(x, y, K, b) {
	b11 = b[1]; b21 = b[2];
	err1 = x^5 + y^5 - 5*b11*K*x^3 - 5*b21*K*y^3 + 5*b11^2*K^2*x + 5*b21^2*K^2*y - 2*K^4 - (b11^5 + b21^5)*K;
	err2 = x*y*(x^3 + y^3) + 2*x^2*y^2*(x+y) - (b11+2*b21)*K*x^3 - (2*b11+b21)*K*y^3 - 6*(b11+b21)*K*x*y*(x+y) +
		+ (3*b11^2 + 4*b21^2 + 8*b11*b21)*K^2*x + (4*b11^2 + 3*b21^2 + 8*b11*b21)*K^2*y +
		- 6*K^4 - b11*b21*((b11+b21)^2 - b11*b21)*(b11+b21)*K;
	err = rbind(err1, err2);
	err = round0(err);
	return(err);
}

### Examples
b = c(-1, 2)
K = 3
#
n = 5
b11 = b[1]; b21 = b[2];
k = rootn(K, n);
x = k^4 + b[1]*k;
y = k^4 + b[2]*k;
S = x + y;
### concrete Example:
x^5 + y^5 + 5*K*x^3 - 10*K*y^3 + 5*K^2*x + 20*K^2*y - 2*K^4 - 31*K # = 0
x*y*(x^3 + y^3) + 2*x^2*y^2*(x+y) - 3*K*x^3 - 6*K*x*y*(x+y) + 3*K^2*x - 6*K^4 + 6*K # = 0


### Ex 1:
b = c(-1, 2)
K = 3
sol = solve.DP5(K, b);
x = sol[,1]; y = sol[,2];

### Test
test.DP5(x, y, K, b)


### Classic Polynomial
p1 = toPoly.pm("x^5 + y^5 + 5*K*x^3 - 10*K*y^3 + 5*K^2*x + 20*K^2*y - 2*K^4 - 31*K");
p2 = toPoly.pm("x^4*y + x*y^4 + 2*x^3*y^2 + 2*x^2*y^3 - 3*K*x^3 - 6*K*x^2*y - 6*K*x*y^2 + 3*K^2*x - 6*K^4 + 6*K");
# Big-Integers
library(gmp)
p1$coeff = as.bigz(p1$coeff);
p2$coeff = as.bigz(p2$coeff);
pR = solve.pm(p1, p2, "y", asBigNum=TRUE);
str(pR) # pR$Rez = polynomial with 221 terms;
max(pR$Rez$x) # the question to the universe!
# Q: Are there also other roots?
# Or: Are the remaining "roots" False-roots?
# Solver: only 20 roots;
# TODO: factor the big polynomial;


round0.p(poly.calc(x) * (32*K^3 + 1))
round0(poly.calc.mpfr(x, tol=1E-4) * (32*K^3 + 1))


#####################
#####################

#####################
### Base: Class 2 ###
#####################

# Binomial expansions of Class 2 polynomials;

###############
### Order 4 ###
###############

### Base:
# x = s2*m^2 + s1*m
# y = s3*m^3 - s1*m
# where m^5 = 1;

### Simple System:
# x^4 + y^4 = ...;
# (x + y)^4 = ...;

### Entangled Variant:
# x^4 + 4*x*y*(x^2+y^2) + 6*x^2*y^2 = ...;
# y^4 + 4*x*y*(x^2+y^2) + 6*x^2*y^2 = ...;


### Solution:

### Simple system:
# Step 1: compute S = x + y;
# Step 2:
# - solve:
#   x + y = S;
#   the 1st Eq;


### Examples:

m = unity(5, all=FALSE);

s = c(3, 2, -1)
s1 = s[1]; s2 = s[2]; s3 = s[3];
x = s2*m^2 + s1*m;
y = s3*m^3 - s1*m;

### Test
### Eq 1:
x^4 + (s1 + s2)*x^3 + (s1^2 + 2*s1*s2 + s2^2)*x^2 + (s1^3 + 3*s1^2*s2 - 2*s1*s2^2 + s2^3)*x +
	+ s1^4 - s1^3*s2 + s1^2*s2^2 - s1*s2^3 + s2^4 +
y^4 - (s1 - s3)*y^3 + (s1^2 - 2*s1*s3 + s3^2)*y^2 - (s1^3 + 2*s1^2*s3 + 3*s1*s3^2 - s3^3)*y +
	+ s1^4 + s1^3*s3 + s1^2*s3^2 + s1*s3^3 + s3^4 # = 0
### Eq 2:
x^4 + y^4 + 4*x^3*y + 4*x*y^3 + 6*x^2*y^2 + 3*(s2+s3)*x*y*(x+y) + 2*(s2^2 - 3*s2*s3 + s3^2)*y*x +
	+ (s2 + s3) * (x^3 + y^3) + (s2^2 - 3*s2*s3 + s3^2) * (x^2 + y^2) +
	+ (s2^3 - 2*s2^2*s3 - 2*s2*s3^2 + s3^3) * (x + y) +
	+ s2^4 + s3^4 - s2^3*s3 + s2^2*s3^2 - s2*s3^3 # = 0


### Derivation:
p1 = toPoly.Class2.pm(4, s.id=c(1,2));
print.p(p1, "x")

p2 = toPoly.Class2.pm(4, s.id=c(1,3), xn="y");
p2 = replace.pm(p2, toPoly.pm("-s1"), "s1", 1);
print.p(p2, "y")

p3 = toPoly.Class2.pm(4, s.id=c(2,3));
p3 = replace.pm(p3, toPoly.pm("x+y"), "x", 1);
print.p(p3, c("x", "y"))


###################
###################

###############
### Class 3 ###
###############

### Base:
# x = s2*m[2] + s1*m[1] + s0;
# y = s2*m[2] - s1*m[1] - s0;
# where m = 2*cos(2*id*pi/7), id = seq(3);

### System:
x^3 + 3*x^2*y + 3*x*y^2 - (3*s0*y^2 - s1*y^2 + s2*y^2 + 3*s0^2*y - 2*s0*s1*y - 2*s1^2*y + 2*s0*s2*y - 3*s1*s2*y - 2*s2^2*y +
	+ s0^3 - s0^2*s1 - 2*s0*s1^2 + s1^3 + s0^2*s2 - 3*s0*s1*s2 - 3*s1^2*s2 - 2*s0*s2^2 - 4*s1*s2^2 - s2^3) +
	+ 2*s2*x^2 + 4*s2*x*y + 2*s2*y^2 - 8*s2^2*x - 8*s2^2*y - 8*s2^3 # = 0
y^3 + 3*x^2*y + 3*x*y^2 - (- 3*s0*x^2 + s1*x^2 + s2*x^2 + 3*s0^2*x - 2*s0*s1*x - 2*s1^2*x - 2*s0*s2*x + 3*s1*s2*x - 2*s2^2*x +
	- s0^3 + s0^2*s1 + 2*s0*s1^2 - s1^3 + s0^2*s2 - 3*s0*s1*s2 - 3*s1^2*s2 + 2*s0*s2^2 + 4*s1*s2^2 - s2^3) +
	+ 2*s2*x^2 + 4*s2*x*y + 2*s2*y^2 - 8*s2^2*x - 8*s2^2*y - 8*s2^3 # = 0


### Examples:
m = 2*cos(2*pi/7 * (1:3));

s = c(1, -2, 3, 0);
s0 = s[1]; s1 = s[2]; s2 = s[3];
x = roots.Cl3(s, n=3);
y = roots.Cl3(s * c(-1,-1,1,0), n=3);

### concrete Example:
x^3 + 3*x^2*y + 3*x*y^2 + 6*x^2 + 12*x*y - 2*y^2 - 72*x - 77*y - 215 # = 0
y^3 + 3*x^2*y + 3*x*y^2 + 8*x^2 + 12*x*y + 6*y^2 - 29*x - 72*y - 133 # = 0


### Derivation
s = c(1, -2, 3, 0);
p = system.S2Cl3Ht(s, s*c(-1,-1,1,0), n=3)
print.p(p[[1]], c("x","y"))
print.p(p[[2]], c("x","y"))

### Base-Test:
x^3 - 3*s0*x^2 + s1*x^2 + s2*x^2 +
	+ 3*s0^2*x - 2*s0*s1*x - 2*s1^2*x - 2*s0*s2*x + 3*s1*s2*x - 2*s2^2*x +
	- s0^3 + s0^2*s1 + 2*s0*s1^2 - s1^3 + s0^2*s2 - 3*s0*s1*s2 - 3*s1^2*s2 + 2*s0*s2^2 + 4*s1*s2^2 - s2^3


### Ex 2:
s1 = c(1, -2, 3, 0);
s2 = c(0, -1,-2, 0);
x = roots.Cl3(s1, n=3);
y = roots.Cl3(s2, n=3);
p = system.S2Cl3Ht(s1, s2, n=3)
print.p(p[[1]], c("x","y"))
print.p(p[[2]], c("x","y"))
### Test:
x^3 + 3*x^2*y + 3*x*y^2 - 5*x^2 - 10*x*y - 2*y^2 - 22*x - 18*y + 14 # = 0
y^3 + 3*x^2*y + 3*x*y^2 - 3*x^2 - 10*x*y - 5*y^2 + 21*x - 22*y + 96 # = 0


### Ex 3:
n = 3
s1 = c(1, -2, 3);
s2 = c(0, -1,-2);
x = roots.Cl3P(s1, n=n);
y = roots.Cl3P(s2, n=n);
p = system.S2Cl3Ht(s1, s2, n=n, type="Powers")
print.p(p[[1]], c("x","y"))
print.p(p[[2]], c("x","y"))
### Test:
x^3 + 3*x^2*y + 3*x*y^2 - 11*x^2 - 22*x*y - 20*y^2 + 10*x - 10*y + 30 # = 0
y^3 + 3*x^2*y + 3*x*y^2 + 9*x^2 - 22*x*y - 11*y^2 - 79*x + 10*y + 142 # = 0


### Ex 4:
n = 5
s1 = c(1, -2, 3,-1);
s2 = c(0, -1,-2, 1);
x = roots.Cl3P(s1, n=n);
y = roots.Cl3P(s2, n=n);
p = system.S2Cl3Ht(s1, s2, n=n, type="Powers")
print.p(p[[1]], c("x","y"))
print.p(p[[2]], c("x","y"))
### Test:
x^5 + 5*x^4*y + 10*x^3*y^2 + 10*x^2*y^3 + 5*x*y^4 - 17*x^4 - 68*x^3*y - 102*x^2*y^2 - 68*x*y^3 - 38*y^4 +
	+ 65*x^3 + 195*x^2*y + 195*x*y^2 - 63*y^3 + 54*x^2 + 108*x*y - 235*y^2 - 171*x - 372*y - 108 # = 0
y^5 + 5*x^4*y + 10*x^3*y^2 + 10*x^2*y^3 + 5*x*y^4 + 21*x^4 - 68*x^3*y - 102*x^2*y^2 - 68*x*y^3 - 17*y^4 +
	- 330*x^3 + 195*x^2*y + 195*x*y^2 + 65*y^3 + 1252*x^2 + 108*x*y + 54*y^2 - 1551*x - 171*y + 398 # = 0


### Derivation
# [old] [explicit]
p1 = toPoly.Class3.pm(3)
print.p(p1, "x")

p2 = toPoly.Class3.pm(3, xn="y")
p2 = replace.pm(p2, data.frame(s0=1, coeff=-1), "s0")
p2 = replace.pm(p2, data.frame(s1=1, coeff=-1), "s1")
print.p(p2, "y")

p3 = toPoly.Class3.pm(3, s.id=c(1,2));
p3 = replace.pm(p3, data.frame(x=1:0, y=0:1, coeff=1), "x");
p3 = replace.pm(p3, data.frame(s2=1, coeff=2), "s2");
# simplify: s11 + s21 = 0
p3 = replace.pm(p3, 0, "s1");
print.p(p3, c("x", "y"))

x^3 + 3*x^2*y + 3*x*y^2 + y^3 + 2*s2*x^2 + 4*s2*x*y + 2*s2*y^2 - 8*s2^2*x - 8*s2^2*y - 8*s2^3 # = 0

pS1 = diff.pm(p3, p2)
pS1 = replace.pm(pS1, s[-4], paste0("s", 0:2));
print.p(pS1, c("x","y"))

pS2 = diff.pm(p3, p1)
pS2 = replace.pm(pS2, s[-4], paste0("s", 0:2));
print.p(pS2, c("x","y"))

