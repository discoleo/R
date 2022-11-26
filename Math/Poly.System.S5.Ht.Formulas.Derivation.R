########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S5: Hetero-Symmetric
### Useful Formulas
###
### draft v.0.1f


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

#######################
#######################

### Base-System

### Numerical Solution:

### List of Initial values
source("Poly.System.S5.Ht.Formulas.Derivation.x0.R")

### Base-Function
solve.S5HtMixed.Num = function(x, R=c(0,1,0,0,1)) {
	x = matrix(x, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
	x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4]; x5 = x[5];
	s1 = x1 + x3; s2 = x2 + x4; S = s1 + s2 + x5;
	p1 = x1 * x3; p2 = x2 * x4; E5 = p1 * p2 * x5;
	ps = s1 * s2; sp = p1 + p2;
	# E2 = x1*(S - x1) + x2*(x3 + x4 + x5) + x3*(x4 + x5) + x4*x5;
	# E2 = sp + ps + x5*(S - x5);
	E3 = p1*s2 + p2*s1 + x5*(sp + ps);
	E4 = p1*p2 + x5*(p1*s2 + p2*s1);

	### E2:
	E11a = x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1;
	#
	y = c(S, E11a, E3, E4, E5) - R;
	y = rbind(Re(y), Im(y));
	return(y);
}
solve.path.S5 = function(R, R0, x0, steps=6, start.at=0, check=TRUE, verbose=FALSE) {
	path = expand.path(R0, R, steps=steps, start.at=start.at);
	if(is.character(x0)) x0 = x0All[[x0]];
	x.all = solve.path(solve.S5HtMixed.Num, x0, path=path, debug=verbose);
	if(check) {
		idDuplicates = which.perm.S5(x.all, verbose=verbose);
		if(ncol(idDuplicates) > 0) {
			cat("Warning: Duplicate roots!\n");
			print(idDuplicates);
		}
	}
	return(x.all);
}
polyS = function(R, x0, sol.rm=NULL, debug=FALSE, tol=5E-7) {
	if(is.character(x0)) x0 = x0All[[x0]];
	x.all = solve.all(solve.S5HtMixed.Num, x0, R=R, debug=debug);
	if( ! is.null(sol.rm)) x.all = x.all[ - sol.rm, ];
	if(nrow(x.all) != 7) {
		# stop("Wrong number of solutions!");
		cat("ERROR: Wrong number of solutions!\n");
		return(x.all);
	}
	p = poly.calc0(apply(x.all, 1, function(x) sum(x * x[c(3,4,5,1,2)])), tol=tol) * 27;
	print(p);
	return(x.all);
}
poly.calc.S5 = function(x, tol=5E-7) {
	if(is.character(x)) x = x0All[[x]];
	if(length(x) != 35) stop("Wrong roots!");
	p = poly.calc0(apply(x, 1, function(x) sum(x * x[c(3,4,5,1,2)])), tol=tol) * 27;
	return(p);
}
### Test solutions
test.sol = function(x, R=c(0,0,0,0,0)) {
	R = R;
	test.f = function(x) {
		err = solve.S5HtMixed.Num(rbind(Re(x), Im(x)), R=R);
		err = err[1,] + 1i * err[2,];
	}
	round0(apply(x.all, 1, test.f));
}
# Solve Eqs for coefficients
# - mini-solve for 2 unknown coefficients:
#   eval(EXPc) %*% B = ValCoeff, where
#   EXPc = "c(..., ...)" and
#   ValCoeff = b0[i] + FUN(R"i"), i = 1:2;
solve.coeff = function(R1, R2, b0, EXPc, FUN) {
	if(is.character(EXPc)) EXPc = parse(text=EXPc);
	coeff = function(R, b) {
		S = R[1]; E11a = R[2]; E3 = R[3]; E4 = R[4]; E5 = R[5];
		m = eval(EXPc, list(S=S, E11a=E11a, E3=E3, E4=E4, E5=E5));
		v = b - FUN(R);
		list(m=m, v=v);
	}
	l1 = coeff(R1, b0[1]);
	l2 = coeff(R2, b0[2]);
	r  = solve(rbind(l1$m, l2$m), c(l1$v, l2$v));
	return(r);
}

# Estimate Maximal Power of a Parameter
# npos = which variable in vector R;
# pow = which power in the polynomial;
# v = value used to vary parameter R;
max.pow.S = function(R, x0, pow, FUN, npos=1, v=2, R2=NULL, skip.path=FALSE, check=TRUE,
		steps=6, debug=FALSE) {
	if(is.character(x0)) x0 = x0All[x0];
	Rn = R; Rn[npos] = - Rn[npos];
	x01 = polyS(R,  x0=x0[[1]], debug=debug);
	x02 = polyS(Rn, x0=x0[[2]], debug=debug);
	#
	if(is.null(R2)) {
		R2  = R; R2[npos] = v;
		R2n = R2; R2n[npos] = - v;
	} else {
		v = R2[npos];
		R2n = R2; R2n[npos] = - v;
	}
	if(skip.path) {
		x03 = polyS(R2,  x0=x0[[3]], debug=debug);
		x04 = polyS(R2n, x0=x0[[4]], debug=debug);
	} else {
		path = expand.path(R, R2, steps=steps);
		x03  = solve.path(solve.S5HtMixed.Num, x0[[1]], path=path, debug=debug);
		idDuplicates = which.perm.S5(x03, verbose=FALSE);
		if(ncol(idDuplicates) != 0) {
			cat("Error in roots: step 3!\n");
			print(idDuplicates);
			if(check) return(x03);
		}
		#
		path = expand.path(Rn, R2n, steps=steps);
		x04  = solve.path(solve.S5HtMixed.Num, x0[[2]], path=path, debug=debug);
		idDuplicates = which.perm.S5(x04, verbose=FALSE);
		if(ncol(idDuplicates) != 0) {
			cat("Error in roots: step 4!\n");
			print(idDuplicates);
			if(check) return(x04);
		}
	}
	#
	pow = pow + 1;
	p1 = poly.calc.S5(x01)[[pow]];
	p2 = poly.calc.S5(x02)[[pow]];
	p3 = poly.calc.S5(x03);
	p4 = poly.calc.S5(x04);
	if( ! skip.path) { print(p3); print(p4); }
	p3 = p3[[pow]]; p4 = p4[[pow]];
	b0 = c(FUN(R), FUN(Rn), FUN(R2), FUN(R2n));
	powPlus  = (p3 + p4 - b0[3] - b0[4]) / (p1 + p2 - b0[1] - b0[2]);
	powMinus = (p3 - p4 - b0[3] + b0[4]) / (p1 - p2 - b0[1] + b0[2]);
	div = log(abs(v / R[npos]));
	r = c(powPlus, powMinus,
		log(abs(powPlus)) / div, log(abs(powMinus)) / div);
	return(r);
}

plot.path.S5 = function(R, R0, x0, steps=20, col=seq(length(R0)), ...,
		subset=NULL, p0.pch=5, debug=FALSE) {
	plot.path(R, R0=R0, x0=x0, FUN=solve.S5HtMixed.Num, steps=steps, col=col, ...,
		subset=subset, p0.pch=p0.pch, debug=debug);
}
plot.path = function(R, R0, x0, FUN, steps=20, col=seq(length(R0)), ...,
		subset=NULL, p0.pch=5, debug=FALSE) {
	if(is.character(x0)) x0 = x0All[[x0]];
	isMatrix = is.matrix(x0);
	len = if(isMatrix) nrow(x0) else 1;
	# Path:
	path = expand.path(R0, R, steps=steps);
	x = array(0, c(length(R0), 0));
	for(id in seq(steps)) {
		x0 = solve.all(FUN, x0=x0, R=path[[id]], ..., debug=debug);
		x  = cbind(x, t(x0));
	}
	nr = nrow(x);
	isSubset = ! is.null(subset);
	if(nr > 1) {
		xlim = c(min(Re(x)), max(Re(x)));
		ylim = c(min(Im(x)), max(Im(x)));
		plot(Re(x[1,]), Im(x[1,]), col=col[1], xlim=xlim, ylim=ylim);
		idAll = if(isSubset) subset else seq(2, nr);
		for(id in idAll) {
			points(Re(x[id,]), Im(x[id,]), col=col[id]);
		}
		# Initial points
		if(p0.pch > 0) {
			p0 = function(x, col) {
				points(jitter(Re(x)), Im(x), col=col, pch=p0.pch, cex=1.5, lwd=1.5);
			}
			if(isSubset) {
				p0(x[- subset, seq(len)], col="blue");
				p0(x[  subset, seq(len)], col="red");
			} else {
				p0(x[, seq(len)], col="red");
			}
		}
	}
	invisible(x);
}

###################
###################

# Note:
# - the double permutation (x1, x3) & (x4, x5) like c(x3, x2, x1, x5, x4),
#   and equivalent permutations (e.g. (x1, x2) & (x3, x5)),
#   are also a solution;

### Set 1:
x0 = c(-0.70449+0.64i, 0.8913, -0.70449-0.64i, 0.2589 + 1.08i, 0.2589 - 1.08i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x)


### Set 2:
x0 = c(-1.70449-0.64i, -0.8913, -1.70449+0.64i, -0.2589 + 1.08i, -0.2589 - 1.08i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 3:
x0 = c(0.05-0.94i, -2.8913, 0.05+0.04i, -0.2589 + 1.08i, -0.2589 - 1.08i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 4:
x0 = c(-1.274+0.729i, 1.23-0.47i, 0.58-0.67i, 0.178+0.725i, -0.71 - 0.31i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 5:
x0 = c(-1.274-0.729i, 1.23+0.47i, 0.58+0.67i, 0.178-0.725i, -0.71 + 0.31i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 6:
x0 = c(-0.8108+1.5014i, 0.7763-1.6039i, 0.6605-0.1778i, -0.5008-0.4342i, -0.1252+0.7145i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 7:
x0 = c(-0.8108-1.5014i, 0.7763+1.6039i, 0.6605+0.1778i, -0.5008+0.4342i, -0.1252-0.7145i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)
x.all = matrix(x.all, nc=5, byrow=T)

round0(poly.calc(x.all)) * 27
poly.calc(apply(x.all, 1, function(x) sum(x * x[c(3,4,5,1,2)]))) * 27

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

max.pow.S(c(1,0,1,0,2), c("E3V101", "E3Vn101"), pow=4, FUN=f4, v=3)


max.pow.S(c("E3V131", "E3Vn131", "E3V531", "E3Vn531"), pow=4, FUN=f4, skip.path=T, R=c(1,2.8,1,0,2), R2=c(5,3,1,0,2))


solve.coeff(c(1,2.8,1,0,2), c(5,3,1,0,2), c(- 217.8506 - 1100.902, - 19610.25 + 65015.25) / 2,
"c(...)", function(R) { Rn = R; Rn[1] = - Rn[1]; (f4(R) - f4(Rn))/2; })


27*(E11a^7 + E11b^7)*E5^2 +
	# x^6
	- 27*(E11a^6 + E11b^6)*E5^2*S^2 - (E11a*E11b)^6 +
	+ 81*E11a*E11b*(E11a^5 + E11b^5)*E5^2 + 9*(E11a*E11b)^3*(E11a^3 + E11b^3)*E5*S +
	# x^5:
	+ (E11a^5 + E11b^5)*(9*E5^2*S^4 + 4*E3^3*E5 - 90*E3*E5^2*S) +
	+ (E11a*E11b)*(E11a^4 + E11b^4)*(9*E3^2*E5*S) +
	+ (E11a*E11b)^2*(E11a^3 + E11b^3)*(27*E5^2 - 18*E3*E5*S^2) +
	- (E11a*E11b)^3*(E11a^2 + E11b^2)*(2*E5*S^3 + 15*E3*E5 + 0*E3^2*S^2) +
	- (E11a*E11b)^4*(E11a + E11b)*(3*E5*S + 2*E3^2) +
	+ 4*(E11a*E11b)^5*E3*S +
	# x^4: Note: 1 term from x^5 contributes as well;
	- (E11a^4 + E11b^4)*E5^2*S^6 + 18*(E11a*E11b)^2*(E11a^2 + E11b^2)*E5^2*S^2 +
	+ 4*(E11a*E11b)^4*E5*S^3 - 10*(E11a*E11b)^3*(E11a + E11b)*E5^2 +
	- 21*(E11a*E11b)*(E11a^3 + E11b^3)*E5^2*S^4 + 150*(E11a^4 + E11b^4)*E3^2*E5^2 +
	- 11*(E11a*E11b)*(E11a^3 + E11b^3)*E3^3*E5 - (E11a*E11b)^2*(E11a^2 + E11b^2)*E3^4 +
	+ 20*(E11a*E11b)^4*E3*E5 +
	+ (E11a^4 + E11b^4)*(48*E3*E5^2*S^3 - 10*E3^3*E5*S^2) +
	# x^3:
	+ (E11a^3 + E11b^3)*(6*E3^4*E5*S - 6*E3*E5^2*S^5 +
		+ 2*E3^3*E5*S^4 - 50*E5^3*S^3 + 5*E3^2*E5^2*S^2) +
	- 375*(E11a*E11b)*(E11a^2 + E11b^2)*E5^3*S + 4*(E11a*E11b)*(E11a^2 + E11b^2)*E5^2*S^6 +
	+ 12*(E11a*E11b)^2*(E11a + E11b)*E5^2*S^4 - 68*(E11a*E11b)^3*E5^2*S^2 +
	+ 275*(E11a*E11b)*(E11a^2 + E11b^2)*E3^2*E5^2 - 4*(E11a*E11b)^3*E3^4 +
	- 3*(E11a*E11b)^2*(E11a + E11b)*E3^3*E5 +
	# x^2:
	+ (E11a^2 + E11b^2)*(14*E5^3*S^5 + 500*E3*E5^3*S^2 - 21*E3^2*E5^2*S^4 - E3^6*S^2 +
		- 200*E3^3*E5^2*S + 8*E3^4*E5*S^3) +
	+ 200*(E11a*E11b)*(E11a + E11b)*E5^3*S^3 +
	- (E11a*E11b)^2*(6*E5^2*S^6 + 750*E5^3*S) +
	- 5^5*E5^4*(E11a^2 + 3*E11a*E11b + E11b^2) +
	+ 12*(E11a^2 + E11b^2)*E3^5*E5 + 5^4*(E11a*E11b)*(E11a + E11b)*E3*E5^3 +
	- 2*(E11a*E11b)*(E11a + E11b)*E3^6 - 375*(E11a*E11b)^2*E3^2*E5^2 +
	# x^1:
	+ (E11a + E11b)*(5^5*E5^4*S^2 + 125*E3^4*E5^2 + 2*E3^7*S - 1250*E3^2*E5^3*S +
		+ 110*E3^3*E5^2*S^3 - 26*E3^5*E5*S^2 - 150*E3*E5^3*S^4) +
	- E11a*E11b*(28*E5^3*S^5 + 34*E3^5*E5 + ...*E3^4*E5*S^3 + ...*E3^3*E5^2*S) +
	# B0:
	- 5^4*E5^4*S^4 - E3^8 - 150*E3^4*E5^2*S^2 + 20*E3^6*E5*S + 500*E3^2*E5^3*S^3 # = 0


### Coefficients:
f0 = function(R) {
	S = R[1]; E11a = R[2]; E3 = R[3]; E4 = R[4]; E5 = R[5];
	cc = 27*E5^2*E11a^7 - 27*S^2*E5^2*E11a^6 +
		+ E11a^5*(9*E5^2*S^4 + 4*E3^3*E5 - 90*E3*E5^2*S) +
		- E11a^4*(E5^2*S^6 - 150*E3^2*E5^2 - 48*E3*E5^2*S^3 + 10*E3^3*E5*S^2) +
		+ E11a^3*(- 6*E3*E5^2*S^5 - 50*E5^3*S^3 + 6*E3^4*E5*S + 2*E3^3*E5*S^4 + 5*E3^2*E5^2*S^2) +
		+ E11a^2*(14*E5^3*S^5 + 500*E3*E5^3*S^2 - 21*E3^2*E5^2*S^4 - E3^6*S^2 +
			- 200*E3^3*E5^2*S + 8*E3^4*E5*S^3 - 3125*E5^4 + 12*E3^5*E5) +
		+ E11a*(3125*E5^4*S^2 + 125*E3^4*E5^2 + 2*E3^7*S - 1250*E3^2*E5^3*S + 110*E3^3*E5^2*S^3 +
			- 26*E3^5*E5*S^2 - 150*E3*E5^3*S^4) +
		- 625*S^4*E5^4 - E3^8 - 150*E3^4*E5^2*S^2 + 20*E3^6*E5*S + 500*E3^2*E5^3*S^3;
	cc / E5^2;
}
f1 = function(R) {
	S = R[1]; E11a = R[2]; E3 = R[3]; E4 = R[4]; E5 = R[5];
	cc = 3125*E5^4*S^2 + 125*E3^4*E5^2 + 2*E3^7*S - 1250*E3^2*E5^3*S + 110*E3^3*E5^2*S^3 +
		- 26*E3^5*E5*S^2 - 150*E3*E5^3*S^4 +
		+ 81*E5^2*E11a^6 - E11a^4*(21*S^4*E5^2 + 11*E3^3*E5) +
		+ E11a^3*(4*E5^2*S^6 - 375*E5^3*S + 275*E3^2*E5^2) +
		+ E11a^2*(200*E5^3*S^3 - 2*E3^6 + 5^4*E3*E5^3) +
		- E11a*(9375*E5^4 + 28*E5^3*S^5 + 34*E3^5*E5);
	cc / E5^2;
}
f2 = function(R) {
	S = R[1]; E11a = R[2]; E3 = R[3]; E4 = R[4]; E5 = R[5];
	14*S^5*E5 - 3125*E5^2 + 12*S^4*E11a^3 + 18*S^2*E11a^4 + 27*E11a^5 +
	- 6*E11a^2*S^6 - 750*E11a^2*E5*S + 200*E11a*E5*S^3 + 12*E3^5/E5 +
	- E11a^4*E3^4/E5^2 - 3*E11a^3*E3^3/E5 + 5^4*E11a*E3*E5 +
	- 2*E11a*E3^6/E5^2 - 375*E11a^2*E3^2 + 500*E3*E5*S^2 - 21*E3^2*S^4 - E3^6*S^2/E5^2 +
	- 200*E3^3*S + 8*E3^4*S^3/E5;
}
f3 = function(R) {
	S = R[1]; E11a = R[2]; E3 = R[3]; E4 = R[4]; E5 = R[5];
	- 50*S^3*E5 - 10*E11a^4 - 2*S^3*E11a^5/E5 + 9*S*E11a^6/E5 +
	+ 12*E11a^2*S^4 - 375*E11a*E5*S - 68*E11a^3*S^2 + 4*E11a*S^6 +
	- 15*E11a^5*E3/E5 - 3*E11a^2*E3^3/E5 +
	+ 275*E11a*E3^2 - 4*E11a^3*E3^4/E5^2 + 6*E3^4*S/E5 - 6*E3*S^5 +
	+ 2*E3^3*S^4/E5 + 5*E3^2*S^2;
}
f4 = function(R) {
	S = R[1]; E11a = R[2]; E3 = R[3]; E4 = R[4]; E5 = R[5];
	150*E3^2 - S^6 + 18*E11a^2*S^2 + 4*E11a^4*S^3/E5 - 10*E11a^3 - 21*E11a*S^4 - 3*E11a^5*S/E5 +
	- 2*E11a^5*E3^2/E5^2 - E11a^2*E3^4/E5^2 - 11*E11a*E3^3/E5 + 20*E11a^4*E3/E5 +
	+ 48*E3*S^3 - 10*E3^3*S^2/E5;
}
f5 = function(R) {
	S = R[1]; E11a = R[2]; E3 = R[3]; E4 = R[4]; E5 = R[5];
	9*S^4 + 4*E3^3/E5 - 90*E3*S +
	+ 9*E11a*E3^2*S/E5 + E11a^2*(27 - 18*E3*S^2/E5) +
	- E11a^3*(2*S^3/E5 + 15*E3/E5) - E11a^4*(3*S/E5 + 2*E3^2/E5^2) +
	+ 4*E11a^5*E3*S/E5^2;
}
f6 = function(R) {
	S = R[1]; E11a = R[2]; E3 = R[3]; E4 = R[4]; E5 = R[5];
	- 27*S^2 - E11a^6/E5^2 + 81*E11a + 9*E11a^3*S/E5;
}


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

