########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S5: Hetero-Symmetric
### Helper Functions
###   for Baseline System
###
### draft v.0.1a


### Derivation:
# - Helper functions:
#   were used to solve numerically the Baseline system;

# this file:
# source("Poly.System.S5.Ht.Formulas.Derivation.HS0.R")


#######################

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
test.S5Ht = function(x, R=c(0,0,0,0,0)) {
	R = R;
	test.f = function(x) {
		err = solve.S5HtMixed.Num(rbind(Re(x), Im(x)), R=R);
		err = err[1,] + 1i * err[2,];
	}
	round0(apply(x, 1, test.f));
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
	if(is.character(x0)) x0 = x0All[[x0]];
	plot.path(R, R0=R0, x0=x0, FUN=solve.S5HtMixed.Num, steps=steps, col=col, ...,
		subset=subset, p0.pch=p0.pch, debug=debug);
}

