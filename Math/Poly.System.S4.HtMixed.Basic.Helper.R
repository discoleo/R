########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Hetero-Symmetric S4: Mixed
### Basic Types: Helper Functions
###
### draft v.0.1a


### this file
# source("Poly.System.S4.HtMixed.Basic.Helper.R")


########################

### Helper Functions

source("Polynomials.Helper.R")


### Other

test.S4HtMixed = function(sol, n=2, nE2 = 1, R = NULL) {
	x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];
	err1 = x1^n + x2^n + x3^n + x4^n;
	# Ht
	ht = cbind(x1*x2, x2*x3, x3*x4, x4*x1);
	ht = if(nE2 == 1) ht else ht^nE2;
	err2 = apply(ht, 1, sum);
	err3 = x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4;
	err4 = x1*x2*x3*x4;
	err = rbind(err1, err2, err3, err4);
	if( ! is.null(R)) {
		for(id in 1:4) err[id,] = err[id,] - R[id];
	}
	err = round0(err);
	return(err);
}

test.S4HtMixed.En3 = function(sol, R=NULL, n=2, nE=c(1,2,1)) {
	# Ht
	ht.f = function(x) {
		sum(x^nE[1] * (x^nE[2])[c(2,3,4,1)] * (x^nE[3])[c(3,4,1,2)]);
	}
	err2 = apply(sol, 1, ht.f);
	#
	x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];
	err1 = x1^n + x2^n + x3^n + x4^n;
	err3 = x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4;
	err4 = x1*x2*x3*x4;
	err = rbind(err1, err2, err3, err4);
	if( ! is.null(R)) {
		for(id in 1:4) err[id,] = err[id,] - R[id];
	}
	err = round0(err);
	return(err);
}

#######################

### Solvers: Base Cases

### E22a:
solve.S4HtM.Ord2Base = function(R, E2, sort=TRUE, all.sol=TRUE) {
	S = R[1]; E22a = R[2]; E3 = R[3]; E4 = R[4];
	#
	len = length(E2);
	x1 = sapply(seq(len), function(id) roots(c(1, -S, E2[id], -E3, E4)));
	x1 = as.vector(x1);
	E2 = rep(E2, each=4);
	len = length(x1);
	# fully robust:
	x13T1 = x1^2*(E2^2 - E22a - 2*S*E3 + 2*E4);
	x3sq = x1^6*(x1^2 + 2*E2 - S^2) + x1^4*E22a - E4^2;
	div  = x1^4*(2*x1^2 + 2*E2 - S^2) + x13T1;
	x3sq = - x3sq / div;
	x3 = (x3sq*(S - x1) + E4/x1) / (x3sq + E2 - x1*(S - x1));
	x3 = as.vector(x3);
	# x2, x4:
	xs = S - x1 - x3; x24 = E4 / (x1*x3);
	# Note: both sqrt() values are valid; (option: all.sol)
	xd = sqrt(xs^2 - 4*x24 + 0i);
	x2 = (xs + xd)/2; x4 = (xs - xd)/2;
	sol = cbind(x1=x1, x2=x2, x3=x3, x4=x4);
	if(all.sol) sol = rbind(sol, sol[, c(1,4,3,2)]); # all roots
	#
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}

### E121a:
solve.S4HtM.E121Base = function(R, E2, sort=TRUE, all.sol=FALSE) {
	len = length(E2);
	# robust:
	S = R[1]; E121a = R[2]; E3 = R[3]; E4 = R[4];
	x1 = sapply(seq(len), function(id) roots(c(1, -S, E2[id], -E3, E4)));
	E2 = rep(E2, each=4); x1 = as.vector(x1);
	len = length(x1);
	#
	x3 = sapply(seq(len), function(id) solve.x3.S4HtM.E121P1(R, E2[id], x1[id]));
	#
	xs = S - x1 - x3; x24 = E4 / (x1*x3);
	xd = sqrt(xs^2 - 4*x24 + 0i);
	x2 = (xs + xd)/2;
	x4 = (xs - xd)/2;
	sol = cbind(x1=x1, x2=x2, x3=x3, x4=x4);
	if(all.sol) sol = rbind(sol, sol[, c(1,4,3,2)]);
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}
solve.x3.S4HtM.E121P1 = function(R, E2, x1) {
	S = R[1]; E121 = R[2]; E3 = R[3]; R4 = R[4]; # E4
	E2d = E2 - x1*(S-x1);
	x0 = x1^2*R4^2*S^2 - 4*x1^3*R4^2*S - x1^3*R4*E121*S + 4*x1^4*R4^2 + 4*R4^3 +
		+ 2*x1^4*R4*E121 - 4*x1*R4^2*E3 + x1^2*R4*E3^2;
	div = - x1*R4^2*S^2 + x1^2*R4*E3*S^2 + 4*x1^2*R4^2*S - 4*x1^3*R4*E3*S - x1^3*E121*E3*S - 4*x1^3*R4^2 +
		+ 4*x1*E2d*R4^2 + 2*x1^3*R4*E121 + x1^3*E121^2 + 4*x1^4*R4*E3 - 4*x1^2*E2d*R4*E3 + x1^4*E121*E3 + x1^3*E2d*E3^2;
	# TODO
	if(round0(div) == 0) warning("Div by 0!");
	return(x0 / div);
}

### Compute E2

e2.f = function(x) {
	e2.f0 = function(x) x[1]*sum(x, -x[1]) + x[2]*(x[3]+x[4]) + x[3]*x[4];
	e2 = if(is.matrix(x)) apply(x, 1, e2.f0) else e2.f0(x);
	sort.sol(matrix(e2, ncol=1), useRe=TRUE);
}
e3.f = function(x) {
	e3.f0 = function(x) {
		x34 = x[3]*x[4];
		x[1]*(x[2]*(x[3] + x[4]) + x34) + x[2]*x34;
	}
	e3 = if(is.matrix(x)) apply(x, 1, e3.f0) else e3.f0(x);
	sort.sol(matrix(e3, ncol=1), useRe=TRUE);
}
e2a.f = function(x) {
	e2.f0 = function(x) sum(x * x[c(2,3,4,1)]);
	e2 = if(is.matrix(x)) apply(x, 1, e2.f0) else e2.f0(x);
	sort.sol(matrix(e2, ncol=1), useRe=TRUE);
}
e3a.f = function(x, pow=c(1,2,1)) {
	e2.f0 = function(x) sum(x^pow[1] * x[c(2,3,4,1)]^pow[2] * x[c(3,4,1,2)]^pow[3]);
	e2 = if(is.matrix(x)) apply(x, 1, e2.f0) else e2.f0(x);
	sort.sol(matrix(e2, ncol=1), useRe=TRUE);
}

### Solve Coefficients:
# - hack the formulas;

which.sq = function(x, sq=2, iter=1000, digits=6, pow=2) {
	if(is.na(x)) return(NA);
	if(round(x) == round(x, digits)) return(0);
	i = seq(iter);
	sq = if(pow == 2) sqrt(sq) else rootn(sq, n=pow);
	d = round(i * sq - x, digits);
	id = which(d == round(d));
	if(length(id) > 0) return(id);
	d = round(i * sq + x, digits);
	id = which(d == round(d));
	if(length(id) == 0) return(NA);
	return(- id);
}

# - hack higher powers;
whichHasPower = function(R, id=2, type=1, FUN=NULL, print=FALSE, digits=5, iter=1000) {
	if(length(R) == 1) {len = R; R0 = rep(1, len); }
	else {len = length(R); R0 = R; }
	vals = c(2,3,5); vsqrt = sqrt(vals);
	if(type == 1) {
		vsqrt = c(vsqrt, - vsqrt); vals = c(vals, vals);
	} else if(type == 2) {
		vsqrt = c(vsqrt, 2*vsqrt); vals = c(vals, vals);
	} else if(type == 3) {
		vsqrt = c(vsqrt, - vsqrt) + 1; vals = c(vals, vals);
	} else if(type == 4) {
		vsqrt = c(vsqrt, 2*vsqrt) - 1; vals = c(vals, vals);
	}
	VLEN = length(vals);
	m = array(NA, c(len, VLEN));
	f0 = if( ! is.null(FUN)) {
		function(vid, nr) {
			R = R0;
			R[nr] = vsqrt[vid];
			which.sq(FUN(R), sq=vals[vid], pow=2, digits=digits, iter=iter)
		}
	} else function(vid, nr) {
		R = R0;
		R[nr] = vsqrt[vid];
		which.coeff(R, sq=vals[vid], id=id, pow=2, digits=digits, iter=iter, print=print)
	}
	for(nr in seq(len)) {
		tmp = sapply(seq(VLEN), f0, nr);
		m[nr, ] = tmp;
	}
	return(m);
}
which.coeff = function(R, FUN, sq=2, id=3, pow=2, DIFF=NULL, print=TRUE, digits=6, iter=1000) {
	p = FUN(R);
	if(print) print(p);
	x = p[id]; # which coefficient;
	if( ! is.null(DIFF)) {
		x = x - DIFF(R);
	}
	return(which.sq(x, sq=sq, pow=pow, digits=digits, iter=iter))
}
which.coeff.gen = function(FUN) {
	f = function(R, sq=2, id=3, pow=2, DIFF=NULL, print=TRUE, digits=6, iter=1000, debug=FALSE) {
		p = FUN(R, debug=debug);
		if(print) print(p);
		x = p[id]; # which coefficient;
		if( ! is.null(DIFF)) {
			x = x - DIFF(R);
		}
		return(which.sq(x, sq=sq, pow=pow, digits=digits, iter=iter))
	}
}

### Formulas

### Simple: E2a
polyE2Ord1 = function() {
	p = toPoly.pm("E2a*E2^2 - (S*E3 + 2*E2a^2)*E2 +
		+ S^2*E4 + S*E2a*E3 + E2a^3 - 4*E2a*E4 + E3^2");
	return(p);
}
polyE2Ord2 = function() {
	p = toPoly.pm("(4*E4 - E22a)*E2^4 - 2*(E3^2 + S^2*E4)*E2^3 +
	+ (4*(S*E3 - 2*E4)*E22a - 8*S*E3*E4 + S^2*E3^2 + 2*E22a^2)*E2^2 +
	+ (4*S^3*E3*E4 + 4*S*E3^3 + 2*(S^2*E4 + E3^2)*E22a)*E2 +
	- S^4*E4^2 - 2*S^3*E3^3 - E22a^3 - E3^4 + 4*E22a^2*E4 + 2*S^2*E3^2*E4 +
		+ 8*S*E22a*E3*E4 - 4*S*E22a^2*E3 - 5*S^2*E22a*E3^2");
	return(p);
}
### Type: E121a
### Order 1:
polyE2_E121P1 = function() {
	p = toPoly.pm("4*(4*E4 + E121a)*E4*E2^2 - (E121a + 8*E4)*(S^2*E4 + E3^2)*E2 +
		+ S^4*E4^2 + S*E121a^2*E3 - E121a^3 + E3^4 - 4*E121a^2*E4 + 2*S^2*E3^2*E4");
	return(p);
}
