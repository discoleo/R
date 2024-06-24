########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Multi-Variable Polynomials


### fast load:
# source("Polynomials.Helper.R")

### Note:
# - Substitutions moved to new file:
#   Polynomials.Helper.Substitution.R;
# - Division moved to new file:
#   Polynomials.Helper.Div.R;
# - mpfr-Functions moved to new file:
#   Polynomials.Helper.mpfr.R;
# - Factorization moved to new file:
#   Polynomials.Helper.Factorize.R;
# - Format & Print moved to new file:
#   Polynomials.Helper.Format.R;
# - Parser moved to new file:
#   Polynomials.Helper.Parser.R;
# - Derivation moved to new file:
#   Polynomials.Helper.D.R;
# - Generic Solvers moved to new file:
#   Polynomials.Helper.Solvers.R;


#######################

library(polynom)
library(pracma)

### Helper Functions

# required:
source("Polynomials.Helper.Parser.R")
source("Polynomials.Helper.Format.R")
source("Polynomials.Helper.Substitution.R")
source("Polynomials.Helper.Div.R")
source("Polynomials.Helper.D.R")
source("Polynomials.Helper.Factorize.R")

### Note:
# D-Functions: required in factorize.p();


### Big Numbers
startBigNumbers = function() {
	# is loaded automatically
	# library(gmp)
	source("Polynomials.Helper.BigNumbers.R");
}


### Basic Algebra

### Round to 0
round0 = function(x, ...) UseMethod("round0")
round0.default = function(m, tol=1E-7) {
	m[abs(Re(m)) < tol & abs(Im(m)) < tol] = 0
	isZero = (Re(m) != 0) & (abs(Re(m)) < tol)
	if(any(isZero)) {
		m[isZero] = complex(re=0, im=Im(m[isZero]))
	}
	isZero = (Im(m) != 0) & (abs(Im(m)) < tol)
	if(any(isZero)) {
		m[isZero] = Re(m[isZero])
	}
	return(m)
}
round0.p = function(p, tol=1E-7) round0.polynomial(p, tol=tol)
round0.polynomial = function(p, tol=1E-7) {
	p = round0(as.vector(p), tol=tol)
	class(p) = "polynomial"
	return(p)
}
round0.pm = function(p, tol=1E-7) {
	p$coeff = round0(p$coeff, tol=tol);
	return(p);
}
### Round:
round.pm = function(p, digits=4) {
	p$coeff = round(p$coeff, digits=digits);
	return(p);
}

### Root
rootn = function(r, n) {
	if(n %% 2 == 0) {
		if(all(Im(r) == 0) && all(Re(r) >= 0)) return(Re(r)^(1/n));
		return((r + 0i)^(1/n));
	}
	ifelse( (Im(r) == 0 & Re(r) < 0), - (-r)^(1/n), r^(1/n) )
}
unity = function(n=3, all=TRUE) {
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	if(all) {
		m = m^(0:(n-1))
	}
	return(m)
}
# Specific Roots:
roots.f = function(K, s, n=length(s)) {
	# roots for Class 1 polynomials;
	# s = includes s0;
	m = unity(n=n, all=T)
	k = rootn(K, n);
	order1 = n - 1;
	r = sapply(seq(n), function(id) sum(s * (k*m[id])^(0:order1)))
	return(r)
}
roots.cl2.f = function(s, n = length(s)) {
	# roots for basic Class 2 polynomials;
	m = unity(n+1, all=T)[-1]; # exclude 1;
	r = sapply(seq(n), function(id) sum(s * m[id]^(0:n)))
	r = round0(r)
}
### Compute the roots:
roots.pm = function(p, ..., xn="x", sort=TRUE, split=TRUE, tol=1E-8) {
	r = roots(evalCoeff(p, xn=xn, ...));
	if(sort) r = sort(r);
	# move real roots at the end;
	if(split) r = split.roots.complex(r, as.list=FALSE, tol=tol)
	return(r);
}
### Split into complex and pure real values:
# - move real roots at the end of the vector;
# - useful for automatic access to conjugated roots;
split.roots.complex = function(x, as.list=FALSE, tol=1E-8, f=NULL, ...) {
	if( ! is.null(f)) {
		return(split.default(x, f=f, ...));
	}
	len = length(x);
	# move real values at the end of the vector;
	isRe = (round0(Im(x), tol=tol) == 0);
	x = split.default(x, f = isRe, ...);
	names(x) = NULL;
	if( ! as.list) x = unlist(x, recursive=FALSE);
	lenRe = sum(isRe);
	attr(x, "dim.z") = c(len - lenRe, lenRe);
	return(x);
}
multiplicity.pm = function(p, val, tol=1E-10, doGCD=TRUE) {
	if( ! is.pm(p)) stop("Not a polynomial!");
	nc = ncol(p);
	if(nc > 2) stop("p must be a Uni-variate polynomial!");
	idx = which( ! names(p) %in% "coeff");
	if(length(idx) != 1) stop("Not a polynomial!");
	xn = names(p)[idx];
	len = maxPow.pm(p, xn=xn);
	if(len == 0) return(0);
	#
	if(round0(eval.pm(p, val), tol=tol) != 0) return(0);
	dp = p;
	for(id in seq(len)) {
		dp = dp.pm(dp, xn=xn);
		if(doGCD) {
			xgcd = gcd.vpm(dp);
			if(xgcd > 1) dp$coeff = dp$coeff / xgcd;
		}
		if(round0(eval.pm(dp, val), tol=tol) != 0) return(id);
	}
	return(len);
}

### Sort
sort.sol = function(sol, useRe=TRUE, mod.first=TRUE, ncol=1, digits=5) {
	if(useRe) {
		if(length(digits) == 1) digits = rep(digits, 3);
		mod = round(abs(sol[,ncol]), digits[1]);
		re = round(Re(sol[,ncol]), digits[2]);
		im = round(Im(sol[,ncol]), digits[3]);
		if(mod.first) {
			id = order(mod, re, im);
		} else {
			id = order(re, mod, im);
		}
	} else {
		id = order(abs(round(sol[,ncol], digits)) );
	}
	return(sol[id,, drop=FALSE]);
}
isConj.f = function(x, y, tol=1E-3) {
	isConj = (abs(Re(x) - Re(y)) < tol) & (abs(Im(x) + Im(y)) < tol);
	return(isConj);
}

as.conj = function(x, nrow=0) {
	if(nrow > 0) {
		x[nrow,] = abs(x[nrow,])^2 / x[nrow,];
		return(x);
	}
	abs(x)^2 / x;
}

###################
### Polynomials ###
###################

# Is polynomial:
is.pm = function(p, strict=FALSE) {
	if(strict) {
		inherits(p, "pm");
	} else inherits(p, "data.frame");
}
# isNonZero:
isNZ.pm = function(p) return(is.NZ.pm(p));
is.NZ.pm = function(p) {
	is.data.frame(p) && (nrow(p) > 0);
}
is.zero.pm = function(p) {
	if(is.numeric(p)) return(p == 0);
	if(is.complex(p)) return(p == 0);
	if(nrow(p) == 0) return(TRUE);
	isZero = all(p$coeff == 0);
	return(isZero);
}

### Simple Multiplication
# p1, p2 = simple vectors of coefficients;
mult.p = function(p1, p2) {
	p.m = outer(p1, p2)
    p = as.vector(tapply(p.m, row(p.m) + col(p.m), sum))
	return(p)
}

# Replace "poly.calc" & "polynomial" constructor in package polynom;
# - only basic construction:
poly.calc0 = function(x, toReal=TRUE, tol=1E-8, digits=0, as.pm=TRUE, warn=TRUE) {
	p = 1;
	for (xi in x) {
		p = c(0, p) - c(xi * p, 0);
	}
	return(polynomial(p, tol=tol, digits=digits, toReal=toReal, as.pm=as.pm, warn=warn));
}
polynomial = function(coef = c(0, 1), tol=1E-8, digits=0, toReal=TRUE,
		simplify = TRUE, as.pm=TRUE, warn=TRUE) {
	b   = round0(coef, tol=tol);
	isC = Im(b) != 0;
	if(warn && any(isC)) {
		Res = sum(abs(Im(b)));
		warning("Some coefficients are complex:\n  Residual = ", Res, "!");
	}
	if(toReal) b = Re(b);
	if(digits > 0) b = round(b, digits);
	last = length(b);
	if(last > 0) {
		for(npos in seq(last, 1)) {
			if(b[npos] != 0) break;
		}
		b = b[seq(1, npos)];
	}
	if(as.pm) {
		return(as.pm.polynomial(b, tol=0, simplify=simplify));
	}
	structure(b, class = "polynomial");
}

##################

### Multi-variable Polynomials:
# p = multi-variable polynomial;
#  => data.frame with exponents of variables;
#  => coeff = column with the coefficients;
# Note:
# - initial idea was to allow also basic lists,
#   but most functions would work only with data frames!

### Workhorse for SUM & Multiplication:
# - a hacked version is needed for Big Numbers (gmp, mpfr);
aggregate0.pm = function(p) {
	p.r = aggregate(coeff~., p, sum);
	return(p.r);
}

### Utils
dimnames.pm = function(p) {
	nms = names(p);
	# excludes the Coefficients:
	id = match("coeff", nms);
	return(nms[ - id]);
}

### names.pm():
# - interferes with 'names<-';
# - would require an explicit setter as well:
#   possible using explicit attr(p, "names") = ...;
#
# names.pm = function(p) {
#	nms = NextMethod();
#	# excludes the Coefficients:
#	id = match("coeff", nms);
#	return(nms[ - id]);
# }

### Leading Monomials
top.pm = function(p, xn="x", exclude=FALSE) {
	# exclude = exclude variable xn;
	if(nrow(p) == 0) return(data.frame(coeff=numeric(0)));
	# check names:
	idn = if(is.numeric(xn)) xn else match(xn, names(p));
	if(all(is.na(idn))) return(data.frame(coeff=numeric(0)));
	idn = idn[ ! is.na(idn)];
	#
	if(length(idn) == 1) {
		pow.max = max(p[, idn]);
		p = if(exclude) { p[p[,idn] == pow.max, -idn, drop=FALSE]; }
			else p[p[,idn] == pow.max, , drop=FALSE];
		return(p);
	} else {
		stop("Not yet implemented!")
	}
}
### Free Term
B0.pm = function(p, xn="x", warn=TRUE) {
	# Global B0:
	if(is.null(xn)) {	
		idc = match("coeff", names(p));
		px  = p[, - idc, drop=FALSE]; class(px) = "data.frame";
		isB0 = sapply(seq(nrow(px)), function(nr) all(px[nr,] == 0));
		return(p[isB0, , drop=FALSE]);
	}
	#
	idx = match(xn, names(p));
	if(is.na(idx)) {
		if(warn) warning("Variable not found!");
		return(p);
	}
	b0 = p[p[, xn] == 0, , drop=FALSE];
	return(b0);
}

maxPow.pm = function(p, xn) {
	max(p[, xn]);
}

### Basic Operations

### Multiplication
# TODO:
# - check consequences!
Ops.pm = function(e1, e2) {
	r = switch(.Generic,
		'-' = { diff.pm(e1, e2); },
		'*' = { mult.pm(e1, e2); },
		'%%' = {
			if(e2 < 2) stop("Modulus must be >= 2!");
			e1$coeff = e1$coeff %% e2;
			reduce.pm(e1);
		}
	)
	return(r);
}

prod.pm = function(...) {
	p = list(...);
	len = length(p);
	if(len < 1) return(data.frame()); # TODO
	isList = function(x) is.list(x) && ! is.pm(x);
	r = if(isList(p[[1]])) mult.lpm(p[[1]]) else p[[1]];
	if(len == 1) return(r);
	# len > 1
	for(id in seq(2, length(p))) {
		p0 = p[[id]];
		p0 = if(isList(p0)) mult.lpm(p0) else p0;
		r = mult.pm(r, p0);
	}
	return(r);
}
mult.all.pm = function(p) return(mult.lpm(p));
mult.lpm = function(p) {
	if( ! is.list(p) || is.pm(p, strict=TRUE))
		stop("p must be a list of polynomials!");
	len = length(p);
	if(len == 1) return(p[[1]]);
	pR = p[[1]];
	isNum = is.numeric(pR) || is.complex(pR);
	for(id in seq(2, len)) {
		p2 = p[[id]];
		if(is.numeric(p2) || is.complex(p2)) {
			if(isNum) { pR = pR * p2; }
			else pR = mult.sc.pm(pR, p2);
		} else {
			if(isNum) { isNum = FALSE; pR = mult.sc.pm(p2, pR); }
			else pR = mult.pm(pR, p2);
		}
	}
	return(pR);
}
mult.m = function(p, m) {
	# optimized for multiplication with monomial;
	if(nrow(m) > 1) stop("Cannot multiply with a polynomial!")
	id = match("coeff", names(m));
	if( ! is.na(id)) {
		p$coeff = p$coeff * m[,id];
		m = m[, - id, drop=FALSE];
	}
	xn.all = intersect(names(p), names(m));
	for(xn in xn.all) {
		p[, xn] = p[, xn] + m[1,xn];
	}
	isNew = ! (names(m) %in% xn.all);
	xn.all = names(m)[isNew];
	for(xn in xn.all) {
		p[, xn] = m[1,xn];
	}
	return(p);
}
mult.pm = function(p1, p2, sc=1) {
	if(is.character(p2)) {
		# multiply by a single variable:
		if(length(p2) > 1) stop("Not yet implemented!");
		idx = match(p2, names(p1));
		if(is.na(idx)) {
			p1[, p2] = 1;
		} else {
			p1[, idx] = p1[, idx] + 1;
		}
		return(p1);
	}
	# P1
	split.df = function(p.df) {
		p.l = lapply(colnames(p.df), function(name) p.df[,name]);
		names(p.l) = colnames(p.df);
		return(p.l);
	}
	# Zero:
	if(inherits(p1, "data.frame") && nrow(p1) == 0) return(p1);
	if( ! missing(p2) &&
		inherits(p2, "data.frame") && nrow(p2) == 0) return(p2);
	if(is.numeric(p1) || is.complex(p1) || ncol(p1) == 1) {
		if(inherits(p1, "data.frame")) p1 = p1$coeff;
		if(missing(p2)) p2 = p1; # Note: computes p1^2; (but may be dropped)
		if(is.list(p2)) return(mult.sc.pm(p2, p1));
		return(data.frame(coeff=p1*p2));
	}
	if( ! missing(p2)) {
	if(is.numeric(p2) || is.complex(p2) || ncol(p2) == 1) {
		if(inherits(p2, "data.frame")) p2 = p2$coeff;
		return(mult.sc.pm(p1, p2));
	}
	}
	# TODO: drop support for simple lists
	if(is.data.frame(p1)) {
		p1 = split.df(p1);
	}
	p1.b0 = p1$coeff;
	p1 = p1[ ! names(p1) %in% "coeff"];
	# P2
	if(missing(p2)) {
		p2.b0 = p1.b0;
		p2 = p1;
	} else {
		if(is.data.frame(p2)) {
			p2 = split.df(p2);
		}
		p2.b0 = p2$coeff;
		p2 = p2[ ! names(p2) %in% "coeff"];
	}
	# helper
	prod.b0 = function(p1, p2=p1) outer(p1, p2);
	# Adjust Vars
	vars = unique(c(names(p1), names(p2)));
	len  = if(length(p1) > 0) length(p1[[1]]) else 0; # ???
	for(v in vars[ ! vars %in% names(p1)]) p1[[v]] = rep(0, len);
	len  = if(length(p2) > 0) length(p2[[1]]) else 0;
	for(v in vars[ ! vars %in% names(p2)]) p2[[v]] = rep(0, len);
	# print(p1); print(p2);
	# Multiply
	p.m = lapply(vars, function(name) outer(p1[[name]], p2[[name]], function(i, j) i+j))
	p.b0 = as.vector(prod.b0(p1.b0, p2.b0));
	p.l = lapply(p.m, as.vector);
	p.v = as.data.frame(do.call(cbind, p.l));
	# p.v = cbind(p.v, coeff = p.b0);
	p.v$coeff = p.b0;
	p.r = aggregate0.pm(p.v);
	colnames(p.r) = c(vars, "coeff");
	if(sc != 1) p.r$coeff = p.r$coeff * sc;
	if( ! inherits(p.r$coeff, c("bigz", "bigq"))) {
		p.r = p.r[round0(p.r$coeff) != 0, ];
	}
	return(p.r);
}

### Power:
'^.pm' = function(p, n) {
	pow.pm(p, n=n);
}
pow.pm = function(p, n=2, do.order=TRUE, debug=TRUE) {
	if(n == 1) return(p);
	if(is.double(n) && (n == round(n))) n = as.integer(n);
	if( ! is.integer(n)) stop("n must be integer!")
	# Multiply
	# TODO: vectorize: as.integer(intToBits(3));
	p.r = NULL;
	p.pow = p;
	while (n > 0) {
		if(debug) print(paste0("Pow = ", n));
		if (n %% 2 == 1) {
			if(is.null(p.r)) p.r = p.pow else p.r = mult.pm(p.r, p.pow);
		}
		if(n == 1) break;
        p.pow = mult.pm(p.pow, p.pow);
        n = n %/% 2;
    }
	if(do.order) {
		x.name = names(p)[1];
		id = order(p.r[, x.name], decreasing=TRUE);
		p.r = p.r[id,];
	}
	return(p.r);
}

### All Powers:
pow.all = function(x, n, ...) {
	UseMethod("pow.all");
}
powAll.pm = function(p, n=2, asList=TRUE) {
	# TODO: Deprecate;
	pow.all.pm(p, n=n, asList=asList);
}
pow.all.pm = function(p, n, asList=TRUE) {
	if(n == 1) return(if(asList) list(p) else p);
	if(is.double(n) && (n == round(n))) n = as.integer(n);
	if( ! is.integer(n)) stop("n must be integer!")
	# Multiply
	p.r = list(p);
	p.pow = p;
	for(i in seq(2, n)) {
        p.pow = mult.pm(p.pow, p);
		p.r[[i]] = p.pow;
    }
	return(p.r);
}
pow.all.numeric = function(x, n, start.zero = FALSE) {
	x = rep(x, n);
	if(start.zero) x = c(1, x);
	cumprod(x);
}
pow.all.binom = function(n, coeff = c(1,1), names = c("x", "y")) {
	if(n < 0) stop("Invalid n!");
	if( ! inherits(n, "integer")) n = as.integer(n);
	lst = pascal(n);
	if(length(lst[[1]]) == 1) return(data.frame(coeff = 1));
	#
	c1 = rev(pow.all.numeric(coeff[1], n=n, start.zero = TRUE));
	c2 = pow.all.numeric(coeff[2], n=n, start.zero = TRUE);
	cc = c1 * c2;
	nms = c(names, "coeff");
	lst = lapply(seq(n), function(id) {
		coeff = lst[[id]] * cc[seq(1, id + 1)];
		tmp = data.frame(seq(id, 0, by = -1), seq(0, id), coeff = coeff);
		names(tmp) = nms;
		class(tmp) = c("pm", class(tmp));
		return(tmp);
	});
	return(lst);
}
mult.sc.pm = function(p, s, div=1, coeff.name="coeff") {
	# Multiplication by scalar
	# div = for numerical stability;
	if(is.data.frame(p)) {
		p[ , coeff.name] = p[ , coeff.name] * s;
		if(div != 1) p[ , coeff.name] = p[ , coeff.name] / div;
	} else if(is.list(p)) {
		# TODO: remove support for simple lists & reuse list;
		p[[coeff.name]] = p[[coeff.name]] * s;
		if(div != 1) p[[coeff.name]] = p[[coeff.name]] / div;
	} else stop("p must be a polynomial!")
	return(p);
}

# Pascal's Triangle
pascal = function(n) {
	if(n < 0) stop("Invalid n!");
	if( ! inherits(n, "integer")) n = as.integer(n);
	if(n == 0) return(list(1));
	#
	lst = list(c(1,1));
	if(n == 1) return(lst);
	for(i in seq(2, n)) {
		tmp = lst[[i - 1]];
		lst[[i]] = c(0, tmp) + c(tmp, 0);
	}
	return(lst);
}

### Simplify functions

### Simplify p: Powers & Coefficients
# - useful when solving: p = 0;
simplify.spm = function(p1, do.gcd=FALSE) {
	return(simplify.pm.pow(p1, do.gcd=do.gcd))
}
simplify.pm.pow = function(p, do.gcd=FALSE) {
	if(nrow(p) == 0) return(p);
	nms = names(p);
	nms = nms[ ! nms %in% "coeff"];
	for(nm in nms) {
		v.pow = min(p[,nm]);
		if(v.pow > 0) {
			p[,nm] = p[,nm] - v.pow;
		}
	}
	if(do.gcd && (xgcd <- gcd.vpm(p)) > 1) p$coeff = p$coeff / xgcd;
	return(p);
}
simplify.pm = function(p1, p2) {
	# simplify fractions: x^n1 / x^n2;
	if(is.null(p2)) return(simplify.spm(p1));
	com.nm = intersect(names(p1), names(p2));
	com.nm = com.nm[ ! com.nm %in% "coeff"];
	for(nm in com.nm) {
		v.pow = min(p1[,nm], p2[,nm]);
		if(v.pow > 0) {
			p1[,nm] = p1[,nm] - v.pow;
			p2[,nm] = p2[,nm] - v.pow;
		}
	}
	return(list(p1=p1, p2=p2));
}
### Simplify functions
reduce.pm = function(p) {
	# remove Monomials with coeff == 0;
	if( ! is.pm(p)) stop("p must be a Polynomial!")
	id = which(p$coeff != 0);
	return(p[id, , drop=FALSE]);
}
reduce0.pm = function(p) {
	return(p[p$coeff != 0, , drop=FALSE]);
}
reduce.var.pm  = function(p) drop.pm(p);
drop.pm = function(p) {
	# remove Vars with power == 0;
	id = match("coeff", names(p));
	if(ncol(p) <= 1) return(p);
	nc = rep(TRUE, ncol(p));
	nc[-id] = sapply(seq(ncol(p))[-id], function(id) any(p[,id] != 0));
	return(p[, nc, drop=FALSE]);
}
reduce.cpm = function(p, asBigNum=FALSE) {
	warning("Obsolete function: use reduce.coef.pm!");
	reduce.coef.pm(p, asBigNum=asBigNum);
}
reduce.coef.pm = function(p, asBigNum=FALSE) {
	# simplify coefficients
	xgcd = if(asBigNum) as.bigz(0) else 0;
	if(is.complex(p$coeff)) {
		xgcd = gcd.vpm.complex(p, xgcd);
	} else {
		xgcd = gcd.vpm(p, xgcd);
	}
	if(xgcd != 1) {
		p$coeff = p$coeff / xgcd;
		if(asBigNum) p$coeff = as.bigz(p$coeff);
	}
	return(p);
}

### Convert Coefficients
toDouble.pm = function(p, scale=1) {
	# convert Big Integers to some sensible Double;
	isZero = (p$coeff == 0);
	div = min(abs(p$coeff[ ! isZero]));
	div = div * scale;
	p$coeff = as.double(p$coeff / div);
	return(p);
}
toDouble.lpm = function(lp) {
	div = lapply(lp, function(p) {
		isZero = (p$coeff == 0);
		div = min(abs(p$coeff[ ! isZero]));
		matrix(div, ncol=1, nrow=1);
	});
	div = do.call(cbind, div); # workaround for bigz;
	div = max(div);
	for(id in seq(length(lp))) {
		lp[[id]]$coeff = as.double(lp[[id]]$coeff / div);
	}
	return(lp);
}
toBigz.pm = function(p) {
	if(is.complex(p$coeff)) {
		if(any(Im(p$coeff) != 0)) stop("Complex coefficients!");
		p$coeff = Re(p$coeff);
	}
	p$coeff = as.bigz(p$coeff);
	return(p);
}
as.bigz.pm = function(p) toBigz.pm(p);
as.numeric.pm = function(p) {
	p$coeff = as.numeric(p$coeff);
	return(p);
}

### Sum
### Helper functions
align.pm = function(p1, p2, align.names=TRUE, doReduce=TRUE) {
	# align columns of 2 data.frames for sum.pm();
	if(doReduce) {p1 = reduce.pm(p1); p2 = reduce.pm(p2);}
	n1 = names(p1); n2 = names(p2);
	### Coefficients
	n1 = n1[ ! n1 %in% "coeff"];
	n2 = n2[ ! n2 %in% "coeff"];
	xc = intersect(n1, n2);
	xall = union(n1, n2);
	pad.pm = function(p, vnew) {
		if(is.data.frame(p)) {
			p[, vnew] = 0;
		} else {
			len = length(p$coeff);
			zero = rep(0, len);
			for(nn in vnew) {
				p[[nn]] = zero;
			}
		}
		return(p);
	}
	# add missing variables
	if(length(xall) != length(xc)) {
		### p1
		n1new = n2[ ! n2 %in% n1];
		p1 = pad.pm(p1, n1new);
		### p2
		n2new = n1[ ! n1 %in% n2];
		p2 = pad.pm(p2, n2new);
	}
	#
	if(align.names) {
		id = match(names(p1), names(p2));
		list(p1=p1, p2=p2[id]);
	} else {
		list(p1=p1, p2=p2);
	}
}
sum.sc.pm = function(p, sc) {
	if(sc == 0) return(p);
	idCoeff = match("coeff", names(p));
	B0.pm = function() {
		sapply(seq(nrow(p)), function(nr) all(p[nr, - idCoeff] == 0));
	}
	isB0 = B0.pm();
	if(any(isB0)) {
		p$coeff[isB0][1] = p$coeff[isB0][1] + sc;
	} else {
		nr = nrow(p) + 1;
		p[nr, ] = rep(0, ncol(p));
		p[nr, "coeff"] = sc;
	}
	return(p);
}
sum.pm = function(p1, p2, doReduce=FALSE) {
	isDF1 = is.data.frame(p1); isDF2 = is.data.frame(p2);
	if(isDF1 && nrow(p1) == 0) {
		# TODO: if(doReduce);
		p2 = if(isDF2) reduce.pm(p2) else p2;
		return(p2);
	} else if(isDF2 && nrow(p2) == 0) {
		p1 = if(isDF1) reduce.pm(p1) else p1;
		return(p1);
	}
	if(is.numeric(p1) || is.complex(p1)) {
		if(isDF2) return(sum.sc.pm(p2, p1)) else return(p1 + p2);
	} else if(is.numeric(p2) || is.complex(p2)) {
		return(sum.sc.pm(p1, p2))
	}
	#
	l = align.pm(p1, p2, doReduce=doReduce); # no need to pre-reduce;
	p1 = l[[1]]; p2 = l[[2]];
	n1 = names(p1); n2 = names(p2);
	### to DF
	id = match(n2, n1);
	p = rbind(as.data.frame(p1), as.data.frame(p2)[,id]);
	### Sum
	p.r = aggregate0.pm(p);
	return(reduce0.pm(p.r));
}
sum.lpm = function(lp) {
	pR = data.frame();
	for(pd in lp) {
		pR = sum.pm(pR, pd);
	}
	return(pR);
}
add.pm = function(p1, p2) return(sum.pm(p1, p2));
add.lpm = function(lp) return(sum.lpm(lp));

### Diff
diff.pm = function(p1, p2) {
	if(is.pm(p2)) { p2$coeff = - p2$coeff; }
	else p2 = - p2;
	return(sum.pm(p1, p2));
}
diff.lpm = function(p1, lp) {
	for(pd in lp) {
		p1 = diff.pm(p1, pd);
	}
	return(p1);
}

########################

### Evaluations

### Evaluate Polynomial:
eval.pm = function(p, x, progress=FALSE) {
	# x = c(values of variables) OR
	# x = list(values of variables);
	# progress = currently not used;
	pP = p[, - which(names(p) == "coeff"), drop=FALSE];
	if(is.list(x) || ! is.null(names(x))) {
		len = sapply(x, length);
		if(any(len != 1)) stop("Each list element must have 1 entry!")
		nms = names(x);
		if( ! is.null(nms)) {
			nmsP = names(pP);
			id = pmatch(nmsP, nms);
			if(any(is.na(id))) stop(paste0("Variables missing: ", nmsP[is.na(id)]));
			x = x[id];
		}
		if(is.list(x)) {
			isBigz = sapply(x, function(x) inherits(x, c("bigz", "bigq")));
			if(any(isBigz)) {
				x = do.call(c, x);
			} else
				x = unlist(x, recursive=FALSE);
		}
	}
	eval.p = function(id) {
		idx = which(unlist(pP[id,]) != 0);
		if(length(idx) == 0) return(p$coeff[id]);
		prod(p$coeff[id], x[idx]^unlist(pP[id, idx]));
	}
	tmp = sapply(seq(nrow(p)), eval.p);
	if(inherits(p$coeff, c("bigz", "bigq"))) {
		tmp = do.call(c, tmp);
		# tmp = as.vector.bigz(tmp);
		print(tmp)
	}
	sum(tmp);
}

### Evaluate a list of polynomials
# - numeric values in the list are interpreted as the corresponding value;
eval.lpm = function(p, vals) {
	len = length(p);
	rez = sapply(seq(len), function(id) {
		p1 = p[[id]];
		if(is.numeric(p1)) return(p1);
		rez = eval.pm(p1, vals);
		return(rez);
	});
	return(rez);
}
### Evaluate with a set of solutions
eval.pm.sol = function(p, sol, ..., round0=TRUE, tol=1E-7) {
	vals = list(...);
	if(is.null(dim(sol)) || nrow(sol) == 1) {
		vals = c(as.list(sol), vals);
		rez  = eval.pm(p, vals);
		if(round0) rez = round0(rez, tol=tol);
		return(rez);
	}
	#
	len = nrow(sol);
	rez = sapply(seq(len), function(id) {
		iVal = c(as.list(sol[id,]), vals);
		eval.pm(p, iVal);
	})
	if(round0) rez = round0(rez, tol=tol);
	return(rez);
}

##################

### Solve Variable
# - list of polynomials is usually in ascending order of powers of substituted variables;
# - xn = vector of variables which are sequentially eliminated;
#   -- can be shorter than (number of polynomials - 1);
solve.lpm = function(..., xn, stop.at=NULL, asBigNum=FALSE) {
	pL  = list(...);
	len = length(pL);
	if(len < 2) {
		if(inherits(pL[[1]], c("pm", "data.frame"))) {
			return(pL);
		}
		if(inherits(pL[[1]], c("list"))) {
			pL  = pL[[1]];
			len = length(pL);
		} else stop("Not supported!");
	}
	#
	pR = list();
	lenX = min(len - 1, length(xn));
	for(id in seq(lenX)) {
		cat(paste0("\nStarting step: ", id, "\n"));
		if( ! is.null(stop.at) && id == lenX) {
			# Stop only during the last elimination;
			tmp = solve.pm(pL[[id+1]], pL[[id]], xn=xn[[id]], stop.at=stop.at, asBigNum=asBigNum);
		} else {
			tmp = solve.pm(pL[[id+1]], pL[[id]], xn=xn[[id]], asBigNum=asBigNum);
		}
		pR[[id]] = tmp;
		pL[[id+1]] = tmp$Rez;
		idS = id + 2;
		if(idS > len) next;
		for(id2 in seq(idS, len)) {
			if(is.na(match(xn[[id]], names(pL[[id2]])))) {
				warning(paste0("Missing Variable: ", xn[[id]], "; step = ", id));
				next;
			}
			# TODO: reuse powers of tmp$x0, tmp$div;
			pL[[id2]] = replace.fr.pm(pL[[id2]], tmp$x0, tmp$div, xn=xn[[id]]);
		}
	}
	lenX = lenX + 2; # Note: + 2 avoids repeating (lenX + 1);
	if(lenX <= len) {
		pR = c(pR, pL[seq(lenX, len)]);
	}
	return(pR);
}
solve.pm = function(p1, p2, by = xn, xn = NULL, stop.at=NULL, simplify=TRUE,
		asBigNum=FALSE, verbose=TRUE) {
	if(is.null(by)) stop("Missing variable name!");
	if(is.pm(by)) stop("Invalid variables: Did you mean to use solve.lpm()?");
	xn = by;
	#
	max1 = max(p1[,xn]); max2 = max(p2[,xn]);
	if(max2 == 0) stop("No variable!");
	if(max1 == 0) {
		# stop("No variable!");
		warning("No variable!");
		# Result may still be useful!
		return(list(Rez = p1, x0 = p2, div=NULL, xn=xn));
	}
	if(max1 < max2) {
		tmp = p1; p1 = p2; p2 = tmp;
		tmp = max1; max1 = max2; max2 = tmp;
	} else if(max1 == max2 && nrow(p1) < nrow(p2)) {
		tmp = p1; p1 = p2; p2 = tmp;
	}
	# for debugging:
	if( ! is.null(stop.at) && max2 == stop.at) return(list(p1, p2));
	split.pm = function(p, pow) {
		px = p[,xn];
		p2 = p[, - match(xn, names(p)), drop=FALSE];
		p2x = p2[px == pow, , drop=FALSE];
		p20 = p2[px != pow, , drop=FALSE]; # changed to (- coeff) in max2-section;
		return(list(p20, p2x));
	}
	if(max2 == 1) {
		lp2 = split.pm(p2, max2);
		if(nrow(lp2[[1]]) == 0) {
			print("Warning: x == 0!");
			p1 = p1[p1[,xn] == 0, , drop=FALSE];
			return(p1);
		}
		print(paste0("Substituting: Len = ", nrow(p1),
			"; Len = ", nrow(lp2[[1]]), " + ", nrow(lp2[[2]])));
		lp2[[2]]$coeff = - lp2[[2]]$coeff; # "-" !!!
		if(asBigNum) {
			if(verbose) print("Computing GCD!")
			xgcd = as.bigz(0);
			for(i in seq(nrow(lp2[[1]]))) xgcd = gcd(xgcd, lp2[[1]]$coeff[i]);
			for(i in seq(nrow(lp2[[2]]))) xgcd = gcd(xgcd, lp2[[2]]$coeff[i]);
			if(xgcd > 1) {
				if(verbose) print("Simplifying by GCD!")
				lp2[[1]]$coeff = as.bigz(lp2[[1]]$coeff / xgcd);
				lp2[[2]]$coeff = as.bigz(lp2[[2]]$coeff / xgcd);
			}
		}
		p1 = replace.fr.pm(p1, lp2[[1]], lp2[[2]], x=xn, pow=1, verbose=verbose);
		if(simplify) p1 = simplify.spm(p1, do.gcd=TRUE);
		return(list(Rez=p1, x0=lp2[[1]], div=lp2[[2]], xn=xn));
	}
	leading.pm = function(p, pow) {
		px = p[,xn];
		p2 = p[, - match(xn, names(p)), drop=FALSE];
		p2x = p2[px == pow,]; # TODO: drop = FALSE;
		return(p2x);
	};
	p1cf = leading.pm(p1, max1);
	p2cf = leading.pm(p2, max2);
	dmax = max1 - max2;
	if(dmax < 0) p2cf[,xn] = -dmax;
	if(dmax > 0) p1cf[,xn] = dmax;
	# TODO: gcd of coefficients & polynomials;
	p1 = sum.pm(mult.pm(p1, p2cf), mult.pm(p2, p1cf, -1));
	if(verbose) print(paste0("Max pow: ", max1, "; Len = ", nrow(p1)));
	if(simplify) { p1 = simplify.spm(p1); p2 = simplify.spm(p2); }
	return(solve.pm(p1, p2, xn=xn, stop.at=stop.at, simplify=simplify, asBigNum=asBigNum));
}

# p(r) => p(r^2), where r = root of p(x)
square.pm = function(p, xn) {
	xp = p[, xn];
	xp2 = xp %% 2;
	xp0 = (xp2 == 0);
	xp1 = ! xp0;
	if(all(xp0) || all(xp1)) {
		pR = mult.pm(p, p);
	} else {
		p1 = p[xp0, , drop=FALSE];
		p1 = mult.pm(p1, p1);
		p2 = p[xp1, , drop=FALSE];
		p2 = mult.pm(p2, p2);
		pR = diff.pm(p1, p2);
	}
	pR[, xn] = pR[, xn] / 2;
	maxPow = max(pR[, xn]);
	isMax  = which(pR[, xn] == maxPow)
	if(length(isMax) == 1) {
		if(pR$coeff[isMax] < 0) pR$coeff = - pR$coeff;
	}
	return(pR);
}

#######################
#######################

### Differentiation:
# - moved to file:
#   Polynomials.Helper.D.R;


### Polynomial Generators:
# - moved Generators to file:
#   Polynomials.Helper.Generators.R;


#############
### Tests ###
#############

# - moved to file:
#   Polynomials.Helper.Tests.R;

