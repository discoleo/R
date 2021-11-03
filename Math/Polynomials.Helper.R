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
# - Parser moved to new file:
#   Polynomials.Helper.Parser.R;
# - Derivation moved to new file:
#   Polynomials.Helper.D.R;
# - Generic Solvers moved to new file:
#   Polynomials.Helper.Solvers.R;


#######################

library(polynom)
library(pracma)

### helper Functions

# required:
source("Polynomials.Helper.Parser.R")
source("Polynomials.Helper.D.R")
# D: required in factorize.p();


### Basic Algebra

### Round to 0
round0 = function(m, tol=1E-7) {
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
round0.p = function(p, tol=1E-7) {
	p = round0(as.vector(p), tol=tol)
	class(p) = "polynomial"
	return(p)
}
round0.pm = function(p, tol=1E-7) {
	p$coeff = round0(p$coeff, tol=tol);
	return(p);
}
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
	# roots for Class 2 polynomials;
	m = unity(n+1, all=T)[-1]; # exclude 1;
	r = sapply(seq(n), function(id) sum(s * m[id]^(0:n)))
	r = round0(r)
}

sort.sol = function(sol, useRe=TRUE, ncol=1, digits=5) {
	if(useRe) {
		id = order(
			abs(round(sol[,ncol], digits)),
			round(Re(sol[,ncol]), digits) );
	} else {
		id = order(abs(round(sol[,ncol], digits)) );
	}
	return(sol[id,]);
}
isConj.f = function(x, y, tol=1E-3) {
	isConj = (abs(Re(x) - Re(y)) < tol) & (abs(Im(x) + Im(y)) < tol);
	return(isConj);
}

###################
### Polynomials ###
###################

# isNonZero:
isNZ.pm = function(p) {
	is.data.frame(p) && (nrow(p) > 0);
}

### Simple Multiplication
# p1, p2 = simple vectors of coefficients;
mult.p = function(p1, p2) {
	p.m = outer(p1, p2)
    p = as.vector(tapply(p.m, row(p.m) + col(p.m), sum))
	return(p)
}
# Compute polynomial from Roots:
poly.calc.mpfr = function(x, bits=120, tol=1E-7) {
	one  = mpfr(1, precBits=bits);
	zero = mpfr(0, precBits=bits);
	p  = mpfr2array(c(Re=one, Im=zero), c(1,2));
	p0 = mpfr2array(c(Re=zero, Im=zero), c(1,2));
	xRe = Re(x); xIm = Im(x);
	for (i in seq(length(x))) {
		px = mult.mpfr(p[,1], p[,2], xRe[i], xIm[i]);
		p = rbind(p0, p);
		p  = p - rbind(px, p0);
	}
	p = round0(p, tol=tol);
	return(p);
}

##################

### Multi-variable Polynomials:
# p = multi-variable polynomial;
#  => data.frame with exponents of variables;
#  => coeff = column with the coefficients;
# Note:
# - initial idea was to allow also basic lists,
#   but most functions would work only with data frames!
aggregate0.pm = function(p) {
	p.r = aggregate(coeff~., p, sum);
	return(p.r);
}
### Multiplication
mult.all.pm = function(p) return(mult.lpm(p));
mult.lpm = function(p) {
	if( ! is.list(p)) stop("p must be a list of polynomials");
	len = length(p);
	pR = p[[1]];
	for(id in seq(2, len)) {
		p2 = p[[id]];
		if(is.numeric(p2)) {
			if(is.numeric(pR)) pR = pR * p2
			else pR = mult.sc.pm(pR, p2);
		} else {
			if(is.numeric(pR)) pR = mult.sc.pm(p2, pR)
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
	# P1
	split.df = function(p.df) {
		p.l = lapply(colnames(p.df), function(name) p.df[,name]);
		names(p.l) = colnames(p.df);
		return(p.l);
	}
	if(is.numeric(p1) || is.complex(p1) || ncol(p1) == 1) {
		if( ! (is.numeric(p1) || is.complex(p1)) ) p1 = p1$coeff;
		if(missing(p2)) p2 = p1;
		if(is.list(p2)) return(mult.sc.pm(p2, p1));
		return(data.frame(coeff=p1*p2));
	}
	if( ! missing(p2)) {
	if(is.numeric(p2) || is.complex(p2) || ncol(p2) == 1) {
		if( ! (is.numeric(p2) || is.complex(p2)) ) p2 = p2$coeff;
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
		if(debug) print(n);
		if (n %% 2 == 1) {
			if(is.null(p.r)) p.r = p.pow else p.r = mult.pm(p.r, p.pow);
		}
		if(n == 1) break;
        p.pow = mult.pm(p.pow);
        n = n %/% 2;
    }
	if(do.order) {
		x.name = names(p)[1];
		id = order(p.r[,x.name], decreasing=TRUE);
		p.r = p.r[id,];
	}
	return(p.r);
}
powAll.pm = function(p, n=2, asList=TRUE) {
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
simplify.spm = function(p1, do.gcd=FALSE) {
	nms = names(p1);
	nms = nms[ ! nms %in% "coeff"];
	for(nm in nms) {
		v.pow = min(p1[,nm]);
		if(v.pow > 0) {
			p1[,nm] = p1[,nm] - v.pow;
		}
	}
	if(do.gcd && (xgcd <- gcd.vpm(p1)) > 1) p1$coeff = p1$coeff / xgcd;
	return(p1);
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
	id = which(p$coeff != 0)
	if(is.data.frame(p)) {
		return(p[id, , drop=FALSE]);
	}
	# old list code: but still assumes equal length;
	p = lapply(p, function(m) m[id]);
	return(p);
}
reduce.var.pm = function(p) {
	# remove Vars with power == 0;
	id = match("coeff", names(p));
	nc = rep(TRUE, ncol(p));
	nc[-id] = sapply(seq(ncol(p))[-id], function(id) any(p[,id] != 0));
	return(p[, nc, drop=FALSE]);
}
reduce.cpm = function(p, asBigNum=FALSE) {
	# simplify coefficients
	xgcd = if(asBigNum) as.bigz(0) else 0;
	xgcd = gcd.vpm(p, xgcd);
	if(xgcd != 1) {
		p$coeff = p$coeff / xgcd;
		if(asBigNum) p$coeff = as.bigz(p$coeff);
	}
	return(p);
}
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
	p$coeff = as.bigz(p$coeff);
	return(p);
}
as.numeric.pm = function(p) {
	p$coeff = as.numeric(p$coeff);
	return(p);
}
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
	return(reduce.pm(p.r));
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
	p2$coeff = - p2$coeff;
	return(add.pm(p1, p2));
}
diff.lpm = function(p1, lp) {
	for(pd in lp) {
		p1 = diff.pm(p1, pd);
	}
	return(p1);
}

### Shift variables:
shift.pm = function(p, val, x="x", tol=1E-10) {
	len = length(val);
	if(len > 1) {
		if(length(x) == 1) {
			warning("Same variable used!");
			return(shift.pm(p, round0(sum(val), tol=tol), x=x, tol=tol))
		}
		for(i in seq(len)) {
			p = shift.pm(p, val = val[i], x=x[i], tol=tol);
		}
		return(p);
	}
	# TODO: efficient algorithm;
	x.new = paste0(x, "__sh");
	new.x = function(x)	{
		x.df = data.frame(x=1:0, x.new=0:1, coeff=1);
		names(x.df)[1:2] = c(x, x.new);
		return(x.df);
	}
	p = replace.pm(p, new.x(x), x, pow=1);
	p = replace.pm(p, val, x=x.new, pow=1);
	idv = match(x.new, names(p));
	if( ! is.na(idv) && all(p[, idv] == 0)) p = p[, - idv];
	if( ! is.null(tol)) p$coeff = round0(p$coeff, tol=tol);
	return(p);
}
replace.withVal.pm = function(p, x, pow=1, val, simplify=TRUE) {
	if(length(val) > 1) {
		len = length(val);
		if(length(pow) == 1) pow = rep(pow, len);
		if(length(x) == 1) {x = rep(x, len); warning("Same variable used!");}
		for(i in seq(len)) {
			p = replace.withVal.pm(p, x=x[i], pow=pow[i], val=val[i], simplify=simplify);
		}
		return(p);
	}
	if(val == 0) {
		if(pow != 1) {
			warning("Only some terms will be replaced with 0!");
			isZero = p[,x] >= pow;
			p = p[ ! isZero, , drop=FALSE];
		} else p = p[p[,x] == 0, , drop=FALSE];
		return(p);
	}
	xpow = p[,x];
	p[,x] = if(pow == 1) 0 else xpow %% pow;
	hasX  = xpow >= pow; xpow = xpow[hasX];
	# Replace with value:
	x.pow = if(pow == 1) xpow else xpow %/% pow;
	xpow.unq = sort(unique(x.pow));
	id = match(x.pow, xpow.unq);
	xval = val^xpow.unq;
	p[hasX, "coeff"] = p[hasX, "coeff"] * xval[id];
	p = aggregate0.pm(p);
	p = reduce.pm(p);
	if(simplify) p = reduce.var.pm(p);
	return(p)
}
replace.pm = function(p1, p2, x, pow=1) {
	# replace x^pow by p2;
	idx = match(x, names(p1));
	if(any(is.na(idx))) {
		warning(paste0("Polynomial does NOT contain variable: ", x));
		return(p1);
	}
	if(is.numeric(p2) || is.complex(p2)) return(replace.withVal.pm(p1, x=x, pow=pow, val=p2));
	# xPow
	rpow = if(pow == 1) p1[,idx] else p1[,idx] %/% pow;
	p1[,idx] = if(pow == 1) 0 else p1[,idx] %% pow;
	max.pow = max(rpow);
	p2.pows = list(p2);
	if(max.pow > 1) {
		if(ncol(p2) > 1) {
			for(ipow in seq(2, max.pow)) {
				p2.pows[[ipow]] = mult.pm(p2.pows[[ipow - 1]], p2);
			}
		} else {
			p2.pows = sum(p2$coeff)^seq(max.pow);
		}
	}
	pR = data.frame();
	for(nr in seq(nrow(p1))) {
		pR = if(rpow[nr] == 0) add.pm(pR, p1[nr,])
			else add.pm(pR, mult.pm(p2.pows[[rpow[nr]]], p1[nr,]));
	}
	return(reduce.var.pm(pR));
}
replace.fr.pm = function(p1, p2, p2fr, x, pow=1) {
	# replace x^pow by p2/p2fr;
	idx = match(x, names(p1));
	if(is.na(idx)) stop(paste0("Polynomial does NOT contain variable: ", x));
	# xPow
	rpow = if(pow == 1) p1[,idx] else p1[,idx] %/% pow;
	p1[,idx] = if(pow == 1) 0 else p1[,idx] %% pow;
	max.pow = max(rpow);
	p2.pows = list(p2); # powers of p2
	p2fr.pows = list(p2fr); # powers of p2fr
	if(max.pow > 1) {
		for(ipow in seq(2, max.pow)) {
			print(paste0("Pow = ", ipow))
			p2.pows[[ipow]] = mult.pm(p2.pows[[ipow - 1]], p2);
			p2fr.pows[[ipow]] = mult.pm(p2fr.pows[[ipow - 1]], p2fr);
		}
		print("Finished Powers!");
	}
	p2m = c(tail(p2fr.pows, 1), p2.pows);
	print("Starting cross-multiplication:");
	# pre-align p1 & p2: verify if useful?
	pall = align.pm(p1, p2);
	p1 = pall[[1]];
	if(max.pow > 1) {
	for(ipow in seq(1, max.pow - 1)) {
		p2m[[ipow + 1]] = mult.pm(p2m[[ipow + 1]], p2fr.pows[[max.pow - ipow]]);
		p2m[[ipow + 1]] = align.pm(p1, p2m[[ipow + 1]], doReduce=FALSE)[[2]];
	}
	}
	print("Finished cross-multiplication!");
	pR = data.frame();
	for(nr in seq(nrow(p1))) {
		ipow = rpow[nr];
		tmp = mult.pm(p2m[[ipow + 1]], p1[nr,]);
		# TODO: separate sum
		pR = sum.pm(pR, tmp);
	};
	return(reduce.var.pm(pR));
}
eval.pm = function(p, x, progress=FALSE) {
	# x = c(values of variables)
	# TODO: list!
	pP = p[, - which(names(p) == "coeff"), drop=FALSE];
	eval.p = function(id) {
		idx = which(pP[id,] != 0);
		if(length(idx) == 0) return(p$coeff[id]);
		prod(x[idx]^pP[id, idx], p$coeff[id]);
	}
	sum(sapply(seq(nrow(p)), eval.p))
}
toPolar.mpfr = function(x, len=1, bits=120, tol=1E-10) {
	re = mpfr(Re(x), bits); im = mpfr(Im(x), bits);
	r = sqrt(re^2 + im^2);
	pib = Const("pi", bits); pih = pib / 2;
	# seems NO difference between asin & atan versions;
	th  = asin(abs(im/r));
	if(re == 0) {th = pih; if(im < 0) th = - th;}
	else if(re < 0) {
		th = pib + if(im > 0) - th else th;
	} else if(im < 0) {
		th = -th;
	}
	if(len > 1) {
		r  = r^seq(len);
		th = th * seq(len);
	}
	re = r * cos(th); im = r * sin(th);
	return(cbind(Re=re, Im=im));
}
mult.mpfr = function(re1, im1, re2, im2) {
	reN = re1 * re2 - im1 * im2;
	imN = re1 * im2 + im1 * re2;
	return(cbind(Re=reN, Im=imN));
}
eval.cpm = function(p, x, bits=120, tol=1E-10, doPolar=TRUE, progress=FALSE) {
	# uses the Rmpfr package;
	# currently assumes that only coeffs are big and
	# have an impact on numeric stability;
	pP = p[, - which(names(p) == "coeff"), drop=FALSE];
	# pow = lapply(seq(ncol(pP)), function(nc) sort(unique(pP[,nc])));
	# currently only max:
	pow = lapply(seq(ncol(pP)), function(nc) max(pP[,nc]));
	xpows = lapply(seq(length(x)), function(id) {
		x0 = round0(x[id], tol=tol);
		len = tail(pow[[id]], 1);
		if(is.complex(x0) && Im(x0) != 0) {
			div = 1;
			# polar coordinates:
			# - but less accuracy with certain complex numbers;
			# - needed when r^max.pow overflows;
			if(doPolar) {
				re = mpfr(Re(x0), bits); im = mpfr(Im(x0), bits);
				r = sqrt(re^2 + im^2);
				pib = Const("pi", bits); pih = pib / 2;
				# seems NO difference between asin & atan versions;
				th  = asin(abs(im/r));
				if(re == 0) {th = pih; if(im < 0) th = - th;}
				else if(re < 0) {
					th = pib + if(im > 0) - th else th;
				} else if(im < 0) {
					th = -th;
				}
				r  = r^seq(len);
				th = th * seq(len);
				re = r * cos(th); im = r * sin(th);
			} else {
				x = x0^seq(len);
				re = Re(x); im = Im(x);
				re = mpfr(re * div, bits);
				im = mpfr(im * div, bits);
			}
			return(cbind(Re=re, Im=im, Div=div));
		} else {
			x0 = mpfr(Re(x0), bits);
			x = x0^seq(len); # power 0 NOT needed;
			if(x0 == 0) {
				return(cbind(Re=x, Im=0, Div=1));
			} else {
				div = 1; # 12 - round(log(abs(x)) / log(10));
				return(cbind(Re=x, Im=0, Div=div));
			}
		}
	})
	if(progress) cat("Processing row:\n");
	eval.p = function(id) {
		if(progress && id %% 16 == 1) cat(paste0(id, if(id %% 96 == 1) "\n" else ", "));
		idx = which(pP[id,] != 0);
		xpows = xpows[idx]; lenv = length(idx);
		re = mpfr2array(sapply(seq(lenv), function(id2) xpows[[id2]][pP[id, idx[id2]], 1]), lenv);
		im = mpfr2array(sapply(seq(lenv), function(id2) xpows[[id2]][pP[id, idx[id2]], 2]), lenv);
		if(length(re) > 0) {
		while(TRUE) {
			len = length(re);
			if(len == 1) break;
			iend = (len %% 2);
			i  = seq(1, len - iend, by=2);
			re1 = re[i] * re[i+1] - im[i] * im[i+1];
			im1 = re[i] * im[i+1] + im[i] * re[i+1];
			if(iend > 0) {
				re.last = re[len]; im.last = im[len];
				re2 = re1[1] * re.last - im1[1] * im.last;
				im2 = re1[1] * im.last + im1[1] * re.last;
				re1[1] = re2; im1[1] = im2;
			}
			re = re1; im = im1;
		}
		} else {re = 1; im = 0;}
		re = re * p$coeff[id]; im = im * p$coeff[id];
		sol = mpfr2array(c(Re=re, Im=im, Div=1), c(3));
		return(sol);
	}
	if(progress) cat("\n");
	sol = sapply(seq(nrow(p)), eval.p);
	sdim = attr(sol, "dim"); sol = mpfr2array(t(sol), rev(sdim));
	sol = apply(sol, 2, sum);
	return(sol);
}
## === Div ===
div.pm = function(p1, p2, by="x", debug=TRUE) {
	# very simple division
	xn = by[1];
	idx2 = match(xn, names(p2));
	if(is.na(idx2)) stop(paste0("P2 must contain the variable: ", xn));
	idx1 = match(xn, names(p1));
	if(is.na(idx1)) stop(paste0("P1 must contain the variable: ", xn));
	if( ! is.data.frame(p2)) p2 = as.data.frame(p2);
	#
	xpow2 = max(p2[,idx2]);
	pDx = p2[p2[,idx2] == xpow2, ];
	idc2 = match("coeff", names(p2));
	idc1 = match("coeff", names(p1));
	c2 = pDx[,idc2];
	pRez = as.data.frame(array(0, c(0,2)));
	names(pRez) = c(xn, "coeff");
	#
	idn = match(names(pDx)[-idc2], names(p1));
	print(idn);
	if(any(is.na(idn))) stop(paste0("No matching variables: ", names(pDx)[is.na(idn)]));
	if(nrow(pDx) == 1) {
		while(TRUE) {
			if(nrow(p1) == 0) break;
			xpow1 = max(p1[,xn]);
			if(xpow1 < xpow2) break;
			px1 = p1[p1[,xn] == xpow1, ];
			for(nc in seq_along(idn)) {
				px1[, idn[nc]] = px1[, idn[nc]] - pDx[, nc];
			}
			px1[, idc1] = px1[, idc1] / c2;
			pRez = add.pm(pRez, px1);
			p1 = diff.pm(p1, mult.pm(px1, p2));
		}
	} else {
		if(length(by) == 1) stop("Not yet implemented!")
		while(TRUE) {
			if(nrow(p1) == 0) break;
			xpow1 = max(p1[,xn]);
			if(xpow1 < xpow2) break;
			pCoeff = div.pm(p1[p1[,xn] == xpow1, ], p2[p2[,xn] == xpow2,], by=by[2])$Rez;
			print(pCoeff)
			pRez = add.pm(pRez, pCoeff);
			p1 = diff.pm(p1, mult.pm(pCoeff, p2));
		}
	}
	if(debug) {
		if(nrow(p1) > 0) print("Not divisible!")
		else print("Divisible!");
	}
	return(list(Rez=pRez, Rem=p1));
}
divByZero.pm = function(p1, pDiv, xn="x", asBigNum=TRUE) {
	return(gcd.exact.p(p1, pDiv, xn=xn, asBigNum=asBigNum, doGCD=FALSE))
}
divOK.pm = function(p, div, xn="x", warn=TRUE) {
	pR = div.pm(p, div, xn);
	if(nrow(pR$Rem) == 0) {
		return(list(Rez=pR$Rez, isDiv=TRUE));
	} else if(warn) warning("Division went wrong!");
	return(list(Rez=p, isDiv=FALSE));
}
gcd.vpm = function(p, xgcd=0) {
	for(i in seq(nrow(p))) xgcd = gcd(xgcd, p$coeff[i]);
	return(xgcd);
}
gcd.pm = function(p1, p2, by="x", div.sc=1) {
	if(missing(p2)) return(gcd.vpm(p1)); # scalar gcd on coeff;
	# basic implementation without much thought!
	if(max(p2[,by]) > max(p1[,by])) { tmp = p1; p1 = p2; p2 = tmp; }
	pR = div.pm(p1, p2, by=by);
	if(nrow(pR$Rem) == 0) return(pR);
	while(TRUE) {
		p1 = diff.pm(p1, mult.pm(p2, pR$Rez));
		p1$coeff = round0(p1$coeff);
		p1 = p1[p1$coeff != 0, ];
		print(mult.sc.pm(p1, 1, div.sc));
		n1 = max(p1[,by]); n2 = max(p2[,by]);
		if(n1 == 0 || n2 == 0) return(pR);
		if(n2 > n1) { tmp = p1; p1 = p2; p2 = tmp; }
		else if(n1 == n2 && max(abs(p1$coeff[p1[,by] == n1])) < max(abs(p2$coeff[p2[,by] == n1]))) {
			# non-robust for multi-variate polynomials!
			tmp = p1; p1 = p2; p2 = tmp;
		}
		pR = div.pm(p1, p2, by=by);
		if(nrow(pR$Rem) == 0) return(pR);
	}
	return(pR);
}
gcd.exact.p = function(p1, p2, xn="x", asBigNum=TRUE, doGCD=TRUE) {
	# exact implementation: only univariate polynomials;
	if( ! doGCD) fact = if(asBigNum) as.bigz(1) else 1;
	while(TRUE) {
		n1 = max(p1[,xn]); n2 = max(p2[,xn]);
		if( ! doGCD && n1 < n2) return(list(p=p1, f=fact));
		if(doGCD && (n1 < n2)) {
			tmp = p1; p1 = p2; p2 = tmp;
			tmp = n1; n1 = n2; n2 = tmp;
		}
		c1 = p1$coeff[p1[,xn] == n1];
		c2 = p2$coeff[p2[,xn] == n2];
		div = gcd(c1, c2);
		if(div != 1) {
			c1 = c1 / div; c2 = c2 / div;
			if(asBigNum) {c1 = as.bigz(c1); c2 = as.bigz(c2);}
		}
		p1m = p1; p1m$coeff = p1m$coeff * c2;
		p2m = p2; p2m$coeff = p2m$coeff * c1;
		dn = n1 - n2;
		if(dn > 0) p2m[,xn] = p2m[,xn] + dn
		else if(dn < 0) p1m[,xn] = p1m[,xn] - dn;
		#
		dp = diff.pm(p1m, p2m);
		if(nrow(dp) == 0) {
			print("Factor found!");
			if(doGCD) return(p2) else return(list(p=0, f=0));
		}
		# simplify the coefficients: robust for BigNumbers;
		xgcd = gcd.vpm(dp, xgcd=dp$coeff[1]);
		if(xgcd != 1) {
			dp$coeff = dp$coeff / xgcd;
			if(asBigNum) dp$coeff = as.bigz(dp$coeff);
		}
		if( ! doGCD) fact = fact * c2 / xgcd;
		n0 = max(dp[, xn, drop=TRUE]);
		print(paste0("Pow = ", n0, ", Len = ", nrow(dp)));
		if(n0 == 0) {
			print("Not divisible!");
			if(doGCD) return(dp);
			return(list(p=dp, f=fact));
		}
		p1 = dp;
	}
}
### Solve Variable
solve.pm = function(p1, p2, xn, stop.at=NULL, simplify=TRUE, asBigNum=FALSE) {
	max1 = max(p1[,xn]); max2 = max(p2[,xn]);
	if(max1 == 0) stop("No variable!")
	if(max2 == 0) stop("No variable!")
	if(max1 < max2) {
		tmp = p1; p1 = p2; p2 = tmp;
		tmp = max1; max1 = max2; max2 = tmp;
	} else if(max1 == max2 && nrow(p1) < nrow(p2)) {
		tmp = p1; p1 = p2; p2 = tmp;
	}
	if( ! is.null(stop.at) && max2 == stop.at) return(list(p1, p2));
	split.pm = function(p, pow) {
		px = p[,xn];
		p2 = p[, - match(xn, names(p)), drop=FALSE];
		p2x = p2[px == pow,];
		p20 = p2[px != pow,]; # changed to (- coeff) in max2-section;
		return(list(p20, p2x));
	}
	if(max2 == 1) {
		lp2 = split.pm(p2, max2);
		if(nrow(lp2[[1]]) == 0) {
			print("Warning: x == 0!");
			p1 = p1[p1[,xn] == 0,];
			return(p1);
		}
		print(paste0("Substituting: Len = ", nrow(p1),
			"; Len = ", nrow(lp2[[1]]), " + ", nrow(lp2[[2]])));
		lp2[[2]]$coeff = - lp2[[2]]$coeff; # "-" !!!
		if(asBigNum) {
			print("Computing GCD!")
			xgcd = as.bigz(0);
			for(i in seq(nrow(lp2[[1]]))) xgcd = gcd(xgcd, lp2[[1]]$coeff[i]);
			for(i in seq(nrow(lp2[[2]]))) xgcd = gcd(xgcd, lp2[[2]]$coeff[i]);
			if(xgcd > 1) {
				print("Simplifying by GCD!")
				lp2[[1]]$coeff = as.bigz(lp2[[1]]$coeff / xgcd);
				lp2[[2]]$coeff = as.bigz(lp2[[2]]$coeff / xgcd);
			}
		}
		p1 = replace.fr.pm(p1, lp2[[1]], lp2[[2]], x=xn, pow=1);
		if(simplify) p1 = simplify.spm(p1);
		return(list(Rez=p1, x0=lp2[[1]], div=lp2[[2]], xn=xn));
	}
	leading.pm = function(p, pow) {
		px = p[,xn];
		p2 = p[, - match(xn, names(p)), drop=FALSE];
		p2x = p2[px == pow,];
		return(p2x);
	};
	p1cf = leading.pm(p1, max1);
	p2cf = leading.pm(p2, max2);
	dmax = max1 - max2;
	if(dmax < 0) p2cf[,xn] = -dmax;
	if(dmax > 0) p1cf[,xn] = dmax;
	# TODO: gcd of coefficients & polynomials;
	p1 = sum.pm(mult.pm(p1, p2cf), mult.pm(p2, p1cf, -1));
	print(paste0("Max pow: ", max1, "; Len = ", nrow(p1)));
	if(simplify) { p1 = simplify.spm(p1); p2 = simplify.spm(p2); }
	return(solve.pm(p1, p2, xn=xn, stop.at=stop.at, simplify=simplify, asBigNum=asBigNum));
}

### Factorize
factorize.p = function(p, xn="x", f.all=FALSE, asBigNum=TRUE, file="_R.Temp.") {
	# factorize.all = FALSE
	# - p1 is usually sufficient;
	# - dos NOT handle: p1 * p2^3*p3^4 or p1^2*p2^3;
	id = match(xn, names(p));
	if(is.na(id)) stop("Variable NOT present!");
	lvl = 1; # level of factorization;
	doSave = ! (is.null(file) || is.na(file));
	rez = list();
	dp = dp.pm(p, xn);
	if(nrow(dp) == 0 || ncol(dp) < 2) return(list(list(GCD=NULL, p1=p)));
	while(TRUE) {
		# Step 1: GCD
		cat("\n"); print(paste0("Level = ", lvl));
		pGCD = gcd.exact.p(p, dp, xn, asBigNum=asBigNum);
		pGCD = reduce.var.pm(pGCD);
		if(nrow(pGCD) < 1 || ncol(pGCD) < 2) break;
		id = match(xn, names(pGCD));
		if(is.na(id)) break;
		# Leading Sign:
		isMaxPow = (pGCD[,id] == max(pGCD[,id]));
		if(pGCD$coeff[isMaxPow][[1]] < 0) pGCD$coeff = - pGCD$coeff;
		if(doSave) write.csv(pGCD, file=paste0(file, "GCD.", lvl, ".csv"), row.names=FALSE);
		# Step 2:
		p.all = div.pm(p, pGCD, xn)$Rez;
		if(asBigNum) {
			if(all(denominator(p.all$coeff) == 1)) p.all$coeff = as.bigz(p.all$coeff)
			else print("Warning: some Denominators != 1!")
		}
		if(doSave) write.csv(p.all, file=paste0(file, "ALL.", lvl, ".csv"), row.names=FALSE);
		p.minus1 = gcd.exact.p(pGCD, p.all, xn, asBigNum=asBigNum);
		# TODO: IF(p.minus1 == p.all) => multiplicity!
		p1 = div.pm(p.all, p.minus1, xn)$Rez;
		if(doSave) write.csv(p1, file=paste0(file, "p1.", lvl, ".csv"), row.names=FALSE);
		#
		rez[[lvl]] = list();
		rez[[lvl]][["GCD"]] = pGCD;
		rez[[lvl]][["p1"]]  = p1;
		rez[[lvl]][["All"]] = p.all;
		if( ! f.all) break;
		lvl = lvl + 1;
		p = pGCD;
		dp = dp.pm(pGCD, xn)
		if(nrow(dp) < 1 || ncol(dp) < 2) break;
		if(is.na(match(xn, names(dp)))) break;
	}
	return(rez);
}

### Derivation:
# - moved to file:
#   Polynomials.Helper.D.R;


#############

#############
### Print ###

order.df = function(x, decreasing=TRUE) {
	order.s = function(...) order(..., decreasing=decreasing);
	id = do.call(order.s, x);
	return(id);
}
### TODO: update use of sort.pm() everywhere!
sort.pm = function(p, xn=NULL, sort.coeff, xn2=NULL) {
	### Special Cols:
	# TODO: different approach;
	# - over xn: c(1,2,3,4) = Sum, Max, Min, MinNZ; (IF length(xn) > 1)
	# - over all: c(5,6,7,8,9) = SumAll, MaxAll, MinAll, MinNZAll, Coeff;
	isM = ( ! is.null(xn) && length(xn) > 1); # isMultiple
	if( ! is.null(xn2)) xn = c(xn, xn2);
	if(missing(sort.coeff)) {
		sort.coeff = if(isM) c(1,2, seq(10, length.out=length(xn)), 5,6)
			else if(is.null(xn)) c(1,2) else c(seq(6, length.out=length(xn)), 1,2);
	}
	pP = p[, - which(names(p) == "coeff"), drop=FALSE];
	summary.sort = function(p, FUN=sum) sapply(seq(nrow(p)), function(id) FUN(unlist(p[id, , drop=TRUE])));
	to.df = function(i, p, FUN) if(any(sort.coeff == i)) summary.sort(p, FUN) else rep(0, nrow(p));
	if(isM) {
		pP.pp = pP[, xn, drop=FALSE];
		s.df = to.df(1, pP.pp, sum);
		s.df = cbind(s.df, to.df(2, pP.pp, max));
		s.df = cbind(s.df, to.df(3, pP.pp, min));
		s.df = cbind(s.df, to.df(4, pP.pp, function(x) min(x[x != 0])));
		id0 = 4;
		s.df = cbind(s.df, to.df(1 + id0, pP, sum));
	} else {
		id0 = 0;
		s.df = to.df(1, pP, sum);
	}
	s.df = cbind(s.df, to.df(2 + id0, pP, max));
	s.df = cbind(s.df, to.df(3 + id0, pP, min));
	s.df = cbind(s.df, to.df(4 + id0, pP, function(x) min(x[x != 0]) ));
	s.df = cbind(s.df, to.df(5 + id0, p[, "coeff", drop=FALSE], function(x) - abs(x) ));
	if( ! is.null(xn)) s.df = cbind(s.df, pP[, xn, drop=FALSE]);
	s.df = s.df[, sort.coeff, drop=FALSE];
	id = order.df(s.df, decreasing=TRUE);
	p = p[id,]; rownames(p) = seq(nrow(p));
	return(p)
}
# Print multi-variable Poly
print.monome = function(name, p) {
	v = p[,name];
	v.r = rep("", length(v));
	isCoeff = v != 1 & v != 0;
	v.r[isCoeff] = paste0(name, "^", v[isCoeff]);
	v.r[v == 1] = name;
	return(v.r);
}
print.pm = function(...) {
	ch = as.character.pm(...);
	print(ch);
	invisible(ch);
}
print.p = function(...) {
	ch = as.character.pm(...);
	print(ch);
	invisible(ch);
}
# TODO: format.complex();
sort.simple.pm = function(p, leading=1, do.rev=FALSE, sort.order=TRUE) {
	if(length(leading) == 1) {
		p = p[order(p[, leading], decreasing=sort.order), , drop=FALSE];
	} else {
		coeff.df = data.frame(vs=apply(p[, leading, drop=FALSE], 1, sum));
		coeff.df = cbind(coeff.df, (p[, leading]));
		order.s = function(...) order(..., decreasing=sort.order);
		id = do.call(order.s, coeff.df);
		p = p[id,];
		if(do.rev) leading = rev(leading);
	}
	p = cbind(p[,-leading, drop=FALSE], p[,leading, drop=FALSE]);
	return(p)
}
as.character.pm = function(p, leading=NA, do.sort=TRUE, do.rev=FALSE, sort.order=TRUE, simplify.complex=TRUE) {
	### Var order
	isNA = all(is.na(leading));
	if( ! isNA && ! is.numeric(leading)) leading = match(leading, names(p));
	if(any(is.na(leading))) {
		if( ! isNA) warning("Sort var does NOT exist!");
		leading = leading[ ! is.na(leading)];
	}
	if(length(leading) > 0) {
		if(do.sort) {
			p = sort.simple.pm(p, leading=leading, do.rev=do.rev);
		} else {
			if(do.rev) leading = rev(leading);
			p = cbind(p[,-leading, drop=FALSE], p[,leading, drop=FALSE]);
		}
	}
	###
	id.coeff = match("coeff", colnames(p));
	coeff = p[,id.coeff]; p = p[, - id.coeff, drop=FALSE];
	p.str = sapply(colnames(p), print.monome, p=p);
	# print(p.str)
	paste.nonempty = function(str, collapse="*") {
		str = str[nchar(str) > 0]
		paste(str, collapse=collapse)
	}
	if( ! is.null(dim(p.str))) p.str = apply(p.str, 1, paste.nonempty)
	else p.str = paste.nonempty(p.str);
	# Sign: 0 treated as "+"
	isPlus = (Re(coeff) > 0) | (Re(coeff) == 0 & Im(coeff) >= 0);
	sign.str = ifelse(isPlus, " + ", " - ");
	sign.str[1] = if(isPlus[1]) "" else "- ";
	# Complex numbers
	hasNegativ = (Re(coeff) < 0) | (Re(coeff) == 0 & Im(coeff) < 0);
	coeffPlus = coeff; coeffPlus[hasNegativ] = - coeffPlus[hasNegativ];
	coeff.str = as.character(coeffPlus);
	if(simplify.complex) {
		isPureImaginary = (Re(coeff) == 0) & (Im(coeff) != 0);
		coeff.str[isPureImaginary] = paste0(as.character(Im(coeffPlus[isPureImaginary])), "i");
	}
	# Coeff == 1
	hasCoeff = (coeffPlus != 1 & nchar(p.str) > 0); # TODO: verify if fixed!
	isB0 = (nchar(p.str) == 0);
	p.str[hasCoeff] = paste(coeff.str[hasCoeff], p.str[hasCoeff], sep = "*");
	p.str[isB0] = coeff.str[isB0]; # [should be fixed]: ERROR "+ b0*";
	return(paste(sign.str, p.str, sep="", collapse=""));
}
toCoeff = function(p, x="x") {
	idx = match(x, names(p));
	if(idx < 0) stop(paste0("No variable ", x));
	px = p[,x]; p = p[, - idx, drop=FALSE];
	str = tapply(seq(nrow(p)), px, function(nr) print.p(p[nr,, drop=FALSE], leading=NA))
	str[nchar(str) == 0] = "1";
	# missing powers
	x.all = seq(0, max(px));
	p.all = rep("0", length(x.all));
	p.all[1 + sort(unique(px))] = str;
	return(p.all)
}
evalCoeff = function(p, x="x", ...) {
	idx = match(x, names(p));
	if(idx < 0) stop(paste0("No variable ", x));
	px = p[,x]; p = p[, - idx, drop=FALSE];
	coeff = tapply(seq(nrow(p)), px, function(nr) eval.pm(p[nr,, drop=FALSE], ...));
	# missing powers
	x.all = seq(0, max(px));
	p.all = rep(0, length(x.all));
	p.all[1 + sort(unique(px))] = coeff;
	p.all = rev(p.all);
	return(p.all);
}
print.coeff = function(p, x="x") {
	p = rev(toCoeff(p, x));
	last = tail(p, 1);
	sapply(head(p, -1), function(p) cat(paste(p, ",\n", sep="")));
	cat(paste(last, "\n", sep=""));
	invisible(p);
}
print.pcoeff = function(l, print=TRUE, strip=NULL, len=10) {
	nlast = length(l);
	lsep = rep(", ", nlast);
	lsep[nlast] = ""; # tail(lsep, 1) = ""; # DOES NOT function!
	lsep[seq(len+1, nlast-1, by=len)] = ",\n";
	l.str = paste0(l, lsep, collapse="");
	if( ! is.null(strip)) {
		l.str = gsub(paste0("[ *]*+", strip, "\\^[0-9]++"), "", l.str, perl=TRUE);
	}
	if(print) { cat(l.str); cat("\n"); }
	return(invisible(l.str));
}


### Polynomial Generators:
# - moved Generators to file:
#   Polynomials.Helper.Generators.R;


#######################
#######################

#############
### Tests ###
#############

# - moved to file:
#   Polynomials.Helper.Tests.R;

