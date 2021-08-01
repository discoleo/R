########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions


#######################

library(polynom)
library(pracma)

### helper Functions

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

### Root
rootn = function(r, n) {
	if(n %% 2 == 0) return(r^(1/n)); # TODO: complex?
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

sort.sol = function(sol, useRe=TRUE, ncol=1) {
	if(useRe) {
		id = order(abs(sol[,ncol]), Re(sol[,ncol]));
	} else {
		id = order(abs(sol[,ncol]));
	}
	return(sol[id,]);
}
isConj.f = function(x, y, tol=1E-3) {
	isConj = (abs(Re(x) - Re(y)) < tol) & (abs(Im(x) + Im(y)) < tol);
	return(isConj);
}

### Polynomials

### Multiplication
mult.p = function(p1, p2) {
	p.m = outer(p1, p2)
    p = as.vector(tapply(p.m, row(p.m) + col(p.m), sum))
	return(p)
}

### Multi-variable Polynomials:
# p = multi-variable polynomial;
#  => data.frame with exponents of variables;
#  => coeff = column with the coefficients;
# Note:
# - initial idea was to allow also basic lists,
#   but most functions work only with data frames!
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
	if(is.numeric(p1)) {
		if(is.list(p2)) return(mult.sc.pm(p2, p1));
		return(data.frame(coeff=p1*p2));
	}
	if(is.data.frame(p1)) {
		if(ncol(p1) == 1) {
			if(is.numeric(p2)) return(data.frame(coeff=p1$coeff * p2));
			if(ncol(p2) == 1) return(data.frame(coeff=p1$coeff * p2$coeff));
			return(mult.sc.pm(p2, p1$coeff));
		}
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
pow.pm = function(p, n=2) {
	if(n == 1) return(p);
	if(is.double(n) && (n == round(n))) n = as.integer(n);
	if( ! is.integer(n)) stop("n must be integer!")
	# Multiply
	# TODO: vectorize: as.integer(intToBits(3));
	p.r = NULL;
	p.pow = p;
	while (n > 0) {
		print(n)
		if (n %% 2 == 1) {
			if(is.null(p.r)) p.r = p.pow else p.r = mult.pm(p.r, p.pow);
		}
		if(n == 1) break;
        p.pow = mult.pm(p.pow);
        n = n %/% 2;
    }
	x.name = names(p)[1];
	id = order(p.r[,x.name], decreasing=TRUE);
	p.r = p.r[id,];
	return(p.r);
}
mult.sc.pm = function(p, s, div=1, coeff.name="coeff") {
	# Multiplication by scalar
	# div = for numerical stability;
	if(is.data.frame(p)) {
		p[ , coeff.name] = p[ , coeff.name] * s;
		if(div != 1) p[ , coeff.name] = p[ , coeff.name] / div;
	} else if(is.list(p)) {
		p[[coeff.name]] = p[[coeff.name]] * s;
		if(div != 1) p[[coeff.name]] = p[[coeff.name]] / div;
	} else stop("p must be a polynomial!")
	return(p);
}
simplify.spm = function(p1) {
	nms = names(p1);
	nms = nms[ ! nms %in% "coeff"];
	for(nm in nms) {
		v.pow = min(p1[,nm]);
		if(v.pow > 0) {
			p1[,nm] = p1[,nm] - v.pow;
		}
	}
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
### Helper functions
align.pm = function(p1, p2, align.names=TRUE) {
	# align columns of 2 data.frames for sum.pm();
	p1 = reduce.pm(p1); p2 = reduce.pm(p2);
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
sum.pm = function(p1, p2) {
	if(is.data.frame(p1) && nrow(p1) == 0) return(reduce.pm(p2));
	if(is.data.frame(p2) && nrow(p2) == 0) return(reduce.pm(p1));
	l = align.pm(p1, p2);
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
replace.withVal.pm = function(p, x, pow=1, val, simplify=TRUE) {
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
	if(is.na(idx)) {
		warning(paste0("Polynomial does NOT contain variable: ", x));
		return(p1);
	}
	if(is.numeric(p2)) return(replace.withVal.pm(p1, x=x, pow=pow, val=p2));
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
	pR = data.frame();
	for(nr in seq(nrow(p1))) {
		ipow = rpow[nr];
		lp = if(ipow == 0) list(p2fr.pows[[max.pow]])
			else if(max.pow == ipow) list(p2.pows[[max.pow]])
			else list(p2.pows[[ipow]], p2fr.pows[[max.pow - ipow]]);
		lp = c(lp, list(p1[nr,]));
		tmp = mult.all.pm(lp);
		pR = sum.pm(pR, tmp);
	};
	return(reduce.var.pm(pR));
}
sort.pm = function(p, sort.coeff=1, xn=NULL) {
	pP = p[, - which(names(p) == "coeff"), drop=FALSE];
	pow.tot = sapply(seq(nrow(p)), function(id) sum(pP[id, ]));
	pow.max = sapply(seq(nrow(p)), function(id) max(pP[id, ]));
	if(length(sort.coeff) == 1) {
		id = order(abs(p$coeff), pow.tot, pow.max);
	} else {
		coeff.df = data.frame(abs(p$coeff), -pow.tot, -pow.max);
		if( ! is.null(xn)) coeff.df = cbind(coeff.df, -pP[,xn]);
		if(length(sort.coeff) > 3 + length(xn)) {
			# minimum power of Monome: may be 0;
			coeff.df$min = sapply(seq(nrow(p)), function(id) -min(pP[id, ]));
		}
		coeff.df = coeff.df[, sort.coeff];
		id = do.call(order, coeff.df)
	}
	return(p[id,])
}
eval.pm = function(p, x, progress=FALSE) {
	pP = p[, - which(names(p) == "coeff")];
	eval.p = function(id) {
		idx = which(pP[id,] != 0);
		prod(x[idx]^pP[id, idx], p$coeff[id]);
	}
	sum(sapply(seq(nrow(p)), eval.p))
}
eval.cpm = function(p, x, bits=120, tol=1E-12, progress=FALSE) {
	# uses the Rmpfr package;
	# currently assumes that only coeffs are big and
	# have an impact on numeric stability;
	pP = p[, - which(names(p) == "coeff")];
	# pow = lapply(seq(ncol(pP)), function(nc) sort(unique(pP[,nc])));
	# currently only max:
	pow = lapply(seq(ncol(pP)), function(nc) max(pP[,nc]));
	xpows = lapply(seq(length(x)), function(id) {
		x0 = x[id];
		len = tail(pow[[id]], 1);
		if(is.complex(x0) && Im(x0) != 0) {
			# TODO: explore polar coordinates;
			x = x0^seq(len);
			re = Re(x); im = Im(x);
			div = 1;
			re = mpfr(re * div, bits);
			im = mpfr(im * div, bits);
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
		if(progress && id %% 16 == 1) cat(paste0(id, if(id %% 64 == 1) "\n" else ", "));
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
	sol = sapply(seq(nrow(p)), eval.p);
	sdim = attr(sol, "dim"); sol = mpfr2array(t(sol), rev(sdim));
	sol = apply(sol, 2, sum);
	return(sol);
}
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
gcd.vpm = function(p, xgcd=0) {
	for(i in seq(nrow(p))) xgcd = gcd(xgcd, p$coeff[i]);
	return(xgcd);
}
gcd.pm = function(p1, p2, by="x", div.sc=1) {
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
dp.pm = function(p, xn="x") {
	p = p[(p[,xn] != 0),];
	p$coeff = p$coeff * p[,xn];
	p[,xn] = p[,xn] - 1;
	return(p);
}


#############
### Other ###

perm.gen = function(x) {
	len = length(x)
	id = seq(len)
	id.m = outer(id, id, function(i, j) ((i+j+1) %% len + 1))
	p.m = x[id.m]
	dim(p.m) = dim(id.m)
	p.m
}

### Solvers:

### Simple systems:
solve.En = function(x, max.perm=0, n=4, duplicates=FALSE) {
	id = 1:length(x)
	if(max.perm == 1) {
		id.gr = matrix(id, nrow=1)
	} else {
		id.l = rep(list(id), n)
		id.gr = expand.grid(id.l)
		if( ! duplicates) {
			isDuplic = apply(id.gr, 1, function(id.val) any(duplicated(id.val)))
			id.gr = id.gr[ ! isDuplic , ]
		}
		if(max.perm > 0) {
			max.perm = min(max.perm, nrow(id.gr));
			id.gr = head(id.gr, n=max.perm);
		}
	}
	sol = if(n == 4) cbind(
			x1=x[id.gr[,1]], x2=x[id.gr[,2]], x3=x[id.gr[,3]], x4=x[id.gr[,4]])
		else cbind(x[id.gr[,1]], x[id.gr[,2]], x[id.gr[,3]])
	sol.names = if(n >= 4) paste0("x", seq(n)) else
		if(n == 3) c("x", "y", "z") else c("x", "y");
	colnames(sol) = sol.names;
	return(sol);
}
solve.EnAll = function(m, max.perm=0, n=4) {
	# generates ncol(m) * (nrow(m)!) root combinations/permutations!
	l = lapply(seq(ncol(m)), function(id) solve.En(as.vector(m[,id]), max.perm=max.perm, n=n));
	do.call(rbind, l)
}

### decomposed polynomial systems
solve.S = function(S, R, b=0) {
	# generic solver (based on existing S = x+y+z)
	b2 = if(length(b) > 1) b[2] else 0; # Ext A2;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A3;
	x = sapply(S, function(x) roots(c(1, -x, R[2] - b2*x, - R[3] + b3*x)))
	len = length(S)
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	yz = R[3]/x - b3
	yz.s = S - x
	# TODO: robust (when necessary)
	# Note: this simple method is NOT robust!
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

solve.mS = function(S, b=0) {
	# generic solver (based on existing S = x+y+z)
	# S = cbind(S, E2, E3)
	b2 = if(length(b) > 1) b[2] else 0; # Ext A2;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A3;
	x = sapply(seq(nrow(S)), function(id) roots(c(1, -S[id, 1], S[id, 2] - b2*S[id, 1], - S[id, 3] + b3*S[id, 1])))
	E2 = S[,2]; E3 = S[,3]; S = S[,1]; len = length(S)
	S  = matrix(S, ncol=len, nrow=3, byrow=T)
	E3 = matrix(E3, ncol=len, nrow=3, byrow=T)
	yz = E3/x - b3
	yz.s = S - x
	# TODO: robust (when necessary)
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

#########

### Print

# Print multi-variable Poly
print.monome = function(name, p) {
	v = p[,name];
	v.r = rep("", length(v));
	v.r[v > 1] = paste0(name, "^", v[v > 1]);
	v.r[v == 1] = name;
	return(v.r);
}
print.p = function(p, leading=1, order=TRUE, sort.order=TRUE) {
	### Var order
	if( ! is.numeric(leading)) leading = match(leading, names(p));
	if( ! is.na(leading)) {
		if(order) p = p[order(p[, leading], decreasing=sort.order), ];
		p = cbind(p[,-leading, drop=FALSE], p[,leading, drop=FALSE]);
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
	sign.str = ifelse(coeff > 0, " + ", " - ");
	sign.str[1] = if(coeff[1] > 0) "" else "- ";
	coeff.str = as.character(abs(coeff));
	hasCoeff = (abs(coeff) != 1); # TODO: ERROR "+ b0*";
	p.str[hasCoeff] = paste(coeff.str[hasCoeff], p.str[hasCoeff], sep = "*");
	return(paste(sign.str, p.str, sep="", collapse=""));
}
toCoeff = function(p, x="x") {
	idx = match(x, names(p));
	if(idx < 0) stop(paste0("No variable ", x));
	px = p[,x]; p = p[, - idx, drop=FALSE];
	str = tapply(seq(nrow(p)), px, function(nr) print.p(p[nr,], leading=NA))
	str[nchar(str) == 0] = "1";
	# missing powers
	x.all = seq(0, max(px));
	p.all = rep("0", length(x.all));
	p.all[1 + sort(unique(px))] = str;
	return(p.all)
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
### Parse expressions / polynomials
toPoly.pm = function(e) {
	if(is.character(e)) {
		if(length(e) > 1) {
			pl = lapply(e, function(e) toPoly.pm(e));
			return(pl);
		}
		e = parse(text=e);
	}
	if( ! is.expression(e)) stop("Input must be an expression!")
	if( ! is.language(e[[1]])) return(NULL);
	e = e[[1]];
	p = data.frame();
	while(TRUE) {
		isSymbol = is.symbol(e);
		if(isSymbol || is.symbol(e[[1]])) {
			op = if(isSymbol) e else e[[1]];
			if(op == "+") {
				m = toMonom.pm(e[[3]]);
				p = if(nrow(p) == 0) m else sum.pm(p, m);
				e = e[[2]];
			} else if(op == "-") {
				if(length(e) > 2) {
					m = toMonom.pm(e[[3]], xsign=-1);
					p = if(nrow(p) == 0) m else sum.pm(p, m);
					e = e[[2]];
				} else {
					m = toMonom.pm(e[[2]], xsign=-1);
					p = if(nrow(p) == 0) m else sum.pm(p, m);
					break;
				}
			} else {
				m = toMonom.pm(e);
				p = if(nrow(p) == 0) m else sum.pm(p, m);
				break;
			}
		} else break;
	}
	return(p);
}
### Parse expression
parse.pm = function(e) {
	if( ! is.expression(e)) stop("Input must be an expression!")
	if( ! is.language(e[[1]])) return(NULL);
	e = e[[1]];
	e.txt = character(0);
	c.e = function(e, x.sign) {
		xi = if(nchar(x.sign) == 0) format(e[[3]]) else paste(x.sign, format(e[[3]]));
		c(e.txt, xi);
	}
	while(TRUE) {
		if(is.symbol(e[[1]])) {
			x.sign = paste0(e[[1]]);
			if(x.sign == "+") x.sign = ""
			else if(x.sign != "-") {
				e.txt = c(e.txt, format(e));
				break;
			}
		} else {
			print(e); break;
		}
		e.txt = c.e(e, x.sign);
		if(is.language(e[[2]])) e = e[[2]]
		else break;
		
	}
	return(e.txt);
}

toMonom.pm = function(e, xsign = 1) {
	m = data.frame(coeff=xsign);
	acc = list();
	while(TRUE) {
		if(length(e) == 1) {
			if(is.symbol(e)) {
				if(e == "-") {
					m$coeff = - m$coeff; # NO effect?
				} else if(e == "+") {
					# an extra "+"; # NO effect?
				} else {
					vn1 = as.character(e); # a variable name;
					m[, vn1] = 1;
				}
			} else if(is.numeric(e)) {
				m[, "coeff"] = m[, "coeff"] * e;
			} else print(paste0("Error: ", e));
		} else {
			op = e[[1]];
			if(is.symbol(op)) {
				if(op == "*") {
					acc = c(acc, e[[2]]);
					e = e[[3]]; next;
				}
				if(op == "^") {
					vn1 = as.character(e[[2]]); # TODO: 8^8
					pow = e[[3]];
					if( ! is.numeric(pow)) {
						warning(paste0("Power = ", pow, " is NOT numeric!"));
						pow = NA;
					}
					m[, vn1] = pow;
				} else if(op == "-") {
					m$coeff = - m$coeff;
					e = e[[2]]; next;
				} else if(op == "+") {
					e = e[[2]]; next;
				} else if(op == "/") {
					m[, "coeff"] = m[, "coeff"] / e[[3]];
					e = e[[2]]; next;
				} else {
					vn1 = as.character(op); # a variable name;
					m[, vn1] = 1;
				}
			} else if(is.numeric(op)) {
				m[, "coeff"] = m[, "coeff"] * op;
			}
		}
		if(length(acc) == 0) break;
		e = acc[[length(acc)]];
		acc = head(acc, -1);
	}
	return(m);
}

### Classic Polynomials
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

### Derived Polynomials
roots.derived = function(n, pow=seq(n-1), rn="r", sn="s", all.roots=TRUE) {
	slen = length(pow);
	S = diag(slen);
	s = lapply(seq(nrow(S)), function(nr) S[nr,]);
	s = data.frame(s); names(s) = paste0(sn, pow);
	p1 = data.frame(x=rep(0, slen), r1=pow, s, coeff=-1);
	names(p1)[2] = paste0(rn, 1);
	p1 = rbind(p1, c(1, rep(0, slen+1), 1))
	if(all.roots) {
		p.list = list(p1);
		for(id in seq(2, n)) {
			pT = p1;
			names(pT)[2] = paste0(rn, id);
			p.list[[id]] = pT;
		}
		p1 = p.list;
	}
	return(p1)
}

### Extensions
extend.spm = function(p, n=2, vb="be", vR="R", vS="S", sort=TRUE) {
	pS = data.frame(R=0, S=seq(n), coeff=-1);
	pS = rbind(pS, c(1,0,1));
	names(pS)[1:2] = c(vR, vS);
	b.all = paste0(vb, seq(n));
	for(id in seq(n)) {
		pB = data.frame(rep(0, n+1));
		pB[id,1] = id; names(pB) = b.all[id];
		pS = cbind(pS, pB);
	}
	# substitute in p
	p = replace.pm(p, pS, vR, pow=1);
	if(sort) p = sort.pm(p, c(4,3), xn=vS);
	return(p);
}

#######################
#######################

#############
### Tests ###
#############

### Multi-variable Multiplication

# (x^3 + b1*x - R)^3
p = list(
	x = c(3,1,0),
	b1 = c(0,1,0),
	R = c(0,0,1),
	coeff = c(1,1,-1)
)

### Test
mult.pm(p)

p.v = pow.pm(p, 3)
p.v

print.p(p.v[,c(2,3,4,1)])

