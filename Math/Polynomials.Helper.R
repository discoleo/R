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
### Multi-variable Multiplication
# p = multi-variable polynomial;
#  => data.frame with exponents of variables;
#  => coeff = column with the coefficients;
# Note:
# - initial idea was to allow also basic lists,
#   but most functions work only with data frames!
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
	p.b0 = prod.b0(p1.b0, p2.b0);
	p.l = lapply(p.m, as.vector);
	p.v = do.call(cbind, p.l)
	p.v = cbind(p.v, b0=as.vector(p.b0));
	p.r = aggregate(b0~., p.v, sum);
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
reduce.pm = function(p) {
	# remove Monomials with coeff == 0;
	id = which(p$coeff != 0)
	if(is.data.frame(p)) {
		return(p[id, ]);
	}
	p = lapply(p, function(m) m[id]);
	return(p);
}
reduce.var.pm = function(p) {
	# remove Vars with power == 0;
	id = match("coeff", names(p));
	nc = rep(TRUE, ncol(p));
	nc[-id] = sapply(seq(ncol(p))[-id], function(id) any(p[,id] != 0));
	return(p[, nc]);
}
align.pm = function(p1, p2, align.names=TRUE) {
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
	p.r = aggregate(coeff~., p, sum);
	return(reduce.pm(p.r));
}
sum.lpm = function(lp) {
	pR = data.frame();
	for(pd in lp) {
		pR = add.pm(pR, pd);
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
replace.pm = function(p1, p2, x, pow=1) {
	# replace x^pow by p2;
	idx = match(x, names(p1));
	if(is.na(idx)) stop(paste0("Polynomial does NOT contain variable: ", x));
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
			p2.pows[[ipow]] = mult.pm(p2.pows[[ipow - 1]], p2);
			p2fr.pows[[ipow]] = mult.pm(p2fr.pows[[ipow - 1]], p2fr);
		}
	}
	pR = data.frame();
	for(nr in seq(nrow(p1))) {
		ipow = rpow[nr];
		lp = if(ipow == 0) list(p2fr.pows[[max.pow]])
			else if(max.pow == ipow) list(p2.pows[[max.pow]])
			else list(p2.pows[[ipow]], p2fr.pows[[max.pow - ipow]]);
		lp = c(lp, list(p1[nr,]));
		tmp = mult.all.pm(lp);
		pR = add.pm(pR, tmp);
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
eval.pm = function(p, x) {
	pP = p[, - which(names(p) == "coeff")];
	eval.p = function(id) {
		idx = which(pP[id,] != 0);
		prod(x[idx]^pP[id, idx], p$coeff[id]);
	}
	sum(sapply(seq(nrow(p)), eval.p))
}
div.pm = function(p1, p2, by="x", debug=TRUE) {
	# very simple division
	xn = by;
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
	if(nrow(pDx) == 1) {
		idn = match(names(pDx)[-idc2], names(p1));
		print(idn);
		if(any(is.na(idn))) stop(paste0("No matching variables: ", names(pDx)[is.na(idn)]));
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
	}
	if(debug) {
		if(nrow(p1) > 0) print("Not divisible!")
		else print("Divisible!");
	}
	return(list(Rez=pRez, Rem=p1));
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
solve.pm = function(p1, p2, xn, stop.at=NULL, simplify=TRUE) {
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
		lp2[[2]]$coeff = - lp2[[2]]$coeff; # !
		p1 = replace.fr.bigpm(p1, lp2[[1]], lp2[[2]], x=xn, pow=1);
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
	p1 = sum.bigpm(mult.bigpm(p1, p2cf), mult.bigpm(p2, p1cf, -1));
	print(paste0("Max pow: ", max1, "; Len = ", nrow(p1)));
	if(simplify) { p1 = simplify.spm(p1); p2 = simplify.spm(p2); }
	return(solve.pm(p1, p2, xn=xn, stop.at=stop.at));
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
	px = p[,x]; p = p[, - idx];
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
	sapply(p, function(p) cat(paste(p, ",\n", sep="")));
	invisible(p);
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
### Poly Calculations
rotate = function(p, n, val0=0, asPoly=TRUE) {
	m = matrix(val0, nrow=n, ncol=n);
	len = length(p);
	max.len = n - len + 1;
	for(nc in seq(n)) {
		if(nc <= max.len) {
			m[seq(nc, length.out=len), nc] = p;
		} else {
			dn = nc - max.len;
			m[seq(1, dn), nc] = p[seq(len - dn + 1, len)];
			m[seq(nc, n), nc] = p[seq(1, len - dn)];
		}
	}
	m = t(m);
	if(asPoly) {
		m = data.frame(m);
		names(m) = paste0("x", seq(n));
		m$coeff = 1;
	}
	return(m)
}
perm1 = function(n, p=1, val0=0) {
	m = matrix(val0, nrow=n, ncol=n);
	diag(m) = p;
	return(m);
}
perm2 = function(n, p=c(1,1), val0=0) {
	# all 2 permutations
	n1 = n - 1;
	len = n*n1/2;
	m = matrix(val0, nrow=len, ncol=n);
	ioff = 0;
	for(i in seq(1, n1)) {
		m[seq(ioff+i, ioff+n1), i] = p[1];
		ioff = ioff + n1 - i;
	}
	nC = unlist(sapply(seq(2, n), function(id) seq(id, n)));
	for(i in seq(1, len)) {
		m[i, nC[i]] = p[2];
	}
	return(m)
}
perm3 = function(n, p=c(1,1,1), val0=0) {
	if(min(p) == max(p)) {
		if(length(p) == n-1) {
			# works for: length(p) = n-1;
			return(perm1(n, p=val0, val0=p[1]));
		} else if(length(p) == n-2) {
			return(perm2(n, p=c(val0, val0), val0=p[1]));
		} else stop("Not yet implemented!")
	}
	if(any(p[3] == p[1:2])) {
		pP = perm2(n, p=c(p[3], p[3]), val0=val0);
		p3 = p[p != p[3]];
	} else {
		pP = perm2(n, p=p[1:2], val0=val0);
		if(p[1] != p[2]) pP = rbind(pP, array(rev(pP), dim(pP)));
		p3 = p[3];
	}
	{
		pP = t(pP);
		id = sapply(seq(ncol(pP)), function(id) which(pP[,id] == val0));
		# TODO: c(T,F) => c(T, rep(F, ...));
		putVar = function(pos=c(T,F)) {
			p21 = pP;
			id1 = id[rep(pos, length(id) %/% 2)];
			p21[seq(0, ncol(pP)-1)*n + id1] = p3;
			return(p21);
		}
		pP = cbind(putVar(c(T,F)), putVar(c(F,T)));
		pP = t(pP);
	}
	return(pP);
}
prod.perm.poly = function(n, pow=c(1,1)) {
	# all 2 permutations
	m = perm2(n, p=pow);
	xn = paste0("x", seq(n));
	toPoly = function(mr) {
		id = which(mr != 0);
		pP = list(c(pow[1],0), c(0,pow[2]));
		names(pP) = xn[id];
		pP$coeff = c(1,1);
		return(pP);
	}
	if(n == 2) return(toPoly(m[1,]));
	pR = mult.pm(toPoly(m[1,]), toPoly(m[2,]));
	for(id in seq(3, nrow(m))) {
		pR = mult.pm(pR, toPoly(m[id,]));
	}
	return(pR)
}
perm.poly = function(n, p=c(1,1), val0=0) {
	if(length(p) == 1) {
		m = as.data.frame(perm1(n, p=p, val0=val0))
	} else if(length(p) == 2) {
		m = as.data.frame(perm2(n, p=p, val0=val0))
	} else if(length(p) == 3) {
		m = as.data.frame(perm3(n, p=p, val0=val0))
	} else if(length(p) == 4) {
		# may work sometimes (see perm3());
		m = as.data.frame(perm3(n, p=p, val0=val0))
	}
	names(m) = paste0("x", seq(n));
	m$coeff = rep(1, nrow(m))
	return(m);
}
sym.poly = function(p, var="x") {
	n = length(p)
	l = list(p);
	l = rep(l, n);
	pP = 0; # TODO!
	names(pP) = paste0(var, seq(n));
	pP$coeff = 1;
	return(pP);
}
perm2.pm = function(p, xn) {
	# permute the variables with the set in xn;
	xn0 = names(p);
	idn = match(xn0, xn);
	idVar = which( ! is.na(idn)); # TODO: exclude coeff;
	m = perm2(length(xn));
	m = (m == 1); m = t(m);
	permute = function(nc) {
		pp = p;
		names(pp)[idVar] = xn[m[ , nc]];
		return(pp);
	}
	pAll = lapply(seq(ncol(m)), permute);
	return(pAll);
}
Eprod.pm = function(n, p=1, xn="x") {
	p = data.frame(array(p, c(1, n)));
	names(p) = if(length(xn) == n) xn else paste0(xn, seq(n));
	p$coeff = 1;
	return(p);
}
Esum.pm = function(n, p=1, xn="x") {
	p = as.data.frame(perm1(n, p=p));
	names(p) = if(length(xn) == n) xn else paste0(xn, seq(n));
	p$coeff = 1;
	return(p);
}
Epoly.base = function(n, v=4, E=NULL) {
	vn = c("S", paste0("E", seq(2, v)));
	to.df = function(v) {
		x.df = data.frame(1, 1);
		names(x.df) = c(v, "coeff");
		return(x.df);
	}
	if(is.null(E)) {
		E = list(
			data.frame(S=1, coeff=1),
			data.frame(S=c(2,0), E2=c(0,1), coeff=c(1,-2)),
			data.frame(S=c(3,1,0), E2=c(0,1,0), E3=c(0,0,1), coeff=c(1,-3,3)));
		ipow = 4;
		while(ipow <= v) {
			p = mult.pm(E[[ipow - 1]], to.df(vn[1]));
			signE = -1;
			for(i2 in 2:ipow) {
				if(ipow == i2) {
					p = add.pm(p, mult.sc.pm(to.df(vn[i2]), ipow*signE));
				} else {
					p = add.pm(p, mult.pm(E[[ipow - i2]], to.df(vn[i2]), sc=signE));
				}
				signE = -signE;
			}
			E[[ipow]] = p;
			ipow = ipow + 1;
		}
	} else {
		ipow = length(E) + 1;
	}
	# compute remaining
	while(ipow <= n) {
		p = mult.pm(E[[ipow - 1]], to.df(vn[1]));
		signE = -1;
		for(i2 in 2:v) {
			p = add.pm(p, mult.pm(E[[ipow - i2]], to.df(vn[i2]), sc=signE));
			signE = -signE;
		}
		E[[ipow]] = p;
		ipow = ipow + 1;
	}
	return(E); # TODO
}
Epoly.gen = function(n, v=4, e=1, E=NULL, full=FALSE) {
	if(e > v) return(data.frame(S=0, coeff=0));
	if(e == v) {
		p = data.frame(n, 1);
		names(p) = c(paste0("E", v), "coeff");
		return(p);
	}
	if(e == 1) {
		if(is.null(E) || length(E) < n) {
			E = Epoly.base(n, v, E=E);
		}
		if(full) return(E); # TODO: list;
		return(E[[n]]);
	}
	if(e == 2) {
		if(is.null(E) || length(E) < 2*n) {
			E = Epoly.base(2*n, v=v, E=E);
		}
		p = diff.pm(pow.pm(E[[n]], 2), E[[2*n]]);
		p = mult.sc.pm(p, 1, 2);
		if(full) return(list(E=E, p=p));
		return(p);
	}
	if(e == 3) {
		if(is.null(E) || length(E) < 3*n) {
			E = Epoly.base(3*n, v=v, E=E);
		}
		p = diff.pm(E[[3*n]], pow.pm(E[[n]], 3));
		p = mult.sc.pm(p, 1, 3);
		p = add.pm(p, mult.pm(Epoly.gen(n, v=v, e=2, E=E), E[[n]]))
		if(full) return(list(E=E, p=p));
		return(p);
	}
	# e > 3
	return(Epoly.adv(n=n, v=v, e=e, E=E, full=full));
}
Epoly.adv = function(n, v=4, e=1, E=NULL, full=FALSE) {
	if(is.null(E) || length(E) < e*n) {
		E = Epoly.base(e*n, v=v, E=E);
	}
	p = E[[e*n]];
	signE = -1;
	iPow = 1;
	while(iPow < e) {
		p = add.pm(p, mult.pm(Epoly.gen(n, v=v, e=iPow, E=E), E[[(e - iPow)*n]], sc=signE));
		iPow = iPow + 1;
		signE = -signE;
	}
	#
	p = mult.sc.pm(p, 1, -e);
	if(full) return(list(E=E, p=p));
	return(p);
}
Epoly.distinct = function(pow, v=3, E=NULL, full=FALSE) {
	pow = pow[pow != 0];
	if(length(pow) > v) stop("Error! Longer than no. of variables."); # TODO
	if(length(pow) == v) {
		# Case: Prod[over_v]
		p.min = min(pow);
		pow = pow - p.min; pow = pow[pow != 0];
		p = Epoly.distinct(pow, v=v, E=E);
		Ep.nm = paste0("E", v);
		id = match(Ep.nm, names(p));
		if(is.na(id)) {
			p[, Ep.nm] = p.min;
		} else p[, id] = p[, id] + p.min;
		if(full) return(list(E=E, p=p));
		return(p);
	}
	if(length(pow) == 1) return(Epoly.gen(p.rg[1], v=v, e=1, E=E, full=full));
	# Composite Cases:
	p.rg = range(pow);
	if(p.rg[1] == p.rg[2]) return(Epoly.gen(p.rg[1], v=v, e=length(pow), E=E, full=full));
	len = sum(pow);
	E = Epoly.base(len, v=v, E=E); # TODO: more accurate power;
	if(length(pow) == 2) {
		print(p.rg);
		### Case: p.rg[1] == p.rg[2]: already handled;
		p = mult.pm(E[[p.rg[2]]], E[[p.rg[1]]]);
		p = diff.pm(p, E[[len]]);
		if(full) return(list(E=E, p=p));
		return(p);
	}
	if(length(pow) == 3) {
		p = Epoly.distinct(pow[1:2], v=v, E=E);
		p = mult.pm(p, E[[pow[3]]]);
		sc = if(pow[1] + pow[3] != pow[2]) 1 else 2;
		p = diff.pm(p, mult.sc.pm(Epoly.distinct(c(pow[1]+pow[3], pow[2]), v=v, E=E), sc));
		if(pow[1] != pow[2]) {
			sc = if(pow[1] != (pow[2] + pow[3])) 1 else 2;
			p = diff.pm(p, mult.sc.pm(Epoly.distinct(c(pow[1], pow[2]+pow[3]), v=v, E=E), sc));
		}
		if(any(pow[1:2] == pow[3])) p = mult.sc.pm(p, 1, 2);
		if(full) return(list(E=E, p=p));
		return(p);
	}
	# TODO: generalize & check;
	tbl = table(pow);
	tbl.max = max(tbl);
	if(tbl.max > 1) {
		stop("Not yet implemented!");
		id = which(tbl == tbl.max)[1];
		pow.unq = as.integer(names(tbl));
		last.pow = pow.unq[id];
		isLastDuplic = (pow == last.pow);
		p  = Epoly.distinct(pow[ ! isLastDuplic], v=v, E=E);
		p2 = Epoly.distinct(pow[isLastDuplic], v=v, E=E);
		p = mult.pm(p, p2);
		# TODO!
	} else {
		last.pow = tail(pow, 1);
		other.pow = head(pow, -1);
		p = Epoly.distinct(head(pow, -1), v=v, E=E);
		p = mult.pm(p, E[[last.pow]]);
		for(idx in seq_along(other.pow)) {
			xpow = other.pow[idx];
			sc = if(any((xpow + last.pow) == other.pow)) 2 else 1; # can match mostly 1!
			cmb.pow = c(xpow + last.pow, other.pow[-idx]);
			p = diff.pm(p, mult.sc.pm(Epoly.distinct(cmb.pow, v=v, E=E), sc));
		}
	}
	#
	if(full) return(list(E=E, p=p));
	return(p);
}



#######################
#######################

#############
### Tests ###
#############

### E Polynomials:
### TODO:
# - thoroughly check!
Epoly.gen(5, 4)

###
x1 = sqrt(2); x2 = sqrt(3); x3 = sqrt(5); x4 = sqrt(7); x5 = sqrt(11);
xx = c(x1, x2, x3, x4, x5);
E2 = eval.pm(perm.poly(5, c(1,1)), xx);
E3 = x1*x2*(x3+x4+x5) + (x1+x2)*x3*(x4+x5) + (x1+x2+x3)*x4*x5;
E4 = prod(xx) * sum(1/xx); E5 = prod(xx); S = sum(xx);
#
eval.pm(perm.poly(5, c(3,3,3)), xx)
eval.pm(Epoly.gen(3, 5, 3), c(S,E2,E3,E4,E5))

### TODO:
# check for larger n:
Epoly.gen(4, 5, 4)
#
eval.pm(perm.poly(5, rep(4, 4)), xx)
eval.pm(Epoly.gen(4, 5, 4), c(S,E2,E3,E4,E5))

###
N = 5;
pseq.all = 0:4;
# pseq.all = c(0, 1,3,4,6);
pseq = pseq.all[pseq.all != 0]
s = sapply(1:10000, function(id) sample(pseq.all, N))
s = t(s)
s = unique(s)
nrow(s)
sum(sapply(seq(nrow(s)), function(id) prod(xx^s[id,])))
eval.pm(Epoly.distinct(pseq, N), c(S, E2, E3, E4, E5))


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

################
### Permutations
n = 4
p = prod.perm.poly(n)
p = sort.pm(p);
p

print.p(p)

apply(perm3(4, p=c(3,2,1)), 1, sum)
table(duplicated(perm3(4, p=c(3,2,1))))

