########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Elementary Polynomials
###
### draft v.0.2a


#######################


### Helper Functions

### Polynomial Tools
source("Polynomials.Helper.R")

### this file:
# source("Polynomials.Helper.EP.R")


########################
########################

### Processing Symmetric Polynomials

### Extract Base-Terms of Symmetric Polynomials:
unique.rpm =  function(p, xn="r", v=5) {
	vn = if(length(xn) > 1) xn else paste0(xn, seq(v));
	r = p[, vn, drop=FALSE];
	pr = t(sapply(seq(nrow(p)), function(id) sort(r[id,], decreasing=TRUE)));
	p = cbind(pr, p[, - match(vn, names(p)) ]);
	p = unique(p); rownames(p) = seq(nrow(p));
	return(p);
}
filter.byCol = function(E, idCol, val=0) {
	# returns ids of rows
	apply(E[ , idCol, drop=FALSE] == val, 1, all)
}
### Decompose Symmetric Polynomial
coef.rpm = function(p, flt=c("S", "E2", "E3"), rpl=NULL, pow=NULL, xn="x", sort=TRUE, debug=TRUE) {
	xpow = if(is.null(pow)) sort(unique(p$x), decreasing=FALSE) else pow;
	l = lapply(xpow, function(pow) {
		if(debug) print(paste0("Processing power = ", pow));
		pp = replace.rpm(p[p$x == pow, ], flt=flt);
		if( ! is.null(rpl)) {
			for(i in seq(length(rpl))) {
				pp = replace.pm(pp, rpl[[i]], names(rpl)[i]);
			}
		}
		return(pp);
	});
	pp = sum.lpm(l);
	if(sort) pp = sort.pm(pp, c(4,3), xn);
	rownames(pp) = seq(nrow(pp));
	return(pp);
}
replace.rpm = function(p, xn="r", v=5, flt=c("S", "E2", "E3")) {
	# extract Base-Terms: does NOT check that poly is Symmetric!!!
	p = unique.rpm(p, xn=xn, v=v);
	vn = if(length(xn) > 1) xn else paste0(xn, seq(v));
	pr = p[, vn];
	ps = p[, - match(vn, names(p))];
	r = lapply(seq(nrow(p)), function(id) {
		pow = unlist(pr[id, pr[id,] > 0]);
		if(is.null(pow)) return(ps[id,]);
		E = Epoly.distinct(pow, v);
		# filter: S = 0; E2 = 0; E3 = 0
		idCol = match(flt, names(E));
		idCol = idCol[ ! is.na(idCol)];
		if(length(idCol) > 0) E = E[filter.byCol(E, idCol, val=0), , drop=FALSE];
		if(nrow(E) == 0) return(data.frame());
		E = reduce.var.pm(E);
		return(mult.pm(ps[id,], E));
	})
	r = r[sapply(r, function(x) nrow(x) > 0)]
	return(sum.lpm(r));
}
sort.rpm = function(p, xn="r", v=5, sortPowers=TRUE) {
	# sortPowers = sort individual Powers (partially);
	vn = if(length(xn) > 1) xn else paste0("r", seq(v, 1));
	pr.max = sapply(seq(nrow(p)), function(id) max(p[id, vn]));
	pr.min = sapply(seq(nrow(p)), function(id) {
		mP = p[id, vn];
		min(mP[mP > 0]); # extract only non-0;
	});
	pr.nmax = sapply(seq(nrow(p)), function(id) {
		mP = p[id, vn];
		max.p = max(mP);
		sum(mP == max.p); # how many vars have max power?
	});
	pr.id = if( ! sortPowers) 0
		else sapply(seq(nrow(p)), function(id) {
			mP = p[id, vn]; max.p = max(mP);
			sum(which(mP == max.p));
		});
	id = if(sortPowers) order(-pr.max, -pr.nmax, -pr.min, - pr.id) else order(-pr.max, -pr.nmax, -pr.min);
	p = p[id,];
	rownames(p) = seq(nrow(p));
	return(p)
}

### Poly Generators
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
### ???
sym.poly = function(p, var="x") {
	n = length(p)
	l = list(p);
	l = rep(l, n);
	pP = 0; # TODO!
	names(pP) = paste0(var, seq(n));
	pP$coeff = 1;
	return(pP);
}

########################

### E2: Hetero-Symmetric

### Sum & Diff:
E2.pm = function(n) {
	epow = powAll.pm(toPoly.pm("S^2 - 4*E2"), n %/% 2, asList=TRUE);
	pEDiff = lapply(seq(n), function(n) diff.E2.pm(n, epow=epow));
	pESum  = lapply(seq(n), function(n) sum.E2.pm(n, epow=epow));
	return(list(pEDiff=pEDiff, pESum=pESum));
}
### Diff: x^n - y^n
diff.E2.pm = function(n, epow=NULL) {
	if(n == 1) return(data.frame(S=0, D=1, coeff=1));
	if(n == 2) return(data.frame(S=1, D=1, coeff=1));
	if(n == 0) return(0); # as.data.frame?
	# n >= 3
	np = (n-1) %/% 2;
	if(is.null(epow)) {
		epow = powAll.pm(toPoly.pm("S^2 - 4*E2"), np, asList=TRUE);
	}
	SF = function(n, c) data.frame(S=n, E2=0, coeff=c);
	# Result:
	r = SF(n-1, n);
	for(id in seq(1, np)) {
		tmp = epow[[id]];
		tmp = mult.pm(tmp, SF(n - 1 - 2*id, choose(n, 2*id+1)));
		r = sum.pm(r, tmp);
	}
	r$D = 1;
	r$coeff = r$coeff / 2^(n-1);
	return(r);
}
### Sum: x^n + y^n
sum.E2.pm = function(n, epow=NULL) {
	if(n == 1) return(data.frame(S=1, coeff=1));
	if(n == 2) return(data.frame(S=c(2,0), E2=c(0,1), coeff=c(1,-2)));
	if(n == 0) return(2); # as.data.frame?
	# n >= 3
	np = n %/% 2;
	if(is.null(epow)) {
		epow = powAll.pm(toPoly.pm("S^2 - 4*E2"), np, asList=TRUE);
	}
	SF = function(n, c) data.frame(S=n, E2=0, coeff=c);
	# Result:
	r = SF(n, 1);
	for(id in seq(1, np)) {
		tmp = epow[[id]];
		tmp = mult.pm(tmp, SF(n - 2*id, choose(n, 2*id)));
		r = sum.pm(r, tmp);
	}
	r$coeff = r$coeff / 2^(n-1);
	return(r);
}

### Permutations

# TODO: use this function;
perm = function(x, unique=FALSE, sort=TRUE) {
	# adapted from package DescTools;
	# other variants:
	# https://stackoverflow.com/questions/11095992/generating-all-distinct-permutations-of-a-list-in-r
	n = length(x)
	if (n == 1L) return(matrix(x));
	z = matrix(1L);
	for (i in 2L:n) {
		y = cbind(z, i);
		a = c(1L:i, 1:(i - 1L));
		z = matrix(0L, ncol = ncol(y), nrow = i * nrow(y));
		z[seq(nrow(y)), ] = y;
		for (j in seq(i - 1L)) {
			z[j * nrow(y) + seq(nrow(y)), ] = y[, a[seq(i) + j]]
		}
	}
    dimnames(z) = NULL;
	m = apply(z, 2L, function(i) x[i]);
	if (unique && any(duplicated(x))) m = unique(m);
	if (sort) {
		order.f = function(...) order(..., decreasing=TRUE);
		id = do.call(order.f, data.frame(m));
		m = m[id,];
	}
	return(m);
}

# countDuplicates(perm(c(3,1,1,0,0)))

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
		pP = data.frame(c(pow[1],0), c(0,pow[2]));
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
### used by Epoly.distinct()
permSum.simple = function(p1, p2, n1, n2) {
	# rep(p1, n1) + permutations of(rep(p2, n2));
	if(n1 < n2) {
		tmp = p1; p1 = p2; p2 = tmp;
		tmp = n1; n1 = n2; n2 = tmp;
	}
	p0 = c(rep(p1, n1), rep(0, n2-1));
	p = p0;
	p[seq(n2)] = p1 + p2;
	if(n2 == 1) return(matrix(p, nrow=1));
	# partially overlapping
	for(npos in seq(n2 - 1)) {
		pp = p0;
		pp[seq(n2 - npos)] = p1 + p2;
		pp[n1 + seq(npos)] = p2;
		p = cbind(p, pp);
	}
	p = t(p);
	return(p);
}
permSum.all = function(p1, p2) {
	# TODO
	if(length(p1) < length(p2)) {
		tmp = p1; p1 = p2; p2 = tmp;
	}
	n2 = length(p2);
	if(length(p1) > 2 || n2 > 2) stop("Not yet implemented!")
	p0 = c(p1, rep(0, n2 - 1));
	p = p0;
	p[seq(n2)] = p[seq(n2)] + p2;
	#
	pp = p0;
	pp[seq(n2)] = pp[seq(n2)] + p2[c(2,1)];
	p = cbind(p, pp);
	#
	pp = p0;
	pp[1] = pp[1] + p2[1]; pp[3] = p2[2];
	p = cbind(p, pp);
	pp = p0;
	pp[2] = pp[2] + p2[1]; pp[3] = p2[2];
	p = cbind(p, pp);
	#
	pp = p0;
	pp[1] = pp[1] + p2[2]; pp[3] = p2[1];
	p = cbind(p, pp);
	pp = p0;
	pp[2] = pp[2] + p2[2]; pp[3] = p2[1];
	p = cbind(p, pp);
	p = unique(t(p));
	isD2 = (p[,1] == p[,2]) & (p[,3] > 0);
	pD2 = p[isD2,]; p = rbind(p, pD2);
	return(p);
}

### Elementary Polynomials

### Basic EP
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
### Derived EP
### Sum(x[i]^n)
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
		if(length(pow) == 0) { p = data.frame(coeff=1); }
		else p = Epoly.distinct(pow, v=v, E=E);
		Ep.nm = paste0("E", v);
		id = match(Ep.nm, names(p));
		if(is.na(id)) {
			p[, Ep.nm] = p.min;
		} else p[, id] = p[, id] + p.min;
		if(full) return(list(E=E, p=p));
		return(p);
	}
	if(length(pow) == 1) return(Epoly.gen(pow[1], v=v, e=1, E=E, full=full));
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
		id = which(tbl == tbl.max)[1];
		pow.unq = as.integer(names(tbl));
		last.pow = pow.unq[id];
		isLastDuplic = (pow == last.pow);
		p  = Epoly.distinct(pow[ ! isLastDuplic], v=v, E=E);
		p2 = Epoly.distinct(pow[isLastDuplic], v=v, E=E);
		p = mult.pm(p2, p);
		if(length(pow[ ! isLastDuplic]) == 1) {
			cmb.pow = pow[isLastDuplic]; cmb.pow[1] = cmb.pow[1] + pow[ ! isLastDuplic];
			p = diff.pm(p, Epoly.distinct(cmb.pow, v=v, E=E));
		} else if(length(tbl) == 2) {
			len1 = sum(isLastDuplic); len2 = length(pow) - len1;
			pp = permSum.simple(last.pow, pow[ ! isLastDuplic][1], n1=len1, n2=len2);
			for(nr in seq(nrow(pp))) {
				p = diff.pm(p, Epoly.distinct(pp[nr,], v=v, E=E));
			}
		} else {
			pp = permSum.all(pow[ ! isLastDuplic], pow[isLastDuplic]);
			for(nr in seq(nrow(pp))) {
				p = diff.pm(p, Epoly.distinct(pp[nr,], v=v, E=E));
			}
		}
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

