########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions


#######################

library(polynom)
library(pracma)
library(gmp) # BigNumbers

### helper Functions

# - hack to work with BigNumbers;
# - workaround for: aggregate();
# - a lot of duplicated code!


mult.all.bigpm = function(p) {
	if( ! is.list(p)) stop("p must be a list of polynomials");
	len = length(p);
	pR = p[[1]];
	for(id in seq(2, len)) {
		p2 = p[[id]];
		if(is.numeric(p2)) {
			pR = mult.sc.pm(pR, p2);
		} else {
			pR = mult.bigpm(pR, p2);
		}
	}
	return(pR);
}
mult.bigpm = function(p1, p2, sc=1) {
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
	b0 = as.vector(p.b0);
	p.v = cbind(p.v, b0=seq(nrow(p.v)));
	aggr.bignum = function(id) {
		x.df = data.frame(b0=0);
		s = sum(b0[id]);
		x.df$b0 = s;
	}
	p.r = aggregate(b0~., p.v, aggr.bignum);
	colnames(p.r) = c(vars, "coeff");
	if(sc != 1) p.r$coeff = p.r$coeff * sc;
	return(p.r);
}
sum.bigpm = function(p1, p2) {
	if(is.data.frame(p1) && nrow(p1) == 0) return(reduce.pm(p2));
	if(is.data.frame(p2) && nrow(p2) == 0) return(reduce.pm(p1));
	l = align.pm(p1, p2);
	p1 = l[[1]]; p2 = l[[2]];
	n1 = names(p1); n2 = names(p2);
	### to DF
	id = match(n2, n1);
	p = rbind(as.data.frame(p1), as.data.frame(p2)[,id]);
	### Sum
	coeff = p$coeff; p$coeff = seq(nrow(p));
	aggr.bignum = function(id) {
		x.df = data.frame(coeff=0);
		s = sum(coeff[id]);
		x.df$coeff = s;
	}
	p.r = aggregate(coeff~., p, aggr.bignum);
	return(reduce.pm(p.r));
}
replace.fr.bigpm = function(p1, p2, p2fr, x, pow=1) {
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
			p2.pows[[ipow]] = mult.bigpm(p2.pows[[ipow - 1]], p2);
			p2fr.pows[[ipow]] = mult.bigpm(p2fr.pows[[ipow - 1]], p2fr);
		}
	}
	pR = data.frame();
	for(nr in seq(nrow(p1))) {
		ipow = rpow[nr];
		lp = if(ipow == 0) list(p2fr.pows[[max.pow]])
			else if(max.pow == ipow) list(p2.pows[[max.pow]])
			else list(p2.pows[[ipow]], p2fr.pows[[max.pow - ipow]]);
		lp = c(lp, list(p1[nr,]));
		pR = sum.bigpm(pR, mult.all.bigpm(lp));
	}
	return(reduce.var.pm(pR));
}
div.bigpm = function(p1, p2, by="x", debug=TRUE) {
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
			pRez = sum.bigpm(pRez, px1);
			tmp = mult.bigpm(px1, p2); tmp$coeff = - tmp$coeff;
			p1 = sum.bigpm(p1, tmp);
		}
	}
	if(debug) {
		if(nrow(p1) > 0) print("Not divisible!")
		else print("Divisible!");
	}
	return(list(Rez=pRez, Rem=p1));
}
# solve.bigpm()
# - is basically the same:
#   only uses mult.bigpm(), sum.bigpm() & replace.fr.bigpm();
#   [currently hard-coded in solve.pm()]

