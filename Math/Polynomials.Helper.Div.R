########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Multi-Variable Polynomials
### Polynomial Division


### fast load:
# - is automatically loaded in Polynomials.Helper.R;
# source("Polynomials.Helper.Div.R")


################
################

################
### Division ###
################

div.pm = function(p1, p2, by="x", NF.stop=TRUE, debug=TRUE) {
	if(is.character(p2) || inherits(p2, "formula")) p2 = toPoly.pm(p2, env=parent.frame());
	if(is.character(p1) || inherits(p1, "formula")) p1 = toPoly.pm(p1, env=parent.frame());
	# move Coeffs to last column:
	idc1 = match("coeff", names(p1));
	if(idc1 < ncol(p1)) {
		p1 = cbind(p1[ , -idc1, drop=FALSE], p1[ , idc1, drop=FALSE]);
		idc1 = ncol(p1);
	}
	# very simple division
	xn = by[1];
	idx  = checkDiv.pm(p1, p2, xn=xn);
	idx1 = idx[[1]]; idx2 = idx[[2]];
	# should NOT be necessary anymore; [DEPRECATED]
	# if( ! is.data.frame(p2)) p2 = as.data.frame(p2);
	#
	xpow2 = max(p2[,idx2]);
	pDx = p2[p2[,idx2] == xpow2, , drop=FALSE];
	pDx = drop.pm(pDx); # only vars from Leading monomial;
	idcDx = match("coeff", names(pDx));
	# idc1  = match("coeff", names(p1));
	c2 = pDx[, idcDx];
	pRez = as.data.frame(array(0, c(0,2)));
	names(pRez) = c(xn, "coeff");
	#
	idn = match(names(pDx)[-idcDx], names(p1));
	allMatched = ! any(is.na(idn));
	# print(paste0("Leading Variables: all matched = ", allMatched));
	if( ! allMatched) {
		if(NF.stop) {
			msg = "No matching variables for Leading Divisor:\n";
			nms = names(pDx)[-idcDx][is.na(idn)];
			if(length(nms) > 1) nms = paste0(nms, collapse=", ");
			nms[1] = paste0("  ", nms[1]); # nicer formatting
			stop(msg, nms);
		}
	}
	if(nrow(pDx) == 1) {
		if(allMatched) {
		while(TRUE) {
			if(nrow(p1) == 0) break;
			# Leading Monomials:
			xpow1 = max(p1[,xn]);
			if(xpow1 < xpow2) break;
			px1 = p1[p1[,xn] == xpow1, , drop=FALSE]; # TODO: check;
			# Diff Powers:
			for(nc in seq_along(idn)) {
				px1[, idn[nc]] = px1[, idn[nc]] - pDx[, nc];
			}
			px1[, idc1] = px1[, idc1] / c2;
			pRez = add.pm(pRez, px1);
			p1 = diff.pm(p1, mult.pm(px1, p2));
		}
		} else {
			# NOT all Leading variables matched
			# stop("Not yet implemented: NF-option!");
			# only Monomial:
			mDx = pDx[1, - idcDx, drop=FALSE];
			pDivL = NULL;
			while(TRUE) {
				if(nrow(p1) == 0) break;
				# Leading Monomials:
				xpow1 = max(p1[, xn]);
				if(xpow1 < xpow2) break;
				# check leading variables:
				pLead = p1[p1[, xn] == xpow1, , drop=FALSE];
				pLead = drop.pm(pLead);
				cat("Lead:\n"); print(pLead); # DEBUG
				# TODO: may be > 1 Leading rows;
				idn = match(names(mDx), names(pLead));
				isNot = is.na(idn);
				idn = idn[ ! isNot];
				isNot = isNot | sapply(idn, function(nc) any(pLead[ , nc] == 0));
				if(any(isNot)) {
					# Add missing variables
					m0 = mDx[, isNot, drop=FALSE];
					p1 = mult.m(p1, m0);
					if(nrow(pRez) > 0) pRez = mult.m(pRez, m0);
					idn = match(names(mDx), names(p1));
					pDivL = if(is.null(pDivL)) {
							m0$coeff = 1; m0;
						} else mult.m(pDivL, m0);
				}
				px1 = p1[p1[,xn] == xpow1, , drop=FALSE]; # TODO: check;
				# Diff Powers:
				for(nc in seq_along(idn)) {
					px1[, idn[nc]] = px1[, idn[nc]] - pDx[, nc];
				}
				px1[, idc1] = px1[, idc1] / c2;
				pRez = add.pm(pRez, px1);
				p1 = diff.pm(p1, mult.pm(px1, p2));
			}
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
	pRez = list(Rez=pRez, Rem=p1);
	if( ! allMatched) pRez$pDiv = pDivL;
	class(pRez) = c("pm.div", "list");
	return(pRez);
}
# only Remainder & eval(pDiv) == 0;
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

### Helper Functions:
checkDiv.pm = function(p1, pDiv, xn="x") {
	idx2 = match(xn, names(pDiv));
	if(is.na(idx2)) stop(paste0("P2 must contain the variable: ", xn));
	idx1 = match(xn, names(p1));
	if(is.na(idx1)) stop(paste0("P1 must contain the variable: ", xn));
	return(c(idx1, idx2));
}

###########

###########
### GCD ###
###########

gcd.vpm = function(p, xgcd=0) {
	### [Hack] mpfr
	# - skip the gcd;
	# - Note: coefficients may grow larger than they should be!
	if(inherits(p$coeff, "mpfr")) return(1);
	#
	if(xgcd == 0 && inherits(p$coeff, c("bigz", "bigq"))) {
		xgcd = as.bigz(0);
	}
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
gcd.exact.p = function(p1, p2, xn="x", asBigNum=TRUE, doGCD=TRUE, debug=FALSE) {
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
		if(dn > 0) { p2m[,xn] = p2m[,xn] + dn; }
		else if(dn < 0) p1m[,xn] = p1m[,xn] - dn;
		#
		dp = diff.pm(p1m, p2m);
		if(debug) print(toPoly.pm(dp)); # e.g. overflows massively;
		if(nrow(dp) == 0) {
			print("Factor found!");
			if(doGCD) return(p2) else return(list(p=0, f=0));
		}
		# simplify the coefficients: robust for BigNumbers;
		xgcd = gcd.vpm(dp, xgcd=dp$coeff[1]);
		if( ! is.na(xgcd)) {
			if(xgcd != 1) {
				dp$coeff = dp$coeff / xgcd;
				if(asBigNum) dp$coeff = as.bigz(dp$coeff);
			}
			if( ! doGCD) fact = fact * c2 / xgcd;
		} else if( ! doGCD) fact = fact * c2;
		if(debug) print(toPoly.pm(dp)); # e.g. overflows massively;
		# Remaining x:
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

### Multivariate Polynomials
# - initial attempt;
# - but far more challenging;
gcd.pm.exact = function(p1, p2, xn="x", asBigNum=NULL, doGCD=TRUE, multi.stop=FALSE,
			debug=FALSE, MAX.ITER=10) {
	if(is.null(asBigNum)) asBigNum = inherits(p1$coeff, c("bigz", "bigq"));
	if( ! doGCD) fact = if(asBigNum) as.bigz(1) else 1;
	lead.f = function(p, n, id=1) p[ p[, xn[[id]]] == n, , drop=FALSE];
	hasManyRows.f = function(p, n, id=1) sum(p[, xn[[id]]] == n) > 1;
	#
	iter = 0;
	while(TRUE) {
		if(MAX.ITER > 0) {
			iter = iter + 1;
			if(iter > MAX.ITER) {
				return(list(p1=p1, p2=p2));
			}
		}
		n1 = max(p1[, xn[[1]]]); n2 = max(p2[, xn[[1]]]);
		if( ! doGCD && n1 < n2) return(list(p=p1, f=fact));
		if(doGCD) {
			hasOneRow1 = ! hasManyRows.f(p1, n1);
			if((n1 < n2) || (n1 == n2 && hasOneRow1)) {
				tmp = p1; p1 = p2; p2 = tmp;
				tmp = n1; n1 = n2; n2 = tmp;
			}
		}
		# Leading Monomial:
		mL2 = lead.f(p2, n2);
		idMulti = 1;
		doElimination = FALSE;
		if(nrow(mL2) > 1) {
			# stop("Multi-Lead: Not yet supported!");
			if(length(xn) == 1 && multi.stop[[1]]) {
				warning("Multi-Lead: Not yet supported!");
				return(list(p1=p1, p2=p2));
			}
			if(length(xn) > 1) {
				# Different approach:
				# - sequential processing: benefit still not known;
				if(length(multi.stop) > 1) multi.stop = multi.stop[-1];
				return(gcd.pm.exact(p1, p2, xn=xn[-1], asBigNum=asBigNum, doGCD=doGCD,
					multi.stop=multi.stop, debug=debug, MAX.ITER=MAX.ITER));
				# [old approach]
				idMulti = idMulti + 1;
				n2  = max(mL2[, xn[[idMulti]]]);
				mL2 = lead.f(mL2, n2, id=idMulti);
				if(nrow(mL2) > 1) {
					warning("Even bigger Multi-Lead: Not yet supported!");
					return(list(p1=p1, p2=p2));
				}
			} else {
				doElimination = TRUE;
			}
		}
		if(doElimination) {
			mL1 = lead.f(p1, n1);
			div = gcd.vpm(mL1, xgcd = mL1$coeff[1]);
			if(div > 1) div = gcd.vpm(mL2, xgcd = div);
			if(div > 1) {
				mL1$coeff = mL1$coeff / div;
				mL2$coeff = mL2$coeff / div;
				if(asBigNum) { mL1$coeff = as.bigz(mL1$coeff); mL2$coeff = as.bigz(mL2$coeff); }
			}
			dp = diff.pm(mult.pm(p1, mL2), mult.pm(p2, mL1));
		} else {
			# piece-wise elimination
		mL1 = lead.f(p1, n1);
		if(idMulti > 1) {
			n1 = max(mL1[, xn[[idMulti]]]);
			if(n1 == 0) warning("Variable NOT present in Lead!");
			mL1 = lead.f(mL1, n1, id=idMulti);
		}
		# TODO: optimize selection of monomial;
		mL1 = mL1[1, , drop=FALSE];
		# Coefficients:
		c1 = mL1$coeff;
		c2 = mL2$coeff;
		div = gcd(c1, c2);
		if(div != 1) {
			c1 = c1 / div; c2 = c2 / div;
			if(asBigNum) {c1 = as.bigz(c1); c2 = as.bigz(c2);}
		}
		p1m = p1; p1m$coeff = p1m$coeff * c2;
		p2m = p2; p2m$coeff = p2m$coeff * c1;
		# Powers:
		mLL = align.pm(mL1, mL2, align.names=TRUE, doReduce=FALSE);
		mL1 = mLL[[1]]; mL2 = mLL[[2]];
		nms = names(mL1);
		idc = match("coeff", nms);
		for(nc in seq(ncol(mL1))) {
			if(nc == idc) next;
			dn = mL1[[nc]] - mL2[[nc]];
			nm = nms[nc];
			if(is.na(match(nm, names(p1m)))) p1m[, nm] = 0;
			if(is.na(match(nm, names(p2m)))) p2m[, nm] = 0;
			#
			if(dn > 0) { p2m[,nm] = p2m[,nm] + dn; }
			else if(dn < 0) p1m[,nm] = p1m[,nm] - dn;
		}
		#
		dp = diff.pm(p1m, p2m);
		}
		### Simplifications:
		if(doGCD && nrow(dp) > 0) {
			nms = names(dp);
			idc = match("coeff", nms);
			nms = nms[ - idc];
			# Reduce the Powers:
			for(nc in nms) {
				minPow = min(dp[ , nc]);
				# print(nc); print(dp);
				if(minPow > 0) dp[ , nc] = dp[ , nc] - minPow;
			}
			dp = drop.pm(dp);
		}
		# if(debug) print(toPoly.pm(dp)); # e.g. overflows massively;
		if(nrow(dp) == 0) {
			print("Factor found!");
			if(doGCD) return(p2) else return(list(p=0, f=0));
		}
		# simplify the coefficients: robust for BigNumbers;
		xgcd = gcd.vpm(dp, xgcd=dp$coeff[1]);
		if( ! is.na(xgcd)) {
			if(xgcd != 1) {
				dp$coeff = dp$coeff / xgcd;
				if(asBigNum) dp$coeff = as.bigz(dp$coeff);
			}
			if( ! doGCD) fact = fact * c2 / xgcd;
		} else if( ! doGCD) fact = fact * c2;
		if(debug) print(toPoly.pm(dp)); # e.g. overflows massively;
		# Remaining x:
		idVar = match(xn[[1]], names(dp));
		if(is.na(idVar)) {
			# can still have factor: e.g. by sequential processing and
			# a variable NOT present in the factor;
			warning("Something went wrong!");
			return(list(p1 = dp, p2 = p2));
		}
		n0 = max(dp[, xn[[1]], drop=TRUE]);
		print(paste0("Pow = ", n0, ", Len = ", nrow(dp)));
		if(n0 == 0) {
			print("Not divisible!");
			if(doGCD) return(dp);
			return(list(p=dp, f=fact));
		}
		p1 = dp;
		print(dp); # more DEBUG;
	}
}

############################

### Fraction Rationalization

rational.pm = function(p, p0, by=NULL, debug=TRUE) {
	if(is.null(by)) {
		if(ncol(p0) > 2) stop("Variable NOT specified!");
		idc = match("coeff", names(p0));
		by  = names(p0)[ - idc];
	}
	# TODO: check name in p;
	powMax = max(p0[, by]);
	if(max(p[, by]) >= powMax) {
		p = div.pm(p, p0, by=by, debug=debug)$Rem;
		# NO roots anymore in Denominator:
		if(max(p[, by]) == 0) return(list(Rez=as.pm(1), Div=p));
	}
	pR = data.frame(coeff=numeric(0));
	pD = p;
	fs = 1;
	while(TRUE) {
		pQ = div.pm(p0, pD, by=by, debug=debug);
		if(nrow(pR) == 0) { pR = pQ$Rez; }
		else {
			pR = mult.pm(pR, pQ$Rez);
			# Reduce polynomial:
			if(max(pR[, by]) > powMax)
				pR = div.pm(pR, p0, by=by, debug=debug)$Rem;
		}
		if(nrow(pQ$Rem) == 0) {
			warning("Division by 0!");
			return(list(Rez=NULL, Div=pD));
		}
		#
		fs = - fs;
		pD = pQ$Rem;
		if(debug) print(as.character(pD$coeff));
		isRez = is.na(match(by, names(pD))) || (max(pD[, by]) == 0);
		if(isRez) {
			if(fs < 0) pR$coeff = - pR$coeff;
			return(list(Rez = pR, Div = pD));
		}
	}
}