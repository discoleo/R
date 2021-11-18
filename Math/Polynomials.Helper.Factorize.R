########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Multi-Variable Polynomials
### Factorize
###
### draft v.0.1b


### Factorize Multi-Variable Polynomials


### fast load:
# - is automatically loaded in: Polynomials.Helper.R;
# source("Polynomials.Helper.Factorize.R")

### may need BigNumbers
# library(gmp)

######################

# TODO: a lot of work;
# - reduce variables;
# - multi-variable support;
# - advanced factorization;


### Factorize
# - only univariate Polynomial;
factorizeExt.p = function(p, xn="x", asBigNum=FALSE, file=NULL, debug=FALSE) {
	pR = factorize.p(p, xn=xn, asBigNum=asBigNum, file=file, f.all=FALSE, debug=debug);
	if(length(pR) > 0 && ! is.null(pR[[1]]$GCD)) return(pR);
	### Level 2: P(1/x)
	pinv = rev.pm(p, xn=xn);
	pR = factorize0.p(p, pinv, xn=xn, asBigNum=asBigNum, file=file, f.all=FALSE, debug=debug);
	if(length(pR) > 0 && ! is.null(pR[[1]]$GCD)) return(pR);
	### Level 2: P( - 1/x)
	isOdd = (pinv[,xn] %% 2 == 1);
	pinv$coeff[isOdd] = - pinv$coeff[isOdd];
	pR = factorize0.p(p, pinv, xn=xn, asBigNum=asBigNum, file=file, f.all=FALSE, debug=debug);
	#
	return(pR);
}
factorize.p = function(p, xn="x", f.all=FALSE, asBigNum=TRUE,
		file="_R.Temp.", debug=FALSE) {
	id = match(xn, names(p));
	if(is.na(id)) stop("Variable NOT present!");
	# gcd(p, D(p))
	dp = dp.pm(p, xn);
	if(nrow(dp) == 0 || ncol(dp) < 2) return(list(list(GCD=NULL, p1=p)));
	pRez = factorize0.p(p, dp, xn=xn, f.all=f.all, asBigNum=asBigNum, file=file, debug=debug);
	return(pRez)
}
factorize0.p = function(p, dp, xn="x", f.all=FALSE, asBigNum=TRUE,
		file="_R.Temp.", debug=FALSE) {
	# factorize.all = FALSE
	# - p1 is usually sufficient;
	# - dos NOT fully handle:
	#   p1 * p2^3*p3^4 or p1^2*p2^3 => p1*p2^2 * (...);
	### Checks:
	id = match(xn, names(p));
	if(is.na(id)) stop("Variable NOT present!");
	if(nrow(dp) == 0 || ncol(dp) < 2) return(list(list(GCD=NULL, p1=p)));
	### Factorization:
	lvl = 1; # level of factorization; (TODO)
	doSave = ! (is.null(file) || is.na(file));
	rez = list();
	while(TRUE) {
		# Step 1: GCD
		cat("\n"); print(paste0("Level = ", lvl));
		pGCD = gcd.exact.p(p, dp, xn, asBigNum=asBigNum, debug=debug);
		pGCD = reduce.var.pm(pGCD);
		if(nrow(pGCD) < 1 || ncol(pGCD) < 2) break;
		id = match(xn, names(pGCD));
		if(is.na(id)) break;
		# Leading Sign:
		isMaxPow = (pGCD[,id] == max(pGCD[,id]));
		if(pGCD$coeff[isMaxPow][[1]] < 0) pGCD$coeff = - pGCD$coeff;
		if(doSave) write.csv(pGCD, file=paste0(file, "GCD.", lvl, ".csv"), row.names=FALSE);
		# Step 2:
		cat("\n"); print("Starting division:");
		p.all = div.pm(p, pGCD, xn)$Rez;
		if(asBigNum) {
			if(all(denominator(p.all$coeff) == 1)) p.all$coeff = as.bigz(p.all$coeff)
			else print("Warning: some Denominators != 1!")
		}
		if(doSave) write.csv(p.all, file=paste0(file, "ALL.", lvl, ".csv"), row.names=FALSE);
		cat("\n"); print("Starting Step 3:");
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

# TODO:
# - change sign in Result: if (x^max) < 0;


#######################

rev.pm = function(p, xn="x", sort=TRUE) {
	idx = match(xn, names(p));
	if(is.na(idx)) {
		warning("Variable name not found!");
		return(p);
	}
	p[ , xn] = max(p[ , xn]) - p[ , xn];
	if(sort) p = sort.pm(p, xn);
	return(p);
}

rPoly = function(n, coeff=c(0,1,-1), p=NULL, b0 = TRUE) {
	coeffs = sample(coeff, n, replace=TRUE, prob=p);
	if(b0 && coeffs[n] == 0) {
		coeffs[n] = sample(coeff[coeff != 0], 1);
	}
	coeffs = c(1, coeffs);
	isNZero = (coeffs != 0);
	pR = data.frame(x=seq(n, 0)[isNZero], coeff=coeffs[isNZero]);
	return(toPoly.pm(pR));
}

### Rules:

# if(b0 == 1)
#  => (... + 1) * (... + 1) OR
#  => (... - 1) * (... - 1);
#  => Test: P(1), P(-1);

factorSimple.pm = function(p) {
	if(ncol(p) > 2) stop("Only univariate polynomials!");
	xn = "x"; # TODO
	p = sort.pm(p, xn);
	len = nrow(p);
	if(p[len, xn] != 0) stop("Divisible by x!"); # TODO
	#
	b0 = p$coeff[len];
	if(abs(b0) != 1) warning("Not yet implemented!");
	# Trivial Factors:
	pval1 = c(eval.pm(p, 1), eval.pm(p, -1));
	isP1  = (pval1 == 0);
	if(any(isP1)) {
		pval1 = pval1[isP1];
		if(length(pval1) == 1) {
			pDiv = data.frame(x=1:0, coeff=c(1, pval1));
		} else {
			pDiv = data.frame(x=2:0, coeff=c(1, sum(pval1), prod(pval1)));
			pDiv = reduce.pm(pDiv);
		}
		names(pDiv)[1] = xn;
		pR = div.pm(p, pDiv, by=xn);
		return(pR$Rez);
	}
	### Primes:
	pp = c(3,5,7)
	for(pr in pp) {
		print(paste0("mod ", pr, ": ", pval1 %% pr));
	}
	pp = c(2,3,5)
	b1 = p$coeff[p[, xn] == 1];
	if(length(b1) == 0) b1 = 0;
	print("P^2:")
	for(pr in pp) {
		p2 = pr*pr;
		print(paste0("mod ", p2, ": ", (pval1*pr + b0) %% p2));
	}
	# TODO
}

