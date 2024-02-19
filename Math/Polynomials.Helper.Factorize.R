########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Multi-Variable Polynomials
### Factorize
###
### draft v.0.1h-refactor


### Factorize Multi-Variable Polynomials


### fast load:
# - is automatically loaded in: Polynomials.Helper.R;
# source("Polynomials.Helper.R")
# source("Polynomials.Helper.Factorize.R")

### may need BigNumbers
# library(gmp)

# needed to Factorize Symmetric Polynomials:
source("Polynomials.Helper.Categories.R");

# needed to Factorize other Polynomials:
source("Polynomials.Helper.Mod.R");

######################

# TODO: a lot of work;
# - reduce variables;
# - multi-variable support;
# - advanced factorization;


### Factorize
# - only univariate Polynomial;
factorizeExt.p = function(p, by = xn, xn = "x", asBigNum = FALSE, skip.squares = FALSE,
		file=NULL, debug=FALSE) {
	# Level 1:
	if( ! skip.squares) {
		pR = factorize.p(p, by=by, asBigNum=asBigNum, f.all=FALSE, file=file, debug=debug);
		if(length(pR) > 0 && ! is.null(pR[[1]]$GCD)) return(pR);
	}
	### Level 2: P(1/x)
	pinv = rev.pm(p, xn=by);
	pR = factorize0.p(p, pinv, by=by, asBigNum=asBigNum, asSquares=FALSE, f.all=FALSE,
		file=file, debug=debug);
	if(length(pR) > 0 && ! is.null(pR[[1]]$GCD)) return(pR);
	### Level 3: P( - 1/x)
	isOdd = (pinv[, by] %% 2 == 1);
	pinv$coeff[isOdd] = - pinv$coeff[isOdd];
	pR = factorize0.p(p, pinv, by=by, asBigNum=asBigNum, asSquares=FALSE, f.all=FALSE,
		file=file, debug=debug);
	#
	return(pR);
}
### Factorizing Squares
factorize.p = function(p, by = xn, xn = "x", f.all=FALSE, asBigNum=TRUE,
		file="_R.Temp.", debug=FALSE) {
	id = match(by, names(p));
	if(is.na(id)) stop("Variable NOT present!");
	# gcd(p, D(p))
	dp = dp.pm(p, by);
	if(nrow(dp) == 0 || ncol(dp) < 2) return(list(list(GCD=NULL, p1=p)));
	pRez = factorize0.p(p, dp, by=by, f.all=f.all, asBigNum=asBigNum, file=file, debug=debug);
	return(pRez)
}
### Factorizes using: GCD(p, dp)
factorize0.p = function(p, dp, by = xn, xn = "x", f.all=FALSE, asBigNum=TRUE, asSquares=TRUE,
		file="_R.Temp.", debug=FALSE) {
	# factorize.all = FALSE
	# - p1 is usually sufficient;
	# - does NOT fully handle:
	#   p1 * p2^3*p3^4 or p1^2*p2^3 => p1*p2^2 * (...);
	### Checks:
	id = match(by, names(p));
	if(is.na(id)) stop("Variable NOT present!");
	if(nrow(dp) == 0 || ncol(dp) < 2) return(list(list(GCD=NULL, p1=p)));
	if(asBigNum && ! inherits(p$coeff, c("bigz", "bigq")))
		warning("Polynomial is NOT of type BigNumber!");
	### Factorization:
	lvl = 1; # level of factorization; (TODO)
	doSave = ! (is.null(file) || is.na(file));
	rez = list();
	while(TRUE) {
		# Step 1: GCD
		cat("\n"); print(paste0("Level = ", lvl));
		pGCD = gcd.exact.p(p, dp, by, asBigNum=asBigNum, debug=debug);
		pGCD = reduce.var.pm(pGCD);
		if(nrow(pGCD) < 1 || ncol(pGCD) < 2) break;
		id = match(by, names(pGCD));
		if(is.na(id)) break;
		# Leading Sign:
		isMaxPow = (pGCD[,id] == max(pGCD[,id]));
		if(pGCD$coeff[isMaxPow][[1]] < 0) pGCD$coeff = - pGCD$coeff;
		if(doSave) write.csv(pGCD, file=paste0(file, "GCD.", lvl, ".csv"), row.names=FALSE);
		# Step 2:
		cat("\n"); print("Starting division:");
		p.all = div.pm(p, pGCD, by)$Rez;
		if(asBigNum) {
			if(all(denominator(p.all$coeff) == 1)) p.all$coeff = as.bigz(p.all$coeff)
			else print("Warning: some Denominators != 1!")
		}
		if(doSave) write.csv(p.all, file=paste0(file, "ALL.", lvl, ".csv"), row.names=FALSE);
		# Data
		rez[[lvl]] = list();
		rez[[lvl]][["GCD"]] = pGCD;
		rez[[lvl]][["All"]] = p.all;
		# Step 3: additional Step useful for factorizing Squares
		if(asSquares) {
			cat("\n"); print("Starting Step 3:");
			p.minus1 = gcd.exact.p(pGCD, p.all, by, asBigNum=asBigNum);
			# TODO: IF(p.minus1 == p.all) => multiplicity!
			p1 = div.pm(p.all, p.minus1, by)$Rez;
			if(doSave) write.csv(p1, file=paste0(file, "p1.", lvl, ".csv"), row.names=FALSE);
			rez[[lvl]][["p1"]]  = p1;
		}
		#
		if( ! f.all) break;
		lvl = lvl + 1;
		p  = pGCD;
		dp = dp.pm(pGCD, by);
		if(nrow(dp) < 1 || ncol(dp) < 2) break;
		if(is.na(match(by, names(dp)))) break;
	}
	return(rez);
}

# TODO:
# - change sign in Result: if (x^max) < 0;

# proof of concept:
factorizeByB0.p = function(p, xn=by, by="x", pow=2, max.rows=100, digits=6, debug=TRUE) {
	b0 = B0.pm(p, xn=xn);
	if(nrow(b0) == 0) return(p); # TODO
	if(nrow(b0) > 1 || ncol(b0) > 2) stop("Only univariate polynomials supported!");
	scaleBack = function(pF, sc) {
		pF$coeff = pF$coeff / sc^pF[, xn];
		pF$coeff = pF$coeff * sc^pow;
		if(digits > 0) pF$coeff = round(pF$coeff, digits=digits);
		pF = reduce.cpm(pF);
		pF = toPoly.pm(pF);
		return(pF);
	}
	f.grid = expand.primes(b0$coeff, pow=pow, max.rows=max.rows);
	for(nr in seq(nrow(f.grid))) {
		### Scale x:
		ptmp = p;
		sc   = prod(f.grid[nr,]);
		if(debug) cat(paste0("\nScaling = ", sc, "^", pow));
		ptmp$coeff = ptmp$coeff * sc^ptmp[, xn];
		pR = factorizeExt.p(ptmp, xn=xn, skip.squares=TRUE, file=NULL);
		if(length(pR) > 0 && ! is.null(pR[[1]]$GCD)) {
			return(scaleBack(pR[[1]]$GCD, sc));
		}
	}
	print("No factors found!");
	invisible();
}

#######################

rev.pm = function(p, xn=by, by="x", sort=TRUE) {
	idx = match(xn, names(p));
	if(is.na(idx)) {
		warning("Variable name not found!");
		return(p);
	}
	p[ , xn] = max(p[ , xn]) - p[ , xn];
	if(sort) p = sort.pm(p, xn);
	return(p);
}

### Factorize x:
factors.numeric = function(x) {
	# requires: library(pracma)
	x = abs(x);
	f = data.frame(F = pracma::factors(x), count=1);
	f = aggregate(count ~ F, f, length);
	return(f)
}
# - expand factors with power >= pow;
expand.primes = function(x, pow=2, max.rows=100, one.rm=TRUE) {
	f = factors.numeric(x);
	if(pow > 1) {
		f$count = f$count %/% pow; # multiple of pow;
		f = f[f$count > 0, , drop=FALSE];
	}
	f = lapply(seq(nrow(f)), function(id) f$F[id]^seq(0, f$count[id]));
	f.grid = expand.grid(f);
	if(one.rm) f.grid = f.grid[-1, , drop=FALSE];
	if( ! is.null(max.rows) && nrow(f.grid) > max.rows) {
		f.grid = f.grid[seq(max.rows), , drop=FALSE];
		warning("There are more factors than max.rows!\nThe additional factors have been removed.")
	}
	return(f.grid)
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

# Analysis for Factoring Algorithms
factorSimple.pm = function(p, by = "x") {
	if(ncol(p) > 2) stop("Only univariate polynomials!");
	xn = by; # TODO
	p = sort.pm(p, xn);
	len = nrow(p);
	if(p[len, xn] != 0) stop("Divisible by x!"); # TODO
	#
	b0 = p$coeff[len];
	if(abs(b0) != 1) warning("Factorize: Not yet implemented!");
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

###################
###################

### Special Cases

### V1 P9: Almost-Symmetric


checkE21_S3 = function(e, e21) {
	S = e[1]; E2 = e[2]; E3 = e[3];
	E21 = e21[1]; E12 = e21[2];
	# Checks:
	chk1 = E21 + E12 - E2*S + 3*E3;
	if(round0(chk1) != 0) return(FALSE);
	chk2 = E21*E12 - (E2^3 + 3*E3^2 - 3*E3*E2*S) +
		- E3*(S^3 - 3*E2*S + 3*E3) - 3*E3^2;
	if(round0(chk2) != 0) return(FALSE);
	#
	return(TRUE);
}

factorize.V1P9.QuasiSym = function(p, digits=8) {
	if( ! is.pm(p)) stop("Not a polynomial!");
	len = ncol(p);
	if(len > 2) stop("p must be a Univariate polynomial!");
	idc = match("coeff", names(p));
	if(is.na(idc)) stop("Not a polynomial: Missing coefficients!");
	idx = seq(2)[ - idc];
	xn  = names(p)[idx];
	# P[9]
	pow = max(p[, idx]);
	if(pow != 9) {
		warning("Polynomial must be of Order 9!");
		return(list(p=p, Factor=NULL));
	}
	# Symmetric Polynomial:
	msg = "Not an Almost Symmetric Polynomial!";
	wf = function(msg) {
		warning(msg);
		return(list(p, Factor=NULL));
	}
	bx = coef.pm(p, xn=xn, pow=0:9);
	isSymm = isSymmetric.numeric(bx, len=4);
	if( ! isSymm) return(wf(msg));
	#
	if(bx[1] != 1) bx = bx / bx[1];
	S = bx[2]; E2 = bx[3] - S; E3 = bx[4] - S^2 + E2 - 3;
	b4 = bx[5]; b5 = bx[6];
	if(b4 == b5) {
		# Fully Symmetric:
		m = multiplicity.pm(p, -1);
		if(m > 0) {
			pDiv = toPoly.pm(paste0(xn, "+ 1"));
			# TODO: arg "n" => "pow";
			pRez = div.pm(p, pow.pm(pDiv, n=m, debug = FALSE), by=xn)$Rez;
			pRez = toPoly.pm(pRez);
			return( list(p=pRez, Factor=pDiv, Multiple=m) );
		}
	}
	# Checks:
	if( ! checkE21_S3(c(S,E2,E3), c(b4, b5) - (2*S + E2))) return(wf(msg));
	# Factors
	b = roots(c(1, -S, E2, -E3));
	# Round to Integer:
	b = sapply(b, function(b) {
		if(round(b) == round(b, digits=digits)) round(b) else b;
	})
	asP = function(b) {
		pR = data.frame(x=3:0, coeff=c(1, b[1], b[2], 1));
		names(pR)[1] = xn;
		return(toPoly.pm(pR));
	}
	# TODO: E21a => order of roots matters!
	b = b[c(1,3,2)];
	pR = list(p = asP(b[1:2]), Factor = asP(b[2:3]), Factor2 = asP(b[c(3,1)]));
	return(pR);
}

