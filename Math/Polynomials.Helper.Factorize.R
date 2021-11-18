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

# TODO:
# - change sign in Result: if (x^max) < 0;


#######################

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
#  => NOT: P[1] * P[n-1];

