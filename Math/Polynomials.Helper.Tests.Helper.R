########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: Helper Functions
### Tests


### fast load:
source("Polynomials.Helper.R")


### this file:
# source("Polynomials.Helper.Tests.Helper.R")


####################

### Helper Functions

# Check value of a specific coefficient
checkCoeff.pm = function(p, val, pow=1, xn="x") {
	if(is.null(xn)) {
		coeff = p$coeff;
	} else {
		coeff = p$coeff[p[, xn] == pow];
	}
	if(is.na(val)) {
		if(length(coeff) >= 1) stop("Wrong number of Monoms!");
	} else {
		if(length(coeff) == 0) {
			if(val != 0) stop("Missing Monom!");
		} else if(length(coeff) > 1) { stop("Wrong number of Monoms!"); }
		else if(coeff != val) stop("Wrong value!");
	}
	cat("Coeff: Success!\n");
	invisible(TRUE);
}
checkMaxPow.pm = function(p, val, xn="x") {
	idx = match(xn, names(p));
	if(is.na(idx)) {
		stopifnot(val == 0);
		cat("Power == 0: Success!\n");
		return(invisible());
	}
	pow = max(p[, xn]);
	if(pow != val) stop("Wrong power!");
	cat("Power: Success!\n");
	invisible(TRUE);
}
checkVal.pm = function(pval, val) {
	print(pval);
	if(inherits(pval, "data.frame")) {
		if(nrow(pval) > 1 || ncol(pval) > 1) stop("Value is a Data frame!");
		pval = pval$coeff;
	}
	stopifnot(pval == val);
	cat("Value: Success!\n");
	invisible(TRUE);
}
checkEmpty.pm = function(p) {
	stopifnot(nrow(p) == 0);
	cat("Empty: Success!\n");
	invisible(TRUE);
}
checkLength.pm = function(p, nr) {
	stopifnot(nrow(p) == nr);
	cat("Length: Success!\n");
	invisible(TRUE);
}
checkWarning.pm = function(wrn, len=1, txt=NULL) {
	stopifnot( ! is.null(wrn) && length(wrn) == len);
	if( ! is.null(txt)) {
		len = length(txt);
		stopifnot(names(wrn)[seq(len)] == txt);
	}
	cat("Warning: Success!\n");
	invisible(TRUE);
}

