########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Multi-Variable Polynomials
### Categories
###
### draft v.0.1a


### Categories of Multi-Variable Polynomials


### fast load:
# source("Polynomials.Helper.Categories.R")


####################
####################

### Helper Functions

# - assumes that Polynomials.Helper.R is loaded;
# source("Polynomials.Helper.R")


####################
####################


isSymmetric.numeric = function(x, len=NULL) {
	if(len < 0) stop("Invalid length parameter!");
	len.half = length(x) %/% 2;
	isOdd    = (length(x) %% 2) == 1;
	posStart = if(isOdd) len.half + 2 else len.half + 1;
	if(is.null(len) || len >= length(x)) {
		posOff = 0;
	} else {
		posOff = if(len >= len.half) 0 else len.half - len;
	}
	isEq = (x[seq(len.half - posOff)] == x[seq(length(x), posStart + posOff)]);
	return(all(isEq));
}

isSymmetric.pm = function(p, pow.max=NULL) {
	# Strictly Symmetric: b0 == bn;
	# Assumes: reduced polynomial;
	nc = ncol(p);
	if(nc > 2) stop("p must be a Univariate polynomial!");
	idc = match("coeff", names(p));
	if(is.na(idc)) stop("Not a polynomial: Missing coefficients!");
	len = nrow(p);
	# 0-length Polynomial;
	if(len == 0) return(TRUE);
	if(nc == 1) {
		if(any(p[, idc] != 0)) return(FALSE);
		return(TRUE); # Coeff = 0;
	}
	# Variable:
	idx = seq(nc)[ - idc];
	xn  = names(p)[idx];
	#
	id = order(p[, xn]);
	p  = p[id, , drop=FALSE];
	maxPow = p[len, xn];
	### Checks
	if(len == 1 && p[1, idc] != 0) return(FALSE);
	if(p[1, xn] != 0 || p[1, idc] != p[len, idc]) return(FALSE);
	len.half = ((len + 1) %/% 2);
	if(is.null(pow.max)) {
		pow.len = len.half;
	} else {
		# TODO: How to encode powers?
		if(pow.max > maxPow) {
			warning(paste0("Checking only up to power", maxPow));
			pow.len = len.half;
		} else if(pow.max >= ((maxPow + 1) %/% 2)) {
			pow.len = len.half;
		} else {
			pow.max = min(pow.max, (maxPow + 1) %/% 2);
			pow.len = sum(p[, xn] <= pow.max);
		}
	}
	# "x^n + 1"
	if(pow.len <= 1) return(TRUE);
	if(pow.len >= len.half && (len %% 2 == 1)) {
		isEven = (maxPow %% 2 == 0);
		if(isEven) {
			if(p[pow.len, xn] != (maxPow %/% 2)) return(FALSE);
		} else {
			return(FALSE); # NOT possible;
		}
		if(pow.len <= 2) return(TRUE);
		pow.len = pow.len - 1;
	}
	for(nr in seq(2, pow.len)) {
		nr.end = len - nr + 1;
		if(p[nr, xn] != (maxPow - p[nr.end, xn])) return(FALSE);
		if(p[nr, idc] != p[nr.end, idc]) return(FALSE);
	}
	return(TRUE);
}

