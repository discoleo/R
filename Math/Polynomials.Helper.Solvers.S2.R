########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Solvers: S2 Systems



#######################

# library(polynom)
# library(pracma)

### helper Functions
source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


### fast load:
# source("Polynomials.Helper.Solvers.S2.R")


###############
###############

###############
### Solvers ###
###############

### Hetero-Symmetric

### Decompose into Symmetric Polynomials
decompose.S2Ht = function(p, vars=c("x", "y"), drop=TRUE) {
	if(length(vars) != 2) stop("Only S2 can be solved!");
	id = match(vars, names(p));
	if(all(is.na(id))) stop("All variables are missing!");
	### Decompose system: x*y
	xy = "E2"; # paste0(vars, collapse="");
	p = replace.pm.character.pm(p, xy, toPoly.pm(paste0(vars, collapse="*")));
	id = match(vars, names(p));
	id = id[ ! is.na(id)];
	if(length(id) == 0 || (max.pow <- maxPow.pm(p, id)) == 0) {
		# pure polynomial in "x*y";
		return(list(P0 = p));
	}
	### x^n + y^n & x^n - y^n
	# Powers used for Decomposition:
	epow = E2.pm(max.pow);
	pEDiff = epow$pEDiff;
	pESum  = epow$pESum;
	### Decomposition:
	diff.f = function(p, idv) {
		pDiff = data.frame();
		for(nr in seq(nrow(p))) {
			pTmp = mult.pm(pEDiff[[p[nr, idv]]], p[nr, -id, drop=F]);
			pDiff = sum.pm(pDiff, pTmp); # TODO: do not aggregate;
		}
		if(nrow(pDiff > 0)) pDiff = aggregate0.pm(pDiff);
		return(pDiff);
	}
	sum.f = function(p, idv) {
		pSum = data.frame();
		for(nr in seq(nrow(p))) {
			pTmp = mult.pm(pESum[[p[nr, idv]]], p[nr, -id, drop=F]);
			pSum = sum.pm(pSum, pTmp); # TODO: do not aggregate;
		}
		if(nrow(pSum) > 0) pSum = aggregate0.pm(pSum);
		return(pSum);
	}
	# remaining variables:
	isX = (p[, id[1]] > 0);
	pDiff = diff.f(p[isX, , drop=F], id[1]);
	pSum  = sum.f(p[isX, , drop=F], id[1]);
	isY = (p[, id[2]] > 0);
	pTmp = diff.f(p[isY, , drop=F], id[2]);
	pDiff = diff.pm(pDiff, pTmp);
	pSum  = sum.pm(pSum, sum.f(p[isY, , drop=F], id[2]));
	isNot = ! (isX | isY);
	if(any(isNot)) {
		tmp = p[isNot, , drop=F];
		tmp$coeff = tmp$coeff * 2;
		pSum = sum.pm(pSum, tmp);
	}
	if(drop) { pDiff = drop.pm(pDiff); pSum = drop.pm(pSum); }
	return(list(pDiff=pDiff, pSum=pSum));
}

