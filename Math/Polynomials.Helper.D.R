########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Derivation


### fast load:
# source("Polynomials.Helper.D.R")

### requires:
# - but is already loaded inside;
# source("Polynomials.Helper.R")

######################

###############
### History ###
###############


### draft v.0.1a:
# - moved Section on Derivation
#   from Polynomials.Helper.R;


#####################
#####################


# D( p(x) )
dp.pm = function(p, xn="x") {
	p = p[(p[,xn] != 0), , drop=FALSE];
	if(nrow(p) == 0) return(0);
	p$coeff = p$coeff * p[,xn];
	p[,xn] = p[,xn] - 1;
	return(p);
}
# D( p1 * exp(p2) )
dp.exp.pm = function(p, xn="x") {
	pr = mult.pm(dp.pm(p$Exp, xn), p$Poly);
	pr = sum.pm(pr, dp.pm(p$Poly, xn));
	p$Poly = pr;
	return(p);
}
# D( p1 / pdiv ])
dp.div.pm = function(p1, pdiv, xn="x") {
	if(is.numeric(p1) || is.complex(p1)) {
		r = mult.pm(dp.pm(pdiv, xn=xn), - p1);
		div = mult.pm(pdiv, pdiv);
		return(list(p = r, div = div));
	}
	r  = mult.pm(pdiv, dp.pm(p1, xn=xn));
	r2 = mult.pm(p1, dp.pm(pdiv, xn=xn));
	r = diff.pm(r, r2);
	r = list(p=r, div = mult.pm(pdiv, pdiv));
	return(r)
}

### D( ODE )

# D( dny ) => d(n+1)y;
dy.pm = function(p, yn="y", xn="x", starts.with=FALSE) {
	regEnd = if(starts.with) "" else "$";
	reg = paste0("^(?:d[0-9]*+|)", yn, regEnd);
	ynms = grepl(reg, names(p), perl=TRUE);
	ynms = names(p)[ynms];
	if(length(ynms) == 0) {
		# Note: NOT 0 ???
		if(is.null(xn)) return(p);
		return(dp.pm(p, xn=xn));
	}
	# Dx
	r = if(is.null(xn)) 0 else dp.pm(p, xn);
	# Dy
	dy.f = function(yn, dyn) {
		dp = p[p[,yn] != 0, , drop=FALSE];
		if(nrow(dp) == 0) return(0);
		dp$coeff = dp$coeff * dp[,yn];
		dp[,yn]  = dp[,yn] - 1;
		if( ! dyn %in% names(dp)) dp[,dyn] = 0;
		dp[,dyn] = dp[,dyn] + 1;
		return(dp);
	}
	y0 = if(starts.with) {
		ynms[startsWith(ynms, yn)];
	} else {
		ynms[ynms %in% yn];
	}
	if(length(y0) > 0) {
		for(y00 in y0) {
			dy = paste0("d", y00);
			r = sum.pm(r, dy.f(y00, dy));
		}
		ynms = ynms[ ! ynms %in% y0];
		if(length(ynms) == 0) return(r);
	}
	#
	dy = Dy.names(ynms);
	for(id in seq_along(ynms)) {
		r = sum.pm(r, dy.f(ynms[[id]], dy[[id]]));
	}
	#
	return(r);
}
Dy.names = function(nms) {
	# processes only "dy" & "dny";
	# dy.pm: "(?=.)" should NOT be necessary!
	npos = regexpr("^d[0-9]*+(?=.)", nms, perl=TRUE);
	len = attr(npos, "match.length");
	if(any(len < 0)) stop("Invalid dy name!");
	#
	n = rep(1, length(len));
	isDn = (len > 1);
	s = substr(nms[isDn], 2, len[isDn]); # pmax(2, len[isDn])
	n[isDn] = as.integer(s);
	n = n + 1;
	nms = substr(nms, len + 1, nchar(nms));
	nms = paste0("d", n, nms);
	return(nms)
}

