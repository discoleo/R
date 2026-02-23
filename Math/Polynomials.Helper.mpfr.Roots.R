########################
##
## Leonard Mada
## [the one and only]
##
## Polynomials: Helper Functions
## mpfr Functions: Polynomial Roots
##
## draft v.0.2g


### Roots

roots.Ext.mpfr = function(p, which = 1, by = "x", digits = 240, iter = 3,
		cpl.names = c("a", "b", "i"), verbose = 1) {
	rr = roots.pm(p);
	if(which < 0) stop("Invalid id of Root!");
	if(iter  < 0) stop ("Invalid number of Iterations!");
	if(iter == 0) return(rr); # Not yet implemented
	r  = rr[which];
	# Verbose:
	if(is.logical(verbose)) {
		verbose = if(verbose) 1 else 0;
	}
	isVerbose = verbose >= 1;
	# Complex root:
	pC = data.frame(a = c(1,0), b = c(0,1), i = c(0,1), coeff = c(1,1));
	names(pC)[1:3] = cpl.names;
	nm.i = cpl.names[3];
	pz = replace.pm(p, pC, by);
	pz = replace.pm(pz, -1, nm.i, pow = 2)
	p1 = drop.pm(pz[pz[nm.i] == 0, , drop = FALSE]);
	p2 = pz[pz[nm.i] != 0, , drop = FALSE];
	p2[nm.i] = 0; p2 = drop.pm(p2);
	#
	p1a = dp.pm(p1, by=cpl.names[1]);
	p1b = dp.pm(p1, by=cpl.names[2]);
	p2a = dp.pm(p2, by=cpl.names[1]);
	p2b = dp.pm(p2, by=cpl.names[2]);
	p1a$coeff = mpfr(p1a$coeff, digits);
	p1b$coeff = mpfr(p1b$coeff, digits);
	p2a$coeff = mpfr(p2a$coeff, digits);
	p2b$coeff = mpfr(p2b$coeff, digits);
	#
	p1$coeff = mpfr(p1$coeff, digits);
	p2$coeff = mpfr(p2$coeff, digits);
	# Compute improved Roots: Iterate
	a = mpfr(Re(r), digits);
	b = mpfr(Im(r), digits);
	#
	for(step_i in seq(iter)) {
		vals = list(a=a, b=b);
		names(vals) = cpl.names[1:2];
		v1a = eval.pm(p1a, vals);
		v1b = eval.pm(p1b, vals);
		v2a = eval.pm(p2a, vals);
		v2b = eval.pm(p2b, vals);
		v1  = eval.pm(p1,  vals);
		v2  = eval.pm(p2,  vals);
		#
		div = (v1a*v2b - v2a*v1b);
		dx  = (v2b*v1 - v1b*v2) / - div;
		dy  = (v2a*v1 - v1a*v2) / div;
		# Verbose:
		if(isVerbose) {
			err1 = v1a*dx + v1b*dy;
			err2 = v2a*dx + v2b*dy;
			if(verbose == 1) {
				err = sqrt(err1^2 + err2^2);
				cat("Residual = "); print(err);
			} else {
				cat("Residual Re = "); print(err1);
				cat("Residual Im = "); print(err2);
			}
		}
		# Update:
		a = a + dx; b = b + dy;
	}
	r.mpfr = list(Re = a, Im = b);
	if(isVerbose) cat("\n");
	return(r.mpfr);
}

# Complex Root:
# - only 1 Root!
roots.mpfr = function(p, x0, x0i, precBits = 200, iter = 24,
		multiplicity = 1, warn.multiplicity = TRUE, verbose = TRUE) {
	if(multiplicity < 1) stop("Invalid multiplicity!")
	if(ncol(p) > 2) stop("Polynomial must be univariate!");
	if(is.complex(p$coeff)) stop("Complex polynomials not yet implemented!");
	# Precision:
	isMpfr = inherits(x0, "mpfr");
	if(isMpfr) {
		precBits = getPrec(x0)[1];
	}
	if(is.pm(p)) {
		if(inherits(p$coeff, "mpfr")) {
			precBits = getPrec(p$coeff)[1];
		} else {
			p$coeff = mpfr(p$coeff, precBits=precBits);
		}
	}
	if( ! isMpfr) {
		x0  = mpfr(x0,  precBits=precBits);
		x0i = mpfr(x0i, precBits=precBits);
	}
	ddh = mpfr(paste0("1E-", round(precBits/6)), precBits=precBits);
	dd0 = ddh*ddh;
	idX = which(names(p) != "coeff");
	xn = names(p)[idX];
	dp = dp.pm(p, by = xn);
	nn = max(p[, xn]);
	id = p[, xn] + 1; idD = dp[, xn] + 1;
	if(multiplicity > 1) p$coeff = multiplicity * p$coeff;
	# Iterations:
	for(i in seq(iter)) {
		xpow = pow.all.complex.mpfr(x0, x0i, n = nn, start.zero = TRUE);
		pval = sum(xpow$Re[id] * p$coeff);
		pvi  = sum(xpow$Im[id] * p$coeff);
		dval = sum(xpow$Re[idD] * dp$coeff);
		dvi  = sum(xpow$Im[idD] * dp$coeff);
		### Div:
		# TODO: dp == 0
		div  = dval*dval + dvi*dvi;
		if(div <= dd0) {
			ddh = 2*ddh;
			if(abs(pval) <= ddh && abs(pvi) <= ddh) {
				if(warn.multiplicity) warning("Multiplicity!");
				if(verbose) cat("Iteration: ", i, "\n");
				if(multiplicity > 1) multiplicity = multiplicity - 1;
				sol = roots.mpfr(dp, x0, x0i, iter = max(24, iter - 8),
					multiplicity=multiplicity);
				return(sol);
			}
			warning("Division by 0!");
			sol = list(Re = x0, Im = x0i);
			return(sol);
		}
		dval = dval / div;
		dvi  = - dvi / div;
		# Newton-Raphson:
		x0  = x0 - pval*dval + pvi*dvi;
		x0i = x0i - pval*dvi - pvi*dval;
	}
	# TODO: test residual / convergence
	sol = list(Re = x0, Im = x0i);
	return(sol);
}

######################
### Specific Roots ###


### Class 1 Polynomials
roots.Class1.mpfr = function(K, s, n=5, bits=120) {
	if( ! inherits(K, "mpfr")) K = mpfr(K, precBits=bits)
	k = rootn.mpfr(K, n=n);
	m = unity.mpfr(n=n, all=TRUE, bits=bits);
	len = length(s);
	if(len > n) stop("Sequence too long!");
	s0 = 0; # TODO
	pows = seq(len);
	# TODO: complex k;
	k = mpfr2array(k[1]^pows, c(len));
	r = sapply(seq(n), function(id) {
		# TODO
		mp = if(id == 1) m
			else if(id == n) matrix(c(1,0), nrow=len, nc=2, byrow=TRUE)
			else m[ ((id*pows) %% n), ]; # TODO: non-primes + 1;
		re = sum(s*k*mp[,1]) + s0;
		im = sum(s*k*mp[,2]);
		return(c(re, im));
	})
	r = t(mpfr2array(r, c(2, n)));
}

