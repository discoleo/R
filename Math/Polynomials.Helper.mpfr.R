########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: Helper Functions
### mpfr Functions
###
### draft v.0.2f


### fast load:
source("Polynomials.Helper.R")
source("Polynomials.Helper.Matrix.mpfr.R")


# Required libraries:
library(Rmpfr)


### this file:
# source("Polynomials.Helper.mpfr.R")


#######################
#######################

###############
### History ###
###############

### draft v.0.2:
# - det.mpfr, det.complex.mpfr;
# - generate mpfr Vandermonde determinants;
# - [0.2.f] moved Matrix code to new file:
#   Polynomials.Helper.Matrix.mpfr.R;
### draft v.0.1c:
# - various fixes;
### draft v.0.1a - v.0.1b:
# - moved mpfr-specific functions
#   to this file, from files:
#   Polynomials.Helper.R;
#   Polynomials.Class1.R;

# TODO: complex.mpfr vs cmpfr;


########################
########################

format.mpfr = function (x, digits = NULL, drop0trailing = TRUE, right = TRUE, 
		max.digits = getOption("Rmpfr.print.max.digits", 999L),
		exponent.plus = getOption("Rmpfr.print.exponent.plus", TRUE), ...) {
	stopifnot(is.mpfr(x), is.null(digits) || digits >= 1);
	n = length(x);
	if (n >= 1) {
		lFormat = function(x, na.print, print.gap, max, useSource, ...)
			format(x, digits = digits, max.digits = max.digits, 
				drop0trailing = drop0trailing, exponent.plus = exponent.plus, 
				...);
		x = lFormat(x, ...);
    }
    invisible(x);
}

prod.complex.mpfr = function(Re, Im) {
	x = Re; xi = Im;
	n = length(x);
	m = length(xi);
	if(n != m) {
		if(m == 1) {
			xi = rep(xi, n);
		} else stop("Wrong length!");
	}
	# Prod:
	rezr = x[1];
	rezi = xi[1];
	for(i in seq(2, n)) {
		tmp  = rezr * x[i] - rezi * xi[i];
		rezi = rezr * xi[i] + rezi * x[i];
		rezr = tmp;
	}
	return(list(Re = rezr, Im = rezi));
}

# Note: should be numerically more robust;
pow.all.complex.mpfr = function(Re, Im, n, start.zero = FALSE) {
	x = Re; xi = Im;
	prec = getPrec(x);
	z0 = mpfr(0, prec);
	z1 = mpfr(1, prec);
	if(n == 0) {
		if(start.zero) return(list(Re = z1, Im = z0));
		return(NULL);
	}
	y  = rep(x, n);
	yi = rep(xi, n);
	if(n == 1) {
		if(start.zero) { y = c(z1, y); yi = c(z0, yi); }
		return(list(Re = y, Im = yi));
	}
	# Powers:
	i1 = 1L;
	i2 = integer(0);
	# Powers of 2:
	while((tmp <- i1 + i1) <= n) {
		y[tmp]  = y[i1] * y[i1] - yi[i1] * yi[i1];
		yi[tmp] = 2 * y[i1] * yi[i1];
		i1 = tmp;
		i2 = c(i2, tmp);
	}
	for(i in i2) {
		id = seq(min(i-1, n - i));
		y[i + id]  = y[i] * y[id] - yi[i] * yi[id];
		yi[i + id] = y[i] * yi[id] + yi[i] * y[id];
	}
	if(start.zero) { y = c(z1, y); yi = c(z0, yi); }
	return(list(Re = y, Im = yi));
}
pow.all.complex.old.mpfr = function(Re, Im, n, start.zero = FALSE) {
	x = Re; xi = Im;
	prec = getPrec(x);
	z0 = mpfr(0, prec);
	z1 = mpfr(1, prec);
	if(n == 0) {
		if(start.zero) return(list(Re = z1, Im = z0));
		return(NULL);
	}
	y  = rep(x, n);
	yi = rep(xi, n);
	if(n == 1) {
		if(start.zero) { y = c(z1, y); yi = c(z0, yi); }
		return(list(Re = y, Im = yi));
	}
	for(i in seq(2, n)) {
		y[i]  = y[i-1] * Re - yi[i-1]*Im;
		yi[i] = y[i-1] * Im + yi[i-1]*Re;
	}
	if(start.zero) { y = c(z1, y); yi = c(z0, yi); }
	return(list(Re = y, Im = yi));
}

# (a + 1i*b)^n
pascal.complex.mpfr = function(n, precBits, as.pm  = TRUE) {
	b  = pascal.mpfr(n, precBits);
	if(as.pm) {
		lst = split.cpm.pascal(b);
	} else {
		lst = split.complex.pascal(b);
	}
	return(lst);
}
split.complex.pascal = function(x) {
	re = lapply(x, function(x) {
		re  = x[seq(1, length(x), by = 2)];
		len = length(re);
		if(len < 2) return(re);
		id = seq(2, len, by = 2);
		re[id] = - re[id];
		return(re);
	});
	im = lapply(x, function(x) {
		im  = x[seq(2, length(x), by = 2)];
		len = length(im);
		if(len < 2) return(im);
		id = seq(2, len, by = 2);
		im[id] = - im[id];
		return(im);
	});
	lst = list(Re = re, Im = im);
	return(lst);
}
split.cpm.pascal = function(x) {
	re = lapply(x, function(x) {
		re  = x[seq(1, length(x), by = 2)];
		len = length(re);
		if(len < 2) {
			# TODO: len == 0;
			tmp = data.frame(a = 1, b = 0);
			tmp$coeff = re;
			return(tmp);
		}
		id = seq(2, len, by = 2);
		re[id] = - re[id];
		len = length(x) - 1;
		tmp = data.frame(
			a = seq(len, 0, by = -2),
			b = seq(0, len, by=2));
		tmp$coeff = re;
		return(tmp);
	});
	im = lapply(x, function(x) {
		len = length(x);
		im  = x[seq(2, len, by = 2)];
		if(len <= 3) {
			# TODO: len == 0;
			if(len == 2) {
				tmp = data.frame(a = 0, b = 1);
			} else if(len == 3) {
				tmp = data.frame(a = 1, b = 1);
			}
			tmp$coeff = im;
			return(tmp);
		}
		len = length(im);
		id = seq(2, len, by = 2);
		im[id] = - im[id];
		len = length(x) - 1;
		tmp = data.frame(
			a = seq(len - 1, 0, by = -2),
			b = seq(1, len, by=2));
		tmp$coeff = im;
		return(tmp);
	});
	lst = list(Re = re, Im = im);
	return(lst);
}


###################
### Polynomials ###

### Roots

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


### Generators

# Compute polynomial from Roots:
poly.calc.mpfr = function(x, bits=120, tol=1E-7) {
	one  = mpfr(1, precBits=bits);
	zero = mpfr(0, precBits=bits);
	p  = mpfr2array(c(Re=one, Im=zero), c(1,2));
	p0 = mpfr2array(c(Re=zero, Im=zero), c(1,2));
	if(inherits(x, "mpfrMatrix")) {
		len = nrow(x);
		xRe = mpfr2array(x[,1], nrow(x));
		xIm = mpfr2array(x[,2], nrow(x));
	} else {
		xRe = Re(x); xIm = Im(x); len = length(x);
	}
	for (i in seq(len)) {
		px = mult.mpfr(p[,1], p[,2], xRe[i], xIm[i]);
		p = rbind(p0, p);
		p  = p - rbind(px, p0);
	}
	p = round0(p, tol=tol);
	return(p);
}

split.cmpfr.pm = function(p, precBits = NULL) {
	nc = ncol(p);
	if(nc > 2) warning("Only univariate polynomials!");
	if(is.numeric(p$coeff)) {
		p$coeff = mpfr(p$coeff, precBits=precBits);
		prec = precBits;
	} else prec = getPrec(p$coeff);
	#
	z0 = mpfr(0, prec);
	z1 = mpfr(1, prec);
	idX = which(names(p) != "coeff")[1];
	pMax = max(p[idX]);
	if(pMax == 0) return(list(Re = p, Im = data.frame(coeff = z0)));
	if(pMax == 1) {
		pC = p[p[idX] == 1,, drop = FALSE];
		names(p)[idX] = "a";
		names(pC)[idX] = "b";
		lst = list(Re = p, Im = pC);
		return(lst);
	}
	# Binomial Expansion:
	b = pascal.complex.mpfr(pMax, precBits = prec);
	# TODO
}

### Generate specific Roots

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


########################

### Evaluate Polynomials

eval.pm.mpfr = function(p, x, bits=120, tol=10^(10 - bits %/% 3)) {
	# TODO
}
eval.cpm = function(p, x, bits=120, tol=1E-10, doPolar=TRUE, progress=FALSE) {
	# uses the Rmpfr package;
	# currently assumes that only coeffs are big and
	# have an impact on numeric stability;
	pP = p[, - which(names(p) == "coeff"), drop=FALSE];
	### Powers:
	# currently only max:
	pow  = lapply(seq(ncol(pP)), function(nc) max(pP[,nc]));
	# pow = lapply(seq(ncol(pP)), function(nc) sort(unique(pP[,nc])));
	#
	xpows = pow.cpm.mpfr(x, pow=pow, bits=bits, tol=tol, doPolar=doPolar);
	#
	if(progress) cat("Processing row:\n");
	eval.p = function(id) {
		# Progress bar:
		if(progress && id %% 16 == 1) cat(paste0(id, if(id %% 96 == 1) "\n" else ", "));
		#
		idx = which(pP[id,] != 0);
		xpows = xpows[idx]; lenv = length(idx);
		### Free Term:
		if(lenv == 0) {
			re = mpfr(Re(p$coeff[id]), precBits=bits);
			im = mpfr(Im(p$coeff[id]), precBits=bits);
			return(mpfr2array(c(Re=re, Im=im, Div=1), c(3)));
		}
		#
		re = mpfr2array(sapply(seq(lenv), function(id2) xpows[[id2]][pP[id, idx[id2]], 1]), lenv);
		im = mpfr2array(sapply(seq(lenv), function(id2) xpows[[id2]][pP[id, idx[id2]], 2]), lenv);
		if(length(re) > 0) {
		while(TRUE) {
			len = length(re);
			if(len == 1) break;
			iend = (len %% 2);
			i  = seq(1, len - iend, by=2);
			re1 = re[i] * re[i+1] - im[i] * im[i+1];
			im1 = re[i] * im[i+1] + im[i] * re[i+1];
			if(iend > 0) {
				re.last = re[len]; im.last = im[len];
				re2 = re1[1] * re.last - im1[1] * im.last;
				im2 = re1[1] * im.last + im1[1] * re.last;
				re1[1] = re2; im1[1] = im2;
			}
			re = re1; im = im1;
		}
		} else {re = 1; im = 0;}
		re = re * p$coeff[id]; im = im * p$coeff[id];
		sol = mpfr2array(c(Re=re, Im=im, Div=1), c(3));
		return(sol);
	}
	if(progress) cat("\n");
	sol = sapply(seq(nrow(p)), eval.p);
	sdim = attr(sol, "dim"); sol = mpfr2array(t(sol), rev(sdim));
	sol = apply(sol, 2, sum);
	return(sol);
}

pow.cpm.mpfr = function(x, pow, bits=120, tol=1E-10, doPolar=TRUE) {
	# x = c(v1, v2, ..., vn); # multi-variable polynomials: n variables!
	# pow = list(pow.v1, pow.v2, ..., pow.vn);
	# - mpfr-quirk: vector is actually list!
	isList = ( ! inherits(x, "mpfr")) && inherits(x, "list");
	lenx = if(isList) seq(length(x)) else 1;
	# pow = is always list;
	if( ! inherits(pow, "list")) stop("pow must be a list!")
	#
	xpows = lapply(lenx, function(id) {
		x = if(isList) x[[id]] else x;
		isMpfr = inherits(x, "mpfr");
		# round small values to 0:
		x0 = if( ! isMpfr && tol > 0) round0(x, tol=tol) else x;
		len = tail(pow[[id]], 1);
		#
		isComplex = (is.complex(x0) && Im(x0) != 0) ||
			(isMpfr && length(x) >= 2 && x[2] != 0);
		if(isComplex) {
			div = 1; # defunct
			# polar coordinates:
			# - but less accuracy with certain complex numbers;
			# - needed when r^max.pow overflows;
			if(doPolar) {
				if(isMpfr) {
					re = x0[1]; im = x0[2];
				} else {
					re = mpfr(Re(x0), bits); im = mpfr(Im(x0), bits);
				}
				r = sqrt(re^2 + im^2);
				pib = Const("pi", bits); pih = pib / 2;
				# seems NO difference between asin & atan versions;
				th  = asin(abs(im/r));
				if(re == 0) {th = pih; if(im < 0) th = - th;}
				else if(re < 0) {
					th = pib + if(im > 0) - th else th;
				} else if(im < 0) {
					th = -th;
				}
				r  = r^seq(len);
				th = th * seq(len);
				re = r * cos(th); im = r * sin(th);
			} else {
				if(isMpfr) {
					# TODO
					stop("Not yet implemented!")
				}
				x = x0^seq(len);
				re = Re(x); im = Im(x);
				re = mpfr(re, bits);
				im = mpfr(im, bits);
			}
			return(cbind(Re=re, Im=im, Div=div));
		} else {
			x0 = if(isMpfr) x0[1] else mpfr(Re(x0), bits);
			x = x0^seq(len); # power 0 NOT needed;
			return(cbind(Re=x, Im=0, Div=1));
		}
	})
	return(xpows);
}


####################

### Helper Functions

### Roots of Unity
unity.mpfr = function(n=5, all=TRUE, bits=120, include1=FALSE) {
	pib = Const("pi", bits); pin = 2*pib / n;
	if( ! all) {
		m = c(cos(pin), sin(pin));
		m = mpfr2array(m, c(2));
	} else {
		from = if(include1) 0 else 1;
		m = sapply(seq(from, n-1), function(id) {
				c(cos(id*pin), sin(id*pin)) });
		m = t(mpfr2array(m, c(2, n - from)));
	}
	return(m);
}

### Conversions

as.double.cmpfr = function(x) {
	len = nrow(x);
	x = sapply(seq(len), function(nr) {
		complex(re=as.double(x[nr,1]), im=as.double(x[nr,2]));
	})
	return(x);
}

toPolar.lmpfr = function(x, bits=120) {
	pib = Const("pi", bits);
	xpol = sapply(seq(nrow(x)), function(nr) toPolar.mpfr(as.list(x[nr,]), bits=bits, piConst=pib));
	xpol = mpfr2array(xpol, c(2, nrow(x)));
	return(t(xpol));
}
toPolar.mpfr = function(x, bits=120, piConst=NULL) {
	if(inherits(x, "mpfrMatrix")) {
		re = x[,1]; im = x[,2];
	} else if(inherits(x, "mpfrArray")) {
		# already mpfr:
		re = x[1]; im = x[2];
	} else if(inherits(x, "mpfr")) {
		# already mpfr:
		re = x[[1]]; im = x[[2]];
	} else {
		re = mpfr(Re(x), bits); im = mpfr(Im(x), bits);
	}
	### Polar coordinates:
	r = sqrt(re^2 + im^2); # TODO: use function;
	pib = if( ! is.null(piConst)) piConst
		else Const("pi", bits);
	pih = pib / 2;
	# seems NO difference between asin & atan versions;
	th  = mpfr2array(asin(abs(im/r)), length(re));
	# Correct quadrant:
	isZero = (re == 0);
	th[isZero] = ifelse(im[isZero] < 0, pih, - pih);
	isQ23 = (re < 0);
	th[isQ23] = ifelse(im[isQ23] > 0, pib - th[isQ23], pib + th[isQ23]);
	isQ34 = ( ! isQ23) & (im < 0);
	th[isQ34] = -th[isQ34];
	rez = mpfr2array(c(r, th), c(length(th), 2));
	colnames(rez) = c("M", "Theta");
	return(rez);
}

### Math: Complex Numbers

### Multiplication
mult.mpfr = function(re1, im1, re2, im2) {
	reN = re1 * re2 - im1 * im2;
	imN = re1 * im2 + im1 * re2;
	return(cbind(Re=reN, Im=imN));
}

### Power
pow.mpfr = function(x, len=1, bits=120) {
	x = toPolar.mpfr(x, bits=bits);
	r = x[,1]; th = x[,2]; # TODO: Bug in mpfr selectors;
	doPow = FALSE;
	if(length(len) > 1) { pow = len; doPow = TRUE; }
	else if(len > 1) { pow = seq(len); doPow = TRUE; }
	if(doPow) {
		r  = r^pow;
		th = th * pow;
	}
	re = r * cos(th); im = r * sin(th);
	return(cbind(Re=re, Im=im));
}

### Radicals
rootn.mpfr = function(x, n) {
	if(n %% 2 == 0) {
		# TODO:
		stop("Not yet implemented!")
		if(is.list(x) || any(x < 0)) {
		} else {
		}
	} else {
		ninv = 1/mpfr(n, precBits=120);
		if(is.matrix(x) && dim(x)[2] >= 2) {
			# Complex numbers:
			bits = getPrec(x[1,1]);
			xpol = toPolar.lmpfr(x, bits=bits);
			r = sapply(seq(nrow(x)), function(id) {
				x = x[id, ];
				if(x[[2]] == 0) {
					r = if(x[[1]] >=0) x[[1]]^(1/n)
						else - (-x[[1]])^(1/n);
					return(c(r, 0));
				}
				r = xpol[id, 1]; th = xpol[id, 2];
				re = r^ninv; th = th*ninv;
				im = re * sin(th); re = re * cos(th);
				return(c(re, im));
			})
			r = t(mpfr2array(r, c(2, nrow(x))));
		} else {
			r = ifelse(x >= 0, x^ninv, - (-x)^ninv);
			r = mpfr2array(r, c(length(r)));
		}
	}
	return(r);
}

log.cmpfr = function(x, bits=120) {
	x = toPolar.mpfr(x, bits=bits);
	r1 = log(x[,1]);
	r2 = x[,2];
	r = cbind(r1, r2);
	len = if(inherits(x, "mpfrMatrix")) c(nrow(x), 2) else c(2);
	r = mpfr2array(r, len);
	return(r);
}

#################

