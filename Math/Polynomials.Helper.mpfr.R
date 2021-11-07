########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: Helper Functions
### mpfr Functions
###
### draft v.0.1a


### fast load:
source("Polynomials.Helper.R")


# required libraries
library(Rmpfr)


### this file:
# source("Polynomials.Helper.mpfr.R")


#######################
#######################

###############
### History ###
###############


### draft v.0.1a:
# - moved mpfr-specific functions
#   to this file, from file:
#   Polynomials.Helper.R;


########################
########################


### Compute Polynomials

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

### Evaluate Polynomials


eval.cpm = function(p, x, bits=120, tol=1E-10, doPolar=TRUE, progress=FALSE) {
	# uses the Rmpfr package;
	# currently assumes that only coeffs are big and
	# have an impact on numeric stability;
	pP = p[, - which(names(p) == "coeff"), drop=FALSE];
	# pow = lapply(seq(ncol(pP)), function(nc) sort(unique(pP[,nc])));
	# currently only max:
	pow = lapply(seq(ncol(pP)), function(nc) max(pP[,nc]));
	xpows = lapply(seq(length(x)), function(id) {
		x0 = round0(x[id], tol=tol);
		len = tail(pow[[id]], 1);
		if(is.complex(x0) && Im(x0) != 0) {
			div = 1;
			# polar coordinates:
			# - but less accuracy with certain complex numbers;
			# - needed when r^max.pow overflows;
			if(doPolar) {
				re = mpfr(Re(x0), bits); im = mpfr(Im(x0), bits);
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
				x = x0^seq(len);
				re = Re(x); im = Im(x);
				re = mpfr(re * div, bits);
				im = mpfr(im * div, bits);
			}
			return(cbind(Re=re, Im=im, Div=div));
		} else {
			x0 = mpfr(Re(x0), bits);
			x = x0^seq(len); # power 0 NOT needed;
			if(x0 == 0) {
				return(cbind(Re=x, Im=0, Div=1));
			} else {
				div = 1; # 12 - round(log(abs(x)) / log(10));
				return(cbind(Re=x, Im=0, Div=div));
			}
		}
	})
	if(progress) cat("Processing row:\n");
	eval.p = function(id) {
		if(progress && id %% 16 == 1) cat(paste0(id, if(id %% 96 == 1) "\n" else ", "));
		idx = which(pP[id,] != 0);
		xpows = xpows[idx]; lenv = length(idx);
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


####################

### Helper Functions

toPolar.lmpfr = function(x, bits=120) {
	pib = Const("pi", bits);
	xpol = sapply(seq(nrow(x)), function(nr) toPolar.mpfr(as.list(x[nr,]), bits=bits, piConst=pib));
	xpol = mpfr2array(xpol, c(2, nrow(x)));
	return(t(xpol));
}
toPolar.mpfr = function(x, bits=120, piConst=NULL) {
	if(is.list(x)) {
		# already mpfr:
		re = x[[1]]; im = x[[2]];
	} else {
		re = mpfr(Re(x), bits); im = mpfr(Im(x), bits);
	}
	r = sqrt(re^2 + im^2); # TODO: use function;
	pib = if( ! is.null(piConst)) piConst
		else Const("pi", bits);
	pih = pib / 2;
	# seems NO difference between asin & atan versions;
	th  = asin(abs(im/r));
	if(re == 0) {th = pih; if(im < 0) th = - th;}
	else if(re < 0) {
		th = pib + if(im > 0) - th else th;
	} else if(im < 0) {
		th = -th;
	}
	return(cbind(M=r, Theta=th));
}

### Multiplication
mult.mpfr = function(re1, im1, re2, im2) {
	reN = re1 * re2 - im1 * im2;
	imN = re1 * im2 + im1 * re2;
	return(cbind(Re=reN, Im=imN));
}

### Power
pow.pmfr = function(x, len=1, bits=120) {
	x = toPolar(x, bits=bits);
	r = x$M; th = x$Theta;
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

