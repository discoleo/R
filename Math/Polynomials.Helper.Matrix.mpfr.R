########################
###
### Leonard Mada
### [the one and only]
###
### Matrix Functions: mpfr
###
### draft v.0.2g

# Matrix Operations using mpfr;

# Note:
# - code is independent of code for Polynomials;


### Load ALL files:
# source("Polynomials.Helper.mpfr.R")
# this file:
# source("Polynomials.Helper.Matrix.mpfr.R")

####################

### Helper Functions

# Required libraries:
# library(Rmpfr)

as.matrix.mpfrMatrix = function(x) {
	n = dim(x);
	x = matrix(as.numeric(x), nrow = n[1], ncol = n[2]);
	return(x);
}

as.pow.mpfr = function(x, precBits, pow = 1, pow.div = 2) {
	z = mpfr(x, precBits=precBits);
	p = mpfr(pow, precBits=precBits) / pow.div;
	isOdd = pow.div %% 2 == 1;
	isNeg = FALSE;
	if(isOdd) {
		idNeg = which(x < 0);
		if(length(idNeg) > 0) isNeg = TRUE;
		if(isNeg) z[idNeg] = - z[idNeg];
		if(pow %% 2 == 0) isNeg = FALSE;
	}
	z = z^p;
	if(isNeg) z[idNeg] = - z[idNeg];
	return(z);
}

# Note:
# - also defined in file:
#   Polynomials.Helper.mpfr.R;
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


####################
####################

##############
### Matrix ###

# Solve:
# - rather naive, but probably sufficiently robust
#   by increasing slightly the number of bits;
solve.mpfr = function(b, y, transpose = TRUE) {
	if(transpose) b = t(b);
	mdim = dim(b);
	nn = mdim[1];
	if(nn != mdim[2]) stop("Please provide a Square Matrix!");
	if(nn == 0) return(NULL);
	#
	prec = getPrec(b[1,1]);
	if(is.numeric(y)) {
		y = mpfr(y, precBits = prec);
	}
	if(nn == 1) {
		sol = y / b;
		return(sol);
	}
	# Upper:
	z0 = mpfr(0, precBits = prec);
	for(nr in seq(nn - 1)) {
		b1 = b[nr, nr];
		if(b1 == z0) {
			# Swap columns:
			isZero = TRUE;
			for(i in seq(nr + 1, nn)) {
				if(x[nr, i] != 0) {
					isZero = FALSE; break;
				}
			}
			if(isZero) return(NULL);
			# sg  = - sg;
			tmp = x[, nr]; x[, nr] = x[, i]; x[, i] = tmp;
			tmp = y[nr]; y[nr] = y[i]; y[i] = tmp;
			b1  = x[nr, nr];
		}
		for(nc in seq(nr + 1, nn)) {
			b2 = b[nr, nc];
			if(b2 != z0) {
				ff = - b2 / b1;
				y[nc] = y[nc] + ff * y[nr];
				b[, nc] = b[, nc] + ff * b[, nr];
			}
		}
	}
	# det = prod(diag(b));
	# Lower: Backwards!
	for(nr in seq(nn, 2)) {
		b1 = b[nr, nr];
		for(nc in seq(nr - 1, 1)) {
			b2 = b[nr, nc];
			ff = - b2 / b1;
			y[nc] = y[nc] + ff*y[nr];
		}
	}
	sol = y / diag(b);
	return(sol);
}
solve.old.mpfr = function(b, y, transpose = TRUE) {
	if(transpose) b = t(b);
	mdim = dim(b);
	if(mdim[1] != mdim[2]) stop("Please provide a Square Matrix!");
	if(mdim[1] == 0) return(NULL);
	prec = getPrec(b[1,1]);
	b0 = b;
	z0 = mpfr(0, precBits = prec);
	LU = list();
	# Determinant:
	# - but NOT really needed;
	# TODO: move to new function;
	for(nr in seq(mdim[1] - 1)) {
		Tr = rep(z0, dim = mdim[2] - nr);
		b1 = b[nr, nr];
		for(nc in seq(nr + 1, mdim[2])) {
			b2 = b[nr, nc];
			if(b2 != z0) {
				ff = - b2 / b1;
				Tr[nc - nr] = ff;
				b[, nc] = b[, nc] + ff * b[, nr];
			}
		}
		LU[[nr]] = Tr;
	}
	det = prod(diag(b));
	# Solution:
	for(nr in seq(2, mdim[1])) {
		ff = LU[[nr - 1]];
		y[seq(nr, mdim[1])]  = y[seq(nr, mdim[1])] + ff * y[nr - 1];
	}
	# Lower: Backwards!
	for(nr in seq(mdim[1], 2)) {
		b1 = b[nr, nr];
		for(nc in seq(nr - 1, 1)) {
			b2 = b[nr, nc];
			ff = - b2 / b1;
			y[nc] = y[nc] + ff*y[nr];
		}
	}
	Sol = y / diag(b);
	#
	sol = list(Sol = Sol, Tr = LU, det = det);
	return(sol);
}

####################
### Determinants ###

# O(N^3)
determinant.complex = function(x, ...) {
	nn = dim(x);
	if(nn[1] != nn[2]) stop("Please provide a Square Matrix!");
	# Special Cases:
	nn = nn[1];
	if(nn == 0) return(NULL);
	if(nn == 1) return(x[1,1]);
	if(nn == 2) return(x[1,1]*x[2,2] - x[2,1]*x[1,2]);
	sg = 1;
	# Determinant:
	for(nr in seq(nn - 1)) {
		b1 = x[nr, nr];
		if(b1 == 0) {
			# Swap columns:
			isZero = TRUE;
			for(i in seq(nr + 1, nn)) {
				if(x[nr, i] != 0) {
					isZero = FALSE; break;
				}
			}
			if(isZero) return(0);
			sg  = - sg;
			tmp = x[, nr]; x[, nr] = x[, i]; x[, i] = tmp;
			b1  = x[nr, nr];
		}
		for(nc in seq(nr + 1, nn)) {
			b2 = x[nr, nc];
			if(b2 != 0) {
				ff = - b2 / b1;
				# another O(N) here;
				x[, nc] = x[, nc] + ff * x[, nr];
			}
		}
	}
	rez = prod(diag(x));
	if(sg < 0) rez = - rez;
	return(rez);
}

# O(N^3 * bits)
# normalize = may sometimes improve the precision,
#  but only slightly and inconsistently;
det.mpfr = function(x, normalize = FALSE) {
	return(determinant.mpfr(x, normalize=normalize));
}
determinant.mpfr = function(x, normalize = FALSE, ...) {
	nn = dim(x);
	if(nn[1] != nn[2]) stop("Please provide a Square Matrix!");
	# Special Cases:
	nn = nn[1];
	if(nn == 0) return(NULL);
	if(nn == 1) return(x[1,1]);
	if(nn == 2) return(x[1,1]*x[2,2] - x[2,1]*x[1,2]);
	#
	prec = getPrec(x[1,1]);
	z0 = mpfr(0, precBits = prec);
	z1 = mpfr(1, precBits = prec);
	sg = 1;
	# Normalization
	if(normalize) {
		# mean, median: do NOT work very well;
		# TODO: optimal strategy?
		minz = apply(x, 2, function(x) min(abs(x)));
		idBig = which(minz > 2);
		norm = z1;
		if(length(idBig) > 0) {
			norm = prod(minz[idBig]);
			for(id in idBig) {
				x[, id] = x[, id] / minz[id];
			}
		}
	}
	# Determinant:
	for(nr in seq(nn - 1)) {
		b1 = x[nr, nr];
		if(b1 == z0) {
			# Swap columns:
			isZero = TRUE;
			for(i in seq(nr + 1, nn)) {
				if(x[nr, i] != 0) {
					isZero = FALSE; break;
				}
			}
			if(isZero) return(z0);
			sg  = - sg;
			tmp = x[, nr]; x[, nr] = x[, i]; x[, i] = tmp;
			b1  = x[nr, nr];
		}
		for(nc in seq(nr + 1, nn)) {
			b2 = x[nr, nc];
			if(b2 != z0) {
				ff = - b2 / b1;
				# another O(N) here;
				x[, nc] = x[, nc] + ff * x[, nr];
			}
		}
	}
	rez = prod(diag(x));
	if(sg < 0) rez = - rez;
	if(normalize) rez = rez * norm;
	return(rez);
}

# Backwards Sequential processing of Columns:
# - should behave better with Vandermonde determinants;
determinant.seq.mpfr = function(x, ...) {
	nn = dim(x);
	if(nn[1] != nn[2]) stop("Please provide a Square Matrix!");
	# Special Cases:
	nn = nn[1];
	if(nn == 0) return(NULL);
	if(nn == 1) return(x[1,1]);
	if(nn == 2) return(x[1,1]*x[2,2] - x[2,1]*x[1,2]);
	#
	prec = getPrec(x[1,1]);
	z0 = mpfr(0, precBits = prec);
	sg = 1; # not needed;
	# Determinant:
	for(nr in seq(nn - 1)) {
		for(nc in seq(nn, nr + 1)) {
			nc.prev = nc - 1;
			b1 = x[nr, nc.prev];
			if(b1 == z0) {
				# Pseudo-Swap columns:
				isZero = TRUE;
				for(i in seq(nc - 1, nr)) {
					if(x[nr, i] != 0) {
						isZero = FALSE; break;
					}
				}
				# TODO: check;
				if(isZero) {
					if(x[nr, nc] == 0) return(z0);
					tmp = x[, nr]; x[, nr] = x[, nc]; x[, nc] = tmp;
					next;
				}
				b1 = x[nr, i];
				nc.prev = i;
			}
			b2 = x[nr, nc];
			if(b2 != z0) {
				ff = - b2 / b1;
				# another O(N) here;
				x[, nc] = x[, nc] + ff * x[, nc.prev];
			}
		}
	}
	rez = prod(diag(x));
	if(sg < 0) rez = - rez;
	return(rez);
}

# Complex: det(x + 1i*xi)
# TODO: Test thoroughly!
det.complex.mpfr = function(x, xi) {
	nn = dim(x);
	ni = dim(xi);
	if(nn[1] != nn[2] || ni[1] != ni[2]) stop("Please provide a Square Matrix!");
	if(nn[1] != ni[1]) stop("The 2 matrixes must have the same dimensions!");
	if(nn[1] == 0) return(NULL);
	prec = getPrec(x[1,1]);
	z0 = mpfr(0, precBits = prec);
	nn = nn[1];
	sg = 1;
	# Determinant:
	for(nr in seq(nn - 1)) {
		b1 = x[nr, nr];
		bi = xi[nr, nr];
		if(b1 == z0 && bi == z0) {
			# Swap columns:
			isZero = TRUE;
			for(i in seq(nr + 1, nn)) {
				if(x[nr, i] != z0 || xi[nr, i] != z0) {
					isZero = FALSE; break;
				}
			}
			if(isZero) return(z0);
			sg  = - sg;
			tmp = x[, nr]; x[, nr] = x[, i]; x[, i] = tmp;
			tmp = xi[, nr]; xi[, nr] = xi[, i]; xi[, i] = tmp;
			b1  = x[nr, nr];
			bi  = xi[nr, nr];
		}
		div  = b1^2 + bi^2;
		invr = b1 / div;
		invi = - bi / div;
		for(nc in seq(nr + 1, nn)) {
			b2 = x[nr, nc];
			b2i = xi[nr, nc];
			ff = - b2 * invr + b2i * invi;
			fi = - b2 * invi - b2i * invr;
			# another O(N) here;
			tmpr = x[, nc] + ff * x[, nr] - fi * xi[, nr];
			tmpi = xi[, nc] + ff * xi[, nr] + fi * x[, nr];
			x[, nc]  = tmpr;
			xi[, nc] = tmpi;
		}
	}
	# rez = prod(diag(x));
	rr = x[1,1]; ri = xi[1,1];
	for(i in seq(2, nn)) {
		b = x[i,i]; bi = xi[i,i];
		tmpr = rr*b - ri*bi;
		tmpi = rr*bi + ri*b;
		rr = tmpr; ri = tmpi;
	}
	if(sg < 0) { rr = - rr; ri = - ri; }
	return(list(Re = rr, Im = ri));
}

### Matrix Generators
# e.g. to test Determinants

vandermond.mpfr = function(x) {
	warning("The correct name is: vandermonde!");
	vandermonde.mpfr(x);
}
vandermonde.mpfr = function(x) {
	n = length(x);
	prec = getPrec(x[1]);
	m = mpfrArray(1, prec, dim = c(n,n));
	tmp = x;
	for(id in seq(2, n)) {
		m[,id] = tmp;
		tmp = tmp * x;
	}
	return(m);
}
vandermonde.complex.mpfr = function(Re, Im) {
	x = Re; xi = Im;
	n = length(x);
	prec = getPrec(x[1]);
	mr = mpfrArray(1, prec, dim = c(n,n));
	mi = mpfrArray(0, prec, dim = c(n,n));
	tmpr = x; tmpi = xi;
	for(id in seq(2, n)) {
		mr[,id] = tmpr;
		mi[,id] = tmpi;
		tmp  = tmpr * x - tmpi * xi;
		tmpi = tmpr * xi + tmpi * x;
		tmpr = tmp;
	}
	return(list(Re = mr, Im = mi));
}
# x = Square Matrix or
#     Vector corresponding to the pow = 1 column;
det.vandermonde.mpfr = function(x, byRow = FALSE) {
	# Checks:
	dimx = dim(x);
	if( ! is.null(dimx)) {
		if(dimx[1] != dimx[2]) stop("Not a square matrix!");
		if(dimx[1] == 1) {
			if(x[1,1] != 1) warning("Not a Vandermonde matrix!");
			return(x[1,1]);
		}
		x = if(byRow) x[2,] else x[,2];
	}
	n = length(x);
	if(n == 1) {
		stop("Not a Vandermonde matrix!");
	}
	#
	rez = mpfr(1, precBits = getPrec(x[1]));
	for(i in seq(1, n-1)) {
		tmp = x[i] - x[seq(i+1, n)];
		rez = rez * prod(tmp);
	}
	# Square free:
	n4 = n %% 4;
	if(n4 != 0 && n4 != 1) rez = - rez;
	return(rez);
}
# Complex Vandermonde:
det.vandermonde.complex.mpfr = function(Re, Im) {
	x = Re; y = Im;
	if(is.null(y)) return(det.vandermonde.mpfr(x));
	# Checks:
	dimx = dim(x); dimy = dim(y);
	if( ! is.null(dimx)) {
		if(dimx[1] != dimx[2]) stop("Not a square matrix!");
		if(dimx[1] == 1) {
			if(x[1,1] != 1) warning("Not a Vandermonde matrix!");
			if(y[1,1] != 0) warning("Not a Vandermonde matrix!");
			return(list(Re = x[1,1], Im = y[1,1]));
		}
		x = x[,2];
	}
	if( ! is.null(dimy)) {
		y = y[,2];
	}
	n = length(x);
	m = length(y);
	if(n != m) stop("Non-equal dimensions!");
	#
	rezr = mpfr(1, precBits = getPrec(x[1]));
	rezi = mpfr(0, precBits = getPrec(x[1]));
	for(i in seq(1, n-1)) {
		tmpr = x[i] - x[seq(i+1, n)];
		tmpi = y[i] - y[seq(i+1, n)];
		rezr = c(rezr, tmpr);
		rezi = c(rezi, tmpi);
		tmp  = prod.complex.mpfr(rezr, rezi);
		rezr = tmp$Re; rezi = tmp$Im;
	}
	# Square free:
	n4 = n %% 4;
	if(n4 != 0 && n4 != 1) { rezr = - rezr; rezi = - rezi; }
	return(list(Re = rezr, Im = rezi));
}
