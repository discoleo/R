

### Eval d(P(x)) at the roots
eval.dpAtRoots.mpfr = function(x) {
	# base-function in file: Polynomials.Helper.D.R
	n = length(x);
	if(n <= 1) return(0);
	tmp = lapply(seq(n), function(id) {
		tmp = prod(x[id] - x[-id]);
	});
	tmp = do.call(Rmpfr:::c.mpfr, tmp);
	return(tmp)
}

### Gauss-Legendre Quadrature
intw = function (lower, upper, n = 81, prec = 128) {
	# Gauss-Legendre: adapted from package pracma;
	if(missing(upper)) {
		if(length(lower) == 2) {
			upper = lower[2]; lower = lower[1];
		} else stop("Missing upper Limit!");
	}
	if(inherits(lower, "mpfr")) {
		if(is.null(prec)) prec = getPrec(lower);
	} else {
		if(is.null(prec)) stop("Precision must be specified!");
		lower = mpfr(lower, prec);
		upper = mpfr(upper, prec);
	}
	#
	i = seq(1, n - 1, by = 1);
	d = i / sqrt(4*i^2 - 1);
	E = eigen(pracma::Diag(d, 1) + pracma::Diag(d, -1),
		symmetric = TRUE);
	L = E$values; ids = order(L);
	L = mpfr(L, prec);
	d = (upper - lower);
	x = L[ids];
	w = 1 / (1 - x^2) / eval.dpAtRoots.mpfr(x)^2;
	w = w / sum(w); w = d * w;
	x = 0.5 * (d*x + upper + lower);
	# Low-precision weights:
	V = E$vectors; # has been recomputed;
	V = t(V[, ids]); ww = d * V[, 1]^2;
	return(list(x = x, w = w, ww = ww))
}

