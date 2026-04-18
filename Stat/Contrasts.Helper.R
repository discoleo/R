

### Contrasts:
### Helper Functions

# Inner function in contr.poly;
make.poly <- function(n, scores = seq(n)) {
    y = scores - mean(scores);
    X = outer(y, seq_len(n) - 1, '^');
    QR = qr(X);
    z  = QR$qr;
    z  = z * (row(z) == col(z)); # z = diag(z)
    raw = qr.qy(QR, z);
    Z <- sweep(raw, 2L, apply(raw, 2L, function(x) sqrt(sum(x^2))), 
            '/', check.margin = FALSE)
    colnames(Z) <- paste0("^", 1L:n - 1L)
    Z
}

### Exact Formulas

# Gen: Component 3
poly.cp3 = function(n) {
	nm  = (n+1)/2; nc = seq(1, nm, by = 1);
	# div = sqrt(nm * (nm-1) * (8*nm^3 - 12*nm^2 - 2*nm + 3) / 5);
	div = sqrt((8*nm^5 - 20*nm^4 + 10*nm^3 + 5*nm^2 - 3*nm) / 5);
	(3*nc^2 - 3*(n+1)*nc + 2*nm^2 + nm) / div;
}
poly.cp3.mpfr = function(n, prec = 240) {
	if(inherits(n, "mpfr")) {
		v1 = Rmpfr::mpfr(1, getPrec(n));
	} else {
		n  = Rmpfr::mpfr(n, prec);
		v1 = Rmpfr::mpfr(1, prec);
	}
	nm  = (n+1)/2; nc = seq(v1, nm, by = 1);
	# div = sqrt(nm * (nm-1) * (8*nm^3 - 12*nm^2 - 2*nm + 3) / 5);
	div = sqrt((8*nm^5 - 20*nm^4 + 10*nm^3 + 5*nm^2 - 3*nm) / 5);
	(3*nc^2 - 3*(n+1)*nc + 2*nm^2 + nm) / div;
}


### Component 4: Cubic

poly.cp4 = function(n) {
	dn  = if(n %% 2 == 1) -1 else 1;
	cc  = poly.c4c(n, dn=dn);
	div = (n^7 - 14*n^5 + n^3*49 - 36*n) / 4032;
	div = sqrt(div) * 6;
	return(cc / div);
}
poly.c4c = function(n, dn = -1) {
	nc = seq(1, (n+dn)/2, by=1);
	cc = (5*nc^3 - 15*(n+1)/2*nc^2 + (6*n^2 + 15*n + 11)/2*nc +
		- (n^3 + 6*n^2 + 11*n + 6)/4);
	return(cc);
}
poly.cp4.mpfr = function(n, prec = 240) {
	if(inherits(n, "mpfr")) {
		prec = Rmpfr::getPrec(n);
		v1 = Rmpfr::mpfr(1, prec);
	} else {
		n  = Rmpfr::mpfr(n, prec);
		v1 = Rmpfr::mpfr(1, prec);
	}
	isOdd = (n %% 2 == 1);
	nm = if(isOdd) (n+1)/2 else n / 2;
	nc = seq(v1, nm, by = 1);
	cc = (5*nc^3 - 15*(n+1)/2*nc^2 + (6*n^2 + 15*n + 11)/2*nc +
		- (n^3 + 6*n^2 + 11*n + 6)/4);
	div = (n^7 - 14*n^5 + n^3*49 - 36*n) / 4032;
	div = sqrt(div) * 6;
	if(isOdd) cc[length(cc)] = mpfr(0, prec);
	return(cc / div);
}

poly.cp5 = function(n) {
	nc  = seq(1, (n+1)/2, by = 1);
	div = (n^9 - 30 * n^7 + 273*n^5 - 820*n^3 + 576*n) / 5184;
	cc  = (70*nc^4 - 140*(n+1) * nc^3 + 10*(9*n^2 + 21*n + 17) * nc^2 +
		- 10*(2*n^3 + 9*n^2 + 17*n + 10) * nc +
		+ (n^4 + 10*n^3 + 35*n^2 + 50*n + 24)) / 24;
	return(cc / sqrt(div));
}
poly.cp5.mpfr = function(n, prec = 240) {
	if(inherits(n, "mpfr")) {
		prec = Rmpfr::getPrec(n);
		v1 = Rmpfr::mpfr(1, prec);
	} else {
		n  = Rmpfr::mpfr(n, prec);
		v1 = Rmpfr::mpfr(1, prec);
	}
	isOdd = (n %% 2 == 1);
	nm = if(isOdd) (n+1)/2 else n / 2;
	nc = seq(v1, nm, by = 1);
	div = (n^9 - 30 * n^7 + 273*n^5 - 820*n^3 + 576*n) / 5184;
	cc  = (70*nc^4 - 140*(n+1) * nc^3 + 10*(9*n^2 + 21*n + 17) * nc^2 +
		- 10*(2*n^3 + 9*n^2 + 17*n + 10) * nc +
		+ (n^4 + 10*n^3 + 35*n^2 + 50*n + 24)) / 24;
	return(cc / sqrt(div));
}

