

### Vandermonde Matrix

vandermonde = function(x) {
	n = length(x);
	if(n < 2) return(matrix(1, ncol=n, nrow=n));
	m = cbind(1, x);
	if(n == 2) return(m);
	xn = x;
	for(id in seq(2, n-1)) {
		xn = xn * x;
		m  = cbind(m, xn)
	}
	colnames(m) = paste0("x", seq(0, n-1));
	return(m);
}

