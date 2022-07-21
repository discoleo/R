

### Roots of Unity
unity = function(n=3, all=TRUE) {
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	if(all) {
		m = m^(0:(n-1))
	}
	return(m)
}


### Real roots
isRe.f = function(x) {
	lapply(x, function(x) (Im(x) == 0));
}
Re.f = function(x, isRe=NULL) {
	if( ! is.null(isRe) && is.list(isRe) && ! is.list(x)) {
		x = lapply(seq_along(isRe), function(id) x);
	}
	if(is.null(isRe)) isRe = isRe.f(x);
	x = lapply(seq_along(x), function(id) Re(x[[id]][isRe[[id]]]));
	return(x);
}
range.c = function(x, isRe=NULL) {
	if(is.null(isRe)) isRe = isRe.f(x);
	x1 = Re(x[[1]][isRe[[1]]]);
	x2 = Re(x[[2]][isRe[[2]]]);
	xmax = max(x1, x2)
	xmin = min(x1, x2)
	return(list(rg=c(xmin, xmax), isRe=isRe));
}

