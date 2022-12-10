########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Multi-Variable Polynomials
### Polynomial Matrices


# this file:
# source("Polynomials.Helper.Matrix.R")


#########################

print.mpm = function(p) {
	ch = sapply(p, function(p) if(is.numeric(p)) p else as.character.pm(p));
	n  = attr(p, "dim");
	matrix(ch, nrow=n[[1]], ncol=n[[2]]);
}


### Upper-Multi-Diagonal matrix
diag.lpm = function(p, n) UseMethod("diag.lpm");
diag.lpm.default = function(p, n) {
	len = n*n;
	pR  = as.list(rep(0, len));
	m = matrix(seq(len), nrow=n, ncol=n);
	for(id in seq(length(p))) {
		for(nr in seq(n)) {
			nc = ((nr + id - 2) %% n) + 1;
			pR[[m[nr, nc]]] = p[[id]];
		}
	}
	attr(pR, "dim") = c(n, n);
	class(pR) = c("mpm", class(pR));
	return(pR);
}
diag.lpm.character = function(p, n, ...) {
	pR = lapply(p, toPoly.pm, ...);
	return(diag.lpm.default(pR, n=n));
}

# - naive implementation;
det.mpm = function(p, verbose=FALSE) {
	n = attr(p, "dim");
	if(is.null(n)) stop("Not a polynomial matrix!");
	if(length(n) != 2 && n[1] != n[2]) stop("Improper Matrix!");
	#
	n = n[1];
	if(n == 1) return(p);
	#
	if(n < 5) {
		det2 = function(p) {
			diff.pm(mult.pm(p[[1]], p[[4]]), mult.pm(p[[2]], p[[3]]));
		}
		if(n == 2) return(det2(p));
		#
		det3 = function(p) {
			if(verbose) print("Started D3");
			pR = if(is.zero.pm(p[[1]])) 0 else mult.pm(p[[1]], det2(p[c(5,6,8,9)]));
			if( ! is.zero.pm(p[[2]])) {
				tmp = det2(p[c(4,5,7,8)]);
				pR = diff.pm(pR, mult.pm(p[[2]], det2(p[c(4,6,7,9)])));
			}
			if( ! is.zero.pm(p[[3]])) {
				pR = sum.pm(pR, mult.pm(p[[3]], det2(p[c(4,5,7,8)])));
			}
			return(pR);
		}
		if(n == 3) return(det3(p));
		#
		len = n*n;
		id.m = matrix(seq(len), ncol=n);
		pR = 0;
		for(nr in seq(4)) {
			if( ! is.zero.pm(p[[nr]])) {
				tmp = mult.pm(p[[nr]], det3(p[id.m[-nr, -1]]));
				if(nr %% 2 == 0) { pR = diff.pm(pR, tmp); }
				else pR = sum.pm(pR, tmp);
			}
		}
		return(pR);
	}
	len = n*n;
	id.m = matrix(seq(len), nrow=n, ncol=n);
	pR = 0;
	# TODO: find optimal column;
	for(nr in seq(n)) {
		if( ! is.zero.pm(p[[nr]])) {
			if(verbose) print(paste0("M", n, ", Row: ", nr));
			tmp = p[id.m[-nr, -1]];
			attr(tmp, "dim") = c(n-1, n-1);
			tmp = mult.pm(p[[nr]], det.mpm(tmp, verbose=verbose));
			if(nr %% 2 == 0) { pR = diff.pm(pR, tmp); }
			else pR = sum.pm(pR, tmp);
		}
	}
	return(pR);
}

