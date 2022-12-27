########################
###
### Leonard Mada
### [the one and only]
###
### Matrices:
### Tools for Matrices
###
### draft v.0.1b


### Tools for Matrices


### Symmetric Diagonal / Band

# Cyclic Permutation:
diag.band = function(x, n=length(x)) {
	len = length(x);
	if(n < len) stop("Invalid dimension!");
	if(n > len) {
		x = c(x, rep(0, n - len));
	}
	m = matrix(0, nrow=n, ncol=n);
	xr = rev(x);
	for(nc in seq(n - 1)) {
		v = c(xr[seq(n - nc + 1, n)], xr[seq(1, n - nc)]);
		m[, nc] = v;
	}
	m[, n] = xr;
	return(m);
}

# Class 2 Poly => Eigenvalues
roots.class2 = function(x, n=length(x)) {
	u = cos(2*pi/n) + 1i*sin(2*pi/n);
	idNZ = which(x[-1] != 0);
	u = c(1, u^idNZ);
	r = sapply(seq(0, n-1), function(id) {
		sum(x * u^id);
	})
	return(r);
}


### Complex Matrices

reduce.cm = function(m, mult=1, div=1) {
	nonzero = which(m[,1] != 0);
	if(length(nonzero) == 0) return(list(m=m, mult=mult, div=div));
	if(nonzero[1] > 1) {
		tmp = m[1,]; m[1,] = m[nonzero[1],]; m[nonzero[1],] = tmp;
		div = if(nonzero[1] %% 2) div else - div;
	}
	nonzero = nonzero[-1];
	for(nr in nonzero) {
		e = m[nr, 1];
		div = div * m[1,1];
		m[nr,] = m[nr,]*m[1,1] - e*m[1,];
	}
	return(list(m=m, mult=mult, div=div));
}

det.cm = function(m, mult=1, div=1) {
	nc = ncol(m);
	if(nc <= 3) {
		if(nc == 1) {
			d = m[1];
		} else if(nc == 2) {
			d = m[1]*m[4] - m[2]*m[3];
		} else {
			d = m[1]*m[5]*m[9] + m[2]*m[6]*m[7] + m[4]*m[8]*m[3] +
				- m[3]*m[5]*m[7] - m[2]*m[4]*m[9] - m[1]*m[6]*m[8];
		}
		return(list(det = d, div = div/mult));
	}
	m.all = reduce.cm(m, mult=mult, div=div);
	m = m.all$m; m1 = m[1,1];
	if(m1 == 0) return(list(det=0, mult=mult, div=div));
	mult = m.all$mult * m1; div = m.all$div;
	return(det.cm(m[-1, -1], mult = mult, div=div));
}

#########################
#########################

#############
### Eigen ###
#############

# - using roots of Class 2 Polynomials;

###
x = c(1,2,3)
m = diag.band(x, n=5)
roots.class2(x, 5)
eigen(m)$values

###
x = c(1,2,3,2)
m = diag.band(x, n=5)
roots.class2(x, 5)
eigen(m)$values

###
x = c(1,2,3,2)
m = diag.band(x, n=7)
roots.class2(x, 7)
eigen(m)$values

###
x = c(1,2,4,3,1)
m = diag.band(x, n=7)
roots.class2(x, 7)
eigen(m)$values


####################

### Test Complex Det

n = 5
x = sample(seq(-5, 5), n^2, T)
# m = matrix(as.numeric(x), ncol=n)
m = matrix(x, ncol=n)
det(m)
d = det.cm(m)
d$det / d$div


###
n = 5
m = matrix(as.numeric(sample(seq(-5, 5), n^2, T)), ncol=n)
x = sample(seq(-5, 5), n^2, T);
m = m + 1i * matrix(x, ncol=n)
d = det.cm(m)
d$det / d$div
# workaround
det(m) # problem in R
prod(eigen(m)$values)


###
n = 5
m = matrix(as.numeric(sample(seq(-5, 5), n^2, T)), ncol=n)
m = m + 1i * matrix(as.numeric(sample(seq(-5, 5), n^2, T)), ncol=n)
d = det.cm(m)
d$det / d$div
# workaround
det(m) # problem in R
prod(eigen(m)$values)

