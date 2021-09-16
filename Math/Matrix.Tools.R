########################
###
### Leonard Mada
### [the one and only]
###
### Matrices:
### Tools for Matrices
###
### draft v.0.1a


### Tools for Matrices


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

###############

### Test

n = 5
m = matrix(as.numeric(sample(seq(-5, 5), n^2, T)), ncol=n)
det(m)
d = det.cm(m)
d$det / d$div


###
n = 5
m = matrix(as.numeric(sample(seq(-5, 5), n^2, T)), ncol=n)
m = m + 1i * matrix(as.numeric(sample(seq(-5, 5), n^2, T)), ncol=n)
det(m) # problem in R
d = det.cm(m)
d$det / d$div
