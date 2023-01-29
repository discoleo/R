########################
###
### Leonard Mada
### [the one and only]
###
### Special Matrices
### Quasi-Vandermonde
###
### draft v.0.1a


### Quasi-Vandermonde Matrix
# - same determinant as a Vandermonde matrix;

# TODO: proper sign;
matrix.quasiVand = function(x) {
	allCoef = function(x) {
		len = length(x);
		r = sapply(seq(2, len), function(k) {
			sum(combn(x, k, FUN=prod));
		})
		r = c(1, sum(x), r);
	}
	len = length(x);
	r = array(0, c(len, len));
	for(id in seq(len)) {
		r[,id] = allCoef(x[-id]);
	}
	return(r);
}
det.Vand = function(x) {
	prod(combn(x, 2, FUN=diff));
}


### Examples:

###
x = sqrt(c(2,3,5,7))
#
m = matrix.quasiVand(x)
det(m)
det.Vand(x)


###
x = sqrt(c(2,3,5,7,11))
#
m = matrix.quasiVand(x)
det(m)
det.Vand(x)


###
x = sqrt(c(2,3,5,7,11,13))
# TODO: sign
m = matrix.quasiVand(x)
det(m)
det.Vand(x)


###
x = sqrt(c(2,3,5,7,11,13, 17))
# TODO: sign
m = matrix.quasiVand(x)
det(m)
det.Vand(x)

