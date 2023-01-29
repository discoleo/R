########################
###
### Leonard Mada
### [the one and only]
###
### Special Matrices
### Quasi-Vandermonde
###
### draft v.0.1b-sign


### Quasi-Vandermonde Matrix
# - same determinant as a Vandermonde matrix;
# - useful for Fraction Decomposition of polynomial fractions;


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
# Direct computation:
det.Vand = function(x) {
	r = prod(combn(x, 2, FUN=diff));
	return(r);
}

# Fraction Decomposition:
coef.fr = function(x) {
	r = sapply(seq(length(x)), function(id) {
		d = det.Vand(x[-id]);
		if(id %% 2 == 0) d = -d;
		return(d);
	});
	# TODO: proper sign;
	d = det.Vand(x);
	if(length(x) %% 2 == 1) d = -d;
	return(r/d)
}


### Examples:

###
x = sqrt(c(2,3,5,7))
#
m = matrix.quasiVand(x)
det(m)
det.Vand(x)
# Fraction Decomposition
xx = 3^(1/3)
a  = coef.fr(x)
1/prod(xx - x) # ==
sum(a / (xx - x))


###
x = sqrt(c(2,3,5,7,11))
#
m = matrix.quasiVand(x)
det(m)
det.Vand(x)
# Fraction Decomposition
xx = 3^(1/3)
a  = coef.fr(x)
1/prod(xx - x) # ==
sum(a / (xx - x))


###
x = sqrt(c(2,3,5,7,11,13))
# TODO: sign
m = matrix.quasiVand(x)
det(m)
det.Vand(x)
# Fraction Decomposition
xx = 3^(1/3)
a  = coef.fr(x)
1/prod(xx - x) # ==
sum(a / (xx - x))


###
x = sqrt(c(2,3,5,7,11,13, 17))
# TODO: sign
m = matrix.quasiVand(x)
det(m)
det.Vand(x)
# Fraction Decomposition
xx = 3^(1/3)
a  = coef.fr(x)
1/prod(xx - x) # ==
sum(a / (xx - x))

