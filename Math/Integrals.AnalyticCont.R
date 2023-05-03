

### Analytic Continuation


# For an Introduction, see link to Digamma Function:
# - Lines That Connect: Extending the Harmonic Numbers to the Reals
#   https://www.youtube.com/watch?v=9p_U_o1pMKo



####################

### Helper Functions

### sum( n / exp(n) )
invExpSum = function(x, N = 30) {
	id = seq(0, N);
	sapply(x, \(x) sum(id/exp(id)) - sum((x+1+id)/exp(x+1+id)));
}
invExpSum0 = function(n) {
	id = seq(0, n);
	sum(id/exp(id));
}


### sum( 1 / (n*exp(n) + 1) )
invXExpSum = function(x, k=1, N = 30) {
	id = seq(0, N);
	sapply(x, function(x) {
		sum(1/(id*exp(id) + k)) - sum(1/((x+1+id)*exp(x+1+id) + k))
	});
}
invXExpSum0 = function(n, k=1) {
	id = seq(0, n);
	sum(1 / (id*exp(id) + k) );
}

#####################

### Examples

### sum( n / exp(n) )
curve(invExpSum(x), xlim=c(0,10))
tmp = sapply(seq(0, 5), \(x) points(x, invExpSum0(x), col="red"))


#####################

### sum( 1 / (n*exp(n) + 1) )
curve(invXExpSum(x), xlim=c(0,10))
tmp = sapply(seq(0, 6), \(x) points(x, invXExpSum0(x), col="red"))
#
k = 1 + 1/16
curve(invXExpSum(x, k=k), add = TRUE)
tmp = sapply(seq(0, 6), \(x) points(x, invXExpSum0(x, k=k), col="red"))


# NOT identical: still some diff?
invXExpSum0(50)
# exp(1)/2 # close, but ???
integrate(\(x) 1/(x*exp(x) + 1), 0, Inf)

