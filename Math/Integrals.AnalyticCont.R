

### Analytic Continuation


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
invXExpSum = function(x, N = 30) {
	id = seq(0, N);
	sapply(x, function(x) {
		sum(1/(id*exp(id) + 1)) - sum(1/((x+1+id)*exp(x+1+id) + 1))
	});
}
invXExpSum0 = function(n) {
	id = seq(0, n);
	sum(1 / (id*exp(id) + 1) );
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

# NOT identical: still some diff?
invXExpSum0(50)
# exp(1)/2 # close, but ???
integrate(\(x) 1/(x*exp(x) + 1), 0, Inf)

