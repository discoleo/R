

### Analytic Continuation


# For an Introduction, see link to Digamma Function:
# - Lines That Connect: Extending the Harmonic Numbers to the Reals
#   https://www.youtube.com/watch?v=9p_U_o1pMKo



####################

### Helper Functions

### Simple ###

### Based on the Zeta-function
# - Zeta Partial-Sums;

### sum( 1 / n^p )
zetaSum = function(x, p=2, N = 50) {
	id = seq(1, N);
	sapply(x, \(x) sum(1/id^p) - sum(1/(x+id)^p));
}
zetaSum0 = function(n, p=2) {
	id = seq(1, n);
	sum(1/id^p);
}


### Exponentials ###

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
#####################

################
### Examples ###
################

### Simple

### sum( 1 / n^p )
curve(zetaSum(x), xlim=c(0,10))
tmp = sapply(seq(1, 8), \(x) points(x, zetaSum0(x), col="red"))
#
p = 3
curve(zetaSum(x, p=p), add=T, col="blue")
tmp = sapply(seq(1, 8), \(x) points(x, zetaSum0(x, p=p), col="red"))

# Q: Can we compute a closed form for zetaSum(3/2)?
zetaSum(1/2, N=200)
4 - pi^2/3
zetaSum(3/2, N=200)
4 - pi^2/3 + 4/9
zetaSum(5/2, N=200)
4 - pi^2/3 + 4/9 + 4/25


################
################

### Exponentials

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

