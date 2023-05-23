

### Analytic Continuation


# For an Introduction, see link to Digamma Function:
# - Lines That Connect: Extending the Harmonic Numbers to the Reals
#   https://www.youtube.com/watch?v=9p_U_o1pMKo



####################

### Helper Functions

### Simple ###

### Based on H(n)

### H(3*n) - H(n)
# - accuracy/convergence problems possible with some of the formulas;
#   (I may have had a bug in the code as well)
H3n = function(x, N = 30) {
	n = seq(0, N);
	sapply(x, \(x) {
		m = n + x;
		sum( (9*n+5) / ((3*n+1)*(3*n+2)*(3*n+3)) ) +
			- sum( (9*m+5) / ((3*m+1)*(3*m+2)*(3*m+3)));
	})
}
H3n0 = function(n) {
	sum(1/seq(3*n)) - sum(1/seq(n));
}
H3n0.alt = function(n) {
	n = seq(0, n-1);
	sum( (9*n+5) / ((3*n+1)*(3*n+2)*(3*n+3)) );
}


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


### Quasi-Zeta

### sum( 1 / n^n )
zetaQSum = function(x, N = 30) {
	id = seq(1, N);
	sapply(x, \(x) {
		m = id + x;
		sum(1/id^id) - sum(1/m^m);
	})
}
zetaQSum0 = function(n) {
	id = seq(1, n);
	sum(1/id^id);
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

### H(3*n) - H(n)

px = c(1,5,10,15);
curve(H3n(x), xlim=c(0, 30), ylim=c(0, 1.2))
points(px, sapply(px, H3n0), col="red")
abline(h = log(3), col="red")
# for derivation of log(3), see file:
# Sums.Trig.R;

###
H3n(1/2, N=10000)
2/3


##################
##################

### Zeta-based

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


##############

### Quasi-Zeta
# sum( 1 / n^n )

curve(zetaQSum(x), xlim=c(0, 6))
tmp = sapply(seq(1, 6), \(x) points(x, zetaQSum0(x), col="red"))

###
curve(zetaQSum(x), xlim=c(1.95, 7))
tmp = sapply(seq(2, 6), \(x) points(x, zetaQSum0(x), col="red"))

###
zetaQSum(1/2)
# TODO: ???


################

### Bernoulli Integral
# TODO:


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

