

### Composition / "Convolution"

### Sum(runif)
rsum = function(n, scale = 1, from = 0, to = 1) {
	n2 = 2*n;
	r = runif(n2, from, to);
	r = (r[1:n] + scale * r[(n+1):n2]) / (scale + 1);
	return(r);
}
# alternative: fraction * runif + (1-fraction)*runif;

hist.rsum = function(n, scale = 1, bins = 32, from = 0, to = 1) {
	x = rsum(n, scale=scale, from=from, to=to);
	hist(x, breaks = bins);
	invisible(x);
}

# Note:
# - see also Moments.Stat.R;

############

### Examples

# Triangular Distribution
hist.rsum(10000)

hist.rsum(10000, scale = 1/2)

hist.rsum(10000, scale = 1/4)

# - Resembles more and more a unif distribution
#   as scale -> 0 or scale -> Inf;


###############

# TODO:

n = 10000
#
x = runif(n, 0, 3)
x = (x + runif(n, 1, 2)) / 5;
hist(x, breaks = 32)

###
n = 10000
#
x = runif(n, 0, 3)
x = (x + 2*runif(n, 1, 2)) / 7;
hist(x, breaks = 32)

###
n = 10000
#
x = runif(n, 0, 3)
x = (2*x + runif(n, 1, 2)) / 8;
hist(x, breaks = 32)

###
n = 10000
#
x = runif(n, 0, 3)
x = (x + 2*runif(n, 1, 2) - 2) / 5;
hist(x, breaks = 32)

