########################
###
### Leonard Mada
### [the one and only]
###
### Clustering: Tools & Simulations
###
### draft v.0.1a



###############
### History ###
###############


### draft v.0.1a:
# - generate random clusters;


####################
####################

### Helper functions

library(LaplacesDemon)
library(ggplot2)

### Other:

rmeans = function(n, dim=2, lim=c(0, 10)) {
	len = n * dim;
	m = runif(len, lim[1], lim[2]);
	matrix(m, nrow=dim, ncol=n);
}
rcluster = function(n, cl, mu=NULL, dim=2, sigma=NULL, add.id=TRUE, seed=NULL) {
	len = length(n);
	ncl = if(len > 1) len else cl;
	if( ! missing(cl) && ncl != cl) stop("Number of clusters must match!");
	if(len == 1 && ncl > 1) n = rep(n, ncl);
	### Parameters
	if(is.null(mu)) {
		mu = rmeans(ncl, dim=dim);
	}
	# Covariance:
	if(is.null(sigma)) sigma = diag(dim);
	if( ! is.list(sigma)) {
		sigma = list(sigma);
		sigma = rep(sigma, ncl);
	}
	### Names
	names.f = function(x) {
		nms = paste0("v", seq(dim));
		if(add.id) nms = c(nms, "ID");
		names(x) = nms;
		return(x);
	}
	### Samples
	# TODO: seed;
	x = rmvn(n[[1]], mu[,1], sigma[[1]]);
	if(add.id) x = cbind(x, ID = rep(1, n[[1]]));
	if(ncl == 1) return(names.f(as.data.frame(x)));
	# many clusters:
	x = list(x);
	for(id in seq(2, ncl)) {
		x0 = rmvn(n[[id]], mu[,id], sigma[[id]]);
		if(add.id) x0 = cbind(x0, ID = rep(id, n[[id]]));
		x[[id]] = x0;
	}
	x = do.call(rbind, x);
	x = as.data.frame(x);
	return(names.f(x));
}

### Matrix Operations
# Generate Covariance Matrix
matrix.sigma = function(triang) {
	len = length(triang);
	n = (-1 + sqrt(1 + 8*len))/2;
	n = round(n);
	n = n + 1;
	### Matrix
	m = diag(n);
	m[lower.tri(m)] = triang;
	m[upper.tri(m)] = t(m)[upper.tri(m)];
	return(m)
}

####################
####################

### Ref:
# 1.) Yudong Chen. Structures Of Local Minima In K-Means And Mixture Models.
#     Simons Institute.
#     https://www.youtube.com/watch?v=4Smq8JX12-8


### Examples

### Ex 1:
x = rcluster(100, cl=5)
x$ID = as.factor(x$ID)

ggplot(x, aes(x=v1, y=v2, fill=ID, col=ID)) +
	geom_point()


### Ex 2:
sigma = matrix(c(1, 0.7, 0.7, 1), nrow=2, ncol=2)
x = rcluster(100, cl=5, sigma=sigma)
x$ID = as.factor(x$ID)

ggplot(x, aes(x=v1, y=v2, fill=ID, col=ID)) +
	geom_point()


### Ex 3:
cl = 5
sigma = lapply(seq(cl), function(id) matrix.sigma(runif(1, 0, 0.9)))
x = rcluster(100, cl=cl, sigma=sigma)
x$ID = as.factor(x$ID)

ggplot(x, aes(x=v1, y=v2, fill=ID, col=ID)) +
	geom_point()


################
################

### TODO:
# - Clustering;

