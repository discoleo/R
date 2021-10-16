########################
###
### Leonard Mada
### [the one and only]
###
### Clustering: Tools & Simulations
###
### draft v.0.1d



###############
### History ###
###############


### draft v.0.1d:
# - rmatrix.sigma(): generate random sigma matrices;
### draft v.0.1c:
# - plot.cluster.3D();
### draft v.0.1a - v.0.1b-ex2:
# - generate random clusters;
# - parameter: set variance;
# - more examples; [v.0.1b-ex & v.0.1b-ex2]


####################
####################

### Helper functions

library(LaplacesDemon)
library(ggplot2)
# 3D scatterplot
library(rgl)

### Other:

rmeans = function(n, dim=2, lim=c(0, 10)) {
	len = n * dim;
	m = runif(len, lim[1], lim[2]);
	matrix(m, nrow=dim, ncol=n);
}
rcluster = function(n, cl, mu=NULL, dim=2, sigma=NULL, add.id=TRUE,
		seed=NULL, id.offset=0) {
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
	if( ! is.null(seed)) print("Warning: seed not yet implemented!")
	x = rmvn(n[[1]], mu[,1], sigma[[1]]);
	if(add.id) x = cbind(x, ID = rep(1 + id.offset, n[[1]]));
	if(ncl == 1) return(names.f(as.data.frame(x)));
	# many clusters:
	x = list(x);
	for(id in seq(2, ncl)) {
		x0 = rmvn(n[[id]], mu[,id], sigma[[id]]);
		if(add.id) x0 = cbind(x0, ID = rep(id + id.offset, n[[id]]));
		x[[id]] = x0;
	}
	x = do.call(rbind, x);
	x = as.data.frame(x);
	return(names.f(x));
}

### Matrix Operations
# Generate Covariance Matrix
matrix.sigma = function(triang, diag=1) {
	len = length(triang);
	n = (-1 + sqrt(1 + 8*len))/2;
	n = round(n);
	n = n + 1;
	### Matrix
	m = diag(diag, n);
	m[lower.tri(m)] = triang;
	m[upper.tri(m)] = t(m)[upper.tri(m)];
	return(m)
}
rmatrix.sigma = function(d=c(1), dim=2, sc=0.9) {
	if(is.list(sc)) {
		# assumes correct dimensions
	} else {
		len.up = (dim-1)*dim/2;
		if(length(sc) == 1) {
			sc = rep(list(c(-sc, sc)), len.up);
		} else if(length(sc) == len.up) {
			sc = lapply(sc, function(sc) c(-sc, sc));
		} else stop("Length of scaling-parameter NOT supported!")
	}
	# Matrix Sigma:
	rcov.f = function(sc, id) runif(1, sc[1]*d[[id]], sc[2]*d[[id]]);
	lapply(seq_along(d), function(id) {
		rndcov = sapply(sc, rcov.f, id=id);
		matrix.sigma(rndcov, diag=d[[id]]);
	})
}

### Graphic
plot.cluster.2D = function(x) {
	x$ID = as.factor(x$ID);
	ggplot(x, aes(x=v1, y=v2, fill=ID, col=ID)) +
		geom_point();
}
plot.cluster.3D = function(x, radius=0.2, col=NULL, add=FALSE) {
	if(is.null(col)) {
		col = length(unique(x$ID));
		col = colorRampPalette(c("#F0F020", "#F080F0", "#20F0F0"), space = "rgb")(col);
	}
	col = col[x$ID];
	plot3d(x[, 1:3], type="s", radius=radius, col=col, add=add);
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

plot.cluster.2D(x)


### Ex 2:
sigma = matrix.sigma(0.7, diag=1)
x = rcluster(100, cl=5, sigma=sigma)

plot.cluster.2D(x)


### Ex 3:
cl = 5
sigma = lapply(seq(cl), function(id) matrix.sigma(runif(1, -0.5, 0.9)))
x = rcluster(100, cl=cl, sigma=sigma)

plot.cluster.2D(x)


### Ex 4:
cl = 5
sdsq = 2;
set.seed(31); # 31, 35
sigma = lapply(seq(cl), function(id) matrix.sigma(runif(1, 0, 1.9), diag=sdsq))
x = rcluster(100, cl=cl, sigma=sigma)

plot.cluster.2D(x)


### Ex 5:
cl = 5
sdsq = 2;
set.seed(31); # 31, 35
sigma = lapply(seq(cl), function(id) matrix.sigma(runif(1, -1.9, 1.9), diag=sdsq))
x = rcluster(100, cl=cl, sigma=sigma)

plot.cluster.2D(x)


### Ex 6:
cl = 10
sdsq = seq(0.25, by=0.25, length.out=cl)
set.seed(35); # 35, 7
sigma = rmatrix.sigma(d=sdsq, sc=0.9, dim=2)
x = rcluster(100, cl=cl, sigma=sigma)

plot.cluster.2D(x)


################

cl = 8
sdsq = seq(0.25, by=0.25, length.out=cl)
set.seed(7);
sigma = rmatrix.sigma(d=sdsq, sc=c(0.9, 0.5, 0.5), dim=3)
x = rcluster(100, cl=cl, dim=3, sigma=sigma)

plot.cluster.3D(x)


################
################

### TODO:
# - Clustering;

