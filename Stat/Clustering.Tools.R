########################
###
### Leonard Mada
### [the one and only]
###
### Clustering: Tools & Simulations
###
### draft v.0.1g



###############
### History ###
###############


### draft v.0.1g:
# - density contours;
### draft v.0.1e - v.0.1f-3D:
# - cluster around a polygon;
# - more complicated examples; [v.0.1f & v.0.1f-ex2]
# - 3D example; [v.0.1f-3D]
### draft v.0.1d:
# - rmatrix.sigma(): generate random cov (sigma) matrices;
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
		# Group ID:
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
	if(is.matrix(d)) d = lapply(seq(ncol(d)), function(nc) d[,nc]);
	isList = FALSE; # big list with scalings per d;
	if(is.list(sc)) {
		# assumes correct dimensions
		isList = TRUE;
	} else {
		len.up = (dim-1)*dim/2;
		len = length(d);
		if(length(sc) == 1) {
			sc = rep(list(c(-sc, sc)), len.up);
		} else if(length(sc) == len.up) {
			sc = lapply(sc, function(sc) c(-sc, sc));
		} else if(length(sc) == len.up*len) {
			# TODO: requires another lapply; or different concept ???
			sc = lapply(seq(len), function(id) sort(c(-sc[id], sc[id])));
			isList = TRUE;
		} else stop("Length of scaling-parameter NOT supported!")
	}
	# Matrix Sigma:
	rcov.f = function(sc, id) runif(1, sc[1]*d[[id]][1], sc[2]*d[[id]][1]);
	lapply(seq_along(d), function(id) {
		sc = if(isList) sc[id] else sc;
		rndcov = sapply(sc, rcov.f, id=id);
		matrix.sigma(rndcov, diag=d[[id]]);
	})
}

### Graphic
plot.cluster.2D = function(x, contour.add=TRUE) {
	x$ID = as.factor(x$ID);
	g = ggplot(x, aes(x=v1, y=v2, fill=ID, col=ID)) +
		geom_point();
	# TODO: Gaussian-contours;
	if(contour.add) g = g + geom_density_2d();
	return(g);
}
plot.cluster.3D = function(x, radius=0.2, col=NULL, add=FALSE) {
	if(is.null(col)) {
		col = length(unique(x$ID));
		col = colorRampPalette(c("#F0F020", "#F080F0", "#20F0F0"), space = "rgb")(col);
	}
	col = col[x$ID];
	plot3d(x[, 1:3], type="s", radius=radius, col=col, add=add);
}

### Shapes / Geometry

polygon.reg = function(n, r=1, a.offset = 0, clockwise=FALSE, closed=FALSE) {
	p = 2*pi / n;
	p = if(clockwise) p * c(0, seq(n-1, 1, by=-1)) else p * seq(0, n-1);
	if(closed) p = c(p, 0);
	p = p + a.offset;
	px = cos(p); py = sin(p);
	px = r * px; py = r * py;
	pxy = cbind(x=px, y=py);
	# TODO: round0();
	return(pxy);
}

####################
####################

### Ref:
# 1) Yudong Chen. Structures Of Local Minima In K-Means And Mixture Models.
#    Simons Institute. https://www.youtube.com/watch?v=4Smq8JX12-8
# 2) Braxton Osting. Archetypal Analysis.
#    Simons Institute. https://www.youtube.com/watch?v=FbP568yWp7U
#  - regularization; image clustering;


### Examples

### Ex 1:
x = rcluster(100, cl=5)

plot.cluster.2D(x)


### Density-Contours:

library(car)
library(cluster)
library(ellipse)

# Note: cannot be mixed with ggplot!

### car-package:
# [Base-R Graphics]
dataEllipse(x$v1, x$v2, as.factor(x$ID), levels=c(0.5, 0.75, 0.95))

cov.m = cov.wt(data.frame(x$v1, x$v2)[x$ID == 1, ]);
# MASS::cov.rob(data.frame(x$v1, x$v2)[x$ID == 1, ]);

### cluster-package:
n.obs = sum(x$ID == 1);
qVal = qchisq(0.99, df=2); # qf(0.99, 2, n.obs - 1)
lines(ellipsoidPoints(cov.m$cov, qVal, loc=cov.m$center), col="green")

### ellipse-package
sc = diag(cov.m$cov); # sc = c(1,1);
lines(ellipse(cov.m$cov, scale=sc, centre=cov.m$center, level=c(0.90)), col="red")


# explicitly construct the Ellipse:
Q = chol(cov.m$cov, pivot=TRUE);
id = order(attr(Q, "pivot"));
Q[, id]
# TODO

#########

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

### Polygons

### Ex 7:
cl = 6
r = 3.5
mu = t(polygon.reg(cl, r=r));
sdsq  = rep(1, cl)
sigma = rmatrix.sigma(d=sdsq, sc=0, dim=2)
x = rcluster(100, cl=cl, mu=mu, sigma=sigma)

plot.cluster.2D(x)


rspikes = function(n, cl, r=1, id.offset=0, dim=2) {
	mu = t(polygon.reg(cl, r=r));
	if(dim == 3) mu = rbind(mu, rep(0, cl));
	# only regular Hexagon;
	# TODO: generalize;
	scale.h = exp((1+1)/2); # exp(sqrt(2))
	sdsq = rbind(
		c(scale.h,1/2,1/2, scale.h,1/2,1/2),
		c(0.075,3/2,3/2, 0.075,3/2,3/2));
	sc = c(0,0.7,-0.7, 0,0.7,-0.7);
	if(dim == 3) {
		sdsq = rbind(sdsq, rep(0.1, cl));
		sc   = lapply(sc, function(sc) c(sc, 0, 0));
	}
	sigma = lapply(seq(cl), function(id) matrix.sigma(sc[[id]], d=sdsq[,id]))
	x = rcluster(n, cl=cl, mu=mu, sigma=sigma, dim=dim, id.offset=id.offset);
	return(x);
}
### Ex 8:
cl = 6
r = 3.5
# Set 1:
mu = t(polygon.reg(cl, r=r, a.offset = pi/cl));
sdsq  = rep(1, cl)
sigma = rmatrix.sigma(d=sdsq, sc=0, dim=2)
x = rcluster(100, cl=cl, mu=mu, sigma=sigma)
# Set 2:
x2 = rspikes(100, cl=cl, r=2*r, id.offset=cl)
#
x = rbind(x, x2);

plot.cluster.2D(x)


### Ex 9:
cl = 6
r = 4.25
# Set 1:
x = rspikes(100, cl=cl, r=r, id.offset=cl)
# Set 2:
mu = t(polygon.reg(cl, r=1.5*r, a.offset = pi/cl));
sdsq  = rep(1, cl)
sigma = rmatrix.sigma(d=sdsq, sc=0, dim=2)
x2 = rcluster(100, cl=cl, mu=mu, sigma=sigma)
#
x = rbind(x, x2);

plot.cluster.2D(x)


################
################

### 3D

cl = 8
sdsq = seq(0.25, by=0.25, length.out=cl)
set.seed(7);
sigma = rmatrix.sigma(d=sdsq, sc=c(0.9, 0.5, 0.5), dim=3)
x = rcluster(100, cl=cl, dim=3, sigma=sigma)

plot.cluster.3D(x)


### Ex 2:
cl = 6
r = 4.25
# Set 1:
x = rspikes(100, cl=cl, r=r, dim=3, id.offset=0)

plot.cluster.3D(x)


################
################

### TODO:
# - Clustering;
# - Mixture Models;

### packages:
# - mixreg, mixdist;

