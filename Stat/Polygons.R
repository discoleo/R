########################
###
### Leonard Mada
### [the one and only]
###
### Polygon Process
###
### draft v.0.1n


### "Polygon"-Process


# this file:
# source("Polygons.R")


# TODO:
# - distribute randomly using a point process, see:
# https://search.r-project.org/CRAN/refmans/spatstat.random/html/00Index.html


### Libraries

library(shape)
library(rootSolve)

### Plot

circle = function(r, mid, col=1, fill=NULL, ...) {
	filledcircle(r1=r, mid=mid, lcol=col, col=fill, ...);
}

plot.dp = function(x, a, type="l", ...) {
	xy = coords.dp(x, a=a)
	plot(xy[,1], xy[,2], type=type, ...);
	invisible(xy);
}
lines.dp = function(x, a, type="l", ...) {
	xy = coords.dp(x, a=a)
	lines(xy[,1], xy[,2], type=type, ...);
	invisible(xy);
}

seq.mult = function(x, length) {
	lenX = length(x);
	if(lenX == 1) {
		v = c(0, x[1] * seq(length - 1));
	} else {
		each = diff(c(1, ceiling(quantile(seq(length), seq(lenX)/lenX))));
		A0 = x[- lenX] * each[- lenX];
		S0 = cumsum(A0);
		v  = unlist(lapply(seq(2, lenX), function(id) {
			seq(each[id]) * x[id] + S0[id - 1];
		}));
		v = c(0, seq(each[1]) * x[1], v);
	}
	return(v);
}
coords.dp = function(x, a) {
	d = x; len = length(d);
	x = rep(0, len);
	y = rep(0, len)
	x[2] = d[1]; y[2] = 0;
	a = seq.mult(a, len);
	for(i in seq(2, length(d))) {
		x[i + 1] = x[i] + d[i]*cos(a[i]);
		y[i + 1] = y[i] + d[i]*sin(a[i]);
	};
	return(cbind(x, y));
}

plot.ini = function(xlim, ylim=xlim, ...) {
	plot.new();
	plot.window(xlim=xlim, ylim=ylim, ...);
	Axis(side=1); Axis(side=2);
}

### Generators

# TODO:
# rtriangle:
# - .sas, .dist, .area, .circle, .incircle;

rtriangle.circle = function(n, ..., r=1,
		type=c("sequential", "random", "half", "phalf", "eqpart"),
		center=c(0,0), asX = FALSE, tol=1E-4) {
	type = match.arg(type);
	if(type == "sequential") {
		# if a1 = rnorm(tol, 2*pi-3*tol)
		# => skewed towards almost-degenerate triangles;
		a1 = rnorm(n, tol, pi);
		dA = 2*pi - a1 - 2*tol;
		dA = pmax(tol, dA);
		a2 = rnorm(n, tol, dA);
		if(length(r) == 1) r = rep(r, n);
		xy = lapply(seq(n), function(id) {
			as.triangle.circle(c(0, a1[id], a2[id]), r=r[id], type="sequential",
				center=center, asX=asX);
		});
		return(xy);
	}
	if(type == "half" || type == "phalf") {
		a1 = rnorm(n, tol, pi);
		a2 = rnorm(n, tol, pi - tol);
		if(type == "phalf") a2 = a2 + pi;
		if(length(r) == 1) r = rep(r, n);
		type = if(type == "half") "sequential" else "random";
		xy = lapply(seq(n), function(id) {
			as.triangle.circle(c(0, a1[id], a2[id]), r=r[id], type=type,
				center=center, asX=asX);
		});
		return(xy);
	}
	if(type == "eqpart") {
		opt = list(...);
		mult = opt$mult;
		if(is.null(mult)) mult = 2;
		id = sample(seq(mult*n), n);
		a1 = id * pi / (mult*n);
		# 2nd Partition:
		part = opt$part;
		if(is.null(part)) part = max(5, round(sqrt(n)));
		id = sample(seq(part), n, replace=TRUE);
		a2 = pi * (1 - 1/(2*sqrt(n)) + id / (part + 1.5));
		if(length(r) == 1) r = rep(r, n);
		xy = lapply(seq(n), function(id) {
			as.triangle.circle(c(0, a1[id], a2[id]), r=r[id], type="random",
				center=center, asX=asX);
		});
		return(xy);
	}
	# Note: probably skews towards "larger" triangles,
	# including close-to-right angles;
	if(type == "random") {
		upper = 2*pi - tol;
		a1 = rnorm(n, tol, upper);
		a2 = rnorm(n, tol, upper);
		isEq = abs(a1 - a2) <= tol;
		n2 = sum(isEq);
		if(n2 > 0) {
			a2[isEq] = rnorm(n2, tol, upper);
		}
		if(length(r) == 1) r = rep(r, n);
		xy = lapply(seq(n), function(id) {
			as.triangle.circle(c(0, a1[id], a2[id]), r=r[id], type="random",
				center=center, asX=asX);
		});
		return(xy);
	}
}


# Side 1 = along OX axis;
# Note: generates only 1 triangle;
# - Issue with multiple triangles:
#   either list with coordinates for each triangle,
#   or matrix/df with an additional id-column;
as.triangle.dist = function(d, tol=1E-8) {
	# Note: does NOT check if triangle is valid;
	len = length(d);
	if(len == 0) return(array(0, c(0,2)));
	if(len != 3) stop("Incorrect number of distances!");
	dd23 = (d[2]^2 + d[3]^2 - d[1]^2);
	### Special Case: right A
	# numerical instability unlikely:
	# => probably unnecessary;
	if(abs(dd23) <= tol) {
		# d1 is along OX:
		h  = d[2]*d[3]/d[1];
		xA = (d[1]^2 + d[3]^2 - d[2]^2) / (2*d[1]);
		x = c(0, d[1], xA);
		y = c(0, 0, h);
		return(cbind(x, y));
	}
	#   A
	# B___C
	xA = (d[1]^2 + d[3]^2 - d[2]^2) / (2*d[1]);
	sq = d^2;
	yA = 2*(sq[1]*sq[2] + sq[1]*sq[3] + sq[2]*sq[3]) - sq[1]^2 - sq[2]^2 - sq[3]^2;
	yA = sqrt(yA) / (2*d[1]);
	x = c(0, d[1], xA);
	y = c(0,0,yA);
	return(cbind(x,y));
}

### Circumscribed Circle
# a = central angles: can be 2 angles or
#     3 angles with sum = 2*pi for type sequential, but NOT checked;
# r = radius of circle;
# mid = c(0,0): convenience parameter;
# asX = normalize side 1 parallel to OX;
as.triangle.circle = function(a, r=1, type=c("sequential", "random"),
		center=c(0,0), asX = FALSE) {
	type = match.arg(type);
	if(length(a) < 3) a = c(0, a);
	if(type == "sequential") {
		a = cumsum(a[1:3]);
	}
	if(asX) {
		dA = (pi - a[2]) / 2 - a[1];
		a = a + dA;
	}
	x = r*cos(a) + center[1];
	y = r*sin(a) + center[2];
	xy = cbind(x, y);
	return(xy);
}

# Based on the incircle
# d = side 1 (along OX);
# r = radius of incircle;
# prop = proportion of side 1 determined by incircle;
as.triangle.incircle = function(d, r, prop, tol=1E-8) {
	dB = d*prop; dC = d - dB;
	sinB = 2*dB*r / (r^2 + dB^2);
	sinC = 2*dC*r / (r^2 + dC^2);
	ds = sinB - sinC;
	if(abs(ds) < tol) {
		r2 = r^2;
		dA = 4*d*r2 / (d^2 - 4*r^2);
		x  = c(0,d,d/2);
		y  = c(0,0, r + dA);
		return(cbind(x, y));
	}
	dA = (dC*sinC - dB*sinB) / ds;
	# TODO: optimize?
	as.triangle.dist(c(d, dC + dA, dB + dA));
}


### Solve/Optimize: Polygon-Angles

# based on LastPx == Origin = c(0,0);
# d = side lengths of polygon;
# a = pi - angles;
polygonOrigin = function(a, d, id=NULL) {
	len = length(d); lenA = length(a) - 2;
	x1  = d[1]; y1 = 0;
	aa  = seq.mult(a, len);
	for(i in seq(2, len)) {
		x1 = x1 + d[i]*cos(aa[i]);
		y1 = y1 + d[i]*sin(aa[i]);
	};
	val = c(x1, y1);
	# hack:
	if(lenA > 0) {
		a = a[-c(1,2)];
		valD = x1 * a^seq(lenA) + y1 * a^seq(lenA, 1);
		val = c(val, valD);
		if( ! is.null(id)) val = val[id];
	}
	return(val);
}

# based on Dist(LastPx, c(0,0));
polygonOptim = function(a, d) {
	len = length(d);
	x1 = d[1]; y1 = 0;
	# len1 = len - 1; len2 = len1 %/% 2;
	# a = c(0, a[1] * seq(1, len2),
	#	a[1]*len2 + a[2] * seq(1, len2 + 1));
	a = seq.mult(a, len);
	for(i in seq(2, len)) {
		x1 = x1 + d[i]*cos(a[i]);
		y1 = y1 + d[i]*sin(a[i]);
	};
	d0 = x1^2 + y1^2;
	return(d0);
}

### Helper

# - checks only length of sides;
# Note:
# - a different polygon construction mechanism
#   can be implemented using this method;
is.valid.poly = function(x, degenerate=FALSE) {
	x = sort(x);
	len = length(x);
	s = sum(x[- len]);
	top = x[len];
	if(s < top) return(FALSE);
	if(s > top || (degenerate && s == top)) return(TRUE);
	return(FALSE);
}

### Triangles: Incircle
# tol = tolerance for degenerate triangle;
is.valid.incircle = function(d, r, prop, degenerate=FALSE, tol=1E-8) {
	d2 = d/2;
	if(r > d2) return(FALSE);
	if(r == d2) {
		# TODO:
		# (abs(prop - 1/2) < tol) OR (abs(r - d2) < tol); ???
		if(prop != 1/2) return(FALSE);
		return(degenerate);
	}
	dB = d*prop; dC = d - dB;
	sinB = 2*dB*r / (r^2 + dB^2);
	sinC = 2*dC*r / (r^2 + dC^2);
	ds = sinB - sinC;
	if((ds < 0 && prop < 1/2) || (ds > 0 && prop > 1/2)) {
		return(FALSE);
	} else if(abs(ds) < tol) {
		return(degenerate);
	}
	return(TRUE);
}

### Analysis

### Side Lengths
dist.triangle = function(x, sort=TRUE) {
	distf = if(sort) {
			function(xy) {
				x = xy[,1]; y = xy[,2];
				sort(c(sqrt( (x[c(1,1,2)] - x[c(2,3,3)])^2 +
					+ (y[c(1,1,2)] - y[c(2,3,3)])^2)));
			}
		} else {
			function(xy) {
				x = xy[,1]; y = xy[,2];
				c(sqrt( (x[c(1,1,2)] - x[c(2,3,3)])^2 +
					+ (y[c(1,1,2)] - y[c(2,3,3)])^2));
			}
		}
	if(inherits(x, "matrix")) {
		d = distf(x);
		class(d) = c("dist", class(d));
		return(d);
	} else if(inherits(x, "list")) {
		d = sapply(x, distf);
		d = t(d);
		class(d) = c("pdist", class(d));
		return(d);
	}
}

'[.pdist' = function(x, op1, op2) {
	tmp = unclass(x);
	tmp = tmp[op1, op2];
	isDist = missing(op2) ||
		(length(op2) == 3 && all(c(1,2,3) %in% op2));
	if(isDist) class(tmp) = c("pdist", class(tmp));
	return(tmp);
}


### Area
area.triangle = function(x, ...) {
	UseMethod("area.triangle");
}
area.triangle.dist = function(d) {
	return(area.triangle.pdist(d));
}
area.triangle.pdist = function(d) {
	if(is.null(dim(d))) d = matrix(d, nrow=1);
	s = (d[,1] + d[,2] + d[,3]) / 2;
	area = s*(s - d[,1])*(s - d[,2])*(s - d[,3]);
	area = sqrt(area);
	return(area);
}

analyse.triangle = function(n, r=1, bw = NULL) {
	par.old = par(mfrow = c(2,2));
	for(type in c("sequential", "random", "phalf", "eqpart")) {
		p = rtriangle.circle(n, r=r, type=type, asX=TRUE);
		d = dist.triangle(p);
		A = area.triangle(d);
		if(is.null(bw)) {
			Density = density(A);
		} else {
			Density = density(A, bw=bw);
		}
		plot(Density, main=type);
	}
	par(par.old);
	invisible();
}


### Transformations

### TODO:
# - rotate:
#   around one vertex, center of incircle,
#   center of circumscribed circle;

as.convex = function(x, y) {
	if(missing(y)) {
		dim = dim(x);
		if(is.null(dim)) stop("Missing y!");
		y = x[,2]; x = x[,1];
	}
	len = length(x);
	for(i in seq(len)) {
		# TODO:
	}
}


# t = proportion of sides (similar to a quasi-Bezier curve);
# r, phi, center = convenience parameters to generate regular n-gon;
# - transform can be applied twice for a Bezier-like effect;
transform.Bezier = function(t, xy=NULL, r=1, phi=0, center=c(0,0)) {
	n = length(t);
	if(is.null(xy)) {
		x = r*cos(2*seq(0, n-1)*pi/n + phi) + center[1];
		y = r*sin(2*seq(0, n-1)*pi/n + phi) + center[2];
	} else {
		x = xy[,1]; y = xy[,2];
	}
	xp = rep(0, n); yp = rep(0, n);
	id = c(seq(2, n), 1);
	xp = (1-t)*x + t*x[id];
	yp = (1-t)*y + t*y[id];
	p = cbind(x=xp, y=yp);
	class(p) = c("polygon", class(p));
	return(p);
}

# TODO:
# - transform.lattice();


###################
###################

### Examples:
# - moved to file:
#   Polygons.Examples.R;

