########################
###
### Leonard Mada
### [the one and only]
###
### Polygon Process
###
### draft v.0.1a


### "Polygon"-Process

# TODO:
# https://search.r-project.org/CRAN/refmans/spatstat.random/html/00Index.html


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

### Solve/Optimize: Polygon-Angles

# based on LastPx == Origin = c(0,0);
polygonOrigin = function(a, d) {
	d = x; len = length(d);
	x1 = d[1]; y1 = 0;
	# len1 = len - 1; len2 = len1 %/% 2;
	# a = c(0, a[1] * seq(1, len2),
	#	a[1]*len2 + a[2] * seq(1, len2 + 1));
	a = seq.mult(a, len);
	for(i in seq(2, length(d))) {
		x1 = x1 + d[i]*cos(a[i]);
		y1 = y1 + d[i]*sin(a[i]);
	};
	return(c(x1, y1));
}
# based on Dist(LastPx, c(0,0));
polygonOptim = function(a, d) {
	d = x; len = length(d);
	x1 = d[1]; y1 = 0;
	# len1 = len - 1; len2 = len1 %/% 2;
	# a = c(0, a[1] * seq(1, len2),
	#	a[1]*len2 + a[2] * seq(1, len2 + 1));
	a = seq.mult(a, len);
	for(i in seq(2, len)) {
		x1 = x1 + d[i]*cos(a[i]);
		y1 = y1 + d[i]*sin(a[i]);
	};
	d = x1^2 + y1^2;
	return(d);
}

### Transformations

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

###################

### Ex 1:
x = runif(10, 1, 3)
plot.dp(x, c(0.9, 0.5))

a1 = multiroot(polygonOrigin, c(0.9, 0.5), d=x)
a2 = optim(c(0.9, 0.5), polygonOptim)

print(a1); cat("=========\n"); print(a2);

plot.dp(x, a1$root)
lines.dp(x, a2$par, col="red", lty=4, lwd=2)


### Ex 2:
x = runif(20, 1, 3)
x0 = c(0.3, 0.3)
plot.dp(x, x0)

# multiroot: sensible to x0!
a1 = multiroot(polygonOrigin, x0, d=x)
a2 = optim(x0, polygonOptim)

print(a1); cat("=========\n"); print(a2);

plot.dp(x, a2$par)
lines.dp(x, a1$root, col="red", lty=4, lwd=2)

