########################
###
### Leonard Mada
### [the one and only]
###
### Polygon Process
###
### draft v.0.1b


### "Polygon"-Process

# TODO:
# https://search.r-project.org/CRAN/refmans/spatstat.random/html/00Index.html


###

library(rootSolve)

### Plot

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
# d = side lengths of polygon;
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
a2 = optim(c(0.9, 0.5), polygonOptim, d=x)

print(a1); cat("=========\n"); print(a2);

plot.dp(x, a1$root)
lines.dp(x, a2$par, col="red", lty=4, lwd=2)


### Ex 2:
x = runif(20, 1, 3)
x0 = c(0.3, 0.3)
plot.dp(x, x0)

# multiroot: sensible to x0!
a1 = multiroot(polygonOrigin, x0, d=x)
a2 = optim(x0, polygonOptim, d=x)

print(a1); cat("=========\n"); print(a2);

plot.dp(x, a2$par)
lines.dp(x, a1$root, col="red", lty=4, lwd=2)


### Ex 3:
x = runif(7, 1, 4)
x0 = c(1, 0.7)
plot.dp(x, x0)

# multiroot: sensible to x0!
a1 = multiroot(polygonOrigin, x0, d=x)
a2 = optim(x0, polygonOptim, d=x)

print(a1); cat("=========\n"); print(a2);

plot.dp(x, a2$par)
lines.dp(x, a1$root, col="red", lty=4, lwd=2)


### Ex 4:
x = runif(4, 1, 5)
x0 = c(1.3, 1.3)
plot.dp(x, x0)

# multiroot: sensible to x0!
a1 = multiroot(polygonOrigin, x0, d=x)
a2 = optim(x0, polygonOptim, d=x)

print(a1); cat("=========\n"); print(a2);

plot.dp(x, a2$par)
lines.dp(x, a1$root, col="red", lty=4, lwd=2)


### Ex 5:
x = runif(4, 1, 5)
# - 3 values of x0 are necessary in some situations;
# - but only optim handles properly 3 values;
# - multiroot: very sensible to x0 with 3 values!
x0 = c(0.327, 2.78, 1.72)
plot.dp(x, x0)

a1 = multiroot(polygonOrigin, x0, d=x)
a2 = optim(x0, polygonOptim, d=x)

print(a1); cat("=========\n"); print(a2);

plot.dp(x, a2$par)
lines.dp(x, a1$root, col="red", lty=4, lwd=2)

