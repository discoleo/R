
####################
### Random Walks ###
####################

## Leonard Mada
##
## draft v.0.1d

# 1. Initial code for this module started as an issue:
#    https://github.com/jmzobitz/ModelingWithR/issues/22
# 2. Related Topics: Percolation
#    https://github.com/discoleo/Percolation


####################

### 2D

library(rgl)

### 8 Neighbours
# - can move on Diagonals as well;
rwalk.2D8N = function(n, x, y, xprob = NULL, yprob = NULL, iter = 100) {
	n = n + 1; # for Time = 0;
	r = lapply(seq(iter), function(id) {
		x = sample(x, n, prob = xprob, replace = TRUE);
		y = sample(y, n, prob = yprob, replace = TRUE);
		x[1] = 0; y[1] = 0; # Start Position
		x = cumsum(x); y = cumsum(y);
		tt = seq(0, by=1, length.out = n)
		data.frame(t = tt, x=x, y=y, id=id);
	});
	r = do.call(rbind, r);
	attr(r, "n") = n;
	class(r) = c("rw", class(r));
	invisible(r);
}

### 4 Neighbours
# - will NOT move on Diagonals;
rwalk.2D4N = function(n, x, y, xprob = NULL, yprob = NULL, iter = 100) {
	xw = x; yw = y;
	n = n + 1; # for Time = 0;
	r = lapply(seq(iter), function(id) {
		isX = sample(c(TRUE, FALSE), n, replace = TRUE);
		nX = sum(isX);
		x  = rep(0, n);
		x[isX] = sample(xw, nX, prob = xprob, replace = TRUE);
		nY = n - nX;
		y  = rep(0, n);
		y[! isX] = sample(yw, nY, prob = yprob, replace = TRUE);
		x[1] = 0; y[1] = 0; # Start Position
		x = cumsum(x); y = cumsum(y);
		tt = seq(0, by=1, length.out = n)
		data.frame(t = tt, x=x, y=y, id=id);
	});
	r = do.call(rbind, r);
	attr(r, "n") = n;
	class(r) = c("rw", class(r));
	invisible(r);
}

### Analysis
atan.rw = function(x) {
	th = atan2(x$y, x$x);
	th = data.frame(th=th, id = x$id, t = x$t);
	attr(th, "n") = n;
	class(th) = c("rwth", class(th));
	invisible(th);
}
quantile.rwth = function(x, quantiles = c(0.25, 0.5, 0.75), na.rm = TRUE) {
	r = tapply(x$th, x$t, function(x) {
		as.data.frame(as.list(
			quantile(x, probs = quantiles, na.rm=na.rm)));
	});
	r = do.call(rbind, r);
	names(r) = paste0("q", quantiles);
	r$t = sort(unique(x$t));
	attr(r, "quantiles") = quantiles;
	invisible(r);
}

### Plot:
# sc = scale Base-Square by sc * sqrt(n);
plot.rw = function(x, col = 1:4, alpha = 0.75, b.alpha = 0.6,
		type = c("4", "2", "A4", "A2", "Custom"),
		base = TRUE, ..., sc = 2) {
	n  = attr(x, "n");
	n2 = n %/% 2;
	type  = match.arg(type);
	doCol = type != "Custom";
	if(doCol) {
		if(length(col) == 1) { col = rep(col, 4); }
		else if(length(col) == 2) {col = rep(col, 2); }
	}
	### Plot:
	# Base:
	if(base) {
		sc = sc * sqrt(n);
		polygon3d(rep(0,4), c(1,1,-1,-1)*sc, c(1,-1,-1,1)*sc,
			alpha = b.alpha, col = "red");
	}
	# Walk:
	x   = split(x, x$id);
	tmp = lapply(x, function(x) {
		if(doCol) {
			if(type == "4") {
				idCol = (x$x[n2] < 0) * 2 + (x$y[n2] < 0) + 1;
			} else if(type == "2") {
				idCol = if(x$x[n2] < 0) "black" else "red";
			} else if(type == "A4") {
				# Alternating 4 cols:
				idCol = (x$id %% 4) + 1;
			} else {
				# Alternating 2 cols:
				idCol = (x$id %% 2) + 1;
			}
			col = col[idCol];
		}
		lines3d(x$t, x$x, x$y, col=col, alpha=alpha, ...);
	});
	invisible();
}

# TODO: plot using more colours;

################

### Examples
n = 121;
x = c(-2,-1,0,1,2); xprob = c(1,2,1,2,1)/7;
y = c(-1,1); yprob = c(1,1)/2;
walk = rwalk.2D8N(n, x, y, xprob, yprob)
plot(walk, alpha = 0.6)


###
close3d()
n = 121;
x = c(-2,-1,0,1,2); xprob = c(1,2,1,2,1)/7;
y = c(-1,1); yprob = c(1,1)/2;
walk = rwalk.2D4N(n, x, y, xprob, yprob)
plot(walk, alpha = 0.6)


### Difference: 4N vs 8N
# Walk in 4N is delayed compared to 8N;
close3d()
n = 120;
x = c(-2,-1,0,1,2); xprob = c(1,2,1,2,1)/7;
y = c(-1,1); yprob = c(1,1)/2;
walk4 = rwalk.2D4N(n, x, y, xprob, yprob)
walk8 = rwalk.2D8N(n, x, y, xprob, yprob)
plot(walk4, col = c("#12E844", "#B2FF44"), alpha = 0.6)
plot(walk8, col = c("red", "#FF4496"), alpha = 0.6, base = FALSE)


### Phase:
n = 120;
x = c(-2,-1,0,1); xprob = c(2,1,1,5)/9;
# x = c(-3,-2,-1,0,1); xprob = c(1,1,2,1,7)/12;
# x = c(-3,-2,-1,0,1,2,3); xprob = c(1,1,2,1,2,1,1)/9;
y = c(-1,2); yprob = c(2,1)/3;
walk = rwalk.2D8N(n, x, y, xprob, yprob, iter = 400);
th = quantile.rwth(atan.rw(walk));
matplot(th$t, th[names(th) != "t"], type = "l", lty = 1, ylab = "Phase")


### Note: EEG
# - Phase resembles the fast waves (e.g. beta waves)
#   in an Electroencephalogram (EEG);
# - Other wave types may be possible to generate
#   using different values for the random walk;
# - EEG may be actually the ensemble average
#   of hundreds of basic neuron-units;
# - Random walk may represent either the activation
#   of neurons at a distance proportional to the walk,
#   or the number of neurons activated by one input neuron;
# - Both values may be actually correlated: if the signal
#   travels further outwards, more neurons get activated;

