

### Lecture: Over-Fitting

# We want to predict a binary class:
# - Data: 2 predictor variables;
# - Outcome: binary class;
# - Problem: Can we obtain a perfect classifier?


##############

library(shape)


### Data

# Helper Functions:
plot.data0 = function(x, y, slope = 3/4, r=1/6, col = c("#D03264", "#6432D0")) {
	plot.new();
	plot.window(xlim = c(0.5,10.5), ylim=c(0.5,10.5));
	for(yy in y) {
		isUp = (x*slope) >= yy;
		for(xx in x[isUp]) {
			plotcircle(mid=c(xx, yy), r=r, col=col[1]);
		}
		for(xx in x[ ! isUp]) {
			plotcircle(mid=c(xx, yy), r=r, col=col[2]);
		}
	}
}
points.data = function(x, y, col, r=1/6) {
	for(id in seq_along(x)) {
		plotcircle(mid=c(x[id], y[id]), r=r, col=col[1]);
	}
}
plot.data = function(x, y, out = NULL, main, col, FUN=NULL) {
	# out = outliers;
	par.old = par(mar = c(0.5,0,1,0))
	plot.data0(px, py, col=col)
	if( ! is.null(out)) {
		points.data(out[[1]]$x, out[[1]]$y, col=col[2]);
		points.data(out[[2]]$x, out[[2]]$y, col=col[1]);
	}
	text(median(x), max(y) + 1.75, main, cex = 2, col="#16D048")
	if( ! is.null(FUN)) FUN();
	par(par.old);
}


# Data Points
px = seq(10)
py = seq(8)
out = list(
	list(x = c(7.5, 7.75), y = c(4.25, 3.5)),
	list(x = c(3.5, 3), y = c(3.5, 4)))
col = c("#D03264", "#6432D0")
plot.data(px, py, out=out, col=col, main = "The Data")


########################

### Which Model is Best?

par.old = par(mfrow = c(2,1))

FUN = function() lines(c(0.5, 10.25), c(0.4, 7.75), col="green", lwd=2);
plot.data(px, py, out=out, col=col, main = "Model 1", FUN = FUN)

FUN = function() {
	lines(
		c(0.5, 3.45, 2.5, 3, 4, 7.1, 7.35, 7.7, 8, 7.5, 10.25),
		c(0.4, 2.8,  4,  4.75,  3.25, 5.6, 3.75, 3.1, 3.5, 5.8, 7.75),
		col="green", lwd=2);
}
plot.data(px, py, out=out, col=col, main = "Model 2", FUN = FUN)

par(par.old)

##########################

### Over-Fitting

# Model 2: is *NOT* better!
# - it overfits the data;
# - the 2 pockets with different-type data are *anomalies";
# - they arise because of the limited sample size and
#   because of the projection of the *complete/complex* data on a 2D space;

# Real Data:
# - the entire genome at birth +
#   all the environmental data since birth!
# - it is a trillion-dimensional data set;
# - if we project this data set on a 2D space,
#   then the 2 classes will inevitably overlap;
# - this overlap would be obvious if we plot the entire population;

### Entire population:
# Entire population = 8 billion humans + past population + future population;

px = seq(1, 10, length.out = 30);
py = seq(1, 8, length.out = 24);
out = list(
	list(
		x = c(seq(2, 9.5, by=0.5) + 0.25, seq(3, 9) + 0.5),
		y = c(seq(2, 9.5, by=0.5) - 0.5, seq(3, 9) - 1) * 3/4),
	list(
		x = c(seq(1, 9, by=0.5) + 0.25, seq(1, 9) + 0.5),
		y = c(seq(1, 9, by=0.5) + 2/3, seq(1,  9) + 5/3) * 3/4) )
FUN = function() lines(c(0.5, 10.25), c(0.4, 7.75), col="green", lwd=2);
plot.data(px, py, out=out, col=col, main = "Entire Population", FUN = FUN)

### Problem:
# - there is *NO* way to correctly fit/classify the data
#   using the 2 input variables;
#   (or any limited-dimensional data sets that happen to be available)
# - the real population will always overlap inside this limited data-space:
#   the initial pockets were spurious anomalies;
# - Fundamental limitation: we do not have access to the trillion-dimensional data-set,
#   and therefore will not be able to correctly classify
#   *ALL* individuals of the population;
# - any such attempts will over-fit the data: applying the over-fitted model
#   to the entire population will induce massive miss-classifications;
# => the real reason to avoid over-fitting!


#########################
#########################

### Future Topics:
# - Gradient descent vs Metaheuristic algorithms
#   (including Monte Carlo / Gibbs sampling);
# - Spin glasses: long term "memory" & local minima;
# - Long range correlations;

