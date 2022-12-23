

library(biogeom)


####################
### Enhancements ###
####################


### 1. Return Invisibly

# inside code
return(invisible(list(x = x.coordi, y = y.coordi)));
# last line
invisible(list(x = x.coordi, y = y.coordi));


### 2. Plot Offset

# Add one of the following arguments to the curve/plot functions:
# offset = c(0,0);
# origin = c(0,0);
# xy0 = c(0,0);
# center = c(5,5);

### Note:
# - the center depends on the graphic;
# - the preferred argument is one of offset or origin;


### 3. Plot Rings

plot.rings = function(x, which=0, id0 = 10,
		ub.np = 200, len.pro=1/10, ub0.np=2000, len0.pro=1/20,
		col=1, col0="grey73", lwd=1, lwd0=3, xlim=c(3, 13), ylim=xlim, add=FALSE, ...) {
	xy = x; # initial data;
	# All rings
	uni.C <- sort( unique(xy$Code) );
	if(length(which) > 1 || which[1] != 0) {
		# TODO: possible to save subset in new variable;
		uni.C = intersect(uni.C, which);
	}
	if( ! add) {
		plot.new();
		plot.window(xlim=xlim, ylim=ylim, asp=1, cex.lab=1.5, cex.axis=1.5,
			xlab=expression(italic("x")), ylab=expression(italic("y")));
	}
	plotf = function(code, col, lwd) {
		Data <- xy[xy$Code == code, ];
		x0 <- Data$x;
		y0 <- Data$y;
		Res <- adjdata(x0, y0, ub.np=ub.np, len.pro=len.pro);
		# close the circle
		Res$x = c(Res$x, Res$x[1]);
		Res$y = c(Res$y, Res$y[1]);
		lines(Res$x, Res$y, col=col, lwd=lwd, type="l");
	}
	for(i in 1:length(uni.C)) {
		plotf(uni.C[i], col=col, lwd=lwd);
	}
	# Base ring:
	if(id0 > 0 && id0 %in% uni.C) {
		plotf(id0, col=col0, lwd=lwd0);
	}
}

# Plot all rings:
plot.rings(whitespruce)

# Plot a specific subset of rings:
plot.rings(whitespruce, which=c(10, seq(15,1000, by=2)))

# Advanced:
class(whitespruce) = c("rings", "data.frame");
# the plot command works now automatically:
# (Note: plot.rings is defined above)
plot(whitespruce, which=seq(2,1000, by=2));

