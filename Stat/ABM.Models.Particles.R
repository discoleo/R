########################
###
### ABM
### Turtle Dynamics
###
### UVT Team Project 2021
### Simulations in Epidemiology
###
### supervised by:
### Leonard Mada
### [the one and only]
###
### draft v.0.1e

# - based on:
#   https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.04516
# - the original used package: lcmix;
# - but this package is NOT maintained anymore;
# install.packages("lcmix", repos="http://R-Forge.R-project.org")
# library(lcmix) # multivariate distributions

################

###############
### History ###
###############

### draft v.0.1e:
# - exploring custom defined methods;
### draft v.0.1b - v.0.1d:
# - added lines.turtles();
# - improved colors;
# - path inherits colors;
### draft v.0.1a:
# - replaced package lcmix with LaplacesDemon;
# - TODO: multivariate rgamma;


################


# install.packages("NetLogoR")
# install.packages("LaplacesDemon")

library(NetLogoR)
library(LaplacesDemon) # multivariate distributions

######################
######################

### helper Functions

rWorld.gen = function(size, ntypes) {
	createWorld(minPxcor=1, maxPxcor=size, minPycor=1, maxPycor=size,
		sample(seq(ntypes), size*size, replace = TRUE));
}
run.model = function(turtles, land, distRate, iter=10, plot=FALSE, pch=16) {
	# t.df = data.frame(xcor=numeric(), ycor=numeric(), who=numeric());
	t.df = data.frame(turtles@.Data[ , c("xcor", "ycor", "who")]);
	for(i in seq(iter)) {
		# Identify the cells the turtles are on
		cellTurtle = patchHere(world=land, turtles=turtles)
		### Distance
		distMove = of(world=land, agents=cellTurtle)
		distShape = distMove * distRate
		ntrt = nrow(turtles);
		rho = matrix(rep(0.8, length=ntrt * ntrt), ncol=ntrt)
		diag(rho) = 1 # used for rvgamma & TODO;
		# distMoveRan = rvgamma(2, distShape, distRate, rho)[1, ]
		distMoveRan = rgamma(ntrt, distShape, distRate)
		# set Distance
		turtles = fd(turtles=turtles, dist=distMoveRan, world=land, torus=FALSE, out=FALSE)
		### Direction
		meanHeading = mean(of(agents=turtles, var="heading"))
		Sigma = matrix(rep(0.8 * meanHeading, length=ntrt * ntrt), ncol=ntrt)
		diag(Sigma) = meanHeading
		angleInd = rmvn(n=1, mu=rep(meanHeading, ntrt), Sigma=Sigma)
		# set Direction
		turtles = right(turtles=turtles, angle=as.vector(angleInd))
		if(plot) points(turtles, pch=pch, col= of(agents=turtles, var="color"))
		t.df = rbind(t.df, turtles@.Data[,c("xcor", "ycor", "who")]);
	}
	attr(t.df, "col") = of(agents=turtles, var="color");
	rez = list(t=turtles, path=t.df);
	class(rez) = c("agentsWithPath", class(rez));
	invisible(rez);
}
plot.turtles = function(turtles, pch=16, col, ...) {
	if(missing(col)) col = of(agents=turtles, var="color");
	points(turtles@.Data, pch=pch, col=col, ...);
}
# overwriting the default function:
plot.agentMatrix = function(x, pch=16, col, ...) {
	plot.turtles(x, pch=pch, col=col, ...);
}
setClass("agentsWithPath", contains = c(t="agentMatrix"))
setMethod("plot", signature(x="agentsWithPath", y="missing"),
	function(x, pch=16, col, ...) {
		plot.turtles(x$t, pch=pch, col=col, ...);
	} )
# plot the Path
lines.turtles = function(path, col) {
	if(inherits(path, "agentsWithPath")) path = t.all$path;
	if(missing(col)) col = attr(path, "col");
	who = unique(path$who);
	sapply(seq(length(who)),
		function(id) lines(path[path$who == (id-1), 1:2], col=col[id]))
	invisible();
}

######################

### Options
size = 20
agents = 4
# dynamics
distRate = 0.5

### Grid
# create raster-grid
land = rWorld.gen(size, ntypes=2)

### Agents
t1 = createTurtles(n=agents, world=land)

# Visualize the turtles on the landscape with their respective color
plot(land, c(1,2), col=c("#FFFFFF", "#00FF0064"))
plot(t1)
# points(t1, pch=16, col=of(agents=t1, var="color"))


### MODEL

# plot(land)

# Note:
# - the code below can be repeated a few times,
#   allowing the turtles to move a larger distance;

t.all = run.model(t1, land, distRate);
t1 = t.all$t;
plot(t.all)
lines.turtles(t.all)

