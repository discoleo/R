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
### draft v.0.2b

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

### draft v.0.2a - v.0.2b:
# - [refactoring] S4 classes;
# - additional parameters passed to lines();
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

### S4 Class: agentsWithPath

.agentsWithPath = setClass("agentsWithPath",
	slots = c(Agents="agentMatrix", Path="data.frame"))

setMethod(
	"show",
	signature(object = "agentsWithPath"),
	definition = function(object) {
		tmp = object@Agents;
		print(class(tmp))
		show(tmp)
})

setGeneric("agentsWithPath",
	function(agents, path, ...) {
		standardGeneric("agentsWithPath")
	}
)
setMethod("agentsWithPath",
	signature = c(agents = "agentMatrix"),
	definition = function(agents, path, ...) {
		obj = new("agentsWithPath", Agents = agents, ...);
		if( ! is.null(path)) obj@Path = path;
		return(obj);
	}
)

### Methods

### Plot:
plot.agentsWithPath = function(x, pch=16, col, ...) {
	x = x@Agents;
	if(missing(col)) col = of(agents=x, var="color");
	points(x@.Data, pch=pch, col=col, ...);
}
# overwriting the default function:
plot.agentMatrix = function(x, pch=16, col, ...) {
	if(missing(col)) col = of(agents=x, var="color");
	points(x@.Data, pch=pch, col=col, ...);
}

### Lines:
lines.agentsWithPath = function(x, col, ...) {
	if(missing(col)) col = of(agents=x@Agents, var="color");
	path = x@Path;
	who  = unique(path$who); len = length(who);
	if(length(col) < len) col = rep(col, len);
	sapply(seq(len),
		function(id) lines(path[path$who == (id-1), 1:2], col=col[id], ...))
	invisible();
}

#####################

### Helper Functions:

rWorld.gen = function(size, ntypes) {
	createWorld(minPxcor=1, maxPxcor=size, minPycor=1, maxPycor=size,
		sample(seq(ntypes), size*size, replace = TRUE));
}
run.model = function(turtles, land, distRate, iter=10, add=FALSE, plot=FALSE, pch=16) {
	# t.df = data.frame(xcor=numeric(), ycor=numeric(), who=numeric());
	if(inherits(turtles, "agentsWithPath")) {
		t.df = if(add) turtles@Path
			else {
				data.frame(turtles@Agents@.Data[ , c("xcor", "ycor", "who")]);
			}
		turtles = turtles@Agents;
	} else if(inherits(turtles, "agentMatrix")) {
		t.df = data.frame(turtles@.Data[ , c("xcor", "ycor", "who")]);
	} else stop("Agents NOT supported!")
	#
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
		# Plot:
		if(plot) points(turtles, pch=pch, col= of(agents=turtles, var="color"))
		t.df = rbind(t.df, turtles@.Data[,c("xcor", "ycor", "who")]);
	}
	rez = agentsWithPath(turtles, t.df);
	invisible(rez);
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
t.all = createTurtles(n=agents, world=land)

# Visualize the turtles on the landscape with their respective color
plot(land, c(1,2), col=c("#FFFFFF", "#00FF0064"))
plot(t.all)


### MODEL

# Note:
# - the code below can be repeated a few times,
#   allowing the turtles to move a larger distance;

t.all = run.model(t.all, land, distRate);
plot(t.all)
lines(t.all, lty=3)


##############
