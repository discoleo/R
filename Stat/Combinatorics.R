########################
###
### Leonard Mada
### [the one and only]
###
### Combinatorics
###
### draft v.0.1a


####################
### Introduction ###

### The famous Evening Combinatorics Provlems
# Best attacked before sleep, ... at least 3-4 hours before sleep.

# see Combinatorics.wiki


####################
####################

### helper functions

plot.cross = function(x1, x2, jitter=TRUE, col=c("red", "green"), ...) {
	ncols = length(col)
	plot(NA, NA, xlim=range(x1, x2), ylim = c(0, 1))
	id.j = if(jitter) jitter(seq_along(x1)) else 0;
	r = sapply(seq_along(x1), function(id)
		lines(c(if(jitter) id.j[id] else id, match(x1[id], x2)), c(1, 0), col=col[1 + (id %% ncols)]))
}

crossings = function(x1, x2, all=TRUE) {
	# assumes x1 = 1:n;
	less.f = function(id) {
		if(id == 1) return(0);
		sum(x2[1:(id - 1)] < x2[id])
	}
	interx = x2 - 1 - sapply(seq_along(x2), less.f)
	if(all) {
		return(rbind(x1, x2, interx))
	} else {
		return(interx);
	}
}

mean.crossings = function(n=10, iter=1000) {
	x1 = 1:n
	sum = 0
	for(i in seq(iter)) {
		x2 = sample(x1, n)
		sum = sum + sum(crossings(x1, x2, all=FALSE))
		if(i %% 100 == 1) {
			cat(i); cat(", ")
		}
	}
	cat("\n")
	return(sum / iter)
}

###############

n = 11

x1 = 1:n
x2 = sample(x1, n)

cross = crossings(x1, x2)
cross
sum(cross[3,])

plot.cross(x1, x2, jitter=F)

#########################


mean.crossings(n=11)

### TODO:
# - MC vs Quasi-Monte Carlo;


#########################

n = seq(7, 36)
cross.mean = sapply(n, mean.crossings, iter=6000)

plot(n, cross.mean)

summary(n^2 / cross.mean)
plot(n, n^2 / cross.mean)


cross.mean
n^2 / (4 + 4/n)

####################

####################
### Find Formula ###

library(gramEvol)

# - attempt to find a valid formula;

###########

### Grammer
ruleDef <- list(
	expr  = grule(op(expr, expr), func(expr), var),
	func  = grule(sqrt, log, exp, cosh),
	op    = grule('+', '-', '*', "/"),
	var   = grule(x, x^pow, div),
	pow   = grule(2, 1/3, 3),
	div   = grule(4, 3, 2, 1))

grammarDef <- CreateGrammar(ruleDef)

grammarDef


### Data
# from above
x = n;

### Fitness
fitness <- function(expr) {
	rmse(eval(expr), cross.mean)
}
rmse <- function(error, val) {
	if (any(is.nan(error))) return(Inf)
	r = sqrt(mean((error - val)^2, na.rm = TRUE))
}

### Evolution
ge = GrammaticalEvolution(grammarDef, fitness, terminationCost = 0.005, iterations=3000)

ge

# most likely: x^2 / 4 or similar formulas;


