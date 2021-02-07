########################
###
### Leonard Mada
### [the one and only]
###
### Combinatorics
###
### draft v.0.1c


####################
### Introduction ###

### The famous Evening Combinatorics Problems
# Best attacked before sleep, ... at least 3-4 hours before sleep.

# see Combinatorics.wiki


####################
####################

### helper functions

plot.cross = function(x1, x2, jitter=TRUE, col=c("red", "green"), add=FALSE, ...) {
	ncols = length(col)
	ylim = c(0, 1)
	if( ! add) {
		plot(NA, NA, xlim=range(x1, x2), ylim = ylim)
	}
	id.j = if(jitter) jitter(seq_along(x1)) else 0;
	r = sapply(seq_along(x1), function(id)
		lines(c(if(jitter) id.j[id] else id, match(x1[id], x2)), c(ylim[2], ylim[1]), col=col[1 + (id %% ncols)]))
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


plot(n, n^2 / cross.mean - (3.85 + 1/log(n)))


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


#####################
#####################

### Famous Square ###

n = 11

x1 = 1:n
x2 = sample(x1, n)
x3 = sample(x1, n)
x4 = sample(x1, n)

plot.crossy = function(x1, x2, jitter=TRUE, col=c("purple", "blue")) {
	ncols = length(col)
	len = length(x1) + 0.1
	id.j = if(jitter) jitter(seq_along(x1)) else 0;
	r = sapply(seq_along(x1), function(id)
		lines(c(1, length(x1)), c(if(jitter) id.j[id] else id, match(x1[id], x2))/len, col=col[1 + (id %% ncols)]))
}

### Opposite sides:
# - the number of inter-crossings is fixed: + n^2;
plot.cross(x1, x2)
plot.crossy(x3, x4)


######################
######################

### Cyclic sequence

plot.cycle = function(s1, s2, r=c(1, 2), col=c("red", "green"), f=1.25, add=FALSE) {
	len1 = length(s1); len2 = length(s2);
	plot.cyc.text = function(s, id=1) {
		len = length(s); arc = 2*pi/len;
		x = r[id] * cos( (0:(len-1)) * arc - pi/2);
		y = r[id] * sin( (0:(len-1)) * arc - pi/2);
		text(x, y, s, col=col[id])
	}
	# plot(x1, y1)
	r.max = max(r) * f; lim = c(-r.max, r.max);
	old.par = par(mar=c(1,1,1,1))
	plot.new(); plot.window(xlim=lim, ylim=lim);
	plot.cyc.text(s1, id=1)
	plot.cyc.text(s2, id=2)
	par(old.par)
}
shift.seq = function(s, first=1) {
	id = match(first, s)
	if(is.na(id)) stop("NO such element!")
	if(id > 1) {
		s = c(tail(s, 1-id), head(s, id-1))
	}
	return(s)
}


### Test
n = 6
s1 = seq(n)
s2 = sample(s1, n)
s2 = shift.seq(s2)

### Plot:
plot.cycle(s1, s2)


### Save image
SAVE=FALSE
if(SAVE) {
	id = 2
	png(file=paste0("img/Combinatorics.Cyclic.", id, ".png"))
	plot.cycle(s1, s2)
	dev.off()
}


