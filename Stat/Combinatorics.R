########################
###
### Leonard Mada
### [the one and only]
###
### Combinatorics
###
### draft v.0.1h


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
# approx n^2 / 4
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

plot.cycle = function(s1, s2, r=c(1, 2), col=c("red", "green"), f=1.25, plot=TRUE, add=FALSE) {
	len1 = length(s1); len2 = length(s2);
	plot.cyc.text = function(s, id=1) {
		len = length(s); arc = 2*pi/len;
		x = r[id] * cos( (0:(len-1)) * arc - pi/2);
		y = r[id] * sin( (0:(len-1)) * arc - pi/2);
		text(x, y, s, col=col[id])
		return(cbind(x,y))
	}
	# plot(x1, y1)
	r.max = max(r) * f; lim = c(-r.max, r.max);
	old.par = par(mar=c(1,1,1,1))
	plot.new(); plot.window(xlim=lim, ylim=lim);
	xy1 = plot.cyc.text(s1, id=1)
	xy2 = plot.cyc.text(s2, id=2)
	par(old.par)
	return(list(xy1, xy2))
}
shift.by = function(s, by=1) {
	if(by >= length(s)) by = by %% length(s);
	if(by == 0) return(s);
	c(tail(s, -by), head(s, by))
}
shift.seq = function(s, by, first) {
	if(missing(by)) {
		id = match(first, s)
		if(is.na(id)) stop("NO such element!")
		if(id > 1) {
			s = c(tail(s, 1-id), head(s, id-1))
		}
	} else {
		s = shift.by(s, by)
	}
	return(s)
}
shift.max = function(s, ref) {
	by = match.max(s, ref)
	return(shift.by(s, by))
}
match.max = function(s1, s2) {
	len = length(s1) - 1
	id = 0:len;
	matches = sapply(id, function(id) sum(shift.by(s1, id) == s2))
	max.id = match(max(matches), matches);
	return(max.id - 1);
}
connect = function(s1, s2, xy, col=c("orange", "purple"), lwd=2, delta=0.1, add=TRUE) {
	if(missing(xy)) xy = plot.cycle(s1, s2);
	delta.f = function(x) {
		dl = round(x[1] - x[2], 4)
		if(dl > 0) {
			x = x + delta * c(-1, 1);
		} else if(dl < 0) {
			x = x + delta * c(1, -1);
		}
		return(x)
	}
	lines.match = function(id, diff=0, id2=id - diff, col) {
		x = c(xy[[1]][id,1], xy[[2]][id2,1])
		y = c(xy[[1]][id,2], xy[[2]][id2,2])
		# shrink
		x = delta.f(x); y = delta.f(y);
		lines(x, y, col=col, lwd=lwd)
	}
	#
	isEq = s1 == s2;
	if(any(isEq)) {
		id = seq(length(s1))[isEq]
		sapply(id, lines.match, col=col[1])
	}
	# shifted
	con.shift = function(sh=1, ...) {
		z.sh = rep(0, sh);
		# shifted R;
		isEq = c(z.sh, s1) == c(s2, z.sh);
		if(any(isEq)) {
			isEq = head(isEq, -sh)
			id = seq(length(s1))[isEq]
			sapply(id-sh, lines.match, diff=-sh, ...)
		}
		# shifted L: TODO: check!
		isEq = c(s1, z.sh) == c(z.sh, s2);
		if(any(isEq)) {
			isEq = head(isEq, -sh)
			id = seq(length(s1))[isEq]
			sapply(id, lines.match, diff=sh, ...)
		}
	}
	con.shift(1, col=col[1]); con.shift(2, col=col[2]);
}

crossings.cyc = function(s) {
	len = length(s)
	len2 = length(s) / 2
	s0 = seq(len)
	# from L to R
	typeL1 = (s[s0] >= s0) & (s[s0] <= s0 + len2);
	typeL2 = (s[s0] < s0) & (s[s0] + len2 < s0);
	typeL = typeL1 | typeL2;
	# TODO: check R-R crossings;
	#
	cr = sapply(s0, function(id) {
			# avoid double counting!
			if( ! typeL[id]) {
				# R-R crossings;
				sp.id = c(); sn.id = c();
				if(id == len) {
					sn.id = 1:(s[id] - len2);
				} else if(s[id] < id) {
					sp.id = (id+1):min(len, s[id] + len2, id+len2);
				} else {
					# R-R crossings across START pos;
					sn.id = (id+1):(s[id] - len2);
				}
				cr = sum((s[sp.id] < s[id]) & (s[sp.id] + len2 > sp.id))
				# R-R crossings across START pos;
				cr = cr + sum((s[sn.id] < s[id]) & (s[sn.id] > sn.id + len2))
				return(cr);
			}
			# upper bound
			if(typeL1[id]) {
				sn.id = c()
				if(id == len) {
					sp.id = 1:(len2 - 1)
				} else if(s[id] <= len2) {
					sp.id = (id+1):(s[id]+len2)
				} else if(s[id] > len2) {
					sp.id = (id+1):len
					sn.id = 1:(s[id] - len2)
				}
				cr = sum((s[sp.id] < s[id]) & (s[sp.id] + len2 > sp.id))
				cr = cr + sum(s[sp.id] > sp.id + len2)
				if(length(sn.id) > 0) {
					cr = cr + sum((s[sn.id] < s[id]) & (s[sn.id] > sn.id + len2))
				}
			} else {
				# type L2: crosses over [1]
				sn.id = c()
				if(id == len) {
					sp.id = 1:(s[id] + len2 - 1)
				} else {
					sn.id = (id+1):len
					sp.id = c(sn.id, 1:(s[id] + len2 - 1))
				}
				cr = sum((s[sp.id] < s[id]) & (s[sp.id] + len2 < sp.id))
				cr = cr + sum((s[sp.id] > s[id]) & (sp.id + len2 < s[sp.id]))
				if(length(sn.id) > 0) {
					cr = cr + sum(s[sn.id] + len2 > id)
				}
			}
			return(cr)
		})
	return(cr)
}


##############

### Test
n = 20
s1 = seq(n)
s2 = sample(s1, n)
s2 = shift.max(s2, s1)

### Plot:
xy = plot.cycle(s1, s2)
connect(s1, s2, xy)


cr = crossings.cyc(s2)
sum(cr); cr;


### Save image
SAVE=FALSE
if(SAVE) {
	id = 1
	png(file=paste0("img/Combinatorics.Cyclic.", n, ".", id, ".png"))
	xy = plot.cycle(s1, s2)
	connect(s1, s2, xy)
	dev.off()
}



### Upper bound to Average:
n = 60
it = 1000
s1 = 1:n
sum( sapply(1:it, function(id) sum(crossings.cyc(shift.max(sample(s1, n), s1)))) ) / it / n



### Average Number of Crossings

### Liniar case:
# - aprox n^2 / 4;

### Cyclic case:
# - there are 2 pathways for the matching;
#  -- the max number is therefore at most 1/2 of the liniar case;
#  -- especially the extreme outliers are therefore reduced;
# - in addition: we can rotate one of the sequences
#   to better match the 2 sequences;
#  -- the "rotation" can be virtual,
#     but a true rotation aids visualisation;
# - the average may be far less than n^2 / 8;
# - the upper bound seems to converge to aprox n^2 / 8;
### Q: 
# Can it converge to some O(n^(1.5))?
# Or does if converge also to some O(n^2)?

