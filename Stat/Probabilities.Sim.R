#######################
###
### Leonard Mada
### [the one and only]
###
### Coupled / Recursive Poisson Processes
###
### draft v.0.1c


### Coupled / Recursive Poisson Processes
# - some basic experiments in Probabilities;


######################

rpois.heavy = function(n, lambda, rm.zero=FALSE) {
	# heavy-tailed pseudo-Poisson distribution
	x = rpois(n, lambda)
	if(rm.zero) x = x[x != 0];
	y = sapply(x, function(n) rpois(1, n))
	invisible(y);
}
rpois.cor = function(n, lambda, dL, scale.by=0, rm.zero=TRUE,
		type=c("recursive", "binomial")) {
	# correlated pseudo-Poisson processes
	# rm.zero = remove the 0s in the initial process;
	x = rpois(n, lambda1)
	if(rm.zero) x = x[x != 0];
	# Type
	type = match.arg(type);
	if(is.na(type)) stop("Type NOT supported!")
	if(type == "recursive") {
		x1 = x - dL; x1[x1 < 0] = 0;
		if(scale.by != 0) {
			isGr = (x > lambda); isLs = (x < lambda);
			x1[isGr] = x1[isGr] - scale.by;
			x1[isLs] = x1[isLs] + scale.by;
			x1[x1 < 0] = 0;
		}
		y = sapply(x1, function(n) rpois(1, n));
	} else if(type == "binomial") {
		p = 0.5 - dL/lambda / 2;
		y = sapply(x, function(n) rbinom(1, 2*n, p));
	}
	invisible(list(x=x, y=y));
}
simulate = function(n, FUN, FUN2, iter=100, ...) {
	# FUN = generator function;
	m = numeric(iter);
	for(id in seq(iter)) {
		r = FUN(n, ...);
		m[id] = FUN2(r);
	}
	return(m)
}
sum.rpois.cor = function(l) {
	sum(l$x, -l$y);
}

###########

### Example
dL = 1.5
lambda1 = 3

n = 50

### Simulation
# Central Limit Theorem: holds;
iter = 1000
m = simulate(n, rpois.cor, sum.rpois.cor, iter=iter, lambda=lambda1, dL=dL)

mean(m); sd(m);
hist(m, breaks=20)

plot(density(m))
gs.seq = seq(min(m), max(m))
y = dnorm(gs.seq, mean(m), sd(m));
lines(gs.seq, y, col="red")

############

############
### Test ###

scale.by = 0
r = rpois.cor(n, lambda1, dL=dL, scale.by=scale.by)
x = r$x; y = r$y;

### Analysis
length(x); mean(x); table(x); table(y)
table(y <= x)
table(y - x)

v = c(x, -y)
hist(v, breaks=12)
sum(v) / n; cor(x, y);
plot(jitter(x), jitter(y))

### Heavy-tailed (negative) Poisson
x = rpois(1000, lambda1*2 - dL) -  lambda1 - dL
plot(density(v))
lines(density(x), col="red")

###
n = 1000
lambda = 5
x = rpois.heavy(n, lambda)
mean(x)

plot(density(x))
gs.seq = seq(min(x), max(x))
y1 = dpois(gs.seq, lambda);
y2 = dpois(gs.seq, lambda + 0.5);
lines(gs.seq, y1, col="red")
lines(gs.seq, y2, col="green")

### Heatmap
cor(x, y)
tbl = table(x, y)
heatmap(tbl, Rowv=NA, Colv=NA)


############

### Test 2:
### Correlated events
n = 1000
lambda1 = 3
dL = 0.25
#
r = rpois.cor(n, lambda1, dL=dL, type="binomial")
x = r$x; y = r$y;

### Analysis
length(x); mean(x); mean(y);
table(x); table(y)
table(y <= x)
table(y - x)

v = c(x, -y)
hist(v, breaks=12)
sum(v) / n;

### Heatmap
cor(x, y)
tbl = table(x, y)
heatmap(tbl, Rowv=NA, Colv=NA)


######################
######################

### TODO

lambda2 = lambda1 - dL

nx = rpois(n, lambda1)
ny = rpois(n, lambda2)

x = sapply(nx, function(n) rpois(n, n))
y = sapply(ny, function(n) rpois(n, n))

table(sapply(x, length))
table(sapply(y, length))

v = c(unlist(x), -unlist(y))
hist(v)
sum(v)


isMiss = sapply(x, function(x) length(x) == 0)
table(unlist(y[isMiss]))
tmp = x[isMiss]
x[isMiss] = y[isMiss]
y[isMiss] = tmp;
v = c(unlist(x), -unlist(y))
hist(v)
sum(v)

