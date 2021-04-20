
### Coupled / Recursive Poisson Processes

# some basic experiments in Probabilities;


######################

rpois.heavy = function(n, lambda, rm.zero=FALSE) {
	# heavy-tailed Poisson distribution
	x = rpois(n, lambda)
	if(rm.zero) x = x[x != 0];
	y = sapply(x, function(n) rpois(1, n))
	invisible(y);
}
rpois.cor = function(n, lambda, dL, rm.zero=TRUE) {
	# correlated Poisson processes
	x = rpois(n, lambda1)
	if(rm.zero) x = x[x != 0];
	y = sapply(x - dL, function(n) rpois(1, max(0, n)))
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

r = rpois.cor(n, lambda1, dL=dL)
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

