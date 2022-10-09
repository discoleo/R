########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Generic Solvers



#######################

# library(polynom)
# library(pracma)

### Helper Functions
source("Polynomials.Helper.R")


### fast load:
# source("Polynomials.Helper.Solvers.R")


#######################

### Solvers:

# TODO:
# - the code that uses these functions needs a lot of refactoring & cleanup!


### Simple systems:
# - creates a full set of permutations of x;
# - valid only for totally symmetric systems!
# - TODO: rename to permute.sol();
solve.En = function(x, max.perm=0, n=4, duplicates=FALSE) {
	warning("Deprecated! Use permute.sol instead!");
	return(permute.sol(x, max.perm=max.perm, n=n, duplicates=duplicates));
}
permute.sol = function(x, max.perm=0, n=4, duplicates=FALSE) {
	id = 1:length(x)
	if(max.perm == 1) {
		id.gr = matrix(id, nrow=1)
	} else {
		id.l = rep(list(id), n)
		id.gr = expand.grid(id.l)
		if( ! duplicates) {
			isDuplic = apply(id.gr, 1, function(id.val) any(duplicated(id.val)))
			id.gr = id.gr[ ! isDuplic , ]
		}
		if(max.perm > 0) {
			max.perm = min(max.perm, nrow(id.gr));
			id.gr = head(id.gr, n=max.perm);
		}
	}
	sol = if(n == 4) cbind(
			x1=x[id.gr[,1]], x2=x[id.gr[,2]], x3=x[id.gr[,3]], x4=x[id.gr[,4]])
		else cbind(x[id.gr[,1]], x[id.gr[,2]], x[id.gr[,3]])
	sol.names = if(n >= 4) paste0("x", seq(n)) else
		if(n == 3) c("x", "y", "z") else c("x", "y");
	colnames(sol) = sol.names;
	return(sol);
}
solve.EnAll = function(m, max.perm=0, n=4) {
	# generates ncol(m) * (nrow(m)!) root combinations/permutations!
	l = lapply(seq(ncol(m)), function(id) solve.En(as.vector(m[,id]), max.perm=max.perm, n=n));
	do.call(rbind, l)
}

### S3: decomposed polynomial systems
solve.S = function(S, R, b=0) {
	# generic solver (based on existing S = x+y+z)
	b2 = if(length(b) > 1) b[2] else 0; # Ext A2;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A3;
	x = sapply(S, function(x) roots(c(1, -x, R[2] - b2*x, - R[3] + b3*x)))
	len = length(S)
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	yz = R[3]/x - b3
	yz.s = S - x
	# TODO: robust (when necessary)
	# Note: this simple method is NOT robust!
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

### S3: Solve Step 2
# - but uses old/non-robust approach;
solve.mS = function(S, b=0) {
	# generic solver (based on existing S = x+y+z)
	# S = cbind(S, E2, E3)
	b2 = if(length(b) > 1) b[2] else 0; # Ext A2;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A3;
	x = sapply(seq(nrow(S)), function(id) roots(c(1, -S[id, 1], S[id, 2] - b2*S[id, 1], - S[id, 3] + b3*S[id, 1])))
	E2 = S[,2]; E3 = S[,3]; S = S[,1]; len = length(S)
	S  = matrix(S, ncol=len, nrow=3, byrow=T)
	E3 = matrix(E3, ncol=len, nrow=3, byrow=T)
	yz = E3/x - b3
	yz.s = S - x
	# TODO: robust (when necessary)
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

#############
### Other ###

perm.gen = function(x) {
	len = length(x)
	id = seq(len)
	# What was the purpose ???
	# id.m = outer(id, id, function(i, j) ((i-j) %% len + 1))
	id.m = outer(id, id, function(i, j) ((i+j+1) %% len + 1))
	p.m = x[id.m]
	dim(p.m) = dim(id.m)
	p.m
}
