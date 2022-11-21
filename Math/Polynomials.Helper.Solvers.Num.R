########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Numerical Solvers


### fast load:
# source("Polynomials.Helper.Solvers.Num.R")


### Solvers: Generic Tools
# - Numerical approaches;


####################

### Helper Functions

source("Polynomials.Helper.R")

library(rootSolve)


# Solve for all initial tuples in x0
solve.all = function(FUN, x0, ..., debug=TRUE) {
	if(is.null(dim(x0))) x0 = matrix(x0, nrow=1);
	nr = nrow(x0);
	x.all = array(0, c(ncol(x0), 0));
	#
	for(id in seq(nr)) {
		xi = rbind(Re(x0[id,]), Im(x0[id,]));
		xx = multiroot(FUN, start=xi, ...);
		#
		x = xx$root;
		x = matrix(x, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
		if(debug) {
			cat(paste0("ID = ", id, "; Iter = ", xx$iter, "; Prec = ", xx$estim.precis));
			cat("\n   "); cat.sol(xx);
		}
		x.all = cbind(x.all, x);
	}
	x.all = t(x.all);
	colnames(x.all) = paste0("x", seq(ncol(x0)));
	rownames(x.all) = NULL;
	return(x.all);
}
# Solve using given Path:
solve.path = function(FUN, x0, path, debug=TRUE) {
	n = length(path);
	for(i in seq(n)) {
		if(debug) cat(paste0("Step: ", i, "\n"));
		x0 = solve.all(FUN, x0=x0, R=path[[i]], debug=debug);
	}
	return(x0)
}

# create a seq from Xstart to Xend;
expand.path = function(xs, xe, steps=6, start.at=0) {
	if(steps <=0 ) stop("Invalid number of steps!");
	if(steps == 1) {
		warning("Only final step!");
		return(xe);
	}
	if(length(xe) == 1) xe = rep(xe, length(xs));
	#
	i  = seq(start.at, steps, by=1) / steps;
	sq = lapply(i, function(i) (1 - i)*xs + i*xe );
	return(sq);
}

### Print/Format
cat.sol = function(x, digits=4, sep="\n") {
	if(inherits(x, "list")) {
		x = x$root;
		x = matrix(x, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
	}
	cat("c(");
	cat(paste0(round(x, digits=digits), collapse=", ")); cat(")"); cat(sep);
}
cat.sol.m = function(x, digits=4, sep=",\n") {
	if( ! inherits(x, "matrix")) {
		return(cat.sol(x, digits=digits));
	}
	apply(x, 1, cat.sol, digits=digits, sep=sep);
	invisible();
}
# Print also variable names:
cat.sol.var = function(x, digits=8, sep=";\n") {
	if(inherits(x, "list")) {
		x = x$root;
		x = matrix(x, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
	}
	#
	len = length(x);
	if(len == 2) { vn = c("x", "y"); }
	else if(len == 3) { vn = c("x", "y", "z"); }
	else vn = paste0("x", seq(len));
	#
	cat(paste0(vn, " = ", round(x, digits=digits), sep), sep=""); cat("\n");
}
