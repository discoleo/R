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
solve.all = function(FUN, x0, ..., debug=TRUE, maxiter=100) {
	if(is.null(dim(x0))) x0 = matrix(x0, nrow=1);
	nr = nrow(x0);
	x.all = array(0, c(ncol(x0), 0));
	#
	for(id in seq(nr)) {
		xi = rbind(Re(x0[id,]), Im(x0[id,]));
		xx = multiroot(FUN, start=xi, ..., maxiter=maxiter);
		#
		x = xx$root;
		x = matrix(x, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
		if(debug) {
			cat(paste0("ID = ", id, "; Iter = ", xx$iter, "; Prec = ", xx$estim.precis));
			cat("\n   "); cat.sol0(x);
		}
		x.all = cbind(x.all, x);
	}
	x.all = t(x.all);
	colnames(x.all) = paste0("x", seq(ncol(x0)));
	rownames(x.all) = NULL;
	if(nr == 1) {
		# TODO: new class?
		class(x.all) = c("sol", class(x.all));
	} else {
		class(x.all) = c("sol", class(x.all));
	}
	return(x.all);
}
# Solve using given Path:
solve.path = function(FUN, x0, path, ..., debug=TRUE, maxiter=100) {
	n = length(path);
	for(i in seq(n)) {
		if(debug) cat(paste0("Step: ", i, "\n"));
		x0 = solve.all(FUN, x0=x0, R=path[[i]], ..., debug=debug, maxiter=maxiter);
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
print.sol = function(x, digits=4, sep=",\n") {
	# cat() is NOT generic!
	cat.sol(x, digits=digits, sep=sep);
}
cat.sol0 = function(x, digits=4, sep="\n") {
	if(inherits(x, "list")) {
		x = x$root;
		x = matrix(x, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
	}
	cat("c(");
	cat(paste0(round(x, digits=digits), collapse=", ")); cat(")"); cat(sep);
}
cat.sol = function(x, digits=4, sep=",\n") {
	if( ! inherits(x, "matrix")) {
		return(cat.sol0(x, digits=digits));
	}
	apply(x, 1, cat.sol0, digits=digits, sep=sep);
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

##################

### Graphics

### Plot path of roots: The Leo-Diagram
# R = parameter of function FUN;
# FUN = function with the system of NLEs;
plot.path = function(R, R0, x0, FUN, steps=20, col=seq(length(R0)), ...,
		subset=NULL, p0.pch=5, debug=FALSE) {
	isMatrix = is.matrix(x0);
	len = if(isMatrix) nrow(x0) else 1;
	# Path:
	path = expand.path(R0, R, steps=steps);
	x = array(0, c(length(R0), 0));
	for(id in seq(steps)) {
		x0 = solve.all(FUN, x0=x0, R=path[[id]], ..., debug=debug);
		x  = cbind(x, t(x0));
	}
	nr = nrow(x);
	isSubset = ! is.null(subset);
	if(nr > 1) {
		xlim = c(min(Re(x)), max(Re(x)));
		ylim = c(min(Im(x)), max(Im(x)));
		plot(Re(x[1,]), Im(x[1,]), col=col[1], xlim=xlim, ylim=ylim);
		idAll = if(isSubset) subset else seq(2, nr);
		for(id in idAll) {
			points(Re(x[id,]), Im(x[id,]), col=col[id]);
		}
		# Initial points
		if(p0.pch > 0) {
			p0 = function(x, col) {
				points(jitter(Re(x)), Im(x), col=col, pch=p0.pch, cex=1.5, lwd=1.5);
			}
			if(isSubset) {
				p0(x[- subset, seq(len)], col="blue");
				p0(x[  subset, seq(len)], col="red");
			} else {
				p0(x[, seq(len)], col="red");
			}
		}
	}
	invisible(x);
}
