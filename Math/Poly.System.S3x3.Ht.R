

### Complex Symmetries
### S 3x3 Ht


####################

### Helper Functions

# library(rootSolve)

source("Polynomials.Helper.Solvers.Num.R")

###
test.S3x3.Simple = function(sol, R = NULL, n = c(1,2), round0 = FALSE, warn = TRUE) {
	if(warn) {
		warn.S3x3(n=n, R=R);
	}
	if( ! inherits(sol, "matrix")) {
		sol = matrix(sol, nrow = 3, ncol = 3, byrow = FALSE);
	}
	#
	err1 = apply(sol^n[1], 1, sum);
	err2 = apply(sol^n[2], 2, sum);
	#
	tmp = lapply(seq(3), function(i) {
		nc = if(i == 1) c(1,2,3)
			else if(i == 2) c(2,3,1) else c(3,1,2);
		sol[i, nc];
	});
	tmp = do.call(rbind, tmp);
	err3 = apply(tmp, 2, function(x) {
		x[1]*x[2] + x[3]*(x[1] + x[2]);
	});
	if(! is.null(R)) {
		err1 = err1 - R[1];
		err2 = err2 - R[2];
		err3 = err3 - R[3];
	}
	err = cbind(err1, err2, err3);
	if(round0) err = round0(err);
	return(err);
}
warn.S3x3 = function(n, R = null) {
	if(n[1] == n[2]) {
		if( ! is.null(R)) {
			if(R[1] != R[2]) warning("No solution!");
		}
	}
}

### Complex Solutions
wrap.D2 = function(FUN) {
	FUN = FUN;
	wF = function(x, start, R = NULL, ...) {
		if(is.null(R)) R = start;
		x = matrix(x, nrow=2); xc = x[2,]; x = x[1,] + 1i * xc;
		R = matrix(R, nrow=2); Rc = R[2,]; R = R[1,] + 1i * Rc;
		y = FUN(x, R=R, ...);
		y = rbind(Re(y), Im(y));
		return(y);
	}
	return(wF);
}

as.variables = function(x, n=4) {
	tmp = which(x != 0);
	id1 = (tmp - 1) %/% n + 1;
	id2 = (tmp - 1) %% n + 1;
	ss = letters[seq(n)];
	sn = ss[id1];
	sn = paste0(sn, id2);
	return(sn);
}


#####################

### System 3x3 Ht

# x1^n + x2^n + x3^n = R1
# y1^n + y2^n + y3^n = R1
# z1^n + z2^n + z3^n = R1
#
# x1^2 + y1^2 + z1^2 = R2
# x2^2 + y2^2 + z2^2 = R2
# x3^2 + y3^2 + z3^2 = R2
#
# x1*y2 + y2*z3 + z3*x1 = R3
# x2*y3 + y3*z1 + z1*x2 = R3
# x3*y1 + y1*z2 + z2*x3 = R3

### Note:
# If ((x1,x2,x3), (y1,y2,y3), (z1,z2,z3)) is a solution:
# - then any simultaneous cyclic permutation of
#   the rows or the columns or the diagonals are also a solution;
# => 9 permutations;
# - this particular type & n = even: (-1) * solutions is also solution;

# If n = 2: => R1 = R2, but system is indeterminate;
# Aux eq: x1*y1*z1 + x2*y2*z2 + x3*y3*z3 = R4

# TODO: solve;

### Debug

R = c(1,2,3)
x0 = seq(9) / 11 + seq(9,1) * 1i/3;
#
xm = matrix(x0, ncol=3)
R0 = test.S3x3.Simple(xm)
Rp = expand.path(R0, rep(R, each=3))

# Solver: fails directly;
sol = solve.path(wrap(test.S3x3.Simple), x0, path=Rp)
x0 = c(1.0758-0.6126i, 1.2275+0.5367i, 2e-04+1.0003i,
	-1.1516+0i, -0.2277+0.4636i, -0.2277-0.4636i,
	1.0758+0.6126i, 2e-04-1.0003i, 1.2275-0.5367i);
#
xc = rbind(Re(x0), Im(x0));

maxiter = 32;
sol = multiroot(wrap(test.S3x3.Simple), start = xc, R=R, n = c(1,2), maxiter=maxiter);
x0 = matrix(sol$root, nrow = 2);
x0 = x0[1,] + 1i*x0[2,];

test.S3x3.Simple(x0, n = c(1,2), round0 = TRUE)


#####################
#####################

### Matrix
### Complex Symmetries

### Size = 3
m = matrix(0L, nrow=9, ncol=9)
for(i in seq(3)) m[i, c(1,2,3) + (i-1)*3] = 1L;
for(i in seq(3)) m[i + 3, c(0,3,6) + i] = 1L;
for(i in seq(3)) m[i + 6, c(i, i %% 3 + 4, (i+1) %% 3 + 7)] = 1L;

det(m)


### Size = 4
m = matrix(0L, nrow=16, ncol=16)
for(i in seq(4)) m[i, c(1,2,3,4) + (i-1)*4] = 1L;
for(i in seq(4)) m[i + 4, c(0,4,8,12) + i] = 1L;
for(i in seq(4)) m[i + 8, c(i, i %% 4 + 5, (i+1) %% 4 + 9, (i+2) %% 4 + 13)] = 1L;
for(i in seq(4)) m[i + 12, c((i+2) %% 4 + 1, (i+1) %% 4 + 5, i %% 4 + 9, (i-1) %% 4 + 13)] = 1L;

det(m)

tmp = apply(m, 1, function(x) paste0(as.variables(x), collapse = " + "))

