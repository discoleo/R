

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
		if(length(R) == 9) {
			err1 = err1 - R[1:3];
			err2 = err2 - R[4:6];
			err3 = err3 - R[7:9];
		} else {
			err1 = err1 - R[1];
			err2 = err2 - R[2];
			err3 = err3 - R[3];
		}
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

solve.S3x3.byPath = function(x0, R, n = c(1,2), maxiter = 32) {
	xm = matrix(x0, ncol=3)
	R0 = test.S3x3.Simple(xm)
	Rp = expand.path(R0, rep(R, each=3))
	# Solver: fails directly for some starting values;
	sol = solve.path(wrap(test.S3x3.Simple), x0, path=Rp, n=n)
	x0 = round(sol[1,], 4);
	xc = rbind(Re(x0), Im(x0));
	#
	sol = multiroot(wrap(test.S3x3.Simple), start = xc, R=R, n=n, maxiter=maxiter);
	x0 = matrix(sol$root, nrow = 2);
	x0 = x0[1,] + 1i*x0[2,];
	x0 = matrix(x0, nrow=3);
	return(x0)
}

n = c(1,2)
R = c(1,2,3)
#
x0 = c(1.0758-0.6126i, 1.2275+0.5367i, 2e-04+1.0003i,
	-1.1516+0i, -0.2277+0.4636i, -0.2277-0.4636i,
	1.0758+0.6126i, 2e-04-1.0003i, 1.2275-0.5367i);
# x0 = seq(9) / 11 + seq(9,1) * 1i/3;
# x0 = c(-5,-1, seq(6), -6) / 3 + seq(-3,5) * 1i/3;
#
x = solve.S3x3.byPath(x0, R=R, n=n)

test.S3x3.Simple(x, n=n, round0 = TRUE)
prod(x)

# there seems to be only 2 distinct solution-sets?

# xall = c(xall, x);
# poly.calc0(xall) * 7*8

56*x^18 - 336*x^17 + 672*x^16 - 163*x^15 +
 + 129*x^14 - 6672*x^13 + 19268.72*x^12 # ??? + ...

### Derivation:

n = 1;
p1 = as.pm("x1^n + x2^n + x3^n - R1")
p2 = as.pm("y1^n + y2^n + y3^n - R1")
p3 = as.pm("z1^n + z2^n + z3^n - R1")
#
p4 = as.pm("x1^2 + y1^2 + z1^2 - R2")
p5 = as.pm("x2^2 + y2^2 + z2^2 - R2")
p6 = as.pm("x3^2 + y3^2 + z3^2 - R2")
#
p7 = as.pm("x1*y2 + y2*z3 + z3*x1 - R3")
p8 = as.pm("x2*y3 + y3*z1 + z1*x2 - R3")
p9 = as.pm("x3*y1 + y1*z2 + z2*x3 - R3")

pR = solve.lpm(p1,p2,p3, p4,p5,p6, p7,p8,p9, by = c("x3","y3","z3"))
pR = pR[3:8]
pR[[1]] = pR[[1]]$Rez;
pR = pR[c(6,5,4,1,2,3)]
pR$by = c("z2","z1")
pR = do.call(solve.lpm, pR)
pR = pR[2:5]
pR[[1]] = pR[[1]]$Rez

sapply(pR, function(p) max(p$y2))
# but last poly has already 272 monomials;


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

