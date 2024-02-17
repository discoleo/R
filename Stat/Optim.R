

### Tests using Package optimx

# Vignette: Explaining Gradient Minimizers in R
# https://cran.r-project.org/web/packages/optimx/vignettes/ExplainGradMinR.pdf


library(optimx)

# Requires also package lbfgs

# Note:
#   Modified functions;
#   for-loop => vectorized sum;
bfn = function(par) {
	n <- length(par)
	if (n < 1) stop("Variably Dimensioned: n must be positive");
	# fsum <- 0; fn1 <- 0;
	fd   = par - 1;
	fsum = sum(fd * fd);
	fn1  = sum(seq(n) * fd);
	# for (j in 1:n) {
	#	fj <- par[j] - 1; fsum <- fsum + fj * fj; fn1 <- fn1 + j * fj
	# }
	fn1_fn1 <- fn1 * fn1 # f_n+1 and f_n+2
	fsum <- fsum + fn1_fn1 + fn1_fn1 * fn1_fn1
	fsum
}

bgr = function(par) {
	n <- length(par)
	# if (n < 1) {stop("Variably Dimensioned: n must be positive") }
	# fsum <- 0; grad <- rep(0, n); fn1 <- 0;
	fd   = par - 1;
	fn1  = sum(fd * seq(n));
	grad = 2 * fd;
	# for (j in 1:n) {
	#	fj <- par[j] - 1; fn1 <- fn1 + j * fj; grad[j] <- grad[j] + 2 * fj
	# }
	fn1_2 = fn1 * 2; fn13_4 <- fn1_2 * fn1_2 * fn1;
	grad  = grad + 1:n * (fn1_2 + fn13_4)
	grad
}

n  = 100; # Number of variables
# Initial guess:
x0 = rep(pi, n)
defaultW = getOption("warn")
options(warn = -1)
mm   = c("ncg","Rcgmin","lbfgs","L-BFGS-B","Rtnmin")
res1 = opm(x0, bfn, bgr, method=mm)
options(warn = defaultW)
res1[ , seq(n, n +8)]

### Test
# All Solutions == 1
# Error of solutions:
summary(unlist(res1[1, 1:n]) - 1)


### Note: Loss of precision!
# - Methods "bcg" and "Rcgmin" loose precision with these new functions;
# - Original precision: x closer to 1 & function "value" (residual error) also better;
#  -- value for ncg = 4.247862e-24;
#  -- value for Rcgmin = 2.320238e-26;
#  -- How is this precision achieved?
#  -- Value for ncg drops massively for n = 101 with both versions of the function!
#     (to 5E-12 ... 5E-14)
# - Issue with sum(...)?

