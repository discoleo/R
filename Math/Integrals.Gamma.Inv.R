

### NLE: Gamma Function


# Strategy to solve:
# gamma(x) = y


### Further reading:
# 1) Package optimx
# - Vignette: Explaining Gradient Minimizers in R
# https://cran.r-project.org/web/packages/optimx/vignettes/ExplainGradMinR.pdf
# 2) Package nlstools
# - focused on fitting NLS models;
# https://cran.r-project.org/web/packages/nlstools/vignettes/vignetteJSS.pdf
#
# TODO


####################

### Helper Functions

library(rootSolve)


# used as Example:
Euler   = 0.57721566490153286060651209008240243079;


gamma.inv = function(x, start = 0.5, positive = TRUE, atol = 1E-10) {
	r = rootSolve::multiroot(\(x, parms) { gamma(x) - parms$v; },
		start = start, parms = list(v = x), positive=positive, atol=atol, rtol=atol);
	cat("Root = ", format(r$root, digits = 10), "\n");
	cat("Diff = ", r$f.root, "\n");
	return(r$root);
}
gamma.invN = function(x, x0, lim, ...,
		ndeps = 1E-5, rtol = 1E-9, gtol = 0.0001, digits = 10) {
	v = x;
	if(length(lim) == 1) lim = c(lim - 1 + gtol, lim - gtol);
	cntr = c(list(...), ndeps = ndeps, pgtol = rtol);
	# abs(): should provide greater precision than ()^2 when computing the gradients;
	# param ndeps: needed to improve accuracy;
	r = optim(x0, \(x) { abs(gamma(x) - v); }, lower = lim[1], upper =lim[2],
		method = "L-BFGS-B", control = cntr);
	if( ! is.null(digits)) print(r$par, digits);
	return(r$par);
}

### Plot
plot.gamma = function(xlim = c(-6, -1), ylim = c(-1,3),
		hline = NULL, col.hline = "#00FF00A0", n = 1000) {
	if(length(xlim) < 2) xlim = c(xlim, xlim + 1);
	# Disconnected segments: manual processing
	# TODO: x[1] > x[2];
	x  = floor(xlim);
	xS = seq(x[1], x[2]);
	xE = xS + 1;
	xS[1] = xlim[1];
	xE[length(xE)] = xlim[2];
	#
	len = length(xS);
	if(xS[len] >= xE[len]) {
		xS = xS[-len]; xE = xE[-len];
		len = len - 1;
	}
	n = round(n / len);
	#
	curve(gamma(x), from = xS[1], to = xE[1], xlim=xlim, ylim=ylim, n=n);
	#
	for (i in seq_along(xS)[-1]) {
		curve(gamma(x), from = xS[i], to = xE[i], add = TRUE, n=n);
	}
	if( ! is.null(hline)) abline(h = hline, col = col.hline);
}

# see also prototype function by Bill Dunlap:
# https://stat.ethz.ch/pipermail/r-help/2024-February/478897.html


#####################
#####################

### Example 1:
# gamma(x) = pi

x = gamma.inv(pi)
# Test:
gamma(x) - pi


##############

### Example 2:
# gamma(x) = Euler

# Note:
# - is more problematic, as Euler's constant < min(gamma( on (0, Inf) ))!
# - gamma.inv fails even with a good starting value!
# gamma.inv(Euler)
# gamma.inv(Euler, start = - 3.3)

# there are 2 real solutions in the interval (-4, -3):
plot.gamma(hline = Euler)

### Solution using optim:
# - explicitly specifying x0 and the valid interval;
x = gamma.invN(Euler, -3.5, lim = -3)
gamma(x) - Euler
x = gamma.invN(Euler, -3.9, lim = -3)
gamma(x) - Euler


### Example 3:

### Min value of gamma() on (-4, -3)
x = gamma.invN(0, -3.635, c(-4,-3))
gamma(x)

### Min on (0, Inf)
x = gamma.invN(0, 1/2, c(0,2))
gamma(x)
gamma(x) - gamma(x - 2E-8)
gamma(x) - gamma(x + 2E-8)

plot.gamma(xlim = c(0,2), hline = gamma(x))
abline(v = x, col = "red")

