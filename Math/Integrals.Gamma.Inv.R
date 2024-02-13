

### NLE: Gamma Function


# Strategy to solve:
# gamma(x) = y


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
plot.gamma = function(xlim = c(-6, -1), ylim = c(-1,3), hline = NULL, n = 1000) {
	curve(gamma(x), from = xlim[1], to = xlim[2], ylim=ylim, n=n);
	if( ! is.null(hline)) abline(h = hline, col = "green");
}


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
# - is more problematic, as Euler's constant < 1!
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

