

### Simulation of "Noisy" ODEs


####################

### Helper Functions

library(bvpSolve)


####################

stoch.gen = function(n, lim, scale = 1, ...) {
	x = seq(lim[1], lim[2], length.out = n);
	y = rnorm(n, ...) * scale;
	splinefun(x, y, method = "fmm");
}

Idny = function(x, y, parms) {
	# a relatively simple ODE:
	d2y = y[2] + x * y[1] + parms$FUN(x);
	list(c(y[2], d2y));
}
solve.ODE = function(x, init, ODE, FUN, guess = 0.5) {
	bvpshoot(
		yini = c(init[1], NA),
		yend = c(init[2], NA),
		x = x, func = ODE, guess = guess, parms = list(FUN = FUN));
}

solve.ODE.stoch = function(x, ODE, xlim, n = 10, iter = 10, ...) {
	sol = solve.ODE(x, xlim, ODE = ODE, FUN = stoch.gen(n, xlim, ...));
	# Plot:
	ylim = range(sol[,2]) + c(-2,2);
	plot(sol[, 1:2], type="l", ylim = ylim, col = "#FF0000");
	for(i in seq(iter)) {
		sol = solve.ODE(x, xlim, ODE = ODE, FUN = stoch.gen(n, xlim, ...));
		val = 15 * i;
		hex = as.hexmode(val);
		if(val < 16) hex = paste0("0", hex);
		lines(sol[, 1:2], type="l", lty = 2, col = paste0("#F032", hex));
	}
	FUN = function(x) 0; # or return the mean?
	sol = solve.ODE(x, xlim, ODE = ODE, FUN = FUN);
	lines(sol[, 1:2], type="l", lwd = 2, col = "green");
	return(invisible(sol));
}


###
init = c(10, 1);
x.start = 0; x.end = 5;
x = seq(x.start, x.end, by = 0.01)
n = 10;

sol = solve.ODE.stoch(x, ODE = Idny, init, n=n, sd = 1.5)

sol = solve.ODE.stoch(x, ODE = Idny, init, n=n, sd = 1.5, mean = 2)


### Test

plot(sol)


