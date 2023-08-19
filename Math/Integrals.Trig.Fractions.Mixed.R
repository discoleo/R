

### Trig: Mixed Fractions

int = function(FUN, upper, lower, rel.tol=1E-8, subdivisions=4097, ...) {
	r = integrate(FUN, upper, lower, subdivisions=subdivisions, rel.tol=rel.tol, ...);
	cat(paste0("Err = ", r$abs.error, "; subd = ", r$subdivisions, "\n"));
	return(r$value);
}

###############

### I( 1 / (x - sin(x)) )
# Fails:
integrate(\(x) 1 / (x - sin(x)) - 6/x^3 - 3/10/x, 0, 1)
# pracma: did not finish for > 5 minutes
# pracma::integral(\(x) 1 / (x - sin(x)) - 6/x^3 - 3/10/x, 0, 1)

# Interval split:
mid = 0.0075; up = 1;
integrate(\(x) 1 / (x - sin(x)) - 6/x^3 - 3/10/x, mid, 10*mid, rel.tol=1E-7, subdivisions=65)$value +
	+ integrate(\(x) 1 / (x - sin(x)) - 6/x^3 - 3/10/x, 10*mid, up, rel.tol=1E-8)$value +
	+ 11/2800 * mid^2 + 17/(4*126000) * mid^4;

# Note:
# Base Int correction: - 3/x^2 + 3/10*log(x);
# - Case: up = 1 => - 3;

# TODO: ???


###############

### on [0, Inf]
FUN = \(x) 1 / (x - sin(x)) - 6/x^3 - 3/10/x - 7/10/(x+1);
mid = 0.01;
int(FUN, mid, 1, rel.tol=1E-7, subdiv=8025) +
int(FUN, 1, 1E+5, rel.tol=1E-7, subdiv=8025) +
int(FUN, 1E+5, 1E+6, rel.tol=1E-6, subd=1025) +
int(FUN, 1E+6, 5E+6, rel.tol=1E-6, subd=1025) +
+ 11/2800 * mid^2 + 17/(4*126000) * mid^4;

# TODO: ???


###
integrate(\(x) 1 / (x - sin(pi/2*x)) - 1/(1-pi/2)/x - 1/(x-1), 0, 1)

# TODO: ???

