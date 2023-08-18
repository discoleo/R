

### Trig: Mixed Fractions


###############

### I( 1 / (x - sin(x)) )
# Fails:
integrate(\(x) 1 / (x - sin(x)) - 6/x^3 - 3/10/x, 0, 1)
# pracma: did not finish for > 5 minutes
# pracma::integral(\(x) 1 / (x - sin(x)) - 6/x^3 - 3/10/x, 0.1, 1)

# Interval split:
mid = 0.0075; up = 1;
integrate(\(x) 1 / (x - sin(x)) - 6/x^3 - 3/10/x, mid, 10*mid, rel.tol=1E-7, subdivisions=65)$value +
	+ integrate(\(x) 1 / (x - sin(x)) - 6/x^3 - 3/10/x, 10*mid, up, rel.tol=1E-8)$value +
	+ 11/2800 * mid^2 + 17/(4*126000) * mid^4;

# Note:
# Base Int correction: - 3/x^2 + 3/10*log(x);
# - Case: up = 1 => - 3;

# TODO: ???

