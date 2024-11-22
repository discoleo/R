
### Helper Functions

polylog2 = function (z, n = 2) {
	# Bug in pracma::polylog;
	stopifnot(is.numeric(n));
	if(is.complex(z)) {
		if(Im(z) == 0) { z = Re(z); }
		else {
			warning("Complex not yet implemented!");
			if(abs(abs(z) - 1) > 1E-8) warning("Result will be incorrect!");
			if(n == 2) {
				y  = - (pi^2/3 + Re(log(-z)^2)) / 4;
				# TODO: complex-part properly;
				thz = Im(log(z));
				yi = integrate(\(x) Re(log(1- exp(1i*x))), 0, thz)$value;
				# thz = Im(log(z))/2;
				# yi = integrate(\(x) 2*log(abs(2*sin(x))), 0, thz)$value;
				y  = y - 1i*yi;
			} else if(n == 3) {
				# TODO: compute real-part;
				x = log(z) / (2*pi);
				y = (4*pi^3/3 * (- x^3 + 3i/2*x^2 + x/2)) / 2;
			} else stop("Not yet implemented!");
			return(y);
		}
	}
	if(z == 1) {
		return(pracma::zeta(n));
	}
	if(z == -1) {
		return(pracma::zeta(n) * (1/2^(n-1) - 1));
	}
	if(z < 0) {
		if(n == 2) {
			y = polylog2(z^2, n=n)/2 - polylog2(-z, n=n);
		} else {
			y = polylog2(z^2, n=n) * 2^(1-n) - polylog2(-z, n=n);
		}
		return(y);
	}
	if(z > 1) {
		if(n == 2) {
			y = - pi^2/6 - log(0i - z)^2/2 - polylog2(1/z, n=n);
		} else if(n == 3) {
			x = log(z) / (2*pi);
			y = 4*pi^3/3 * (- x^3 + 3i/2*x^2 + x/2) + polylog2(1/z, n=n);
		} else stop("Not yet implemented!");
		return(y);
	} else if(z >= 0.55) {
		# pracma::polylog FAILS!
		if(n == 2) {
			y = pi^2/6 - log(z)*log(1-z) - pracma::polylog(1 - z, n=n);
		} else if(n == 3) {
			# see Ref;
			y = log(1-z)^3 / 6 - log(z)*log(1-z)^2/2 + pi^2/6 * log(1-z) +
				+ pracma::zeta(3) - polylog(1-z, 3) - polylog2(z/(z-1), n=3);
		} else stop("Not yet implemented!");
		return(y);
	}
	return(pracma::polylog(z, n=n));
}
# Ref:
# https://math.stackexchange.com/questions/942796/compute-polylog-of-order-3-at-frac12


### Tests

if(FALSE) {

# Test: 1 / (0.6 / (1-0.6)) == 2/3 > 0.55;
# - requires 2 iterations;
li3_06 = 0.6560025136329806832346611928113322291802755981380034005592;
polylog2(0.6, 3) - li3_06;


### Test:
x = 1/3; # value behaves badly!
# Note:
# - in practice, one may know only the z value and NOT x;
z = exp(1i*x);
integrate(\(x) Re(log(1- exp(1i*x))), 0, x) # OK
# Approximations: both fail!
thz = Im(log(z)); # recompute x;
integrate(\(x) Re(log(1- exp(1i*x))), 0, thz)
#
thz = Im(log(z))/2;
integrate(\(x) 2*log(abs(2*sin(x))), 0, thz)
}

