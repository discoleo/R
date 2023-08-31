


### Helper Functions


library(bvpSolve)


#####################

### Y = G(x) * I[0, x](...) + F(x)

### y = x^p * I( exp(-k*z^n) dz ) + F(x)

### D =>
# [not run]
# dy = p*x^(p-1)*I(...) + x^p * exp(-k*x^n) + df;
# dy = p*(y - f)/x + x^p * exp(-k*x^n) + df;
x*dy - p*y - x^(p+1) * exp(-k*x^n) + p*f - x*df # = 0

### D2 =>
x*d2y + dy - p*dy - ((p+1)*x^p - k*n*x^(p+n)) * exp(-k*x^n) + (p-1)*df - x*d2f # = 0
x^2*d2y - (p-1)*x*dy + (k*n*x^n - p - 1) * x^(p+1) * exp(-k*x^n) + (p-1)*x*df - x^2*d2f # = 0
x^2*d2y - (p-1)*x*dy + (k*n*x^n - p - 1) * (x*dy - p*y  + p*f - x*df) + (p-1)*x*df - x^2*d2f # = 0
# ODE:
x^2*d2y + (k*n*x^n - 2*p)*x*dy - p*(k*n*x^n - p - 1)*y +
	 + (k*n*x^n - p - 1) * (p*f - x*df) + (p-1)*x*df - x^2*d2f # = 0

### Special Cases:

# TODO


### Tests:

# F(x) = x^2 + x + 1;

Ipk = function(x, parms) {
	p = parms$p; n = parms$n; k = parms$k; up = x;
	(x^2 + x + 1) + x^p * integrate(\(x) exp(-k * x^n), 0, up)$value;
}
dyf = function(x, parms) {
	p = parms$p; n = parms$n; k = parms$k;
	f = x^2 + x + 1; df = 2*x + 1;
	y = Ipk(x, parms);
	dy = p*(y - f)/x + x^p * exp(-k*x^n) + df;
	return(dy);
}

Ip = function(x, y, parms) {
	p  = parms$p; n = parms$n; k = parms$k;
	dy = y[2]; y = y[1];
	f  = x^2 + x + 1; df = 2*x + 1; d2f = 2;
	d2y = (k*n*x^n - 2*p)*x*dy - p*(k*n*x^n - p - 1)*y +
		+ (k*n*x^n - p - 1) * (p*f - x*df) + (p-1)*x*df - x^2*d2f;
	d2y = - d2y / x^2;
	list(c(dy, d2y));
}


###
parms = list(p = 1/2, k = 1/3, n = 4/3);
x.start = 0.1; x.end = 3; # x.end = 1.5;
x = seq(x.start, x.end, by = 0.01)

sol <- bvpshoot(
	yini = c(Ipk(x.start, parms), NA),
	yend = c(Ipk(x.end, parms), NA),
	x = x, func = Ip, guess = dyf(x.start, parms),
	parms = parms)

### Test

plot(sol)

# perfect match
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(x) Ipk(x, parms))
lines(x, y, col="red", lty=2)

