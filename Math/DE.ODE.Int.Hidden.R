


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
	- x^2*d2f - (k*n*x^n - 2*p)*x*df + p*(k*n*x^n - p - 1)*f # = 0

### Special Cases:

# TODO


### Tests:

# F(x) = x^2 + x + 1;
# F(x) = - 1/3 * x - x^(1/2) + 1;

Ipk = function(x, parms) {
	p = parms$p; n = parms$n; k = parms$k; up = x;
	- 1/3 * x - x^(1/2) + 1 + # (x^2 + x + 1) +
	x^p * integrate(\(x) exp(-k * x^n), 0, up)$value;
}
dyf = function(x, parms) {
	p = parms$p; n = parms$n; k = parms$k;
	# f = x^2 + x + 1; df = 2*x + 1;
	f = - 1/3 * x - x^(1/2) + 1; df = - 1/3 - 1/(2*x^(1/2));
	y = Ipk(x, parms);
	dy = p*(y - f)/x + x^p * exp(-k*x^n) + df;
	return(dy);
}

Ip = function(x, y, parms) {
	p  = parms$p; n = parms$n; k = parms$k;
	dy = y[2]; y = y[1];
	# f  = x^2 + x + 1; df = 2*x + 1; d2f = 2;
	f  = - 1/3 * x - x^(1/2) + 1; df = - 1/3 - 1/(2*x^(1/2)); d2f = 1/(4*x^(3/2));
	d2y = (k*n*x^n - 2*p)*x*dy - p*(k*n*x^n - p - 1)*y +
		+ (k*n*x^n - p - 1) * (p*f - x*df) + (p-1)*x*df - x^2*d2f;
	d2y = - d2y / x^2;
	list(c(dy, d2y));
}


###
parms = list(p = 1/2, k = 1/3, n = 4/3);
x.start = 0.1; x.end = 10; # x.end = 1.5;
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


#####################

### y = x^p * I( log(z^n + 1) / (z + 1) dz ) + F(x)

### D =>
# [not run]
# dy = p*x^(p-1)*I(...) + x^p * log(x^n + 1) / (x + 1) + df;
# dy = p*(y - f)/x + x^p * log(x^n + 1) / (x + 1) + df;
x*dy - p*y - x^(p+1) * log(x^n + 1) / (x + 1) + p*f - x*df # = 0
(x^2 + x)*dy - p*(x + 1)*y - x^(p+1) * log(x^n + 1) + p*(x+1)*f - (x^2 + x)*df # = 0

# TODO


######################

### y = exp(k1*x) * I( exp(-k2*z^n) dz ) + F(x)

### D =>
# [not run]
# dy = k1*exp(k1*x) * I(...) + exp(-k2*x^n + k1*x) + df;
# dy = k1*(y - f) + exp(-k2*x^n + k1*x) + df;
dy - k1*y - exp(-k2*x^n + k1*x) + k1*f - df # = 0


### D2 =>
d2y - k1*dy + (n*k2*x^(n-1) - k1) * exp(-k2*x^n + k1*x) + k1*df - d2f # = 0
d2y - k1*dy + (n*k2*x^(n-1) - k1) * (dy - k1*y + k1*f - df) + k1*df - d2f # = 0
d2y + (n*k2*x^(n-1) - 2*k1)*dy - k1*(n*k2*x^(n-1) - k1)*y +
	+ k1*(n*k2*x^(n-1) - k1)*f - (n*k2*x^(n-1) - 2*k1)*df - d2f # = 0

### Special Cases:

### F(x) = b*x
d2y + (n*k2*x^(n-1) - 2*k1)*dy - k1*(n*k2*x^(n-1) - k1)*y +
	+ b*n*k1*k2*x^n - b*n*k2*x^(n-1) - b*k1^2*x + 2*b*k1 # = 0

# TODO: check;

