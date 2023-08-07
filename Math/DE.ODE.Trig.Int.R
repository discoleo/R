

library(bvpSolve)


##############


### y(k) = I( z^p * sin(k*z) )
# Note: on [0, pi/2]
# d2y = (p+1)*(p+2) * x^2 * y +
#    - (p+2)*(pi/2)^(p+1) * sin(pi/2 * x) / x^2 +
#    + (pi/2)^(p+2) * cos(pi/2 * x) / x;


Ipk = function(p, k) integrate(\(x) x^p * sin(k*x), 0, pi/2)$value;

Ip = function(x, y, pars) {
	d2y = (p+2)*((p+1)*y[1] - (pi/2)^(p+1) * sin(x*pi/2)) / x^2 +
		+ (pi/2)^(p+2) * cos(x*pi/2) / x;
	list(c(y[2], d2y));
}


###
p = 1/2
k.start = 0.1; k.end = 1;
x = seq(k.start, k.end, by = 0.01)

sol <- bvpshoot(
	yini = c(Ipk(p, k.start), NA),
	yend = c(Ipk(p, k.end), NA),
	x = x, func = Ip, guess = 0)

### Test

plot(sol)

# perfect match
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(p, k))
lines(x, y, col="red", lty=2)


######################
######################

### y = I( sin(k*z) / (z^2 + 1) )
# d2y = y + cos(k*x)/x - 1/x;

Ipk = function(k, lim=pi/2) integrate(\(x) sin(k*x) / (x^2 + 1), 0, lim)$value;

Ip = function(x, y, pars) {
	d2y = y[1] + cos(x * pars$lim)/x - 1/x;
	list(c(y[2], d2y));
}


###
lim = pi;
k.start = 0.1; k.end = 1;
x = seq(k.start, k.end, by = 0.01)

sol <- bvpshoot(
	yini = c(Ipk(k.start, lim=lim), NA),
	yend = c(Ipk(k.end, lim=lim), NA),
	x = x, func = Ip, guess = 0, parms = list(lim=lim))

### Test

plot(sol)

# perfect match
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, lim=lim))
lines(x, y, col="red", lty=2)


######################

### y = I( sin(k*z)^2 / (z^2 + 1) )
# d2y = 4*y + sin(2*c*x)/x - 2*atan(c);

Ipk = function(k, lim=pi/2) integrate(\(x) sin(k*x)^2 / (x^2 + 1), 0, lim)$value;

Ip = function(x, y, pars) {
	c = pars$lim;
	# d2y = 4*y[1] - 4*c*sin(c*x)^2 - 2*c*cos(2*c*x) + sin(2*c*x)/x + 2*c - 2*atan(c);
	d2y = 4*y[1] + sin(2*c*x)/x - 2*atan(c);
	list(c(y[2], d2y));
}


###
lim = 5/4 * pi;
k.start = 0.1; k.end = 1;
x = seq(k.start, k.end, by = 0.01)

sol <- bvpshoot(
	yini = c(Ipk(k.start, lim=lim), NA),
	yend = c(Ipk(k.end, lim=lim), NA),
	x = x, func = Ip, guess = 0, parms = list(lim=lim))

### Test

plot(sol)

# perfect match
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, lim=lim))
lines(x, y, col="red", lty=2)


######################

### y = I( sin(k*z^(3/2)) / (z^3 + 1) )
# d2y = y + ???
# TODO: solve vs proof of concept;

Ipk = function(k, lim=pi/2) integrate(\(x) sin(k*x^(3/2)) / (x^3 + 1), 0, lim)$value;

Ip = function(x, y, pars) {
	int = integrate(\(z) sin(x*z^(3/2)), 0, pars$lim)$value;
	d2y = y[1] - int;
	list(c(y[2], d2y));
}

# TODO: complicated;
# upper = Inf may work;


###
lim = 5/4 * pi;
k.start = 0.1; k.end = 1;
x = seq(k.start, k.end, by = 0.01)

sol <- bvpshoot(
	yini = c(Ipk(k.start, lim=lim), NA),
	yend = c(Ipk(k.end, lim=lim), NA),
	x = x, func = Ip, guess = 0, parms = list(lim=lim))

### Test

plot(sol)

# perfect match
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, lim=lim))
lines(x, y, col="red", lty=2)


#################

lim = Inf;
Ip = function(x, y, pars) {
	d2y = y[1] - 2/3 * gamma(2/3) * sin(pi/3) / x^(2/3);
	list(c(y[2], d2y));
}


#################

### Composite Example
# Free Term: Polynomial

Ipk = function(k, lim=Inf) {
	if(! is.infinite(lim)) warning("Range must be infinite!");
	# numerical issues:
	r1 = pracma::quadinf(\(x) sin(k*x^(4/3)) / (x^(8/3) + 1), 0, Inf)$Q;
	r2 = integrate(\(x) sin(k*x^(5)) / (x^10 + 1), 0, Inf)$value;
	r3 = pracma::quadinf(\(x) sin(k*x^(3/4)) / (x^(3/2) + 1), 0, Inf)$Q;
	# Normalization:
	r1 = r1 / (3/4 * gamma(3/4) * sin(pi*3/8));
	r2 = r2 / (1/5 * gamma(1/5) * sin(pi/10));
	r3 = r3 / (4/3 * gamma(4/3) * sin(pi*2/3));
	# just some factor:
	b = c(2, 5, -2);
	r = sum(b * c(r1, r2, r3));
	return(r);
}

Ip = function(x, y, pars) {
	# Note: same factor as above:
	b = c(2, 5, -2)
	d2y = y[1] - b[1] / x^(3/4) - b[2] / x^(1/5) - b[3] / x^(4/3);
	list(c(y[2], d2y));
}
lim = Inf;


###
k.start = 0.1; k.end = 2;
x = seq(k.start, k.end, by = 0.01)

sol <- bvpshoot(
	yini = c(Ipk(k.start, lim=lim), NA),
	yend = c(Ipk(k.end, lim=lim), NA),
	x = x, func = Ip, guess = 0, parms = list(lim=lim))

### Test

plot(sol)

# almost "perfect" match
# - possible numerical instabilities in computing the quasi-exact integral;
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, lim=lim))
lines(x, y, col="red", lty=2)


######################
######################

### y(k) = k * y0

### y = k * I( sin(k*z^(5/3)) / (z^(10/3) + 1) )
# x*d2y = 2*dy + (x^2 - 2)*y - x^(1/3);

# includes dy term;

Ipk = function(k, lim=Inf, subdivisions=4097) {
	r = integrate(\(x) k * sin(k*x^(5/3)) / (x^(10/3) + 1), 0, lim,
		subdivisions=subdivisions)$value;
	r = r / (3/5 * gamma(3/5) * sin(pi/2*3/5));
	return(r);
}

Ip = function(x, y, pars) {
	d2y = 2*y[2]/x + (x - 2/x)*y[1] - 1/x^(2/3);
	list(c(y[2], d2y));
}
lim = Inf;


###
lim = Inf;
k.start = 0.1; k.end = 1.5;
x = seq(k.start, k.end, by = 0.01)

sol <- bvpshoot(
	yini = c(Ipk(k.start, lim=lim), NA),
	yend = c(Ipk(k.end, lim=lim), NA),
	x = x, func = Ip, guess = 0, parms = list(lim=lim))

### Test

plot(sol)

# perfect match
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, lim=lim))
lines(x, y, col="red", lty=2)


######################
######################

### y(k) = 1/k * y0

# TODO


######################
######################

### d3y = y + P(x)

Ipk = function(k, lim=Inf, n=3) {
	r = integrate(\(x) exp(-k*x^n) / (x^(3*n) + 1), 0, lim)$value;
	# Normalization:
	r = r * n / gamma(1/n);
	return(r);
}
dyIpk = function(k, lim=Inf, n=3) {
	r = integrate(\(x) - x^n * exp(-k*x^n) / (x^(3*n) + 1), 0, lim)$value;
	# Normalization:
	r = r * n / gamma(1/n);
	return(r);
}

Ip = function(x, y, pars) {
	n = pars$n;
	d3y = y[1] - 1 / x^(1/n);
	list(c(y[2], y[3], d3y));
}
lim = Inf;

n = 3;
k.start = 0.1; k.end = 1.5;
x = seq(k.start, k.end, by = 0.01)

sol <- bvpshoot(
	yini = c(Ipk(k.start, lim=lim, n=n), dyIpk(k.start, lim=lim, n=n), NA),
	yend = c(Ipk(k.end, lim=lim, n=n), NA, NA),
	x = x, func = Ip, guess = c(0.1), parms = list(lim=lim, n=n))

### Test

plot(sol)

# perfect match
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, lim=lim, n=n))
lines(x, y, col="red", lty=2)

