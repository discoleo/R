


### Examples:

# d2y = y + b/x^n
# d2y = y + sum( b[i] / x^n[i] )
# d2y = y + 1/(n*x+1)
# d2y = y + sum( b[i] / (n[i]*x+1) )
# d2y = - y + 1/(n*x+1)
# d2y = - y + sum( b[i] / (n[i]*x+1) )


####################

### Helper Functions

library(bvpSolve)


####################
####################

### ODE:
# d2y = y - b/x^(2/3);


Ipk = function(k, lim=Inf, b=1) {
	y = integrate(\(x) sin(k*x^(3/2)) / (x^3 + 1), 0, lim)$value;
	# Normalization:
	y = b * y / (2/3 * gamma(2/3) * sin(pi/3));
	return(y);
}

Ip = function(x, y, pars) {
	b = pars$b;
	d2y = y[1] - b / x^(2/3);
	list(c(y[2], d2y));
}


lim = Inf; # fixed
b = -1/3;
k.start = 0.1; k.end = 1;
x = seq(k.start, k.end, by = 0.01)

sol <- bvpshoot(
	yini = c(Ipk(k.start, lim=lim, b=b), NA),
	yend = c(Ipk(k.end, lim=lim, b=b), NA),
	x = x, func = Ip, guess = 0, parms = list(lim=lim, b=b))

### Test

plot(sol)

# perfect match
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, lim=lim, b=b))
lines(x, y, col="red", lty=2)


#####################
#####################

### Composite Example
# Free Term: Polynomial

### Base:
# y(k) = I( sin(k*z^p) / (z^(2*p) + 1) )

### ODE:
# d2y = y - b1/x^(3/4) - b2/x^(1/5) - b3/x^(4/3);


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

### d3y = y + P(1/x)

### Example:
# d3y = y - b/x^(1/n);

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
	# Note: y[1] = y; y[2] = dy; y[3] = d2y;
	d3y = y[1] - 1 / x^(1/n);
	list(c(y[2], y[3], d3y));
}
lim = Inf;

n = 3; # n = 6;
k.start = 0.1; k.end = 1.5;
x = seq(k.start, k.end, by = 0.01)

sol <- bvpshoot(
	yini = c(Ipk(k.start, lim=lim, n=n), dyIpk(k.start, lim=lim, n=n), NA),
	yend = c(Ipk(k.end, lim=lim, n=n), NA, NA),
	x = x, func = Ip, guess = c(0), parms = list(lim=lim, n=n))

### Test

plot(sol)

# perfect match
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, lim=lim, n=n))
lines(x, y, col="red", lty=2)


######################
######################

### d3y = - y + P(1/x)

Ipk = function(k, lim=Inf, n=3, tol=1E-6) {
	r = integrate(\(x) exp(-k*x^n) / (x^(3*n) - 1), 0, 1 - tol)$value +
		integrate(\(x) exp(-k*x^n) / (x^(3*n) - 1), 1 + tol, lim)$value;
	# Normalization:
	r = r * n / gamma(1/n);
	return(r);
}
dyIpk = function(k, lim=Inf, n=3, tol=1E-6) {
	r = integrate(\(x) - x^n * exp(-k*x^n) / (x^(3*n) - 1), 0, 1 - tol)$value +
		integrate(\(x) - x^n * exp(-k*x^n) / (x^(3*n) - 1), 1 + tol, lim)$value;
	# Normalization:
	r = r * n / gamma(1/n);
	return(r);
}

Ip = function(x, y, pars) {
	n = pars$n;
	d3y = - y[1] - 1 / x^(1/n);
	list(c(y[2], y[3], d3y));
}
lim = Inf;

n = 3; # n = 2/3;
k.start = 0.1; k.end = 1.5;
x = seq(k.start, k.end, by = 0.01)

sol <- bvpshoot(
	yini = c(Ipk(k.start, lim=lim, n=n), dyIpk(k.start, lim=lim, n=n), NA),
	yend = c(Ipk(k.end, lim=lim, n=n), NA, NA),
	x = x, func = Ip, guess = c(0), parms = list(lim=lim, n=n))

### Test

plot(sol)

# perfect match
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, lim=lim, n=n))
lines(x, y, col="red", lty=2)


######################
######################

### d2y = y - 1/(n*x+1) - n/(n*x+1)^2


Ipk = function(k, n=1, lim=1) {
	r = integrate(\(x) (x^(n*k) - exp(-k)) / (n*log(x) + 1), 0, lim, rel.tol=1E-8)$value;
	return(r);
}
# dy = 1/(n*k + 1) - y;
dyIpk = function(k, n=1, lim=1) {
	r = 1/(n*k+1) - Ipk(k=k, n=n, lim=lim);
	return(r);
}

Ip = function(x, y, parms) {
	n = parms$n;
	d2y = y[1] - (1/(n*x + 1) + n/(n*x + 1)^2);
	list(c(y[2], d2y));
}
lim = 1;


###
n = 2; # n = 1/3;
k.start = 0.1; k.end = 1.5;
x = seq(k.start, k.end, by = 0.005)
# Guess:
dyIpk(k.start, n=n, lim=lim)


sol <- bvpshoot(
	yini = c(Ipk(k.start, n=n, lim=lim), NA),
	yend = c(Ipk(k.end, n=n, lim=lim), NA),
	x = x, func = Ip, guess = 0.7, parms = list(n=n, lim=lim))


### Test

plot(sol)

#
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, n=n, lim=lim))
lines(x, y, col="red", lty=2)


### Test: dy
n = 2;
k.start = 0.1; k.end = 1.5;
xguess = seq(k.start, k.end, by = 0.0005);
yguess = sapply(xguess, \(k) Ipk(k, n=n, lim=lim))
# seems a good match:
plot(xguess[-1], diff(yguess)/diff(xguess), type="l")
lines(xguess, 1/(n*xguess+1) - yguess, col="red", lty=2)
abline(h=0); abline(v=2/3, col="blue");
#
dy = 1/(n*xguess+1) - yguess;
plot(xguess[- 1], diff(dy)/diff(xguess), type="l")
lines(xguess, - 1/(n*xguess+1) - n/(n*xguess+1)^2 + yguess, col="red", lty=2)


### Varia:

# [old]
Ipk.old = function(k, n=1, lim=1) {
	r = integrate(\(x) (x^(n*k) - exp(-n*k)) / (log(x) + 1), 0, lim, rel.tol=1E-8)$value;
	r = r / n;
	# Lim: n -> 0 => I/n = k;
	return(r);
}


######################
######################

### d2y = - y + 1/(n*x+1)

# y(k) = I( x^(n*k) / (n^2*log(x)^2 + 1) )

Ipk = function(k, n=1, lim=1) {
	r = integrate(\(x) x^(n*k) / (n^2*log(x)^2 + 1), 0, lim, rel.tol=1E-8)$value;
	return(r);
}

Ip = function(x, y, parms) {
	n = parms$n;
	d2y = - y[1] + 1/(n*x + 1);
	list(c(y[2], d2y));
}
lim = 1;


###
n = 2; # n = 1/3;
k.start = 0.1; k.end = 1.5;
x = seq(k.start, k.end, by = 0.005)


sol <- bvpshoot(
	yini = c(Ipk(k.start, n=n, lim=lim), NA),
	yend = c(Ipk(k.end, n=n, lim=lim), NA),
	x = x, func = Ip, guess = 0.7, parms = list(n=n, lim=lim))


### Test

plot(sol)

#
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, n=n, lim=lim))
lines(x, y, col="red", lty=2)


####################

### Composite:
### d2y = - y + sum( b[i] / (n[i]*x+1) )

Ipk = function(k, n = NULL, b = NULL, lim=1) {
	r = sapply(n, \(n) {
		integrate(\(x) x^(n*k) / (n^2*log(x)^2 + 1), 0, lim, rel.tol=1E-8)$value;
	});
	r = sum(b * r);
	return(r);
}

Ip = function(x, y, parms) {
	n = parms$n;
	b = parms$b;
	d2y = - y[1] + sum(b/(n*x + 1));
	list(c(y[2], d2y));
}
lim = 1;


###
n = c(2,3,5); b = c(2.5, 2, -7);
k.start = 0.1; k.end = 3;
x = seq(k.start, k.end, by = 0.005)


sol <- bvpshoot(
	yini = c(Ipk(k.start, n=n, b=b, lim=lim), NA),
	yend = c(Ipk(k.end, n=n, b=b, lim=lim), NA),
	x = x, func = Ip, guess = 0.7, parms = list(n=n, b=b, lim=lim))


### Test

plot(sol)

#
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, n=n, b=b, lim=lim))
lines(x, y, col="red", lty=2)



######################
######################

### d2y = y + 1/(n*x+1)

# y(k) = I( (x^(n*k) - exp(-k)) / (n^2*log(x)^2 - 1) )

Ipk = function(k, n=1, lim=1) {
	r = integrate(\(x) (x^(n*k) - exp(-k)) / (n^2*log(x)^2 - 1), 0, lim, rel.tol=1E-8)$value;
	return(r);
}

Ip = function(x, y, parms) {
	n = parms$n;
	d2y = y[1] + 1/(n*x + 1);
	list(c(y[2], d2y));
}
lim = 1;


###
n = 2; # n = 1/3;
k.start = 0.1; k.end = 1.5;
x = seq(k.start, k.end, by = 0.005)


sol <- bvpshoot(
	yini = c(Ipk(k.start, n=n, lim=lim), NA),
	yend = c(Ipk(k.end, n=n, lim=lim), NA),
	x = x, func = Ip, guess = 0.7, parms = list(n=n, lim=lim))


### Test

plot(sol)

#
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, n=n, lim=lim))
lines(x, y, col="red", lty=2)


######################
######################

### Composite:
### d2y = y + sum( b[i] / (n[i]*x+1) )

# - theoretically its possible to add also
#   the terms: 1/(n*x + 1)^2 & 1/x^p;

Ipk = function(k, n = NULL, b = NULL, lim=1) {
	r = sapply(n, \(n) {
		integrate(\(x) (x^(n*k) - exp(-k)) / (n^2*log(x)^2 - 1), 0, lim, rel.tol=1E-8)$value;
	});
	r = sum(b * r);
	return(r);
}

Ip = function(x, y, parms) {
	n = parms$n;
	b = parms$b;
	d2y = y[1] + sum(b/(n*x + 1));
	list(c(y[2], d2y));
}
lim = 1;


###
n = c(2,3,5); b = c(-1, -3, 14);
k.start = 0.1; k.end = 2;
x = seq(k.start, k.end, by = 0.005)


sol <- bvpshoot(
	yini = c(Ipk(k.start, n=n, b=b, lim=lim), NA),
	yend = c(Ipk(k.end, n=n, b=b, lim=lim), NA),
	x = x, func = Ip, guess = 0.7, parms = list(n=n, b=b, lim=lim))


### Test

plot(sol)

#
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, n=n, b=b, lim=lim))
lines(x, y, col="red", lty=2)

