#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals:
## Derived from I( Trig-Fractions )
##
## v.0.2b


### Examples:

# d2y = y + b/x^n
# d2y = y + sum( b[i] / x^n[i] )
# d2y = y + 1/(n*x+1)
# d2y = y + sum( b[i] / (n[i]*x+1) )
# d2y = - y + 1/(n*x+1)
# d2y = - y + sum( b[i] / (n[i]*x+1) )

### Exp:
# d2y = 2*dy - b*exp(x) / x^(2/3)


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

### y(x) = x * y0

### y = x * I( sin(x*z^(5/3)) / (z^(10/3) + 1) )
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

### y(x) = 1/x * y0

# TODO


######################

### y(x) = exp(x) * y0

### ODE:
# d2y = 2*dy - b*exp(x) / x^(2/3);

# TODO: incomplete Gamma function;

# y(k)
Ipk = function(k, lim=Inf, b=1) {
	y = integrate(\(x) sin(k*x^(3/2)) / (x^3 + 1), 0, lim)$value;
	y = exp(k) * y;
	# Normalization:
	y = b * y / (2/3 * gamma(2/3) * sin(pi/3));
	return(y);
}
# dy(k)
Idy = function(k, lim=Inf, b=1) {
	# Note: various numerical issues!
	if(lim == Inf) lim = 2000 * pi;
	dy = integrate(\(x) {
		xx = x^(2/3); # x = z^(2/3); dx = 2/3 * z^(-1/3) dz;
		xx * cos(k*x) / (x*x + 1);
	}, 0, lim, subdivisions = 43570)$value;
	# Normalization:
	# dy = 2/3 * b * dy / (2/3 * gamma(2/3) * sin(pi/3));
	dy = b * dy / (gamma(2/3) * sin(pi/3));
	dy = dy * exp(k) + Ipk(k, lim=lim, b=b);
	return(dy);
}

Ip = function(x, y, pars) {
	b = pars$b;
	d2y = 2*y[2] - b * exp(x) / x^(2/3);
	list(c(y[2], d2y));
}


lim = Inf; # fixed
b = 1/3;
k.start = 0.4; k.end = 1.90;
x = seq(k.start, k.end, by = 0.01)

sol <- bvpshoot(
	yini = c(Ipk(k.start, lim=lim, b=b), NA),
	yend = c(Ipk(k.end, lim=lim, b=b), NA),
	x = x, func = Ip, guess = 0, parms = list(lim=lim, b=b))

### Test

plot(sol)

# Plot y: perfect match
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, lim=lim, b=b))
lines(x, y, col="red", lty=2)

# Plot dy:
par(mfrow = c(1, 1))
plot(sol[, c(1,3)], type="l", col="green")
dy = sapply(x, \(k) Idy(k, lim=lim, b=b))
lines(x, dy, col="red", lty=2)


#########################

### y(x) = exp(-1/x) * y0

### ODE:
# x^4*d2y = 2*x^2*dy + (x^4-2*x-1)*y - b*x^(10/3)*exp(-1/x)

# dy  = (dy0 + 1/x^2*y0) * exp(-1/x)
# d2y = (d2y0 + 2/x^2*dy0 + (1/x^4-2/x^3)*y0)*exp(-1/x)
# d2y = 2/x^2*dy + (1-2/x^3-1/x^4)*y - b/x^(2/3)*exp(-1/x)
# where d2y0 = y0 - b/x^(2/3)


# y(k)
Ipk = function(k, b=1, lim=Inf) {
	FUN = function(k)
		integrate(\(x) sin(k*x^(3/2)) / (x^3 + 1), 0, lim,
			subdivisions = 1025)$value;
	y = try( FUN(k) );
	if(inherits(y, "try-error")) {
		cat("y(x): x = ", k, "\n");
		y = FUN(k + 0.001);
	}
	y = exp(-1/k) * y;
	# Normalization:
	y = b * y / (2/3 * gamma(2/3) * sin(pi/3));
	return(y);
}
# dy(k)
Idy = function(k, b=1, lim=Inf, up.scale = 2000) {
	# Note: various numerical issues!
	if(lim == Inf) lim = up.scale * pi;
	FUN = function(k) integrate(\(x) {
		xx = x^(2/3); # x = z^(2/3); dx = 2/3 * z^(-1/3) dz;
		xx * cos(k*x) / (x*x + 1);
	}, 0, lim, subdivisions = 43570)$value;
	dy = try( FUN(k) ); # print(str(dy))
	if(inherits(dy, "try-error")) {
		cat("x = ", k, "\n");
		dy = FUN(k + 0.001);
	}
	# Normalization:
	# dy = 2/3 * b * dy / (2/3 * gamma(2/3) * sin(pi/3));
	dy = b * dy / (gamma(2/3) * sin(pi/3));
	dy = dy * exp(-1/k) + Ipk(k, lim=lim, b=b) / k^2;
	return(dy);
}

Ip = function(x, y, pars) {
	b = pars$b;
	# y[2] = dy;
	d2y = 2/x^2*y[2] + (x^4-2*x-1)/x^4 * y[1] - b*exp(-1/x)/x^(2/3);
	list(c(y[2], d2y));
}


lim = Inf; # fixed
b = 4/3;
k.start = 0.4; k.end = 2.15;
x = seq(k.start, k.end, by = 0.01)

# bvptwp # bvpshoot(guess = 0)
sol <- bvptwp(
	yini = c(Ipk(k.start, b=b, lim=lim), NA),
	yend = c(Ipk(k.end, b=b, lim=lim), NA),
	x = x, func = Ip, parms = list(b=b, lim=lim))

### Test

plot(sol)

# Plot y: perfect match;
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(k) Ipk(k, lim=lim, b=b))
lines(x, y, col="red", lty=2)

# Plot dy:
par(mfrow = c(1, 1)); ylim = range(sol[,3]) + c(-0.075, 0.025);
plot(sol[, c(1,3)], type="l", col="green", ylim=ylim)
dy = sapply(x, \(k) Idy(k, lim=lim, b=b))
lines(x, dy, col="red", lty=2)


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


###################
###################

### Solutions

# - using Reductions / Integrating factor;
# - for more examples, see file:
#   DE.ODE.Transforms.Reductions.R;

### d2y = y + f0(x)
d2y - y - f0 # = 0
d2y + dy - dy - y - f0 # = 0
exp(x)*(d2y + dy - dy - y) - f0*exp(x) # = 0
# =>
# D(exp(x)*(dy - y)) = f0*exp(x);
# exp(x) * (dy - y)  = I( f0 * exp(x) );
# exp(-x) * (dy - y) = exp(-2*x) * I( f0 * exp(x) );
# D(exp(-x)*y) = exp(-2*x) * I( f0 * exp(x) );
# y = exp(x) * I( exp(-2*x) * I( f0 * exp(x) ) );

# TODO:
# - simplification using Integration by parts;



### d2y = k^2 * y + f0(x)
d2y - k^2 * y - f0 # = 0
d2y + k*dy - k*dy - k^2*y - f0 # = 0
exp(k*x)*(d2y + k*dy - k*(dy + k*y)) - f0*exp(k*x) # = 0
# =>
# D(exp(k*x)*(dy - k*y)) = f0*exp(k*x);
# exp(k*x) * (dy - k*y)  = I( f0 * exp(k*x) );
# D(exp(-k*x)*y) = exp(-2*k*x) * I( f0 * exp(k*x) );
# y = exp(k*x) * I( exp(-2*k*x) * I( f0 * exp(k*x) ) );

