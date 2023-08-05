

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

