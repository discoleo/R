########################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Transforms
## Order Reduction
##
## v.0.1e


### Order Reduction of ODEs


####################

### Helper Functions

library(bvpSolve)


####################
####################

### Integrating Factors

### Note:
# - for the direct approach using Laplace transforms, see:
#   Dr Peyam: Convolution and ODE
#   https://www.youtube.com/watch?v=OWxZ7frYiUo
# - but see below for generalization of this approach, e.g.:
#   d2y = 4*k^2*x^2*y + f0;


### Base ODE:
# d2y = y + f0

### Transform: Exp(x)
# y = z * exp(x)
# dy = (dz + z) * exp(x)
# d2y = (d2z + 2*dz + z) * exp(x)

### Order Reduction
# =>
# d2z + 2*dz = f0 * exp(-x)



### Transform: Exp(-x)
# y = z * exp(-x)
# dy = (dz - z) * exp(-x)
# d2y = (d2z - 2*dz + z) * exp(-x)

### Order Reduction
# =>
# d2z - 2*dz = f0 * exp(x)


###################
###################

### Base ODE:
# d2y = 4*k^2*x^2*y + f0


### Transform: Exp(k*x^2)
# y = z * exp(k*x^2)
# dy = (dz + 2*k*x*z) * exp(k*x^2)
# d2y = (d2z + 4*k*x*dz + 4*k^2*x^2*z) * exp(k*x^2)

### Order Reduction
# =>
# d2z + 4*k*x*dz = f0 * exp(-k*x^2)
# D(exp(2*k*x^2) * dz) = f0 * exp(k*x^2);


###################

### Gen: x^n * y

### Base ODE:
# d2y = (n*k)^2 * x^(2*n-2)*y + f0


### Transform: Exp(k*x^n)
# y = z * exp(k*x^n)
# dy = (dz + n*k*x^(n-1)*z) * exp(k*x^n)
# d2y = (d2z + 2*n*k*x^(n-1)*dz + n^2*k^2*x^(2*n-2)*z) * exp(k*x^n)

### Order Reduction
# =>
# d2z + 2*n*k*x^(n-1)*dz = f0 * exp(-k*x^n)


#################
### Composite ###

### Base ODE:
# d2z = (4*k2^2*x^2 + 4*k1*k2*x + k1^2 - 2*k2) * z + f0(x);

### Transform: Exp(k2*x^2 + k1*x)
# y  =
z * exp(k2*x^2 + k1*x);
# dy =
dz * exp(k2*x^2 + k1*x) + (2*k2*x + k1) * y;
# d2y =
d2z * exp(k2*x^2 + k1*x) +
	+ 2*(2*k2*x + k1) * dy +
	- (4*k2^2*x^2 + 4*k1*k2*x + k1^2 - 2*k2) * y;

# Order reduction:
# d2y = 2*(2*k2*x + k1) * dy + f0 * exp(k2*x^2 + k1*x);


### Example

# Base ODE:
Ip0 = function(x, y, pars) {
	ff = pars$FUN;
	k1 = pars$k1; k2 = pars$k2;
	z = y[1];
	d2y = (4*k2^2*x^2 + 4*k1*k2*x + k1^2 - 2*k2) * z + ff(x);
	list(c(y[2], d2y));
}
# Transformed ODE:
Ip = function(x, y, pars) {
	ff = pars$FUN;
	k1 = pars$k1; k2 = pars$k2;
	dy = y[2];
	d2y = 2*(2*k2*x + k1) * dy + ff(x) * exp(k2*x^2 + k1*x);
	list(c(dy, d2y));
}

###
x.start = 0.2; x.end = 2;
x = seq(x.start, x.end, by = 0.01)
y = c(1, 5)
k = c(4/3, -1/5); # k = c(7/3, -6/5);
FUN = function(x) x * exp(-x^2/3);
parms = list(FUN=FUN, k1=k[1], k2=k[2]);

# bvptwp # bvpshoot(guess = 0)
sol0 <- bvptwp(
	yini = c(y[1], NA),
	yend = c(y[2], NA),
	x = x, func = Ip0, parms = parms)
fs = exp(k[2]*x[1]^2 + k[1]*x[1]);
xe = x[length(x)]; fe = exp(k[2]*xe^2 + k[1]*xe);
sol <- bvptwp(
	yini = c(fs * y[1], NA),
	yend = c(fe * y[2], NA),
	x = x, func = Ip, parms = parms)

### Test

plot(sol)
par(mfrow = c(1, 1))

# Plot y vs y0: perfect match;
y.tmp = sol0[,2] * exp(k[2]*x^2 + k[1]*x);
plot(sol0[, 1], y.tmp, type="l", col="green")
lines(sol[, 1:2], col="red", lty=2)


########
### Gen:

### Base ODE:
# x^2 * d2z = (k1*n1*x^n1 + k2*n2*x^n2 + 1) * (k1*n1*x^n1 + k2*n2*x^n2) * z +
#	- (k1*n1^2*x^n1 + k2*n2^2*x^n2) * z + x^2 * f0(x);


### Transform: Exp(k1*x^n1 + k2*x^n2)
# y  =
z * exp(k1*x^n1 + k2*x^n2);
# x*dy =
x * dz * exp(k1*x^n1 + k2*x^n2) + (k1*n1*x^n1 + k2*n2*x^n2) * y;
# x^2*d2y =
x^2 * d2z * exp(k1*x^n1 + k2*x^n2) +
	+ 2*x * (k1*n1*x^n1 + k2*n2*x^n2) * dy +
	- (k1*n1*x^n1 + k2*n2*x^n2 + 1) * (k1*n1*x^n1 + k2*n2*x^n2) * y +
	+ (k1*n1^2*x^n1 + k2*n2^2*x^n2) * y;

### Order Reduction
# =>
d2y - 2*(k1*n1*x^n1 + k2*n2*x^n2) / x * dy +
	- f0 * exp(k1*x^n1 + k2*x^n2) # = 0


### Example

# Base ODE:
Ip0 = function(x, y, pars) {
	ff = pars$FUN;
	n1 = pars$n[1]; n2 = pars$n[2];
	k1 = pars$k[1]; k2 = pars$k[2];
	z = y[1];
	d2y = (k1*n1*x^n1 + k2*n2*x^n2 + 1) * (k1*n1*x^n1 + k2*n2*x^n2) * z +
		- (k1*n1^2*x^n1 + k2*n2^2*x^n2) * z + x^2 * ff(x);
	d2y = d2y / x^2;
	list(c(y[2], d2y));
}
# Transformed ODE:
Ip = function(x, y, pars) {
	ff = pars$FUN;
	n1 = pars$n[1]; n2 = pars$n[2];
	k1 = pars$k[1]; k2 = pars$k[2];
	dy = y[2];
	d2y = 2*(k1*n1*x^n1 + k2*n2*x^n2) / x * dy +
		+ ff(x) * exp(k1*x^n1 + k2*x^n2) # = 0
	list(c(dy, d2y));
}

###
x.start = 0.2; x.end = 3;
x = seq(x.start, x.end, by = 0.01)
y = c(1, 5)
n = c(1/2, 4/3);
k = c(4/3, -2/5); # k = c(7/3, -6/5);
FUN = function(x) x * exp(-x^2/3);
parms = list(FUN=FUN, n=n, k=k);

# bvptwp # bvpshoot(guess = 0)
sol0 <- bvptwp(
	yini = c(y[1], NA),
	yend = c(y[2], NA),
	x = x, func = Ip0, parms = parms)
fs = exp(k[2]*x[1]^n[2] + k[1]*x[1]^n[1]);
xe = x[length(x)]; fe = exp(k[2]*xe^n[2] + k[1]*xe^n[1]);
sol <- bvptwp(
	yini = c(fs * y[1], NA),
	yend = c(fe * y[2], NA),
	x = x, func = Ip, parms = parms)

### Test

plot(sol)
par(mfrow = c(1, 1))

# Plot y vs y0: perfect match;
y.tmp = sol0[,2] * exp(k[2]*x^n[2] + k[1]*x^n[1]);
plot(sol0[, 1], y.tmp, type="l", col="green")
lines(sol[, 1:2], col="red", lty=2)

