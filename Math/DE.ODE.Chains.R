

###################
### ODE: Chains ###
###################

### Helper Functions

library(bvpSolve)


#########################

### Fundamental Solution:
### d2y = y

### Mult: dy =>
# dy * d2y = y * dy
# 1/2 * (dy)^2 = 1/2 * y^2 + 1/2 * ct
# (dy)^2 = y^2 + ct;

### Classic Solutions:
# y = b1*exp(x) + b2*exp(-x);
# where 4*b1*b2 = - ct;

# TODO:
# - Are these solutions unique?
# - Does the uniqueness generalize to higher orders?


##############

##############
### Chains ###
##############

### Exp

# y = y0 + b/2 * x * exp(x)
# y0 satisfies the fundamental solution
# =>
# d2y = d2y0 + b/2 * x * exp(x) + b * exp(x)
# d2y = y0 + b/2 * x * exp(x) + b * exp(x)
# =>
# ODE: d2y = y + b * exp(x);

# TODO: check & expand;


#############

### d2y = y + sum(b[i] * exp(k[i] * x))

# y = y0 + b1/(k1^2 - 1) * exp(k1*x) + ...
# where k[i] != 1;
# y0 satisfies the fundamental solution;

# =>
# d2y = d2y0 + b1*k1^2/(k1^2 - 1) * exp(k1*x) + ...
# d2y = y0 + b1/(k1^2 - 1) * exp(k1*x) + b1*exp(k1*x) + ...
# ODE: d2y = y + b1*exp(k1*x) + ...

### Examples

# F.S.: trivial y0 = 0;

Iyk = function(x, k, b) sum(b / (k^2 - 1) * exp(k*x));

Idny = function(x, y, parms) {
	b = parms$b; k = parms$k;
	d2y = y[1] + sum(b * exp(k*x));
	list(c(y[2], d2y));
}


###
b = c(2, -3, -1/5, 1/7)
k = c(1/2,1/3, 2, 2.5)
x.start = 0; x.end = 1.5;
x = seq(x.start, x.end, by = 0.01)

sol <- bvpshoot(
	yini = c(Iyk(x.start, k=k, b=b), NA),
	yend = c(Iyk(x.end,   k=k, b=b), NA),
	x = x, func = Idny, guess = 0, parms = list(k=k, b=b))

### Test

plot(sol)

# perfect match
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(x) Iyk(x, k=k, b=b))
lines(x, y, col="red", lty=2)


#############

### d2y = y + be * x^2 * exp(ke * x)

# y = y0 + be/(ke^2 - 1) * x^2 * exp(ke*x) +
#     - 4*be*ke / (ke^2 - 1)^2 * x * exp(ke*x) +
#     - 2*be*(1 / (ke^2 - 1)^2 - 4*ke^2 / (ke^2 - 1)^3) * exp(ke*x);
# where ke[i] != 1;
# y0 satisfies the fundamental solution;

# =>
# d2y # =
d2y0 + be * x^2 * exp(ke*x) + be/(ke^2 - 1) * x^2 * exp(ke*x) +
	- 4*be*ke / (ke^2 - 1)^2 * x * exp(ke*x) +
	- 2*be*(1 / (ke^2 - 1) - 4*ke^2 / (ke^2 - 1)^2) * exp(ke*x);
# d2y = y + be * x^2 * exp(ke*x);


### Examples

# F.S.: trivial y0 = 0;
# - includes also previous example;

Iyk = function(x, ke, be, k, b) {
	div = (ke^2 - 1);
	be/div * x^2 * exp(ke*x) +
		- 4*be*ke / div^2 * x * exp(ke*x) +
		- 2*be*(1 / div^2 - 4*ke^2 / div^3) * exp(ke*x) +
	+ sum(b / (k^2 - 1) * exp(k*x));
}

Idny = function(x, y, parms) {
	b = parms$b; k = parms$k;
	be = parms$be; ke = parms$ke;
	d2y = y[1] + be * x^2 * exp(ke*x) + sum(b * exp(k*x));
	list(c(y[2], d2y));
}


###
b = c(2, -3, -1/5, 1/7)
k = c(1/2,1/3, 2, 2.5)
ke = -1/2; be = 3/2;
x.start = 0; x.end = 2.5;
x = seq(x.start, x.end, by = 0.01)

sol <- bvpshoot(
	yini = c(Iyk(x.start, ke=ke, be=be, k=k, b=b), NA),
	yend = c(Iyk(x.end,   ke=ke, be=be, k=k, b=b), NA),
	x = x, func = Idny, guess = 0,
	parms = list(ke=ke, be=be, k=k, b=b))

### Test

plot(sol)

# perfect match
par(mfrow = c(1, 1))
plot(sol[, 1:2], type="l", col="green")
y = sapply(x, \(x) Iyk(x, ke=ke, be=be, k=k, b=b))
lines(x, y, col="red", lty=2)


# TODO: chains with x^p;

