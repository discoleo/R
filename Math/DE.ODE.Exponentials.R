
########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Exponentials
###
### draft v.0.1d


### ODEs Derived from Exponentials

### TODO:
# - move relevant sections from
#   DE.ODE.Fractions.Lambert.R
#   & DE.ODE.Gaussian.R to this file;



#########################

### Helper functions

library(pracma)
# - may be needed to solve various equations;
#   (e.g. Lambert W)


# include: DE.ODE.Helper.R;
source("DE.ODE.Helper.R")

#########################
#########################


### Section B: Simple
### Non-Linear ODEs

##############
### e^y = P(x)

### D =>
e^y * dy = dp
### D2 =>
e^y * d2y + e^y * dy^2 - d2p # = 0

### ODE:
p*d2y + p*dy^2 - d2p # = 0
# - can be converted trivially to order 1 (but still non-liniar);

### Examples:
### P(x) = x^2
x^2*d2y + x^2*dy^2 - 2 # = 0
### P(x) = sin(a*x)
d2y + dy^2 + a^2 # = 0
### P(x) = cos(x)^3
d2y + dy^2 - 6*tan(x)^2 + 3 # = 0


##################

### Section C: High-Power
### Non-Linear ODEs

##################
### e^(y^2) = P(x)

### D =>
2*e^(y^2) * y*dy = dp
### D2 =>
2*e^(y^2) * (y*d2y + dy^2) + 4*e^(y^2) * y^2*dy^2 - d2p # = 0
2*p*y*d2y + 4*p*y^2*dy^2 + 2*p*dy^2 - d2p # = 0

### ODE:
2*p*y*d2y + 4*p*y^2*dy^2 + 2*p*dy^2 - d2p # = 0

### Examples:
### P(x) = x
y*d2y + 2*y^2*dy^2 + dy^2 # = 0
### P(x) = x^2
x^2*y*d2y + 2*x^2*y^2*dy^2 + x^2*dy^2 - 1 # = 0
### P(x) = sqrt(x^3)
8*x^2*y*d2y + 16*x^2*y^2*dy^2 + 8*x^2*dy^2 - 3 # = 0
### P(x) = sqrt(x^5)
8*x^2*y*d2y + 16*x^2*y^2*dy^2 + 8*x^2*dy^2 - 15 # = 0
### P(x) = sin(a*x)
2*y*d2y + 4*y^2*dy^2 + 2*dy^2 + a^2 # = 0
### P(x) = sin(sqrt(3)/2 * x)
8*y*d2y + 16*y^2*dy^2 + 8*dy^2 + 3 # = 0


### Solution & Plot:
y.gen = function(n) {
	n1=n;
	y1.f = function(x, n=n1) { x^n }
	y1.d1.df = function(x, n=n1) { n*x^(n-1) }
	y1.d2.df = function(x, n=n1) { if(n == 2) 2  else n*(n-1)*x^(n-2) }
	list(y1.f, y1.d1.df, y1.d2.df)
}
y1.lst = y.gen(2)
y = function(x, PFUN, posRoot=TRUE) {
	val = sqrt(log(PFUN[[1]](x)));
	if( ! posRoot) val = - val;
	return(val)
}
dy = function(x, PFUN, posRoot=TRUE) {
	y.x = y(x, PFUN=PFUN, posRoot=posRoot)
	dp = PFUN[[2]](x) / y.x / exp(y.x^2) / 2;
	return(dp)
}
d2y = function(x, PFUN, posRoot=TRUE) {
	p.v = PFUN[[1]](x);
	y.x  = y(x, PFUN=PFUN, posRoot=posRoot);
	dy.x = dy(x, PFUN=PFUN, posRoot=posRoot);
	dp = PFUN[[3]](x) - 4*p.v*y.x^2*dy.x^2 - 2*p.v*dy.x^2
	div = 2 * p.v * y.x;
	dp = ifelse(div != 0, dp/ div, 1) # TODO
	return(dp)
}
### Plot:
curve(y(x, PFUN=y1.lst), from= 1, to = 3)
sapply(c(3/5 + (1:3)*3/5), line.tan, dx=3, p=y, dp=dy, PFUN=y1.lst)
# inverse-exp-like:
curve(dy(x, PFUN=y1.lst), add=T, col="green")
sapply(c(3/5 + (1:3)*3/5), line.tan, dx=3, p=dy, dp=d2y, PFUN=y1.lst, col="orange")


### Ex 2:
n = 3/2;
y1.lst = y.gen(n=n)
curve(y(x, PFUN=y1.lst), from= 1, to = 3)
sapply(c(3/5 + (1:3)*3/5), line.tan, dx=3, p=y, dp=dy, PFUN=y1.lst)
# inverse-exp-like:
curve(dy(x, PFUN=y1.lst), add=T, col="green")
sapply(c(3/5 + (1:3)*3/5), line.tan, dx=3, p=dy, dp=d2y, PFUN=y1.lst, col="orange")


########################
########################

########################
### Section D:
### Non-Linear ODEs
### Product-Type
########################


########################
### Integral-Product ###
########################

### I(e^y) dx * I(e^(-y)) dx = P(x)

# - trivial: y^1;
#   [but interesting approach to some integrals]

### D =>
e^y * In + e^(-y) * Ip - dp # = 0
### D2 =>
e^y * dy * In - e^(-y) * dy * Ip - d2p + 2 # = 0
### Solve liniar system =>
e^y * In = (dp*dy + d2p - 2) / (2*dy)
e^(-y) * Ip = (dp*dy - d2p + 2) / (2*dy)

### ODE:
(dp*dy + d2p - 2)*(dp*dy - d2p + 2) - 4*p*dy^2 # = 0
(dp^2 - 4*p)*dy^2 - (d2p - 2)^2 # = 0

### Examples:
### P(x) = x
(1 - 4*x)*dy^2 - 4 # = 0
# y = sqrt(1 - 4*x)

### Plot
# TODO: proper sign pairing;
y = function(x, posRoot=TRUE) {
	val = sqrt(1 - 4*x)
	if( ! posRoot) val = - val;
	return(val)
}
y.exp = function(x, posInt=TRUE) {
	val = sqrt(1 - 4*x)
	if( ! posInt) val = - val;
	val = exp(val);
	return(val)
}
dy = function(x, posRoot=TRUE) {
	dp = 2 / sqrt(1 - 4*x)
	if( ! posRoot) dp = - dp;
	return(dp)
}
Ip.f = function(x) {
	dy.x = dy(x)
	(dy.x + 2) * y.exp(x) / (2*dy.x)
}
Ip.int = function(x, upper=0, diff=3/4) {
	val = sapply(x - diff,
		function(lower) {
		integrate(y.exp, lower=lower, upper=upper)$value
		})
}
y.int = function(lower, upper=0) {
	val = sapply(lower,
		function(lower) {
		integrate(y.exp, lower=lower, upper=upper)$value *
		integrate(y.exp, lower=lower, upper=upper, posInt=FALSE)$value # - upper
	} )
}
###
lim = c(-7, 1/4)
# interesting integral
# D ( 1/2*(1-sqrt(1-4*x))*exp(sqrt(1-4*x))) = exp(sqrt(1-4*x));
# Note: switched signs; [TODO]
curve(Ip.f(x), from=lim[1], to=lim[2])
curve(Ip.int(x, upper=1/4, diff=3/4), add=T, col="green")
# TODO:
# curve(y.int(x), from=lim[1], to=lim[2])



########################


###################
### Log-Product ###
###################

### log(y + x) * log(y - x) = P(x)

### D =>
(dy + 1)*(y - x)*log(y - x) + (dy - 1)*(y + x)*log(y + x) - (y^2 - x^2)*dp # = 0

### D2 =>
(d2y*(y - x) + dy^2 - 1)*log(y - x) + (d2y*(y + x) + dy^2 - 1)*log(y + x) +
	+ 2*dy^2 - 2 - 2*(y*dy - 2*x)*dp - (y^2 - x^2)*d2p # = 0

### Solve liniar system:
log(y + x) = - ... /2 / (y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)
log(y - x) = ... /2 / (y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)


### ODE:
# TODO

### Example:
### P(x) = x
(dy + 1)*(y - x)*log(y - x) + (dy - 1)*(y + x)*log(y + x) - (y^2 - x^2) # = 0
(d2y*(y - x) + dy^2 - 1)*log(y - x) + (d2y*(y + x) + dy^2 - 1)*log(y + x) +
	+ 2*(dy^2 - y*dy + x - 1) # = 0
# =>
2*log(y + x) = - (y - x) * (2*(dy + 1)*(dy^2 - y*dy + x - 1) + (d2y*(y - x) + dy^2 - 1)*(y + x)) /
	(y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)
2*log(y - x) = - (y + x) * (2*(dy - 1)*(dy^2 - y*dy + x - 1) + (d2y*(y + x) + dy^2 - 1)*(y - x)) /
	(y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)
# =>
2*log(y + x) = - (y - x)*(y^2*d2y - x^2*d2y + 2*dy^3 - y*dy^2 + x*dy^2 + 2*dy^2 - 2*y*dy + 2*x*dy - 2*dy - y + x - 2) /
	(y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)
2*log(y - x) = - (y + x)*(...) /
	(y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)

### ODE:
(y^2 - x^2)*(...)*(...) +
	- 4*x*(y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)^2 # = 0

### TODO: check!


#########################
#########################

### I(f1(y) dy) = f2(y, x)

### I( exp(y) / y dy) = P1(x)*P2(y)

### Example:
# I( exp(y) / y dy) = x*y

### D =>
exp(y) / y - x*dy - y # = 0
exp(y) - x*y*dy - y^2 # = 0

### D2 =>
exp(y)*dy - x*y*d2y - x*dy^2 - y*dy - 2*y*dy # = 0
(x*y*dy + y^2)*dy - x*y*d2y - x*dy^2 - 3*y*dy # = 0
x*y*d2y - x*y*dy^2 + x*dy^2 - y^2*dy + 3*y*dy # = 0

# TODO: check!

