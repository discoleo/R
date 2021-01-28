
########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Exponentials
###
### draft v.0.1a-ex2


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

