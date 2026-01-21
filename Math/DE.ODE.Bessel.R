########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs: Bessel Functions
##
## draft v.0.1d


### Bessel Functions

# - see e.g. Wikipedia:
#   https://en.wikipedia.org/wiki/Bessel_function


####################

### Helper Functions


dBesselJ = function(x, n) {
	# (besselJ(x+eps, n) - besselJ(x, n)) / eps;
	id  = 0:49; id1 = id+1;
	idn = 2*id + n;
	sg = rep(c(1,-1), 25);
	sum(sg * idn / gamma(id1) / gamma(id1+n) * (x/2)^(idn-1)) / 2;
}
d2BesselJ = function(x, n) {
	# (besselJ(x+2*eps, n) + besselJ(x, n) - 2*besselJ(x+eps, n)) / eps^2;
	id  = 0:49; id1 = id+1;
	sg = rep(c(1,-1), 25);
	idn = 2*id + n;
	sum(sg * idn*(idn-1) / gamma(id1) / gamma(id1+n) * (x/2)^(idn-2)) / 4;
}
d2BesselJ.all = function(x, n, eps = 1E-5) {
	# eps: is not used;
	y  = besselJ(x, n);
	dy = dBesselJ(x, n);
	d2y = d2BesselJ(x, n);
	list(y=y, dy=dy, d2y=d2y);
}
d2BesselJ.all.old = function(x, n, eps = 1E-5) {
	y  = besselJ(x, n);
	y1 = besselJ(x+eps, n);
	y2 = besselJ(x+2*eps, n);
	dy  = (y1 - y) / eps;
	d2y = (y2 + y - 2*y1) / eps^2;
	list(y=y, dy=dy, d2y=d2y);
}


####################
####################

### Basic:
n = sqrt(2); x = sqrt(3);
tmp = d2BesselJ.all(x, n=n);
y = tmp$y; dy = tmp$dy; d2y = tmp$d2y;

### ODE:
x^2*d2y + x*dy + (x^2 - n^2)*y # = 0


##################

##################
### Transforms ###

### EXP

### Ex 1:
n = sqrt(2); x = sqrt(3);
tmp = d2BesselJ.all(x, n=n);
y  = exp(x) * tmp$y;
dy = exp(x) * tmp$dy + y;
d2y = exp(x) * (tmp$d2y + tmp$dy) + dy;

### ODE:
x^2*d2y - (2*x^2 - x)*dy + (2*x^2 - x - n^2)*y # = 0


#########
### Ex 2:
n = sqrt(2); x = sqrt(3);
tmp = d2BesselJ.all(x, n=n);
y  = exp(2*x) * tmp$y;
dy = exp(2*x) * tmp$dy + 2*y;
d2y = exp(2*x) * (tmp$d2y + 2*tmp$dy) + 2*dy;

### ODE:
x^2*d2y - (4*x^2 - x)*dy + (5*x^2 - 2*x - n^2)*y # = 0


#########
### Ex 3:
n = sqrt(2); k = 1/3;
x = sqrt(3);
tmp = d2BesselJ.all(x, n=n);
y  = exp(k*x) * tmp$y;
dy = exp(k*x) * tmp$dy + k*y;
d2y = exp(k*x) * (tmp$d2y + k*tmp$dy) + k*dy;

### ODE:
x^2*d2y - (2*k*x^2 - x)*dy + ((k^2+1)*x^2 - k*x - n^2)*y # = 0



### EXP(x^n)

### Ex 1: Exp(x^2) * Bessel(x)
n = sqrt(2); x = sqrt(3);
tmp = d2BesselJ.all(x, n=n);
y  = exp(x^2) * tmp$y;
dy = exp(x^2) * tmp$dy + 2*x*y;
d2y = exp(x^2) * (tmp$d2y + 2*x*tmp$dy) + 2*x*dy + 2*y;

### ODE:
x^2*d2y - x*(4*x^2 - 1)*dy + (4*x^4 - 3*x^2 - n^2)*y # = 0


#########
### Ex 2: Exp(2*x^2) * Bessel(x)
n = sqrt(2); x = sqrt(3);
tmp = d2BesselJ.all(x, n=n);
y  = exp(2*x^2) * tmp$y;
dy = exp(2*x^2) * tmp$dy + 4*x*y;
d2y = exp(2*x^2) * (tmp$d2y + 4*x*tmp$dy) + 4*x*dy + 4*y;

### ODE:
x^2*d2y - x*(8*x^2 - 1)*dy + (16*x^4 - 7*x^2 - n^2)*y # = 0


#########
### Ex 3: Exp(k*x^2) * Bessel(x)
n = sqrt(2); x = sqrt(3); k = 1/3;
tmp = d2BesselJ.all(x, n=n);
y  = exp(k*x^2) * tmp$y;
dy = exp(k*x^2) * tmp$dy + 2*k*x*y;
d2y = exp(k*x^2) * (tmp$d2y + 2*k*x*tmp$dy) + 2*k*x*dy + 2*k*y;

### ODE:
x^2*d2y - x*(4*k*x^2 - 1)*dy + (4*k^2*x^4 - (4*k-1)*x^2 - n^2)*y # = 0


# Special Case: k = 1/4
n = sqrt(2); x = sqrt(3);
tmp = d2BesselJ.all(x, n=n);
y  = exp(x^2 / 4) * tmp$y;
dy = exp(x^2 / 4) * tmp$dy + x*y/2;
d2y = exp(x^2 / 4) * (tmp$d2y + 1/2*x*tmp$dy) + x*dy/2 + y/2;

### ODE:
4*x^2*d2y - 4*x*(x^2 - 1)*dy + (x^4 - 4*n^2)*y # = 0


##################
##################

### Bessel-Like

# Exploration of Bessel-Like Functions

### y = x^p * sin(k*x^n)
# Check:
x = sqrt(3); k = 1/5; p = sqrt(2); n = sqrt(3); c0 = -1/3;
# k = 1/n; p = - n/2;
params = list(x=x, p=p, n=n, c0=c0);
e = expression(x^p * sin(k*x^n) + c0)[[1]];
#
y   = eval(e, params); dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*d2y - (2*p+n-1)*x*dy +
	+ (k^2*n^2*x^(2*n) + p*(p+n)) * (y - c0) # = 0


# D =>
x*dy - p*(y - c0) - k*n*x^(p+n) * cos(k*x^n) # = 0

# D2 =>
x^2*d2y - (p-1)*x*dy +
	- k*n*(p+n)*x^(p+n) * cos(k*x^n) +
	+ k^2*n^2*x^(p+2*n) * sin(k*x^n) # = 0
x^2*d2y - (2*p+n-1)*x*dy +
	+ (k^2*n^2*x^(2*n) + p*(p+n)) * (y - c0) # = 0


### Special Cases:

### Case: p = - n/2; k = 1/n;
m = n/2;
x^2*d2y + x*dy + (x^(4*m) - m^2) * (y - c0) # = 0


###########
### y = x^p * sin(k*log(x))
# Check:
x = sqrt(3); k = 1/5; p = sqrt(2); n = sqrt(3); c0 = -1/3;
# k = 1/n; p = - n/2;
params = list(x=x, k=k, p=p, c0=c0);
e = expression(x^p * sin(k*log(x)) + c0)[[1]];
#
y   = eval(e, params); dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*d2y - (2*p-1)*x*dy + (k^2 + p^2) * (y - c0) # = 0


# D =>
x*dy - p*(y - c0) - k*x^p*cos(k*log(x)) # = 0

# D2 =>
x^2*d2y - (p-1)*x*dy - k*p*x^p*cos(k*log(x)) +
	+ k^2*x^p * sin(k*log(x)) # = 0
x^2*d2y - (2*p-1)*x*dy + (k^2 + p^2) * (y - c0) # = 0


##############
### y = exp(k1*x) * sin(k2*log(x))

# Check:
x = sqrt(3); k1 = sqrt(2/3); k2 = 1/5; c0 = -1/3;
params = list(x=x, k1=k1, k2=k2, c0=c0);
e = expression(exp(k1*x) * sin(k2*log(x)) + c0)[[1]];
#
y   = eval(e, params); dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*d2y - x*(2*k1*x-1)*dy + (k1^2*x^2 - k1*x + k2^2) * (y - c0) # = 0


# D =>
x*dy - k1*x*(y - c0) - k2*exp(k1*x)*cos(k2*log(x)) # = 0

# D2 =>
x^2*d2y - x*(k1*x-1)*dy - k1*x*(y - c0) +
	- k1*k2*x*exp(k1*x)*cos(k2*log(x)) +
	+ k2^2*exp(k1*x)*sin(k2*log(x)) # = 0
x^2*d2y - x*(2*k1*x-1)*dy + (k1^2*x^2 - k1*x + k2^2) * (y - c0) # = 0


###############
### y = sqrt(x) * exp(k*sqrt(x))

# Check:
x = sqrt(3); k = 1/5; c0 = -1/3;
params = list(x=x, k=k, c0=c0);
e = expression(sqrt(x) * exp(k*sqrt(x)) + c0)[[1]];
#
y   = eval(e, params); dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
4*x^2*d2y - 2*x*dy - (k^2*x - 2) * (y - c0) # = 0


# D =>
2*x*dy - (y - c0) - k * x * exp(k*sqrt(x)) # = 0

# D2 =>
4*x*d2y + 2*dy - 2*k * exp(k*sqrt(x)) +
	- k^2 * sqrt(x) * exp(k*sqrt(x)) # = 0
4*x^2*d2y - 2*x*dy - (k^2*x - 2) * (y - c0) # = 0

########
### Gen:
# y = x^n * sqrt(x) * exp(k*sqrt(x))

# Check:
x = sqrt(3); n = 1/3; k = 1/5; c0 = -1/3;
# n = 1/4;
params = list(x=x, k=k, c0=c0);
e = expression(x^n * sqrt(x) * exp(2*k*sqrt(x)) + c0)[[1]];
#
y   = eval(e, params); dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*d2y - (2*n+1/2)*x*dy +
	- (k^2*x - (n+1)*(n+1/2)) * (y - c0) # = 0

# D =>
2*x*dy - (2*n+1)*(y - c0) - 2*k * x^(n+1) * exp(2*k*sqrt(x)) # = 0

# D2 =>
4*x*d2y - 2*(2*n-1)*dy - 4*k*(n+1) * x^n * exp(2*k*sqrt(x)) +
	- 4*k^2 * x^n * sqrt(x) * exp(2*k*sqrt(x)) # = 0
x^2*d2y - (2*n+1/2)*x*dy +
	- (k^2*x - (n+1)*(n+1/2)) * (y - c0) # = 0


### Special Cases:

### Case: n = 1/4
x^2*d2y - x*dy - (k^2*x - 15/16) * (y - c0) # = 0

