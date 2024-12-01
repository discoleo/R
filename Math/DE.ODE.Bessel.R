########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs: Bessel Functions
##
## draft v.0.1c


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

