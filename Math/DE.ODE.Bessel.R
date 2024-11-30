########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs: Bessel Functions
##
## draft v.0.1a


### Bessel Functions


####################

### Helper Functions


dBesselJ = function(x, n, eps = 1E-5) {
	(besselJ(x+eps, n) - besselJ(x, n)) / eps;
}
d2BesselJ = function(x, n, eps = 1E-5) {
	(besselJ(x+2*eps, n) + besselJ(x, n) - 2*besselJ(x+eps, n)) / eps^2;
}
d2BesselJ.all = function(x, n, eps = 1E-5) {
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
n = sqrt(2); x = sqrt(3); eps = 1E-5;
tmp = d2BesselJ.all(x, n=n, eps=eps);
y = tmp$y; dy = tmp$dy; d2y = tmp$d2y;

### ODE:
x^2*d2y + x*dy + (x^2 - n^2)*y # = 0


##################
### Transforms ###

### EXP

### Ex 1:
n = sqrt(2); x = sqrt(3);
tmp = d2BesselJ.all(x, n=n, eps=eps);
y  = exp(x) * tmp$y;
dy = exp(x) * tmp$dy + y;
d2y = exp(x) * (tmp$d2y + tmp$dy) + dy;

### ODE:
x^2*d2y - (2*x^2 - x)*dy + (2*x^2 - x - n^2)*y # = 0


#########
### Ex 2:
n = sqrt(2); x = sqrt(3);
tmp = d2BesselJ.all(x, n=n, eps=eps);
y  = exp(2*x) * tmp$y;
dy = exp(2*x) * tmp$dy + 2*y;
d2y = exp(2*x) * (tmp$d2y + 2*tmp$dy) + 2*dy;

### ODE:
x^2*d2y - (4*x^2 - x)*dy + (5*x^2 - 2*x - n^2)*y # = 0


#########
### Ex 3:
n = sqrt(2); k = 1/3;
x = sqrt(3);
tmp = d2BesselJ.all(x, n=n, eps=eps);
y  = exp(k*x) * tmp$y;
dy = exp(k*x) * tmp$dy + k*y;
d2y = exp(k*x) * (tmp$d2y + k*tmp$dy) + k*dy;

### ODE:
x^2*d2y - (2*k*x^2 - x)*dy + ((k^2+1)*x^2 - k*x - n^2)*y # = 0

