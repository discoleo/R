#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Mixed Exp & Log
##
## draft v.0.1b


### Linear ODE of Order 2
#   with Polynomial coefficients;


### Derived from:
# y = G(x) * Log(P1(x)) * Exp(P2(x)) + F0(x)

# G, P1, P2 = Polynomials (or Polynomial Fractions);


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#####################

#####################
### Mixed Exp-Log ###

### Product:
### y = log(P1(x)) * exp(P2(x)) + F(x);

### Example 1:
### y = log(x^2 + b0) * exp(x^2)

### Check:
ye = expression(log(x^2 + b0) * exp(x^2))[[1]]
x  = sqrt(3); b0 = sqrt(2);
params = list(x=x, b0=b0);
#
y = eval(ye, params); dy = eval(D(ye, "x"), params);
d2y = eval(D(D(ye, "x"), "x"), params);


### ODE:
x*(x^2 + b0)*d2y - (4*x^4 + (4*b0 - 1)*x^2 + b0)*dy +
	+ 4*x^3*(x^2 + (b0 - 1))*y # = 0


### Derivation:

### D =>
dy - 2*x*y - 2*x*exp(x^2)/(x^2 + b0) # = 0
(x^2 + b0)*dy - 2*x*(x^2 + b0)*y - 2*x*exp(x^2) # = 0

### D2 =>
(x^2 + b0)*d2y - 2*x*(x^2 + b0 - 1)*dy - (6*x^2 + 2*b0)*y +
	- 2*(2*x^2 + 1)*exp(x^2) # = 0
x*(x^2 + b0)*d2y - 2*x^2*(x^2 + b0 - 1)*dy - x*(6*x^2 + 2*b0)*y +
	- (2*x^2 + 1)*((x^2 + b0)*dy - 2*x*(x^2 + b0)*y) # = 0


### Solution & Plot:
y = function(x, b=1) {
	x2 = x^2;
	val = log(x2 + b[1]) * exp(x2);
	return(val)
}
dy = function(x, b=1) {
	yx = y(x, b=b);
	x2 = x^2; x2b = x2 + b[1];
	dyx = 2*x*(x2b*yx + exp(x2));
	div = x2b;
	dyx = ifelse(div != 0, dyx/div, 1); # TODO: check;
	return(dyx)
}
d2y = function(x, b=1) {
	yx  =  y(x, b=b);
	dyx = dy(x, b=b);
	#
	x2 = x^2; x2b = x2 + b[1];
	d2p = (4*x2*x2b - x2 + b[1])*dyx - 4*x2*x*(x2b - 1)*yx;
	div = x*x2b;
	d2p = ifelse(div != 0, d2p/div, 1); # TODO: check;
	return(d2p);
}
### Plot:
b = 1
px = (1:3)*2/7 + 1; px = c(0, -px, px);
curve(y(x, b=b), from = -2, to = 2, ylim=c(-30, 80))
line.tan(px, dx=3, p=y, dp=dy, b=b)
# global minimum:
curve(dy(x, b=b), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, b=b, col="orange")

### Ex 2:
b = -1
px = (1:4)*2/7 + 1;
curve(y(x, b=b), from = 1, to = 2.5, ylim=c(-10, 200))
line.tan(px, dx=3, p=y, dp=dy, b=b)
# global minimum:
curve(dy(x, b=b), add=T, col="green")
line.tan(px, dx=3, p=dy, dp=d2y, b=b, col="orange")


###################
###################

### y = log(x + b0) * exp(k/x)

### Check:
ye = expression(log(x + b0) * exp(k/x))[[1]]
x  = sqrt(3); b0 = sqrt(2); k = -1/sqrt(5);
params = list(x=x, k=k, b0=b0);
#
y = eval(ye, params); dy = eval(D(ye, "x"), params);
d2y = eval(D(D(ye, "x"), "x"), params);

### ODE:
x^4*(x + b0)*d2y + x^2*(x^2 + 2*k*x + 2*k*b0)*dy +
	- k*(x^2 + (2*b0 - k)*x - k*b0)*y # = 0


### Derivation:

# D =>
x^2*(x + b0)*dy + k*(x + b0)*y - x^2*exp(k/x) # = 0

# D2 =>
x^2*(x + b0)*d2y + (3*x^2 + 2*b0*x + k*x + k*b0)*dy + k*y +
	- (2*x - k)*exp(k/x) # = 0
x^4*(x + b0)*d2y + x^2*(3*x^2 + 2*b0*x + k*x + k*b0)*dy + k*x^2*y +
	- (2*x - k)*(x^2*(x + b0)*dy + k*(x + b0)*y) # = 0
x^4*(x + b0)*d2y + x^2*(x^2 + 2*k*x + 2*k*b0)*dy +
	- k*(x^2 + (2*b0 - k)*x - k*b0)*y # = 0


### Special Cases:

### k = 2*b0
b0 = 1/sqrt(5); k = 2*b0;
params = list(x=x, k=k, b0=b0);
#
x^4*(x + b0)*d2y + x^2*(x + 2*b0)^2*dy - 2*b0*(x^2 - 2*b0^2)*y # = 0

