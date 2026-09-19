########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Lambert W
##
## draft v.0.1a



#########################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


# Evaluate Derivatives:
eval.d = function(ye, fe, params) {
	fe = substitute(fe, list(fe = fe));
	fe = fe[[1]]; ye = ye[[1]];
	fe = do.call(substitute, list(fe, list(y = ye)));
	#
	y   = eval(ye, params); ye = D(ye, "x");
	dy  = eval(ye, params); ye = D(ye, "x");
	d2y = eval(ye, params);
	#
	f0  = eval(fe, params); fe = D(fe, "x");
	df0 = eval(fe, params); fe = D(fe, "x");
	d2f = eval(fe, params);
	lst = list(y=y, dy=dy, d2y=d2y,
		f0=f0, df0=df0, d2f=d2f);
	return(lst);
}


#########################
#########################

### y * log(y)^2 = F0(x)

# y = exp(lambertWp(sqrt(x)/2))^2;

# Check:
ye = expression(x * log(x+1));
fe = expression(y * log(y)^2);
#
x = sqrt(3);
params = list(x = x);
z = eval.d(ye, fe, params=params);
y  = z$y;  dy  = z$dy;  d2y = z$d2y;
f0 = z$f0; df0 = z$df0; d2f = z$d2f;


# D =>
log(y)^2 * dy + 2*log(y)*dy - df0 # = 0
2*y*log(y)*dy + f0 * dy - df0*y # = 0

# D2 =>
2*y*log(y)*d2y + 2*dy^2 + 2*log(y)*dy^2 +
	+ f0 * d2y - d2f*y # = 0
2*y^2*log(y)*d2y + 2*y*dy^2 - (f0*dy - df0*y)*dy +
	+ f0 * y*d2y - d2f*y^2 # = 0
- y*(f0 * dy - df0*y)*d2y + 2*y*dy^3 - (f0*dy - df0*y)*dy^2 +
	+ f0 * y*dy*d2y - d2f*y^2*dy # = 0

###  ODE:
df0*y^2*d2y + 2*y*dy^3 - f0*dy^3 + df0*y*dy^2 - d2f*y^2*dy # = 0


### Special Cases:

# Case: f0 = x;
y^2*d2y + 2*y*dy^3 - x*dy^3 + y*dy^2 # = 0

# Case: f0 = 1/x;
x*y^2*d2y - 2*x^3 * y*dy^3 + x^2 * dy^3 + x*y*dy^2 + y^2*dy # = 0

# Check:
x = sqrt(3);
f0 = x; df0 = 1; d2f = 0;
zW = pracma::lambertWp(sqrt(x)/2);
y  = exp(2*zW);
dy = y / sqrt(x) / (2*(zW + 1) * exp(zW));
# TODO: d2y = dy^2 / y - dy / (2*x) - ...;

