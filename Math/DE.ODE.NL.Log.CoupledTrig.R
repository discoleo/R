########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs - Log w. Coupled Trig
##
## draft v.0.1a


### Examples:

# y*d2y - k*dy + 4*k^2 * y^2 + k^2 = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y = sin(k*x)*cos(k*x) * log(tan(k*x))
# - Simple example;

# Check:
k = sqrt(3); # k = 1;
x = 3^(3/5);
params = list(x=x, k=k);
e = expression(sin(k*x) * cos(k*x) * log(tan(k*x)))[[1]];
# e = expression(1/(tan(k*x) + 1/tan(k*x)) * log(tan(k*x)))[[1]];
# y = sin(k*x) * cos(k*x) * log(tan(k*x));
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
y*d2y - k*dy + 4*k^2 * y^2 + k^2 # = 0

# D =>
dy - k*(cos(k*x)^2 - sin(k*x)^2) * log(tan(k*x)) - k;

# D2 =>
d2y + 4*k^2*sin(k*x)*cos(k*x) * log(tan(k*x)) +
	- k^2*(cos(k*x)^2 - sin(k*x)^2) / (cos(k*x)*sin(k*x)) # = 0
d2y + 4*k^2 * y +
	- k^2*(cos(k*x)^2 - sin(k*x)^2) / (cos(k*x)*sin(k*x)) # = 0
d2y + 4*k^2 * y - k*(dy - k) / y # = 0

