########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log w. Coupled Trig
##
## draft v.0.1b


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

#########################

### y = x^p * sin(k*x)*cos(k*x) * log(tan(k*x))
# - Simple example;

# Check:
k = sqrt(3); # k = 1;
p = 1/2^(2/5); # p = -1;
x = 3^(3/5);
params = list(x=x, k=k, p=p);
e = expression(x^p * sin(k*x) * cos(k*x) * log(tan(k*x)))[[1]];
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - p*y - k*x^(p+1) * (cos(k*x)^2 - sin(k*x)^2) * log(tan(k*x)) - k*x^(p+1) # = 0

# D2 =>
x*d2y - (p-1)*dy +
	- (p+1)*(x*dy - p*y - k*x^(p+1)) / x +
	+ 4*k^2*x^(p+1) * sin(k*x)*cos(k*x) * log(tan(k*x)) +
	- k^2*x^(p+1) * (cos(k*x)^2 - sin(k*x)^2) / (sin(k*x)*cos(k*x)) +
	- k*(p+1)*x^p # = 0
x*d2y - 2*p*dy +
	+ 4*k^2*x * y + p*(p+1)*y / x +
	- k*x^p * (x*dy - p*y - k*x^(p+1)) / y # = 0

### ODE:
x^2*y*d2y - 2*p*x*y*dy - k*x^(p+2) * dy +
	+ (4*k^2*x^2 + p*(p+1)) * y^2 +
	+ k*p*x^(p+1) * y + k^2*x^(2*p+2) # = 0

### Special Cases:
# - For case p = 0, see above;

### Case: p = -1;
x^2 * y*d2y + 2*x*y*dy - k*x*dy + 4*k^2*x^2 * y^2 - k*y + k^2 # = 0

