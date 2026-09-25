########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs - Trig: Other Variants
##
## draft v.0.1b


### Examples:

# x^2*d2y + x*dy + k^2*y = 0;
# x^2*d2y + k*y = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### Aka Complex Powers

### y = sin(k*log(x))
# - Simple example;

# Check:
k = sqrt(3); # k = 1;
x = 3^(3/5);
params = list(x=x, k=k);
e = expression(sin(k*log(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*d2y + x*dy + k^2*y # = 0

# D =>
x*dy - k*cos(k*log(x)) # = 0

# D2 =>
x*d2y + dy + k^2/x*sin(k*log(x)) # = 0


###########################

### y = x^p * sin(k*log(x))

# Check:
k = sqrt(3); # k = 1;
p = sqrt(2/5); # p = 1/2;
x = 3^(3/5);
params = list(x=x, k=k, p=p);
e = expression(x^p * sin(k*log(x)))[[1]];
# e = expression(x^p * sin(sqrt(k-1/4)*log(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
x^2*d2y - (2*p-1)*x*dy + (p^2+k^2)*y # = 0

# D =>
x*dy - p*y - k*x^p * cos(k*log(x)) # = 0

# D2 =>
x*d2y - (p-1)*dy +
	- p*k*x^(p-1) * cos(k*log(x)) +
	+ k^2*x^(p-1) * sin(k*log(x)) # = 0

### Special Cases:

### Case: p = 1/2;
x^2*d2y + (k^2+1/4)*y # = 0
# k => sqrt(k - 1/4);
x^2*d2y + k*y # = 0


###########################

### y = x^p * (sin(k*log(x)) + c1*cos(k*log(x)))
# Same Polynomial coefficient for Sin & Cos;
# Note: only c1 differs, but c1 has NO impact on ODE;

# Check:
x = 3^(3/5); c1 = 1/sqrt(2);
k = sqrt(3); # k = 1;
p = sqrt(2/5); # p = 1/2;
params = list(x=x, k=k, p=p);
e = expression(x^p * (sin(k*log(x)) + c1*cos(k*log(x))))[[1]];
# e = expression(x^p * sin(sqrt(k-1/4)*log(x)))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - p*y - k*x^p * (cos(k*log(x)) - c1*sin(k*log(x))) # = 0

# D2 =>
x^2*d2y - (p-1)*x*dy +
	+ k^2*x^p * (sin(k*log(x)) + c1*cos(k*log(x))) +
	- p*(x*dy - p*y) # = 0

### ODE:
x^2*d2y - (2*p-1)*x*dy + (p^2+k^2)*y # = 0

