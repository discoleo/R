########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Trig: Tan
##
## draft v.0.1b


### NL ODEs: Trigonometric
#   Trivial Tan


### Examples

# x*d2y - 2*K*x^n * y*dy - (n-1)*dy = 0;
# x^2*d2y - 2*K*x^(n-p+1) * y*dy - (n-1)*x*dy + p*(n-p)*y - 2*p*K*x^(n+p) = 0;
# d2y - K*x^(n-2) * (2*K*x^n * y^3 + (n-1)*y^2 + 2*K*x^n * y + n-1) = 0;


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


########################

### Simple Components

### y = tan(k*x^n)

# Check:
e = expression(tan(k*x^n));
k = sqrt(5); n = 1/sqrt(2);
x = sqrt(3);
params = list(x=x, k=k);
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
dy - k*n*x^(n-1) * (1 + y^2) # = 0

# D2 =>
d2y - k*n*(n-1)*x^(n-2) * (1 + y^2) +
	- 2*k*n*x^(n-1) * y*dy # = 0

### Variant 1:
d2y - k*n*(n-1)*x^(n-2) * (1 + y^2) +
	- 2*k^2*n^2*x^(2*n-2) * y*(1 + y^2) # = 0

### ODE:
d2y - k*n*x^(n-2) * (2*k*n*x^n * y^3 + (n-1)*y^2 + 2*k*n*x^n * y + n-1) # = 0


### Variant 2:
d2y - k*n*(n-1)*x^(n-2) * (1 + y^2) +
	- 2*k^2*n^2*x^(2*n-2) * y*(1 + y^2) # = 0

### ODE:
x*d2y - 2*k*n*x^n * y*dy - (n-1)*dy # = 0


#############################

### Gen: y = x^p * tan(k*x^n)

# Check:
e = expression(x^p * tan(k*x^n));
k = sqrt(5);
n = 1/sqrt(2); p = 2/5; # p = n;
x = sqrt(3);
params = list(x=x, k=k);
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*dy - p*y - k*n*x^(n-p) * (y^2 + x^(2*p)) # = 0

# D2 =>
x*d2y - (p-1)*dy - k*n*(n-p)*x^(n-p-1) * (y^2 + x^(2*p)) +
	- 2*k*n*x^(n-p) * (y*dy + p*x^(2*p-1)) # = 0

### Variant 1:
x^2 * d2y - (p-1)*x*dy - (n-p) * (x*dy - p*y) +
	- 2*k*n*x^(n-p+1) * (y*dy + p*x^(2*p-1)) # = 0

### ODE:
x^2 * d2y - 2*k*n*x^(n-p+1) * y*dy +
	- (n-1)*x * dy + p*(n-p) * y - 2*p*k*n*x^(n+p) # = 0


### Special Cases:

### Case: p = n;
x*d2y - 2*k*n*y*dy - (n-1)*dy - 2*k*n^2*x^(2*n-1) # = 0

