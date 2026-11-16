########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Lambert W
##
## draft v.0.1d

## NL ODEs from Coupled/Persisting Terms:
# y^y = ... OR y * Exp(y) = ...,
# where D(RHS) = RHS * ...;


### Examples:

# k*y*d2y + dy^3 = 0;
# k*x*y*d2y - x^3 * dy^3 + 2*k*y*dy = 0;
# x*y*d2y + (x^2+1) * dy^3 + (x^2-1)/(x^2+1) * y*dy = 0;
# x*y*(y+1) * d2y - x*dy^2 - k*y^2 = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")

#####################
#####################

### y^y = exp(k*x)

# Check:
k = 1/sqrt(3);
x = sqrt(5); params = list(x=x, k=k);
w = pracma::lambertWp(k*x);
y = exp(w); params$y = y;
dy = y/x * w/(w+1);
ye = expression(k / (log(y) + 1))[[1]];
ye = D(ye, "y"); d2y = eval(ye, params) * dy;

# D =>
(log(y) + 1)*dy * y^y - k*exp(k*x) # = 0
(log(y) + 1)*dy - k # = 0

# D2 =>
y*(log(y) + 1)*d2y + dy^2 # = 0

### ODE:
k*y*d2y + dy^3 # = 0

# Alternative Formulation:
Ki = 1/k;
y*d2y + Ki*dy^3 # = 0


#####################

### y^y = exp(k*x^n)

# Check:
k = 1/sqrt(3); n = exp(-1); # n = -1;
x = sqrt(5); params = list(x=x, k=k, n=n);
w = pracma::lambertWp(k*x^n);
y = exp(w); params$y = y;
dy = n*y/x * w/(w+1);
ye = expression(n*k*x^(n-1) / (log(y) + 1))[[1]];
yy = D(ye, "y"); yx = D(ye, "x");
d2y = eval(yx, params) + eval(yy, params) * dy;

# D =>
(log(y) + 1)*dy * y^y - n*k*x^(n-1) * exp(k*x^n) # = 0
(log(y) + 1)*dy - n*k*x^(n-1) # = 0

# D2 =>
y*(log(y) + 1)*d2y + dy^2 - k*n*(n-1)*x^(n-2) * y # = 0

### ODE:
k*n*x^(n-1) * y*d2y + dy^3 - k*n*(n-1)*x^(n-2) * y*dy # = 0

### Special Cases:

### Case: n = -1
k*x*y*d2y - x^3 * dy^3 + 2*k*y*dy # = 0


########################
########################

### y*exp(y) = x^(k*x)

# Check:
k = 1/sqrt(3);
x = sqrt(5); params = list(x=x, k=k);
y = pracma::lambertWp(x^(k*x)); params$y = y;
dy = k*(log(x)+1) * y/(y+1);
ye = expression(k*(log(x)+1)*x^(k*x) / ((y+1)*exp(y)))[[1]];
yy = D(ye, "y"); yx = D(ye, "x");
d2y = eval(yx, params) + eval(yy, params) * dy;

# D =>
(y+1)*exp(y) * dy - k*(log(x)+1)*x^(k*x) # = 0
(y+1) * dy - k*(log(x)+1)*y # = 0

# D2 =>
x*(y+1) * d2y + x*dy^2 - k*x*(log(x)+1)*dy - k*y # = 0
x * y*(y+1) * d2y + x*y*dy^2 - x*(y+1) * dy^2 - k*y^2 # = 0

### ODE:
x*y*(y+1) * d2y - x*dy^2 - k*y^2 # = 0


########################
########################

### Radicals

### y^y = sqrt(x^2 + 1)

# Check:
x = sqrt(5); params = list(x=x);
w = pracma::lambertWp(log(x^2+1)/2);
y = exp(w); params$y = y;
dy = x / (x^2+1) * 2/log(x^2+1) * y * w/(w+1);
ye = expression(x / (x^2+1) / (log(y) + 1))[[1]];
yy = D(ye, "y"); yx = D(ye, "x");
d2y = eval(yx, params) + eval(yy, params) * dy;

# D =>
(log(y) + 1)*dy * y^y - x*sqrt(x^2 + 1) / (x^2+1) # = 0
(x^2+1)*(log(y) + 1)*dy - x # = 0

# D2 =>
(x^2+1)*(log(y) + 1)*y*d2y + (x^2+1) * dy^2 + 2*x*(log(y) + 1)*y*dy - y # = 0

### ODE:
x * y*d2y + (x^2+1) * dy^3 + (x^2-1)/(x^2+1) * y*dy # = 0


########################

### sqrt(y^2 + 1) = x^x

# Check:
x = sqrt(5);
params = list(x=x);
ye = expression(sqrt(x^(2*x) - 1))[[1]];
y  = eval(ye, params);
ye = D(ye, "x"); dy  = eval(ye, params);
ye = D(ye, "x"); d2y = eval(ye, params);

# D =>
y*dy * sqrt(y^2 + 1) / (y^2 + 1) - (log(x)+1)*x^x # = 0
y*dy - (log(x)+1) * (y^2 + 1) # = 0

# D2 =>
x*y*d2y + x*dy^2 - 2*x*(log(x)+1) * y*dy - (y^2 + 1) # = 0
x*y*d2y + x*dy^2 - 2*x * y^2*dy^2 / (y^2 + 1) - (y^2 + 1) # = 0

### ODE:
x*y*(y^2 + 1)*d2y - x*y^2*dy^2 + x*dy^2 - (y^2 + 1)^2 # = 0

