########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Lambert W
##
## draft v.0.1h

## NL ODEs from Coupled/Persisting Terms:
# y^y = ... OR y * Exp(y) = ...,
# where D(RHS) = RHS * ...;


### Examples:

# k*y*d2y + dy^3 = 0;
# k*x*y*d2y - x^3 * dy^3 + 2*k*y*dy = 0;
# x*(k*x+1)*y*d2y + x^2 * dy^3 + y*dy = 0;
# x*y*d2y + (x^2+1) * dy^3 + (x^2-1)/(x^2+1) * y*dy = 0;
# x*y*(y+1)*d2y - x*dy^2 - k*y^2 = 0;
# n*x*y*(y+1)*d2y + x*(y^2 + 2*y - n+1) * dy^2 + n*y*(y+1)*dy = 0;
# (x^2+k^2) * (y+1)*d2y + (x^2+k^2) * (y+2)*dy^2 + 2*x * (y+1)*dy = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")

# Note:
# Case: y^y = F(x)
# => w = pracma::lambertWp(log(F)); y = exp(w);
#    D(y) = y * D(F)/F * w/(w+1);
#
# Case: y*exp(y) = F(x)
# => w = pracma::lambertWp(F); y = w;
#    D(y) = D(F)/F * w/(w+1);

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


######################

### y^y = x * exp(k*x)

# Check:
k = 1/sqrt(3);
x = sqrt(5); params = list(x=x, k=k);
w = pracma::lambertWp(k*x + log(x));
y = exp(w); params$y = y;
dy = y * (k*x+1)/x / (k*x + log(x)) * w/(w+1);
ye = expression((k*x+1)/x / (log(y) + 1))[[1]];
yy = D(ye, "y"); yx = D(ye, "x");
d2y = eval(yx, params) + eval(yy, params) * dy;

# D =>
(log(y) + 1) * dy * y^y - (k*x+1) * exp(k*x) # = 0
x*(log(y) + 1) * dy - (k*x+1) # = 0

# D2 =>
x*(log(y) + 1) * y*d2y + (log(y) + 1) * y*dy + x*dy^2 - k*y # = 0

### ODE:
x*(k*x+1) * y*d2y + x^2 * dy^3 + y*dy # = 0


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


########################
########################

### Other

### y*exp(y) = log(x)^n

# Check:
n = sqrt(3);
x = sqrt(5); params = list(x=x, n=n);
y = pracma::lambertWp(log(x)^n); params$y = y;
dy = n/x / log(x) * y/(y+1);
ye = expression(n*y / (x*log(x) * (y+1)))[[1]];
yy = D(ye, "y"); yx = D(ye, "x");
d2y = eval(yx, params) + eval(yy, params) * dy;

# D =>
x*(y+1)*exp(y) * dy - n*log(x)^(n-1) # = 0
x*log(x) * (y+1) * dy - n*y # = 0

# D2 =>
x*log(x) * (y+1) * d2y + x*log(x)*dy^2 + (log(x)+1) * (y+1) * dy - n*dy # = 0
n*x * y*(y+1)*d2y + n*x * y*dy^2 + x*(y+1)*(y - (n-1)) * dy^2 + n*y*(y+1)*dy # = 0

### ODE:
n*x * y*(y+1)*d2y + x*(y^2 + 2*y - n+1) * dy^2 + n*y*(y+1)*dy # = 0


##########################

### y*exp(y) = atan(x/k)^n

# Check:
k = sqrt(2); n = sqrt(3); # n = 1;
x = sqrt(5); params = list(x=x, n=n, k=k);
y = pracma::lambertWp(atan(x/k)^n); params$y = y;
dy = n*k/(x^2+k^2) / atan(x/k) * y/(y+1);
ye = expression(n*k*y / ((x^2+k^2)*atan(x/k) * (y+1)))[[1]];
yy = D(ye, "y"); yx = D(ye, "x");
d2y = eval(yx, params) + eval(yy, params) * dy;

# D =>
(x^2+k^2)*(y+1)*exp(y) * dy - n*k * atan(x/k)^(n-1) # = 0
(x^2+k^2)*atan(x/k) * (y+1)*dy - n*k * y # = 0

# D2 =>
(x^2+k^2)*atan(x/k)*(y+1) * d2y + (x^2+k^2)*atan(x/k)*dy^2 +
	+ 2*x*atan(x/k)*(y+1) * dy + k*(y - n+1) * dy # = 0
n*k*(x^2+k^2) * y*(y+1)*d2y + n*k*(x^2+k^2) * y*dy^2 +
	+ k*(x^2+k^2) * (y+1)*(y - n+1) * dy^2 + 2*n*k*x * y*(y+1)*dy # = 0

### ODE:
n*(x^2+k^2) * y*(y+1)*d2y + (x^2+k^2) * (y^2 + 2*y - n+1) * dy^2 +
	+ 2*n*x * y*(y+1)*dy # = 0


### Special Cases:

### Case: n = 1;
(x^2+k^2) * (y+1)*d2y + (x^2+k^2) * (y+2)*dy^2 + 2*x * (y+1)*dy # = 0


#######################
#######################

### y*exp(y) = sqrt(log(x)^2 + 1)

# Check:
n = sqrt(3);
x = sqrt(5); params = list(x=x, n=n);
y = pracma::lambertWp(sqrt(log(x)^2 + 1)); params$y = y;
dy = log(x)/x / (log(x)^2+1) * y/(y+1);
ye = expression(log(x)*y / (x*(log(x)^2 + 1) * (y+1)))[[1]];
yy = D(ye, "y"); yx = D(ye, "x");
d2y = eval(yx, params) + eval(yy, params) * dy;

# D =>
x*(y+1)*exp(y) * dy - log(x) / sqrt(log(x)^2 + 1) # = 0
x*(log(x)^2 + 1) * (y+1)*dy - log(x) * y # = 0

# D2 =>
x*log(x) * y*d2y + x^2*(log(x)^2 + 1) * dy * dy^2 +
	+ x*log(x) * (2*y+1)*dy^2 + log(x) * y*dy - y*dy # = 0
x*log(x) * y*(y+1)*d2y + x*log(x) * (2*y^2 + 4*y + 1)*dy^2 +
	+ log(x) * y*(y+1)*dy - y*(y+1)*dy # = 0

### ODE:
x^3 * y^2*(y+1)^2*d2y^2 +
	+ 2*x^3 * y*(y+1)*(2*y^2 + 4*y + 1) * dy^2*d2y +
	+ 2*x^2 * y^2*(y+1)^2 * dy*d2y  -  x*y^3*(y+1)*d2y +
	+ x^3 * (2*y^2 + 4*y +1)^2 * dy^4 +
	+ 2*x^2 * y*(y+1)*(2*y^2 + 4*y + 1) * dy^3 +
	+ x*y^2 * dy^2 - y^3*(y+1)*dy # = 0

