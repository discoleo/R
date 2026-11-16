########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Lambert W
##
## draft v.0.1c



####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")

#####################
#####################

### y^y = exp(k*x)
# TODO: Check;

# D =>
(log(y) + 1)*dy * y^y - k*exp(k*x) # = 0
(log(y) + 1)*dy - k # = 0

# D2 =>
y*(log(y) + 1)*d2y + dy^2 # = 0

### ODE:
k*y*d2y + dy^3 # = 0


#####################

### y^y = exp(k*x^n)
# TODO: Check;

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

### Radicals

### y^y = sqrt(x^2 + 1)
# TODO: Check;

# D =>
(log(y) + 1)*dy * y^y - x*sqrt(x^2 + 1) / (x^2+1) # = 0
(x^2+1)*(log(y) + 1)*dy - x # = 0

# D2 =>
(x^2+1)*(log(y) + 1)*y*d2y + (x^2+1) * dy^2 + 2*x*(log(y) + 1)*y*dy - y # = 0
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

