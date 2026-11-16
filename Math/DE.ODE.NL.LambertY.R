########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Lambert W
##
## draft v.0.1b



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

