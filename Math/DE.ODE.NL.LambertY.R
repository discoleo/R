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

