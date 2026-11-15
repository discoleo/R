########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Lambert W
##
## draft v.0.1a



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

