########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Lambert W
##
## draft v.0.1a



#########################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y * log(y)^2 = F0(x)

# y = LambertW(1/2 * sqrt(F0))^2

# Check:
# TODO


# D =>
log(y)^2 * dy + 2*log(y)*dy - df0 # = 0
2*y*log(y)*dy + f0 * dy - df0*y # = 0

# D2 =>
2*y*log(y)*d2y + 2*dy^2 + 2*log(y)*dy^2 +
	+ f0 * d2y - d2f*y # = 0
2*y^2*log(y)*d2y + 2*y*dy^2 - (f0*dy - df0*y)*dy +
	+ f0 * y*d2y - d2f*y^2 # = 0
- y*(f0 * dy - df0*y)*d2y + 2*y*dy^3 - (f0*dy - df0*y)*dy^2 +
	+ f0 * y*dy*d2y - d2f*y^2*dy # = 0

###  ODE:
df0*y^2*d2y + 2*y*dy^3 - f0*dy^3 + df0*y*dy^2 - d2f*y^2*dy # = 0

