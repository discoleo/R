########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Tests for
### Solvers: S2 Systems



####################

### helper Functions

source("Polynomials.Helper.Solvers.S2.R")


### fast load:
# source("Polynomials.Helper.Solvers.S2.Tests.R")


###################

#############
### Tests ###
#############

### S2Ht P3 Simple:
decompose.S2Ht(toPoly.pm("x^3 + b1*y - R"))

### S2Ht P3 + xy:
decompose.S2Ht(toPoly.pm("x^3 + b2*x*y + b1*y - R"))

