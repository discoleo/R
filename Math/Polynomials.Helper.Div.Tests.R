########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Multi-Variable Polynomials
### Polynomial Division: Tests
###
### draft v.0.1a


### Tests:
### Division in Multi-Variable Polynomials


######################

### Helper Functions

### fast load:
# source("Polynomials.Helper.Div.Tests.R")


source("Polynomials.Helper.R")
# - is automatically loaded in: Polynomials.Helper.R;
# source("Polynomials.Helper.Div.R")


########################
########################

b1 = 3;
p1   = toPoly.pm("x^3 + b1()*x + 1")
pDiv = toPoly.pm("a*x + 1")
pR = div.pm(p1, pDiv, by="x", NF.stop=FALSE)
pR # TODO: fix bugs;
