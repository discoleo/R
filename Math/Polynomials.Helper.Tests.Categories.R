########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Multi-Variable Polynomials
### Categories: Tests
###
### draft v.0.1a


### Tests
### Multi-Variable Polynomials: Categories


### this file:
# options(warn=1); source("Polynomials.Helper.Tests.Categories.R")


####################

### Helper Functions

source("Polynomials.Helper.Categories.R");


####################
####################

### Symmetric Polynomials

### Tests:
checkVal.pm(isSymmetric.pm(toPoly.pm("x^5 + 1")), TRUE)
checkVal.pm(isSymmetric.pm(toPoly.pm("x^4 + 1")), TRUE)
checkVal.pm(isSymmetric.pm(toPoly.pm("1")), FALSE)
checkVal.pm(isSymmetric.pm(toPoly.pm("0")), TRUE)
#
checkVal.pm(isSymmetric.pm(toPoly.pm("x^5")), FALSE)
checkVal.pm(isSymmetric.pm(toPoly.pm("x^5 + 2")), FALSE)
checkVal.pm(isSymmetric.pm(toPoly.pm("x^4 + 2")), FALSE)

###
checkVal.pm(isSymmetric.pm(toPoly.pm("x^4 + 3*x^2 + 1")), TRUE)
checkVal.pm(isSymmetric.pm(toPoly.pm("x^5 + 3*x^2 + 1")), FALSE)
checkVal.pm(isSymmetric.pm(toPoly.pm("x^4 + 3*x + 1")), FALSE)

###
checkVal.pm(isSymmetric.pm(toPoly.pm("x^5 + 3*x^3 + 3*x^2 + 1")), TRUE)
checkVal.pm(isSymmetric.pm(toPoly.pm("x^4 + 3*x^3 + 3*x + 1")), TRUE)
checkVal.pm(isSymmetric.pm(toPoly.pm("x^4 + 3*x^3 + 4*x^2 + 3*x + 1")), TRUE)
checkVal.pm(isSymmetric.pm(toPoly.pm("x^8 + 3*x^5 + 4*x^4 + 3*x^3 + 1")), TRUE)
checkVal.pm(isSymmetric.pm(toPoly.pm("x^9 + 3*x^6 + 4*x^5 + 4*x^4 + 3*x^3 + 1")), TRUE)

###
checkVal.pm(isSymmetric.pm(toPoly.pm("x^8 + 3*x^5 + 4*x^2 + 3*x^3 + 1")), FALSE)
checkVal.pm(isSymmetric.pm(toPoly.pm("x^9 + 3*x^6 + 4*x^5 + 5*x^4 + 3*x^3 + 1")), FALSE)


###################
###################

cat("\n\nAll Tests: Success!\n");

