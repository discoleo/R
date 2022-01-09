########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: Helper Functions
### Tests: Roots



### this file:
# options(warn=1); source("Polynomials.Helper.Tests.Roots.R")


####################

### Helper Functions

### fast load:
source("Polynomials.Helper.Tests.Helper.R")


#####################
#####################

###
n = 3; val = 1;
p = toPoly.pm("(x - val[1])^n")
m = multiplicity.pm(p, val=val)
checkVal.pm(m, n)

###
n = 4; val = 1;
p = toPoly.pm("(x - val[1])^n")
m = multiplicity.pm(p, val=val)
checkVal.pm(m, n)

###
n = 3; val = 2;
p = toPoly.pm("(x - val[1])^n")
m = multiplicity.pm(p, val=val)
checkVal.pm(m, n)

###
n = 4; val = 3;
p = toPoly.pm("(x - val[1])^n")
m = multiplicity.pm(p, val=val)
checkVal.pm(m, n)


### Composite:
n = 7; val = 1;
p = toPoly.pm("(x^3 - val[1])^n * (x^5 - x - 1)")
m = multiplicity.pm(p, val=rootn(val, 3))
checkVal.pm(m, n)


### Composite:
n = 7; val = 2;
p = toPoly.pm("(x^3 - val[1])^n * (x^5 - x - 1)")
# FAILS due to Tolerance!
# TODO: mpfr; (roots cannot be BigZ)
m = multiplicity.pm(p, val=rootn(val, 3), tol=1E-5)
checkVal.pm(m, n)

