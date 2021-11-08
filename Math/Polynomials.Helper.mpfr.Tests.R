########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: Helper Functions
### mpfr Functions: Tests
###
### draft v.0.1a


### fast load:
# source("Polynomials.Helper.mpfr.Tests.R")


### Requirements:

source("Polynomials.Helper.mpfr.R")


### Libraries:
# library(Rmpfr)

### Other requirements:
# source("Polynomials.Helper.R")


#######################
#######################

#############
### Tests ###
#############

###
K = 3
x = roots.Class1.mpfr(K, c(1,-K,0,1))
poly.calc.mpfr(x)
# x^5 - 15*x^3 + 1305*x + 1293

###
K = 2
x = roots.Class1.mpfr(K, c(1,-K,0,1))
poly.calc.mpfr(x)
# x^5 - 10*x^3 + 200*x - 50
p = toPoly.pm("x^5 - 10*x^3 + 200*x - 50")
eval.cpm(p, as.double.cmpfr(x)[1])
eval.cpm(p, x[1, ]) # better accuracy
eval.cpm(p, x[x[,2] == 0, ]) # better accuracy

