########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: Helper Functions
### mpfr Functions: Tests
###
### draft v.0.1b


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


##############
### Matrix ###
##############

# Note: Scale by k[j]
# - should NOT have any relevant impact on precision,
#   as all entries affecting column j get scaled by k[j];
# => x[i,j] = k[j]*(x[i,j] - x[i,1]*x[1,j]/x[1,1]);

###
pow = 0.5
x = c(2,3,5,7,13,11,4,6,8,9,10,12)
x = (mpfr(x, 200))^pow

mv = vandermonde.mpfr(x)
det.mpfr(mv, normalize = FALSE)
det.mpfr(mv, normalize = TRUE) # OK: but why?
determinant.seq.mpfr(mv)
det.vandermonde.mpfr(mv)


###
prec = 200
x = 2:35
x = (mpfr(x, prec))^0.5

mv = vandermonde.mpfr(x)
det.mpfr(mv, normalize = FALSE)
det.mpfr(mv, normalize = TRUE) # minimal
determinant.seq.mpfr(mv) # slight improvement;
det.vandermonde.mpfr(mv)


###
x = 2:35
x = (mpfr(x, 240))^0.5

mv = vandermonde.mpfr(x)
det.mpfr(mv, normalize = FALSE)
det.mpfr(mv, normalize = TRUE) # NO effect
det.vandermonde.mpfr(mv)


### Complex

###
x0 = c(2,3,5,7,11,7,2,13,15,14)
xi0 = c(5,3,7,2,2,3,4,5,0,2)

prec = 200
x  = mpfr(x0, prec)
xi = mpfr(xi0, prec)
x  = sqrt(x); xi = sqrt(xi)

m  = vandermonde.complex.mpfr(x, xi)
m0 = as.numeric(m$Re) + 1i*as.numeric(m$Im);
dim(m0) = dim(m$Re)

determinant.complex(m0)
det.complex.mpfr(m$Re, m$Im)
det.vandermonde.complex.mpfr(x, xi)


###################
### Polynomials ###
###################

###
K = 3
x = roots.Class1.mpfr(K, c(1,-K,0,1))
poly.calc.mpfr(x)
# x^5 - 15*x^3 + 1305*x + 1293

### Polar coordinates:
toPolar.mpfr(x)

###
K = 2
x = roots.Class1.mpfr(K, c(1,-K,0,1))
poly.calc.mpfr(x)
# x^5 - 10*x^3 + 200*x - 50
p = toPoly.pm("x^5 - 10*x^3 + 200*x - 50")
eval.cpm(p, as.double.cmpfr(x)[1])
eval.cpm(p, x[1, ]) # better accuracy
eval.cpm(p, x[x[,2] == 0, ]) # better accuracy


####################

### m^3 = -1
n = 3
m = unity.mpfr(2*n, all=F)
toPolar.mpfr(m)
pow.mpfr(m, 3)
log.cmpfr(pow.mpfr(m, 3))

###
m1 = mpfr(1, 120) / 10^35;
m1 = mpfr2array(c(m1, mpfr(0, 120)), 2)
m1 = m + m1;
m1;
# TODO:
log.cmpfr(pow.mpfr(m1, 3)[3,] + c(1,0))
