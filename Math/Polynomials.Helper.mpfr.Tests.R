########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: Helper Functions
### mpfr Functions: Tests
###
### draft v.0.2a


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
pd = 2
x0 = c(2,3,5,7,13,11,4,6,8,9,10,12)
x = as.pow.mpfr(x0, 200, pow.div = pd)

mv = vandermonde.mpfr(x)
det.mpfr(mv, normalize = FALSE)
det.mpfr(mv, normalize = TRUE) # OK: but why?
determinant.seq.mpfr(mv)
det.vandermonde.mpfr(mv)


###
pd = 2; prec = 200;
x0 = 2:35;
x = as.pow.mpfr(x0, prec, pow.div = pd)

mv = vandermonde.mpfr(x)
det.mpfr(mv, normalize = FALSE)
det.mpfr(mv, normalize = TRUE) # minimal
determinant.seq.mpfr(mv) # slight improvement;
det.vandermonde.mpfr(mv)


###
pd = 2; prec = 240;
x0 = 2:35
x = as.pow.mpfr(x, prec, pow.div = pd)

mv = vandermonde.mpfr(x)
det.mpfr(mv, normalize = FALSE)
det.mpfr(mv, normalize = TRUE) # NO effect
det.vandermonde.mpfr(mv)


###
pd = 3; prec = 220;
x0 = 2:35 * rep(c(1,-1), 17)
x = as.pow.mpfr(x0, prec, pow.div = pd)

mv = vandermonde.mpfr(x)
det.mpfr(mv, normalize = FALSE)
det.mpfr(mv, normalize = TRUE) # NO effect
determinant.seq.mpfr(mv) # NO effect
det.vandermonde.mpfr(mv) # Maths formula
det(as.matrix(mv)) # leading digit wrong;


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


#############
### Solve ###

prec = 200;
xup = 8;
div = 2; # div = 3; # x^2 vs x^3;
x = as.pow.mpfr(2:xup, prec, pow.div = div)
b = vandermonde.mpfr(x)

#
y = seq(length(x));
solve.mpfr(b, y)
solve(as.matrix(b), y)

# Pow: (x^div - 1)^2
solve.mpfr(b, y^2)
solve(as.matrix(b), y^2)

# just some values:
y = c(1,3,7,5,2,2, sample(1:100, length(x) - 6, TRUE));
solve.mpfr(b, y)
solve(as.matrix(b), y)


### Test: Swap columns
prec = 200;
xup = 12;
div = 2; # div = 3;
x = as.pow.mpfr(2:xup, prec, pow.div = div)
b = vandermonde.mpfr(x)
b[1:7, 5] = b[1:7, 4] = b[1:7, 3]; # force swap columns;

y = c(1,3,7,5,2,2, sample(1:100, length(x) - 6, TRUE));
solve.mpfr(b, y)
solve(as.matrix(b), y)


### Complex
prec = 200;
xup = 8;
div = 3;
xr = as.pow.mpfr(2:xup, prec, pow.div = div)
xi = as.pow.mpfr((xup:2) * sample(c(-1,1), length(xr), TRUE), prec, pow.div = div)
b  = vandermonde.complex.mpfr(xr, xi)

#
yr = -2 + 1.5 * b$Re[,2] - 3*b$Re[,3] + b$Re[,4];
yi = -4 + 1.5 * b$Im[,2] - 3*b$Im[,3] + b$Im[,4];
solve.complex.mpfr(b$Re, b$Im, yr, yi)


###################
### Polynomials ###
###################

### Roots


###
b = 2; # b = 1;
n = 17;
p = as.pm("x^n - x - b[1]")

# Note: very sensitive to x0 & x0i!
x0i = - 0.2; # x0i = -0.4; # both for b = 2;
# x0i = - 0.3; # for b = 1;
x = roots.mpfr(p, 1, x0i, it = 64, prec = prec)
xx = pow.all.complex.mpfr(x$Re, x$Im, n);
xx$Re[n] - xx$Re[1]; xx$Im[n] - xx$Im[1]


###
b = c(-1, -2); # b = c(-2, -2)
# b = c(-1,-1)
p = as.pm("x^17 + x^15 + b[2]*x + b[1]")

prec = 200
# Note: very sensitive to x0 & x0i!
x0i = - 0.4; # x0i = - 0.8; # all 3 variants;
x = roots.mpfr(p, 1, x0i, it = 64, prec = prec)
xx = pow.all.complex.mpfr(x$Re, x$Im, 17);
xx$Re[17] + xx$Re[15] + b[2]*xx$Re[1];
xx$Im[17] + xx$Im[15] + b[2]*xx$Im[1];


### Multiplicity
b = 1; # b = 2;
p = as.pm("(x^5 - x - b[1])^3")

#
prec = 200;
x = roots.mpfr(p, -0.1, -0.2, it = 160, prec = prec)
xx = pow.all.complex.mpfr(x$Re, x$Im, 5); xx$Re[5] - xx$Re[1]; xx$Im[5] - xx$Im[1]

#
prec = 240;
x = roots.mpfr(p, -0.1, -0.2, it = 160, prec = prec)
xx = pow.all.complex.mpfr(x$Re, x$Im, 5); xx$Re[5] - xx$Re[1]; xx$Im[5] - xx$Im[1]


#############

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
