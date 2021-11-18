########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Multi-Variable Polynomials
### Factorize: Tests
###
### draft v.0.1a


### Tests:
### Factorize Multi-Variable Polynomials


######################

### Helper Functions

### fast load:
# source("Polynomials.Helper.Factorize.Tests.R")


# - is automatically loaded in: Polynomials.Helper.R;
# source("Polynomials.Helper.Factorize.R")


########################
########################

b1 = 3;
p0 = toPoly.pm("x^2 + b1()*x + 1")
p1 = toPoly.pm("x^3 - 4*x^2 - x + 1")
p2 = toPoly.pm("p0()*p1()")

eval.pm(p2, roots.pm(p1)[1])
p2

### Factors: "strictly" Symmetric Polynomials
### P(x) * P(1/x)
p2
rev(p2)
div.pm(rev(p2), rev(p1), "x")

###
pX = mult.pm(p2, rev(p2))

# Note: overflows massively;
gcd.exact.p(pX, dp.pm(pX, "x"), "x", asBigNum=FALSE, debug=TRUE)

### using BigNumbers
source("Polynomials.Helper.BigNumbers.R")
pR = factorize.p(toBigz.pm(pX), xn="x", f.all=FALSE, asBigNum=TRUE, file=NULL)
str(pR)
pR[[1]]$GCD

### faster version
gcd.exact.p(p2, rev(p2), xn="x", asBigNum=FALSE, debug=TRUE)


### Multiple Techniques:
factorizeExt.p(p2, xn="x", asBigNum=FALSE, debug=T)

### Anti-Symmetric:
b1 = 3;
p0 = toPoly.pm("x^2 + b1()*x - 1")
p1 = toPoly.pm("x^3 - 4*x^2 - x + 1")
p2 = toPoly.pm("p0()*p1()")
#
factorizeExt.p(p2, xn="x", asBigNum=FALSE, debug=F)

### TODO
b1 = 3;
p0 = toPoly.pm("x^2 + b1()*x + 1")
p3 = toPoly.pm("p0() * (x^4 + 5*x^3 + 5*x + 1)")
factorizeExt.p(p3, xn="x", asBigNum=FALSE, debug=F)


################
################

eval.pm(p2, 4) %% 35
eval.pm(p2, 9) %% 35
