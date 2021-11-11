########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Derivation: Tests
###
### draft v.0.1a


### requires:
source("Polynomials.Helper.D.R")


### fast load:
# source("Polynomials.Helper.D.Test.R")

#####################
#####################

###
p1 = toPoly.pm("(x+1)^3")
p2 = toPoly.pm("x+1")

split.pm.fraction(p1, p2)


###
p1 = toPoly.pm("(x+a)^3")
p2 = toPoly.pm("x+a")

split.pm.fraction(p1, p2)


###
p1 = toPoly.pm("(x+a)^3 + b")
p2 = toPoly.pm("x+a")

split.pm.fraction(p1, p2)


###
p1 = toPoly.pm("(x+a1)^2*(x+a2) + x^2 + b")
p2 = toPoly.pm("(x+a1)*(x+a2)")
#
pR = split.pm.fraction(p1, p2)
pR
diff.pm(pR$P0, toPoly.pm("x + a1 + 1"))


######################
######################

###
p = toPoly.pm(data.frame(x=5:0, coeff=1));
#
I.pm(p, xn="x") * 6


###
p1 = mult.pm(p, toPoly.pm("x^2 + b1*x + b0"))
p1 = toPoly.pm(p1)
#
pR = I.pm(p1, xn="x")
pR = sort.pm(pR, "x")
pR * 7

###
p1 = toPoly.pm("x^2 + b1*x + b0 + blog*x^-1")
pR = I.pm(p1, xn="x")
pR = sort.pm(pR, "x")
pR * 3