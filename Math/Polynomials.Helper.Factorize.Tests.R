########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Multi-Variable Polynomials
### Factorize: Tests
###
### draft v.0.1d


### Tests:
### Factorize Multi-Variable Polynomials


######################

### Helper Functions

### fast load:
# source("Polynomials.Helper.Factorize.Tests.R")


source("Polynomials.Helper.R")
# - is automatically loaded in: Polynomials.Helper.R;
# source("Polynomials.Helper.Factorize.R")


########################
########################

b  = 3;
p0 = toPoly.pm("x^2 + b[1]*x + 1")
p1 = toPoly.pm("x^3 - 4*x^2 - x + 1")
p2 = toPoly.pm("p0()*p1()")

err = eval.pm(p2, roots.pm(p1)[1])
err
stopifnot(round0(err) == 0)
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
b  = 3;
p0 = toPoly.pm("x^2 + b[1]*x - 1")
p1 = toPoly.pm("x^3 - 4*x^2 - x + 1")
p2 = toPoly.pm("p0()*p1()")
#
pR = factorizeExt.p(p2, xn="x", asBigNum=FALSE, debug=F)
stopifnot( ! is.null(pR[[1]]$GCD))
print(pR[[1]]$GCD)

### TODO
b  = 3;
p0 = toPoly.pm("x^2 + b[1]*x + 1")
p3 = toPoly.pm("p0() * (x^4 + 5*x^3 + 5*x + 1)")
pR = factorizeExt.p(p3, xn="x", asBigNum=FALSE, debug=F)
stopifnot( ! is.null(pR[[1]]$GCD))
stopifnot( max(pR[[1]]$GCD$x) == 6)
print(pR[[1]]$GCD)


################
################

genPoly = function(b) {
	p0 = toPoly.pm("x^2 + b[1]*x + b[2]");
	p1 = toPoly.pm("p0() * (x^3 + 2*x^2 - 5*x + b[3])");
	return(p1);
}
### Ex 1:
b  = c(3, 4, 1);
p1 = genPoly(b);
factorizeByB0.p(p1, xn="x")

### Ex 2:
b  = c(3, -9, 4);
p1 = genPoly(b);
factorizeByB0.p(p1, xn="x")


###
p = rPoly(10)

### using BigNumbers
source("Polynomials.Helper.BigNumbers.R")
# Note: the code becomes usually incompatible after switching to Bigz; 
p = toBigz.pm(p);
factorizeExt.p(p, xn="x", asBigNum=TRUE, debug=T)


################
################

### Generation of Squares
b = 3
p = toPoly.pm("(x^4 + b[1]*x^2 + 2)*(x^3 - 2*x +3)");
#
p1 = replace.pm(p, toPoly.pm("1i*x"), xn="x")
p2 = replace.pm(p, toPoly.pm("-1i*x"), xn="x")
# with BIGZ:
p3 = mult.pm(p1, p2);
p3 = as.bigz.pm(p3);
pR = factorize.p(p3, "x", asBigNum=TRUE, file=NULL)
pR[[1]]$GCD # TODO: scale back by i^2;

# alternative:
# - but ugly algorithm, due to complex numbers;
# - fails in gcd(a, b);
# gcd.exact.p(p1, p2, asBigNum=FALSE)

# alternative 2:
p1 = rescale.pm(p, -1, "x")
pR = gcd.exact.p(p, p1, asBigNum=FALSE)
pR
div.pm(p, pR, by="x")


################
################

b = c(2,3,5)
plst = lapply(b, function(b) toPoly.pm("x^2 + b[1]*x + 1"))
# a symmetric Polynomial:
p = mult.lpm(plst)


################
################

b = c(8, 2)
p1 = toPoly.pm("x^2 + b[1]*x + b[2]")
p2 = toPoly.pm("x^3 - x^2 + 2*x + 1")
p = toPoly.pm(p1 * p2)

# 3^2 = 2 (mod 7)
pM = rescale.pm(p, 3, mod=7)
pMinv = rev.pm(pM)
#
pGCD1 = toPoly.pm(diff.pm(2*pM, 5*pMinv)) %% 7
pGCD2 = toPoly.pm(diff.pm(2*pM, 5*pGCD1 * "x")) %% 7
pGCD2
pGCD1 = toPoly.pm(diff.pm(3*pGCD1, 2*pGCD2)) %% 7
pGCD2 = toPoly.pm(diff.pm(2*pGCD2, pGCD1 * toPoly.pm("x^2 + x"))) %% 7
pGCD1; pGCD2
# "6 + 2*x + 6*x^2"
# "1 + 5*x + x^2" # 6*5 = 30 = 2 (mod 7)

# TODO: gcd.pm.mod(p1, p2, mod)
# gcd.exact.p(pM, pMinv, asBigNum=F)


################
################

eval.pm(p2, 4) %% 35
eval.pm(p2, 9) %% 35
