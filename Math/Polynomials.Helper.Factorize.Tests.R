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


source("Polynomials.Helper.R")
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

genPoly = function(b01, b02) {
	p0 = toPoly.pm("x^2 + b1()*x + b01()");
	p1 = toPoly.pm("p0() * (x^3 + 2*x^2 - 5*x + b02())");
	return(p1);
}
### Ex 1:
b1 = 3; b01 = 4; b02 = 1;
p1 = genPoly(b01, b02);
factorizeByB0.p(p1, xn="x")

### Ex 2:
b1 = 3; b01 = -9; b02 = 4;
p1 = genPoly(b01, b02);
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
p = toPoly.pm("(x^4 + b()*x^2 + 2)*(x^3 - 2*x +3)");
#
p1 = replace.pm(p, toPoly.pm("1i*x"), xn="x")
p2 = replace.pm(p, toPoly.pm("-1i*x"), xn="x")
#
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

eval.pm(p2, 4) %% 35
eval.pm(p2, 9) %% 35
