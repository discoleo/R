########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Multi-Variable Polynomials
### Factorize: Tests
###
### draft v.0.1i


### Tests:
### Factorize Multi-Variable Polynomials

# Note:
# - existing code covers only univariate polynomials;


######################

### Helper Functions

### fast load:
# source("Polynomials.Helper.Factorize.Tests.R")


source("Polynomials.Helper.Tests.Helper.R")
# - are automatically loaded:
#   source("Polynomials.Helper.R");
# - is automatically loaded in: Polynomials.Helper.R;
#   source("Polynomials.Helper.Factorize.R");

### GMP:
# source("Polynomials.Helper.BigNumbers.R")

########################
########################

b  = 3;
p0 = toPoly.pm("x^2 + b[1]*x + 1")
p1 = toPoly.pm("x^3 - 4*x^2 - x + 1")
p2 = toPoly.pm("p0()*p1()")

# basic Test:
p2
err = eval.pm(p2, roots.pm(p1)[1])
err
stopifnot(round0(err) == 0)

### Factorize
# p0 is symmetric w respect to Inversion: p0(x) == p0(1/x);
# (aka "Strictly symmetric")
pF = gcd.exact.p(p2, rev(p2), asBigNum = FALSE)
print.pm(pF)
stopifnot(nrow(diff.pm(p0, pF)) == 0)


### Simple Test: Divisibility

###
div.pm(p2, pF, by = "x")

###
p2
rev(p2)
div.pm(rev(p2), rev(p1), "x")


### Factors: "strictly" Symmetric Polynomials
### P(x) * P(1/x)
# IF: pF(x) == pF(1/x)
#  => Prod contains pF^2;
# - but it is much simpler to perform gcd(P(x), P(1/x));
pX = mult.pm(p2, rev(p2))

# Note: overflows massively;
gcd.exact.p(pX, dp.pm(pX, "x"), "x", asBigNum=FALSE, debug=TRUE)

### using BigNumbers
if(FALSE) {
source("Polynomials.Helper.BigNumbers.R")
pR = factorize.p(toBigz.pm(pX), xn="x", f.all=FALSE, asBigNum=TRUE, file=NULL)
str(pR)
pR[[1]]$GCD
}

### faster version
gcd.exact.p(p2, rev(p2), by = "x", asBigNum=FALSE, debug=TRUE)


### Multiple Techniques:
factorize.ext.p(p2, xn="x", asBigNum=FALSE, debug=T)


###################

### Anti-Symmetric:
# F(x) = x^2 + b1*x - 1
# F(- 1/x) = 1 - b1*x - x^2 = - F(x);
b  = 3;
p0 = toPoly.pm("x^2 + b[1]*x - 1")
p1 = toPoly.pm("x^3 - 4*x^2 - x + 1")
p2 = toPoly.pm("p0()*p1()")
#
pR = factorize.ext.p(p2, by = "x", asBigNum=FALSE, debug=F)
stopifnot( ! is.null(pR[[1]]$GCD))
print(pR[[1]]$GCD)


################
################

### Asymmetric

genPoly = function(b) {
	p0 = toPoly.pm("x^2 + b[1]*x + b[2]");
	p1 = toPoly.pm("p0() * (x^3 + 2*x^2 - 5*x + b[3])");
	return(p1);
}

### Simple Approach
# Note:
# - only of type x^2 + b1*x + b0^2;
# - more robust approach: in sections below;

### Ex 1:
b  = c(3, 4, 1);
p1 = genPoly(b);
factorizeByB0.p(p1, by="x")

### Ex 2:
b  = c(3, -9, 4);
p1 = genPoly(b);
factorizeByB0.p(p1, by="x")


###################

### Various Methods
### Even Polynomial Factors

### Generation of Squares
# Q = P(1i*x) * P(-1i*x);
# Q = P(x) * P(-x) behaves similarly, but still inefficient;
# Note:
# - very limited and rather inefficient;

b = 3
p = toPoly.pm("(x^4 + b[1]*x^2 + 2)*(x^3 - 2*x + 3)");
#
p1 = replace.pm(p, toPoly.pm("1i*x"), xn="x")
p2 = replace.pm(p, toPoly.pm("-1i*x"), xn="x")

# Note:
# - Direct gcd: more efficient, but still limited usefulness;
# - gcd(P(x), P(-x)) works as well (see Alternative 2);
pR = gcd.pm(p1, p2, by = "x")
mult.pm(pR$Rez, 6)


if(FALSE) {
# with BIGZ:
p3 = mult.pm(p1, p2);
p3 = as.bigz.pm(p3);
pR = factorize.p(p3, "x", asBigNum=TRUE, file=NULL)
pR[[1]]$GCD # TODO: scale back by i^2;
}

# Alternative:
# - but ugly algorithm, due to complex numbers;
# - fails in gcd(a, b);
# gcd.exact.p(p1, p2, asBigNum=FALSE)

# Alternative 2:
p1 = rescale.pm(p, -1, "x")
pR = gcd.exact.p(p, p1, asBigNum=FALSE)
pR
div.pm(p, pR, by="x")


################
################

### Rescale mod p

b = c(8, 2)
p1 = toPoly.pm("x^2 + b[1]*x + b[2]")
p2 = toPoly.pm("x^3 - x^2 + 2*x + 1")
p = toPoly.pm(p1 * p2)


# 3^2 = 2 (mod 7)
# 5^2 = 2 (mod 23)
pR = factorize.mod.p(p, mod = 7, scale = 3, inv.scale = 5)
print.pm(pR)
# (x^2 + x + 2) == (x^2 + 8*x + 2) (mod 7)

pR = factorize.mod.p(p, mod = 23, scale = 5, inv.scale = 14, mult = -1)
print.pm(pR)


###################
###################

### Fully Symmetric

### TODO
b  = 3;
p0 = toPoly.pm("x^2 + b[1]*x + 1")
p3 = toPoly.pm("p0() * (x^4 + 5*x^3 + 5*x + 1)")
pR = factorize.ext.p(p3, by = "x", asBigNum=FALSE, debug=F)
stopifnot( ! is.null(pR[[1]]$GCD))
stopifnot( max(pR[[1]]$GCD$x) == 6)
print(pR[[1]]$GCD)


###
b = c(2,3,5)
plst = lapply(b, function(b) toPoly.pm("x^2 + b[1]*x + 1"))
# Symmetric Polynomial:
p = mult.lpm(plst)

# TODO


###############
### Experiments
p = rPoly(10)

### using BigNumbers
source("Polynomials.Helper.BigNumbers.R")
# Note: the code becomes usually incompatible after switching to Bigz; 
p = toBigz.pm(p);
factorize.ext.p(p, xn="x", asBigNum=TRUE, debug=T)


################
################

### Transform Order 3:
p = toPoly.pm("(x^3 - 3*x + 5)*(x^3 + 4*x + 3)")
r = roots.pm(p)
r

# Note:
# - can be computed directly from the coefficients;
p2 = round(poly.calc(r^3))
p2 = as.pm.polynom(p2)
p2

# r^3 => {3*r - 5, -4*r - 3};
# TODO: explore methods which can benefit from the transformation;
# - still NO factors!
p2 = as.bigz.pm(p2)
factorize.p(p2, "x", asBigNum=TRUE, file=NULL, debug=TRUE)

### only test
p2 = round(poly.calc(r^3 + 5))
p2 = as.pm.polynom(p2)
p2

p2 = as.bigz.pm(p2)
# Note: scaling of "3" is NOT known a priori;
gcd.exact.p(rescale.pm(p, 3, div=TRUE), p2, "x")


#########################
#########################

### Symmetric Polynomials

# Note;
# - all these polynomials are easily decomposable/solvable, see:
#   Polynomials.Derived.P6.Symmetric.R;


### Subtypes of Symmetric Polynomials:

### Subtype Ht:
# b1 + b2 = B[1]
# b1*b2 = B[2] - B[1]
# where: B[id] = coefficient of P[6] polynomial;
# B[3] confirms / rejects sub-type;
factorize.V1P6.F3F3 = function(p) {
	# TODO: check (Partial) Symmetry;
	S  = coef.pm(p, pow=1);
	E2 = coef.pm(p, pow=2) - S;
	# Check:
	b3 = coef.pm(p, pow=3);
	if((b3 - S^2 + 2*E2 - 2) != 0) return(list(isF=FALSE));
	# Solve:
	bd = sqrt(S^2 - 4*E2 + 0i);
	b1 = (S + bd)/2; b2 = (S - bd)/2;
	sol = c(b1=b1, b2=b2);
	return(list(isF=TRUE, F=sol));
}

p = toPoly.pm("(x^3 + b1*x^2 + b2*x + 1) * (x^3 + b2*x^2 + b1*x + 1)")
print.coeff(p)

p = toPoly.pm("(x^2 + b1*x + 1) * (x^2 + b2*x + 1) * (x^2 + b3*x + 1)")
print.coeff(p)

# not symmetric:
# (but easy factorizable)
p = toPoly.pm("(x^3 + b1^2*x^2 + b2*x + 1) * (x^3 + b2^2*x^2 + b1*x + 1)")
print.coeff(p)

# not symmetric:
p = toPoly.pm("(x^3 + b1*x + b2) * (x^3 + b2*x + b1)")
print.coeff(p)

# not symmetric:
p = toPoly.pm("(x^2 + b1*x + b2) * (x^2 + b2*x + b3) * (x^2 + b3*x + b1)")
print.coeff(p)


#############

b = c(3,4,-5)
### Fully Symmetric
p = toPoly.pm("(x^3 + b[1]*x^2 + b[1]*x + 1) * (x^3 + b[2]*x^2 + b[2]*x + 1)")
m = multiplicity.pm(p, -1)
print.pm(p)
checkVal.pm(m, 4)
div.pm(p, toPoly.pm("(x+1)^4"), "x")

### Hidden Ht:
p = toPoly.pm("(x^3 + b[2]*x^2 + b[1]*x + 1) * (x^3 + b[1]*x^2 + b[2]*x + 1)")
m = multiplicity.pm(p, -1)
print.pm(p)
checkVal.pm(m, 0)
# Factorize:
factorize.V1P6.F3F3(p)


### P2-Based:
p = toPoly.pm("(x^2 + b[1]*x + 1)*(x^2 + b[2]*x + 1)*(x^2 + b[3]*x + 1)")
m = multiplicity.pm(p, -1)
print.pm(p)
checkVal.pm(m, 0)
# TODO: factorize


################
################

### Partially-Symmetric

p = toPoly.pm("(x^3+b1*x^2+b2*x+1)*(x^3+b2*x^2+b3*x+1)*(x^3+b3*x^2+b1*x+1)")
print.coeff(p, "x")

### Solution:
# S = B[1]
# E2 = B[2] - B[1]
# E3 = B[3] - S^2 + E2 - 3;
### Consistency Check:
# - Check 1: B[4] + B[5];
#   B[4] + B[5] - 4*S - 2*E2 - E2*S + 3*E3; # == 0!
# - Check 2:
#   (B[4] - 2*S - E2)*(B[5] - 2*S - E2) - (E2^3 + 3*E3^2 - 3*E3*E2*S) - E3*(S^3 - 3*E2*S + 3*E3) - 3*E3^2;

### Ex 1:
pR = replace.pm(p, list(b1=-1, b2=3, b3=4));
pL = factorize.V1P9.QuasiSym(pR)
pL
print.pm(pR)
print.pm(prod.pm(pL))


### Ex 2:
pR = replace.pm(p, list(b1=-1, b2=-3, b3=5));
factorize.V1P9.QuasiSym(pR)

### Ex 3:
pR = replace.pm(p, list(b1=1/2, b2=4, b3=5/2));
factorize.V1P9.QuasiSym(pR)

### Ex 4: S = 0
pR = replace.pm(p, list(b1=2, b2=5, b3=-7));
pL = factorize.V1P9.QuasiSym(pR)
pL
print.pm(pR)
print.pm(prod.pm(pL))

### Ex 5: Fully Symmetric
pR = toPoly.pm("x^9 + 3*x^5 + 3*x^4 + 1");
factorize.V1P9.QuasiSym(pR)



################
################

### P[3] * P[3]

# (x^3 + b2*x^2 + b1*x + 1) * (x^3 + c2*x^2 + c1*x + 1)
p = toPoly.pm("(x^3 + b2*x^2 + b1*x + 1) * (x^3 + c2*x^2 + c1*x + 1)");

b = c(3, -4); c = c(2, 5);
b1 = b[1]; b2 = b[2]; c1 = c[1]; c2 = c[2];
P = function(x) eval.pm(p, list(x=x, b1=b1, b2=b2, c1=c1, c2=c2));

# TODO


################
################

### Experimental

eval.pm(p2, 4) %% 35
eval.pm(p2, 9) %% 35
