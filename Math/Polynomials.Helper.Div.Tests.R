########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Multi-Variable Polynomials
### Polynomial Division: Tests
###
### draft v.0.1a


### Tests:
### Division in Multi-Variable Polynomials


######################

### Helper Functions

### fast load:
# source("Polynomials.Helper.Div.Tests.R")


source("Polynomials.Helper.R")
# - is automatically loaded in: Polynomials.Helper.R;
# source("Polynomials.Helper.Div.R")


########################
########################


###
cat("\nSection 1: Simple\n\n")

###
p = data.frame(coeff=c(1, -1), x=c(3,0), y=c(0,3));
pR = div.pm(p, "x-y", "x");
if(nrow(pR$Rem) != 0) stop("Division went wrong!");


cat("\nFinished: Section 1!\n")

#####################

cat("\nSection 2: Non-compliant Leads\n\n")

###
cat("Example 1: Intermediary Lead w 2 Rows\n")
b1 = 3;
p1   = toPoly.pm("x^3 + b1()*x + 1")
pDiv = toPoly.pm("a*x + 1")
pR = div.pm(p1, pDiv, by="x", NF.stop=FALSE)
pR

if(nrow(pR$Rem) != 3) stop("Division went wrong! Should be non-divisible!");
tmp = diff.pm(p1*pR$pDiv, pDiv * as.pm(pR$Rez))
tmp = diff.pm(tmp, pR$Rem);
if(nrow(tmp) != 0) stop("Division went wrong!");


###
cat("\nExample 2: Divisible w Div Factor\n")
p1   = toPoly.pm("x^3 + y^3")
pDiv = toPoly.pm("a*x + a*y")
pR = div.pm(p1, pDiv, by="x", NF.stop=FALSE)
pR

if(nrow(pR$Rem) != 0) stop("Division went wrong! Should be Divisible!");
tmp = diff.pm(p1*pR$pDiv, pDiv * as.pm(pR$Rez))
if(nrow(tmp) != 0) stop("Division went wrong!")


cat("\nFinished: Section 2!\n")

#####################

cat("\nSection 3: GCD\n\n")


#####################

cat("\nSection 4: Multivariable GCD\n\n")

### Test 1:
cat("Test 1: Simple Multi-Var GCD\n")
p1 = toPoly.pm("a*x^3 + a*y^3")
p2 = toPoly.pm("2*b*x + 2*b*y")

pR = gcd.pm.exact(p1, p2, "x");

if(nrow(pR) != 2) stop("Wrong GCD");


### Test 2:
cat("\nTest 2: Multi-Var GCD\n")
p1 = toPoly.pm("x^3 + y^3")
p2 = toPoly.pm("x^2 - y^2")

pR = gcd.pm.exact(p1, p2, "x");

if(nrow(pR) != 2) stop("Wrong GCD");


### Test 3:
cat("\nTest 3: Multi-Var GCD\n")
p1 = toPoly.pm("x^5 + y^5")
p2 = toPoly.pm("x^3 + y^3")

pR = gcd.pm.exact(p1, p2, "x");

if(nrow(pR) != 2) stop("Wrong GCD");


### Test 4:
cat("\nTest 4: Multi-Var GCD\n")
p1 = toPoly.pm("(x^2 + y^2)*(x^3 + y^3 + 2)")
p2 = toPoly.pm("(x + y)*(x^2 + y^2)")

pR = gcd.pm.exact(p1, p2, "x");

if(nRow(pR) != 2) stop("Wrong GCD");


#############
# Challenges:

cat("\nChallanges:\n")

### Test 1:
# - nice factorization of P1(x) * P2(x, y);
cat("\nTest 3: P1(x) * P2(x, y)\n")
p0 = toPoly.pm("x^3 + 3*x + 1")
p1 = toPoly.pm("p0()*(x^2 + 2*x*y + 5*y^2)")
p2 = dp.pm(p1, "y")

pR = gcd.pm.exact(p1, p2, "x");
pR$coeff = pR$coeff * sign(top.pm(pR, "x")$coeff)

if(nrow(pR) != 3) stop("Wrong GCD");
if(nrow(diff.pm(pR, p0)) != 0) stop("Wrong GCD");


### Test 2:
# - TODO: automatic factorization of P1(x) * P2(x, y);
# - works with higher derivatives;
cat("\nTest 3: P1(x) * P2(x, y)\n")
p0 = toPoly.pm("x^4 - 5*x^2 + 3*x + 1")
p1 = toPoly.pm("p0()*(x^2 + 2*x*y + 5*y^2 - x*y^2 + y^3)")
p2 = dp.pm(p1, "y")

pR = gcd.pm.exact(p1, p2, "x");
pR$coeff = pR$coeff * sign(top.pm(pR, "x")$coeff)

if(max(pR$x) != 4) stop("Wrong GCD");
# if(nrow(pR) != 4) stop("Wrong GCD");
# if(nrow(diff.pm(pR, p0)) != 0) stop("Wrong GCD");

# higher D():
pR = gcd.pm.exact(p1, dnp.pm(p1, n=3, "y"), "x");
pR$coeff = pR$coeff * sign(top.pm(pR, "x")$coeff)
if(nrow(pR) != 4) stop("Wrong GCD");
if(nrow(diff.pm(pR, p0)) != 0) stop("Wrong GCD");


### Test 3:
# - more complex factorization of type P1(x) * P2(x, y);
cat("\nTest 3: P1(x) * P2(x, y)\n")
p1 = toPoly.pm("(x^2 + 2*y^2)*(x^2 + 3*z^2)*(y^2 + 4*z^2)")

p2 = dp.pm(p1, "y")
pR = gcd.pm.exact(p1, p2, "x");

p2 = dp.pm(p1, "z")
pR = gcd.pm.exact(p1, p2, "y");

p2 = dp.pm(p1, "x")
pR = gcd.pm.exact(p1, p2, "z");

# TODO


### Test 4:
cat("\nTest 4: Multi-Var GCD\n")
p1 = toPoly.pm("x^5 + y^5 + 2*(x + y)")
p2 = toPoly.pm("x^3 + y^3")

pR = gcd.pm.exact(p1, p2, "x");

# if(nrow(pR) != 2) stop("Wrong GCD");


###########
### Test 5:
cat("\nTest 5: Multi-Var GCD\n")
p1 = toPoly.pm("a*x^3 + a*y^3")
p2 = toPoly.pm("2*b*(x + y)*(b*x - 3*y)")

pR = gcd.pm.exact(p1, p2, "x");

# if(nrow(pR) != 2) stop("Wrong GCD");


###########
### Test 6:
cat("\nTest 6: Multi-Var GCD\n")
p1 = toPoly.pm("(x^2 + y^2)*(x^3 + y^3 + 3*x*y + 2)")
p2 = toPoly.pm("(x + y + 1)*(x^2 + y^2)")

pR = gcd.pm.exact(p1, p2, "x");

# if(nrow(pR) != 2) stop("Wrong GCD");


###########
### Test 7:
cat("\nTest 7: Multi-Var GCD\n")
p1 = toPoly.pm("(x^2 + a*y^2)^2 + b*x*(x^2 + a*y^2)")
p2 = toPoly.pm("(x^2 + a*y^2)*(x^2 - c*y)")

pR = gcd.pm.exact(p1, p2, "x");

# finds Original: (x^2 + a*y^2)*(x^2 - c*y)
pR = gcd.pm.exact(p1, p2, c("x","a"));

# if(nrow(pR) != 2) stop("Wrong GCD");


###########
### Test 8:
cat("\nTest 8: Multi-Var GCD\n")
p1 = toPoly.pm("(x^3 + a*y^2)^2 + b*x*(x^3 + a*y^2)")
p2 = toPoly.pm("(x^3 + a*y^2)*(x^2 - c*y)")

# Maximum Power in Factor:
pR = gcd.pm.exact(p1, p2, "x");
pR = gcd.pm.exact(p1, p2, "y");
pR = gcd.pm.exact(p1, p2, "a");

# if(nrow(pR) != 2) stop("Wrong GCD");


###########
### Test 9:
cat("\nTest 9: Multi-Var GCD\n")
p1 = toPoly.pm("(x^3*y + a*y^2 + 1)^2 + b*c^2*x*(x^3*y + a*y^2 + 1)")
p2 = toPoly.pm("(x^3*y + a*y^2 + 1)*(x^2 - c*y)")

# Maximum Power in Factor:
pR = gcd.pm.exact(p1, p2, "x");
pR = gcd.pm.exact(p1, p2, "y");
pR = gcd.pm.exact(p1, p2, "a");

# if(nrow(pR) != 2) stop("Wrong GCD");


cat("\nFinished: Section 4!\n")


#################

cat("\nDiv Tests: Finished!\n")
cat("Warning: Very few Tests!\n")

