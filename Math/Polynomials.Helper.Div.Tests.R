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
cat("Example 2: Divisible w Div Factor\n")
p1   = toPoly.pm("x^3 + y^3")
pDiv = toPoly.pm("a*x + a*y")
pR = div.pm(p1, pDiv, by="x", NF.stop=FALSE)
pR

if(nrow(pR$Rem) != 0) stop("Division went wrong! Should be Divisible!");
tmp = diff.pm(p1*pR$pDiv, pDiv * as.pm(pR$Rez))
if(nrow(tmp) != 0) stop("Division went wrong!")


#################

cat("\nDiv Tests: Finished!\n")
cat("Warning: Very few Tests!\n")

