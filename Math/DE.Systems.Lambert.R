########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### DE Systems: Lambert-type
###
### draft v.0.1b


#############
### Types ###
#############

### Symmetric:
# - TODO
### Hetero-Symmetric:
# - Simple, derived from:
#   exp(y1) = b1*y2 + R1;


####################

### Helper functions

library(pracma)
# needed for Lambert W;


# include: DE.ODE.Helper.R;
source("DE.ODE.Helper.R")


#########################
#########################

### Hetero:
### Symmetric & non-Symmetric

# exp(y1) = b1*y2 + R1
# exp(y2) = b2*y1 + R2

# Note:
# - "slightly easier" to solve when:
#   b1 = b2, R1 = R2;
#   => trivial solution: y1 = y2;

### D =>
exp(y1)*dy1 - b1*dy2 - db1*y2 - dR1 # = 0
exp(y2)*dy2 - b2*dy1 - db2*y1 - dR2 # = 0

# DE System:
(b1*y2 + R1)*dy1 - b1*dy2 - db1*y2 - dR1 # = 0
(b2*y1 + R2)*dy2 - b2*dy1 - db2*y1 - dR2 # = 0

### Special Cases:
# b1 = b2 = ct;
(b*y2 + R1)*dy1 - b*dy2 - dR1 # = 0
(b*y1 + R2)*dy2 - b*dy1 - dR2 # = 0

# TODO: check & solve;


#############

### Variants:

# exp(y1+y2) = b1*(y1 - y2) + R1
# exp(y1-y2) = b2*(y1 + y2) + R2

# Note:
# - "slightly easier" to solve when:
#   b1 = b2, R1 = R2;
#   => trivial solution: y2 = 0;

### D =>
exp(y1+y2)*(dy1 + dy2) - b1*dy1 - db1*y1 + b1*dy2 + db1*y2 - dR1 # = 0
exp(y1-y2)*(dy1 - dy2) - b2*dy1 - db2*y1 - b2*dy2 - db2*y2 - dR2 # = 0

### DE System:
(b1*y1 - b1*y2 + R1)*(dy1 + dy2) - b1*dy1 - db1*y1 + b1*dy2 + db1*y2 - dR1 # = 0
(b2*y1 + b2*y2 + R2)*(dy1 - dy2) - b2*dy1 - db2*y1 - b2*dy2 - db2*y2 - dR2 # = 0

### Transforms:
((b1+b2)*y1 - (b1-b2)*y2 + R1 + R2)*dy1 +
	+ ((b1-b2)*y1 - (b1+b2)*y2 + R1 - R2)*dy2 +
	- (b1+b2)*dy1 - (db1+db2)*y1 + (b1-b2)*dy2 + (db1-db2)*y2 - dR1 - dR2 # = 0
((b1-b2)*y1 - (b1+b2)*y2 + R1 - R2)*dy1 +
	+ ((b1+b2)*y1 - (b1-b2)*y2 + R1 + R2)*dy2 +
	- (b1-b2)*dy1 - (db1-db2)*y1 + (b1+b2)*dy2 + (db1+db2)*y2 - dR1 + dR2 # = 0
# bs = b1 + b2; bd = b1 - b2; =>
(bs*y1 - bd*y2 + Rs)*dy1 + (bd*y1 - bs*y2 + Rd)*dy2 +
	- bs*dy1 - dbs*y1 + bd*dy2 + dbd*y2 - dRs # = 0
(bd*y1 - bs*y2 + Rd)*dy1 + (bs*y1 - bd*y2 + Rs)*dy2 +
	- bd*dy1 - dbd*y1 + bs*dy2 + dbs*y2 - dRd # = 0

### Special Case:
# b1 = b2 = ct; bd = 0; dbs = dbd = 0;
  (bs*y1 + Rs - bs)*dy1 - (bs*y2 - Rd)*dy2 - dRs # = 0
- (bs*y2 - Rd)*dy1 + (bs*y1 + Rs + bs)*dy2 - dRd # = 0


#######################
#######################

### Supplementary Info:

### Solver:
# exp(y1) = b*y2 + R;
# exp(y2) = b*y1 + R;

### Case 1:
# - Trivial solution: y1 = y2;

### Case 2:
# - non-Trivial: y1 != y2;

source("Polynomials.Helper.ODE.R")
source("Polynomials.Helper.Solvers.S2.R")

decompose.S2Exp = function(n, pLinear, xn=c("y1", "y2")) {
	pExp = expand.Exp(n = n, xn=xn[[1]], asDiv = TRUE)
	pLin = mult.pm(pLinear, pExp$Div);
	pExp = diff.pm(pExp$P, pLin);
	#
	pR = decompose.S2Ht(pExp, vars=xn);
	pR$pDiff$D = NULL;
	return(pR);
}

### Test

n = 5 # Number of Terms
pLin = toPoly.pm("b*y2 + R")
pR0  = decompose.S2Exp(n=n, pLin)
str(pR0)

### b = ...
b = 1; R = 1
p1 = replace.pm(pR0$pDiff, c(b=b))
p2 = replace.pm(pR0$pSum, c(b=b, R=R))
pR = solve.pm(p1, p2, xn="E2")
str(pR)
# NO overflow yet!
max(abs(pR$Rez$coeff))

# TODO

