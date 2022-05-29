########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### DE Systems: Lambert-type
###
### draft v.0.1f


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

### Initial system:
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
y = log(- b * lambertWp( - exp(-R/b) / b));

### Case 2:
# - non-Trivial: y1 != y2;


source("Polynomials.Helper.ODE.R")
source("Polynomials.Helper.Solvers.S2.R")

decompose.S2Exp = function(n, pLinear, xn=c("y1", "y2"), asDiv=TRUE) {
	pExp = expand.Exp(n = n, xn=xn[[1]], asDiv = asDiv);
	pLin = mult.pm(pLinear, pExp$Div);
	pExp = diff.pm(pExp$P, pLin);
	#
	pR = decompose.S2Ht(pExp, vars=xn);
	pR$pDiff$D = NULL;
	return(pR);
}
solve.pm.S2Basic = function(sol, debug=TRUE) {
	S   = roots(sol$Rez$coeff);
	if(debug) print(S);
	xy  = sapply(S, function(S) eval.pm(sol$x0, c(S=S)));
	div = sapply(S, function(S) eval.pm(sol$div, c(S=S))); print(xy)
	xy = xy / div;
	xy.diff = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.diff) / 2;
	y = (S - xy.diff) / 2;
	sol = cbind(x=x, y=y);
	return(sol);
}

### Test

n = 6 # Number of Terms
pLin = toPoly.pm("b*y2 + R")
pR0  = decompose.S2Exp(n=n, pLin)
str(pR0)

### b = ...
b = 1; R = 1/2
p1 = replace.pm(pR0$pDiff, c(b=b))
p2 = replace.pm(pR0$pSum, c(b=b, R=R))
# p1 = pR0$pDiff; p2 = pR0$pSum;
pR = solve.pm(p1, p2, xn="E2")
pR$Rez$coeff = - pR$Rez$coeff;
if(n == 6) {
	pR$Rez = div.pm(pR$Rez, toPoly.pm("S^2 + 4*S + 4"), "S")$Rez
}
pR$Rez = sort.pm(pR$Rez, "S")
str(pR)
# NO overflow yet!
max(abs(pR$Rez$coeff))

### Solve: explicitly
sol = solve.pm.S2Basic(pR)

# TODO: debug!
# - no valid solutions;

### Test
exp(sol[,1]) - b*sol[,2]
exp(sol[,2]) - b*sol[,1]

id = 5
eval.pm(p1, c(S=sum(sol[id,]), E2=prod(sol[id,])))
eval.pm(p2, c(S=sum(sol[id,]), E2=prod(sol[id,])))

### Debug:
# for n = 6:
poly.calc(sapply(14:15, function(id) sum(sol[id, ])))

b = 1; R = 1/2;
x1 = 1.9398591713 + 1.3497579982i;
x2 = 1.0254412580 + 6.7884907688i;
exp(x1) - b*x2; exp(x2) - b*x1;

eval.pm(pR0$pDiff, c(S=x1+x2, b=1, E2=x1*x2))
eval.pm(expand.Exp(n = 6, xn="x1", asDiv = FALSE), c(x1=x1))

eval.pm(expand.Exp(n = 6, xn="x1", asDiv = FALSE)$P, c(x1=x1))
exp(x1)
# necessitates > 15 terms in the expansion:
eval.pm(expand.Exp(n = 15, xn="x1", asDiv = FALSE)$P, c(x1=x2))
exp(x2)

### Diff =>
exp(x1) + b*x1 - (exp(x2) + b*x2) # = 0
exp(1/b*exp(x1)) * exp(x1) - exp(1/b*exp(x2)) * exp(x2) # = 0
# Lambert W: Wp or Wn =>
exp(x1) - exp(x2) # = 0 *OR*
# different branches of W:
lambertWp(exp(1/b*exp(x1)) * 1/b*exp(x1))
# [but pracma-implementation does NOT accept complex numbers]
# print(exp(1/b*exp(x1)) * 1/b*exp(x1), 12)


### Case 1:
# - Trivial solution: y1 = y2;
# Note: R is actually a function of x;

R = 3
b = 1.2
y = log(- b * lambertWp( - exp(-R/b) / b));

### Test
exp(y) - b*y - R;


##########################
##########################

### Initial Exp-System

# exp(y1) + exp(y2) = R1
# y1*exp(y1) + y2*exp(y2) = R2

# TODO


######################

### Initial Exp-System

# y1*exp(y1) + y2*exp(y2) = R1
# y2*exp(y1) + y1*exp(y2) = R2

# Note:
# - easier to solve when:
#   R1 = R2;
#   => trivial solution: y1 = y2;

### ODE System:

### Solve Linear System:
# exp(y1) = (R1*y1 - R2*y2) / (y1^2 - y2^2)
# exp(y2) = (R2*y1 - R1*y2) / (y1^2 - y2^2)

### D =>
(y1 + 1)*exp(y1)*dy1 + (y2 + 1)*exp(y2)*dy2 - dR1 # = 0
(y2*dy1 + dy2)*exp(y1) + (y1*dy2 + dy1)*exp(y2) - dR2 # = 0
# =>
(y1 + 1)*(R1*y1 - R2*y2)*dy1 + (y2 + 1)*(R2*y1 - R1*y2)*dy2 - dR1*(y1^2 - y2^2) # = 0
(y2*dy1 + dy2)*(R1*y1 - R2*y2) + (y1*dy2 + dy1)*(R2*y1 - R1*y2) - dR2*(y1^2 - y2^2) # = 0

### System:
(y1 + 1)*(R1*y1 - R2*y2)*dy1 + (y2 + 1)*(R2*y1 - R1*y2)*dy2 - dR1*(y1^2 - y2^2) # = 0
(R1*y1*y2 - R2*y2^2 + R2*y1 - R1*y2)*dy1 +
	+ (R2*y1^2 - R1*y1*y2 + R1*y1 - R2*y2)*dy2 - dR2*(y1^2 - y2^2) # = 0


### Special Cases:
# R1 = R2 = R
R*(y1 + 1)*(y1 - y2)*dy1 + R*(y2 + 1)*(y1 - y2)*dy2 - dR*(y1^2 - y2^2) # = 0
R*(y1*y2 - y2^2 + y1 - y2)*dy1 + R*(y1^2 - y1*y2 + y1 - y2)*dy2 - dR*(y1^2 - y2^2) # = 0
# Diff (Eq 2 - Eq 1) =>
- R*(y1 - y2)^2*dy1 + R*(y1 - y2)^2*dy2 # = 0
# =>
dy1 - dy2 # = 0
# [only solution]


# TODO: check;

