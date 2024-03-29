########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### DE Systems: Lambert-type
###
### draft v.0.1h


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
#   => Non-trivial solution: ???

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
# - Solution to original system: see below,
#   Section Supplementary Info;

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


### Non-Trivial: Numeric

### Solver Tools
source("Polynomials.Helper.Solvers.Num.R")


solve.SExp = function(x, R, bb=1) {
	x = matrix(x, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
	y = exp(x) - bb*x[c(2,1)] - R;
	y = rbind(Re(y), Im(y));
	return(y);
}

### Examples
R = 2;
b = 1;
x0 = rbind(
	c(1.8881+5.3416i, 1.8881-5.3416i),
	c(2.5041+11.3727i, 2.5041-11.3727i),
	c(2.000+1.0893i, 1.4221+6.5491i),
	c(2.000-1.0893i, 1.4221-6.5491i),
	c(2.9171+4.9272i, 1.9403-18.0631i),
	c(2.9171-4.9272i, 1.9403+18.0631i),
	c(2.9044+11.245i, 2.507-17.69i),
	c(2.9044-11.245i, 2.507+17.69i)
)
x.all = solve.all(solve.SExp, x0, R=R, bb=b, debug=T)

# Test:
exp(x.all) - b*x.all[, c(2,1)]


### Numerical Approach: Starting Solution
# - just an experiment:
#   but we get another solution;
b = 1
x0 = 1i*sqrt(2)*c(-1, -6); # alternative: c(-1, -10)
R0 = exp(x0) - b*x0[c(2,1)]
R  = 2
path = expand.path(R0, R)
x.all = solve.path(solve.SExp, x0, path=path, bb=b)

### Analyse Path
# - a little bit boring;
tmp = plot.path(c(7,7), c(2,2), x0=x0, solve.SExp, bb=1, steps=37)
# - slightly more going on;
tmp = plot.path(c(-7,-7), c(2,2), x0=x0, solve.SExp, bb=1, steps=101)


### Experimental:
### Case 2.a:
# (x, y) = conjugate roots;
# x = a + d*1i; y = a - d*1i;

id = 1;
a = Re(x.all[id,1]); d = Im(x.all[id,1]);
# Sum =>
exp(a)*cos(d) - b*a - R # = 0
# Diff =>
exp(a)*sin(d) + b*d # = 0

### Sum(squares(...)) =>
exp(2*a) - (b*a + R)^2 - b^2*d^2 # = 0


# Square(Base-Eqs) =>
exp(2*a)*cos(2*d) - R^2 - 2*a*b*R - b^2*(a^2 - d^2) # = 0
# Prod =>
exp(2*a)*sin(2*d) + 2*b*d*(b*a + R) # = 0
# =>
exp(4*a) - (R^2 + 2*a*b*R + b^2*(a^2 - d^2))^2 - 4*b^2*d^2*(b*a + R)^2 # = 0
# unfortunately redundant;


### Case 2.b:
# (x1, y1), (x2, y2) = conjugate set of roots;
# x1 = a1 + d1*1i; y1 = a2 + d2*1i;
# x2 = a1 - d1*1i; y2 = a2 - d2*1i;
# => 4 coupled equations;

# TODO


################
### Non-Trivial: "Brute-Force"

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


#######################

### Extension

# exp(x*y1) = b1*y2 + R1
# exp(x*y2) = b2*y1 + R2

# Note:
# - "slightly easier" to solve when:
#   b1 = b2, R1 = R2;
#   => trivial solution: y1 = y2;
#   => Non-trivial solution: ???

### D =>
exp(x*y1)*(x*dy1 + y1) - b1*dy2 - db1*y2 - dR1 # = 0
exp(x*y2)*(x*dy2 + y2) - b2*dy1 - db2*y1 - dR2 # = 0
# =>
(b1*y2 + R1)*(x*dy1 + y1) - b1*dy2 - db1*y2 - dR1 # = 0
(b2*y1 + R2)*(x*dy2 + y2) - b2*dy1 - db2*y1 - dR2 # = 0

### Special Cases:

### b1 = b2; R1 = R2;
b*x*y2*dy1 + x*R*dy1 - b*dy2 + b*y1*y2 + R*y1 - db*y2 - dR # = 0
b*x*y1*dy2 + x*R*dy2 - b*dy1 + b*y1*y2 + R*y2 - db*y1 - dR # = 0
### b = const;
b*x*y2*dy1 + x*R*dy1 - b*dy2 + b*y1*y2 + R*y1 - dR # = 0
b*x*y1*dy2 + x*R*dy2 - b*dy1 + b*y1*y2 + R*y2 - dR # = 0

### TODO: check;


#########################
#########################

### S2: Both Exponentials

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


###################
###################

### Inv-Trig

# x*atan(y1) = b1*y2 + R1
# x*atan(y2) = b2*y1 + R2

### D =>
(y1^2 + 1)*atan(y1) + x*dy1 - (y1^2 + 1)*(b1*dy2 + db1*y2 + dR1) # = 0
(y2^2 + 1)*atan(y2) + x*dy2 - (y2^2 + 1)*(b2*dy1 + db2*y1 + dR2) # = 0
# Subst =>
(y1^2 + 1)*(b1*y2 + R1) + x^2*dy1 - x*(y1^2 + 1)*(b1*dy2 + db1*y2 + dR1) # = 0
(y2^2 + 1)*(b2*y1 + R2) + x^2*dy2 - x*(y2^2 + 1)*(b2*dy1 + db2*y1 + dR2) # = 0
# =>
x^2*dy1 - x*b1*(y1^2 + 1)*dy2 + (y1^2 + 1)*(b1*y2 - x*db1*y2 + R1 - x*dR1) # = 0
x^2*dy2 - x*b2*(y2^2 + 1)*dy1 + (y2^2 + 1)*(b2*y1 - x*db2*y1 + R2 - x*dR2) # = 0


### Special Cases:

# 1.) b1 = b2 = b; R1 = R2 = R;
x^2*dy1 - x*b*(y1^2 + 1)*dy2 + (y1^2 + 1)*(b*y2 - x*db*y2 + R - x*dR) # = 0
x^2*dy2 - x*b*(y2^2 + 1)*dy1 + (y2^2 + 1)*(b*y1 - x*db*y1 + R - x*dR) # = 0
# b = ct =>
x^2*dy1 - b*x*(y1^2 + 1)*dy2 + (y1^2 + 1)*(b*y2 + R - x*dR) # = 0
x^2*dy2 - b*x*(y2^2 + 1)*dy1 + (y2^2 + 1)*(b*y1 + R - x*dR) # = 0

# 2.) b1 = b2 = x^n; R1 = R2 = x^m;
x^2*dy1 - x^(n+1)*(y1^2 + 1)*dy2 - (y1^2 + 1)*((n-1)*x^n*y2 + (m-1)*x^m) # = 0
x^2*dy2 - x^(n+1)*(y2^2 + 1)*dy1 - (y2^2 + 1)*((n-1)*x^n*y1 + (m-1)*x^m) # = 0

### Example:
# Special case [2] with n = m = 2;
dy1 - x*(y1^2 + 1)*dy2 - (y1^2 + 1)*(y2 + 1) # = 0
dy2 - x*(y2^2 + 1)*dy1 - (y2^2 + 1)*(y1 + 1) # = 0

### TODO: check;

