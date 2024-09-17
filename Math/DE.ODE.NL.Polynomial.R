########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs: Polynomial types
###
### draft v.0.4h


### TODO:
# - clean-up/reorder sections;


###############
### History ###
###############

### Order 1 Non-Linear

### draft v.0.4g - h:
# - moved Section: "Derived from other ODEs"
#   to file: DE.ODE.FromODEs.R;
# - moved Simple Radicals to file:
#   DE.ODE.NL.Poly.Radicals.R;
# - moved ODEs derived from Integrals to file:
#   DE.ODE.NL.Poly.Int.R;
### draft v.0.4f - v.0.4f-exp:
# - derived from other ODEs:
#   x*y*d2y - 1/n*x*dy^2 - y*dy + 1/n*x^2*dy = 0;
# - some experiments with simple base-ODE;
### draft v.0.4e - v.0.4e-fix2:
# - Example based on radicals of y:
#   (n-1)*x*y*dy + x^3*dy - 2*n*y^2 = 0;
### draft v.0.4d:
# - Example based on P[5]:
#   y^3*dy - 2*b0*y*dy + 5*x*dy - y = 0;
### draft v.0.4c:
# - improved formatting; [v.0.4c & v.0.4d-bis]
### draft v.0.4a - v.0.4b:
# - derived from:
#   y = x * I(1 / (x^n + k) dx) + F0(x); [fixed: v.0.4a-fix]
#   y * I(1 / (x^n + k) dx) = F0(x); [v.0.4b]
### draft v.0.3b:
# - Base-solution with Radicals:
#   2*x*y^2*d2y - (n-2)*y^2*dy - n*dy = 0;
### draft v.0.3a - v.0.3a-gen: [16-11-2020]
# - Technique: Integration by parts:
#   c/2*dz^2 + (3*x+b)*dz - 4*z = 0;
# - slight generalization & more examples; [v.0.3a-gen]
### draft v.0.2e: [11-11-2020]
# - solved: x*(x^n - 1)*d2y + x*(x^n - 1)*y*dy + n*dy - b0*x*y = -n*b0;
#   where b0 = constant; [nice mix]
### draft v.0.2d: [11-11-2020]
# - cleanup: moved section with Trigonometric functions
#   to DE.ODE.Trigonometric.R;
### draft v.0.2c: [06-11-2020]
# - solved:
#   x^2*(x^3-1)*dz - c^2*z^3 - c*x^2*z^2 + 2*x*z = 0;
#   [where c = constant]
### draft v.0.2b:
# - improved code;
### draft v.0.2a: [01-11-2020]
# - re-organization of sections;
# - brief theory of polynomial entangling:
#   the introduction of 1 degree of freedom;
### draft v.0.2-pre-a: [31-10-2020]
# - moved section Exponentials/Lambert
#   to separate file: DE.ODE.Fractions.Lambert.R;
### draft v.0.1f - v.0.1f-2: [30-10-2020]
#   [moved now to separate file]
# - Lambert-type equations (not simple polynomial)!
# - solved: x*y*dy + x*(x+1)*dy - 2*y = x;
#   (a Lambert snack)
# - various generalizations, e.g.:
#   (x^2+b)*y*dy + (x^2+b)*(x+1)*dy - 2*x*y = x^2 - b; (v.0.1f-2)

### draft v.0.1e-pre-snack - v.0.1e-snack4:
# - derived: the snack (D(P5));
# - added a simpler snack (snack2): D(P5), where P5 = P((x-1)^(1/5)),
#   and a parametric one (snack3): P5 = P((x-a)^(1/5));
# - basic exploration of (hidden) Trigonometric functions; (snack4/fixed)
### P3 & ODE Transformations:
### draft v.0.1d-tr - v.0.1d-tr2:
# - more polynomial transformations,
#   including nested & combined transformations, e.g.:
#   x*y*dy + x^3*dy - 1/3*x*y^3 - 1/2*y^2 + 2/3*x^4 = 0;
### draft v.0.1d - v.0.1d-dx3:
# - systematic approach to ODE Transformations:
#   -- various substitutions;
#   -- combination of ODE-variants;
# - added order 2 combinations (v.0.1d-dx2 & v.0.1d-dx3);
### draft v.0.1c-tr - v.0.1c-tr4:
# - added transformed base-polynomials:
#   y^3*dy + (x+1)*dy - 3*x*y^2 - y = 0; (v.0.1c-tr)
#   x^2*y*dy + (x+1)*dy - x*y^2 - 1/3*y = 0; (v.0.1c-tr2)
# - transformation of ODE:
#   ### f1(x)*y*dy - f2(x)*y^2 variants:
#   2*x*y*dy + 2*x^3*dy - y^2 - 2*x^2*y = 0;
#   x^3*y*dy + 3*x*log(x)*dy - x^2*y^2 - y = 0; (v.0.1c-tr3)
#   ### higher powers:
#   x*y^2*dy - x^3*dy - 2/3*y^3 + 4*log(x) - 2 = 0; (v.0.1c-tr4)
#   y^4*dy - 6*log(x)*y*dy - 3*x^4*dy - 6*x^3*y - 6*x = 0; (v.0.1c-tr4)

### draft v.0.1c - v.0.1c-sh2:
# - added symmetrically shifted, eg:
#   y^2*dy - 2*y*dy - (x-1)*dy - y = x^2 - 2;
#   y^2*dy + (x^2 - x)*dy + (2*x-1)*y = 9*x^2 - 4*x; (v.0.1c-sh2)
### draft v.0.1b-sh:
# - added classic/full Cardan Polynomials (P3);
# - added shifted version (P3) (v.0.1b-sh);
### draft v.0.1a-plot - v.0.1a-px:
# - added diagnostic plots (+ tangent lines);
# - added more examples (v.0.1a-px);
### draft v.01a:
# - initial draft:
#   some ODEs based on Cardan polynomials;


###################

###################
### Terminology ###

# let f(x), h(x) be 2 functions which are differentiable;
# y(x) = function to be found;
# dy = d(y(x));
# d2y = d(d(y(x)));

### PDE:
# TODO: move to separate file;
# p(x), q(x) = functions to be found;
# dp, dq:
# dp = d(p(x)); dq = d(q(x));
# Note: q(x) is an intermediary function used for various derivations;


#################
###  Theory   ###
###  (brief)  ###
#################

# y^n + P(y) = 0
# where:
# - P(y) = polynomial in y of order < n;
# - the coefficients of P(y) are functions of x;
#   (both polynomial & non-polynomial functions)

### D() =>
# n*y^(n-1)*dy + D(P) = 0

# - unfortunately, the coefficients of D(P) are correlated;
# - we will de-correlate these coefficients by combining/entangling
#   with the original polynomial;

### Multiply with y =>
n*y^n*dy + y*D(P) # = 0
### Substitute: y^n = - P(y)
n*P*dy - y*D(P) # = 0;

### Note:
# - the coefficients of n*P*dy will combine with
#   the corresponding coefficients of y*D(P);
# - this will de-correlate the coefficients;
# - although the coefficients of D(P) are still correlated,
#   the term n*P*dy introduces 1 degree of freedom;
# - the numeric coefficient n from n*P*dy
#   introduces 1 degree of freedom in the equation;


####################
####################

### Helper Functions

# include: DE.ODE.Helper.R;
source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


##########################

########################
###     Order 1      ###
###   Non-Linear     ###
########################


###########################
### Cardano-Polynomials ###
###########################

### Cardano-Polynomials:
### Full Root

### System:
# y = p + q;
# where:
# p^n + q^n = 2*f(x)
# p*q = h(x)
### Examples of such polynomials:
# - see the examples below;

### Solutions
# y = p + q, where:
# p = (f + sqrt(f^2 - h^n))^(1/n)
# q = (f - sqrt(f^2 - h^n))^(1/n)
# Note:
# - it is possible to rotate these solutions using the roots of unity;
# - transformations of the Polynomial & of the ODE: decouple the derivatives;
#   [see later]


##################
### Parametric ###
###  Examples  ###
##################

### Polynomial:
# Simple P3:
# n = 3
# y^3 - 3*h*y - 2*f = 0

### D =>
3*y^2*dy - 3*h*dy - 3*y*dh - 2*df # = 0

### Primary ODE: correlated
y^2*dy - h*dy - dh*y - 2/3*df # = 0
# - see Transformations for de-correlation;

#######################
### Transformations ###
#######################

### T.A.) ODE-Transformations:

### T.A.1: Substitute y^3
# [not run]
y^2*dy - h*dy - y*dh - 2/3*df = 0 # * y
y^3*dy - h*y*dy - y^2*dh - 2/3*df*y = 0
(3*h*y + 2*f)*dy - h*y*dy - dh*y^2 - 2/3*df*y = 0
2*h*y*dy + 2*f*dy - dh*y^2 - 2/3*df*y = 0

### ODE:
h*y*dy + f*dy - 1/2*dh*y^2 - 1/3*df*y = 0

### variant form: (but not reducible)
# 1/2*D(y^2/h) + (f*dy - 1/3*df*y)/h^2 = 0

### T.A.2: Substitute y
# - can replace one occurrence or multiple procurances;
y^2*dy - h*dy - y*dh - 2/3*df = 0 # *h
h*y^2*dy - h^2*dy - h*y*dh - 2/3*h*df = 0
h*y^2*dy - h^2*dy - 1/3*(y^3 - 2*f)*dh - 2/3*h*df = 0
h*y^2*dy - h^2*dy - 1/3*y^3*dh + 2/3*f*dh - 2/3*h*df = 0
# 1/3*(3/h*y^2*dy - y^3*dh/h^2) - dy - 2/3*D(f/h) = 0
# 1/3*D(y^3/h) - dy - 2/3*D(f/h) = 0

### T.A.3: "Linear" Combinations
h*y^2*dy - h^2*dy - 1/3*y^3*dh + 2/3*f*dh - 2/3*h*df +
 + b*(2*h*y*dy + 2*f*dy - y^2*dh - 2/3*df*y) = 0

### T.A.4: y^3 & y:
h*y*dy + f*dy - 1/2*dh*y^2 - 1/3*df*y = 0 # * h
h^2*y*dy + f*h*dy - 1/2*h*dh*y^2 - 1/3*df*h*y = 0
h^2*y*dy + f*h*dy - 1/2*h*dh*y^2 - 1/9*df*(y^3 - 2*f) = 0
h^2*y*dy + f*h*dy - 1/9*df*y^3 - 1/2*h*dh*y^2 + 2/9*f*df = 0


### T.B.) Poly-Transformations:
# - yield similar/equivalent equations;

### T.B.1: * 1/y
# y^3 - 3*h*y - 2*f = 0 # / * 1/y
# y^2 - 3*h - 2*f/y = 0 # D() =>
2*y*dy - 3*dh - 2*df/y + 2*f*dy/y^2 = 0 # * y^2
y^3*dy + f*dy - 3/2*dh*y^2 - df*y = 0

### T.B.1 + T.A.1: y^3
(3*h*y + 2*f)*dy + f*dy - 3/2*dh*y^2 - df*y = 0
3*h*y*dy + 3*f*dy - 3/2*dh*y^2 - df*y = 0
h*y*dy + f*dy - 1/2*dh*y^2 - 1/3*df*y = 0 # same as [T.A.1]

### T.B.1 + T.A.2: y
# - can replace one occurrence or multiple occurrences;
y^3*dy + f*dy - 3/2*dh*y*(y^3 - 2*f)/(3*h) - df*y = 0
y^3*dy + f*dy - 1/2*dh/h*y^4 + f*dh/h*y - df*y = 0 # * h
h*y^3*dy + f*h*dy - 1/2*dh*y^4 + f*dh*y - h*df*y = 0
# ... - D(h*y/f) = 0

# alternative: T.B.1 + T.A.2: y
y^3*dy + f*dy - 3/2*dh*y^2 - df*(y^3 - 2*f)/(3*h) = 0 # * h
h*y^3*dy + f*h*dy - 1/3*df*y^3 - 3/2*h*dh*y^2 + 2/3*f*df = 0

### T.B.1 + T.A.f: f
# - but is trivial: initial D(P3);
# - may be useful with higher orders;
y^3*dy + 1/2 * (y^3 - 3*h*y)*dy - 3/2*dh*y^2 - df*y = 0 # *2
2*y^3*dy + (y^3 - 3*h*y)*dy - 3*dh*y^2 - 2*df*y = 0
y^3*dy - h*y*dy - dh*y^2 - 2/3*df*y = 0
y^2*dy - h*dy - dh*y - 2/3*df = 0 # initial derivative;

### V.B.1 + Back-Sub: f Back into P3
# - may be useful with D of higher Order;
y^3*dy + f*dy - 3/2*dh*y^2 - df*y = 0
# => f*dy = -y^3*dy + 3/2*dh*y^2 + df*y
# y^3 - 3*h*y - 2*f = 0 # * dy
y^3*dy - 3*h*y*dy - 2*f*dy = 0
y^3*dy - 3*h*y*dy - 2*(-y^3*dy + 3/2*dh*y^2 + df*y) = 0
3*y^3*dy - 3*h*y*dy - 3*dh*y^2 - 2*df*y = 0
3*y^2*dy - 3*h*dy - 3*dh*y - 2*df = 0 # initial D(P3)
# D2: 3*y^2*dy^2 + y^3*d2y + df*dy + f*d2y - 3*dh*y*dy - 3/2*d2h*y^2 - df*dy - d2f*y = 0
# -f*d2y = 3*y^2*dy^2 + y^3*d2y - 3*dh*y*dy - 3/2*d2h*y^2 - d2f*y
y^3*dy - 3*h*y*dy - 2*f*dy = 0
y^3*dy - 3*h*y*dy + 2*(3*y^2*dy^2 + y^3*d2y - 3*dh*y*dy - 3/2*d2h*y^2 - d2f*y) = 0
2*y^3*d2y + 6*y^2*dy^2 + y^3*dy - 3*h*y*dy - 6*dh*y*dy - 3*d2h*y^2 - 2*d2f*y = 0
# TODO: substitute dy^2;


### T.B.2: TODO: y^2;


################

################
### Examples ###
################

### Example 1;
# h(x) = x
# f(x) = x^3
# [not run] [NOT decoupled]
y^2*dy - h*dy - y*dh - 2/3*df = 0
y^2*dy - x*dy - y - 2*x^2 = 0
# y^2*dy - D(x*y) = 2*x^2

### Solution & Plot:
# y = (f + sqrt(f^2 - h^3))^(1/3) + (f - sqrt(f^2 - h^3))^(1/3)
y = function(x, n=3) {
	r1 = (x^3 + sqrt(x^6 - x^n + 0i))
	r2 = (x^3 - sqrt(x^6 - x^n + 0i))
	r = round0(rootn(r1, n=n) + rootn(r2, n=n))
	return(r)
}
dy = function(x) {
	y.x = y(x)
	div = (y.x^2 - x)
	dp = (y.x + 2*x^2)
	dp = ifelse(div != 0, dp / div, Inf); # TODO: check limit;
	return(dp)
}
# nice local minimum
curve(y, from=-2, to=2)
sapply(c(-(4:1)/5, (1:4)/5, 1.01), line.tan, dx=3, p=y, dp=dy)


### Variant: [decoupled]
y^2*dy - x*dy - y - 2*x^2 = 0 # *y
y^3*dy - x*y*dy - y^2 - 2*x^2*y = 0
(3*x*y + 2*x^3)*dy - x*y*dy - y^2 - 2*x^2*y = 0
### ODE:
2*x*y*dy + 2*x^3*dy - y^2 - 2*x^2*y # = 0
### variant forms:
# D(y^2/x) + 2*x^2 * D(y/x) = 0 # simple example
# (y(2)^2 + 2*2^2*y(2)) - integrate(function(x) Re(y(x)), lower=0, upper=2)$value * 4*2

### Solution & Plot:
# y = is the same as above;
dy = function(x) {
	y.x = y(x)
	div = (2*x*y.x + 2*x^3)
	dp = (y.x^2 + 2*x^2*y.x)
	dp = ifelse(div != 0, dp / div, Inf);
	return(dp)
}
curve(y, from=-2, to=2)
sapply(c(-(4:1)/5, (1:4)/5, 1.01), line.tan, dx=3, p=y, dp=dy)


### V.A.3: Combinations
# h(x) = x # as above
# f(x) = x^3
# b = 1/2 * 1/x
h*y^2*dy - h^2*dy - 1/3*y^3*dh + 2/3*f*dh - 2/3*h*df +
 + b*(2*h*y*dy + 2*f*dy - y^2*dh - 2/3*df*y) = 0
x*y^2*dy - x^2*dy - 1/3*y^3 + 2/3*x^3 - 2*x^3 +
 + b*(2*x*y*dy + 2*x^3*dy - y^2 - 2*x^2*y) = 0
x*y^2*dy + y*dy - 1/3*y^3 - 1/2/x*y^2 - x*y - 4/3*x^3 = 0
### ODE:
x^2*y^2*dy + x*y*dy - 1/3*x*y^3 - 1/2*y^2 - x^2*y - 4/3*x^4 # = 0

### Solution & Plot:
# y = is the same as above;
dy = function(x) {
	y.x = y(x)
	div = (x*y.x^2 + y.x)
	dp = (1/3*y.x^3 + 1/2/x*y.x^2 + x*y.x + 4/3*x^3)
	dp = ifelse(div != 0, dp / div, Inf);
	return(dp)
}
curve(y, from=-2, to=2)
sapply(c(-(4:1)/5, (1:4)/5, 1.01), line.tan, dx=3, p=y, dp=dy)


### V.A.4: Combinations
# as above;
# [NOT run]
h^2*y*dy + f*h*dy - 1/9*df*y^3 - 1/2*h*dh*y^2 + 2/9*f*df = 0
x^2*y*dy + x^4*dy - 1/3*x^2*y^3 - 1/2*x*y^2 + 2/3*x^5 = 0
x*y*dy + x^3*dy - 1/3*x*y^3 - 1/2*y^2 + 2/3*x^4 = 0
# V.A.3 - V.A.4: [initial Derivative / trivial]
x^2*y^2*dy - x^3*dy - x^2*y - 2*x^4 = 0 # /x^2
### ODE:
y^2*dy - x*dy - y - 2*x^2 # = 0

### Solution & Plot:
# y = is the same as above (previous section);
dy = function(x, alt=FALSE) {
	y.x = y(x)
	if(alt) {
		div = (y.x^2 - x)
		dp = (y.x + 2*x^2) # initial Derivative
	} else {
		div = (x*y.x + x^3)
		dp = (x*y.x^3 + 3/2*y.x^2 - 2*x^4) / 3
	}
	dp = ifelse(div != 0, dp / div, Inf);
	return(dp)
}
curve(y, from=-2, to=2)
sapply(c(-(4:1)/5, (1:4)/5, 1.01), line.tan, dx=3, p=y, dp=dy)


### V.B.1: Tr. Poly
# h(x) = x # as above
# f(x) = x^3
# [not run]
h*y^3*dy + f*h*dy - 1/2*dh*y^4 + f*dh*y - h*df*y # = 0
x*y^3*dy + x^4*dy - 1/2*y^4 + x^3*y - 3*x^3*y # = 0
### ODE:
x*y^3*dy + x^4*dy - 1/2*y^4 - 2*x^3*y # = 0

### Solution & Plot:
# y = is the same as above;
dy = function(x) {
	y.x = y(x)
	div = (x*y.x^3 + x^4)
	dp = (1/2*y.x^4 + 2*x^3*y.x)
	dp = ifelse(div != 0, dp / div, Inf);
	return(dp)
}
curve(y, from=-2, to=2)
sapply(c(-(4:1)/5, (1:4)/5, 1.01), line.tan, dx=3, p=y, dp=dy)

### V.B.1: Tr. Poly / alternative
h*y^3*dy + h*f*dy - 1/3*df*y^3 - 3/2*h*dh*y^2 + 2/3*f*df = 0
x*y^3*dy + x^4*dy - x^2*y^3 - 3/2*x*y^2 + 2*x^5 = 0 # / x
y^3*dy + x^3*dy - x*y^3 - 3/2*y^2 + 2*x^4 = 0
### Solution & Plot:
# y = is the same as above;
dy = function(x) {
	y.x = y(x)
	div = (y.x^3 + x^3)
	dp = (x*y.x^3 + 3/2*y.x^2 - 2*x^4)
	dp = ifelse(div != 0, dp / div, Inf);
	return(dp)
}
curve(y, from=-2, to=2)
sapply(c(-(4:1)/5, (1:4)/5, 1.01), line.tan, dx=3, p=y, dp=dy)


### D2y / Dx2
### D2y / Dx2: Basic
# y^2*dy - h*dy - y*dh - 2/3*df = 0 # =>
# y^2*d2y - h*d2y + 2*y*dy^2 - 2*dh*dy - y*d2h - 2/3*d2f = 0
y^2*d2y - x*d2y + 2*y*dy^2 - 2*dy - 4*x # = 0
# Sum with one of the other ODEs:
y^2*d2y - x*d2y + 2*y*dy^2 + 2*x*y*dy + 2*x^3*dy - 2*dy - y^2 - 2*x^2*y - 4*x # = 0

### Solution & Plot:
# y = is the same as above;
dy = function(x) {
	y.x = y(x)
	div = (2*x*y.x + 2*x^3)
	dp = (y.x^2 + 2*x^2*y.x)
	dp = ifelse((div != 0), dp / div, Inf);
	return(dp)
}
d2y = function(x) {
	y.x = y(x)
	dy.x = dy(x)
	div = (y.x^2 - x)
	dp = - (2*y.x*dy.x^2 + 2*x*y.x*dy.x + 2*x^3*dy.x - 2*dy.x - y.x^2 - 2*x^2*y.x - 4*x)
	dp = ifelse((div != 0), dp / div, Inf);
	return(dp)
}
curve(dy, from=-2, to=2, col="green")
curve(y, from=-2, to=2, add=T, col="grey")
sapply(c(-1.5, -0.5, (1:5)/4.7), line.tan, dx=3, p=dy, dp=d2y)


### D2y / Dx2: Transformed
# [NOT run]
y^2*d2y - x*d2y + 2*y*dy^2 - dy - dy - 4*x = 0
# Replace: dy with one from the other ODEs:
2*x*y*dy + 2*x^3*dy - y^2 - 2*x^2*y = 0 # =>
dy = (y^2 + 2*x^2*y) / (2*x*y + 2*x^3) # =>
(y^2 - x)*d2y + 2*y*((y^2 + 2*x^2*y) / (2*x*y + 2*x^3))^2 - 2*(y^2 + 2*x^2*y) / (2*x*y + 2*x^3) - 4*x = 0

### Solution & Plot:
# y = is the same as above;
dy = function(x) {
	y.x = y(x)
	div = (2*x*y.x + 2*x^3)
	dp = (y.x^2 + 2*x^2*y.x)
	dp = ifelse((div != 0), dp / div, Inf);
	return(dp)
}
d2y = function(x) {
	y.x = y(x)
	dy.x = dy(x)
	div = (y.x^2 - x)
	div2 = (2*x*y.x + 2*x^3)
	val2 = (y.x^2 + 2*x^2*y.x)
	val = 2*y.x*(val2 / div2)^2 - 2*val2 / div2 - 4*x
	dp = - val;
	dp = ifelse((div != 0), dp / div, Inf);
	return(dp)
}
curve(dy, from=-2, to=2, col="green")
curve(y, from=-2, to=2, add=T, col="grey")
sapply(c(-1.5, -0.5, 0.06 + (0:5)/3.9), line.tan, dx=3, p=dy, dp=d2y)


##############
##############
### Example 2:
# h(x) = x^2
# f(x) = 3*log(x)
y^2*dy - h*dy - y*dh - 2/3*df = 0
y^2*dy - x^2*dy - 2*x*y - 2/x = 0
### ODE:
x*y^2*dy - x^3*dy - 2*x^2*y - 2 = 0

### Test & Plot:
# y = (f + sqrt(f^2 - h^3))^(1/3) + (f - sqrt(f^2 - h^3))^(1/3)
y = function(x, n=3) {
	r1 = (3*log(x + 0i) + sqrt(9*log(x + 0i)^2 - x^(2*n) + 0i))
	### imaginary parts do NOT cancel for: x < 0;
	# sign.x = sign(x)
	# sign.x[x >= 0] = 1
	# r2 = (sign.x * 3*log(x + 0i) - sqrt(9*log(x + 0i)^2 - x^(2*n) + 0i))
	r1 = round0(rootn(r1, n=n))
	r2 = x^2 / r1
	return( round0(r1 + r2) )
}
dy = function(x) {
	y.x = y(x)
	div = (x*y.x^2 - x^3)
	dp = (2*x^2*y.x + 2)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=0.01, to=2, ylim=c(-5, 0.5)) # xlim=c(-2, 2)
sapply((1:4)/5, line.tan, dx=3, p=y, dp=dy)
# log(x) is complex for x < 0;
# curve(y, from=-2, to=-0.01, add=T, ylim=c(-5, 0.5))
# sapply(c(-(4:1)/5), line.tan, dx=3, p=y, dp=dy)


### Transformed & Variants:
# h(x) = x^2
# f(x) = 3*log(x)
y^2*dy - h*dy - y*dh - 2/3*df = 0
y^2*dy - x^2*dy - 2*x*y - 2/x = 0
x*y^2*dy - x^3*dy - 2*x^2*y - 2 = 0 # *y
x*y^3*dy - x^3*y*dy - 2*x^2*y^2 - 2*y = 0
x*(3*x^2*y + 6*log(x))*dy - x^3*y*dy - 2*x^2*y^2 - 2*y = 0
x^3*y*dy + 3*x*log(x)*dy - x^2*y^2 - y = 0
# Solution & Plot:
# y = (f + sqrt(f^2 - h^3))^(1/3) + (f - sqrt(f^2 - h^3))^(1/3)
y = function(x, n=3) {
	r1 = (3*log(x + 0i) + sqrt(9*log(x + 0i)^2 - x^(2*n) + 0i))
	### imaginary parts do NOT cancel for: x < 0;
	# sign.x = sign(x)
	# sign.x[x >= 0] = 1
	# r2 = (sign.x * 3*log(x + 0i) - sqrt(9*log(x + 0i)^2 - x^(2*n) + 0i))
	r1 = round0(rootn(r1, n=n))
	r2 = x^2 / r1
	return( round0(r1 + r2) )
}
dy = function(x) {
	y.x = y(x)
	div = (x^3 * y.x + 3*x*log(x))
	dp = (x^2 * y.x^2 + y.x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=0.01, to=3)
sapply((1:8)/5, line.tan, dx=3, p=y, dp=dy)

### Variant 1:
x*y^2*dy - x^3*dy - 2*x^2*y - 2 = 0
x*y^2*dy - x^3*dy - 2/3*(y^3 - 6*log(x)) - 2 = 0
x*y^2*dy - x^3*dy - 2/3*y^3 + 4*log(x) - 2 = 0
### Solution & Plot:
dy = function(x) {
	y.x = y(x)
	div = (x * y.x^2 - x^3)
	dp = (2/3 * y.x^3 - 4*log(x) + 2)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=0.01, to=3)
sapply((1:8)/5, line.tan, dx=3, p=y, dp=dy)

### Variant 2:
x*y^2*dy - x^3*dy - 2*x^2*y - 2 = 0
1/3 * 1/x *y*(y^3 - 6*log(x))*dy - x^3*dy - 2*x^2*y - 2 = 0 # * 3*x
y*(y^3 - 6*log(x))*dy - 3*x^4*dy - 6*x^3*y - 6*x = 0
y^4*dy - 6*log(x)*y*dy - 3*x^4*dy - 6*x^3*y - 6*x = 0
### Solution & Plot:
dy = function(x) {
	y.x = y(x)
	div = (y.x^4 - 6*y.x*log(x) - 3*x^4)
	dp = 6*(y.x*x^3 + x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=0.01, to=3)
sapply((1:8)/5, line.tan, dx=3, p=y, dp=dy)



###############
###############

### Transformed Polynomial
# - "transformation" of initial polynomial;

#########
### n = 3
# y^3 - 3*h*y - 2*f = 0
# y^2 - 3*h - 2*f/y = 0
# 2*y*dy - 3*dh - 2*df/y + 2*f*dy/y^2 = 0
# y^3*dy - 3/2*dh*y^2 - df*y + f*dy = 0

### Examples
# h(x) = x^2
# f(x) = x + 1
# [not run]
2*y^3*dy - 3*dh*y^2 - 2*df*y + 2*f*dy = 0
2*y^3*dy - 2*3*x*y^2 - 2*y + 2*(x+1)*dy = 0
y^3*dy + (x+1)*dy - 3*x*y^2 - y = 0
### Solution:
# y = (f + sqrt(f^2 - h^3))^(1/3) + (f - sqrt(f^2 - h^3))^(1/3)
y = function(x, n=3) {
	det = sqrt((x+1)^2 - x^6 + 0i)
	r1 = round0(rootn(x+1 + det, n=3))
	r2 = round0(rootn(x+1 - det, n=3))
	return( round0(r1 + r2) )
}
dy = function(x) {
	y.x = y(x)
	div = (y.x^3 + x + 1)
	dp = (3*x*y.x^2 + y.x)
	dp = ifelse(div != 0, dp / div, Inf);
	return(dp)
}
curve(y, from=-3, to=3)
sapply(c(-2, c(-(4:1), 1:4)/5, 2), line.tan, dx=3, p=y, dp=dy)

### Variant
# y^3*dy - 3/2*dh*y^2 - df*y + f*dy = 0
# (3*h*y + 2*f)*dy - 3/2*dh*y^2 - df*y + f*dy = 0
# 3*h*y*dy + 3*f*dy - 3/2*dh*y^2 - df*y = 0

### Examples
# h(x) = x^2
# f(x) = x + 1
3*h*y*dy + 3*f*dy - 3/2*dh*y^2 - df*y = 0
3*x^2*y*dy + 3*(x+1)*dy - 3*x*y^2 - y = 0
x^2*y*dy + (x+1)*dy - x*y^2 - 1/3*y = 0
# y = (f + sqrt(f^2 - h^3))^(1/3) + (f - sqrt(f^2 - h^3))^(1/3)
y = function(x, n=3) {
	det = sqrt((x+1)^2 - x^6 + 0i)
	r1 = round0(rootn(x+1 + det, n=3))
	r2 = round0(rootn(x+1 - det, n=3))
	return( round0(r1 + r2) )
}
dy = function(x) {
	y.x = y(x)
	div = (x^2*y.x + x + 1)
	dp = (x*y.x^2 + 1/3 * y.x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-3, to=3)
sapply(c(-2, c(-(4:1), 1:4)/5, 2), line.tan, dx=3, p=y, dp=dy)


### Example [T.B.1]:
# h(x) = 2*x
# f(x) = 1 # Test [trivial]
# [NOT run]
y^3*dy + f*dy - 3/2*dh*y^2 - df*y = 0
y^3*dy + dy - 3*y^2  = 0
### Solution & Plot:
# y = (f + sqrt(f^2 - h^3))^(1/3) + (f - sqrt(f^2 - h^3))^(1/3)
y = function(x, n=3) {
	d = 1; h.x = 2*x
	r1 = (d + sqrt(d^2 - h.x^n + 0i))
	r2 = (d - sqrt(d^2 - h.x^n + 0i))
	r = round0(rootn(r1, n=n) + rootn(r2, n=n))
	return(r)
}
dy = function(x) {
	y.x = y(x)
	div = (y.x^3 + 1)
	dp = (3*y.x^2)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-3, to=3)
sapply(c(-2, -1, -0.4, (2:5)/5), line.tan, dx=3, p=y, dp=dy)

### T.B.1 + A.1:
# [NOT run]
h*y*dy + f*dy - 1/2*dh*y^2 - 1/3*df*y = 0
(2*x*y + 1)*dy - y^2 = 0
### Solution & Plot:
# - same y as above;
dy = function(x) {
	y.x = y(x)
	div = (2*x*y.x + 1)
	dp = (y.x^2)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-3, to=3)
sapply(c(-2, -1, -0.4, (2:5)/5), line.tan, dx=3, p=y, dp=dy)
# same plot as above!


### T.B.1 + A.1: Example 2
# h(x) = x^2 + 1
# f(x) = x^3
# [NOT run]
h*y*dy + f*dy - 1/2*dh*y^2 - 1/3*df*y = 0
(x^2*y + y + x^3)*dy - x*y^2 - x^2*y = 0
### Solution & Plot:
y = function(x, n=3) {
	d = x^3; h.x = x^2 + 1
	r1 = (d + sqrt(d^2 - h.x^n + 0i))
	r2 = (d - sqrt(d^2 - h.x^n + 0i))
	r = round0(rootn(r1, n=n) + rootn(r2, n=n))
	return(r)
}
dy = function(x) {
	y.x = y(x)
	div = (x^2*y.x + y.x + x^3)
	dp = (x*y.x^2 + x^2*y.x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-3, to=3)
sapply(c(-1.5, (0:3)/1.7), line.tan, dx=3, p=y, dp=dy)




#########################
#########################

##################
### P3 Shifted ###

#########
### n = 3
### Shifted
# y^3 + 3*s*y^2 - 3*h*y - 2*f = 0
# 3*y^2*dy + 6*s*y*dy + 3*y^2*ds - 3*h*dy - 3*y*dh - 2*df = 0
# =>
# y^2*dy + 2*s*y*dy - h*dy + y^2*ds - y*dh - 2/3*df = 0
# (y^2 + 2*s*y - h)*dy + y^2*ds - y*dh - 2/3*df = 0

### Solution: reduce shift
# y => y - s =>
# y^3 - 3*(h + s^2)*y - 2*f + 2*s^3 + 3*h*s = 0

###
# h(x) = x
# f(x) = x^3
# s(x) = 1/2
y^2*dy + 2*s*y*dy - h*dy + y^2*ds - y*dh - 2/3*df = 0
y^2*dy + y*dy - x*dy - y - 2*x^2 = 0
# y = - s + (f + sqrt(f^2 - h^3))^(1/3) + (f - sqrt(f^2 - h^3))^(1/3)
y = function(x, n=3) {
	d = x^3 - 3/4*x - 1/8
	det = sqrt(d^2 - (x + 1/4)^3 + 0i)
	r1 = (d + det); r2 = (d - det)
	r = round0(rootn(r1, n=n) + rootn(r2, n=n))
	# shift back
	r = r - 1/2
	return(r)
}
dy = function(x) {
	y.x = y(x)
	div = (y.x^2 + y.x - x)
	dp = (y.x + 2*x^2)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-2, to=2)
sapply(c(-(4:1)/5, (1:6)/5), line.tan, dx=3, p=y, dp=dy)


#########
### n = 3
### y-Shifted
# y^3 + 3*s*y^2 - 3*h*y - 2*f = 0
# 3*y^2*dy + 6*s*y*dy + 3*y^2*ds - 3*h*dy - 3*y*dh - 2*df = 0
# =>
# y^2*dy + 2*s*y*dy - h*dy + y^2*ds - y*dh - 2/3*df = 0
# (y^2 + 2*s*y - h)*dy + y^2*ds - y*dh - 2/3*df = 0

### y-Shift: u = y*m
# (m^2*y^2 + 2*s*m*y - h)*(m*dy + y*dm) + m^2*y^2*ds - m*y*dh - 2/3*df = 0
# (m^3*y^2 + 2*s*m^2*y - h*m)*dy + m^2*y^3*dm + (m^2*ds + 2*m*s*dm)*y^2 - (m*dh + h*dm)*y - 2/3*df = 0
# or
### y-Shift: u = y + n
# ((y+n)^2 + 2*s*(y+n) - h)*(dy + dn) + (y + n)^2*ds - (y + n)*dh - 2/3*df = 0
# (y^2 + 2*(s+n)*y + n^2 + 2*s*n - h)*(dy + dn) + ds*y^2 - (dh - 2*n*ds)*y - 2/3*df + n^2*ds - n*dh = 0
# (y^2 + 2*(s+n)*y + n^2 + 2*s*n - h)*dy + (ds + dn)*y^2 - (dh - 2*n*ds - 2*(s+n))*y - 2/3*df + n^2*ds - n*dh +(n^2 + 2*s*n - h)*dn = 0

### Solution: reduce shift
# y => y - s =>
# y^3 - 3*(h + s^2)*y - 2*f + 2*s^3 + 3*h*s = 0

###
# h(x) = x - 1/4
# f(x) = x^3
# s(x) = 1/2
# n(x) = -1/2
(y^2 + 2*(s+n)*y + n^2 + 2*s*n - h)*dy + (ds + dn)*y^2 - (dh - 2*n*ds - 2*(s+n))*y - 2/3*df + n^2*ds - n*dh +(n^2 + 2*s*n - h)*dn = 0
y^2*dy - x*dy - y - 2*x^2 + 1/2 = 0
# simple: d(y^3 - 3*x*y) = 6*x^2 - 3/2;
# [but useful to test formulas]
# y = n - s + (f.sh + sqrt(f.sh^2 - h.sh^3))^(1/3) + (f.sh - sqrt(f.sh^2 - h.sh^3))^(1/3)
y = function(x, n=3) {
	d = x^3 - 3/4*x + 3/4*1/4 - 1/8
	det = sqrt(d^2 - (x - 1/4 + 1/4)^3 + 0i)
	r1 = (d + det); r2 = (d - det)
	r = round0(rootn(r1, n=n) + rootn(r2, n=n))
	# shift back
	r = r - 1/2 + 1/2
	return(r)
}
dy = function(x) {
	y.x = y(x)
	div = (y.x^2 - x)
	dp = (y.x + 2*x^2 - 1/2)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-2, to=2)
sapply(c(-(4:1)/5, (1:6)/5), line.tan, dx=3, p=y, dp=dy)


# df = from free term;
# y*dy, dy, y^2, y: 4 parameters needed;
# but currently only 3 available: h(x), s(x), n(x);


#############

#############
### n = 5 ###

### Cardano-Polynomial:
y^5 - 5*h*y^3 + 5*h^2*y - 2*f # = 0

### D =>
# 5*y^4*dy - 15*h*y^2*dy - 5*y^3*dh + 5*h^2*dy + 10*h*y*dh - 2*df = 0
### Primary ODE:
y^4*dy - 3*h*y^2*dy + h^2*dy - y^3*dh + 2*h*y*dh - 2/5 * df # = 0

### Entanglement:
y^5*dy - 3*h*y^3*dy + h^2*y*dy - y^4*dh + 2*h*y^2*dh - 2/5 * df*y # = 0
(5*h*y^3 - 5*h^2*y + 2*f)*dy - 3*h*y^3*dy + h^2*y*dy - y^4*dh + 2*h*y^2*dh - 2/5 * df*y # = 0
### ODE:
2*h*y^3*dy - 4*h^2*y*dy + 2*f*dy - y^4*dh + 2*h*y^2*dh - 2/5 * df*y # = 0

### D2 =>
2*h*y^3*d2y + 6*y^2*dy^2 - 4*h^2*y*d2y - 4*h^2*dy^2 - 8*h*dh*y*dy +
	+ 2*f*d2y + 2*df*dy - 4*dh*y^3*dy - d2h*y^4 +
	+ 4*h*dh*y*dy + 2*(dh^2 + h*d2h)*y^2 - 2/5 * df*dy - 2/5 * d2f*y # = 0
# TODO: can substitute y^4 using D1;

### TODO:
# - concrete examples;

### Example 1:
# h(x) = b0; # constant
# f(x) = 5*b0*x;
### ODE:
y^3*dy - 2*b0*y*dy + 5*x*dy - y # = 0


######################
######################

### Shifted Symmetrically

### Base System:
# (p-s)^3 + (q-s)^3 = r
# p*q = c
# y = p + q

y^2*dy - 2*s*y*dy + (s^2 - c)*dy - y^2*ds + (2*s*ds - dc)*y = 1/3*dr - 2*c*ds - 2*s*dc + 2*s^2*ds

### Examples

### Example 1:
# s = 1
# c = x
# r = x^3
y^2*dy - 2*y*dy - (x-1)*dy - y = x^2 - 2

### Solution & Plot:
y = function(x, n=3) {
	s = 1
	d = 1/2 * (x^3 - 3*x + 1)
	det = sqrt(d^2 - x^3 + 0i)
	r1 = (d + det); r2 = (d - det)
	r = round0(rootn(r1, n=n) + rootn(r2, n=n))
	# shift back
	r = r + s
	return(r)
}
dy = function(x) {
	y.x = y(x)
	div = (y.x^2 - 2*y.x - x + 1)
	dp = (y.x + x^2 - 2)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-2, to=2)
line.tan(c(-(4:1)/5, (1:6)/5), dx=3, p=y, dp=dy)

###########

### Shifted Symmetrically

### Base System:
# (p-s1)^3 + (q-s1)^3 + 3*b1*(p+q) = r
# (p-s2)*(q-s2) = c
# y = p + q

### Solution to Polynomial:
# (p-s2)*(q-s2) = c # =>
# p*q - s2*(p+q) + s2^2 - c # = 0
# p*q = s2*y + c - s2^2

p^3 + q^3 - 3*s1*(p^2+q^2) + 3*s1^2*(p+q) - 2*s1^3 + 3*b1*(p+q) = r
y^3 - 3*p*q*y - 3*s1*(y^2 - 2*p*q) + 3*(s1^2+b1)*y - 2*s1^3 - r = 0
y^3 - 3*(s2*y + c - s2^2)*y - 3*s1*(y^2 - 2*(s2*y + c - s2^2)) + 3*(s1^2+b1)*y - 2*s1^3 - r = 0
y^3 - 3*(s1+s2)*y^2 + 3*(s1^2 + s2^2 + 2*s1*s2 + b1 - c)*y + 6*s1*c - 2*s1^3 - 6*s1*s2^2 - r # = 0
# Shift: y => y + (s1+s2)
# Solution is based on this polynomial:
y^3 + (3*b1 - 3*c)*y + (3*b1*s1 + 3*b1*s2 + 3*c*s1 - 3*c*s2 - r - 3*s1*s2^2 + 3*s1^2*s2 - s1^3 + s2^3)

### ODE
y^2*dy - 2*(s1+s2)*y*dy - (ds1+ds2)*y^2 + (s1^2 + s2^2 + 2*s1*s2 + b1 - c)*dy +
 + (2*s1*ds1 + 2*s2*ds2 + 2*ds1*s2 + 2*s1*ds2 + db1 - dc)*y # =
# = - (2*ds1*c + 2*s1*dc - 2*s1^2*ds1 - 2*ds1*s2^2 - 4*s1*s2*ds2 - 1/3 * dr)

### Example 1:
# s1 = x
# s2 = -x
# b1 = x^2
# c = x
# r = x^3
y^2*dy + (x^2 - x)*dy + (2*x-1)*y = 9*x^2 - 4*x

### Solution & Plot:
y = function(x, n=3) {
	s1 = x; s2 = -x;
	b1 = x^2; c = x; r = x^3;
	d = -1/2 * (3*b1*s1 + 3*b1*s2 + 3*c*s1 - 3*c*s2 - r - 3*s1*s2^2 + 3*s1^2*s2 - s1^3 + s2^3)
	det = sqrt(d^2 - (c - b1)^3 + 0i)
	r1 = (d + det); r2 = (d - det)
	r.sol = round0(rootn(r1, n=n) + rootn(r2, n=n))
	# shift back
	r.sol = r.sol + s1 + s2
	return(r.sol)
}
dy = function(x) {
	y.x = y(x)
	div = (y.x^2 + x^2 - x)
	dp =  - y.x * (2*x-1) + 9*x^2 - 4*x
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-3, to=3)
# a nice local minimum
line.tan(c(-(4:1)/5, (1:6)/5), dx=3, p=y, dp=dy)


####################
####################
####################

### Arbitrary Polynomials

### y^5
# y^5 - 5*b1*y - b0 = 0
# 5*y^4*dy - 5*b1*dy - 5*db1*y - db0 = 0
# y^4*dy - b1*dy - db1*y - 1/5*db0 = 0 # *y
# y^5*dy - b1*y*dy - db1*y^2 - 1/5*db0*y = 0
# Transform A.1:
(5*b1*y + b0)*dy - b1*y*dy - db1*y^2 - 1/5*db0*y = 0
(4*b1*y + b0)*dy - db1*y^2 - 1/5*db0*y = 0



### y^n
# y^n - n*b1*y - b0 = 0
# n*y^(n-1)*dy - n*b1*dy - n*db1*y - db0 = 0
# y^(n-1)*dy - b1*dy - db1*y - 1/n*db0 = 0 # *y
# y^n*dy - b1*y*dy - db1*y^2 - 1/n*db0*y = 0
# Transform A.1:
(n*b1*y + b0)*dy - b1*y*dy - db1*y^2 - 1/n*db0*y = 0
((n-1)*b1*y + b0)*dy - db1*y^2 - 1/n*db0*y = 0

### Example:
# b0 = x
# b1 = x
x*((n-1)*y + 1)*dy - y^2 - 1/n*y = 0
# classic:
((n-1)*y + 1)/(y^2 + 1/n*y) * dy = 1/x
(n-1)/2*(2*y + 1/n - 1/n + 2/(n-1))/(y^2 + 1/n*y) * dy = 1/x


### Example:
# b0 = k1 # constant
# b1 = x^p - k2
((n-1)*(x^p - k2)*y + k1)*dy - p*x^(p-1)*y^2 = 0


#######################
#######################

### Class 1 Polynomials

### n = 5

###
x = sqrt(2) # some Test
k = (x^2-x+1)^(1/5)
# root
y = (x^2+x)*k^4 - (x^2+2*x+1)*k^3 + x^2*k^2 + (x^2+x)*k
# P5 Polynomial (in y):
(1 + 7*x + 21*x^2 + 28*x^3 + 6*x^4 - 18*x^5 - 10*x^6 - 69*x^7 - 163*x^8 - 122*x^9 - 122*x^10 - 184*x^11 - 114*x^12 +
  - 50*x^13 - 64*x^14 - 27*x^15 + x^16 - x^17 - x^18) +
+ (5*x + 25*x^2 + 55*x^3 + 85*x^4 + 105*x^5 + 125*x^6 + 155*x^7 + 150*x^8 + 90*x^9 + 85*x^10 + 75*x^11 + 20*x^12 +
  + 10*x^13 + 5*x^14)*y +
+ (- 5*x - 10*x^2 - 5*x^4 - 25*x^5 - 15*x^6 - 10*x^7 - 15*x^8 - 15*x^9 - 10*x^10)*y^2 +
+ y^5

# D(P): The Answer to the Universe and Everything
42*x*y + 50*x*y^2 - 20*x*y^3 + 84*x^2*y + 165*x^2*y^2 + 24*x^3*y + 340*x^3*y^2 - 20*x^3*y^3 - 90*x^4*y +
 + 525*x^4*y^2 - 125*x^4*y^3 - 60*x^5*y + 750*x^5*y^2 - 90*x^5*y^3 - 483*x^6*y + 1085*x^6*y^2 - 70*x^6*y^3 - 1304*x^7*y +
 + 1200*x^7*y^2 - 120*x^7*y^3 - 1098*x^8*y + 810*x^8*y^2 - 135*x^8*y^3 - 1220*x^9*y + 850*x^9*y^2 - 100*x^9*y^3 +
 - 2024*x^10*y + 825*x^10*y^2 - 1368*x^11*y + 240*x^11*y^2 - 650*x^12*y + 130*x^12*y^2 - 896*x^13*y +
 + 70*x^13*y^2 - 405*x^14*y + 16*x^15*y - 17*x^16*y - 18*x^17*y + 7*y + 5*y^2 - 5*y^3 +
(- 5 - 35*x - 105*x^2 - 140*x^3 - 30*x^4 + 90*x^5 + 50*x^6 + 345*x^7 + 815*x^8 + 610*x^9 + 610*x^10 +
 + 920*x^11 + 570*x^12 + 250*x^13 + 320*x^14 + 135*x^15 - 5*x^16 + 5*x^17 + 5*x^18) * dy +
(- 20*x - 100*x^2 - 220*x^3 - 340*x^4 - 420*x^5 - 500*x^6 - 620*x^7 - 600*x^8 - 360*x^9 - 340*x^10 - 300*x^11 +
 - 80*x^12 - 40*x^13 - 20*x^14) * y * dy +
(15*x + 30*x^2 + 15*x^4 + 75*x^5 + 45*x^6 + 30*x^7 + 45*x^8 + 45*x^9 + 30*x^10) * y^2 * dy # = 0

### Solution & Plot:
y = function(x, n=5) {
	# root
	k = rootn(x^2-x+1, n)
	y = (x^2+x)*k^4 - (x^2+2*x+1)*k^3 + x^2*k^2 + (x^2+x)*k
	y = sapply(y, round0)
	return(y)
}
y.free = function(x, y) {
	42*x*y + 50*x*y^2 - 20*x*y^3 + 84*x^2*y + 165*x^2*y^2 + 24*x^3*y + 340*x^3*y^2 - 20*x^3*y^3 - 90*x^4*y +
 + 525*x^4*y^2 - 125*x^4*y^3 - 60*x^5*y + 750*x^5*y^2 - 90*x^5*y^3 - 483*x^6*y + 1085*x^6*y^2 - 70*x^6*y^3 - 1304*x^7*y +
 + 1200*x^7*y^2 - 120*x^7*y^3 - 1098*x^8*y + 810*x^8*y^2 - 135*x^8*y^3 - 1220*x^9*y + 850*x^9*y^2 - 100*x^9*y^3 +
 - 2024*x^10*y + 825*x^10*y^2 - 1368*x^11*y + 240*x^11*y^2 - 650*x^12*y + 130*x^12*y^2 - 896*x^13*y +
 + 70*x^13*y^2 - 405*x^14*y + 16*x^15*y - 17*x^16*y - 18*x^17*y + 7*y + 5*y^2 - 5*y^3
}
dy = function(x) {
	y.x = y(x)
	div = y.x^2 * (15*x + 30*x^2 + 15*x^4 + 75*x^5 + 45*x^6 + 30*x^7 + 45*x^8 + 45*x^9 + 30*x^10) +
		y.x * (- 20*x - 100*x^2 - 220*x^3 - 340*x^4 - 420*x^5 - 500*x^6 - 620*x^7 - 600*x^8 - 360*x^9 - 340*x^10 - 300*x^11 +
			- 80*x^12 - 40*x^13 - 20*x^14) +
		(- 5 - 35*x - 105*x^2 - 140*x^3 - 30*x^4 + 90*x^5 + 50*x^6 + 345*x^7 + 815*x^8 + 610*x^9 + 610*x^10 +
			+ 920*x^11 + 570*x^12 + 250*x^13 + 320*x^14 + 135*x^15 - 5*x^16 + 5*x^17 + 5*x^18)
	dp =  -y.free(x, y.x)
	
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-2, to=2)
# a nice ?global/local? minimum
line.tan(c(-(4:1)/5, (1:6)/5), dx=3, p=y, dp=dy)


#############

### Example: simple n = 5
x = sqrt(2) # some Test
k = (x - 1)^(1/5)
# root
y = k^4 - k^3 + k^2 + k
### P5 Polynomial (in y):
(- 22 + 48*x - 30*x^2 + 5*x^3 - x^4) +
(15 - 20*x + 5*x^3)*y^1 +
(- 10 + 20*x - 10*x^2)*y^2 + y^5
### D(P5)
48*y - 60*x*y + 15*x^2*y - 4*x^3*y - 20*y^2 + 15*x^2*y^2 + 20*y^3 - 20*x*y^3 +
 + (110 - 240*x + 150*x^2 - 25*x^3 + 5*x^4) * dy +
 + (- 60 + 80*x - 20*x^3) * y * dy +
 + (30 - 60*x + 30*x^2) * y^2 * dy

### Solution & Plot:
y = function(x, n=5) {
	# root
	k = rootn(x - 1, n)
	y = k^4 - k^3 + k^2 + k
	y = sapply(y, round0)
	return(y)
}
y.free = function(x, y) {
	48*y - 60*x*y + 15*x^2*y - 4*x^3*y - 20*y^2 + 15*x^2*y^2 + 20*y^3 - 20*x*y^3
}
dy = function(x) {
	y.x = y(x)
	div = y.x^2 * (30 - 60*x + 30*x^2) +
		y.x * (- 60 + 80*x - 20*x^3) +
		(110 - 240*x + 150*x^2 - 25*x^3 + 5*x^4)
	dp =  -y.free(x, y.x)
	
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
### Plot:
curve(y, from=-2, to=3)
# a nice global minimum
line.tan(c(-(4:1)/3, (1:6)/3), dx=3, p=y, dp=dy)


#############

### Example: simple n = 5
a = 4 # Parameter
x = sqrt(2) # some Test
k = rootn(x - a, 5)
# root
y = k^4 - k^3 + k^2 + k
### P5 Polynomial (in y):
# y^5 + b2(a)*y^2 + b1(a)*y + b0(a)
(a + 42*a*x - 3*a*x^2 + 4*a*x^3 - 21*a^2 + 3*a^2*x - 6*a^2*x^2 - a^3 + 4*a^3*x - a^4 - x - 21*x^2 + x^3 - x^4) +
(5*a - 30*a*x - 15*a*x^2 + 15*a^2 + 15*a^2*x - 5*a^3 - 5*x + 15*x^2 + 5*x^3)*y +
(20*a*x - 10*a^2 - 10*x^2)*y^2 + y^5
### ODE based on D(P5)
- 6*a*x*y - 30*a*x*y^2 + 12*a*x^2*y + 42*a*y - 30*a*y^2 + 20*a*y^3 - 12*a^2*x*y + 3*a^2*y + 15*a^2*y^2 +
	+ 4*a^3*y - 42*x*y + 30*x*y^2 - 20*x*y^3 + 3*x^2*y + 15*x^2*y^2 - 4*x^3*y - y - 5*y^2 +
 + (- 5*a - 210*a*x + 15*a*x^2 - 20*a*x^3 + 105*a^2 - 15*a^2*x + 30*a^2*x^2 + 5*a^3 - 20*a^3*x +
	+ 5*a^4 + 5*x + 105*x^2 - 5*x^3 + 5*x^4) * dy +
 + (- 20*a + 120*a*x + 60*a*x^2 - 60*a^2 - 60*a^2*x + 20*a^3 + 20*x - 60*x^2 - 20*x^3) * y * dy +
 + (- 60*a*x + 30*a^2 + 30*x^2) * y^2 * dy # = 0
 
### Solution & Plot:
y = function(x, a, n=5) {
	# root
	k = rootn(x - a, n)
	y = k^4 - k^3 + k^2 + k
	y = sapply(y, round0)
	return(y)
}
y.free = function(x, y, a) {
	- 6*a*x*y - 30*a*x*y^2 + 12*a*x^2*y + 42*a*y - 30*a*y^2 + 20*a*y^3 - 12*a^2*x*y +
	+ 3*a^2*y + 15*a^2*y^2 + 4*a^3*y - 42*x*y + 30*x*y^2 - 20*x*y^3 + 3*x^2*y + 15*x^2*y^2 - 4*x^3*y - y - 5*y^2
}
dy = function(x, a) {
	y.x = y(x, a)
	div = y.x^2 * (- 60*a*x + 30*a^2 + 30*x^2) +
		y.x * (- 20*a + 120*a*x + 60*a*x^2 - 60*a^2 - 60*a^2*x + 20*a^3 + 20*x - 60*x^2 - 20*x^3) +
		(- 5*a - 210*a*x + 15*a*x^2 - 20*a*x^3 + 105*a^2 - 15*a^2*x + 30*a^2*x^2 + 5*a^3 - 20*a^3*x +
		+ 5*a^4 + 5*x + 105*x^2 - 5*x^3 + 5*x^4)
	dp =  -y.free(x, y.x, a)
	
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
### Plot:
curve(y(x, a=a), from=a-2, to=a+3)
# a nice global minimum
line.tan(c(a-(4:1)/3, a+(1:6)/3), dx=3, p=y, dp=dy, a=a)


###########################
### Cardano-Polynomials ###
###########################

### System:
# p^n + q^n = 2*f(x)
# p*q = h(x)

### full P3: see previous section;
# y = p + q;

### P6-partial:
# => q = h / p
# =>
# n*p^(n-1)*dp + n*q^(n-1)*dq = 2 * df
# p^(n-1)*dp + q^(n-1)*dq = 2/n * df
# p^(n-1)*dp + (h/p)^(n-1) * (dh/p - h*dp/p^2) = 2/n * df
# p^(2*n)*dp + h^(n-1) * (p*dh - h*dp) - 2/n * p^(n+1) * df = 0
#
# p^(2*n)*dp - h^n*dp + h^(n-1)*p*dh - 2/n * p^(n+1) * df = 0

### Solutions
# p = (f + sqrt(f^2 - h^n))^(1/n)
# q = (f - sqrt(f^2 - h^n))^(1/n)
# Note:
# - 2 basic solutions are possible;
# - it is possible to rotate these solutions using the roots of unity;


### Transformations:

### Order 3: [redundant]
# same as regular solution;
p^2*dp + q^2*dq = 2/3*df # * p^2
p^4*dp + p^2*q^2*dq - 2/3*df*p^2 = 0
p^4*dp + h^2*dq - 2/3*df*p^2 = 0
p^4*dp + h^2*(dh - q*dp)/p - 2/3*df*p^2 = 0 # * p
p^5*dp + h^2*(dh - q*dp) - 2/3*df*p^3 = 0 # *p
p^6*dp + h^2*(dh*p - h*dp) - 2/3*df*p^4 = 0

### Variant:
p^2*dp + q^2*dq = 2/3*df
p^4*dp^2 + q^4*dq^2 + 2*p^2*q^2*dp*dq = 4/9*df^2
p^4*dp^2 + q^4*((dh - q*dp)/p)^2 + 2*h^2*dp*(dh - q*dp)/p = 4/9*df^2 # *p^2
p^6*dp^2 + q^4*(dh - q*dp)^2 + 2*h^2*dp*(p*dh - p*q*dp) = 4/9*df^2
p^6*dp^2 + q^4*(dh - q*dp)^2 + 2*h^2*dp*(p*dh - h*dp) = 4/9*df^2
p^6*dp^2 + q^4*(dh - q*dp)^2 + 2*h^2*dh*p*dp - 2*h^3*dp^2 = 4/9*df^2
# TODO: ...


################
### Examples ###

#########
### n = 2
p^4*dp - h^2*dp + h*p*dh - p^3 * df = 0

###
# h(x) = x
# f(x) = 1
p^4*dp - x^2*dp + x*p = 0
# p = sqrt(f + sqrt(f^2 - h^2))
p = function(x) {
	sqrt(1 + sqrt(1 - x^2))
}
dp = function(x) {
	p.x = p(x)
	div = (p.x^4 - x^2)
	dp = - x*p.x
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=-1, to=1)
sapply(c((0:5)/6), line.tan, dx=2, p=p, dp=dp)


###
# h(x) = x
# f(x) = 2*x
p^4*dp - x^2*dp - 2*p^3 + x*p = 0
# p = sqrt(f + sqrt(f^2 - h^2))
p = function(x) {
	sqrt(2*x + x*sqrt(3))
}
dp = function(x) {
	p.x = p(x)
	div = (p.x^4 - x^2)
	dp = 2*p.x^3 - x*p.x
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=0, to=3)
sapply(c((0:5)/4.5), line.tan, dx=2, p=p, dp=dp)

#########

#########
### n = 3
p^6*dp - h^3*dp - 2/3 * p^4 * df + h^2*p*dh = 0
# where df, dh = given;

###
# h(x) = x
# f(x) = 1
p^6*dp - x^3*dp + x^2*p = 0
# Solution & Plot:
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (1 + sqrt(1 - x^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
dp = function(x) {
	div = (p(x)^6 - x^3)
	dp = - x^2*p(x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=-3, to=1)
sapply(c((0:5)/6), line.tan, dx=2, p=p, dp=dp)


###
# h(x) = x
# f(x) = 3*x
p^6*dp - x^3*dp - 2*p^4 + x^2*p = 0
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (3*x + sqrt(9*x^2 - x^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
dp = function(x) {
	div = (p(x)^6 - x^3)
	dp = 2*p(x)^4 - x^2*p(x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=-3, to=9)
sapply(c(0.001, 2*(1:3), -1), line.tan, dx=3)

###
# h(x) = x
# f(x) = x^3 + 3*x
p^6*dp - x^3*dp - 2*(x^2+1)*p^4 + x^2*p = 0
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (x^3 + 3*x + sqrt((x^3 + 3*x)^2 - x^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
dp = function(x) {
	div = (p(x)^6 - x^3)
	dp = 2*(x^2+1)*p(x)^4 - x^2*p(x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=-3, to=9)
sapply(c(0.001, 2*(1:3), -1), line.tan, dx=2)


###
# h(x) = x
# f(x) = 3/2 * log(x)
p^6*dp - h^3*dp - 2/3 * p^4 * df + h^2*p*dh = 0
p^6*dp - x^3*dp - 1/x*p^4 + x^2*p = 0
x*p^6*dp - x^4*dp - p^4 + x^3*p = 0
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (3/2 * log(x) + sqrt((3/2 * log(x))^2 - x^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
dp = function(x) {
	div = x*(p(x)^6 - x^3)
	dp = p(x)^4 - x^3*p(x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=0, to=0.8)
sapply(c(0.001, (1:3)/6), line.tan, dx=2)


###
# h(x) = 1/(x^2 + 1)
# f(x) = 3/2 * log(x^2 + 1)
p^6*dp - h^3*dp - 2/3 * p^4 * df + h^2*p*dh = 0
(x^2 + 1)^4 * p^6 * dp - (x^2 + 1)*dp - 2*x*(x^2 + 1)^3 * p^4 - 2*x*p = 0
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (3/2 * log(x^2+1) + sqrt(9/4 * log(x^2+1)^2 - 1/(x^2+1)^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
dp = function(x) {
	x.mult = (x^2+1)^3
	p.x = p(x)
	div = (x^2+1)*(x.mult*p.x^6 - 1)
	dp = 2*x*(x.mult*p.x^4 - p.x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=-5, to=5)
sapply(c(-(1:4), 1:4), line.tan, dx=3)




####################
####################
####################


####################

### F(D1): dz + f(z)

### y^3 - 3*c*y - 2*f = 0
# [not run]
c*y*d2y + f*d2y + c*dy^2 + 2/3*df*dy - 1/2*d2c*y^2 - 1/3*d2f*y # = 0

### Special Cases:
# c = constant; f = x^3 - 1
x^2*(x^3-1)*d2y - c^2*dy^3 - c*x^2*dy^2 + 2*x*dy # = 0
# z = dy; dz = d2y;
x^2*(x^3-1)*dz - c^2*z^3 - c*x^2*z^2 + 2*x*z # = 0

### Solution & Plot:
# y = (f + sqrt(f^2 - c^3))^(1/3) + (f - sqrt(f^2 - c^3))^(1/3)
y = function(x, h, n=3) {
	f.x = x^3 - 1
	det = sqrt(f.x^2 - h^3 + 0i)
	r1 = round0(rootn(f.x + det, n=n))
	r2 = round0(rootn(f.x - det, n=n))
	return( round0(r1 + r2) )
}
dy = function(x, h, y.x) {
	if(missing(y.x)) y.x = y(x, h);
	div = (h*y.x + x^3 - 1)
	dp = x^2*y.x;
	dp = ifelse(div != 0, dp / div, Inf); # TODO
	return(dp)
}
d2y = function(x, h) {
	dy.x = dy(x, h);
	div = x^2*(x^3 - 1)
	dp = h^2*dy.x^3 + h*x^2*dy.x^2 - 2*x*dy.x;
	dp = ifelse(div != 0, dp / div, Inf); # TODO
	return(dp)
}
### Plot:
c = 1;
# a small local minimum
curve(y(x, c), from=-3, to=3)
line.tan(c(-2, c(-(2:1), 1:2)/1.7, 2), dx=3, p=y, dp=dy, h=c)
# D(y)
curve(dy(x, c), from=-3, to=3, add=T, col="green")
line.tan(c(-2, c(-(2:1), 1:2)/1.7, 2), dx=3, p=dy, dp=d2y, h=c, col="orange")


#######################

### (x^n - 1)*dy = f(x);
# [simple Model]
# x*(x^n - 1)*d2y + n*dy = x*df - n*f;
x*(x^n - 1)*d2y + b*(x^n - 1)*y*dy + n*dy - b*f*y - x*df + n*f # = 0
# where b = a constant or a function;
# e.g. b = x; f = b0 (constant);
x*(x^n - 1)*d2y + x*(x^n - 1)*y*dy + n*dy - b0*x*y + n*b0 # = 0

### Solution & Plot:
y = function(x, n=3, b0=1) {
	f.x = x^n - 1
	y.f = function(x) b0/(x^n - 1)
	y = sapply(x, function(x) integrate(y.f, lower=1 + 1E-5, upper=x)$value);
	return( y )
}
dy = function(x, n, b0=1) {
	div = x^n - 1
	dp = b0; # f = b0;
	dp = ifelse(div != 0, dp / div, Inf); # TODO
	return(dp)
}
d2y = function(x, n, b0=1) {
	# uses simple equation;
	dy.x = dy(x, n=n, b0=b0);
	div = x*(x^n - 1);
	dp = - n * (dy.x + b0);
	dp = ifelse(div != 0, dp / div, Inf); # TODO
	return(dp)
}
d2ymixt = function(x, n, b0=1) {
	# checking the mixt equation;
	y.x = y(x, n, b0)
	dy.x = dy(x, n=n, b0=b0);
	xn = x^n - 1
	dp = - (x*xn*y.x*dy.x + n*dy.x - b0*x*y.x + n*b0);
	div = x*xn;
	dp = ifelse(div != 0, dp / div, Inf); # TODO
}
#
c = 1; n = 5;
#
curve(y(x, n, c), from=1 + 1E-5, to=3)
line.tan(c(1.1, 2:5/1.7), dx=3, p=y, dp=dy, n=n, b0=c)
# D(y)
curve(dy(x, n, c), add=T, col="green")
# line.tan(c(1.1, 2:5/1.7), dx=3, p=dy, dp=d2y, n=n, b0=c, col="orange")
line.tan(c(1.1, 2:5/1.7), dx=3, p=dy, dp=d2ymixt, n=n, b0=c, col="orange")


#########################

### Integration by parts:

### Base ODE:
h*y*dy + f*dy - 1/2*dh*y^2 - 1/3*df*y # = 0

### Case: h = c (constant);
c*y*dy + f*dy - 1/3*df*y # = 0

### Ex 1:
# f = b1*x + b0;
c*y*dy + (b1*x+b0)*dy - b1/3 * y # = 0
# I() =>
c/2*y^2 + (b1*x+b0)*y - 4/3*b1*I(y) # = 0
# z = I(y); dz = y;
# c/2*dz^2 + (3*x+b0)*dz - 4/3*b1*z = 0;
### Variant:
c/3*y^3 + 1/2*(b1*x+b0)*y^2 - 2/3*b1*I(y^2) # = 0

### Solution & Plot:
# y = (f + sqrt(f^2 - h^3))^(1/3) + (f - sqrt(f^2 - h^3))^(1/3)
y = function(x, b, h, n=3) {
	d = (b[1]*x + b[2]);
	det = sqrt(d^2 - h^n + 0i)
	r1 = (d + det)
	r2 = (d - det)
	r = round0(rootn(r1, n=n) + rootn(r2, n=n))
	return(r)
}
Iy = function(x, b, h, n=3) {
	yx = y(x, b=b, h=h, n=n)
	yI = (h/2*yx^2 + (b[1]*x+b[2])*yx) / 4 * 3/b[1]
	return(yI)
}
Iy2 = function(x, b, h, n=3) {
	yx = y(x, b=b, h=h, n=n)
	yI = (h/3*yx^3 + 1/2*(b[1]*x+b[2])*yx^2) / 2 * 3/b[1]
	return(yI)
}
dI2 = function(x, b, h, n=3) {
	yx = y(x, b=b, h=h, n=n)
	return(yx^2)
}
### Plot
b = c(3, 1); h = 1;
#
curve(Iy(x, b=b, h=h), from=-3, to=3)
line.tan(c(-2:2), dx=3, p=Iy, dp=y, b=b, h=h)

#
curve(Iy2(x, b=b, h=h), from=-3, to=3)
line.tan(c(-2:2), dx=3, p=Iy2, dp=dI2, b=b, h=h)

###
b = c(3, -1/3); h = 1;
#
curve(Iy(x, b=b, h=h), from=-3, to=3)
line.tan(c(-2:2), dx=3, p=Iy, dp=y, b=b, h=h)

###
b = c(3, -2); h = 2;
#
curve(Iy(x, b=b, h=h), from=-3, to=3)
line.tan(c(-2:2), dx=3, p=Iy, dp=y, b=b, h=h)

###
b = c(2, 2); h = -3;
#
curve(Iy(x, b=b, h=h), from=-3, to=3)
line.tan(c(-2:2), dx=3, p=Iy, dp=y, b=b, h=h)

###
b = c(-2, 2); h = -3;
#
curve(Iy(x, b=b, h=h), from=-3, to=3)
line.tan(c(-2:2), dx=3, p=Iy, dp=y, b=b, h=h)


########################
########################

##################
### y-Radicals ###
##################

### Initial "Polynomial"
# y - G(x)*y^(1/n) = F(x)

### D =>
dy - 1/n*g*y^(1/n-1)*dy - dg*y^(1/n) - df # = 0
g*y*dy - 1/n*g*(y - f)*dy - dg*y*(y - f) - g*df*y # = 0
### ODE:
(n-1)*g*y*dy + g*f*dy - n*dg*y^2 + n*(dg*f - g*df)*y # = 0

### Example:
# g = x^2;
# f = x^2;
(n-1)*x*y*dy + x^3*dy - 2*n*y^2 # = 0

### Solution & Plot:
y = function(x, n=3) {
	coeff0 = c(1, rep(0, n-2));
	r = sapply(x, function(x) {
		coeff = c(coeff0, -x^2, -x^2);
		r = roots(coeff);
		r = r[Im(r) == 0][1];
		return(r);
	})
	r = r^n;
	return(r)
}
dy = function(x, n=3) {
	# (n-1)*x*y*dy + x^3*dy - 2*n*y^2 # = 0
	yx = y(x, n=n); x3 = x^3;
	dy = 2*n*yx^2;
	div = (n-1)*x*yx + x3;
	dy = dy / div;
	return(dy)
}
### Plot
n = 3;
px = -3:3 * 3/5;
curve(y(x, n=n), from=-2, to=2)
line.tan(px, p=y, dp=dy, n=n)


####################

