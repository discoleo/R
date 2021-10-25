########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs: From other ODEs
###
### draft v.0.2e



###############
### History ###
###############


### draft v.0.2e:
# - derived from:
#   dy = exp(G(x)) * (P2(x)*y^2 + P1(x)*y + P0(x));
### draft v.0.2d:
# - derived from a Riccati equation;
# - various variants, e.g.:
#   x*d2y - 2*x^3*y*dy - (x^3 + 2)*dy + 2*k = 0;
### draft v.0.2c:
# - derived from:
#   dy - G(x)*log(y + k) = F(x);
### draft v.0.2b:
# - various examples based on y-Exponentials, eg:
#   x*d2y - x*dy^2 + n*dy = 0;
#   x*d2y - x*dy^2 + (b0*x + n)*dy - n*b0 = 0; [Non-Hom/v.0.2c]
### draft v.0.2a:
# - moved to this file Section with Exponentials:
#   from DE.ODE.Exponentials.R;
### draft v.0.1e:
# - more variants for:
#   y*dy - G(x)*y = F(x);
### draft v.0.1c - v.0.1d:
# - derived from:
#   dy^2 - G(x)*y = F(x);
#   dy - G(x)*log(y) = F(x); [v.0.1d]
### draft v.0.1b:
# - derived from ODEs with higher order polynomials;
### draft v.0.1a:
# - [started] move relevant sections from
#   other files: DE.ODE.Polynomial.R;


####################
####################

### helper functions

# include: DE.ODE.Helper.R;
source("DE.ODE.Helper.R")


########################

########################
### Derived from ODE ###
########################

####################
### Higher Order ###
####################

### dy - G(x)*y^n = F(x)

### D =>
d2y - n*g*y^(n-1)*dy - dg*y^n - df # = 0 # * y =>
y*d2y - n*g*y^(n)*dy - dg*y^(n+1) - df*y # = 0
g*y*d2y - n*g*(dy - f)*dy - dg*(dy - f)*y - g*df*y # = 0
### ODE:
g*y*d2y - n*g*dy^2 - dg*y*dy + n*g*f*dy + (dg*f - g*df)*y # = 0


####################

####################
### Radicals:    ###
### Higher Order ###
####################

### dy - G(x)*y^(1/n) = F(x)
# - similar to the Example with higher powers;

### D =>
d2y - 1/n*g*y^(1/n-1)*dy - dg*y^(1/n) - df # = 0

### Variant 1:
g*y*d2y - 1/n*g*(dy - f)*dy - dg*y*(dy - f) - g*df*y # = 0
g*y*d2y - 1/n*g*dy^2 - dg*y*dy + 1/n*g*f*dy + (dg*f - g*df)*y # = 0

### Example 1:
# f = x; g = x;
x*y*d2y - 1/n*x*dy^2 - y*dy + 1/n*x^2*dy # = 0

### Example 2:
# f = x; g = n;
n*y*d2y - dy^2 + x*dy - n*y # = 0

# TODO: check;
# + concept to check;


########################
########################\

### Riccati Eq:
### dy = p2*y^2 + p1*y + p0

### D =>
d2y - 2*p2*y*dy - p1*dy - dp2*y^2 - dp1*y - dp0 # = 0

### Variant 1:
d2y - 2*p2*y*(p2*y^2 + p1*y + p0) - p1*(p2*y^2 + p1*y + p0) - dp2*y^2 - dp1*y - dp0 # = 0
### ODE:
d2y - 2*p2^2*y^3 - (3*p1*p2 + dp2)*y^2 - (p1^2 + 2*p0*p2 + dp1)*y - p0*p1 - dp0 # = 0

### Variant 2:
p2*d2y - 2*p2^2*y*dy - p1*p2*dy - dp2*(dy - p1*y - p0) - dp1*p2*y - dp0*p2 # = 0
### ODE:
p2*d2y - 2*p2^2*y*dy - (p1*p2 + dp2)*dy + (dp2*p1 - dp1*p2)*y + dp2*p0 - dp0*p2 # = 0

### Solution:
# C2*d2y - C11*y*dy - C10*dy + B1*y + B0 = 0
# 2*(H*C2)^2 = H*C11 => H = C11 / (2*C2^2);

# TODO

### Special Cases:
# p2 = x^2; p0 = k; dp2 = 2*x;
x*d2y - 2*x^3*y*dy - (p1*x + 2)*dy + (2*p1 - dp1*x)*y + 2*k # = 0
# p1 = x^2 =>
x*d2y - 2*x^3*y*dy - (x^3 + 2)*dy + 2*k # = 0


### Generalized Riccati Eq:
### dy = p2*y^n + p1*y + p0

### D =>
d2y - n*p2*y^(n-1)*dy - p1*dy - dp2*y^n - dp1*y - dp0 # = 0 # * y =>
y*d2y - n*p2*y^n*dy - p1*y*dy - dp2*y^(n+1) - dp1*y^2 - dp0*y # = 0

### Variant 3:
p2*y*d2y - n*p2*(dy - p1*y - p0)*dy - p1*p2*y*dy - dp2*y*(dy - p1*y - p0) - dp1*p2*y^2 - dp0*p2*y # = 0


########################
########################

### y*dy Types:

### y*dy - G(x)*y = F(x)
# [very simple eq]

### D =>
y*d2y + dy^2 - g*dy - dg*y - df # = 0

### Variant 1:
# D * y =>
y^2*d2y + y*dy^2 - g*y*dy - dg*y^2 - df*y # = 0
y^2*d2y + (g*y + f)*dy - g*(g*y + f) - dg*y^2 - df*y # = 0
### ODE:
y^2*d2y + g*y*dy + f*dy - dg*y^2 - df*y - g^2*y - g*f # = 0

# Subst 2 =>
# [reducible/redundant]
y^2*d2y + g*(g*y + f) + f*dy - dg*y^2 - df*y - g^2*y - g*f # = 0
### ODE:
y^2*d2y + f*dy - dg*y^2 - df*y # = 0
# => original Eq:
# d2y - dg - D(f/y) # = 0


### Variant 2:
# D * dy =>
y*dy*d2y + dy^3 - g*dy^2 - dg*y*dy - df*dy # = 0
(g*y + f)*d2y + dy^3 - g*dy^2 - dg*(g*y + f) - df*dy # = 0
### ODE:
g*y*d2y + f*d2y + dy^3 - g*dy^2 - df*dy - dg*g*y - dg*f # = 0

### Debug/Check:
x = sqrt(3)
y  = x + log(x^2 + 1);
dy = 2*x/(x^2+1) + 1;
g  = 2*x/(x^2+1);
f  = x + log(x^2 + 1);
d2y = - 2*(x^2-1)/(x^2+1)^2;
dg  = d2y;
df  = 1 + 2*x/(x^2+1);


### [redundant]
### Polynomial Transformations

# * g =>
(y*dy - f)*dy - g^2*y - g*f # = 0
### Derived variant:
y*dy^2 - f*dy - g^2*y - g*f # = 0
(dy + g)*(y*dy - g*y - f) # = 0


#########################
#########################

######################
### y*d2y Variants ###
######################

### D of Higher Power

### dy^2 - G(x)*y = F(x)

### D =>
2*dy*d2y - g*dy - dg*y - df # = 0 # * dy =>
2*dy^2*d2y - g*dy^2 - dg*y*dy - df*dy # = 0
2*(g*y + f)*d2y - g*(g*y + f) - dg*y*dy - df*dy # = 0
### ODE:
2*g*y*d2y + 2*f*d2y - dg*y*dy - df*dy - g^2*y - g*f # = 0

# TODO: check;


########################
########################

### dy^n = y + F(x)

### D =>
n*dy^(n-1)*d2y - dy - df # = 0 # * dy
n*dy^n*d2y - (dy)^2 - df*dy # = 0
n*(y + f)*d2y - (dy)^2 - df*dy # = 0
### ODE:
n*y*d2y + n*f*d2y - (dy)^2 - df*dy # = 0

# TODO: check;


###################

### Generalization:
### dy^n - G(x)*y = F(x)

### D =>
n*dy^(n-1)*d2y - g*dy - dg*y - df # = 0 # * dy =>
n*dy^n*d2y - g*dy^2 - dg*y*dy - df*dy # = 0
n*(g*y + f)*d2y - g*dy^2 - dg*y*dy - df*dy # = 0
### ODE:
n*g*y*d2y + n*f*d2y - g*dy^2 - dg*y*dy - df*dy # = 0

# TODO: check;


#########################
#########################

### D of Logarithm

### dy - G(x)*log(y) = F(x)

### D =>
d2y - g/y * dy - dg*log(y) - df # = 0
g*y*d2y - g^2*dy - dg*g*log(y)*y - g*df*y # = 0
g*y*d2y - g^2*dy - dg*(dy - f)*y - g*df*y # = 0

### ODE:
g*y*d2y - dg*y*dy - g^2*dy + (dg*f - g*df)*y # = 0

# TODO: check;


#########################
### Slight Generalization
### dy - G(x)*log(y + k) = F(x)

### D =>
d2y - g/(y + k) * dy - dg*log(y + k) - df # = 0
g*(y+k)*d2y - g^2*dy - dg*g*log(y+k)*(y+k) - g*df*(y+k) # = 0
g*(y+k)*d2y - g^2*dy - dg*(dy - f)*(y+k) - g*df*(y+k) # = 0
g*(y+k)*d2y - g^2*dy - dg*(y*dy + k*dy - f*y - k*f) - g*df*(y+k) # = 0

### ODE:
g*(y+k)*d2y - dg*y*dy - (g^2 + k*dg)*dy + (dg*f - g*df)*y - k*(g*df - dg*f) # = 0

### Special Case:
# g = x + a; f = x + b;
(x+a)*(y+k)*d2y - y*dy - ((x+a)^2 + k)*dy + (b-a)*y + k*(b-a) # = 0


#########################
#########################

#########################
###   Derived from    ###
### Exponential ODEs  ###
#########################

### Simple:
### dy - G(x)*y = F(x)

### D =>
d2y - g*dy - dg*y - df # = 0

### Entanglements:

### Variant 1:
d2y - g*(g*y + f) - dg*y - df # = 0
d2y - (g^2 + dg)*y - g*f - df # = 0

### Variant 2:
y*d2y - (dy - f)*dy - dg*y^2 - df*y # = 0
y*d2y - dy^2 + f*dy - dg*y^2 - df*y # = 0

### Variant 3:
g*d2y - g^2*(g*y + f) - dg*(dy - f) - g*df # = 0
g*d2y - dg*dy - g^3*y - g^2*f + dg*f - g*df # = 0

### Variant 4:
y*d2y - (dy - f)*(g*y + f) - dg*y^2 - df*y # = 0
y*d2y - g*y*dy - f*dy - dg*y^2 + (g*f - df)*y + f^2 # = 0


### Examples:
# G(x) = exp(x^2) => V2 =>
y*d2y - dy^2 + f*dy - 2*x*exp(x^2)*y^2 - df*y # = 0
### Extra-Variant:
y*d2y - dy^2 + f*dy - 2*x*(dy - f)*y - df*y # = 0
y*d2y - dy^2 - 2*x*y*dy + f*dy + 2*x*f*y - df*y # = 0


### Example 1:
# F(x) = exp(x^3) + b0
# G(x) = -3*x^2
d2y - (9*x^4 - 6*x)*y + 3*b0*x^2 # = 0

### Example 2:
# G(x) = -1/x;
x^2*d2y + x*dy + y - x^2*df # = 0
# F(x) = ln(x^n + b0);
x^2*d2y + x*dy + y - n*x^(n+1)/(x^n + b0) # = 0


########################
########################

### Higher Powers:
### dy - G(x)*y^2 = F(x)
### see also Section: Riccati Eq;

### D =>
d2y - 2*g*y*dy - dg*y^2 - df # = 0

### Variant 1:
g*d2y - 2*g^2*y*dy - dg*(dy - f) - g*df # = 0
g*d2y - 2*g^2*y*dy - dg*dy + dg*f - g*df # = 0


### Other Powers:
### dy - G(x)*y^(1/n) = F(x)

### D =>
d2y - 1/n*g*y^(1/n - 1)*dy - dg*y^(1/n) - df # = 0

### Variant 1:
g*y*d2y - 1/n*g*(dy - f)*dy - dg*y*(dy - f) - g*df*y # = 0
g*y*d2y - 1/n*g*dy^2 - dg*y*dy + 1/n*g*f*dy + dg*f*y - g*df*y # = 0


##################
### Higher Powers:
### Generalization

### dy - G2(x)*y^2 - G1(x)*y = F(x)
### see also Section: Riccati Eq;

### D =>
d2y - 2*g2*y*dy - dg2*y^2 - g1*dy - dg1*y - df # = 0

### Variant 1:
g2*d2y - 2*g2^2*y*dy - dg2*(dy - g1*y - f) - g1*g2*dy - dg1*g2*y - g2*df # = 0
g2*d2y - 2*g2^2*y*dy - (dg2 + g1*g2)*dy + (dg2*g1 - dg1*g2)*y + dg2*f - g2*df # = 0


#########################

### Exponential-Variants
### of Riccati Equation

### dy = exp(G(x)) * (P2(x)*y^2 + P1(x)*y + P0(x))

### D =>
d2y - dg*dy - exp(g)*(2*p2*y*dy + p1*dy + dp2*y^2 + dp1*y + dp0) # = 0
(p2*y^2 + p1*y + p0)*d2y - dg*(p2*y^2 + p1*y + p0)*dy +
	- (2*p2*y*dy + p1*dy + dp2*y^2 + dp1*y + dp0)*dy # = 0
### ODE:
(p2*y^2 + p1*y + p0)*d2y - 2*p2*y*dy^2 - p1*dy^2 +
	- (dp2*y^2 + dg*p2*y^2 + dp1*y + dg*p1*y + dp0 + dg*p0)*dy # = 0


#########################
#########################

######################
### y-Exponentials ###
######################

### Example 1:
x*y*dy - exp(y) # = 0

### D =>
x*y*d2y + x*dy^2 + y*dy - exp(y)*dy # = 0
x*y*d2y - x*y*dy^2 + x*dy^2 + y*dy # = 0


###############
### Example 2a:
x^n*dy - exp(y) # = 0
# trivial: just for check;

### D =>
x^n*d2y + n*x^(n-1)*dy - exp(y)*dy # = 0
x^n*d2y + n*x^(n-1)*dy - x^n*dy^2 # = 0
x*d2y - x*dy^2 + n*dy # = 0


### Example 2a-2:
### Non-Homogeneous:
x^n*dy - exp(y) - x^n*f # = 0

### D =>
x^n*d2y + n*x^(n-1)*dy - exp(y)*dy - x^n*df - n*x^(n-1)*f # = 0
x^n*d2y + n*x^(n-1)*dy - x^n*dy^2 + x^n*f*dy - x^n*df - n*x^(n-1)*f # = 0
x*d2y - x*dy^2 + x*f*dy + n*dy - x*df - n*f # = 0

### Special Case:
# f = b0
x*d2y - x*dy^2 + (b0*x + n)*dy - n*b0 # = 0

### y = ln(z) =>
x^n*dz - z^2 - x^n*f*z # = 0


###############
### Example 2b:
x^n*dy - exp(x^2*y) # = 0

### D =>
x^n*d2y + n*x^(n-1)*dy - (x^2*dy + 2*x*y)*exp(x^2*y) # = 0
x^n*d2y + n*x^(n-1)*dy - (x^2*dy + 2*x*y)*x^n*dy # = 0
x*d2y + n*dy - (x^2*dy + 2*x*y)*x*dy # = 0
x*d2y - x^3*dy^2 - 2*x^2*y*dy + n*dy # = 0


##############
### Example 3:
x*y*dy + y^2 - exp(y) # = 0

### D =>
x*y*d2y + x*dy^2 + y*dy + 2*y*dy - exp(y)*dy # = 0
x*y*d2y + x*dy^2 + 3*y*dy - (x*y*dy + y^2)*dy # = 0
x*y*d2y - x*y*dy^2 + x*dy^2 - y^2*dy + 3*y*dy # = 0

# TODO: check!

