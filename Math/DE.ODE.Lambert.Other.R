########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs: Lambert - Other
###
### draft v.0.1b-ext2


### Lambert W-Like Equations

###############
### History ###
###############


### draft v.0.1b - v.0.1b-ext2:
# - ODE:
#   y*dy - dy - y = 0;
#   x*y*dy - x*dy - x*y + n*y = 0;
### draft v.0.1a:
# - derived from:
#   exp(y) + exp(y^2) = F(x);


#########################

### Helper functions

library(pracma)
# needed for Lambert W;


# include: Polynomials.Helper.R;
# include: DE.ODE.Helper.R;
source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#########################

### exp(y) = exp(x)*y

### D =>
exp(y)*dy - exp(x)*(dy + y) # = 0
# Subst =>
exp(x)*y*dy - exp(x)*(dy + y) # = 0

### ODE:
y*dy - (dy + y) # = 0

# TODO: check;


### Solution: y
x = 3
y = - lambertWp(- exp(-x))

# Test
exp(y) - exp(x)*y


###############
### Extensions:

### x^n * exp(y) = exp(x)*y

### D =>
exp(y)*(x^n*dy + n*x^(n-1)) - exp(x)*(dy + y) # = 0 # * x =>
x^n * exp(y)*(x*dy + n) - exp(x)*x*(dy + y) # = 0
# Subst =>
exp(x)*y*(x*dy + n) - exp(x)*x*(dy + y) # = 0

### ODE:
x*y*dy - x*dy - x*y + n*y # = 0

# TODO: check;

### Solution: y

n = sqrt(2)
x = 3
y = - lambertWp( - x^n * exp(-x) )
# Test
x^n * exp(y) - exp(x)*y


################
### Extension 2:

### x^n * exp(y) = exp(x)*y + x^m*exp(x)

### D =>
exp(y)*(x^n*dy + n*x^(n-1)) - exp(x)*(dy + y) - exp(x)*(x^m + m*x^(m-1)) # = 0
# * x =>
x^n * exp(y)*(x*dy + n) - exp(x)*x*(dy + y) - x^m * exp(x)*(x + m) # = 0
# Subst =>
exp(x)*(y + x^m)*(x*dy + n) - exp(x)*x*(dy + y) - x^m * exp(x)*(x + m) # = 0

### ODE:
x*y*dy + x*(x^m - 1)*dy - x*y + n*y - x^(m+1) + (n-m)*x^m # = 0

# TODO: check;


###########################

### exp(y) + exp(y^2) = F(x)

### D =>
exp(y)*dy + 2*y*exp(y^2)*dy - df # = 0

### Linear System:
# (2*y - 1)*exp(y^2)*dy = (df - f*dy);
# (2*y - 1)*exp(y)*dy = (2*f*dy - df);

### D2 =>
(exp(y) + 2*y*exp(y^2))*d2y +
	+ 4*y^2*exp(y^2)*dy^2 + exp(y)*dy^2 + 2*exp(y^2)*dy^2 - d2f # = 0
df*d2y + (4*y^2*exp(y^2) + exp(y) + 2*exp(y^2))*dy^3 - d2f*dy # = 0
df*(2*y - 1)*d2y + (4*y^2*(df - f*dy) + 2*f*dy - df + 2*(df - f*dy))*dy^2 +
	- d2f*(2*y - 1)*dy # = 0
df*(2*y - 1)*d2y - 4*f*y^2*dy^3 + 4*df*y^2*dy^2 + df*dy^2 - 2*d2f*y*dy + d2f*dy # = 0

# TODO;


#########################

### exp(y) + exp(1/y) = F(x)

### D =>
y^2*exp(y)*dy - exp(1/y)*dy - df*y^2 # = 0

### Linear System:
# (y^2 + 1)*exp(y^2)*dy = (f*y^2*dy - df*y^2);
# (y^2 + 1)*exp(y)*dy = (f*dy + df*y^2);

### D2 =>
# TODO;


#########################

### Generalization: Coefficients
### P1(x)*exp(G1(x)*y) + P2(x)*exp(G2(x)/y) = F(x)

# TODO;

