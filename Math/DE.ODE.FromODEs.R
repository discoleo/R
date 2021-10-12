########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs: From other ODEs
###
### draft v.0.1a



###############
### History ###
###############


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

# TODO


####################
### Radicals:    ###
### Higher Order ###
####################

### dy - G(x)*y^(1/n) = F(x)

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

### y*dy - G(x)*y = F(x)
# [very simple eq]

### D =>
y*d2y + dy^2 - g*dy - dg*y - df # = 0

### Variant 1: * y =>
y^2*d2y + y*dy^2 - g*y*dy - dg*y^2 - df*y # = 0
y^2*d2y + (g*y + f)*dy - g*(g*y + f) - dg*y^2 - df*y # = 0
### ODE:
y^2*d2y + g*y*dy + f*dy - dg*y^2 - df*y - g^2*y - g*f # = 0

# Subst 2 =>
y^2*d2y + g*(g*y + f) + f*dy - dg*y^2 - df*y - g^2*y - g*f # = 0
### ODE:
y^2*d2y + f*dy - dg*y^2 - df*y # = 0
# => original Eq:
# d2y - dg - D(f/y) # = 0

### Check:
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


