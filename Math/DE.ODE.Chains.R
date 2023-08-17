

###################
### ODE: Chains ###
###################


### Fundamental Solution:
### d2y = y

### Mult: dy =>
# dy * d2y = y * dy
# 1/2 * (dy)^2 = 1/2 * y^2 + 1/2 * c
# (dy)^2 = y^2 + c;


##############

##############
### Chains ###
##############

### Exp

# y = y0 + b/2 * x * exp(x)
# y0 satisfies the fundamental solution
# =>
# d2y = d2y0 + b/2 * x * exp(x) + b * exp(x)
# d2y = y0 + b/2 * x * exp(x) + b * exp(x)
# =>
# ODE: d2y = y + b * exp(x);

# TODO: check & expand;
