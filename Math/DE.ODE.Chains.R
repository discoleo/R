

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


#############

### d2y = y + sum(b[i] * exp(k[i] * x))

# y = y0 + b1/(k1^2 - 1) * exp(k1*x) + ...
# where k[i] != 1;
# y0 satisfies the fundamental solution;

# =>
# d2y = d2y0 + b1*k1^2/(k1^2 - 1) * exp(k1*x) + ...
# d2y = y0 + b1/(k1^2 - 1) * exp(k1*x) + b1*exp(k1*x) + ...
# ODE: d2y = y + b1*exp(k1*x) + ...

### Examples

# TODO

