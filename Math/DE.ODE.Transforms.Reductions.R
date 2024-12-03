########################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Transforms
## Order Reduction
##
## v.0.1b


### Order Reduction of ODEs


###################

### Base ODE:
# d2y = y + f0


### Transform: Exp(x)
# y = z * exp(x)
# dy = (dz + z) * exp(x)
# d2y = (d2z + 2*dz + z) * exp(x)

### Order Reduction
# =>
# d2z + 2*dz = f0 * exp(-x)



### Transform: Exp(-x)
# y = z * exp(-x)
# dy = (dz - z) * exp(-x)
# d2y = (d2z - 2*dz + z) * exp(-x)

### Order Reduction
# =>
# d2z - 2*dz = f0 * exp(x)


###################
###################

### Base ODE:
# d2y = 4*k^2*x^2*y + f0


### Transform: Exp(x^2)
# y = z * exp(x^2)
# dy = (dz + 2*k*x*z) * exp(k*x^2)
# d2y = (d2z + 4*k*x*dz + 4*k^2*x^2*z) * exp(k*x^2)

### Order Reduction
# =>
# d2z + 4*k*x*dz = f0 * exp(-k*x^2)

