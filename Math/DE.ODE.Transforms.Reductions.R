############
##
## Leonard Mada
## (the one and only)
##
## Integrals: Transforms
## Order Reduction
##
## v.0.1a


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

