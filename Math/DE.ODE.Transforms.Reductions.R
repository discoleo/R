########################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Transforms
## Order Reduction
##
## v.0.1c


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


### Transform: Exp(k*x^2)
# y = z * exp(k*x^2)
# dy = (dz + 2*k*x*z) * exp(k*x^2)
# d2y = (d2z + 4*k*x*dz + 4*k^2*x^2*z) * exp(k*x^2)

### Order Reduction
# =>
# d2z + 4*k*x*dz = f0 * exp(-k*x^2)


###################

### Gen: x^n * y

### Base ODE:
# d2y = (n*k)^2 * x^(2*n-2)*y + f0


### Transform: Exp(k*x^n)
# y = z * exp(k*x^n)
# dy = (dz + n*k*x^(n-1)*z) * exp(k*x^n)
# d2y = (d2z + 2*n*k*x^(n-1)*dz + n^2*k^2*x^(2*n-2)*z) * exp(k*x^n)

### Order Reduction
# =>
# d2z + 2*n*k*x^(n-1)*dz = f0 * exp(-k*x^n)

