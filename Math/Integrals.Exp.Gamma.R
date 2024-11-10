###############
##
## Integrals: Exp
## Gamma-Variants
##
## Leonard Mada


#################

### on [0,1]

# TODO: find closed formula;


### I( exp(-x^n / n) )
# Mathologer: Ramanujan's easiest hard infinity monster
# (Mathologer Masterclass)
# https://www.youtube.com/watch?v=6iTdNmDHfV0
# Note: very fast-converging way to compute these integrals;


###
integrate(\(x) exp(-x^2/2), 0, 1)
id = seq(1, 31, by=2)
sum(1/cumprod(id)) / exp(1/2)


###
integrate(\(x) exp(-x^3/3), 0, 1)
id = seq(1, 31, by=3)
sum(1/cumprod(id)) / exp(1/3)


###
n = sqrt(3)
integrate(\(x) exp(-x^n/n), 0, 1)
id = seq(1, 31, by=n)
sum(1/cumprod(id)) / exp(1/n)

