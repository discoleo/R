########################
##
## Leonard Mada
## [the one and only]
##
## Exact Integration
## Polynomial Radicals
## Based on Tan^(1/n)
##
## draft v.0.1a


### Radicals: Tan-Based

### Simple Transforms

# Note: I( tan(x)^(p/q) )
# - for (p,q) integers: can be computed exactly on arbitrary intervals;
# - see section "Trigonometric" in file:
#   Integrals.Fractions.Unity.Derived.R;
# - Other Transforms based on TAN are in file:
#   Integrals.Fractions.Unity.Radicals.R;


######################

### Examples:

### I( x / (1-x^5)^(2/5) )
lim = acos(1/2^(5/2)); # gives [0, 1/2]
integrate(\(x) x * (1-x^5)^(-2/5), 0, cos(lim)^(2/5))
integrate(\(x) 2/5 * (1-x^2)^(-2/5) / x^(1/5), 0, cos(lim))
integrate(\(x) 2/5 * tan(x)^(1/5), lim, pi/2)
integrate(\(x) 2 * x^5 / (x^10 + 1), tan(lim)^(1/5), Inf)
integrate(\(x) x^2 / (x^5 + 1), tan(lim)^(2/5), Inf)
# Exact formula available;

# on [0, 1/2]
integrate(\(x) x * (1-x^5)^(-2/5), 0, 1/2)
integrate(\(x) 2/5 * tan(x)^(1/5), acos(1/2^(5/2)), pi/2)

