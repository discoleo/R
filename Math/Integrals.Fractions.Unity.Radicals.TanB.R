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


### Diff-Type

### I( (1 + x^5)^(3/5) / (1-x^10)^(4/5) / x )
lim = 1/2;
integrate(\(x) (1 + x^5)^(3/5) / (1-x^10)^(4/5) / x, lim, 1)
integrate(\(x) 1 / (1 - x^5), 0, tan(acos(lim^5)/2)^(2/5))


# Derivation:
lim = 1/3
integrate(\(x) 1 / (1 - x^5), 0, lim)
integrate(\(x) 2/5 * tan(x)^(-3/5) / (cos(x)^2 - sin(x)^2), 0, atan(lim^(5/2)))
integrate(\(x) 2^(8/5)/5 * cos(x)^(6/5) / sin(2*x)^(3/5) / cos(2*x), 0, atan(lim^(5/2)))
integrate(\(x) 1/5 * (cos(x) + 1)^(3/5) / sin(x)^(3/5) / cos(x), 0, 2*atan(lim^(5/2)))
integrate(\(x) (1 + x^5)^(3/5) / (1-x^10)^(4/5) / x, cos(2*atan(lim^(5/2)))^(1/5), 1)

