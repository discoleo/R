

########################

### I( log(x) / (x + 1)^p ) on [0, 1]

# TODO: truncated Polylog function;

###
integrate(function(x) log(x)/(x+1)^2, 0, 1)
- log(2)

###
integrate(function(x) log(x)/(x+1)^3, 0, 1)
- log(2)/2 + (1/2 - 1)/2

###
integrate(function(x) log(x)/(x+1)^4, 0, 1)
- log(2)/3 + (1/2 + 1/4/2 - 1 - 1/2)/3

###
integrate(function(x) log(x)/(x+1)^5, 0, 1)
- log(2)/4 + (1/2 + 1/4/2 + 1/8/3 - 1 - 1/2 - 1/3)/4

###
integrate(function(x) log(x)/(x+1)^6, 0, 1)
- log(2)/5 + (1/2 + 1/4/2 + 1/8/3 + 1/16/4 - 1 - 1/2 - 1/3 - 1/4)/5


### I( log(x) * sqrt(x + 1) )
integrate(\(x) log(x) * sqrt(x + 1), 0, 1)
- 2/3 * (2*log(2) - 2*asinh(1) + 2*sqrt(2) - 2 + 2/3 * (2^(3/2) - 1))
4/3 * (asinh(1) - log(2) - (1 + 2/3)*sqrt(2) + 1 + 1/3)

# Derivation:
integrate(\(x) (sqrt(x + 1) - 1)/x, 0, 1)
2*log(2) - 2*asinh(1) + 2*sqrt(2) - 2
#
integrate(\(x) ((x+1)*sqrt(x + 1) - 1)/x, 0, 1)
2*log(2) - 2*asinh(1) + 2*sqrt(2) - 2 + 2/3 * (2^(3/2) - 1)


### I( log(x) / sqrt(x + 1) )
integrate(\(x) log(x) / sqrt(x + 1), 0, 1)
4 * (asinh(1) - log(2) - sqrt(2) + 1)

