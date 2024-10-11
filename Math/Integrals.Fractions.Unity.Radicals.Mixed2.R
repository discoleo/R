########################
##
## Leonard Mada
## [the one and only]
##
## Exact Integration:
##   Polynomial Fractions:
##   Mixed Radicals: 2 Radicals
##
## draft v.0.1a


### Mixed Radicals: 2 Radicals



### I( x^4 * (1 - x^4)^(3/4) / (x^2 * (1 - x^4)^(1/2) + 1) )
integrate(\(x) x^4 * (1 - x^4)^(3/4) / (x^2 * (1 - x^4)^(1/2) + 1), 0, 1)
1/4 * (beta(3/4, 2-3/4) - beta(1/4, 3/4)) + pi / tan(pi/6) / 6;


# Derivation:
integrate(\(x) x^4 * (1 - x^4)^(3/4) / (x^2 * (1 - x^4)^(1/2) + 1), 0, 1)
integrate(\(x) 1/4 * (x*(1-x)^3)^(1/4) / ((x*(1-x))^(1/2) + 1), 0, 1)
integrate(\(x) x^4 / (x^4 + 1)^3 / (x^2/(x^4+1) + 1), 0, Inf)
integrate(\(x) x^4 / (x^4 + 1)^2 / (x^4 + x^2 + 1), 0, Inf)
integrate(\(x) x^2 / (x^4 + 1)^2 - 1/(x^4+1) + 1/(x^4 + x^2 + 1), 0, Inf)
1/4 * (beta(3/4, 2-3/4) - beta(1/4, 3/4)) + pi / tan(pi/6) / 6;


#
integrate(\(x) 1/(x^4 + x^2 + 1), 0, Inf)
- pi / tan(3*pi/6) / 6 + pi / tan(pi/6) / 6;
pi / tan(pi/6) / 6

