########################
##
## Leonard Mada
## [the one and only]
##
## Exact Integration:
##   Polynomial Fractions:
##   Mixed Radicals
##
## draft v.0.1a


### Mixed Radicals

### Examples
# I( x^2 * sqrt(1 - x^2) / (x^4 + 1) )
# I( x^2 * sqrt(1 + x^2) / (x^4 + 1) )


### History

# - moved "Mixed Radicals" to this file;
# TODO: move the remaining section;


########################
########################

### I( sqrt(x * (1 - x)) / (x^2 + 1) ) on [0, 1]
# Maths 505: Why this integral is MUCH harder than it looks
# https://www.youtube.com/watch?v=vPWTmwjYM8s
# Euler subst: sqrt(x - x^2) = x*y => x = 1/(1 + y^2);

### I( x^2 * sqrt(1 - x^2) / (x^4 + 1) )
integrate(\(x) 2 * x^2 * sqrt(1 - x^2) / (x^4 + 1), 0, 1)
integrate(\(x) sqrt(x * (1 - x)) / (x^2 + 1), 0, 1)
integrate(\(x) sqrt(tan(x) * (1 - tan(x))), 0, pi/4); # initial I;
- pi + pi * 2^(1/4) * sin(3*pi/8)

### I( x^2 * sqrt(1 + x^2) / (x^4 + 1) ) on [0, 1]
integrate(\(x) x^2 * sqrt(1 + x^2) / (x^4 + 1), 0, 1)
integrate(\(x) 1/2 * sqrt(x + 1) / x / (x^2 + 1), 1, Inf)
integrate(\(x) sin(x)/(1 - cos(x)^2) / cos(x)^2 / (tan(x)^4 + 1), pi/4, pi/2)
integrate(\(x) x^2 / (1 - x^2) / (2*x^4 - 2*x^2 + 1), 0, 1/sqrt(2))
integrate(\(x) x^2 / (x^2 - 1) / (x^4 - 2*x^2 + 2), sqrt(2), Inf)
integrate(\(x) 1/2 * Im((1+1i)/(x^2 - 1) - (1+1i)/(x^2 - 1 + 1i) +
	- (1-1i)/(x^2 - 1) + (1-1i)/(x^2 - 1 - 1i)), sqrt(2), Inf)
integrate(\(x) 1/(x^2 - 1) + 1/2 * Im(
	- (1+1i)/(x^2 - 1 + 1i) +
	+ (1-1i)/(x^2 - 1 - 1i) ), sqrt(2), Inf)
1/2 * log((sqrt(2)+1)/(sqrt(2)-1)) - 1/4 * Im(
	(1+1i)/sqrt(1-1i) * log((sqrt(2) + sqrt(1-1i))/(sqrt(2) - sqrt(1-1i))) +
	+ (1-1i)/sqrt(1+1i) * log((sqrt(2) - sqrt(1+1i))/(sqrt(2) + sqrt(1+1i))));


# Derivation:
integrate(\(x) sqrt(x * (1 - x)) / (x^2 + 1), 0, 1)
integrate(\(x) 2*x^2 / ((x^2 + 1)*(x^4 + 2*x^2 + 2)), 0, Inf)
integrate(\(x) 1 / (x^2 + 1) *
	Re((1-1i)/(x^2 + 1 + 1i) + (1+1i)/(x^2 + 1 - 1i)), 0, Inf)
integrate(\(x) Im((1-1i)/(x^2 + 1) - (1-1i)/(x^2 + 1 + 1i) +
	- (1+1i)/(x^2 + 1) + (1+1i)/(x^2 + 1 - 1i)), 0, Inf)
integrate(\(x) -2/(x^2 + 1) +
	- Im((1-1i)/(x^2 + 1 + 1i) - (1+1i)/(x^2 + 1 - 1i)), 0, Inf)
- pi - pi/2 * Im((1-1i)/sqrt(1+1i) - ((1+1i)/sqrt(1-1i)))
- pi - pi/2 * 2^(1/4) * Im(exp(-3i*pi/8) - exp(3i*pi/8))
- pi + pi * 2^(1/4) * sin(3*pi/8)


### I( x^2 * sqrt(1 - x^2) / (x^4 + 1) ) on [0, 1/2]
integrate(\(x) x^2 * sqrt(1 - x^2) / (x^4 + 1), 0, 1/2)
integrate(\(x) 1/2 * sqrt(x * (1 - x)) / (x^2 + 1), 0, 1/4)
integrate(\(x) x^2 / ((x^2 + 1)*(x^4 + 2*x^2 + 2)), sqrt(3), Inf)
- pi/2 + pi * 2^(-3/4) * sin(3*pi/8) + atan(sqrt(3)) +
	+ 2^(-3/4) * Im(exp(-3i*pi/8) * atan(sqrt(3)/sqrt(1+1i)) +
	- exp(3i*pi/8) * atan(sqrt(3)/sqrt(1-1i)));


### I( sqrt(x * (1 - x)) / (x^2 + 1) ) on [0, 1/2]
integrate(\(x) sqrt(x * (1 - x)) / (x^2 + 1), 0, 1/2)
integrate(\(x) sqrt(tan(x) * (1 - tan(x))), 0, atan(1/2)) # as initial I;
- pi/2 + pi * 2^(1/4) * sin(3*pi/8) +
	+ Im((1-1i)/sqrt(1+1i)*atan(1/sqrt(1+1i)) +
	- (1+1i)/sqrt(1-1i)*atan(1/sqrt(1-1i)));

###
integrate(\(x) sqrt(tan(x) * (1 - tan(x))), atan(1/2), pi/4)
integrate(\(x) sqrt(x * (1 - x)) / (x^2 + 1), 1/2, 1)
- pi/2 - Im((1-1i)/sqrt(1+1i)*atan(1/sqrt(1+1i)) +
	- (1+1i)/sqrt(1-1i)*atan(1/sqrt(1-1i)))
- pi/2 - 2^(1/4) * Im(exp(-3i*pi/8) * atan(1/sqrt(1+1i)) +
	- exp(3i*pi/8) * atan(1/sqrt(1-1i)))
- pi/2 + pi * 2^(1/4) * sin(3*pi/8) +
	- 2^(-3/4) * exp(-3i*pi/8) * log((1 + 1i*sqrt(1+1i)) / (1 - 1i*sqrt(1+1i))) +
	+ 2^(-3/4) * exp(3i*pi/8) * log((1 + 1i*sqrt(1-1i)) / (1 - 1i*sqrt(1-1i)));

