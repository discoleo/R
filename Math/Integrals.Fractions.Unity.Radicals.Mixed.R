########################
##
## Leonard Mada
## [the one and only]
##
## Exact Integration:
##   Polynomial Fractions:
##   Mixed Radicals
##
## draft v.0.1b


### Mixed Radicals

### Examples
# I( x^2 * sqrt(1 - x^2) / (x^4 + 1) )
# I( x^2 * sqrt(1 + x^2) / (x^4 + 1) )
# I( x^3 * (1 + x^3)^(2/3) / (x^6 + 1) )
# I( x^3 * (1 - x^3)^(2/3) / (x^6 + 1) )


### History

# - moved "Mixed Radicals" to this file;
# TODO: move the remaining section;


########################
########################

### I( x^2 * sqrt(1 + x^2) / (x^4 + 1) )
# Maths 505: Why this integral is MUCH harder than it looks
# https://www.youtube.com/watch?v=vPWTmwjYM8s
# Initial: I( sqrt(x * (1 - x)) / (x^2 + 1) ) on [0, 1]
# Euler subst: sqrt(x - x^2) = x*y => x = 1/(1 + y^2);

### I( x^2 * sqrt(1 - x^2) / (x^4 + 1) )
integrate(\(x) x^2 * sqrt(1 - x^2) / (x^4 + 1), 0, 1)
integrate(\(x) 1/2 * sqrt(x * (1 - x)) / (x^2 + 1), 0, 1)
integrate(\(x) 1/2 * sqrt(tan(x) * (1 - tan(x))), 0, pi/4); # initial I;
pi/2 * (-1 + 2^(1/4) * sin(3*pi/8))

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


##################
##################

### Radical Pow: 3

### I( x^3 * (x^3 + 1)^(2/3) / (x^6 + 1) )
integrate(\(x) x^3 * (x^3 + 1)^(2/3) / (x^6 + 1), 0, 1)
integrate(\(x) x^4 / (x^3 - 1) / (x^6 - 2*x^3 + 2), 2^(1/3), Inf)
# Fraction decomposition: see below;


### I( x^3 * (1 - x^3)^(2/3) / (x^6 + 1) )
integrate(\(x) x^3 * (1 - x^3)^(2/3) / (x^6 + 1), 0, 1)
integrate(\(x) 1/3 * (x*(1-x)^2)^(1/3) / (x^2 + 1), 0, 1)
# y = x^3 / (1+x^3) => x^3 = 1/(1-y) - 1;
integrate(\(x) x^3 / (1 + x^3)^3 / ((x^3/(1+x^3))^2 + 1), 0, Inf)
integrate(\(x) x^3 / (1 + x^3) / (2*x^6 + 2*x^3 + 1), 0, Inf)
integrate(\(x) x^4 / (x^3 + 1) / (x^6 + 2*x^3 + 2), 0, Inf)
# Fraction decomposition:
bp = (1+1i)^(1/3); bn = (1-1i)^(1/3);
integrate(\(x) 1/3 * (1/(x+1) - (x+1)/(x^2 - x + 1)) +
	- 1/6 * Re(bp^2 / (x+bp) + bn^2 / (x+bn) +
	- bp^2 * (x+bp) / (x^2 - bp*x + bp^2) +
	- bn^2 * (x+bn) / (x^2 - bn*x + bn^2) ), 0, Inf)
beta(2/3,1/3) / 3 * (-1 + 1/2*(bp^2 + bn^2))
beta(2/3,1/3) / 3 * (-1 + 2^(1/3) * cos(pi/6))


# Derivation:

###
integrate(\(x) x^3 * (x^3 + 1)^(2/3) / (x^6 + 1), 0, 1)
integrate(\(x) 1/3 * (x*(x+1)^2)^(1/3) / (x^2 + 1), 0, 1)
# rev of: y = x^3 / (1-x^3) => x^3 = 1 - 1/(y+1);
integrate(\(x) x^3 / (1 - x^3)^3 / ((x^3/(1-x^3))^2 + 1), 0, 2^(-1/3))
integrate(\(x) x^3 / (1 - x^3) / (2*x^6 - 2*x^3 + 1), 0, 2^(-1/3))
integrate(\(x) x^4 / (x^3 - 1) / (x^6 - 2*x^3 + 2), 2^(1/3), Inf)
# Fraction decomposition;
integrate(\(x) 1/2 * x^4 / (x^3 - 1) *
	Im(1/(x^3 - (1+1i)) - 1/(x^3 - (1-1i))), 2^(1/3), Inf)
integrate(\(x) 1/2 * x^4 * Re(2/(x^3 - 1) +
	- 1/(x^3 - (1+1i)) - 1/(x^3 - (1-1i))), 2^(1/3), Inf)
integrate(\(x) x^4/(x^3 - 1) - 1/2 * x^4 *
	Re(1/(x^3 - (1+1i)) + 1/(x^3 - (1-1i))), 2^(1/3), Inf)
integrate(\(x) x/(x^3 - 1) - 1/2 * x *
	Re((1+1i)/(x^3 - (1+1i)) + (1-1i)/(x^3 - (1-1i))), 2^(1/3), Inf)
bp = (1+1i)^(1/3); bn = (1-1i)^(1/3);
integrate(\(x) 1/3 * (1/(x-1) - (x-1)/(x^2 + x + 1)) +
	- 1/6 * Re((1+1i)/bp / (x-bp) + (1-1i)/bn / (x-bn) +
	- (1+1i)/bp * (x-bp) / (x^2 + bp*x + bp^2) +
	- (1-1i)/bn * (x-bn) / (x^2 + bn*x + bn^2) ), 2^(1/3), Inf)
# TODO


x = sqrt(7)
x / (x^3 - 1)
x * (1/(x-1) - (x+2)/(x^2 + x + 1)) / 3
(1/(x-1) - (x-1)/(x^2 + x + 1)) / 3


x = sqrt(7); b = (sqrt(2) + 1i)^(1/3);
x / (x^3 - b^3)
x * (1/(x-b) - (x + 2*b)/(x^2 + b*x + b^2)) / (3*b^2)
(1/(x-b) - (x-b)/(x^2 + b*x + b^2)) / (3*b)
# Note: x = b*y is an alternative;


###############

### I( x^5 * (1 - x^5)^(4/5) / (x^10 + 1) )
integrate(\(x) x^5 * (1 - x^5)^(4/5) / (x^10 + 1), 0, 1)
integrate(\(x) x^8 / (x^5 + 1) / (x^10 + 2*x^5 + 2), 0, Inf)
beta(1/5, 4/5) / 5 * (-1 + 2^(2/5) * cos(2*pi/10))

