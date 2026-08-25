#####################
##
## Leonard Mada
## (the one and only)
##
## Double Integrals:
## Radicals of Polynomials
## Diff-Type
##
## v.0.1a


### Double Integrals:
# Radicals of Polynomials:
# Diff-Type: ( abs(x - y) * ... )^(1/n)


### Examples:
# I( sqrt( abs(x - y) / (1 - x^2*y^2) ) )


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;


####################

############
### SQRT ###

### Pow = 1

### I( sqrt( abs(x-y) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( abs(x-y) / (1 - x*y) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- (pi + 2*log(2) - 5) * 4/3;


### I( sqrt( x * abs(x-y) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x * abs(x-y) / (1 - x*y) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(pracma::psi(1, 1/4) - 3*pracma::psi(1, 3/4)) / 16 - 1/6;
2*Catalan - pi^2/8 - 1/6;


### I( sqrt( x*y * abs(x-y) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x*y * abs(x-y) / (1 - x*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- gamma(3/4)^3 / gamma(9/4) * sqrt(2) * 3/2 + 3 + 11/15;


### I( sqrt( x/y * abs(x-y) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x/y * abs(x-y) / (1 - x*y) ), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
8/9;

### I( sqrt(x * abs(x-y)) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt(x * abs(x-y)) / (1 - x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pracma::psi(1, 3/4) / 2 - pi * 3/2 + 4;

### I( sqrt(x*y * abs(x-y)) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt(x*y * abs(x-y)) / (1 - x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
gamma(3/4)^2 * gamma(1/2) * sqrt(2) * 4 - pi * 14/3;

### I( x * sqrt(abs(x-y)) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * sqrt(abs(x-y)) / (1 - x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
4/9;

### I( sqrt(abs(x-y) / x) / (1 - x*y) )
# Note: symmetric, but better numerical stability;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt(abs(x-y) / y) / (1 - x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2/4 + pi * (1-log(2)) - 2;


### 2 Components:

### I( sqrt( abs(x-y) * (1-x) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( abs(x-y) * (1-x) / (1 - x*y) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=3E-12)
# TODO


### I( sqrt( abs(x-y) / ((1-x)*(1-y)*(1 - x*y)) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( abs(x-y) / ((1-x)*(1-y)*(1 - x*y)) ), 0, 1, rel.tol=1E-9)$value), 0, 1, rel.tol=3E-10)
4*pi - gamma(3/4)^2 / gamma(1/2) * sqrt(2) * 8;


### I( sqrt( abs(x-y) * (1-x)*(1-y) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( abs(x-y) * (1-x)*(1-y) / (1 - x*y) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=3E-12)
(5*pi - gamma(3/4)^2 / gamma(1/2) * sqrt(2) * 64/5) * 2/3;


### I( sqrt( abs(x-y) * (1 - x*y) * (1-x)/(1-y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( abs(x-y) * (1 - x*y) * (1-x)/(1-y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi*13/2 - gamma(3/4)^2 / gamma(1/2) * sqrt(2) * 64/5) * 2/15;


### I( sqrt( abs(x-y) / (1 - x*y) * (1-x)/(1-y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( abs(x-y) / (1 - x*y) * (1-x)/(1-y) ), 0, 1, rel.tol=5E-12)$value), 0, 1, rel.tol=1E-12)
pi * 7/3 - gamma(3/4)^2 / gamma(1/2) * sqrt(2) * 16/3;


###########
### Pow = 2

### I( sqrt( abs(x^2-y^2) * (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( abs(x^2-y^2) * (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
16/15 - gamma(3/4)^3 / gamma(9/4) * sqrt(2) / 4;
16/15 - gamma(3/4)^4 / gamma(1/2)^2 * 4/5;


### I( sqrt( abs(x^2-y^2) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( abs(x^2-y^2) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
2 - 2*gamma(3/4)^4 / gamma(1/2)^2;

### I( sqrt( x*y * abs(x^2-y^2) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x*y * abs(x^2-y^2) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
8/3 - gamma(3/4)^2 / gamma(1/2) * sqrt(2) * 2;


### I( sqrt( x/y * abs(x^2-y^2) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x/y * abs(x^2-y^2) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
8/3 - beta(3/4, 1/2) * beta(1/2, 1/2) / 4;


### I( x * sqrt( abs(x^2-y^2) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * sqrt( abs(x^2-y^2) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
log(2)/2 - pi/4 + 3/4;


### I( x*y^2 * sqrt( abs(x^2-y^2) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y^2 * sqrt( abs(x^2-y^2) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan/2 - pi^2/32 - 1/24;


### I( x^2*y^2 * sqrt( abs(x^2-y^2) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2*y^2 * sqrt( abs(x^2-y^2) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
14/15 - gamma(3/4)^4 / gamma(1/2)^2 * 6/5;


### Type: abs(x-y)

### I( sqrt( abs(x-y)/(x+y) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( abs(x-y)/(x+y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
2*Catalan - pi^2/8;


### I( sqrt( abs(x-y) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( abs(x-y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO

### I( x * sqrt( abs(x-y) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * sqrt( abs(x-y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(1/3 - 1/5 - 1/2^(3/2)/3 + 1/2^(5/2)/5) * sqrt(2) * 40/3 - 2/3;


####################
####################

### Power = 3

### I( (abs(x^3 - y^3) / (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(abs(x^3 - y^3) / (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
# TODO

### I( x * (abs(x^3 - y^3) / (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * (abs(x^3 - y^3) / (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
3/4 - (digamma(2/3) - digamma(1/3)) * beta(2/3, 2/3) / 9;

### I( x^2 * (abs(x^3 - y^3) / (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * (abs(x^3 - y^3) / (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
(1 - 2/3*log(2)) * 4/9;

### I( (abs(x - y) / (1 - x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(abs(x - y) / (1 - x*y))^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


####################

### Power = 4

### I( (abs(x - y^4) / (1 - x*y^4))^(1/4) )
# Note: same on both intervals;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(abs(x - y^4) / (1 - x*y^4))^(1/4), 0, y^4, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(abs(x - y^4) / (1 - x*y^4))^(1/4), y^4, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)$value;
(digamma(1/8) - digamma(5/8)) / 8 + pi/8 + 5/4;

### I( (abs(x - y^2) / (1 - x*y^2))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(abs(x - y^2) / (1 - x*y^2))^(1/4), 0, y^2, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(abs(x - y^2) / (1 - x*y^2))^(1/4), y^2, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)$value;
# TODO


### I( y^(-1/4) * (abs(x - y) / (1 - x*y))^(1/4) )
# on [0, y]
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y^(-1/4) * (abs(x - y) / (1 - x*y))^(1/4), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi/3 + (log(2) - 1) * 2;
# on [y, 1]
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y^(-1/4) * (abs(x - y) / (1 - x*y))^(1/4), y, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# 4/3 - I( x * sqrt( (1+x) / (1 - x^2*y^2) ) );
# TODO


### I( y^(3/4) * (abs(x - y) / (1 - x*y))^(1/4) )
# on [0, y]
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y^(3/4) * (abs(x - y) / (1 - x*y))^(1/4), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(3 - pi/2 - log(2)) * 2/5;


####################

### Power = 5

### I( (abs(x - y^5) / (1 - x*y^5))^(1/5) )
n = 5;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(abs(x - y^n) / (1 - x*y^n))^(1/n), 0, y^n, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(abs(x - y^n) / (1 - x*y^n))^(1/n), y^n, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)$value;
# TODO

