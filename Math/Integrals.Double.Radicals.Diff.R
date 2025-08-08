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


### I( sqrt( x*y * abs(x-y) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x*y * abs(x-y) / (1 - x*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- gamma(3/4)^3 / gamma(9/4) * sqrt(2) * 3/2 + 3 + 11/15;


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


### I( x^2*y^2 * sqrt( abs(x^2-y^2) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2*y^2 * sqrt( abs(x^2-y^2) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
14/15 - gamma(3/4)^4 / gamma(1/2)^2 * 6/5;


### I( sqrt( abs(x-y)/(x+y) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( abs(x-y)/(x+y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
2*Catalan - pi^2/8;


### I( sqrt( abs(x-y) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( abs(x-y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO

