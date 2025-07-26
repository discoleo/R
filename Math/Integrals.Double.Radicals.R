#####################
##
## Leonard Mada
## (the one and only)
##
## Double Integrals:
## Radicals of Polynomials
##
## v.0.1a


### Double Integrals:
# Radicals of Polynomials

### Examples:
# I( 1 / sqrt( (1-x^2) * (1-y^2) * (1 - x^2*y^2) ) )


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;


########################

### SQRT:

### I( sqrt( (1-x)*(1-y) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x)*(1-y) / (1 - x*y) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
6*Catalan - 5;


### I( 1 / sqrt( (1-x^2) * (1-y^2) * (1 - x^2*y^2) ) )
# - Series started as Trig-Radicals in:
#   Integrals.Double.Trig.R;

### Base:
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / sqrt((1-x^2) * (1-y^2) * (1 - x^2*y^2)), 0, 1)$value), 0, 1)
gamma(1/4)^3 / gamma(3/4) * sqrt(2) / 16;


### I( sqrt( (1 - x^2*y^2) * (1-x^2) / (1-y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1 - x^2*y^2) * (1-x^2) / (1-y^2) ), 0, 1)$value), 0, 1)
gamma(1/4)^3 / gamma(3/4) * sqrt(2) / 48;


### I( sqrt( (1-x^2) / ((1-y^2) * (1 - x^2*y^2)) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^2) / ((1-y^2) * (1 - x^2*y^2)) ), 0, 1)$value), 0, 1)
(gamma(1/4)^3 / gamma(3/4) / 16 - gamma(3/4)^3 / gamma(1/4)) * sqrt(2) / 2;


### I( 1 / sqrt((1-x^2) * (1 - x^2*y^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / sqrt((1-x^2) * (1 - x^2*y^2)), 0, 1)$value), 0, 1)
2*Catalan;


### I( 1 / sqrt(1 - x^2*y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / sqrt(1 - x^2*y^2), 0, 1)$value), 0, 1)
pi*log(2)/2;


### I( sqrt( x / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sqrt( x / (1 - x^2*y^2) ), 0, 1)$value), 0, 1)
pi - gamma(3/4)^2 / gamma(1/2) * 2 * sqrt(2);


### I( sqrt( x/y / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x/y / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
gamma(1/4)^2 / gamma(1/2) * sqrt(2) / 4 - gamma(3/4)^2 / gamma(1/2) * sqrt(2);
(beta(1/4, 1/4) - 2*beta(3/4, 3/4)) * sqrt(2) / 4;


### I( sqrt( (1-x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( sqrt( (1-x)*(1-y) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x)*(1-y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( x / sqrt((1-x^2) * (1-y^2) * (1 - x^2*y^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / sqrt((1-x^2) * (1-y^2) * (1 - x^2*y^2)), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/2 / sqrt((1-x) * (1-y^2) * (1 - x*y^2)), 0, 1)$value), 0, 1);
pi^2 / 4;


### I( x / sqrt( (1-y^2) * (1 - x^2*y^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / sqrt((1-y^2) * (1 - x^2*y^2)), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/2 / sqrt((1-y^2) * (1 - x*y^2)), 0, 1)$value), 0, 1);
1;


### I( x / sqrt((1-x^2) * (1 - x^2*y^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / sqrt((1-x^2) * (1 - x^2*y^2)), 0, 1)$value), 0, 1)
pi^2/8;

