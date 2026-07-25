#####################
##
## Leonard Mada
## (the one and only)
##
## Double Integrals:
## Radicals of Polynomials
## Inhomogeneous Radicals
##
## v.0.1a


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;


####################

### SQRT:

### Type: w. (1-x)

### I( sqrt( (1-x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
#
integrate(\(x) 2 * sapply(x, \(y) integrate(\(x)
	sqrt( (1-x) * (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)$value +
	- (1/3 - 1/5 - 1/2^(3/2)/3 + 1/2^(5/2)/5) * sqrt(2) * 8;
# TODO


### I( sqrt( (1 - x^2*y^2) / (1-x) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1 - x^2*y^2) / (1-x) ), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-13)
# I( sqrt( (1-x) * (1 - x^2*y^2) ) ) +
	+ (1/3 - 1/5 - 1/2^(3/2)/3 + 1/2^(5/2)/5) * sqrt(2) * 16;
# TODO


### I( sqrt( x*(1-x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x*(1-x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( sqrt( y*(1-x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( y*(1-x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( sqrt( (1-x)*(1-y) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x)*(1-y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( sqrt( (1-x) * (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x) * (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/2 * sqrt( (1-x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value
	), 0, 1, rel.tol=1E-13)$value +
	+ (1/3 - 1/5 - 1/2^(3/2)/3 + 1/2^(5/2)/5) * sqrt(2) * 4;
# TODO


### I( sqrt( (1-x)*(1-y) * (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x)*(1-y) * (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# I(sqrt( (1-x)*(1-y) / (1 - x^2*y^2) ) ) * 2/5 +
# I(sqrt( (1-x)/(1-y) * (1 - x^2*y^2) ) ) / 5;
# TODO


### I( 1 / sqrt( (1+x) * (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / sqrt( (1+x) * (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) 8 * x^3 * atan(x) / (x^4 + 1), 0, 1, rel.tol=1E-13);
(digamma(7/8) - digamma(3/8))^2 / 16 +
	- (digamma(5/8) - digamma(1/8))^2 / 16 + 4*Catalan;


### I( sqrt( (1+x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1+x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	2 * sqrt( (1+x) * (1 - x^2*y^2) ), 0, 1)$value), 0, 1)$value - 14/15;
# TODO


### I( y * sqrt( (1+x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * sqrt( (1+x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
log(tan(pi/8)) - sqrt(2) + 3;


### I( x*y * sqrt( (1+x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * sqrt( (1+x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
1.337828528819 - 0.9388152447781;
# see: I( sqrt( (1+x) / (1 - x^2*y^2) ) ) & I( ((1+x) * (1 - x*y) / (1 + x*y))^(1/2) );
# TODO


### Mixed:

### I( sqrt( (1-x)/(1+y) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x)/(1+y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO

### I( sqrt( (1+x)/(1+y) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1+x)/(1+y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### I( sqrt( (1-y) / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-y) / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(sqrt(2) - 1)*6 + 2*log(sqrt(2) - 1);


### I( sqrt( (1-x)*(1-y) / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x)*(1-y) / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


