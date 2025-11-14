#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
## Type: PROD( Atan )
##
## v.0.1a


### Type: PROD( ATAN )

### Examples:
# I( atan(x/y)^2 )
# I( atan(x*y) * atan(x/y) )


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;
dzeta2  = - 0.937548254316;

#####################
#####################

### Powers & Product

### I( atan(x/y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y)^2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi*log(2)/2 + pi^2/16 - Catalan;

### I( atan(x*y) * atan(x/y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * atan(x/y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi/2 * (- pi^2 / 96 + pi/8 - log(2)/4);

### I( atan(x*y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y)^2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) - atan(x)^2 * log(x), 0, 1, rel.tol=1E-12)
# TODO


### Div:

### I( 1/x * atan(x/y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/x * atan(x/y)^2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi*Catalan / 2 - Catalan + pi^2 / 16 + 3/4 * pi*log(2) - 7/8 * pracma::zeta(3);

### I( y/x * atan(x/y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y/x * atan(x/y)^2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi*Catalan / 4 - pi^2 / 32 + pi/8 + log(2)/4 - 7/16 * pracma::zeta(3);


### I( 1/x * atan(x*y) * atan(x/y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * atan(x/y) / x, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(Catalan - pi/4 + log(2)/2) * pi/2 - 7/8 * pracma::zeta(3) + 2*Catalan - pi^2 / 8;
pi*Catalan/2 + 2*Catalan - pi^2/4 + pi*log(2)/4 - 7/8 * pracma::zeta(3);

### I( 1/y * atan(x*y) * atan(x/y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * atan(x/y) / y, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
7/8 * pracma::zeta(3) - 2*Catalan + pi^2 / 8;

### I( y/x * atan(x*y) * atan(x/y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y/x * atan(x*y) * atan(x/y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi*Catalan/4 - pi^2 / 24 + pi/8 - 7/16 * pracma::zeta(3);

### I( x/y * atan(x*y) * atan(x/y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x/y * atan(x*y) * atan(x/y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
7/16 * pracma::zeta(3) - pi^2 / 48;


### I( atan(x*y)^2 / x )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y)^2 / x, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
-7/8 * pracma::zeta(3) + pi*Catalan / 2 - (pi/4)^2 - pi*log(2)/4 + Catalan;

### I( y/x * atan(x*y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y/x * atan(x*y)^2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- 7/16 * pracma::zeta(3) + pi*Catalan / 4 - pi^2 / 32 + pi/8 - log(2)/4;



### Coeff: x or y

### I( x * atan(x/y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x/y)^2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(pi^2/8 + 3/4 * pi*log(2) - pi/4 + log(2)/2 - Catalan) / 3;

### I( y * atan(x/y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(x/y)^2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(pi*log(2)/4 + pi/4 + log(2)/2 - Catalan) / 3;


### I( x * atan(x/y) * atan(x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x/y) * atan(x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- pi^2 / 16 - 5/12 * pi*log(2) + 5/12 * pi + log(2)/2;

### I( y * atan(x/y) * atan(x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(x/y) * atan(x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2 / 16 + pi * log(2)/6 - pi/6 - log(2)/2;


### I( x * atan(x*y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x*y)^2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi*log(2)/4 + pi/4 - log(2)/2 - Catalan;

