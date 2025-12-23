#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
## Type: PROD( Atan )
##
## v.0.1b


### Type: PROD( ATAN )

### Examples:
# I( atan(x/y)^2 )
# I( atan(x*y) * atan(x/y) )
# I( atan(x/y) * log(1-x) )
# I( atan(x/y) * log(1+x) )


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;
dzeta2  = - 0.937548254316;

#####################
#####################

### Powers & Product

### Prod w. Log

### Type: ATAN(x/y)

### I( atan(x/y) * log(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * log(x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2 / 96 - pi/4 + log(2)/4;

### I( atan(x/y) * log(y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * log(y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- pi^2 / 96 - pi/4 - log(2)/4;

### I( x * atan(x/y) * log(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x/y) * log(x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (7/4*pi + log(2) + 3*Catalan - 11/2) / 18;

### I( x * atan(x/y) * log(y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x/y) * log(y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-13)
(- pi^2 / 8 + 2*pi - 4*log(2) + 1) / 18 - pi/4;

### I( y * atan(x/y) * log(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(x/y) * log(x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (- pi^2 / 8 + 2*pi - 4*log(2) + 1) / 18;

### I( y * atan(x/y) * log(y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(x/y) * log(y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-13)
(7/4*pi + log(2) + 3*Catalan - 11/2) / 18 - pi/8;
- (pi/2 - log(2) - 3*Catalan + 11/2) / 18;


### I( atan(x/y) * log(1-x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * log(1-x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi^2 / 12 + pi/8 * log(2) - 3/8 * pi - Catalan;

### I( atan(x/y) * log(1-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * log(1-y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (pi^2 / 12 + pi/8 * log(2) + pi/8 - Catalan);


### I( atan(x/y) * log(1+x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * log(1+x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (pi^2 / 24 - 3/8 * pi*log(2) + pi/8 - log(2)/2);

### I( atan(x/y) * log(1+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * log(1+y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2 / 24 + 5/8 * pi*log(2) - 3/8 * pi - log(2)/2;


### I( atan(x/y) * log(x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * log(x+y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi/4 * (2*log(2) - 3/2);

### I( atan(x/y) * log(y-x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * log(y-x), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(5/96*pi^2 + pi*log(2)/8 - log(2)^2 / 8 - 3/8 * pi + 3/4 * log(2) - Catalan) / 2;
#
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * log(x-y), y, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (5/96*pi^2 + pi*log(2)/8 - log(2)^2 / 8 + 3/8 * pi + 3/4 * log(2) - Catalan) / 2;


### I( atan(x/y) * log(x) * log(y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * log(x) * log(y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi/4;

### I( x * atan(x/y) * log(x) * log(y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x/y) * log(x) * log(y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(3*Catalan + pi^2 / 8 + 17/4 * pi + 5*log(2) - 13/2) / 54;

### I( y * atan(x/y) * log(x) * log(y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(x/y) * log(x) * log(y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (3*Catalan + pi^2 / 8 - 5/2 * pi + 5*log(2) - 13/2) / 54;

### I( x*y * atan(x/y) * log(x) * log(y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * atan(x/y) * log(x) * log(y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi/64;


### ATAN(x*y)

### I( atan(x*y) * log(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * log(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi^2 / 48 - pi/4 + 1/2*log(2) + 3/32 * pracma::zeta(3);


### I( atan(x*y) * log(1-x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * log(1-x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
-23/64 * pracma::zeta(3) + pi * Catalan / 4 - Catalan +
	+ 5/96*pi^2 + pi*log(2)/8 - pi/4 - log(2)^2 / 8 + log(2)/2;

### I( atan(x*y) * log(1+x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * log(1+x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
33/64 * pracma::zeta(3) - pi*Catalan / 4 +
	+ (- pi^2 / 48 + 3/4*pi*log(2) - pi/2 - log(2)^2 / 4 + log(2)) / 2;


### I( atan(x*y) * log(x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * log(x+y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2/48 - (- 9/16 * pracma::zeta(3) + pi * Catalan / 4) +
	- (pi^2/24 + pi*log(2)/2 - pi/2 - 4*log(2) + 3) / 4 +
	- (1/96*pi^2 - pi*log(2)/2 + log(2)^2 / 8 + log(2)/4 + pi/2 - 3/4);
9/16 * pracma::zeta(3) - pi * Catalan / 4 +
	+ (3*pi*log(2) - 3*pi - log(2)^2 + 6*log(2)) / 8;

### I( atan(x*y) * log(y-x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * log(y-x), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (5/16 * pracma::zeta(3) - Catalan * pi / 4 + Catalan) / 2 +
	+ (pi^2 / 16 + pi*log(2)/8 - log(2)^2 / 8 + 3/4*log(2) - 3/8 * pi) / 2;


##################

### w. Simple-Atan

### I( atan(x/y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * atan(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi^2 / 16 + 3/8 * pi*log(2) - pi/8 - log(2)/4 - Catalan/2;

### I( atan(x/y) * atan(y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * atan(y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi/4 * (pi/2 - log(2)) - (pi^2/16 + 3/8 * pi*log(2) - pi/8 - log(2)/4 - Catalan/2);
pi^2 / 16 - 5/8 * pi*log(2) + pi/8 + log(2)/4 + Catalan/2;


### I( atan(x*y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * atan(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) - log(1-x) * log(x) / (x^2+1), 0, 1, rel.tol=1E-13)$value +
	+ pi*log(2)/4 + (pi/4)^3 + (pi/4)^2 - Catalan * log(2)/2 - Catalan;
# TODO

### I( x * atan(x*y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x*y) * atan(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(pi/4 - log(2)/2)^2 / 2;

### I( y * atan(x*y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(x*y) * atan(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi * log(2) / 4 - Catalan / 2;

### I( 1/x * atan(x*y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/x * atan(x*y) * atan(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- 7/8 * pracma::zeta(3) + pi*Catalan / 2 - (5/12 * pi^2 - pi*log(2) - log(2)^2) / 8;

### I( y/x * atan(x*y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y/x * atan(x*y) * atan(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(- 7/4 * pracma::zeta(3) + pi*Catalan) / 4 - pi^2/32 + 1/4;

### I( 1/x^2 * atan(x*y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/x^2 * atan(x*y) * atan(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- pi^2 / 16 + 3/8 * pi * log(2) - pi/8 + log(2)/4 + Catalan / 2;

### I( y/x^2 * atan(x*y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y/x^2 * atan(x*y) * atan(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(2*Catalan - pi^2 / 4 + pi*log(2)/2 + pi/4 + 1/2) / 6;

### I( y^2/x^2 * atan(x*y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y^2/x^2 * atan(x*y) * atan(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(- pi^2 / 16 + pi/4 * log(2) + Catalan) / 3 - (Catalan/2 + (pi/8 - log(2)/2 - 1) / 3 + 1/4) / 6;
(- pi^2 / 4 + pi*log(2) + 3*Catalan) / 12 - (pi/4 - log(2) - 1/2) / 36;


### I( 1/y * atan(x*y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/y * atan(x*y) * atan(x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-12)
# TODO

### I( x/y * atan(x*y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x/y * atan(x*y) * atan(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
7/16 * pracma::zeta(3) - (pi^2 / 16 - pi/2 + log(2) + Catalan) / 2;

### I( 1/(x*y) * atan(x*y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/(x*y) * atan(x*y) * atan(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
Catalan^2 / 2;


### w. Composite-Atan

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
integrate(\(x) 2*x * atan(x) * log(x) / (x^2+1), 0, 1)$value +
	+ (pi/4)^2 + pi*log(2)/4 - Catalan;
# TODO


### I( atan(1-x*y) * atan(x/y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(1-x*y) * atan(x/y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- pi/4 * (5/96*pi^2 + pi*log(2)/8 - log(2)^2 / 8 + log(2)/2 - pi/4 - Catalan);


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

