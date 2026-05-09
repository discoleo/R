#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
## Type: Atan & Diff
##
## v.0.2e


### Type: Diff(Atan) & Atan(Diff(...))


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;
dzeta2  = - 0.937548254316;

#####################


### Div: (x-y)

# Note:
# - many integrals are computed over [0,1] x [0,y]
#   for improved numerical stability;

### I( (atan(x) - atan(y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x) - atan(y)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- pi*log(2)/8 + Catalan + Re(log(1-1i)^2)/2;

### I( x * (atan(x) - atan(y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * (atan(x) - atan(y)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
Catalan / 2 - pi*log(2)/8 - pi/8 - log(2)/4 + 1/2;

### I( x*y * (atan(x) - atan(y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * (atan(x) - atan(y)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(Catalan + pi^2/32 - pi*log(2)/8 - pi/4 - log(2)^2/8 - log(2)/2 + 1/2) / 3;

### I( x^2 * (atan(x) - atan(y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * (atan(x) - atan(y)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(Catalan + pi^2/32 - pi*log(2)/8 + pi/8 - log(2)^2/8 - 5/4 * log(2) - 1/4) / 3;

### I( (x*atan(x) - y*atan(y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(x*atan(x) - y*atan(y)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
Catalan / 2 - pi * log(2) / 8 + pi/8 - log(2)/4;
# Note: same on [y, 1];

### I( (y*atan(x) - x*atan(y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	2 * (y*atan(x) - x*atan(y)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- (pi*log(2)/4 + pi/4 - log(2)/2 - Catalan);

### I( (x^2*atan(x) - y^2*atan(y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(x^2*atan(x) - y^2*atan(y)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(Catalan + pi^2/32 - pi*log(2)/8 + pi/2 - log(2)^2/8 - log(2)/2 - 1) / 3;

### I( x * (x*atan(x) - y*atan(y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * (x*atan(x) - y*atan(y)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(2*Catalan + pi^2/16 - pi*log(2)/4 + pi/2 - log(2)^2/4 - 2*log(2) - 1) / 6;

### I( x * (y*atan(x) - x*atan(y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * (y*atan(x) - x*atan(y)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(4*Catalan + pi^2/8 - pi*log(2)/2 - 3/2 * pi - log(2)^2/2 - 3*log(2) + 3) / 12;


### I( (atan(x)/x - atan(y)/y) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	2 * (atan(x)/x - atan(y)/y) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) 2*x * atan(x) * log(x) / (x^2+1), 0, 1)$value +
integrate(\(x) log(x) * log(x^2+1) / (x^2+1), 0, 1)$value;
# TODO

### I( (atan(x)*y/x - atan(y)*x/y) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	-2 * (atan(x)*y/x - atan(y)*x/y) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
3/48 * pi^2 + pi*log(2)/4 - log(2)^2 / 4;

### I( (atan(x)*y^2/x - atan(y)*x^2/y) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	-2 * (atan(x)*y^2/x - atan(y)*x^2/y) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi*log(2)/4 + pi/4 - log(2)/2;


### I( (atan(x^2) - atan(y^2)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x^2) - atan(y^2)) / (x-y), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# Note: on [0, 1]^2 => 2 * I;
(pracma::psi(1, 7/8) - pracma::psi(1, 3/8)) / 32 +
- (pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 32 +
- (pracma::psi(1, 3/4) - pracma::psi(1, 1/4)) / 32 + # Catalan / 2
	- (digamma(3/4) - digamma(1/4)) * log(2) / 16 +  # pi*log(2) / 16
	- (digamma(7/8) - digamma(3/8))^2 / 64 +
	+ (digamma(5/8) - digamma(1/8))^2 / 64 +
	+ (digamma(7/8) - digamma(3/8)) * log(2) / 16 +
	- (digamma(5/8) - digamma(1/8)) *
		(digamma(3/4) - digamma(1/4)) / 32 +
	+ (digamma(1/4) - digamma(3/4)) / sin(pi/4) * pi/8;


### I( (atan(x/(1-x)) - atan(y/(1-y))) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x/(1-x)) - atan(y/(1-y))) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
Catalan - pi^2 / 16 + pi/4 * log(2);

### I( x * (atan(x/(1-x)) - atan(y/(1-y))) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * (atan(x/(1-x)) - atan(y/(1-y))) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(Catalan - pi^2 / 16 + pi*log(2)/4 - pi/4 + 1/2) / 2;

### I( x*y * (atan(x/(1-x)) - atan(y/(1-y))) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * (atan(x/(1-x)) - atan(y/(1-y))) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(2*Catalan - pi^2 / 16 + pi*log(2)/2 - pi/2 + 1/2) / 6;

### I( (x * atan(x/(1-x)) - y * atan(y/(1-y))) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(x * atan(x/(1-x)) - y * atan(y/(1-y))) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(Catalan - pi^2 / 16 + pi*log(2)/4 + pi/4) / 2;

### I( (y * atan(x/(1-x)) - x * atan(y/(1-y))) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(y * atan(x/(1-x)) - x * atan(y/(1-y))) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(Catalan - pi^2 / 16 + pi*log(2)/4 - pi/4) / 2;

### I( (atan(x/(1-x)) / x - atan(y/(1-y)) / y) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x/(1-x)) / x - atan(y/(1-y)) / y) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi^3 / 96;

### I( (atan(x/(1-x)) * y/x - atan(y/(1-y)) * x/y) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x/(1-x)) * y/x - atan(y/(1-y)) * x/y) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- pi^2 / 16;


#################

### Series: ATAN( x/(1-y) ) / (x - y)

### I( (atan(x/(1-y)) - atan(y/(1-x))) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x/(1-y)) - atan(y/(1-x))) / (x-y), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
0;

### I( atan(y/(1-x)) / (x - y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x/(1-x)) - atan(y/(1-x))) / (x-y), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=3E-13)
(Catalan - pi^2 / 16 + pi*log(2)/4) / 2;

### I( (x*atan(x/(1-y)) - y*atan(y/(1-x))) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(x*atan(x/(1-y)) - y*atan(y/(1-x))) / (x-y), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- Catalan / 2 + pi/4;

### I( (y*atan(x/(1-y)) - x*atan(y/(1-x))) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(y*atan(x/(1-y)) - x*atan(y/(1-x))) / (x-y), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- Catalan / 2;

### I( x * (y*atan(x/(1-y)) - x*atan(y/(1-x))) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * (y*atan(x/(1-y)) - x*atan(y/(1-x))) / (x-y), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (3*Catalan - pi/2 + log(2) - 1/2) / 6;

### I( (atan(x/(1-y)) / x - atan(y/(1-x)) / y) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x/(1-y)) / x - atan(y/(1-x)) / y) / (x-y), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (Catalan * log(2) + pi^3 / 24);

### I( (atan(x/(1-y)) * y/x - atan(y/(1-x)) * x/y) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x/(1-y)) * y/x - atan(y/(1-x)) * x/y) / (x-y), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (Catalan + pi/4 + log(2)/2);

### I( x * (atan(x/(1-y)) * y/x - atan(y/(1-x)) * x/y) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x/(1-y)) * y - atan(y/(1-x)) * x^2/y) / (x-y), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (2*Catalan - pi*log(2)/2 + pi/2 + 2*log(2) - 1) / 4;


### Other

### I( (atan(x/y) - atan(y/x)) / (x-y) )
# I( (2*atan(x/y) - pi/2)  / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x/y) - atan(y/x)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- pi*log(2)/4 + 2*Catalan;


### Series: ATAN( (x-y)/x ) / (x - y)

### I( atan((x-y)/x) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan((x-y)/x) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
Catalan + pi*log(2)/4;

### I( x * atan((x-y)/x) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan((x-y)/x) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
Catalan / 2 + pi*log(2)/8 - pi/8;

### I( y * atan((x-y)/x) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan((x-y)/x) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
Catalan / 2 + pi*log(2)/8;

### I( x*y * atan((x-y)/x) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * atan((x-y)/x) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(Catalan + pi*log(2)/4 - pi/4) / 3;


### ATAN( OTHER )

### I( atan((x-y)/(x+y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan((x-y)/(x+y)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
Catalan - pi * log(2) / 8;

### I( atan((x-y)/(x+y)) * x / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan((x-y)/(x+y)) * x / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(2*Catalan - pi * log(2) / 4 - log(2)) / 4;

### I( atan((x-y)/(x+y)) * y / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan((x-y)/(x+y)) * y / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- pi*log(2)/16 + Catalan/2;


### I( atan((x-y)/(x+y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan((x-y)/(x+y)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- pi*log(2)/8 + Catalan;

### I( x * atan((x-y)/(x+y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan((x-y)/(x+y)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(2*Catalan - pi * log(2) / 4 - log(2)) / 4;

### I( x*y * atan((x-y)/(x+y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * atan((x-y)/(x+y)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(2*Catalan - pi * log(2) / 4 - log(2)) / 6;


### I( atan((x-y)/(1-x*y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan((x-y)/(1-x*y)) / (x-y), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(7/12 * pi^2 - 3*pi*log(2) - log(2)^2) / 8 + Catalan;


### Div: (x^2 - y^2)

### I( (atan(x) - atan(y)) / (x^2-y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x) - atan(y)) / (x^2-y^2), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) log(1-x) * log(x) / (x^2+1), 0, 1)$value + Catalan * log(2) / 2;
# TODO

### I( x * (atan(x) - atan(y)) / (x^2-y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	2 * x * (atan(x) - atan(y)) / (x^2-y^2), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
Catalan - pi^2 / 16 - pi*log(2)/4 + log(2)^2;

### I( (x*atan(x) - y*atan(y)) / (x^2-y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	2 * (x*atan(x) - y*atan(y)) / (x^2-y^2), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
Catalan - pi^2 / 16 + pi*log(2)/4;

### I( (y*atan(x) - x*atan(y)) / (x^2-y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(y*atan(x) - x*atan(y)) / (x^2-y^2), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- (pi*log(2) - log(2)^2 / 2 - 2*Catalan) / 4;


#########################

### ATAN( SQRT(y-x) )

# TODO: some are still in Integrals.Double.Atan.R;


### ATAN( SQRT(y^2-x^2) )

### I( atan(sqrt(y^2-x^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sqrt(y^2-x^2)), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(asinh(1) - 2 + sqrt(2)) * pi/4;

### I( x * atan(sqrt(y^2-x^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(sqrt(y^2-x^2)), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(pi/4 - log(2)/4 - 1/2) * 2/3;

### I( y * atan(sqrt(y^2-x^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(sqrt(y^2-x^2)), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(sqrt(2) - 1 - 1/4) * pi/3;


### I( x*y * atan(sqrt(y^2-x^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * atan(sqrt(y^2-x^2)), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi/8 - 1/3;


### ATAN( sqrt(y^2-x^2) / x )

### I( atan(sqrt(y^2-x^2) / x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sqrt(y^2-x^2) / x), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
1/2;

### I( x * atan(sqrt(y^2-x^2) / x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(sqrt(y^2-x^2) / x), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi/24;

### I( y * atan(sqrt(y^2-x^2) / x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(sqrt(y^2-x^2) / x), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
1/3;

### I( x*y * atan(sqrt(y^2-x^2) / x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * atan(sqrt(y^2-x^2) / x), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi/32;

### on [0,1] x [y,1]

### I( atan(sqrt(y^2-x^2) / x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sqrt(x^2-y^2) / x), y, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi * (sqrt(2) - 1) / 4;

### I( x * atan(sqrt(y^2-x^2) / x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(sqrt(x^2-y^2) / x), y, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(sqrt(2) - 1) * pi/6;

### I( y * atan(sqrt(y^2-x^2) / x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(sqrt(x^2-y^2) / x), y, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi/4 - 1/2) / 3;

### I( x*y * atan(sqrt(y^2-x^2) / x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * atan(sqrt(x^2-y^2) / x), y, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi/16 - 1/8;


### ATAN( sqrt(y^2-x^2) / (x*y) )

### I( atan(sqrt(y^2-x^2) / (x*y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sqrt(y^2-x^2) / (x*y)), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi/2 - 1;

### I( x * atan(sqrt(y^2-x^2) / (x*y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(sqrt(y^2-x^2) / (x*y)), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi*log(2)/4 - pi/8;

### I( y * atan(sqrt(y^2-x^2) / (x*y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(sqrt(y^2-x^2) / (x*y)), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2 / 16 - 1/4;

### I( x*y * atan(sqrt(y^2-x^2) / (x*y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * atan(sqrt(y^2-x^2) / (x*y)), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- pi*log(2)/4 + pi * 5/24;


### ATAN( sqrt(y^2-x^2) / x^2 )

### I( atan(sqrt(y^2-x^2) / x^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sqrt(y^2-x^2) / x^2), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(digamma(5/6) - digamma(1/6)) / 2 - 2/3 * pi;

### on [0,1] x [y,1]
# Note: numerically problematic even with Rmpfr;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sqrt(x^2-y^2) / x^2), y, 1, rel.tol=2E-11)$value), 0, 1, rel.tol=3E-11)
(sqrt(2) - 1) * pi/3;


### I( x * atan(sqrt(y^2-x^2) / x^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(sqrt(y^2-x^2) / x^2), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y/x * atan(sqrt(x^2-1) * x / y) * y/x^2, 1, Inf, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y^2 * cos(x)*sin(x) * atan(sin(x)/cos(x)^2 / y), 0, pi/2, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y^2 * atan(x/(1-x^2) / y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi/12 - 1/2 * integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2*y^3 * (1+x^2) / (x^2 + (1-x^2)^2*y^2), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)$value;
pi/12 - 1/4 + 1/4 * integrate(\(y)
	y^4 * (1+y^2) / (1-y^2)^4 * (log(y^2 + (1-y^2)^2) - log(y^2)) - (3*y^2-1)/(1-y^2)^2, 0, 1)$value;
# TODO

### on [0,1] x [y,1]
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(sqrt(x^2-y^2) / x^2), y, 1, rel.tol=3E-11)$value), 0, 1, rel.tol=1E-11)
(digamma(5/8) - digamma(3/8)) / 4 - (asinh(1) - 2 + sqrt(2)) * pi/16;


### I( y * atan(sqrt(y^2-x^2) / x^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(sqrt(y^2-x^2) / x^2), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
((digamma(5/6) - digamma(1/6)) - 13/12 * pi) / 5;

### on [0,1] x [y,1]
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(sqrt(x^2-y^2) / x^2), y, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(2*pi - log(2) - 2) / 30;

