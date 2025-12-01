#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
## Type: Atan
##
## v.0.2d


### Type: Atan

### Examples:
# I( x*y * atan(x) / (x^2+y^2) )
# I( x * atan(1 + x - x*y) )
# I( atan(x/y) / (x+y) )
# I( atan(x*y) / (1 + x*y) )
# I( atan(x*y) / (1 - x*y) ) [not yet]

### History
# - [refactor] PROD( ATAN ) moved to file:
#   Integrals.Double.Atan.Prod.R;
# - [refactor] ATAN( Trig ) moved to file:
#   Integrals.Double.Atan.Trig.R;


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;
dzeta2  = - 0.937548254316;

#####################
#####################

############
### ATAN ###

# Note: (x+y) does not behave as nicely with ATAN;

### I( atan((x+y)/2) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan((x+y)/2), 0, 1)$value), 0, 1)
3/2*pi + 2*log(5) - 6*log(2) - 3*atan(2)


### I( atan(x/y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y), 0, 1)$value), 0, 1)
pi/4;

### I( atan(x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y), 0, 1)$value), 0, 1)
- pi^2 / 48 + pi/4 - log(2)/2;

### I( x * atan(x/y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x/y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi + log(2) - 1) / 6;

### I( y * atan(x/y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(x/y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi/2 - log(2) + 1) / 6;

### I( atan(1 - x + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(1-x+x*y), 0, 1)$value), 0, 1)
-5/96*pi^2 - pi*log(2)/8 + pi/4 + log(2)^2 / 8 - log(2)/2 + Catalan;

### I( atan(1 + x - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(1+x-x*y), 0, 1)$value), 0, 1)
integrate(\(x) log(x) * x / (x^2 + 2*x + 2) , 0, 1, rel.tol=1E-13)$value +
	+ 2*atan(2) - 1/2*log(5) + 1/2*log(2) - atan(1);
# TODO

### I( x * atan(1 + x - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) x * atan(1+x-x*y), 0, 1)$value), 0, 1)
3/2*atan(2) - log(5) + log(2) - pi/4 + 1/2;

### I( y * atan(1 + x - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) y * atan(1+x-x*y), 0, 1)$value), 0, 1)
# I(atan(1 + x*y)) - I( x * atan(1 + x - x*y) );
# TODO


### I( atan(abs(x-y) / (x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(abs(x-y) / (x+y)), 0, y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
log(2)/4;

### I( atan(2*x / (x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(2*x / (x+y)), 0, 1)$value), 0, 1)
integrate(\(x) atan(2*x/(x+1)) +
	sapply(x, \(y) integrate(\(x) 2*x*y / (4*x^2 + (x+y)^2), 0, 1)$value), 0, 1);
37/100 * pi + 6/5 * log(2) - log(5)/2 + 1/50 * (2*atan(3) - 23*atan(2));


# Helper:
integrate(\(x) atan(2*x/(x+1)), 0, 1)
integrate(\(x) atan(2 - 2/(x+1)), 0, 1)
pi/4 - (log(8)/5 - 1/5 * (atan(3) - atan(1/2)));

#
integrate(\(x) sapply(x, \(y) integrate(\(x) x*y / (4*x^2 + (x+y)^2), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x) y*(x-y) / (4*y^2 + x^2), y, y+1)$value), 0, 1)
integrate(\(x) -1/4 + sapply(x, \(y)
	integrate(\(x) (x*y + x^2/4) / (4*y^2 + x^2), y, y+1)$value), 0, 1)
integrate(\(x) x/2 * (log(5*x^2 + 2*x + 1) - log(5*x^2)) + sapply(x, \(y)
	integrate(\(x) - y^2 / (4*y^2 + x^2), y, y+1)$value), 0, 1)
integrate(\(x) x/2 * (log(5*x^2 + 2*x + 1) - log(5*x^2)) +
	- x/2 * (atan((x+1)/(2*x)) - atan(1/2)), 0, 1)
2/25 * (atan(1/2) - atan(3)) - 11/200 * pi + 1/4*atan(1/2) +
	+ 3/10 * log(8) - log(5)/4;

#
integrate(\(x) x/2 * (log(5*x^2 + 2*x + 1) - log(5*x^2)), 0, 1)
integrate(\(x) -4/25 / (5*x^2 + 2*x + 1), 0, 1)$value +
	+ 7/25*log(8) - log(5)/4 + 1/10;
2/25 * (atan(1/2) - atan(3)) + 7/25*log(8) - log(5)/4 + 1/10;


### I( atan(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(1 - x*y), 0, 1)$value), 0, 1)
- (5/96*pi^2 + pi*log(2)/8 - log(2)^2 / 8 + log(2)/2 - pi/4 - Catalan)


### I( atan(1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(1 + x*y), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	- x*y / ((x*y + 1)^2 + 1), 0, 1)$value), 0, 1)$value +
	+ 2*atan(2) - 1/2*log(5) + 1/2*log(2) - atan(1);
# TODO

# alternative:
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(1 + x*y) + atan(1-x*y), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x^2*y^2/2), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x^2/2) / y, 0, y)$value), 0, 1)
#
integrate(\(x) atan(x^2/2), 0, 1)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x) -4*x^2 / (x^4 + 4) / y, 0, y)$value), 0, 1)$value;
#
integrate(\(x) - 2*sqrt(2) * Re(atan(x / sqrt(2i)) / sqrt(1i)) / x, 0, 1)$value +
	+ log(5)/2 - atan(3/4);
#
integrate(\(x) (pi/4 - atan(x+1)) / x, 0, 1)$value +
integrate(\(x) atan(3*tan(x)), 0, pi/4)$value +
	+ (pi*log(2)/8 - Catalan) + pi^2 / 32 + log(5)/2 - atan(3/4);
# TODO


### Atan( Fraction )
# Composite-Fraction

### I( atan((y-x)/(x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan((y-x)/(x+y)), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
log(2)/4;

### I( x * atan((y-x)/(x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan((y-x)/(x+y)), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(1 - pi/4) / 6;

### I( y * atan((y-x)/(x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan((y-x)/(x+y)), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
log(2)/6;
# Note: on [0, 1]^2
(pi/4 + log(2) - 1) / 6;


### I( atan((y-x)/(1+x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan((y-x)/(1+x+y)), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
-3/8 * (digamma(1/3) + digamma(2/3) + 2*Euler) - 3/8*pi - log(5)/2 + 3/4*atan(2);
9/8 * log(3) - 3/8*pi - log(5)/2 + 3/4*atan(2);


### I( atan( (x+y)/(1+x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan((x+y)/(1+x*y)), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi*log(2)/2 + 2*log(2) - 2*Catalan;

### I( x * atan((x+y)/(1+x*y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan((x+y)/(1+x*y)), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- pi^2 / 12 + pi/8 - log(2)^2/2 + 3/4 * log(2) + 1/2;

### I( 1/x * atan((x+y)/(1+x*y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan((x+y)/(1+x*y)) / x - atan(y) / x, 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2 / 48 + log(2)^2 / 2;

### I( y/x * atan((x+y)/(1+x*y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan((x+y)/(1+x*y)) * y/x - atan(y) * y/x, 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- pi/4 + 2*log(2) - Catalan/2;


### ATAN( OTHER )

### I( atan(x^2 + y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x^2 + y^2), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
integrate(\(x) 2 * Im(sqrt(x^2-1i) * atan(1/sqrt(x^2-1i))), 0, 1)$value +
	+ atan(2) - Im(sqrt(1+1i) * atan(1/sqrt(1+1i)) - sqrt(1-1i)*atan(1/sqrt(1-1i)));
integrate(\(x) 2 * Im(sqrt(x^2+1i) * atan(sqrt(x^2+1i))), 0, 1)$value +
	+ atan(2) - Im(sqrt(1+1i) * atan(1/sqrt(1+1i)) - sqrt(1-1i)*atan(1/sqrt(1-1i))) +
	- pi/2 * Re(1/sqrt(1i) / cos(atan(sqrt(-1i))) - log((1-sin(atan(sqrt(1i))))/cos(atan(sqrt(1i)))));
integrate(\(x) 2 * Im(sqrt(x^2+1i) * atan(sqrt(x^2+1i))), 0, 1)$value +
	+ atan(2) - Im(sqrt(1+1i) * atan(1/sqrt(1+1i)) - sqrt(1-1i)*atan(1/sqrt(1-1i))) +
	- pi/2 * Re(sqrt(-1-1i) - log((1-1/sqrt(1-1i))*sqrt(1+1i)));
# TODO


###################
###################

### Separate ATAN

### I( atan(x) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x) / (x+y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
-1/32 * pi^2 + 3/8*pi*log(2) - log(2)^2 / 8;

### I( x * atan(x) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x) / (x+y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(pi*log(2)/2 + pi/2 - 5*log(2) + 2*Catalan) / 4;

### I( 1/x * atan(x) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x) / x / (x+y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
integrate(\(x) - log(x) * log(x+1) / (x^2+1), 0, 1)$value +
	+ pi^3 * 3/64 - Catalan * log(2)/2;
# TODO

### I( y/x * atan(x) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y/x * atan(x) / (x+y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan + 1/32 * pi^2 - 3/8*pi*log(2) + log(2)^2 / 8;

### I( (atan(y) - atan(x)) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(y) - atan(x)) / (x+y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- log(2)^2 * 3/4 + pi*log(2)/8 - Re(log(1-1i)^2)/2;

### I( (atan(x)/x - atan(y)/y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x)/x - atan(y)/y) / (x+y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) 3 * log(1-x) * log(x) / (x^2+1), 0, 1)$value +
	- (pi/4)^3 - Catalan * log(2)/2;
# TODO


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


### I( (atan(x/y) - atan(y/x)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x/y) - atan(y/x)) / (x-y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- pi*log(2)/4 + 2*Catalan;


### Type: Atan(x/y)

### I( atan(x/y) / x )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) / x, 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan + pi/4 + log(2)/2;

### I( atan(x/y) * y / x )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * y / x, 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan/2 + 1/4;

### I( atan(x/y) * y^2 / x )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * y^2 / x, 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan/3 + (pi/4 - log(2)/2) / 9 + 1/18;


### I( atan(x/y) / (1-x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(1/y) - atan(x/y)) / (1-x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(-11/12*pi^2 - log(2)^2 - pi*log(2) + 8*Catalan + 2*pi + 4*log(2)) / 8;

### I( atan(x/y) * y / (1-x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(1/y) - atan(x/y)) * y / (1-x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan - pi*log(2)/8 - pi/8 - log(2)/4 + 1/4;


### I( atan(x/y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (x+y), 0, 1)$value), 0, 1)
pi*log(2)/2

### I( atan(x/y) * abs(x-y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * abs(x-y) / (x+y), 0, 1, rel.tol=1E-9)$value), 0, 1, rel.tol=1E-9)
# trick for numerical accuracy:
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * abs(x-y) / (x+y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * abs(x-y) / (x+y), y, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)$value;
pi*log(2)/2 - pi/4;


### I( atan(x/y) / (1+x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (1+x+y), 0, 1)$value), 0, 1)
(log(3) * 3/4 - log(2)) * pi;


### I( atan(x/y) / (2 - x - y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (2 - x - y), 0, 1)$value), 0, 1)
pi*log(2)/2;

### I( x * atan(x/y) / (2 - x - y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x/y) / (2 - x - y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13);
(7*pi*log(2) - 3*pi - 10*log(2) + 8*Catalan) / 8;

### I( y * atan(x/y) / (2 - x - y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(x/y) / (2 - x - y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13);
(pi * log(2)/2 + pi/2 + 5*log(2) - 4*Catalan) / 4;

### I( x*y * atan(x/y) / (2 - x - y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * atan(x/y) / (2 - x - y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13);
(log(2) - 1/4) * pi / 3;

### I( x^2 * atan(x/y) / (2 - x - y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * atan(x/y) / (2 - x - y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13);
(17*pi*log(2) - 10*pi - 32*log(2) + 2) / 12 + 2*Catalan;

### I( 1/x * atan(x/y) / (2 - x - y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/x * atan(x/y) / (2 - x - y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13);
# TODO

### I( y/x * atan(x/y) / (2 - x - y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y/x * atan(x/y) / (2 - x - y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13);
# TODO


### Div: (x+y)^2

### I( atan(x/y) * x / (x+y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * x / (x+y)^2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-13)
pi/4 +  pi * log(2)/8 - log(2)/2;

### I( atan(x/y) * y / (x+y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * y / (x+y)^2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-13)
3/8*pi*log(2) - pi/4 + log(2)/2;

### I( atan(x/y) * x*y / (x+y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * x*y / (x+y)^2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-13)
pi*log(2)/4 - pi/8;


### I( atan(x/y) * abs(x-y) / (x+y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * abs(x-y) / (x+y)^2, 0, 1, rel.tol=1E-9)$value), 0, 1, rel.tol=1E-9)
# trick for numerical accuracy:
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * abs(x-y) / (x+y)^2, 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * abs(x-y) / (x+y)^2, y, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)$value;
pi*(1 - log(2)) / 2;


### Pow = 2

### Base: I( atan(x/y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y)^2, 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
(pi*log(2) + pi^2/8 - 2*Catalan) / 2;

### I( atan(x/y)^2 / x )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y)^2 / x, 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pi*Catalan / 2 - Catalan + pi^2 / 16 + 3/4 * pi*log(2) - 7/8 * pracma::zeta(3);

### I( atan(x/y)^2 / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y)^2 / (x+y), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pi^2 * log(2) * 3/16 + pi*Catalan/2 - 21/16*pracma::zeta(3);

### I( x * atan(x/y)^2 / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x/y)^2 / (x+y), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pi^2 * log(2) / 16 + ((pi/4)^2 + pi*log(2)/4 - Catalan) / 2;

### I( y * atan(x/y)^2 / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(x/y)^2 / (x+y), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
(pi*log(2) + pi^2/8 - 2*Catalan) / 2 +
	- (pi^2 * log(2) / 16 + ((pi/4)^2 + pi*log(2)/4 - Catalan) / 2);
- (pi^2 * log(2) / 16 - ((pi/4)^2 + 3/4 * pi*log(2) - Catalan) / 2);


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


###################

### Type: Atan(X*Y)


### I( atan(x*y) / x )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) / x, 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- pi/4 + log(2)/2 + Catalan;

### I( atan(x*y) * y/x )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * y/x, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
Catalan/2 - pi/8 + 1/4;

### I( atan(sqrt(x*y)) / x )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sqrt(x*y)) / x, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
2*Catalan - pi/2 + 1;

### I( atan(sqrt(x*y)) * y/x )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sqrt(x*y)) * y/x, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
Catalan - 1/6;


### I( (atan(x*y) - atan(y)) / (1-x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x*y) - atan(y)) / (1-x), 0, 1)$value), 0, 1)
pi^2/32 + pi*log(2)/8 - Catalan - log(2)^2/8;

### I( atan(sqrt(x*y)) / (1-x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(sqrt(x*y)) - atan(sqrt(y))) / (1-x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi*log(2)/2 - 2*log(2);

### I( atan(x*y) / ((1-x)*(1-y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x*y) - atan(x) - atan(y) + atan(1)) / (1-y) / (1-x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi^3 / 3 + pi*log(2)^2) / 16 - Catalan * log(2)/2;


### Div: (x+y)

### I( atan(x*y) / (x*(x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) / (x*(x+y)), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(pi/4)^3;


### I( atan(x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (x+y), 0, 1)$value), 0, 1)
# Derivation: z = x/y;
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(y^2*x) / (x+1), 0, 1/y)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x) 2 * atan(y^2*x) / (x+1), 0, 1)$value), 0, 1)
# Factor: *2;
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x^2*y) / (y+1), 0, 1)$value), 0, 1)
integrate(\(x) atan(x)/(x+1) - 2/(x*(x+1)) * sapply(x, \(y)
	integrate(\(x) x^2 / (x^4 + 1/y^2), 0, 1)$value), 0, 1);
isq = exp(-1i * pi/4); isq2 = exp(1i * pi/4);
integrate(\(x) atan(x)/(x+1) +
	- sqrt(x)/(x*(x+1)) * Re(atan(sqrt(x)/isq)/isq + atan(sqrt(x)/isq2)/isq2), 0, 1)
integrate(\(x) -4 * Re(atan(x/isq)/isq) / (x^2+1), 0, 1)$value + pi * log(2)/8;
integrate(\(x) 4 * Re(1/(x^2 + 1i)) * atan(x), 0, 1)$value +
	(digamma(3/8) - digamma(7/8)) * pi / 8 + pi * log(2)/8;
integrate(\(x) 4 * x^2/(x^4 + 1) * atan(x), 0, 1)$value +
	(digamma(3/8) - digamma(7/8)) * pi / 8 + pi * log(2)/8;
# TODO

# Helper:
isq = exp(-1i * pi/4); isq2 = exp(1i * pi/4);
integrate(\(x) Re(atan(x/isq)/isq + atan(x/isq2)/isq2) / (x^2+1), 0, 1)
integrate(\(x) Re(atan(x/isq)/isq) / (x^2+1) * 2, 0, 1)
integrate(\(x) (x^2 - 1) * (log(1-x) - Re(log(x + 1i))) / (x^4 + 1), 0, 1)$value +
	- (digamma(7/8) - digamma(3/8)) * pi / 32 +
	- (pracma::psi(1, 3/8) - pracma::psi(1, 7/8)) / 64;
# see Integrals.Log.Fractions.P4.R;
integrate(\(x) (1 - x^2) * Re(log(x + 1i)) / (x^4 + 1), 0, 1)$value +
	- (pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 64 +
	+ (pracma::psi(1, 7/8) - pracma::psi(1, 3/8)) / 32 +
	- (digamma(5/8) - digamma(1/8)) * log(2) / 32 +
	+ (digamma(7/8) - digamma(3/8)) * log(2) / 32 +
	- (digamma(5/8) - digamma(1/8)) * pi / 64 +
	- (digamma(7/8) - digamma(3/8)) * pi * 3 / 64;

# TODO


### I( atan(x*y) * x / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * x / (x+y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- pi^2 / 96 + pi/8 - log(2)/4;

### I( atan(x*y) / x / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) / x / (x + y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi/4)^3;

### I( atan(x*y) * y / x / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * y / x / (x + y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan - pi/4 + log(2)/2 - 0.1936317534825;
# TODO: I( atan(x*y) / (x+y) )

### I( atan(x*y) * y^2 / x / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * y^2 / x / (x + y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan/2 + pi^2 / 96 - pi/4 + log(2)/4  + 1/4;


### I( atan(x*y) * abs(x-y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	2 * atan(x*y) * abs(x-y) / (x+y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
-1/96*pi^2 + pi*log(2)/4 - pi/4 + log(2)^2 * 5/8 + log(2)/2 + Re(log(1-1i)^2)/2;
-1/24*pi^2 + pi*log(2)/4 - pi/4 + log(2)^2 * 3/4 + log(2)/2;


### I( atan(x*y) / (x + y + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (x+y+1), 0, 1)$value), 0, 1)
# TODO


### I( atan(x*y) / (2 - x - y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (2 - x - y), 0, 1)$value), 0, 1)
#
integrate(\(x) sapply(x, \(y) integrate(\(x) 2 * atan(x*(y-x)) / (2 - y), 0, y/2)$value), 0, 1)$value +
integrate(\(x) sapply(x, \(y)
	integrate(\(x) -2*atan(x*(y-x)) / (2 - y), 0, y-1)$value +
	integrate(\(x) +2*atan(x*(y-x)) / (2 - y), 0, y/2)$value ), 1, 2)$value;
#
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*(y-x)) / (2 - y), 0, y)$value), 0, 1)$value +
integrate(\(x) sapply(x, \(y)
	integrate(\(x) -2*atan(x*(y-x)) / (2 - y), 0, y-1)$value +
	integrate(\(x) atan(x*(y-x)) / (2 - y), 0, y)$value ), 1, 2)$value;
#
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*(y-x)) / (2 - y), 0, y)$value), 0, 1)$value +
integrate(\(x) sapply(x, \(y)
	integrate(\(x) 2*atan(x*(y-x)) / (2 - y), 0, 1)$value +
	integrate(\(x) - atan(x*(y-x)) / (2 - y), 0, y)$value ), 1, 2)$value;
# TODO

### I( x * atan(x*y) / (2 - x - y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x*y) / (2 - x - y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
0.5238047434009 + pi^2 / 96 - pi/8 + log(2)/4;
# TODO: I( atan(x*y) / (2 - x - y) )


### Div: includes x*y

### I( atan(x*y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (1 + x*y), 0, 1)$value), 0, 1)
((digamma(1/4) + Euler) * pracma::psi(1, 1/4) - pracma::psi(2, 1/4)/2) / 16 +
	- 7/4 * pracma::zeta(3) - (pi/4)^3 + 3/16 * pi^2 * log(2) +
	+ Catalan * pi/4 + Catalan * log(2);
(pi/4)^3 - Catalan * log(2)/2;


# see also: I( atan(x*y) / (1 - x*y) )
integrate(\(x) -1/3 * atan(x) * log(1-x) / x, 0, 1)$value +
	integrate(\(x) -1/3 * log(1-x) * log(x) / (x^2+1), 0, 1)$value +
	integrate(\(x) 1/3 * log(x) * log(x^2+1) / (x^2+1), 0, 1)$value;
integrate(\(x) -1/3 * atan(x) * log(1-x) / x, 0, 1)$value +
integrate(\(x) +1/3 * log(x) * log(1-x) / (x^2+1), 0, 1, rel.tol=1E-12)$value - Catalan * log(2)/3;

### I( x * atan(x*y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x*y) / (1 + x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi*log(2)/4 - pi/4 + log(2)/2;

### I( 1/x * atan(x*y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/x * atan(x*y) / (1 + x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- pi*log(2)/4 + Catalan;

### I( y/x * atan(x*y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y/x * atan(x*y) / (1 + x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan/2 - pi/8 + log(2)/4;


### I( atan(x*y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (1 - x*y), 0, 1)$value), 0, 1)
integrate(\(x) - atan(x) * log(1-x) / x, 0, 1)$value +
	integrate(\(x) - log(1-x) * log(x) / (x^2+1), 0, 1)$value;
integrate(\(x) - atan(x) * log(x) / (1-x), 0, 1, rel.tol=1E-12); # alternative
#
integrate(\(x) -2 * log(1-x) * log(x) / (x^2+1), 0, 1)$value +
	+ pi^3 * 3/64 - Catalan * log(2)/2;
# TODO

### I( x * atan(x*y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x*y) / (1 - x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi/4 - log(2)/2;

### I( x^2 * atan(x*y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * atan(x*y) / (1 - x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- log(2)/4 + pi/4 - 1/4;

### I( 1/x * atan(x*y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/x * atan(x*y) / (1 - x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan;

### I( y/x * atan(x*y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y/x * atan(x*y) / (1 - x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi/8 - log(2)/4 + Catalan / 2;

### I( y^2/x * atan(x*y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y^2/x * atan(x*y) / (1 - x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(2*Catalan + pi - log(2) - 1) / 6;


### I( atan(x*y) / (1 - x + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) / (1 - x + x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^3 / 32 - Catalan * log(2);

### I( x * atan(x*y) / (1 - x + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(x*y) / (1 - x + x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- pi^2 / 48 + pi/4 - log(2)/2;

### I( y * atan(x*y) / (1 - x + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * atan(x*y) / (1 - x + x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^3 / 32 - pi^2 / 16 - pi/4 + log(2)/2 - Catalan * log(2) + Catalan;

### I( 1/x * atan(x*y) / (1 - x + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/x * atan(x*y) / (1 - x + x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2 / 16;

### I( y/x * atan(x*y) / (1 - x + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y/x * atan(x*y) / (1 - x + x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan/2 + pi^2 / 32 - pi/8;

### I( 1/y * atan(x*y) / (1 - x + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/y * atan(x*y) / (1 - x + x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^3 / 32;

### I( x/y * atan(x*y) / (1 - x + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x/y * atan(x*y) / (1 - x + x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^3 / 32 - pi^2 / 48 + pi/2 - log(2) - Catalan;


### I( atan(x*y) / (1 + x - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) / (1 + x - x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO: same as I( log(1 + x - y) / (x^2 + y^2) )


### Prod

### I( atan(x*y) * log(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * log(x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2 / 48 - pi/4 + log(2)/2 + 3/32 * pracma::zeta(3);

### I( atan(x*y) * log(1-x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * log(1-x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
-23/64 * pracma::zeta(3) + Catalan * pi / 4 - Catalan +
	+ 5/96*pi^2 + pi*log(2)/8 - log(2)^2 / 8 + log(2)/2 - pi/4;

### I( atan(x*y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * atan(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) - log(1-x) * log(x) / (x^2+1), 0, 1)$value +
	+ pi*log(2)/4 + (pi/4)^3 + (pi/4)^2 - Catalan * log(2)/2 - Catalan;
# TODO


### Div

### I( atan(x/y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (1+x*y), 0, 1)$value), 0, 1)
pi^3 / 48

### I( atan(x/y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (1-x*y), 0, 1)$value), 0, 1)
pi^3 / 24;

### I( atan(x/y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (x+y), 0, 1)$value), 0, 1)
pi*log(2)/2

### I( atan(x/y) / (2 - (x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (2-x-y), 0, 1)$value), 0, 1)
pi*log(2)/2

### I( atan(x/y) / (x+y+1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (x+y+1), 0, 1)$value), 0, 1)
(log(3) * 3/4 - log(2)) * pi

### I( atan(x/y) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x/y) - pi/4) / (x-y), 0, y, rel.tol=1E-6)$value), 0, 1, rel.tol=1E-6)
Catalan - pi*log(2)/8;
# Note: on [y, 1] = identical;


### Div: Higher

### I( x*y * atan(x) / (x^2+y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * atan(x) / (x^2+y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi*log(2)/8 - pi/8 - log(2)/4 + Catalan/2;

### I( x*y * atan(x) / (x^2*y^2 + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * atan(x) / (x^2*y^2 + 1), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) log(x) * log(1-x) / (x^2+1), 0, 1)$value +
	- (pi/4)^3 + Catalan * log(2)/2;
# TODO


### Div: SQRT(...)

### I( atan(x/y) / sqrt(1-x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) / sqrt(1-x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi * (1-log(2));


### I( atan(x/y) / sqrt(1+x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) / sqrt(1+x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi * (log(2) - asinh(1)) + (digamma(5/8) - digamma(3/8));


### I( atan(x/y) / sqrt(x^2 + y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) / sqrt(x^2 + y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi/2 * asinh(1);

### I( x*y * atan(x/y) / sqrt(x^2 + y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * atan(x/y) / sqrt(x^2 + y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(sqrt(2) - 1) * pi/6;


### ATAN(x+y)

### I( atan(x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y), 0, 1)$value), 0, 1)
3/2*atan(2) - log(5) + log(2)

### I( atan(x+y) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y) / (x + y), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x) / x, y, y+1)$value), 0, 1)
#
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x) / x, y, 1)$value), 0, 1)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x) / x, 1, y+1)$value), 0, 1)$value;
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x) / x, 1, y+1)$value), 0, 1)$value +
	+ pi/4 - log(2)/2;
# alternative:
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x) / x, y, y+1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x) / x, 0, y)$value), 0, 2)$value +
	- (Catalan - pi/4 + log(2)/2) * 2;
integrate(\(x) atan(x) / x * (2-x), 0, 2)$value - (Catalan - pi/4 + log(2)/2) * 2;
integrate(\(x) 2 * atan(x) / x, 0, 2)$value +
	- (Catalan - pi/4 + log(2)/2) * 2 - 2*atan(2) + log(5)/2;
# TODO

### I( atan(x+y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y) / (1-x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan((x+y)/2) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan((x+y)/2) / (1-x*y), 0, 1)$value), 0, 1)
# TODO

### I( x * atan((x+y)/2) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) x * atan((x+y)/2) / (1-x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y) / (1+x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan((x+y)/2) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan((x+y)/2) / (1+x*y), 0, 1)$value), 0, 1)
# TODO

### I( x * atan((x+y)/2) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) x * atan((x+y)/2) / (1+x*y), 0, 1)$value), 0, 1)
# TODO



### I( atan(x/y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) / (1 - x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi^3 / 24


################

#############
### Varia ###

### I( log(beta(x, y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(beta(x, y)), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
log(2*pi)/2 + 3/4;


### I( atan(beta(x, y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(beta(x, y)), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO

