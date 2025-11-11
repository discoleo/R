#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
## Type: PROD( Log )
##
## v.0.2a

### Double Integrals

### Examples:
# I( log(x) * log(1 - x*y) )
# I( log(1-x) * log(1 - x*y) )
# I( x/y * log(1-x) * log(1 - x*y) )
# I( y/x * log(1-x) * log(1 - x*y) )


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;
dzeta2  = - 0.937548254316;

#####################
#####################

### PROD( LOG )

### I( log(x) * log(y) * log(1 - x*y) )
# Maths 505: A mesmerizing result
# https://www.youtube.com/watch?v=QqVhd_xfjnc
# Note: series expansion of log(1 - x*y);


### Simple: I( log(x) * log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1-x*y) * log(x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
3 - sum(pracma::zeta(2:3));

### Simple: I( log(x) * log(1 - x + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1 - x + x*y) * log(x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
3 - sum(pracma::zeta(2:3));

### Simple: I( log(y) * log(1 - x + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1 - x + x*y) * log(y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
3 - 2*pracma::zeta(3);

### Simple: I( log(x) * log(1 + x - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1 + x - x*y) * log(x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- 3/4 * pracma::zeta(3) - pi^2 / 12 - 2*log(2) + 3;

### Simple: I( log(y) * log(1 + x - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1 + x - x*y) * log(y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (5/8 * pracma::zeta(3) + pi^2 / 6 - log(2)^2 + 2*log(2) - 3);


### Non-Simple:

### I( log(1-x) * log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1-x) * log(1-x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
3 - 2*pracma::zeta(3);

### Type: Plus x: I( log(1+x) * log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1+x) * log(1-x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pracma::zeta(3) * 5/8 - pi^2/6 + log(2)^2 - 4*log(2) + 3;

### Type: Plus xy: I( log(1-x) * log(1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1-x) * log(1+x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (pracma::zeta(3) * 5/8 + pi^2/6 - log(2)^2 + 2*log(2) - 3);

### I( log(1+x) * log(1+x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1+x) * log(1+x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pracma::zeta(3) / 4 + 2*log(2)^2 - 6*log(2) + 3;

# Derivation:
integrate(\(x) (2*log(2) - 1) * log(1+x) - ((1+x)*log(1+x)-x) *
	sapply(x, \(y) integrate(\(x)
	x / (1+x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
integrate(\(x) (2*log(2) - 1) * log(1+x) - ((1+x)*log(1+x)-x) / x +
	+ ((1+x)*log(1+x)-x) * log(1+x) / x^2, 0, 1, rel.tol=1E-13)
integrate(\(x) 2*(log(2) - 1) * log(1+x) - 2*log(1+x) / x + 1 +
	+ (1+x)/x^2 * log(1+x)^2, 0, 1, rel.tol=1E-13);
pracma::zeta(3) / 4 + 2*log(2)^2 - 6*log(2) + 3;

### I( log(1 - x*y) * log(1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1-x*y) * log(1+x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
-21/8 * pracma::zeta(3) + pi^2 * log(2) / 2 - 5/12*pi^2 + log(2)^2 - 4*log(2) + 6;

### Helper: I( x*y * log(1 - x*y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y*log(1-x*y) / (1+x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
13/8 * pracma::zeta(3) - pi^2 * log(2) / 4 + pi^2/6 - 2;


### Variants:

### I( x * log(1-x) * log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * log(1-x) * log(1-x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
1/2;

### I( y * log(1-x) * log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * log(1-x) * log(1-x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
1 + 1/4 - pi^2/12;

### I( log(1-x) * log(1 - x*y) / x )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/x * log(1-x) * log(1-x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
2*pracma::zeta(3) - pi^2/6;

### I( log(1-x) * log(1 - x*y) / y )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/y * log(1-x) * log(1-x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2/6 + 2*pracma::zeta(3) - 3;

### I( x/y * log(1-x) * log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x/y * log(1-x) * log(1-x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2/8 + pracma::zeta(3) - 3/2 - 1/16;

### I( y/x * log(1-x) * log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y/x * log(1-x) * log(1-x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pracma::zeta(3) - pi^2/24 - 1/4;


### Variant: Log(1+x)

### I( x * log(1+x) * log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * log(1+x) * log(1-x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (log(2)^2 - log(2) - pi^2/6 + 2);

### I( y * log(1+x) * log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * log(1+x) * log(1-x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (pi^2/6 - 2*log(2)^2 + 3*log(2) - 5/2) / 2;

### I( x/y * log(1+x) * log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x/y * log(1+x) * log(1-x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
-5/16 * pracma::zeta(3) - pi^2 / 24 + log(2)/2 + 3/16;

### I( y/x * log(1+x) * log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y/x * log(1+x) * log(1-x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
-5/16 * pracma::zeta(3) - pi^2 / 24 + log(2)/2 + 1/4;


### Simple

### I( log(x) * log(y) * log(1 - x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x) * log(y) * log(1 - x*y), 0, 1, rel.tol=1E-12)$value
	), 0, 1, rel.tol=1E-12)
pracma::zeta(2) + pracma::zeta(3) + pracma::zeta(4) - 4;


### I( log(x) * log(y) * log(z) * log(1 - x*y*z) )
# Note: rather slow (5 - 10 s);
integrate(\(x) sapply(x, \(y)
	integrate(\(x) sapply(x, \(z)
	integrate(\(x) log(x) * log(y) * log(z) * log(1 - x*y*z),
		0, 1, rel.tol=1E-12)$value
	), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
6 - sum(pracma::zeta(2:6));

