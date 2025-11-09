#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
##
## v.0.3a

### Double Integrals

### Examples:
# I( 1 / (1 - x - x*y) )
# I( 1 / (1 + x - x*y) )
# I( 1 / (2 - x^2 - y^2) )
# I( x*y / ((x*y - 1)^2 + 1) )


### History

# - [refactor] moved Log-Type to file:
#   Integrals.Double.Log.R;
# - [refactor] moved Prod(Log)-Type to file:
#   Integrals.Double.Log.Prod.R;
# - [refactor] moved Exp-Type to file:
#   Integrals.Double.Exp.R;
# - [refactor] moved Atan-Type to file:
#   Integrals.Double.Atan.R;


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;
dzeta2  = - 0.937548254316;

#####################
#####################

### Simple Fractions

### I( 1 / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / (1-x*y), 0, 1)$value), 0, 1)
pi^2 / 6

### I( 1 / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / (1+x*y), 0, 1)$value), 0, 1)
pi^2/12

### I( 1 / (1 - x^2*y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / (1-x^2*y^2), 0, 1)$value), 0, 1)
pi^2/8

### I( 1 / (1 + x^2*y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / (1+x^2*y^2), 0, 1)$value), 0, 1)
Catalan


### I( 1 / (1 - x + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / (1 - x + x*y), 0, 1)$value), 0, 1)
pi^2 / 6

### I( 1 / (1 + x + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / (1 + x + x*y), 0, 1)$value), 0, 1)
integrate(\(x) log(1+x) / x, 1, 2)
pracma::polylog(1/4, 2) / 2 - log(2)^2;
# TODO

integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / (1 + x - x*y), 0, 1)$value), 0, 1)
pi^2/12


### on [0, 1] x [0, 1/2]
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / (1 - x - x*y), 0, 1/2)$value), 0, 1)
pi^2 / 12 + log(2)^2/2;


### I( x*y / ((x*y + 1)^2 + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) x*y / ((x*y + 1)^2 + 1), 0, 1)$value), 0, 1)
integrate(\(x) - log(x) * x / (x^2 + 2*x + 2) , 0, 1);
# TODO

# Derivation:
# z = x*y =>
integrate(\(x) 1/x * sapply(x, \(y) integrate(\(x) (x-1) / (x^2 + 1), 1, 1+y)$value), 0, 1)
integrate(\(x) 1/x * (log((x+1)^2+1)/2 - log(2)/2 - atan(x+1) + atan(1)), 0, 1)

# [old] [longer way]
integrate(\(x) sapply(x, \(y) integrate(\(x) Im((1+1i)/ (x*y + 1 + 1i)), 0, 1)$value), 0, 1)
integrate(\(x) pi/4/x - atan(x+1) / x, 0, 1)$value +
integrate(\(x)  Re(log(x + 1+1i) - log(1+1i)) / x, 0, 1)$value;
integrate(\(x) (Re(log(x + 1+1i)) - atan(x+1) + pi/4 - log(2)/2) / x, 0, 1);
integrate(\(x) (Re((1-1i)*log(x + 1+1i)) - pi/4 - log(2)/2) / x, 0, 1);
integrate(\(x) 1/2 * (log(x^2 + 2*x + 2) - log(2)) / x, 0, 1)$value +
integrate(\(x) - (atan(x+1) - pi/4) / x, 0, 1)$value;
integrate(\(x) - log(x) * x / (x^2 + 2*x + 2) , 0, 1);

# TODO


### I( x*y / ((x*y - 1)^2 + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) x*y / ((x*y - 1)^2 + 1), 0, 1)$value), 0, 1)
- (5/96*pi^2 + pi*log(2)/8 - log(2)^2 / 8 - Catalan);

# Derivation:
integrate(\(x) sapply(x, \(y) integrate(\(x) Im(-(1-1i)/ (x*y - 1 + 1i)), 0, 1)$value), 0, 1)
integrate(\(x) Im(-(1-1i)*(log(x - 1+1i) - log(-1+1i))) / x, 0, 1)
integrate(\(x) (Im(-(1-1i)*log(x - 1+1i)) + 3/4*pi - log(2)/2) / x, 0, 1)
- (5/96*pi^2 + pi*log(2)/8 - log(2)^2 / 8 - Catalan);


### I( 1 / (2 - x^2 - y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (2 - x^2 - y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
Catalan;

### I( x / (2 - x^2 - y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / (2 - x^2 - y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
log((2+sqrt(2))/(2-sqrt(2))) / sqrt(2) - log(2);
sqrt(2) * log(1+sqrt(2)) - log(2);


################
### Radicals ###

### I( (x^p + y^p) / (1-x*y) / (x*y)^q )
# Maths 505: A surprisingly interesting integral
# https://www.youtube.com/watch?v=BNzVNcZrbfA
# Note: Series expansion of 1/(1 - x*y)

integrate(\(x) sapply(x, \(y)
	integrate(\(x) (sqrt(x) + sqrt(y)) / (1-x*y) / (x*y)^(1/4), 0, 1)$value), 0, 1)
(digamma(5/4) - digamma(3/4)) * 4


### Gen: I( (x^p + y^p) / (1-x*y) / (x*y)^q )
p = sqrt(3); q = sqrt(5) - 2;
integrate(\(x) sapply(x, \(y)
	integrate(\(x) (x^p + y^p) / (1-x*y) / (x*y)^q, 0, 1)$value), 0, 1)
(digamma(1+p-q) - digamma(1-q)) * 2 / p


### Type: SQRT
# - see also file: Integrals.Double.Radicals.R;

### I( 1 / sqrt(x^2 + y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / sqrt(x^2 + y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
2*log(sqrt(2) + 1);

### I( x / sqrt(x^2 + y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / sqrt(x^2 + y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
log(sqrt(2) + 1) / 2 + sqrt(2)/2 - 1/2;


###############
###############

### Varia


### Log( GAMMA )

### I( log(gamma(x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(gamma(x+y)), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
log(2*pi)/2 - 3/4;


### I( log(gamma(x*y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(gamma(x*y)), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) - log(x) * log(gamma(x)), 0, 1, rel.tol=1E-12)
# TODO


### I( log(gamma(sqrt(x*y))) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(gamma(sqrt(x*y))), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( log(gamma(sqrt(x^2 + y^2))) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(gamma(sqrt(x^2 + y^2))), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### GAMMA(...)

### I( gamma(x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	gamma(x*y) - 1/(x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( gamma(sqrt(x^2 + y^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	gamma(sqrt(x^2 + y^2)), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO

