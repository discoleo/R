#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
##
## v.0.2c

### Double Integrals

### Examples:
# I( (x+y)*log(x+y) / (1+x*y) )


### History

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

### Simple

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


###########

###########
### Log ###

### Simple Logs

### I( log(abs(x+y-1)) )
# Note: numerical issues;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(abs(x+y-1)) , 0, y)$value), 0, 1)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x) log(abs(x+y-1)) , y, 1)$value), 0, 1)$value;
- 3/2;


### I( x * log(x) * log(y)^2 / (1 - x*y) )
# Maths 505: A very interesting double integral with a beautiful result
# https://www.youtube.com/watch?v=lWVBnWCz63I
# - series expansion of 1/(1 - x*y);

### I( x * log(x) * log(y)^2 / (1 - x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) x * log(x) * log(y)^2 / (1 - x*y), 0, 1)$value), 0, 1)
# - 2 * sum(1/j^3 * 1/(j+1)^2); j >= 1;
- 2*pracma::zeta(3) + 6*pracma::zeta(2) - 8;


### I( x * log(x)^2 * log(y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) x * log(x)^2 * log(y) / (1 - x*y), 0, 1)$value), 0, 1)
# - 2 * sum(1/j^2 * 1/(j+1)^3); j >= 1;
- 2*pracma::zeta(3) - 6*pracma::zeta(2) + 12;


### I( x^2 * log(x)^2 * log(y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) x^2 * log(x)^2 * log(y) / (1 - x*y), 0, 1)$value), 0, 1)
# - 2 * sum(1/j^2 * 1/(j+2)^3); j >= 1;
- pracma::zeta(3)/2 - 3/4*pracma::zeta(2) + 3/8 * (1+1/2) + 1/2*(2+1/4+1/8);
- pracma::zeta(3)/2 - 3/4*pracma::zeta(2) + 1 + 3/4;


# TODO: more variants;


### I( log(x)*log(y) / (1 - x*y)^3 )
# Maths 505: A lovely surprise awaits at the end of the video
# https://www.youtube.com/watch?v=24CCdCBF5ck
# Series expansion of 1/(1 - x*y)^3;

integrate(\(x) sapply(x, \(y) integrate(\(x) log(x)*log(y) / (1 - x*y)^3, 0, 1)$value), 0, 1)
(pracma::zeta(2) + pracma::zeta(3)) / 2


### Mixed Logs

### I( log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1-x*y), 0, 1)$value), 0, 1)
pi^2/6 - 2

### I( log(1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1+x*y), 0, 1)$value), 0, 1)
pi^2/12 + 2*log(2) - 2


### I( log(1 + x + x*y) )
# Maths 505: 
# https://www.youtube.com/watch?v=xiSw7Ou4RIY

### I( log(1 + x + x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(1 + x + x*y), 0, 1)$value), 0, 1)
3*log(3) - 2*log(2) + polylog2(-1) - polylog2(-2) - 2;
3*log(3) - 2*log(2) - pi^2/12 - polylog2(-2) - 2;


### I( log(1 - x + x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(1 - x + x*y), 0, 1)$value), 0, 1)
pi^2/6 - 2


### I(x, y, z)
# Fails numerically: but trivial extension of above integral;
# library(Rmpfr)
integrate(\(x) sapply(x, \(y)
	integrate(\(x) sapply(x, \(z)
	integrate(\(x) {
		x = mpfr(x, 240); y = mpfr(y, 240);
		as.numeric(log(x) + log(1 + y + y*z));
	}, 0, 1, rel.tol=1E-12)$value),
	0, 1, rel.tol=1E-10)$value), 0, 1)
3*log(3) - 2*log(2) - pi^2/12 - polylog2(-2) - 3;


### I( log(1 - x^2*y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1-x^2*y^2), 0, 1)$value), 0, 1)
pi^2/4 + 2*log(2) - 4

### I( log(1 + x^2*y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1+x^2*y^2), 0, 1)$value), 0, 1)
integrate(\(x) - log(x) * log(x^2+1), 0, 1)
pi/2 + log(2) + 2*Catalan - 4;


### Log-Fractions

### I( log(1 + x+y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1+x+y) / (x+y), 0, 1)$value), 0, 1)
# Note:
# - polylog2 is in file Integrals.Polylog.Helper.R;
integrate(\(x) sapply(x, \(x) - Re(polylog2(-x-1))), 0, 1)$value +
	- pi^2/12 + 2*log(2) - 1;
- 2*polylog2(-2) - pi^2/6 - 3*log(3) + 4*log(2);


### I( log(1-x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1-x*y) / (x+y),
	0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
log(2)^2 - 4*log(2) - pi^2/6 + 4*Catalan;


### I( log(1+x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1+x*y) / (x+y),
	0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi^2 / 4 + log(2)^2 - 4*log(2);


### Div: (1 +/- x*y)

### I( log(x) / (1 - x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x) / (1-x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- pracma::zeta(3)

### I( log(x) / (1 + x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x) / (1+x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
-3/4 * pracma::zeta(3)


### I( log(x+y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(x+y) / (1 - x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-12)
1/8 * pracma::zeta(3);


# I( log((x+y)/2) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log((x+y)/2) / (1 - x*y), 0, 1)$value), 0, 1)
1/8 * pracma::zeta(3) - log(2) * pi^2 / 6;


### I( log(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(x+y) / (1 + x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- pracma::zeta(3) / 8


### I( (x+y) * log(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) {
	integrate(\(x) (x+y)*log(x+y) / (1+x*y), 0, 1, rel.tol=1E-12)$value
	}), 0, 1, rel.tol=1E-12)
- pi^2/2 + 2*log(2)^2 + 4;


### I( (x+y) * log(x+y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) {
	integrate(\(x) (x+y)*log(x+y) / (1-x*y), 0, 1, rel.tol=1E-12)$value
	}), 0, 1, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y) integrate(\(x) 2 * log(1-x*y) / (x+y),
	0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)$value +
	- 2*log(2)^2 + 4*log(2) + pi^2/3 - 4;
8*Catalan - 4*log(2) - 4;


### I( log(abs(x+y-1)) / (1 - x*y) )
# Note: numerical issues; very slow with mpfr;
FUN = \(x, y) log(abs(x+y-1)) / (1-x*y);
# FUN = \(x, y) { x = mpfr(x, 128); y = mpfr(y, 128); as.numeric(log(abs(x+y-1)) / (1-x*y)); }
integrate(\(x) sapply(x, \(y) integrate(FUN, 0, y, y=y, rel.tol=1E-6)$value), 0, 1, rel.tol=1E-6)$value +
integrate(\(x) sapply(x, \(y) integrate(FUN, y, 1, y=y, rel.tol=1E-6)$value), 0, 1, rel.tol=1E-6)$value;
-5/3 * pracma::zeta(3);


### I( log(x+y+1) / (1 + x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x+y+1) / (1 + x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### Div: (x^2 + y^2)

### I( (x+y)*log(x+y) / (x^2 + y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(x+y)*log(x+y) / (x^2 + y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
1/48*pi^2 + pi*log(2)/4 + log(2)^2 / 4 - log(2) - pi/2;


### I( log(x+y+1) / (x^2 + y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(x+y+1) / (x^2 + y^2), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### I( log(1 + x) / (x^2 + y^2) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(1 + x) / (x^2 + y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# I( log(1 + x*y) / (x^2 + y^2) ) + pi^3 / 48;
integrate(\(x) - 3*atan(x) * log(1-x) / x, 0, 1)$value - 5/96 * pi^3;


### I( log(1 + x*y) / (x^2 + y^2) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(1 + x*y) / (x^2 + y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) - 3*atan(x) * log(1-x) / x, 0, 1)$value - 7/96 * pi^3;
integrate(\(x) log(x) * log(x+1) / (x^2+1), 0, 1)$value +
	+ (pi/4)^3 / 3 + Catalan * log(2)/2;
# TODO

### I( log(1 - x*y) / (x^2 + y^2) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(1 - x*y) / (x^2 + y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) - atan(x) * log(1-x) / x, 0, 1)$value - pi^3 / 24;
# TODO


### I( log(x+y) / (1 + x^2 + y^2) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x+y) / (1 + x^2 + y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( log(x) / (2 - x^2 - y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(x) / (2 - x^2 - y^2), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
integrate(\(x) log(1-x) * log(x) / (x^2+1), 0, 1)$value - pi^3 / 32;


### I( log(x+y) / (2 - x^2 - y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(x+y) / (2 - x^2 - y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
#
integrate(\(x) -2 * log(x) * log(1+x) / (x^2+1), 0, 1, rel.tol=1E-12)$value +
	integrate(\(x) - log(x) * log(1-x) / (x^2+1), 0, 1, rel.tol=1E-12)$value +
	((digamma(1/4) + Euler) * pracma::psi(1, 1/4) - pracma::psi(2, 1/4)/2) / 16 +
	- 7/4 * pracma::zeta(3) + 3/16 * pi^2 * log(2) + pi/4 * Catalan;
integrate(\(x) 6 * atan(x) * log(1-x) / x, 0, 1)$value +
	integrate(\(x) - log(x) * log(1-x) / (x^2+1), 0, 1, rel.tol=1E-12)$value +
	((digamma(1/4) + Euler) * pracma::psi(1, 1/4) - pracma::psi(2, 1/4)/2) / 16 +
	- 7/4 * pracma::zeta(3) + 5/32 * pi^3 + 3/16 * pi^2 * log(2) +
	+ Catalan * pi/4 + Catalan * log(2);

# TODO


### Pow = 2 (inside Log)

### I( log(x^2+y^2) / (x+y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x^2+y^2) / (x+y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- pi^2 / 24 + 3/2 * log(2)^2 - 4*log(2);


### I( log(x^2+y^2) / (1 + x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x^2+y^2) / (1+x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
-57/32 * pracma::zeta(3) + Catalan * pi/2;


### I( log(x^2+y^2) / (1 - x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x^2+y^2) / (1-x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
-55/32 * pracma::zeta(3) + Catalan * pi / 2;


### I( y * log(x^2+1) / (x^2+y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * log(x^2+1) / (x^2+y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi*log(2) - pi/2 + log(2)^2/2 - log(2);

### I( y * log(x^2+1) / (x^2*y^2 + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * log(x^2+1) / (x^2*y^2 + 1), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi*log(2) - log(2)^2 / 2 - 2*Catalan;


### I( x * log(x^2+y^2) / (x^2+y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * log(x^2+y^2) / (x^2+y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- (Catalan - pi*log(2)/2 + pi/2 - log(2)^2/4 + log(2));

### I( x^2 * log(x^2+y^2) / (x^2+y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * log(x^2+y^2) / (x^2+y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi/4 + log(2)/2 - 3/2;


### I( log(x^2-x*y+y^2) / (1 - x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x^2-x*y+y^2) / (1-x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) - log(x^2-x+1) * log(1 - x) / x, 0, 1)$value - pracma::zeta(3);
integrate(\(x) - log(1+x^3) * log(1 - x) / x, 0, 1)$value - 13/8 * pracma::zeta(3);
# TODO

# Helper:
id = 0:1000;
FUN3 = \(x) sapply(x, \(x) sum(x^(id+1) / (id+1)^2) * (2*x-1) / (x^2-x+1));
integrate(FUN3, 0, 1);
- sum((digamma((id+4)/6) - digamma((id+1)/6)) / (id+1)^2 ) / 2 +
+ sum((digamma((id+2)/2) - digamma((id+1)/2)) / (id+1)^2 ) / 2 + 2*pracma::zeta(3);
# TODO: closed form?

# Alternative: I( log(x^3+y^3) / .. ) - zeta(3) / 8;

# z = x/y => (not included: - 2*zeta(3))
FUNd = \(x, y) log(x^2-x+1) / (1-x*y^2) * y
integrate(\(x) sapply(x, \(y) integrate(FUNd, 0, 1/y, y=y, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-11)
#
FUNs = \(x, y) 1/2 * log(x^2-x+1) / (x-y) /x;
integrate(\(x) sapply(x, \(y) integrate(FUNs, sqrt(y), Inf, y=y, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-11)
#
FUNs = \(x, y) - log(x^2-x+1) / (x-y) /x;
integrate(\(x) sapply(x, \(y) integrate(FUNs, 1, sqrt(y), y=y, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-11)$value + pracma::zeta(3);
#
FUNs = \(x, y) 1/2 * log(y-sqrt(y)+1) / (sqrt(y)-x) /y;
integrate(\(x) sapply(x, \(y) integrate(FUNs, 0, y, y=y, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-11)$value + pracma::zeta(3);
#
integrate(\(x) -1/2 * log(x-sqrt(x)+1) * log(1 - sqrt(x)) / x,0, 1)$value + pracma::zeta(3);


# z = x+y on [0, 2] =>
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x^2-x*y+y^2) / (1-x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
FUN = \(x, y) 2*log(y^2 - 3*x*(y-x)) / (1-x*(y-x));
integrate(\(x) sapply(x, \(y) integrate(FUN, 0, y/2, y=y,
	rel.tol=1E-12)$value), 0, 2, rel.tol=1E-12)$value +
integrate(\(x) sapply(x, \(y) - integrate(FUN, 0, y-1, y=y,
	rel.tol=1E-12)$value), 1, 2, rel.tol=1E-12)$value;
#
integrate(\(x) sapply(x, \(y) 1/2 * integrate(FUN, 0, y, y=y,
	rel.tol=1E-12)$value), 0, 2, rel.tol=1E-12)$value +
integrate(\(x) sapply(x, \(y) integrate(FUN, y, 1, y=y,
	rel.tol=1E-12)$value), 1, 2, rel.tol=1E-12)$value;
#
integrate(\(x) sapply(x, \(y) 1/3 * integrate(FUN, 0, y, y=y,
	rel.tol=1E-12)$value), 0, 2, rel.tol=1E-12)$value +
integrate(\(x) sapply(x, \(y) 2/3 * integrate(FUN, y, 1, y=y,
	rel.tol=1E-12)$value), 0, 2, rel.tol=1E-12)$value - 4/27 * pracma::zeta(3);
#
integrate(\(x) sapply(x, \(y) 1/3 * integrate(FUN, 0, 1, y=y,
	rel.tol=1E-12)$value), 0, 2, rel.tol=1E-12)$value +
integrate(\(x) sapply(x, \(y) 1/3 * integrate(FUN, y, 1, y=y,
	rel.tol=1E-12)$value), 0, 2, rel.tol=1E-12)$value - 4/27 * pracma::zeta(3);



### I( log(x^2-x*y+y^2) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(x^2-x*y+y^2) / (1+x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	2/3 * log(x^2-x*y+y^2) / (1-x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)$value +
	- pracma::zeta(3) / 18;
# TODO: sub-integral;


### I( log(x^2+x*y+y^2) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(x^2+x*y+y^2) / (1-x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/3 * log(x^2-x*y+y^2) / (1-x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)$value +
	+ pracma::zeta(3) / 3;
# TODO: sub-integral;


### I( log(x^2+x*y+y^2) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(x^2+x*y+y^2) / (1+x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/2 * log(x^2-x*y+y^2) / (1-x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)$value +
	+ 2/9 * pracma::zeta(3);
# TODO: sub-integral;



### Prod( LOG )

### I( log(x) * log(y) * log(1 - x*y) )
# Maths 505: A mesmerizing result
# https://www.youtube.com/watch?v=QqVhd_xfjnc
# Note: series expansion of log(1 - x*y);


### Simple: I( log(x) * log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1-x*y) * log(x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
3 - sum(pracma::zeta(2:3));

### Simple: I( log(1-x) * log(1 - x*y) )
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

