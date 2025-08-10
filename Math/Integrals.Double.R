#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
##
## v.0.1l

### Double Integrals

### Examples:
# I( (x+y)*log(x+y) / (1+x*y) )
# I( x^p * y^q / (exp(x+y) - 1) )


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

### Type: Plus: I( log(1+x) * log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1+x) * log(1-x*y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pracma::zeta(3) * 5/8 - pi^2/6 + log(2)^2 - 4*log(2) + 3;


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

############
### ATAN ###

### I( atan((x+y)/2) )
# Note: (x+y) does not behave as nicely with ATAN;
integrate(\(x) sapply(x, \(y) integrate(\(x) atan((x+y)/2), 0, 1)$value), 0, 1)
3/2*pi + 2*log(5) - 6*log(2) - 3*atan(2)


### I( atan(x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y), 0, 1)$value), 0, 1)
- pi^2 / 48 + pi/4 - log(2)/2


### I( atan(x/y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y), 0, 1)$value), 0, 1)
pi/4;


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


### I( atan(x/y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (x+y), 0, 1)$value), 0, 1)
pi*log(2)/2

### I( atan(x/y) / (1+x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (1+x+y), 0, 1)$value), 0, 1)
(log(3) * 3/4 - log(2)) * pi

### I( atan(x/y) / (2 - x - y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (2 - x - y), 0, 1)$value), 0, 1)
pi*log(2)/2


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


### I( atan(x*y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (1 - x*y), 0, 1)$value), 0, 1)
integrate(\(x) - atan(x) * log(1-x) / x, 0, 1)$value +
	integrate(\(x) - log(1-x) * log(x) / (x^2+1), 0, 1)$value;
integrate(\(x) - atan(x) * log(x) / (1-x), 0, 1, rel.tol=1E-12); # alternative
#
integrate(\(x) -2 * log(1-x) * log(x) / (x^2+1), 0, 1)$value +
	+ pi^3 * 3/64 - Catalan * log(2)/2;
# TODO


### I( atan(x*y) * log(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * log(x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2 / 48 - pi/4 + log(2)/2 + 3/32 * pracma::zeta(3);


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


### I( atan(x/y) / sqrt(1-x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) / sqrt(1-x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi * (1-log(2));


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

### I( atan(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y) / (1+x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan((x+y)/2) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan((x+y)/2) / (1+x*y), 0, 1)$value), 0, 1)
# TODO


### Atan( Fraction )

### I( atan(x/y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) / (1 - x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi^3 / 24


##############
##############

###########
### Exp ###

### Fractions

### I( exp(x+y) / (1+x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) exp(x+y) / (1+x*y), 0, 1)$value), 0, 1)

### I( exp(x+y) / (1-x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) exp(x+y) / (1-x*y), 0, 1)$value), 0, 1)

### I( (exp(1) - exp((x+y)/2)) / (1-x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) (exp(1) - exp((x+y)/2)) / (1-x*y), 0, 1)$value), 0, 1)


### Exp( x*y )

### I( exp(x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) exp(x*y) / (x+y), 0, 1)$value), 0, 1)
# TODO


### Div: Exp

### I( (x+y) / (exp(x+y) - 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) (x+y) / (exp(x+y) - 1), 0, Inf)$value), 0, Inf)
2 * pracma::zeta(3);

### I( x*y / (exp(x+y) - 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) x*y / (exp(x+y) - 1), 0, Inf)$value), 0, Inf)
pracma::zeta(4); pi^4/90;


### Gen: I( (x+y) / (exp(k*(x+y)) - 1) )
k = sqrt(3);
integrate(\(x) sapply(x, \(y) integrate(\(x) (x+y) / (exp(k*(x+y)) - 1), 0, Inf)$value), 0, Inf)
2 * pracma::zeta(3) / k^3;

### Gen: I( x^p / (exp(x+y) - 1) )
p = sqrt(2);
integrate(\(x) sapply(x, \(y) integrate(\(x) x^p / (exp(x+y) - 1), 0, Inf)$value), 0, Inf)
gamma(p+1) * pracma::zeta(p+2);

### Gen: I( x^p * y^q / (exp(x+y) - 1) )
px = sqrt(2); py = sqrt(3);
integrate(\(x) sapply(x, \(y) integrate(\(x) x^px * y^py / (exp(x+y) - 1), 0, Inf)$value), 0, Inf)
gamma(px+1) * gamma(py+1) * pracma::zeta(px+py+2);

### Gen: I( x^p * y^q / (exp(k1*x + k2*y) - 1) )
px = sqrt(2); py = sqrt(3);
k1 = 1 / sqrt(3); k2 = sqrt(5) - sqrt(3);
integrate(\(x) sapply(x, \(y) integrate(\(x) x^px * y^py / (exp(k1*x+k2*y) - 1), 0, Inf)$value), 0, Inf)
gamma(px+1) * gamma(py+1) * pracma::zeta(px+py+2) / k1^(px+1) / k2^(py+1);


### I( (x+y) / (exp(x+y) + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) (x+y) / (exp(x+y) + 1), 0, Inf)$value), 0, Inf)
3/2 * pracma::zeta(3)

### Gen: I( (x+y) / (exp(k*(x+y)) + 1) )
k = sqrt(3);
integrate(\(x) sapply(x, \(y) integrate(\(x) (x+y) / (exp(k*(x+y)) + 1), 0, Inf)$value), 0, Inf)
3/2 * pracma::zeta(3) / k^3;

### Gen: I( x^p * y^q / (exp(k1*x + k2*y) + 1) )
px = sqrt(2); py = sqrt(3);
k1 = 1 / sqrt(3); k2 = sqrt(5) - sqrt(3);
integrate(\(x) sapply(x, \(y) integrate(\(x) x^px * y^py / (exp(k1*x+k2*y) + 1), 0, Inf)$value), 0, Inf)
gamma(px+1) * gamma(py+1) * pracma::zeta(px+py+2) / k1^(px+1) / k2^(py+1) * (1 - 1/2^(px+py+1));


### Other

### I( 1 / ((exp(x)+1) * (exp(y)+1) * (exp(x) + exp(y))) )
# Maths 505: This double integral will make you love calculus
# https://www.youtube.com/watch?v=-Nm5BofLiS4


### on [0, Inf]^2
integrate(\(x) sapply(x, \(y)
	integrate(\(x) 1 / ((exp(x)+1) * (exp(y)+1) * (exp(x) + exp(y))),
		0, Inf, rel.tol=1E-12)$value), 0, Inf, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y)
	integrate(\(x) 1 / ((x+1) * (y+1) * (1/x + 1/y)),
		0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12);
2*log(2) - log(2)^2 - pi^2/12;


### on [0,1]^2
integrate(\(x) sapply(x, \(y)
	integrate(\(x) 1 / ((exp(x)+1) * (exp(y)+1) * (exp(x)+exp(y))),
		0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y)
	integrate(\(x) 1 / ((x+1) * (y+1) * (1/x + 1/y)),
		exp(-1), 1, rel.tol=1E-12)$value), exp(-1), 1, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y)
	integrate(\(x) x / ((x+1) * (x+y)),
		exp(-1), 1, rel.tol=1E-12)$value), exp(-1), 1, rel.tol=1E-12)$value +
	- log(2/(1+exp(-1)))^2 / 2;
integrate(\(x) - x * log(x + exp(-1)) / (x+1), exp(-1), 1)$value +
2*log(2) - log(2)^2 + exp(-1) - 1 +
	- (1+exp(-1) - log(2))*log(1+exp(-1));
polylog2(-1 - 2/(exp(1)-1)) - polylog2(-2/(exp(1)-1)) +
	+ log(exp(2)-1) * log(2) - 2*log(exp(1)+1)*exp(-1) +
	- log(exp(1)-1) * log(exp(1)+1) +
	+ 2*log(2)*exp(-1) - log(2)^2 + 1;
# TODO: closed formula for Li2?


# Helper:
integrate(\(x) x * log(x + 1) / (x+1), exp(-1), 1)
2*log(2) - log(2)^2/2 - 1 - (exp(-1)+1)*log(exp(-1)+1) +
	+ exp(-1) + log(exp(-1)+1)^2/2

# Note:
# - for function polylog2, see file:
#   Integrals.Trig.Tan.R;
integrate(\(x) x * log(x + exp(-1)) / (x+1), exp(-1), 1)
polylog2(-2/(exp(1)-1)) - polylog2(-1 - 2/(exp(1)-1)) +
	+ 2*(1 - log(2))*exp(-1) - 1 + log(exp(1)+1)*exp(-1) +
	+ (log(exp(1)-1) - 1) * log((exp(1)+1)/2);


#####################

### Mixed: Log( Exp )

### I( log(exp(x) + exp(y)) )
# Maths 505: An extremely captivating double integral
# https://www.youtube.com/watch?v=etR1AI3anpU
# - for polylog2, see file: Integrals.Trig.Tan.R;

integrate(\(x) sapply(x, \(y) integrate(\(x) log(exp(x) + exp(y)), 0, 1)$value), 0, 1)
- 3/2*pracma::zeta(3) - pi^2/6 - 2*polylog2(-exp(1), 3) + 1/3;


### I( log(exp(1) - exp(x*y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(exp(1) - exp(x*y)), 0, 1)$value), 0, 1)
# TODO

###
integrate(\(x) sapply(x, \(y) integrate(\(x) log(exp(x*y) - 1), 0, 1)$value), 0, 1)
# TODO


####################

### Mixed: Log & Exp

### I( log(x) / (exp(x+y) - 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(x) / (exp(x+y) - 1), 0, Inf, rel.tol=1E-8)$value), 0, Inf, rel.tol=1E-6)
# gamma(1) * digamma(1) * zeta(2) + gamma(1) * dzeta(2);
- Euler * pi^2 / 6 + dzeta2;


### I( log(x+y) / (exp(x+y) - 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(x+y) / (exp(x+y) - 1), 0, Inf, rel.tol=1E-12)$value), 0, Inf, rel.tol=1E-11)
pi^2/6 - Euler*pi^2/6 + dzeta2;


### I( log(x+y) / (exp(x+y) + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(x+y) / (exp(x+y) + 1), 0, Inf, rel.tol=1E-12)$value), 0, Inf, rel.tol=1E-11)
(pi^2/6 * (log(2) - Euler + 1) + dzeta2) / 2;


### EXP(x*y)

### I( log(x+y) / (exp(x*y) - 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(x+y)*(1/(exp(x*y) - 1) - 1/(x*y)), 0, 1)$value), 0, 1)
# TODO


### I( log(x+y) / (exp(x*y) + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(x+y) / (exp(x*y) + 1), 0, 1)$value), 0, 1)
# TODO

### I( log(x+y+1) / (exp(x*y) + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(x+y+1) / (exp(x*y) + 1), 0, 1)$value), 0, 1)
# TODO

### I( log(1 - x*y) / (exp(x*y) - 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1 - x*y) / (exp(x*y) - 1), 0, 1)$value), 0, 1)
# TODO

### I( log(1 - x*y) / (exp(x*y) + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1 - x*y) / (exp(x*y) + 1), 0, 1)$value), 0, 1)
# TODO


##############

### Other Exp/Base:

### I( log(2^x + 2^y + 1) / (2^(x+y) * (2^x + 2^y)) )
# Dr Peyam: A surprisingly elegant double integral
# https://www.youtube.com/watch?v=q17ezsqKRis
# Note: => triple integral:
# I( 1 / 2^(x+y) / (z*(2^x + 2^y) + 1) ) with z on [0,1] =>
# I( x*y / (x*y + x*z + y*z) ) / log(2)^2 on [0,1]^3;


# Note: numeric instability: upper = Inf;
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(2^x + 2^y + 1) / (2^(x+y) * (2^x + 2^y)),
		0, 100, rel.tol=1E-10)$value), 0, 100, rel.tol=1E-8)
1/3 / log(2)^2


# library(Rmpfr)
integrate(\(x) sapply(x, \(y)
	integrate(\(x) {
		x = mpfr(x, 240); y = mpfr(y, 240);
		v = log(2^x + 2^y + 1) / (2^(x+y) * (2^x + 2^y));
		as.numeric(v);
		}, 0, Inf, rel.tol=1E-10)$value), 0, Inf, rel.tol=1E-10)

