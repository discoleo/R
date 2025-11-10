#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
## Type: Log
##
## v.0.2a

### Double Integrals

### Examples:
# I( x * log(x+y) / (1+x*y) )
# I( log(x^2+y^2) / (1 - x*y) )
# I( log(x^2+y^2) / (1 + x*y) )


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;
dzeta2  = - 0.937548254316;

#####################
#####################


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


####################
### Fractions of Log

### I( log(1 + x + y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1+x+y) / (x+y), 0, 1)$value), 0, 1)
# Note:
# - polylog2 is in file Integrals.Polylog.Helper.R;
integrate(\(x) sapply(x, \(x) - Re(polylog2(-x-1))), 0, 1)$value +
	- pi^2/12 + 2*log(2) - 1;
- 2*polylog2(-2) - pi^2/6 - 3*log(3) + 4*log(2);


### I( log(1-x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1-x*y) / (x+y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
log(2)^2 - 4*log(2) - pi^2/6 + 4*Catalan;

### I( log(1+x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1+x*y) / (x+y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi^2 / 4 + log(2)^2 - 4*log(2);

### I( log(1 + x + y - x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1 + x + y - x*y) / (x+y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi^2 / 12 + log(2)^2 * 3/2 - 4*log(2) + 2*Catalan;

### I( log(1 + x - x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1 + x - x*y) / (x+y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
integrate(\(x) -18 * (log(1+x)/(x^3 - 1) - log(2)/3/(x-1)), 0, 1)$value +
	+ 3*pi^2/6 - 3*log(2)^2 - 4*log(2);
# TODO


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

