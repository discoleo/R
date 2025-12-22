
#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
## Type: ATAN( TRIG )
##
## v.0.1a


### Type: Atan( TRIG )

### Examples:
# I( atan(tan(x) * tan(y)) )
# I( atan(tan(x) / tan(y)) )
# I( atan(sqrt(tan(x)^2 + tan(y)^2)) )


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;
dzeta2  = - 0.937548254316;

#####################
#####################


### Atan( TRIG )

### I( atan(y * tan(x)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(y * tan(x)), 0, pi/2, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi^2 / 12;

### I( atan(tan(x) * tan(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(tan(x)*tan(y)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi^3/16;

### I( atan(tan(x) / tan(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(tan(x) / tan(y)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi^3/16;

### I( x * atan(tan(x)*tan(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * atan(tan(x)*tan(y)), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
integrate(\(x) - pi * atan(x) * log(1-x) / x, 0, 1)$value - pi^4 / 128;
# TODO


### I( atan(sqrt(tan(x) * tan(y))) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sqrt(tan(x) * tan(y))), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi^3/16;


### I( atan(sqrt(tan(x)^2 + tan(y)^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sqrt(tan(x)^2 + tan(y)^2)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
(digamma(7/8) - digamma(3/8)) * (digamma(5/8) - digamma(1/8)) * pi / 16;


### I( atan(sin(x) * sin(y)) / sin(x) )
# Maths 505: A beautiful double integral
# https://www.youtube.com/watch?v=tsxLAsiLyq4
# Note: Series expansion atan(x)
# => I(sin(x)^(2*k)) * I(sin(y)^(2*k+1)) => Beta;

integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
Catalan * pi / 2;


### I( atan(sqrt(sin(x) * sin(y))) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sqrt(sin(x) * sin(y))) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
Catalan * pi;


### I( atan(sin(x) * sin(y)) * cos(y) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) * cos(y) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi/2 * asinh(1) - (digamma(5/8) - digamma(3/8)) / 2;
# Note: asinh(1) == log(tan(3*pi/8));

### I( atan(sin(x) * sin(y)) * sin(y) / tan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) * sin(y) / tan(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi/2 * asinh(1) - (digamma(5/8) - digamma(3/8)) / 2;


### I( atan(sin(x) * sin(y)) * sin(y) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) * sin(y) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
(beta(1/4, 1/4) - 2*beta(3/4, 3/4)) * sqrt(2) * pi/16;


### I( atan(sin(x) * sin(y)) / (sin(x) * sin(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) / (sin(x) * sin(y)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-13)
# TODO

#
id = 0:100;
sum((-1)^id / (2*id+1) * beta(id+1/2, 1/2)^2) / 4;

# Note:
id = seq(0, 20000)
sum((-1)^id / (id+1) * beta(id+1/2, 1/2) * beta(id+3/2, 1/2)) / 4;
((beta(1/4, 1/4) + 2*beta(3/4, 3/4)) * sqrt(2) / 4 - pi) * pi/2;
#
sum((-1)^id / (id+1) / (2*id+1) * beta(id+1/2, 1/2) * beta(id+3/2, 1/2)) / 4
(pi - beta(3/4, 3/4) * sqrt(2)) * pi/2;


### I( atan(sin(x) * sin(y)) / tan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) / tan(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
# TODO


### I( atan(sin(x) * sin(y)) * cos(y) / tan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) * cos(y) / tan(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
Catalan - pi/4 + log(2)/2;


### I( atan(sin(x) * sin(y)) / tan(x) / tan(y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) / tan(x) / tan(y), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi^3 / 32;


### I( atan(sin(x) * sin(y)) / sin(x)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) / sin(x)^2 - sin(y)/x, 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
# TODO


### I( atan(sin(x) * sin(y)) * sin(y) / sin(x)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) * sin(y) / sin(x)^2 - sin(y)^2/x, 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi/4 * (1/2 - log(pi/2));


### Div: Sum

### I( atan(sin(x) * sin(y)) / sqrt(sin(x)^2 + sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) / sqrt(sin(x)^2 + sin(y)^2), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
beta(1/8, 3/4)^2 * beta(3/8, 3/4) *  beta(1/8, 1+1/4) * 2 / 8^4;
beta(1/8, 3/4)^2 * beta(1/4, 3/4) * 4/3 / 8^3;


### I( atan(sin(x) * sin(y)) * sin(x) / (sin(x)^2 + sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) * sin(x) / (sin(x)^2 + sin(y)^2), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
# ==
integrate(\(x) pi/8 * sapply(x, \(y) integrate(\(x)
	sqrt( x / (1-x^2) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( atan(sin(x) * sin(y)) * sin(x)/sin(y) / (sin(x)^2 + sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) * sin(x)/sin(y) / (sin(x)^2 + sin(y)^2), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
2.309240197354 / 2;
# see I( atan(sin(x) * sin(y)) / (sin(x) * sin(y)) )
# TODO


###########
### DIV ###

### I( atan(sin(x) / sin(y)) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) / sin(y)) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
integrate(\(x) pi/2 * atan(sqrt(sin(x))) / sin(x), 0, pi/2, rel.tol=1E-13)
# TODO


### I( atan(sin(x) / sin(y)) * sin(y) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) / sin(y)) * sin(y) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
(pi + 2*log(2)) * pi / 8;


### I( atan(sin(x) / sin(y)) * cos(y) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) / sin(y)) * cos(y) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
((digamma(5/8) - digamma(1/8))^2 - (digamma(7/8) - digamma(3/8))^2) / 16;
log(sqrt(2) + 1) * pi;

### I( atan(sin(x) / sin(y)) * sin(2*y) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) / sin(y)) * sin(2*y) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi/2 * asinh(1) + (digamma(5/8) - digamma(3/8)) / 2;


### I( atan(sin(x) / sin(y)) * sin(y)^2 / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) / sin(y)) * sin(y)^2 / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
3.668202541275/2 - (beta(1/4, 1/4) - 2*beta(3/4, 3/4)) * sqrt(2) * pi/32;
# TODO: I( atan(sin(x) / sin(y)) / sin(x) )


### I( atan(sin(x) / sin(y)) * sin(y) / tan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) / sin(y)) * sin(y) / tan(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi/2 * (asinh(1) - log(2) + 1) - (digamma(5/8) - digamma(3/8)) / 2;

# Note:
# also = pi/2 - 1/2 * I( atan(x/y) / sqrt(1+x*y) ) on [0,1]^2;


### I( atan(sin(x) / sin(y)) / sqrt(sin(x)^2 + sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) / sin(y)) / sqrt(sin(x)^2 + sin(y)^2), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
gamma(1/8)^2 / gamma(7/8)^2 * gamma(3/4)^3 * gamma(1/4) / 8^2 / 2;
beta(1/8, 3/4)^2 * beta(1/4, 3/4) / 8^2 / 2;


### Power = 2

### I( atan(sin(x)^2 / sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x)^2 / sin(y)^2), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
pi^3 / 16;

### I( atan(sin(x)^3 / sin(y)^3) )
p = 3; # any meaningful value;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x)^p / sin(y)^p), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
pi^3 / 16;


### I( atan(sin(x)^2 / sin(y)^2) * sin(y) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x)^2 / sin(y)^2) * sin(y) / sin(x), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
# see I( (x^2 - 1) * Im(log(x + 1i)) / (x^4 + 1) ) on [0,1];
(pracma::psi(1, 1/8) - pracma::psi(1, 5/8)) / 32 - Catalan / 2 +
- (digamma(7/8) - digamma(3/8)) * pi * 3/32 +
- (digamma(5/8) - digamma(1/8)) * pi / 32 +
- (digamma(7/8) - digamma(3/8)) * log(2) / 16 +
- (digamma(5/8) - digamma(1/8)) * log(2) / 16 +
+ sqrt(2)*pi/8 * (pi/2 + log(2)) + pi*log(2) / 8;


######################

### Prod

### I( atan(tan(x)*tan(y)) * atan(x/y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(tan(x)*tan(y)) * atan(x/y), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
pi^4 / 4^3;


######################
######################

### Varia

### I( atan(y^2/x^2) * exp(-x^2) * cos(y^2) / y )
# Hmath: The most difficult integral on the channel!
# https://www.youtube.com/watch?v=mV5nuLF9ocQ
# [in Russian]
# Note: Feynman's trick: dx on F(x, y); [with x instead of x^2]
# => Complex contour => I(dx(F(x, y))) = - pi/2 * Ei(-x);
# see also: I( exp(-x^2) * Ei(-x^2) ) on [0, Inf];

### Pow = 2
# Note: fails numerically;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(y^2/x^2) * exp(-x^2) * cos(y^2) / y, 0, Inf, rel.tol=1E-12)$value), 0, Inf, rel.tol=1E-12)
# upper = Inf;
FUN = \(x, y) atan(y/x^2) * exp(-x^2);
integrate(\(x) 1/2 * cos(x) / x * sapply(x, \(y) integrate(FUN,
	0, Inf, y=y, rel.tol=1E-12)$value), 0, pi*200, rel.tol=1E-12)$value +
integrate(\(x) 1/2 * cos(x) / x * sapply(x, \(y) integrate(FUN,
	0, Inf, y=y, rel.tol=1E-12)$value), pi*200, pi*900, rel.tol=1E-12, subdivisions = 2000)$value;
pi^(3/2) * log(1+sqrt(2)) / 4;


### Pow = 3
# Note: fails numerically;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	3 * atan(y^3/x^3) * exp(-x^3) * cos(y^3) / y, 0, Inf, rel.tol=1E-12)$value), 0, Inf, rel.tol=1E-12)
# upper = Inf;
FUN = \(x, y) atan(y/x^3) * exp(-x^3);
integrate(\(x) cos(x) / x * sapply(x, \(y) integrate(FUN,
	0, Inf, y=y, rel.tol=1E-12)$value), 0, pi*200, rel.tol=1E-12)$value +
integrate(\(x) cos(x) / x * sapply(x, \(y) integrate(FUN,
	0, Inf, y=y, rel.tol=1E-12)$value), pi*200, pi*900, rel.tol=1E-12, subdivisions = 2000)$value;
pi/2 * gamma(1/3) * (
	- (digamma(1/3) + Euler)/3 + 1/2*log((2^(2/3) + 2^(1/3) + 1)/3) +
	- 1/sqrt(3) * atan(1/sqrt(3) * (2^(1/3) - 1) / (2^(1/3) + 1)));

