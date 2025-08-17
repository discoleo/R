#####################
##
## Leonard Mada
## (the one and only)
##
## Double Integrals:
## Radicals of Polynomials
##
## v.0.2b


### Double Integrals:
# Radicals of Polynomials

### Examples:
# I( 1 / sqrt( (1-x^2) * (1-y^2) * (1 - x^2*y^2) ) )


### History:

# - [refactor] moved Diff-Type to file:
#   Integrals.Double.Radicals.Diff.R;

####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;


########################

### SQRT:

### I( sqrt( (1-x) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x) / (1 - x*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-12)
4*log(2) - 2;
2*(digamma(2) - digamma(3/2));


### I( sqrt( x*(1-x) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x*(1-x) / (1 - x*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-12)
digamma(7/4) - digamma(5/4);


### I( sqrt( x^3 * (1-x) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x^3*(1-x) / (1 - x*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-12)
pi/4 + (1/3 - 3)/5;
(digamma(3/4) - digamma(1/4)) / 4 + (1/3 - 3)/5;
(digamma(7/4) - digamma(5/4)) / 4 - 2*(2/3 - 1)/5;
(digamma(7/4) - digamma(9/4)) / 4 + 1/3;


### I( sqrt( (1-x)/x / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x)/x / (1 - x*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-12)
2*(digamma(7/4-1/2) - digamma(5/4-1/2));


### I( sqrt( y*(1-x) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( y*(1-x) / (1 - x*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-12)
- pi^2/4 + 3;

### I( sqrt( x*y * (1-x) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x*y * (1-x) / (1 - x*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-12)
4*Catalan - 3 - 1/3;


### I( sqrt( y/x * (1-x) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( y/x * (1-x) / (1 - x*y) ), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-8)
3 - 2*Catalan


### I( sqrt( x/y * (1-x) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x/y * (1-x) / (1 - x*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-12)
8/9;


### I( sqrt( (1-x)*(1-y) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x)*(1-y) / (1 - x*y) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
6*Catalan - 5;


### I( sqrt( (1-x)/(1-y) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x)/(1-y) / (1 - x*y) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(2*Catalan - 1) * 2;


###########
### Pow = 2

### I( 1 / sqrt( (1-x^2) * (1-y^2) * (1 - x^2*y^2) ) )
# - Series started as Trig-Radicals in:
#   Integrals.Double.Trig.R;

### Base: 1 / sqrt(...)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / sqrt((1-x^2) * (1-y^2) * (1 - x^2*y^2)), 0, 1)$value), 0, 1)
gamma(1/4)^3 / gamma(3/4) * sqrt(2) / 16;


### I( sqrt( (1 - x^2*y^2) * (1-x^2) / (1-y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1 - x^2*y^2) * (1-x^2) / (1-y^2) ), 0, 1)$value), 0, 1)
gamma(1/4)^3 / gamma(3/4) * sqrt(2) / 48;


### I( sqrt( (1-x^2) / ((1-y^2) * (1 - x^2*y^2)) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^2) / ((1-y^2) * (1 - x^2*y^2)) ), 0, 1)$value), 0, 1)
(gamma(1/4)^3 / gamma(3/4) / 16 - gamma(3/4)^3 / gamma(1/4)) * sqrt(2) / 2;


### Prod: I( sqrt( (1 - x^2*y^2) * (1-x^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1 - x^2*y^2) * (1-x^2) ), 0, 1)$value), 0, 1)
Catalan - 1/6;


### Fraction: I( 1 / sqrt( (1 - x^2*y^2) * (1-x^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / sqrt((1-x^2) * (1 - x^2*y^2)), 0, 1)$value), 0, 1)
2*Catalan;

### Fraction: I( x^2 / sqrt( (1 - x^2*y^2) * (1-x^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 / sqrt((1-x^2) * (1 - x^2*y^2)), 0, 1)$value), 0, 1)
1;


### I( sqrt( (1 - x^2*y^2) / (1-x^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1 - x^2*y^2) / (1-x^2) ), 0, 1)$value), 0, 1)
Catalan + 1/2;

### I( sqrt( (1-x^2) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^2) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
2*Catalan - 1


### I( sqrt( x / ((1-x^2) * (1 - x^2*y^2)) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt(x / ((1-x^2) * (1 - x^2*y^2))), 0, 1)$value), 0, 1)
# TODO


### I( sqrt( (1-x^2)/x / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^2)/x / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( sqrt( y / ((1-x^2) * (1 - x^2*y^2)) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt(y / ((1-x^2) * (1 - x^2*y^2))), 0, 1, rel.tol=3E-10)$value), 0, 1, rel.tol=1E-12)
2 - gamma(3/4)^4 / gamma(1/2)^2;
# alternative Formula:
2 - gamma(3/4)^3 / gamma(1/4) * sqrt(2);


### I( sqrt( (1-x^2)/y / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^2)/y / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
gamma(1/4)^3 / gamma(3/4) * sqrt(2) / 24 - 2/3;


### I( sqrt((1-x^2)/(x*y) / (1 - x^2*y^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^2)/(x*y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-7)$value), 0, 1, rel.tol=1E-7)
# I( x / sqrt((1-x^2) * (1 - x^2*y^2)) ) +
	+ beta(1/4, 1/2) * beta(1/2, 1/2) / 4 - 2;
# TODO


### I( 1 / sqrt(1 - x^2*y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / sqrt(1 - x^2*y^2), 0, 1)$value), 0, 1)
pi*log(2)/2;


### I( sqrt( x / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sqrt( x / (1 - x^2*y^2) ), 0, 1)$value), 0, 1)
pi - gamma(3/4)^2 / gamma(1/2) * 2 * sqrt(2);


### I( sqrt( x/y / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x/y / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
gamma(1/4)^2 / gamma(1/2) * sqrt(2) / 4 - gamma(3/4)^2 / gamma(1/2) * sqrt(2);
(beta(1/4, 1/4) - 2*beta(3/4, 3/4)) * sqrt(2) / 4;


### I( sqrt( x * (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sqrt( x * (1 - x^2*y^2) ), 0, 1)$value), 0, 1)
pi/2 - gamma(3/4)^2 / gamma(1/2) * sqrt(2) * 4/5;
pi/2 - gamma(3/4)*gamma(7/4) / gamma(7/2) * sqrt(2) * 2;


### Type: (1-x) * ...

### I( sqrt( (1-x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
#
integrate(\(x) 2 * sapply(x, \(y) integrate(\(x)
	sqrt( (1-x) * (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)$value +
	- (1/3 - 1/5 - 1/2^(3/2)/3 + 1/2^(5/2)/5) * sqrt(2) * 8;
# TODO


### I( sqrt( (1 - x^2*y^2) / (1-x) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1 - x^2*y^2) / (1-x) ), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-13)
# I( sqrt( (1-x) * (1 - x^2*y^2) ) ) +
	+ (1/3 - 1/5 - 1/2^(3/2)/3 + 1/2^(5/2)/5) * sqrt(2) * 16;


### I( sqrt( x*(1-x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x*(1-x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( sqrt( y*(1-x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( y*(1-x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( sqrt( (1-x)*(1-y) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x)*(1-y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( sqrt( (1-x) * (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x) * (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/2 * sqrt( (1-x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value
	), 0, 1, rel.tol=1E-13)$value +
	+ (1/3 - 1/5 - 1/2^(3/2)/3 + 1/2^(5/2)/5) * sqrt(2) * 4;
# TODO


### I( sqrt( (1-x)*(1-y) * (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x)*(1-y) * (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# I(sqrt( (1-x)*(1-y) / (1 - x^2*y^2) ) ) * 2/5 +
# I(sqrt( (1-x)/(1-y) * (1 - x^2*y^2) ) ) / 5;
# TODO


### I( sqrt( (1+x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1+x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	2 * sqrt( (1+x) * (1 - x^2*y^2) ), 0, 1)$value), 0, 1)$value - 14/15;
# TODO


### I( sqrt( (1-x)/(1+x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x)/(1+x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
2*Catalan - pi^2/8;

### I( sqrt( (1+x)/(1-x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1+x)/(1-x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-8)$value), 0, 1, rel.tol=1E-8)
2*Catalan + pi^2/8;

### I( x * sqrt( (1-x)/(1+x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * sqrt( (1-x)/(1+x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2/8 - 1;

### I( x * sqrt( (1+x)/(1-x) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * sqrt( (1+x)/(1-x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2/8 + 1;


### I( sqrt( (1-x)/(1+y) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x)/(1+y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO

### I( sqrt( (1+x)/(1+y) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1+x)/(1+y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### Mixed: Pow = 1 & 2

### I( x / sqrt((1-x^2) * (1-y^2) * (1 - x^2*y^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / sqrt((1-x^2) * (1-y^2) * (1 - x^2*y^2)), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/2 / sqrt((1-x) * (1-y^2) * (1 - x*y^2)), 0, 1)$value), 0, 1);
pi^2 / 4;


### I( x / sqrt( (1-y^2) * (1 - x^2*y^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / sqrt((1-y^2) * (1 - x^2*y^2)), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1/2 / sqrt((1-y^2) * (1 - x*y^2)), 0, 1)$value), 0, 1);
1;


### I( x / sqrt((1-x^2) * (1 - x^2*y^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / sqrt((1-x^2) * (1 - x^2*y^2)), 0, 1)$value), 0, 1)
pi^2/8;


### I( sqrt( x / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi/4;


### I( sqrt( y / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( y / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
gamma(3/4)^2 / gamma(1/2) * sqrt(2) * 4 - 4;


### I( sqrt( x*y / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x*y / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
gamma(3/4)^2 / gamma(1/2) * sqrt(2) * 4/3 - pi/3;
gamma(3/4)^2 / gamma(5/2) * sqrt(2) - pi/3;


### Div: I( sqrt( x/y / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x/y / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
gamma(1/4)^2 / gamma(1/2) * sqrt(2)/5 - pi/5;


### Div: I( sqrt( y/x / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( y/x / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
2*pi - gamma(3/4)^2 / gamma(1/2) * sqrt(2) * 4;


### I( sqrt( (1-x) / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x) / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2/8 - 1/2;


### I( x * sqrt( (1-x) / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * sqrt( (1-x) / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pracma::psi(1, 1/2) / 16; # pi^2 / 32;


### I( y^2 * sqrt( (1-x) / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y^2 * sqrt( (1-x) / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
3/2 - pi^2/8;


### I( sqrt( (1-y) / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-y) / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(sqrt(2) - 1)*6 + 2*log(sqrt(2) - 1);


### I( sqrt( x*(1-x) / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x*(1-x) / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
4/9;


### I( sqrt( y*(1-x) / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( y*(1-x) / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(2*pi - 4*log(2) - 2) / 3;


### I( sqrt( (1-x)*(1-y) / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x)*(1-y) / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### Div: I( sqrt( x/y * (1-x) / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x/y * (1-x) / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### Div: I( sqrt( y/x * (1-x) / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( y/x * (1-x) / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=3E-9)
4 - 4*gamma(3/4)^4 / gamma(1/2)^2;
4 - beta(3/4, 3/4) * beta(3/4, 3/4);


### Higher Order:

### I( 1 / sqrt( (1 - x*y^2) * (1 - x^2*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / sqrt( (1 - x*y^2) * (1 - x^2*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### I( sqrt( (1 - x*y^2) / (1 - x^2*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1 - x*y^2) / (1 - x^2*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


################
################

### Generalisations:

### Gen: I( x^p * y^q / (1 - x^n*y^n)^(1/n) )
n = sqrt(19);
p = sqrt(3) - 1/5; q = sqrt(2) - 1/7;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p*y^q / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta((p+1)/n, 1-1/n) - beta((q+1)/n, 1-1/n)) / ((q-p)*n);


### Gen Full: I( x^p * y^q / (1 - x^n*y^m)^k )
k = 1/sqrt(19);
m = sqrt(5); n = sqrt(7);
p = sqrt(3) - 1/5; q = sqrt(2) - 1/11;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p*y^q / (1 - x^m*y^n)^k, 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta((p+1)/m, 1-k) - beta((q+1)/n, 1-k)) / ((q+1)*m-(p+1)*n);


### Gen: I( x^p*y^q * (1 - x^n*y^n)^(1/n) )
p = sqrt(2); q = -1/sqrt(3); n = 1/sqrt(7);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p*y^q * (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta((p+1)/n, 1+1/n) - beta((q+1)/n, 1+1/n)) / ((q-p)*n);


### Gen: I( x^p / (1 - x^n*y^n)^(1/n) )
n = sqrt(13) + 1/3; p = sqrt(2);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi / sin(pi/n) - beta((p+1)/n, 1-1/n)) / (p*n);


### Gen: I( x^p*y / (1 - x^n*y^n)^(1/n) )
n = sqrt(19); p = sqrt(3) - 1/5;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p*y / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta(2/n, 1-1/n) - beta((p+1)/n, 1-1/n)) / ((p-1)*n);


### Gen: I( 1 / (1 - x^n*y^n)^(1/n) )
# Special Case: p = 0;
n = sqrt(29)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- pi / sin(pi/n) * (digamma(1/n) + Euler) / n^2;


### Gen: I( x / (1 - x^n*y^n)^(1/n) )
n = sqrt(5) + sqrt(6);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi / sin(pi/n) - beta(2/n, 1-1/n)) / n;


### Gen: I( x^2 / (1 - x^n*y^n)^(1/n) )
n = sqrt(13) + 1/3;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi / sin(pi/n) - beta(3/n, 1-1/n)) / (2*n);


### Gen: I( x^n / (1 - x^n*y^n)^(1/n) )
n = sqrt(3) + sqrt(11);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^n / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta(1/n, 1/n) * beta(2/n, 1 - 1/n) - beta(1/n, 1 - 1/n)) / n^3;


### Gen: I( x*y / (1 - x^n*y^n)^(1/n) )
n = sqrt(5) + 2/3;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
beta(2/n, 1-1/n) * (digamma(1+1/n) - digamma(2/n)) / n^2;

### Gen: I( x^2*y / (1 - x^n*y^n)^(1/n) )
n = sqrt(19);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2*y / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta(2/n, 1-1/n) - beta(3/n, 1-1/n)) / n;


### Derived:

### Gen: I( x^p / (1 - x^m*y^n)^(1/n) )
n = sqrt(13) + 1/3; m = sqrt(3) + sqrt(2); p = sqrt(2);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p / (1 - x^m*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi / sin(pi/n) - beta((p+1)/m, 1-1/n)) / ((p+1)*n-m);


### I( x^p * log(x) / (1 - x^n*y^n)^(1/n) )
n = sqrt(13) + 1/3; p = sqrt(2);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p * log(x) / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- beta((p+1)/n, 1-1/n) * (digamma((p+1)/n) - digamma(1+p/n)) / (p*n^2) +
	- (pi / sin(pi/n) - beta((p+1)/n, 1-1/n)) / (p^2*n);


### Series: I( ((1 - x^n) / (1 - x^n*y^n))^(1/n) )

### Gen: I( ((1 - x^n) / (1 - x^n*y^n))^(1/n) )
n = sqrt(5)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x^n) / (1 - x^n*y^n))^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pracma::psi(1, 1/(2*n)) - pracma::psi(1, 1/2 + 1/(2*n))) / (2*n^2) - 1;


### Gen: I( x^n / ((1-x^n) * (1 - x^n*y^n))^(1/n) )
n = sqrt(19);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^n / ((1-x^n) * (1 - x^n*y^n))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(1/2 - 1/(2*n), 1) * beta(1/n, 2-2/n) / n^2 / 2;


### Gen: I( 1 / ((1 - x^n) * (1 - y^n) * (1 - x^n*y^n))^(1/n) )
n = sqrt(3)*2 + 5/3;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1 - x^n) * (1 - y^n) * (1 - x^n*y^n))^(1/n), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-12)
beta(1/(2*n), 1-2/(2*n)) * beta(1/n, 1-3/(2*n)) / (2*n^2);


### Gen: I( ((1-x^n) / (1-y^n) * (1 - x^n*y^n))^(1/n) )
n = sqrt(31)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1-x^n) / (1-y^n) * (1 - x^n*y^n))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(1/(2*n), 1+1/n) * beta(1/(2*n), 1-1/n) / (2*n)^2;


### Type: (1-x) * ...

### Gen: I( ((1 - x) / (1 - x*y^n))^(1/n) )
n = sqrt(7)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x) / (1 - x*y^n))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pracma::psi(1, 1/n) / n^2 - 1/n;


### I( x * ((1 - x) / (1 - x*y^n))^(1/n) )
n = sqrt(11);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x) / (1 - x*y^n))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(pracma::psi(1, 1/n)*(1-1/n) - (n-2)/2) / (2*n^2);


### Gen: I( 1 / ((1-x) * (1-y^n) * (1 - x*y^n))^(1/n) )
n = sqrt(5) + sqrt(7)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1-x) * (1-y^n) * (1 - x*y^n))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(1/n, 2-3/n) * beta(1-1/n, 1-1/n) / n;


### Gen: I( ((1-x) / (1-y^n) * (1 - x*y^n))^(1/n) )
n = 2*sqrt(5) + 2/7;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1-x) / (1-y^n) * (1 - x*y^n))^(1/n), 0, 1, rel.tol=1E-9)$value), 0, 1, rel.tol=1E-9)
beta(1/n, 1-1/n) * beta(1/n, 2+1/n) / n^2;


###########

### Pow = 3

### I( 1 / (1 - x^3*y^3)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 - x^3*y^3)^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
sum(digamma(c(1,2,4,5)/6) * c(-1,-1,1,1)) * diff(digamma(c(1,3)/3)) / 18


### I( x / (1 - x^3*y^3)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / (1 - x^3*y^3)^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi / sin(pi/3) - beta(2/3, 2/3)) / 3;


### I( x^2 / (1 - x^3*y^3)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 / (1 - x^3*y^3)^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi / tan(pi/3) / 3 - 1/4;


### I( x^3 / (1 - x^3*y^3)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^3 / (1 - x^3*y^3)^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi / sin(pi/3) * 2/27;


### I( x*y / (1 - x^3*y^3)^(1/3) )
n = 3;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
beta(2/3, 2/3) / 3 - gamma(2/3)^3 / 6;

### I( x*y^2 / (1 - x^3*y^3)^(1/3) )
n = 3;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y^2 / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
beta(2/3, 2/3) / 3 - 1/2;
beta(2/3, 2/3) / ((2-1)*n) - beta(2/3, 3/3) / 3;

### I( x*y^3 / (1 - x^3*y^3)^(1/3) )
n = 3;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y^n / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
beta(2/3, 2/3) / 6 - pi / sin(pi/3) / 18;


### Gen: I( x*y^p / (1 - x^3*y^3)^(1/3) )
n = 3; # fixed;
p = sqrt(3);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y^p / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta(2/3, 2/3) - beta(2/3, (p+1)/n)) / ((p-1)*n);


### Mixed:

### I( ((1 - x^3) / (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x^3) / (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pracma::psi(1, 1/6) - pracma::psi(1, 4/6)) * 2/6^2 - 1;

# [old]
(pracma::psi(1, 1/3) - pracma::psi(1, 2/3) + pracma::psi(1, 4/3)) / 9;
(2*pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 9 - 1;


### I( ((1 - x^3) / (1 - x^3*y^3)^2)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x^3) / (1 - x^3*y^3)^2)^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(gamma(1/3)^3 - 2*pi * sqrt(3)) / 9;


### Series: 1 / ((1 - x^3) * (1 - x^3*y^3))^(1/3)

### I( 1 / ((1 - x^3) * (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1 - x^3) * (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-11)
# TODO


### I( x / ((1 - x^3) * (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / ((1 - x^3) * (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pi^2 * 2/27;


### I( x^3 / ((1 - x^3) * (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^3 / ((1 - x^3) * (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-11)
beta(1/3, 1/3) / 12;


### I( y^3 / ((1 - x^3) * (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y^3 / ((1 - x^3) * (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-11)
# TODO


### I( ((1-x^3) * (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1-x^3) * (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(gamma(1/3)^3 - pi / sin(pi/3)) / 18;


### I( ((1 - x^3) * (1 - x^3*y^3)^2)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x^3) * (1 - x^3*y^3)^2)^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(2*pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) * 2/27 - 2/3 + 1/4;
(pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 9 + pi^2 * 4/81 - 2/3 + 1/4;


### Full:

### I( 1 / ((1 - x^3) * (1 - y^3) * (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1 - x^3) * (1 - y^3) * (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-10)
gamma(1/6)^3 / gamma(1/2) * sqrt(3) * 2/6^3;
beta(1/6, 4/6) * beta(3/6, 2/6) * 2/6^2;


### Pow 1 & 3:

### I( 1 / ((1-x) * (1-y^3) * (1 - x*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1-x) * (1-y^3) * (1 - x*y^3))^(1/3), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
beta(2/3, 2/3);


### I( ((1 - x) / (1 - x*y^3)^2)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x) / (1 - x*y^3))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pracma::psi(1, 1/3) / 9 - 1/3;

### I( x * ((1 - x) / (1 - x*y^n))^(1/n) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x) / (1 - x*y^3))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pracma::psi(1, 1/3) / 3^3 - 1/3^2 + 1/12;


### I( ((1 - x) * (1 - x*y^3)^2)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x) * (1 - x*y^3)^2)^(1/3), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pracma::psi(1, 1/3) * 2/27 - 1/18;


### I( ((1 - x)^2 * (1 - x*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x)^2 * (1 - x*y^3))^(1/3), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pi * sin(pi/3)/6 + 1/8;
(digamma(5/6) - digamma(1/6)) / 12 + 1/8;


##############

### Pow = 4

### I( 1 / (1 - x^4*y^4)^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 - x^4*y^4)^(1/4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
sum(digamma(c(1:3, 5:7)/8) * rep(c(-1,1), each = 3)) *
	diff(digamma(c(1,4)/4)) / 4^2 / (2 + 1/sqrt(2));
- pi / sin(pi/4) * (digamma(1/4) + Euler) / 4^2;


### I( x / (1 - x^4*y^4)^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / (1 - x^4*y^4)^(1/4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi / sin(pi/4) - beta(2/4, 3/4)) / 4;


### I( x^2 / (1 - x^4*y^4)^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 / (1 - x^4*y^4)^(1/4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi/sin(pi/4) - beta(3/4, 3/4)) / 8;


### I( x^4 / (1 - x^4*y^4)^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^4 / (1 - x^4*y^4)^(1/4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta(1/4, 1/4) * beta(2/4, 3/4) - beta(1/4, 3/4)) / 64;


### I( x*y / (1 - x^4*y^4)^(1/4) )
n = 4;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
beta(3/n, 2/n) * (digamma(5/n) - digamma(2/n)) / n^2;

# Limit:
p = 1 + 1E-4;
(beta(2/n, 3/n) - beta(3/n, (p+1)/n)) / ((p-1)*n);


### I( x^2*y / (1 - x^4*y^4)^(1/4) )
n = 4;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2*y / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta(2/n, 1-1/n) - beta(3/n, 1-1/n)) / n;


### Mixed:

### I( ((1 - x^4) / (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x^4) / (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pracma::psi(1, 1/8) - pracma::psi(1, 5/8)) / 32 - 1;


### I( 1 / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( y / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# I(x^4*y / ((1 - x^4) * (1 - x^4*y^4))^(1/4)) * 2 + beta(1/4, 2/4) / 6;
# TODO


### I( x / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pracma::psi(1, 1/2) / 8;


### I( x*y / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pi/4 - gamma(3/4)^2 * gamma(1/2) * sqrt(2) / 8;


### I( x^2 / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(3/4, 3/4) * (pi/4 - sqrt(2)/2 - 1/2 * log(tan(pi/8))) / 2;


### I( x^2*y / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2*y / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(3/4, 3/4) * sqrt(2)/2 * (1/2 + sqrt(2)/4 * log(tan(pi/8)));
beta(3/4, 2/4) / 4 + beta(3/4, 3/4) / 4 * log(tan(pi/8));
# TODO: digamma?


### I( x^2*y^2 / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2*y^2 / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(3/4, 3/4) * (-pi/4 + sqrt(2)/2 - 1/2 * log(tan(pi/8))) / 4;


### I( x^3 / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^3 / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( x^3*y / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^3*y / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
gamma(3/4)^4 / gamma(1/2)^2 / 4;
beta(3/4, 3/4) * beta(3/4, 3/4) / 16;


### I( x^3*y^3 / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^3*y^3 / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(log(2)/4 - pi/8 + 1/2) / 3;


### I( x^4 / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^4 / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-11)
beta(1/4, 1/4) * sqrt(2) / 36;


### I( x^4*y / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^4*y / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
# I( - y / ((1 - x^4) * (1 - x^4*y^4))^(1/4)) / 2 + beta(1/4, 2/4) / 12;
# TODO

### I( x^4*y^4 / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^4*y^4 / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
# I( 1 / ((1 - x^4) * (1 - x^4*y^4))^(1/4) ) / 4 - beta(1/4, 2/4) / 24;
# TODO


### I( ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### I( ((1 - x^4) / (1 - x^4*y^4)^2)^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x^4) / (1 - x^4*y^4)^2)^(1/4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
gamma(1/4)^2 * gamma(1/2) * sqrt(2) / 16 - pi * sqrt(2) / 4;


### I( ((1 - x^4) * (1 - x^4*y^4)^2)^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x^4) * (1 - x^4*y^4)^2)^(1/4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
gamma(1/4)^2 * gamma(1/2) * sqrt(2) / 24 +
	+ beta(1/4, 3/4) / 16 - pi * sqrt(2) / 6;
# Note: last part = - pi * sqrt(2) * 5/48;


### I( ((1 - x^4) * (1 - x^4*y^4)^3)^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x^4) * (1 - x^4*y^4)^3)^(1/4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pracma::psi(1, 1/8) - pracma::psi(1, 5/8)) * 3/128 - 11/20;


### I( 1 / ((1 - x^4) * (1 - x^4*y^4)^3)^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1 - x^4) * (1 - x^4*y^4)^3)^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-11)
# TODO


### Full:

### I( 1 / ((1 - x^4) * (1 - y^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1 - x^4) * (1 - y^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-12)
beta(1/8, 6/8) * beta(5/8, 2/8) * 2 / 8^2;
# non-generalised:
beta(1/8, 5/8) * beta(3/8, 3/8) * 2 / 8^2;

# [old]
gamma(1/8)^3 / gamma(3/8) / 8^2 / (2*cos(pi/8) + sin(pi/8));
# Note:
(2*cos(pi/8) + sin(pi/8)) # ==
(2*sqrt(2) + 3) * sin(pi/8);


### I( ((1-x^4) / (1-y^4) * (1 - x^4*y^4))^(1/4) )
n = 4;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1-x^n) / (1-y^n) * (1 - x^n*y^n))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(1/(2*n), 1+1/n) * beta(1/(2*n), 1-1/n) / (2*n)^2;


### Pow 1 & 4:

### I( 1 / ((1-x) * (1-y^4) * (1 - x*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1-x) * (1-y^4) * (1 - x*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pi/2;


### I( ((1-x) / (1-y^4) * (1 - x*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1-x) / (1-y^4) * (1 - x*y^4))^(1/4), 0, 1, rel.tol=1E-9)$value), 0, 1, rel.tol=1E-9)
beta(1/4, 3/4) * beta(1/4, 9/4) / 4^2;


### I( x * ((1 - x) / (1 - x*y^4))^(1/4) )
n = 4;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x) / (1 - x*y^n))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(pracma::psi(1, 1/4) * 3/4 - 1) / 32;


### I( y * ((1 - x) / (1 - x*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * ((1 - x) / (1 - x*y^4))^(1/4), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-13)
(pi/2 + log(2) - 1)/3;


### I( ((1 - x)^3 * (1 - x*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x)^3 * (1 - x*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
(digamma(7/4) - digamma(1/4)) / 8;


##############

### Pow = 5

### I( 1 / (1 - x^5*y^5)^(1/5) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 - x^5*y^5)^(1/5), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
sum(digamma(c(1,2,3,4,6,7,8,9)/10) * c(-1,-1,-1,-1,1,1,1,1)) *
	diff(digamma(c(1,5)/5)) / 25 * cos(2*pi/5);
n = 5;
- pi / sin(pi/n) * (digamma(1/n) + Euler) / n^2;


### I( x / (1 - x^5*y^5)^(1/5) )
n = 5;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi / sin(pi/n) - beta(2/n, 1-1/n)) / n;


### I( x^2 / (1 - x^5*y^5)^(1/5) )
n = 5;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pi/sin(pi/n) - beta(3/n, 1-1/n)) / (2*n);


### I( x*y / (1 - x^5*y^5)^(1/5) )
n = 5
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
beta(1-1/n, 2/n) * (digamma(1+1/n) - digamma(2/n)) / n^2;


### Mixed

### I( 1 / ((1 - x^5) * (1 - y^5) * (1 - x^5*y^5))^(1/5) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1 - x^5) * (1 - y^5) * (1 - x^5*y^5))^(1/5), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-12)
gamma(1/10)^3 / gamma(3/10) / 10^2 / (tan(2*pi/5) - sin(pi/5));
beta(1/10, 8/10) * beta(7/10, 2/10) * 2 / 10^2;


### I( 1 / ((1-x) * (1-y^5) * (1 - x*y^5))^(1/5) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1-x) * (1-y^5) * (1 - x*y^5))^(1/5), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(1/5, 2/5) * beta(4/5, 4/5) * 2/15;
beta(1/5, 7/5) * beta(4/5, 4/5) / 5;


### I( ((1-x) / (1-y^5) * (1 - x*y^5))^(1/5) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1-x) / (1-y^5) * (1 - x*y^5))^(1/5), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(1/5, 4/5) * beta(1/5, 11/5) / 5^2;


### I( x * ((1 - x) / (1 - x*y^5))^(1/5) )
n = 5;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x) / (1 - x*y^n))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(pracma::psi(1, 1/n)*(1-1/n) - 3/2) / (2*n^2);
(pracma::psi(1, 1/n)*(1-1/n) - (n-2)/2) / (2*n^2);



##############

### Pow = 6

### I( 1 / (1 - x^6*y^6)^(1/6) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 - x^6*y^6)^(1/6), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
sum(digamma(c(1:5, 7:11)/12) * rep(c(-1,1), each = 5)) *
	diff(digamma(c(1,6)/6)) / 6^2 / (5/2 + 2/sqrt(3));


### I( ((1 - x^6) / (1 - x^6*y^6))^(1/6) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x^6) / (1 - x^6*y^6))^(1/6), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pracma::psi(1, 1/12) - pracma::psi(1, 7/12)) * 2/12^2 - 1;


### I( 1 / ((1-x) * (1-y^6) * (1 - x*y^6))^(1/6) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1-x) * (1-y^6) * (1 - x*y^6))^(1/6), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi * 3/8 * sqrt(3) * 2^(-2/3);
# alternative:
beta(1/6, 1/2) * beta(5/6, 5/6) / 8;
beta(1/6, 9/6) * beta(5/6, 5/6) / 6;


### I( ((1-x) / (1-y^6) * (1 - x*y^6))^(1/6) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1-x) / (1-y^6) * (1 - x*y^6))^(1/6), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(1/6, 5/6) * beta(1/6, 13/6) / 6^2;


##############

### Pow = 7

### I( 1 / (1 - x^7*y^7)^(1/7) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 - x^7*y^7)^(1/7), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
sum(digamma(c(1:6, 8:13)/14) * rep(c(-1,1), each = 6)) *
	diff(digamma(c(1,7)/7)) / 7^2 / 4;


### I( 1 / ((1 - x^7) * (1 - y^7) * (1 - x^7*y^7))^(1/7) )
n = 7; # fixed value;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1 - x^n) * (1 - y^n) * (1 - x^n*y^n))^(1/n), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-12)
gamma(1/(2*n))^3 / gamma(3/(2*n)) / (2*n)^2 / (2*sin(3*pi/7) + sin(2*pi/7));


### I( 1 / ((1-x) * (1-y^7) * (1 - x*y^7))^(1/7) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1-x) * (1-y^7) * (1 - x*y^7))^(1/7), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(1/7, 11/7) * beta(6/7, 6/7) / 7;

