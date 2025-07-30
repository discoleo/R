#####################
##
## Leonard Mada
## (the one and only)
##
## Double Integrals:
## Radicals of Polynomials
##
## v.0.1a


### Double Integrals:
# Radicals of Polynomials

### Examples:
# I( 1 / sqrt( (1-x^2) * (1-y^2) * (1 - x^2*y^2) ) )


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
	sqrt( (1-x) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


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
# TODO


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

### Pow = 3

### I( 1 / (1 - x^3*y^3)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 - x^3*y^3)^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
sum(digamma(c(1,2,4,5)/6) * c(-1,-1,1,1)) * diff(digamma(c(1,3)/3)) / 18


##############

### Pow = 5

### I( 1 / (1 - x^5*y^5)^(1/5) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 - x^5*y^5)^(1/5), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
sum(digamma(c(1,2,3,4,6,7,8,9)/10) * c(-1,-1,-1,-1,1,1,1,1)) *
	diff(digamma(c(1,5)/5)) / 25 * cos(2*pi/5);

