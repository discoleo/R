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


### Gen: I( ((1 - x^n) / (1 - x^n*y^n))^(1/n) )
n = sqrt(5)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x^n) / (1 - x^n*y^n))^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pracma::psi(1, 1/(2*n)) - pracma::psi(1, 1/2 + 1/(2*n))) / (2*n^2) - 1;


### Gen: I( 1 / ((1-x) * (1-y^n) * (1 - x*y^n))^(1/n) )
n = sqrt(5) + sqrt(7)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1-x) * (1-y^n) * (1 - x*y^n))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(1/n, 2-3/n) * beta(1-1/n, 1-1/n) / n;


###########

### Pow = 3

### I( 1 / (1 - x^3*y^3)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 - x^3*y^3)^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
sum(digamma(c(1,2,4,5)/6) * c(-1,-1,1,1)) * diff(digamma(c(1,3)/3)) / 18


### I( ((1 - x^3) / (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x^3) / (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pracma::psi(1, 1/6) - pracma::psi(1, 4/6)) * 2/6^2 - 1;

# [old]
(pracma::psi(1, 1/3) - pracma::psi(1, 2/3) + pracma::psi(1, 4/3)) / 9;
(2*pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 9 - 1;


### I( x^3 / ((1 - x^3) * (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^3 / ((1 - x^3) * (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-11)
beta(1/3, 1/3) / 12;


### I( 1 / ((1 - x^3) * (1 - y^3) * (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1 - x^3) * (1 - y^3) * (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-10)
gamma(1/6)^3 / gamma(1/2) * sqrt(3) * 2/6^3;


### Pow 1 & 3:

### I( 1 / ((1-x) * (1-y^3) * (1 - x*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1-x) * (1-y^3) * (1 - x*y^3))^(1/3), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
beta(2/3, 2/3);


##############

### Pow = 4

### I( 1 / (1 - x^4*y^4)^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 - x^4*y^4)^(1/4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
sum(digamma(c(1:3, 5:7)/8) * rep(c(-1,1), each = 3)) *
	diff(digamma(c(1,4)/4)) / 4^2 / (2 + 1/sqrt(2));


### I( ((1 - x^4) / (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x^4) / (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(pracma::psi(1, 1/8) - pracma::psi(1, 5/8)) / 32 - 1;


### I( x^4 / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^4 / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-11)
beta(1/4, 1/4) * sqrt(2) / 36;


### I( 1 / ((1 - x^4) * (1 - y^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1 - x^4) * (1 - y^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-12)
gamma(1/8)^3 / gamma(3/8) / 8^2 / (2*cos(pi/8) + sin(pi/8));

# Note:
(2*cos(pi/8) + sin(pi/8)) # ==
(2*sqrt(2) + 3) * sin(pi/8);


### I( 1 / ((1-x) * (1-y^4) * (1 - x*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1-x) * (1-y^4) * (1 - x*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pi/2;


##############

### Pow = 5

### I( 1 / (1 - x^5*y^5)^(1/5) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 - x^5*y^5)^(1/5), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
sum(digamma(c(1,2,3,4,6,7,8,9)/10) * c(-1,-1,-1,-1,1,1,1,1)) *
	diff(digamma(c(1,5)/5)) / 25 * cos(2*pi/5);


### I( 1 / ((1 - x^5) * (1 - y^5) * (1 - x^5*y^5))^(1/5) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1 - x^5) * (1 - y^5) * (1 - x^5*y^5))^(1/5), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-12)
gamma(1/10)^3 / gamma(3/10) / 10^2 / (tan(2*pi/5) - sin(pi/5));


### I( 1 / ((1-x) * (1-y^5) * (1 - x*y^5))^(1/5) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / ((1-x) * (1-y^5) * (1 - x*y^5))^(1/5), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(1/5, 2/5) * beta(4/5, 4/5) * 2/15;
beta(1/5, 7/5) * beta(4/5, 4/5) / 5;


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

