#####################
##
## Leonard Mada
## (the one and only)
##
## Double Integrals:
## Radicals of Polynomials: SQRT
##
## v.0.2a


### Double Integrals:
# SQRT of Unit-Polynomials


### Examples:
# I( 1 / sqrt( (1-x^2) * (1-y^2) * (1 - x^2*y^2) ) )
# Gen: I( sqrt( (1-x^n) * (1 - x^n*y) ) )
# Gen: I( sqrt( (1-x^n) / (1 - x^n*y) ) )
# Gen: I( y * sqrt( (1-x^n) * (1 - x^n*y) ) )
# Gen: I( y * sqrt( (1-x^n) / (1 - x^n*y) ) )


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;


####################

### SQRT: 2 Components

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
3 - 2*Catalan;


### I( sqrt( x/y * (1-x) / (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x/y * (1-x) / (1 - x*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-12)
8/9;


### SQRT: x & y Components

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


### I( sqrt( (1 - x^2*y^2) / (1-x^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1 - x^2*y^2) / (1-x^2) ), 0, 1)$value), 0, 1)
Catalan + 1/2;


### I( sqrt( (1-x^2) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^2) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
2*Catalan - 1;

### I( x^2 * sqrt( (1-x^2) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * sqrt( (1-x^2) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
2/9;

### I( x^2*y^2 * sqrt( (1-x^2) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2*y^2 * sqrt( (1-x^2) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan - 5/6;

### I( sqrt( 1/x * (1-x^2) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^2)/x / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# 2 * I( sqrt(x*y / ((1-x^2) * (1 - x^2*y^2))) ) +
	+ beta(3/4, 1/2) * beta(1/2, 1/2) / 2 - 4;
# Other:
integrate(\(x) 2*x / sin(x)^(3/2), 0, pi/2, rel.tol=1E-13)$value - 4;
# TODO


### I( sqrt( 1/y * (1-x^2) / (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^2)/y / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
gamma(1/4)^3 / gamma(3/4) * sqrt(2) / 24 - 2/3;


### I( sqrt((1-x^2)/(x*y) / (1 - x^2*y^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^2)/(x*y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-7)$value), 0, 1, rel.tol=1E-7)
# more robust:
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-y^2)/(x*y) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# I( sqrt(x) / sqrt((1-x^2) * (1 - x^2*y^2)) ) +
	+ beta(1/4, 1/2) * beta(1/2, 1/2) / 4 - 2;
# TODO


### Div: 2 Fractions

### I( 1 / sqrt( (1-x^2) * (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / sqrt((1-x^2) * (1 - x^2*y^2)), 0, 1)$value), 0, 1)
2*Catalan;

### I( x^2 / sqrt( (1-x^2) * (1 - x^2*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 / sqrt((1-x^2) * (1 - x^2*y^2)), 0, 1)$value), 0, 1)
1;

### I( sqrt( x / ((1-x^2) * (1 - x^2*y^2)) ) )
# - better accuray on [0,y] + [y,1];
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt(x / ((1-x^2) * (1 - x^2*y^2))), 0, 1, rel.tol=1E-9)$value), 0, 1, rel.tol=1E-9)
0.1943688252405 * 4 + 2/3;
# see I( x*y * ((1 - x^4) / (1 - x^4*y^4))^(3/4) );
0.4453575428195 * 8 - gamma(1/4)^2 * gamma(1/2) * sqrt(2) / 8 + 2;
# see I( x * ((1 - x^4) / (1 - x^4*y^4))^(1/4) );
# TODO


### I( sqrt( y / ((1-x^2) * (1 - x^2*y^2)) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt(y / ((1-x^2) * (1 - x^2*y^2))), 0, 1, rel.tol=3E-10)$value), 0, 1, rel.tol=1E-12)
2 - gamma(3/4)^4 / gamma(1/2)^2;
# alternative Formula:
2 - gamma(3/4)^3 / gamma(1/4) * sqrt(2);


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


### I( sqrt( (1-x)/(1+x) / (1 - x^2*y^2) ) )
# - are actually homogenous;
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



### Mixed: Pow = 1 & 2
# Special Case:
# - All non-trivial x are/become Pow = 1;
# - y is Pow = 2;

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
pracma::psi(1, 1/2) / 16; # may not be the generalisation;
pi^2 / 32;


### I( y^2 * sqrt( (1-x) / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y^2 * sqrt( (1-x) / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
3/2 - pi^2/8;


### I( sqrt( (1-y) / (1 - x*y^2) ) )
# - see also Inhomogeneous Radicals;
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


### Div: I( sqrt( x/y * (1-x) / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x/y * (1-x) / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(gamma(1/4)^3 / gamma(3/4) * sqrt(2) / 4 + 4) / 21;


### Div: I( sqrt( y/x * (1-x) / (1 - x*y^2) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( y/x * (1-x) / (1 - x*y^2) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=3E-9)
4 - 4*gamma(3/4)^4 / gamma(1/2)^2;
4 - beta(3/4, 3/4) * beta(3/4, 3/4);


### I( sqrt( x*(1 - x) * (1 - x*y) ) )
# Note: numerical issues with the x-variant;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( y*(1 - y) * (1 - x*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(digamma(7/4) - digamma(9/4)) / 3 + 4/9;

### I( sqrt( x*y*(1 - x) * (1 - x*y) ) )
# Note: numerical issues with the x-variant;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x*y*(1 - y) * (1 - x*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan - 7/10;


### I( sqrt( (1 - x)*(1 - y) * (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1 - x)*(1 - y) * (1 - x*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
15/8*Catalan - 21/16;

### I( sqrt( x*(1 - x)*(1 - y) * (1 - x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( y*(1 - x)*(1 - y) * (1 - x*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- pi^2 / 8 + 7/15*pi;

### I( sqrt( x*y * (1 - x)*(1 - y) * (1 - x*y) ) )
# Note: split into [0,y] + [y,1];
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x*y * (1 - x)*(1 - y) * (1 - x*y) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO

### I( sqrt( y/x * (1 - x)*(1 - y) * (1 - x*y) ) )
# Note: split into [0,y] + [y,1];
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( y/x * (1 - x)*(1 - y) * (1 - x*y) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
gamma(3/4)^3 / gamma(9/4) * sqrt(2) / 4;
beta(3/4,1/2)^2 / 10;
beta(3/4,1/2) * beta(3/4,3/2) / 4;


### SQRT: Higher Power

### Gen: I( sqrt( (1-x^n) * (1 - x^n*y) ) )
n = 7; # n = 8; # n = pi;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^n) * (1 - x^n*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- beta(1/n, 1/n) * 2^(2/n) / (6*(n-1)) + 4/3 * n^2 / (n^2-1);

### I( sqrt( (1 - x^3) * (1 - x^3*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^3) * (1 - x^3*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- beta(1/3, 1/3) * 2^(2/3) / 12 + 3/2;

### I( sqrt( (1-x^4) * (1 - x^4*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^4) * (1 - x^4*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- beta(1/4, 1/4) * sqrt(2) / 18 + 64/45;

### I( sqrt( (1-x^5) * (1 - x^5*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^5) * (1 - x^5*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- beta(1/5, 1/5) * 2^(2/5) / 24 + 25/18;


### I( sqrt( (1 - x^3) / (1 - x^3*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^3) / (1 - x^3*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- beta(1/3, 1/3) * 2^(2/3) / 4 + 3;

### Gen: I( sqrt( (1 - x^n) / (1 - x^n*y) ) )
n = 4; # n = 5; # n = 7; # n = pi;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( (1-x^n) / (1 - x^n*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- beta(1/n, 1/n) * 2^(2/n-1) / (n-1) + 2*n/(n-1);


### Free Terms:

### Gen: I( y * sqrt( (1-x^n) * (1 - x^n*y) ) )
n = 7; # n = pi; # n = exp(2);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * sqrt( (1-x^n) * (1 - x^n*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(- beta(1/n, 1/n) * 2^(2/n) * (n-2)*(n+1)/2 + 4*(4*n - 5)*n^2) / (15*(2*n-1)*(n^2-1));

### I( y * sqrt( (1-x^3) * (1 - x^3*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * sqrt( (1-x^3) * (1 - x^3*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(- beta(1/3, 1/3) * 2^(2/3) * 2 + 21*12) / 600;

### I( y * sqrt( (1-x^4) * (1 - x^4*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * sqrt( (1-x^4) * (1 - x^4*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(- beta(1/4, 1/4) * 2^(2/4) * 5 + 704) / (25*7*9);

### I( y * sqrt( (1-x^5) * (1 - x^5*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * sqrt( (1-x^5) * (1 - x^5*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(- beta(1/5, 1/5) * 2^(2/5) * 9 + 1500) / (3*1080);

### I( y * sqrt( (1-x^6) * (1 - x^6*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * sqrt( (1-x^6) * (1 - x^6*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(- beta(1/6, 1/6) * 2^(2/6) * 14 + 2736) / 5775;

### I( y * sqrt( (1-x^n) * (1 - x^n*y) ) )
n = 7;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * sqrt( (1-x^n) * (1 - x^n*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(- beta(1/n, 1/n) * 2^(2/n) * (n-2)*(n+1)/2 + 4*(4*n - 5)*n^2) / (15*(2*n-1)*(n^2-1));


### Gen: I( y * sqrt( (1-x^n) / (1 - x^n*y) ) )
n = 7; # n = pi;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * sqrt( (1-x^n) / (1 - x^n*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(- beta(1/n, 1/n) * 2^(2/n) * (n-2) + 8*n^2 - 12*n) / (6*(2*n-1)*(n-1));

### I( y * sqrt( (1-x^3) / (1 - x^3*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * sqrt( (1-x^3) / (1 - x^3*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(- beta(1/3, 1/3) * 2^(2/3) + 36) / 60;
(- beta(1/3, 1/6) + 36) / 60;

### I( y * sqrt( (1-x^4) / (1 - x^4*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * sqrt( (1-x^4) / (1 - x^4*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(- beta(1/4, 1/4) * 2^(2/4) * 2 + 80) / 126;
(- beta(1/4, 4/8) * 2 + 40) / 63;

### I( y * sqrt( (1-x^5) / (1 - x^5*y) ) )
n = 5;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * sqrt( (1-x^n) / (1 - x^n*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(- beta(1/5, 1/5) * 2^(2/5) * 3 + 140) / 216;

### I( y * sqrt( (1-x^n) / (1 - x^n*y) ) )
n = 6;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * sqrt( (1-x^n) / (1 - x^n*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(- beta(1/n, 1/n) * 2^(2/n) * (n-2) + 216) / 330;

#
n = 7;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * sqrt( (1-x^n) / (1 - x^n*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(- beta(1/n, 1/n) * 2^(2/n) * (n-2) + 308) / 468;

#
n = 8;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * sqrt( (1-x^n) / (1 - x^n*y) ), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(- beta(1/n, 1/n) * 2^(2/n) * (n-2) + 8*n^2 - 12*n) / (6*(2*n-1)*(n-1));


##############
### Triple ###
##############

### I( sqrt( x*y * (1 - x*y*z) ) )
FUN = \(x,y,z) sqrt( x*y * (1 - x*y*z) )
integrate(\(x) sapply(x, \(z) integrate(\(x, z1=z) sapply(x, \(y, z2=z1) integrate(FUN,
	y=y, z=z2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (pi*log(2) + 3/4 * pi - 16/3) / 2;

### I( x * sqrt( y*z * (1 - x*y*z) ) )
FUN = \(x,y,z) x * sqrt( y*z * (1 - x*y*z) )
integrate(\(x) sapply(x, \(z) integrate(\(x, z1=z) sapply(x, \(y, z2=z1) integrate(FUN,
	y=y, z=z2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(log(2) - 1/4) * pi / 2 - (digamma(7/4) - digamma(9/4)) / 2 - 2/3;

### I( sqrt( x*y*z * (1 - x*y*z) ) )
integrate(\(x) sapply(x, \(z) integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt( x*y*z * (1 - x*y*z) ), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
integrate(\(x) sapply(x, \(z) integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt(x * (1-x)) / (y*z), 0, y*z, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
integrate(\(x) 1/4 * sapply(x, \(z) integrate(\(x) sapply(x, \(y)
	# (asin(sqrt(y*z)) - sin(4*asin(sqrt(y*z)))/4) / (y*z)
	(asin(sqrt(y*z)) - sqrt(y*z*(1-y*z)) * (1 - 2*y*z)) / (y*z)
	), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
integrate(\(x) -1/4 * sapply(x, \(z) integrate(\(x) sapply(x, \(y)
	sqrt(y*z*(1-y*z)) * (1 - 2*y*z) / (y*z)
	), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)$value +
	+ (pi^2 / 12 + log(2)^2) * pi/4;
(pi^2 / 3 + 4*log(2)^2 - 2*log(2) - 5/2) * pi/16;

