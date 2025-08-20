#####################
##
## Leonard Mada
## (the one and only)
##
## Double Integrals:
## Radicals of Polynomials
## w. 2 Components
##
## v.0.2d


### Double Integrals:
# Radicals of Polynomials
#   with 2 Components

# Note:
# - Explicit cases moved to this file;


### Series: 2 Fractions

### Gen: I( x*y^p / ((1 - x^n) * (1 - x^n*y^n))^(1/n) )
n = sqrt(7) + 1/5; p = sqrt(3);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y^p / ((1 - x^n) * (1 - x^n*y^n))^(1/n), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-12)
(pi / sin(2*pi/n) - beta((p+1)/n, 1-2/n) * beta(2/n, 1-1/n) / n) / (n*p);


##################
##################

### Explicit Cases


### Power = 3

### I( x / ((1 - x^3) * (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / ((1 - x^3) * (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pi^2 * 2/27;


### I( x*y / ((1 - x^3) * (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y / ((1 - x^3) * (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
(pi / sin(pi/3) - gamma(2/3)^3) / 3;
(pi / sin(pi/3) - beta(2/3, 2/3) * beta(2/3, 1/3) / 3) / 3;


### I( x*y^2 / ((1 - x^3) * (1 - x^3*y^3))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y^2 / ((1 - x^3) * (1 - x^3*y^3))^(1/3), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
(pi / sin(pi/3) - beta(2/3, 2/3) * beta(3/3, 1/3) / 3) / 6;


##############

### Pow = 4

### I( x / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pracma::psi(1, 1/2) / 8; # pi^2 / 16;
beta(1/2, 1/2) * (digamma(3/4) - digamma(1/4)) / 16;


### I( x*y^(1/2) / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*sqrt(y) / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
(pi - beta(3/8, 1/2) * gamma(3/4)^2 / gamma(1/2) * sqrt(2)/2) / 4 * 2;
(pi - beta(3/8, 1/2) * beta(3/4, 1/2) / 4) / 4 * 2;


### I( x*y^(2/3) / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y^(2/3) / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-13)
(pi - beta(2/12, 6/12) * beta(7/12, 7/12) * beta(1/4, 1/4) / beta(1/3, 1/3) / 12) / (4*2/3);
(pi - beta(5/12, 1/2) * beta(3/4, 1/2) / 4) / (4*2/3);


### I( x*y / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pi/4 - gamma(3/4)^2 * gamma(1/2) * sqrt(2) / 8;
(pi - beta(3/4, 3/4) * beta(1/4, 3/4) / 4) / 4;


### I( x*y^(4/3) / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y^(4/3) / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
(pi - beta(7/12, 1/2) * beta(3/4, 1/2) / 4) / (4*4/3);


### I( x*y^1.5 / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y^1.5 / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
(pi - beta(6/8, 3/8) * beta(7/8, 2/8) / 8) / (4*3/2);


### I( x*y^2 / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y^2 / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pi/8 - gamma(3/4)^4 / gamma(3/2)^2 / 16;
(pi + beta(3/4, 3/4) * gamma(3/4)*gamma(-1/4) / gamma(1/2) / 4) / 8;


### I( x*y^3 / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y^3 / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
(pi - gamma(3/4)^2 / gamma(1/2) * sqrt(2)) / 12;
(pi + gamma(7/4)*gamma(-1/2) / gamma(5/4) * beta(3/4, 1) / 4) / 12;


### I( x*y^4 / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y^4 / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
(pi - pi/3) / 16;
(pi + beta(3/4, 3/4) * gamma(5/4) * gamma(-3/4) / gamma(1/2) / 4) / 16;
pi/24;


### I( x*y^5 / ((1 - x^4) * (1 - x^4*y^4))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y^5 / ((1 - x^4) * (1 - x^4*y^4))^(1/4), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
(pi - gamma(3/4)^2 * gamma(1/2) * sqrt(2) / 4) / 20;
(pi - beta(3/4, 3/4) * beta(5/4, 3/4) / 2) / 20;

