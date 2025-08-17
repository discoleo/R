#####################
##
## Leonard Mada
## (the one and only)
##
## Double Integrals:
## Radicals of Polynomials
##
## v.0.2c


### Double Integrals:
# Radicals: Simple Entanglement

### Full Generalisation:
# I( x^p * y^q / (1 - x*y)^k )
# I( x^p * y^q * (1 - x*y)^k )
# & Various Reformulations;


### History:
# - Explicit Cases moved to this file;


#####################

### Gen: I( x^p * y^q / (1 - x^n*y^n)^(1/n) )
n = sqrt(19);
p = sqrt(3) - 1/5; q = sqrt(2) - 1/7;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p*y^q / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta((p+1)/n, 1-1/n) - beta((q+1)/n, 1-1/n)) / ((q-p)*n);


### Gen: I( x^p * y^q * (1 - x^n*y^n)^(1/n) )
p = sqrt(2); q = -1/sqrt(3); n = 1/sqrt(7);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p*y^q * (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta((p+1)/n, 1+1/n) - beta((q+1)/n, 1+1/n)) / ((q-p)*n);


### Various Reformulations:
# see file: Integrals.Double.Radicals.Beta.R;


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


##############

###########
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


##############

###########
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


##############

###########
### Pow = 6

### I( 1 / (1 - x^6*y^6)^(1/6) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 - x^6*y^6)^(1/6), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
sum(digamma(c(1:5, 7:11)/12) * rep(c(-1,1), each = 5)) *
	diff(digamma(c(1,6)/6)) / 6^2 / (5/2 + 2/sqrt(3));
n = 6;
- pi / sin(pi/n) * (digamma(1/n) + Euler) / n^2;


###########
### Pow = 7

### I( 1 / (1 - x^7*y^7)^(1/7) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 - x^7*y^7)^(1/7), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
sum(digamma(c(1:6, 8:13)/14) * rep(c(-1,1), each = 6)) *
	diff(digamma(c(1,7)/7)) / 7^2 / 4;
n = 7;
- pi / sin(pi/n) * (digamma(1/n) + Euler) / n^2;

