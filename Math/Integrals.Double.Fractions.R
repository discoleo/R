#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
## Type: Fractions
##
## v.0.1a

### Double Integrals
### Type: Fractions

### Examples:
# I( 1 / (1 + x^3 + y^3) )



####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;
dzeta2  = - 0.937548254316;

#####################
#####################

### Basic: on [1, Inf]

### I( 1 / (x^n + y^n) )
n = 7.25; # sqrt(...) behaves very badly!
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (x^n + y^n), 1, Inf, rel.tol=1E-13)$value), 1, Inf, rel.tol=1E-13)
(digamma(1-1/(2*n)) - digamma(1/2-1/(2*n))) / (n*(n-2));


### I( 1 / (x^3 + y^3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (x^3 + y^3), 1, Inf, rel.tol=1E-13)$value), 1, Inf, rel.tol=1E-13)
(digamma(5/6) - digamma(1/3)) / 3;

### I( 1 / (x^4 + y^4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (x^4 + y^4), 1, Inf, rel.tol=1E-13)$value), 1, Inf, rel.tol=1E-13)
(digamma(7/8) - digamma(3/8)) / 8;

### I( x / (x^4 + y^4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / (x^4 + y^4), 1, Inf, rel.tol=1E-13)$value), 1, Inf, rel.tol=1E-13)
(digamma(7/8) - digamma(3/8) + pi) / 8;

### I( x / (x^5 + y^5) )
n = 5;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / (x^n + y^n), 1, Inf, rel.tol=1E-13)$value), 1, Inf, rel.tol=1E-13)
# TODO


### Simple: on [0, Inf]

### I( x^p * y^q / (1 + x^n + y^n)^k )
p = 0.35; q = 0.15;
n = 5; k = 2.125;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p * y^q / (1 + x^n + y^n)^k, 0, Inf, rel.tol=1E-13)$value), 0, Inf, rel.tol=1E-13)
# Numeric issues:
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p * y^q / (1 + x^n + y^n)^k, 0, Inf, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-10)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p * y^q / (1 + x^n + y^n)^k, 0, Inf, rel.tol=1E-8)$value), 1, Inf, rel.tol=1E-8)$value;
beta((p+1)/n, (q+1)/n) * beta((p+q+2)/n, k-(p+q+2)/n) / n^2;


### I( x^p / (1 + x^n + y^n) )
p = 0.35; n = 5;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p / (1 + x^n + y^n), 0, Inf, rel.tol=1E-13)$value), 0, Inf, rel.tol=1E-13)
# Numeric issues:
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p / (1 + x^n + y^n), 0, Inf, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-10)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p / (1 + x^n + y^n), 0, Inf, rel.tol=1E-8)$value), 1, Inf, rel.tol=1E-8)$value;
beta(1/n, (p+1)/n) * beta((p+2)/n, 1-(p+2)/n) / n^2;


### I( x^p * y^q / (1 + x^n + y^n) )
p = 0.35; q = 0.15; n = 5;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p * y^q / (1 + x^n + y^n), 0, Inf, rel.tol=1E-13)$value), 0, Inf, rel.tol=1E-13)
# Numeric issues:
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p * y^q / (1 + x^n + y^n), 0, Inf, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-10)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p * y^q / (1 + x^n + y^n), 0, Inf, rel.tol=1E-8)$value), 1, Inf, rel.tol=1E-8)$value;
beta((p+1)/n, (q+1)/n) * beta((p+q+2)/n, 1-(p+q+2)/n) / n^2;


### I( 1 / (1 + x^3 + y^3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 + x^3 + y^3), 0, Inf, rel.tol=1E-13)$value), 0, Inf, rel.tol=1E-13)
gamma(1/3)^3 / 9;


### I( 1 / (1 + x^4 + y^4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 + x^4 + y^4), 0, Inf, rel.tol=1E-13)$value), 0, Inf, rel.tol=1E-13)
gamma(1/4)^2 * gamma(1/2) / 16;


### I( 1 / (1 + x^n + y^n) )
n = 5;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 + x^n + y^n), 0, Inf, rel.tol=1E-10)$value), 0, Inf, rel.tol=1E-10)
beta(1/n, 1/n) * beta(2/n, 1-2/n) / n^2;

# Note: takes long;
# library(Rmpfr)
FUN = \(x, y, n) {
	x = mpfr(x, 160); y = mpfr(y, 160);
	as.numeric(1 / (1 + x^n + y^n));
}
n = sqrt(15); # sqrt() behaves very badly;
# upper = Inf; # done by splitting y over [0,1] & [1, 10^5];
integrate(\(x) sapply(x, \(y) integrate(FUN, 0, Inf, y=y, n=n,
	rel.tol=1E-7)$value), 0, 10^5, rel.tol=1E-7)
beta(1/n, 1/n) * beta(2/n, 1-2/n) / n^2;


### w. Coeff:

### I( x / (1 + x^n + y^n) )
n = 5.25; # behaves very badly w. most other values;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / (1 + x^n + y^n), 0, Inf, rel.tol=1E-13)$value), 0, Inf, rel.tol=1E-13)
beta(1/n, 2/n) * beta(3/n, 1-3/n) / n^2;

### I( x / (1 + x^4 + y^4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / (1 + x^4 + y^4), 0, Inf, rel.tol=1E-13)$value), 0, Inf, rel.tol=1E-13)
gamma(1/4)^2 * gamma(1/2) / 16;

### I( x / (1 + x^5 + y^5) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / (1 + x^5 + y^5), 0, Inf, rel.tol=1E-13)$value), 0, Inf, rel.tol=1E-13)
beta(1/5, 2/5) * beta(2/5, 3/5) / 25;


################

### on [0, 1]

### I( 1 / (1 + x^3 + y^3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 + x^3 + y^3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
gamma(1/3)^3 / 27;

### I( x*y / (1 + x^3 + y^3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y / (1 + x^3 + y^3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO

### I( x^3 / (1 + x^3 + y^3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^3 / (1 + x^3 + y^3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
1/2 - gamma(1/3)^3 / (2*27);


### Pow = 4

### I( 1 / (1 + x^4 + y^4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 + x^4 + y^4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### I( x / (1 + x^4 + y^4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x / (1 + x^4 + y^4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### I( x^2 / (1 + x^4 + y^4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 / (1 + x^4 + y^4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


###################
###################

### on [0, Inf]

### I( (1 - exp(-x^n-y^n)) / (x^n + y^n) )
n = 3.75; # n = 3.92; # n = 4.75;
# n = 5;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(1 - exp(-x^n-y^n)) / (x^n + y^n), 0, Inf, rel.tol=1E-6)$value), 0, Inf, rel.tol=1E-7)
gamma(1/n)^2 / (n*(n-2));

### I( (1 - exp(-x^3-y^3)) / (x^3 + y^3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(1 - exp(-x^3-y^3)) / (x^3 + y^3), 0, Inf, rel.tol=1E-11)$value), 0, Inf, rel.tol=1E-12)
gamma(1/3)^2 / 3;

### I( (1 - exp(-x^4-y^4)) / (x^4 + y^4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(1 - exp(-x^4-y^4)) / (x^4 + y^4), 0, Inf, rel.tol=1E-7)$value), 0, Inf, rel.tol=1E-7)
gamma(1/4)^2 / 8;

