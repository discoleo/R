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
### Type: Log

### Examples:
# I( 1 / (1 + x^3 + y^3) )



####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;
dzeta2  = - 0.937548254316;

#####################
#####################

### Simple: on [0, Inf]

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
n = sqrt(15);
# upper = Inf; # done by splitting y over [0,1] & [1, 10^5];
integrate(\(x) sapply(x, \(y) integrate(FUN, 0, Inf, y=y, n=n,
	rel.tol=1E-7)$value), 0, 10^5, rel.tol=1E-7)
beta(1/n, 1/n) * beta(2/n, 1-2/n) / n^2;


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


### I( 1 / (1 + x^4 + y^4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	1 / (1 + x^4 + y^4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### I( 1 / (1 + x^4 + y^4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 / (1 + x^4 + y^4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)

