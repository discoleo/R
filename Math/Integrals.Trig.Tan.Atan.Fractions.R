###################
##
## Trig: ATAN
## Polynomial Fractions of Atan



####################

### Helper Constants
Catalan = 0.915965594177219015054603514;


####################
####################

### Simple Fractions

### I( atan(x) / (x+1) )
integrate(\(x) atan(x) / (x+1), 0, 1)
pi * log(2)/8

###
integrate(\(x) atan(x) / (x+1)^2, 0, 1)
log(2)/4

###
integrate(\(x) atan(x) / (x+1)^3, 0, 1)
- pi/32 + log(2)/8 + 1/8


### on [0, Inf]
integrate(\(x) atan(1/x) / (x+1), 0, Inf)
pi*log(2)/4 + Catalan


### I( (atan(x) - pi/4) / (x - 1) )
integrate(\(x) (atan(x) - pi/4) / (x - 1), 0, 1)
- pi*log(2)/8 + Catalan

###
integrate(\(x) (atan(x) - pi/4) / (x - 1)^2 + 1/(2*x), 0, 1)
integrate(\(x) 1 / ((x^2 + 1)*(x - 1)) + 1/(2*x) + pi/4 - 1/2, 0, 1)
-log(2)/4 + pi/8 - 1/2


### Other

### I( (atan(x) - pi/4) / (x^2 - 1) )
integrate(\(x) (atan(x) - pi/4) / (x^2 - 1), 0, 1)
Catalan / 2

# Note: Extra-Component: pi/8*log(2);


####################

### I( x * atan(1-x) / (x^2+1) )
integrate(\(x) x*atan(1-x) / (x^2+1), 0, 1)
atan(1/2)*log(2) / 2 - 0.5984967901246/2 + log(5)*(pi/4 + atan(1/3) - atan(1/2)) / 4
# where Im(Li2((3+1i)/5)) = 0.5984967901246/2;
# TODO: solve Li2;


### Composite Fractions

### I( atan(x) / (x*(x^2+1)) )
integrate(\(x) atan(x) / (x*(x^2+1)), 0, Inf)
pi*log(2)/2

### on [0, 1]
integrate(\(x) atan(x) / (x*(x^2+1)), 0, 1)
pi*log(2)/8 + Catalan/2


### I( atan(x) / ((x+1) * (x^2+1)) )
integrate(\(x) atan(x) / ((x+1) * (x^2+1)), 0, Inf)
1/8 * pi*log(2) + pi^2/16 - Catalan/2

### on [0, 1]
integrate(\(x) atan(x) / ((x+1) * (x^2+1)), 0, 1)
pi*log(2)/8 + pi^2/64 - 1/4 * Catalan

# Note: see next section;


### I( atan(x) / (x*(x+1)*(x^2+1)) )
# Maths 505: A MONSTER INTEGRAL!!!
# int 0 to infty arctan(x)/(x(x+1)(x^2+1))
# https://www.youtube.com/watch?v=u-FKjn_83l8

### on [0, Inf]
integrate(\(x) atan(x) / (x*(x+1)*(x^2+1)), 0, Inf)
integrate(\(x) x / (tan(x) * (tan(x)+1)), 0, pi/2)
integrate(\(x) x / tan(x) - x / (tan(x)+1), 0, pi/2)
3/8 * pi*log(2) - pi^2/16 + Catalan/2

### on [0, 1]
integrate(\(x) atan(x) / (x*(x+1)*(x^2+1)), 0, 1)
- pi^2/64 + 3/4 * Catalan


# Helper: I( x * cos(x) / (cos(x) + sin(x)) )
x = pi/5
integrate(\(x) x * cos(x) / (cos(x) + sin(x)), 0, x)
x*log(cos(x) + sin(x)) / 2 + x^2/4 +
	- integrate(\(x) 1/2 * log(cos(x) + sin(x)), 0, x)$value
x*log(cos(x) + sin(x)) / 2 + x^2/4 - x*log(2)/4 +
	- integrate(\(x) 1/2 * log(sin(x)), pi/4, x+pi/4)$value

# on [0, pi/4]
integrate(\(x) x * cos(x) / (cos(x) + sin(x)), 0, pi/4)
pi^2/64 + pi*log(2)/8 - Catalan/4

#
x = pi/5
integrate(\(x) cos(x) / (cos(x) + sin(x)), 0, x)
(log(cos(x) + sin(x)) + x) / 2;


#####################
#####################

### Fractions: Higher Power

### I( atan(x^2) / (x^2 + 1) )
# Maths 505: An outrageous journey of integration: int 0 to pi/4 arctan(cot^2(x))
# https://www.youtube.com/watch?v=VUGlU_dSgPY
# Initial: I( atan(tan(x)^2) )

### I( atan(x^2) / (x^2 + 1) )
integrate(\(x) atan(tan(x)^2), 0, pi/4)
integrate(\(x) atan(x^2) / (x^2 + 1), 0, 1)
pi^2/16 - (digamma(3/8 + 1/2) - digamma(3/8)) *
	(digamma(1/8 + 1/2) - digamma(1/8)) / 32;


### Pow = 4

### I( x * atan(x) / (x^4 + 1) )
integrate(\(x) x * atan(x) / (x^4 + 1), 0, 1)
(digamma(7/8) - digamma(3/8)) * (digamma(5/8) - digamma(1/8)) / 64;

# Note:
# - uses double integral: on [0,1] x [0,1]
#   I(I( u^2 / ((u^2*t^2 + 1) * (u^4 + 1)) ));


### I( x^3 * atan(x) / (x^4 + 1) )
integrate(\(x) x^3 * atan(x) / (x^4 + 1), 0, 1)
(digamma(3/8 + 1/2) - digamma(3/8))^2 / 128 +
	- (digamma(1/8 + 1/2) - digamma(1/8))^2 / 128 +
	+ Catalan/2;

### on [0, Inf]
integrate(\(x) x^3 * atan(x) / (x^4 + 1) - pi/2/(x+1), 0, Inf)
- pi * log(2*cos(pi/8)) / 2

### on [0, Inf]
integrate(\(x) x * atan(x) / (x^4 + 1), 0, Inf)
integrate(\(x) x * (pi/2 - atan(x)) / (x^4 + 1), 0, Inf)
pi^2/16

### I( atan(x) / (x^4 + 1) )
integrate(\(x) atan(x) / (x^4 + 1), 0, Inf)
(pracma::psi(1, 7/8) - pracma::psi(1, 3/8)) / 64 +
	+ pi/32 * (digamma(5/8) - digamma(1/8));

### I( x^2 * atan(x) / (x^4 + 1) )
integrate(\(x) x^2 * atan(x) / (x^4 + 1), 0, Inf)
integrate(\(x) (pi/2 - atan(x)) / (x^4 + 1), 0, Inf)
(pracma::psi(1, 3/8) - pracma::psi(1, 7/8)) / 64 +
	+ pi/32 * (digamma(1/8) - digamma(5/8)) + pi^2/8 / sin(pi/4);


### Pow: 6

###
# x => 1/x;
integrate(\(x) x^2 * atan(x) / (x^6 + 1), 0, Inf)
integrate(\(x) x^2 * (pi/2 - atan(x)) / (x^6 + 1), 0, Inf)
pi^2 / 24;

###
integrate(\(x) (x^4 + 1) * atan(x) / (x^6 + 1), 0, Inf)
pi^2 / 6

###
integrate(\(x) (x^3 + x) * atan(x) / (x^6 + 1), 0, Inf)
pi^2 / 12 / sin(2*pi/3)
pi^2 / (6*sqrt(3))


# Helper:
integrate(\(x) atan(x) / x, 0, 1)
Catalan

# I( atan(x^n) / x )
n = sqrt(3);
integrate(\(x) atan(x^n) / x, 0, 1)
Catalan / n;

### Fraction Decomposition:
# u^4 / ((u^2*t^2 + 1) * (u^4 + 1))
# 1 / (u^2*t^2 + 1) - Fr0;
# =>
# Fr0 + (u^2*t^2 - 1) / ((u^4 + 1) * (t^4 + 1));
# where Fr0 = 1 / ((u^2*t^2 + 1) * (u^4 + 1));
# and (u^2*t^2 - 1) / Prod() is decomposable into 2 separate integrals:
# I( u^2 / (u^4 + 1) )^2 -  I( 1 / (u^4 + 1) )^2;


### Note: Lim p -> 0
p = 1E-4
integrate(\(x) atan(x^(1/p)) / p, 0, 1)
pi/(4*p) - integrate(\(x) 1/p^2 * x^(1/p) / (x^(2/p) + 1), 0, 1)$value
pi/(4*p) - integrate(\(x) 1/p * x^p / (x^2 + 1), 0, 1)$value
- (pracma::psi(1, 3/4) - pracma::psi(1, 1/4)) / 16
Catalan

