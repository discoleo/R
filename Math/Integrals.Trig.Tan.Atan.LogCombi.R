###################
##
## Trig: ATAN
## Combinations with Log



####################

### Helper Constants
Catalan = 0.915965594177219015054603514;


####################
####################

### Prod: ATAN * LOG


### I( x^p * log(x) * atan(x) )
# Maths 505: A cool result: int (0,1) arctan(x)ln(x)
# https://www.youtube.com/watch?v=zL169U85SZc
# Note: uses series expansion of atan,
#   but method based on I( LOG(x) / (x^2+1) ) is more powerful;


###
p = sqrt(3)
integrate(\(x) x^p * log(x) * atan(x), 0, 1)
(pracma::psi(1, (p+2)/4) - pracma::psi(1, p/4 + 1)) / (16*(p+1)) +
	- (digamma((p+2)/4) - digamma(p/4 + 1) + pi) / (4*(p+1)^2);


### I( x^p * log(x)^2 * atan(x) )
p = sqrt(3)
integrate(\(x) x^p * log(x)^2 * atan(x), 0, 1)
(pracma::psi(2, (p+2)/4) - pracma::psi(2, p/4 + 1)) / (64*(p+1)) +
	- (pracma::psi(1, (p+2)/4) - pracma::psi(1, p/4 + 1)) / (8*(p+1)^2) +
	+ (digamma((p+2)/4) - digamma(p/4 + 1) + pi) / (2*(p+1)^3);


### Example:
integrate(\(x) log(x) * atan(x), 0, 1)
pi^2 / 48 - pi/4 + log(2)/2


### Varia: Base Integral
p = sqrt(3)
integrate(\(x) x^p * log(x) / (x^2 + 1), 0, 1)
(pracma::psi(1, (p+1)/4 + 1/2) - pracma::psi(1, (p+1)/4)) / 16
# Note: used for Integration by parts;

###
p = sqrt(3)
integrate(\(x) x^p * log(x)^2 / (x^2 + 1), 0, 1)
(pracma::psi(2, (p+1)/4 + 1/2) - pracma::psi(2, (p+1)/4)) / 64

###
p = sqrt(3)
integrate(\(x) x^p / (x^2 + 1)^2, 0, 1)
(digamma((p+1)/4) - digamma((p-1)/4)) * (p-1)/8 - 1/4


### Full:

### I( x^p * atan(x^n) * log(x) )
n = sqrt(5); p = sqrt(3);
# works only for n > 0;
integrate(\(x) x^p * atan(x^n) * log(x), 0, 1)
- (pracma::psi(1, ((p+1)/n + 3)/4) - pracma::psi(1, ((p+1)/n + 1)/4)) / (16*n*(p+1)) +
	+ (digamma(((p+1)/n + 3)/4) - digamma(((p+1)/n + 1)/4)) / (4*(p+1)^2) +
	- pi/(4*(p+1)^2);


### Special Case:
integrate(\(x) log(x) * atan(x), 0, 1)
pi^2 / 48 - pi/4 + log(2)/2


integrate(\(x) log(x) * atan(x) / x, 0, 1)
- pi^3 / 32

##################

### Variants: Atan * Log

### I( atan(x) * log(x) )
integrate(\(x) atan(x) * log(x), 0, 1)
- pi^2 / 48 + pi/4 - log(2)/2


### Type: log(1-x)

### I( atan(x) * log(1-x) )
integrate(\(x) atan(x) * log(1-x), 0, 1)
integrate(\(x) ((1-x)*log(1-x) + x) / (x^2+1) - pi/4, 0, 1)
integrate(\(x) (1-x)*log(1-x) / (x^2+1) - pi/4 + log(2)/2, 0, 1)
5/96 * pi^2 - log(2)^2 / 8 + pi*log(2)/8 - Catalan - pi/4 + log(2)/2;

### I( atan(x) * log(1-x) / x )
integrate(\(x) atan(x) * log(1-x) / x, 0, 1)
integrate(\(x) log(1-x) * log(x) / (x^2+1), 0, 1)$value +
	- pi^3 * 3/64 + Catalan * log(2)/2;
# TODO

### I( atan(x) * log(x) / (1-x) )
integrate(\(x) atan(x) * log(x) / (1-x), 0, 1)
integrate(\(x) 2 * log(1-x) * log(x) / (x^2+1), 0, 1)$value +
	- pi^3 * 3/64 + Catalan * log(2)/2;


### I( atan(x) / (1-x) )
integrate(\(x) atan(x) / (1-x) - pi/4/(1-x), 0, 1)
pi*log(2)/8 - Catalan

# Derived:
integrate(\(x) log(1-x) * log(x^2+1) / x^2, 0, 1)
pi*log(2)/4 - 2*Catalan - Re(log(1-1i)^2);


### I( log(x^2+1) * log(1-x) )
integrate(\(x) log(x^2+1) * log(1-x), 0, 1)
integrate(\(x) (1-x)*log(1-x) * 2*x/(x^2+1) + 2 - pi/2 - log(2), 0, 1)
integrate(\(x) 2*(x+1)*log(1-x) / (x^2+1) + 4 - pi/2 - log(2), 0, 1)
-5/48 * pi^2 + log(2)^2 / 4 + pi*log(2)/4 - 2*Catalan - pi/2 - log(2) + 4;

# Complex:
integrate(\(x) Re(log(x + 1i)) * log(1-x), 0, 1)
-5/96 * pi^2 + log(2)^2 / 8 + pi*log(2)/8 - Catalan - pi/4 - log(2)/2 + 2;
#
integrate(\(x) Im(log(x + 1i)) * log(1-x), 0, 1)
- (5/96 * pi^2 - log(2)^2 / 8 + pi*log(2)/8 - Catalan + pi/4 + log(2)/2);


### Type: log(1+x)

### I( atan(x) * log(1+x) )
integrate(\(x) atan(x) * log(1+x), 0, 1)
-1/96*pi^2 + 3/8*pi*log(2) - log(2)^2 / 8 + log(2)/2 - pi/4;

### I( atan(x) * log(1+x) / x )
integrate(\(x) atan(x) * log(1+x) / x, 0, 1)
integrate(\(x) - log(x) * log(x+1) / (x^2+1), 0, 1)$value +
	+ (pi/4)^3 - Catalan * log(2)/2;
# TODO

### I( x * atan(x) * log(1+x) )
integrate(\(x) x * atan(x) * log(1+x), 0, 1)
(pi*log(2)/2 - 5*log(2) + 3) / 4;


#################
### Fractions ###

### I( atan(x) * log(x) / (x^2+1) )
integrate(\(x) atan(x) * log(x) / (x^2+1), 0, 1)
7/16 * pracma::zeta(3) - pi * Catalan / 4


### I( x * atan(x) * log(x) / (x^2+1) )
integrate(\(x)  x * atan(x) * log(x) / (x^2+1), 0, 1)
integrate(\(x) -2 * log(1-x) * log(x) / (x^2+1), 0, 1, rel.tol=1E-12)$value + (pi/4)^3;
# TODO


### I( atan(x)^2 * log(x) )
integrate(\(x) atan(x)^2 * log(x), 0, 1)
integrate(\(x) - 2*x * atan(x) * log(x) / (x^2+1), 0, 1)$value +
	- (pi/4)^2 + Catalan - 1/4 * pi * log(2);
integrate(\(x) 2*log(x) * log(x^2+1) / (x^2+1), 0, 1)$value +
	- (pi/4)^2 + Catalan - 1/4 * pi * log(2) + 2*(Catalan * log(2) - (pi/4)^3);
# TODO

###
integrate(\(x) atan(x)^2 * log(1 - x), 0, 1)
# TODO

