########################
##
## Leonard Mada
## [the one and only]
##
## Integrals: Logarithms
## Log: Complex Fractions


### Helper Constants

Euler   = 0.57721566490153286060651209008240243079;
Catalan = 0.915965594177219015054603514;


################
################

### Base:
p = sqrt(2); n = sqrt(5); k = sqrt(3);
integrate(\(x) x^p / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n

###
integrate(\(x) x^p * log(x^n + 1) / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) *
	(digamma(k) - digamma(k - (p+1)/n)) / gamma(k) / n

###
p = sqrt(2); n = sqrt(5); k = sqrt(3);
b = sqrt(2)
integrate(\(x) x^p * log(x^n + b^n) / (x^n + b^n)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) *
	(digamma(k) - digamma(k - (p+1)/n) + n*log(b)) /
	(gamma(k) * n * b^(n*k - p - 1))

# Special Case: k = 1
p = 1/sqrt(2); n = sqrt(5);
b = sqrt(2)
integrate(\(x) x^p * log(x^n + b^n) / (x^n + b^n), 0, Inf)
- gamma((p+1)/n) * gamma(1 - (p+1)/n) *
	(digamma(1 - (p+1)/n) + Euler - n*log(b)) /
	(n * b^(n - p - 1))


###########

###
p = sqrt(2); n = sqrt(5); k = sqrt(3);
a = sqrt(3) - sqrt(2); b = sqrt(2)
integrate(\(x) x^p * log(x^n + a^n) / (x^n + b^n)^k, 0, Inf)
# TODO


#################
#################

### Simple:

### I( log(x + 1i)/(x - 1i)) )
integrate(\(x) Re(log(x + 1i)/(x - 1i)), 0, 1)
integrate(\(x) Im(log(x + 1i)/(x - 1i)), 0, 1)
-3/32*pi^2 + log(2)^2/8 + 1i*(5/8*pi*log(2) - Catalan);

# Helper:

# Li2((1−i)/2) # =
5*pi^2/96 - log(2)^2/8 + 1i*(pi/8*log(2) - Catalan);

# Li2((1+i)/2) # =
5*pi^2/96 - log(2)^2/8 - 1i*(pi/8*log(2) - Catalan);


####################

### Generalised:

### I( x^p * atan(x) / (x^2 + 1) )
p = - 1/sqrt(3);
integrate(\(x) x^p * atan(x) / (x^2 + 1), 0, Inf)
pi/4 / sin(pi*p/2) * (digamma(-p/2) - 2*digamma(-p) - Euler);

# Derivation:
pi/4 * gamma((p+1)/2) * gamma(1 - (p+1)/2) +
	- gamma((p+2)/2) * gamma(1 - (p+2)/2) *
	(digamma(-p/2) + Euler) / 4 +
	+ gamma(p+1) * gamma(-p) / 2 *
	((digamma(-p) + Euler - log(1i)) / (-1i)^p +
	(digamma(-p) + Euler - log(-1i)) / 1i^p );
# from below:

# Derivation:
p = - 1/sqrt(3);
integrate(\(x) x^p * (pi/2 - atan(x)) / (x^2 + 1), 0, Inf)
integrate(\(x) x^p * Re(log((x+1i)/(x-1i)) / 2i) / (x^2 + 1), 0, Inf)
gamma((p+2)/2) * gamma(1 - (p+2)/2) *
	(digamma(-p/2) + Euler) / 4 +
	- gamma(p+1) * gamma(-p) / 2 *
	((digamma(-p) + Euler - log(1i)) / (-1i)^p +
	(digamma(-p) + Euler - log(-1i)) / 1i^p );

###
p = - 1/sqrt(3);
integrate(\(x) Re(x^p * log(x - 1i) / (x^2 + 1)), 0, Inf)
integrate(\(x) Im(x^p * log(x - 1i) / (x^2 + 1)), 0, Inf)
(gamma((p+2)/2) * gamma(1 - (p+2)/2) *
	(digamma(-p/2) + Euler) +
	- 1i*(gamma((p+1)/2) * gamma(1 - (p+1)/2) *
	(digamma(1 - (p+1)/2) + Euler))) / 4i +
	- gamma(p+1) * gamma(-p) / 2i *
	((digamma(-p) + Euler - log(1i)) / (-1i)^p +
	(digamma(-p) + Euler - log(-1i)) / 1i^p );


###
n = 2; k = 1; # !!!
p = - 1/sqrt(3);
integrate(\(x) Re(x^p * log(x - 1i) / (x + 1i)), 0, Inf)
integrate(\(x) Im(x^p * log(x - 1i) / (x + 1i)), 0, Inf)
(gamma((p+1 + k)/n) * gamma(k - (p+1 + k)/n) *
	(digamma(k) - digamma(k - (p+1 + k)/n)) +
	- 1i*(gamma((p+1)/n) * gamma(k - (p+1)/n) *
	(digamma(k) - digamma(k - (p+1)/n)))) / gamma(k) / n +
	+ gamma(p+1) * gamma(-p) *
	(digamma(-p) + Euler - log(1i)) / 1i^(- p);


###
n = 2; k = 1; # !!!
p = - 1/sqrt(3);
integrate(\(x) Re(x^p * log(x^2 + 1) / (x + 1i)), 0, Inf)
integrate(\(x) Im(x^p * log(x^2 + 1) / (x + 1i)), 0, Inf)
(gamma((p+1 + k)/n) * gamma(k - (p+1 + k)/n)*
	(digamma(k) - digamma(k - (p+1 + k)/n)) +
	- 1i*(gamma((p+1)/n) * gamma(k - (p+1)/n)*
	(digamma(k) - digamma(k - (p+1)/n)))) / gamma(k) / n

# TODO: log(x - 1i) / (x + 1i)


###
n = 4; k = 1; # !!!
p = 1/sqrt(3);
integrate(\(x) Re(x^p * log(x^4 + 1) / (x^2 + 1i)), 0, Inf)
integrate(\(x) Im(x^p * log(x^4 + 1) / (x^2 + 1i)), 0, Inf)
(gamma((p+1 + 2*k)/n) * gamma(k - (p+1 + 2*k)/n)*
	(digamma(k) - digamma(k - (p+1 + 2*k)/n)) +
	- 1i*(gamma((p+1)/n) * gamma(k - (p+1)/n)*
	(digamma(k) - digamma(k - (p+1)/n)))) / gamma(k) / n


#############

### I( log(x^2 + 1i) / (x^2 + 1) )
integrate(function(x) Im(log(x^2 + 1i)) / (x^2 + 1), 0, Inf)
pi^2 / 8
integrate(function(x) Re(log(x^2 + 1i)) / (x^2 + 1), 0, Inf)
pi*log(2*cos(pi/8))

# Derivation:
integrate(function(x) atan(x^2) / (x^2 + 1), 0, Inf)
pi^2 / 8
integrate(function(x) Im(log((x^2 + 1i)/(x^2 - 1i))) / (x^2 + 1), 0, Inf)
pi^2 / 4

integrate(function(x) log(x^4 + 1)/(x^2 + 1), 0, Inf)
2*pi*log(2*cos(pi/8))

integrate(function(x) atan(x^4)/(x^2 + 1), 0, Inf)
pi^2 / 8

###############
###############

############
### ATAN ###

# Note:
# - see also Integrals.Trig.Tan.Atan.R;


### I( atan(k * x^2) / (x^2 + 1i) )
k = 5^(1/3)
integrate(\(x) Re(atan(k*x^2) / (x^2 + 1i)), 0, Inf, rel.tol=1E-9)
integrate(\(x) Im(atan(k*x^2) / (x^2 + 1i)), 0, Inf, rel.tol=1E-9)
x = 1/k^(1/2);
pi^2 / sin(3*pi/4) / 8 +
	+ (2*log(x + 1) - log(x^2 + 1) - 2*atan(x)) * pi / sin(pi/4) / 8 +
	- 1i * pi^2 / sin(pi/4) / 8 +
	- 1i * (log(x^2 + 1) - 2*log(x + 1) - 2*atan(x)) * pi / sin(pi/4) / 8;


####################
####################

### on [0, 1]

### I( log(x) * log(1-x) )
integrate(\(x) log(x) * log(1-x), 0, 1)
2 - pi^2/6

###
integrate(\(x) atan(x) * atan(1-x), 0, 1)
# TODO


##############

### Fractions: Pow = 2

### I( log(x+1i) / (x^2+1) )
integrate(\(x) Re(log(x+1i)) / (x^2+1), 0, 1)
integrate(\(x) Im(log(x+1i)) / (x^2+1), 0, 1)
(pracma::psi(1, 3/4) - pracma::psi(1, 1/4)) / 32 +
	+ (digamma(1) - digamma(1/2)) * pi/8 +
	+ pracma::psi(1, 1/2) * 3/16 * 1i;

# Note:
# (pracma::psi(1, 3/4) - pracma::psi(1, 1/4)) / 16 == Catalan;
# (digamma(1) - digamma(1/2)) == 2*log(2);
# I( log(x - 1i) ... ) = complex conjugate;


### I( log(x^2+1) / (x^2+1) )
integrate(\(x) log(x^2+1) / (x^2+1), 0, 1)
(pracma::psi(1, 3/4) - pracma::psi(1, 1/4)) / 16 + pi*log(2)/2;

###
integrate(\(x) log(x+1) / (x^2+1), 0, 1)
pi * log(2)/8

###
integrate(\(x) log(1-x) / (x^2+1), 0, 1)
(pracma::psi(1, 3/4) - pracma::psi(1, 1/4)) / 16 + pi*log(2)/8;


### from log(x+b) / (x^2+1)
b = sqrt(5)
integrate(\(x) 1 / (x+b) / (x^2+1), 0, 1)
integrate(\(x) (1/(x+b) - (x-b)/(x^2+1)) / (b^2+1), 0, 1)
(log(b+1) - log(b) - log(2)/2 + pi/4 * b) / (b^2+1)


##############

### Pow = 3
# see also Integrals.Log.Fractions.Other.R;

### I( log(1-x) / (x^2-x+1) )
integrate(\(x) log(1-x) / (x^2-x+1), 0, 1)
- (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 6

### I( log(1+x) / (x^2-x+1) )
integrate(\(x) log(1+x) / (x^2-x+1), 0, 1)
# TODO


##############

### Pow = 4
# see also Integrals.Log.Fractions.P4.R;

### I( x^p * log(x + 1i) / (x^4 + 1) )

### I( log(x + 1i) / (x^4 + 1) )
integrate(\(x) Re(log(x + 1i)) / (x^4 + 1), 0, 1)
integrate(\(x) Im(log(x + 1i)) / (x^4 + 1), 0, 1)

#
integrate(\(x) Re((1+1i) * log(x + 1i)) / (x^4 + 1), 0, 1)
(pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 64 +
	- (digamma(3/4) - digamma(1/4)) * # == pi
		(digamma(7/8) - digamma(3/8)) / 64 +
	+ (digamma(5/8) - digamma(1/8)) * log(2) / 32;


### I( x^2 * log(x + 1i) / (x^4 + 1) )
integrate(\(x) x^2 * Re((1-1i) * log(x + 1i)) / (x^4 + 1), 0, 1)
(pracma::psi(1, 7/8) - pracma::psi(1, 3/8)) / 64 +
	+ (digamma(5/8) - digamma(1/8)) * pi / 64 +
	+ (digamma(7/8) - digamma(3/8)) * log(2) / 32;

### Im()
# TODO


### I( log(x+1i) / (x^2+1i) )
integrate(\(x) Re(log(x+1i) / (x^2+1i)), 0, 1, rel.tol=1E-13)
integrate(\(x) Im(log(x+1i) / (x^2+1i)), 0, 1, rel.tol=1E-13)
(pracma::psi(1, 7/8) - pracma::psi(1, 3/8)) / 64 +
	- (pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 64 +
	+ (digamma(7/8) - digamma(3/8)) * log(2) / 32 +
	+ (digamma(5/8) - digamma(1/8)) * pi / 16 +
	- (digamma(5/8) - digamma(1/8)) *
		(digamma(3/4) - digamma(1/4)) / 64 +
	+ (digamma(1/4) - digamma(3/4)) / sin(pi/4) * pi/16 +
	# Im:
	+ ((digamma(7/8) - digamma(3/8)) * pi / 64 * 3 +
		- (digamma(5/8) - digamma(1/8)) * log(2) / 32) * 1i;

#
integrate(\(x) Im((1+1i) * log(x + 1i)) / (x^4 + 1), 0, 1)
# TODO

integrate(\(x) x^2 * Im((1-1i) * log(x + 1i)) / (x^4 + 1), 0, 1)
# TODO

integrate(\(x) (x^2 - 1) * Im(log(x + 1i)) / (x^4 + 1), 0, 1)
(pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 64 +
	+ (digamma(7/8) - digamma(3/8)) * pi * 3/64 +
	+ (digamma(5/8) - digamma(1/8)) * pi / 64 +
	+ (digamma(7/8) - digamma(3/8)) * log(2) / 32 +
	+ (digamma(5/8) - digamma(1/8)) * log(2) / 32 +
	- sqrt(2)*pi/16 * (pi/2 + log(2));


### I( x * log(x + 1) / (x^2 - 1i) )
integrate(\(x) x * log(x + 1) * Re(1/(x^2 - 1i)), 0, 1)
integrate(\(x) x * log(x + 1) * Im(1/(x^2 - 1i)), 0, 1)
integrate(\(x) x^3 * log(x+1) / (x^4+1), 0, 1)$value +
integrate(\(x)   x * log(x+1) / (x^4+1), 0, 1)$value * 1i;
(pracma::psi(1, 1) - pracma::psi(1, 1/2)) / 64 +
	+ log(2)^2 / 32 +
	+ (digamma(5/8) - digamma(1/8)) *
		(digamma(7/8) - digamma(3/8)) / 64 +
	- (digamma(3/4) - digamma(1/4))^2 / 128 +
1i * ((pracma::psi(1, 3/4) - pracma::psi(1, 1/4)) / 64 +
	+ (digamma(3/4) - digamma(1/4)) * log(2) / 32 + # pi*log(2)
	- (digamma(7/8) - digamma(3/8))^2 / 128 +
	+ (digamma(5/8) - digamma(1/8))^2 / 128);


### Derivation:

# Base: see files:
# Integrals.Log.Fractions.R;
# Integrals.Log.Fractions.P4.R;


# Note:
# - Real-variants moved to file:
#   Integrals.Log.Fractions.P4.R;


### Helper:

### from x^2 * log(x^2+b^2) / (x^4+1)
b = sqrt(3)
# Note: factor * 2*b NOT included;
integrate(\(x) x^2 / (x^2+b^2) / (x^4+1), 0, 1)
integrate(\(x) (-b^2/(x^2+b^2) + (b^2*x^2 + 1)/(x^4+1)) / (b^4 + 1), 0, 1)
(-b*atan(1/b) + (digamma(5/8) - digamma(1/8))/8 +
	+ b^2*(digamma(7/8) - digamma(3/8))/8) / (b^4 + 1);
(b*atan(b) - b*pi/2 + (digamma(5/8) - digamma(1/8))/8 +
	+ b^2*(digamma(7/8) - digamma(3/8))/8) / (b^4 + 1);

### from log(x^2+b^2) / (x^4+1)
b = sqrt(3)
# Note: factor * 2*b NOT included;
integrate(\(x) 1 / (x^2+b^2) / (x^4+1), 0, 1)
integrate(\(x) (1/(x^2+b^2) - (x^2 - b^2)/(x^4+1)) / (b^4 + 1), 0, 1)
(atan(1/b)/b - (digamma(7/8) - digamma(3/8))/8 +
	+ b^2*(digamma(5/8) - digamma(1/8))/8) / (b^4 + 1);
(pi/2 / b - atan(b)/b - (digamma(7/8) - digamma(3/8))/8 +
	+ b^2*(digamma(5/8) - digamma(1/8))/8) / (b^4 + 1);

# [duplicate]
integrate(\(x) (Re(log(x+1i)) - Im(log(x+1i))) / (x^4 + 1), 0, 1)
integrate(\(x) Re((1+1i) * log(x + 1i)) / (x^4 + 1), 0, 1) # same as previous;
integrate(\(x) (Im(log(x+1i)) - Re(log(x+1i))) / (x^4 + 1), 0, 1) # * 1i
( (pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 64 +
+ (digamma(5/8) - digamma(1/8)) * log(2) / 32 +
- (digamma(7/8) - digamma(3/8)) *
	(digamma(3/4) - digamma(1/4)) / 64) * (1-1i);

### Other
integrate(\(x) Re((1-1i) * log(x + 1i)) / (x^4 + 1), 0, 1)
integrate(\(x) (log(x^2 + 1)/2 - atan(x)) / (x^4 + 1), 0, 1)$value +
	+ (digamma(5/8) - digamma(1/8)) * pi / 16;
# TODO: ?


### from log(x+b) / (x^4 + 1)
b = sqrt(3)
integrate(\(x) 1 / (x+b) / (x^4+1), 0, 1)
integrate(\(x) (1/(x+b) - (x - b)*(x^2 + b^2)/(x^4+1)) / (b^4 + 1), 0, 1)
integrate(\(x) (1/(x+b) - (x^3 - b*x^2 + b^2*x - b^3)/(x^4+1)) / (b^4 + 1), 0, 1)
(log(b+1) - log(b) - log(2)/4 +
	+ (digamma(7/8) - digamma(3/8)) * b / 8 +
	- (digamma(3/4) - digamma(1/4)) * b^2 / 8 +
	+ (digamma(5/8) - digamma(1/8)) * b^3 / 8 ) / (b^4 + 1);
# Note:
# - back-integration over [0, 1i] is redundant;
#   (see duplicated result above)
# - back-integration over [-1, 0] generates variant log(1-x);


### from x^2 * log(x+b) / (x^4 + 1)
b = sqrt(5)
integrate(\(x) x^2 / (x+b) / (x^4+1), 0, 1)
integrate(\(x) x^2 * (1/(x+b) - (x - b)*(x^2 + b^2)/(x^4+1)) / (b^4 + 1), 0, 1)
integrate(\(x) (b^2/(x+b) - (b^2*x^3 - b^3*x^2 - x + b)/(x^4+1)) / (b^4 + 1), 0, 1)
(b^2*log(b+1) - b^2*log(b) - log(2) * b^2 / 4 +
	+ (digamma(7/8) - digamma(3/8)) * b^3 / 8 +
	+ (digamma(3/4) - digamma(1/4)) / 8 +
	- (digamma(5/8) - digamma(1/8)) * b / 8 ) / (b^4 + 1);


### from x * log(x+b) / (x^4 + 1)
b = sqrt(3)
integrate(\(x) x / (x+b) / (x^4+1), 0, 1)
integrate(\(x) (x/(x+b) - x*(x^3 - b*x^2 + b^2*x - b^3)/(x^4+1)) / (b^4 + 1), 0, 1)
integrate(\(x) (-b/(x+b) + (b*x^3 - b^2*x^2 + b^3*x + 1)/(x^4+1)) / (b^4 + 1), 0, 1)
(-b*log(b+1) + b*log(b) + b*log(2)/4 +
	- (digamma(7/8) - digamma(3/8)) * b^2 / 8 +
	+ (digamma(3/4) - digamma(1/4)) * b^3 / 8 +
	+ (digamma(5/8) - digamma(1/8)) / 8 ) / (b^4 + 1);


### from x^3 * log(x+b) / (x^4 + 1)
b = sqrt(3)
integrate(\(x) x^3 / (x+b) / (x^4+1), 0, 1)
integrate(\(x) (x^3/(x+b) +
	- (x^6 - b*x^5 + b^2*x^4 - b^3*x^3)/(x^4+1)) / (b^4 + 1), 0, 1)
integrate(\(x) (- b^3/(x+b) +
	+ (b^3*x^3 + x^2 - b*x + b^2)/(x^4+1)) / (b^4 + 1), 0, 1)
(-b^3*log(b+1) + b^3*log(b) + b^3*log(2)/4 +
	+ (digamma(7/8) - digamma(3/8)) / 8 +
	- (digamma(3/4) - digamma(1/4)) * b / 8 +
	+ (digamma(5/8) - digamma(1/8)) * b^2 / 8) / (b^4 + 1);


#################

### Pow = 6

### I( log(x^2+1) / (x^6+1) )
integrate(\(x) log(x^2+1) / (x^6+1), 0, 1)
integrate(\(x) 2 * (atan(x) - pi/4) / (x^6 - 1), 0, 1)$value +
(pracma::psi(1, 7/12) - pracma::psi(1, 1/12)) / 72 +
	+ (digamma(1/3) - digamma(1/6)) *
		(digamma(11/12) - digamma(5/12)) / 36 +
	- (digamma(2/3) - digamma(1/6)) *
		(digamma(3/4) - digamma(1/4)) / 36 +
	+ (digamma(1) - digamma(1/6)) *
		(digamma(7/12) - digamma(1/12)) / 36;


# Note: I( (1 - x^p) / (1 - x^n) ) on [0,1]
# (digamma((p+1)/n) - digamma(1/n)) / n;


### from log(x^2+b^2) / (x^6 + 1)
# Note: does NOT include the factor * 2*b;
b = sqrt(3)
integrate(\(x) 1 / (x^2+b^2) / (x^6+1), 0, 1)
integrate(\(x) ((x^4 - b^2*x^2 + b^4)/(x^6+1) +
	- 1/(x^2+b^2)) / (b^6 - 1), 0, 1)
(- pi/(2*b) + atan(b)/b +
	+ (digamma(11/12) - digamma(5/12)) / 12 +
	- (digamma(3/4) - digamma(1/4)) * b^2 / 12 +
	+ (digamma(7/12) - digamma(1/12)) * b^4 / 12) / (b^6 - 1)


### ATAN

### I( (atan(x) - pi/4) / (x^2 - 1) )
integrate(\(x) (atan(x) - pi/4) / (x^2 - 1), 0, 1)
Catalan / 2


### I( (atan(x) - pi/4) / (x^4 - 1) )
integrate(\(x) (atan(x) - pi/4) / (x^4 - 1), 0, 1)
pi^2/64 + Catalan/4


### I( (atan(x) - pi/4) / (x^6 - 1) )
integrate(\(x) (atan(x) - pi/4) / (x^6 - 1), 0, 1)
integrate(\(x) Catalan/6 - 1/3 * (atan(x) - pi/4) * (x^2+2) / (x^4 + x^2 + 1), 0, 1)
integrate(\(x) - 1/3 * atan(x) * (x^2+2) / (x^4 + x^2 + 1), 0, 1)$value +
	+ Catalan/6 + pi/12 * (
		(digamma(5/6) - digamma(1/6)) / 6 +
		(digamma(1/2) - digamma(1/6)) / 6 );
# TODO

