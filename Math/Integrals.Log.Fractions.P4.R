########################
##
## Leonard Mada
## [the one and only]
##
## Integrals: Logarithms
## Log-Fractions: Power = 4
##
## draft v.0.1b


##################
### Logarithms ###

### Fractions: Power = 4

# I( x^p * log(1 + x) / (x^4 + 1) )
# I( x^p * log(1 - x) / (x^4 + 1) )


### History:

# [refactor]
# - moved Pow = 4 to this file from file:
#   Integrals.Log.Fractions.Other.R;
# - moved specific variants from:
#   Integrals.Log.Fractions.Complex.R;
# - TODO: move also the Derivation;

####################

### Helper Constants

Catalan = 0.915965594177219015054603514;


##################

##################
### Log(1 - x) ###

### I( log(1-x) / (x^4+1) )
integrate(\(x) log(1-x) / (x^4+1), 0, 1)
(pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 64 +
	+ (digamma(3/4) - digamma(1/4)) * # == pi
		(digamma(7/8) - digamma(3/8)) / 64 +
	+ (digamma(5/8) - digamma(1/8)) * log(2) / 32;

# on [1, Inf]
integrate(\(x) log(x-1) / (x^4+1), 1, Inf)
(pracma::psi(1, 7/8) - pracma::psi(1, 3/8)) / 64 +
- (pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 64 +
	+ (digamma(7/8) - digamma(3/8)) * log(2) / 32 +
	- (digamma(5/8) - digamma(1/8)) *
		(digamma(3/4) - digamma(1/4)) / 64 +
	+ (digamma(1/4) - digamma(3/4)) / sin(pi/4) * pi/16;


### I( x * log(1-x) / (x^4+1) )
integrate(\(x) x * log(1-x) / (x^4+1), 0, 1)
(pracma::psi(1, 3/4) - pracma::psi(1, 1/4)) / 64 + # Catalan / 4
	+ (digamma(3/4) - digamma(1/4)) * log(2) / 32 + # pi*log(2)
	+ (digamma(7/8) - digamma(3/8))^2 / 128 +
	- (digamma(5/8) - digamma(1/8))^2 / 128;

# on [1, Inf]
integrate(\(x) x * log(x-1) / (x^4+1), 1, Inf)
(digamma(3/4) - digamma(1/4)) * log(2) / 32 + # pi*log(2)
	+ (digamma(7/8) - digamma(3/8))^2 / 128 +
	- (digamma(5/8) - digamma(1/8))^2 / 128;


### I( x^2 * log(1-x) / (x^4 + 1) )
integrate(\(x) x^2 * log(1-x) / (x^4+1), 0, 1)
(pracma::psi(1, 7/8) - pracma::psi(1, 3/8)) / 64 +
	+ (digamma(7/8) - digamma(3/8)) * log(2) / 32 +
	- (digamma(5/8) - digamma(1/8)) *
		(digamma(3/4) - digamma(1/4)) / 64;

# on [1, Inf]
integrate(\(x) x^2 * log(x-1) / (x^4+1), 1, Inf)
(pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 64 +
- (pracma::psi(1, 7/8) - pracma::psi(1, 3/8)) / 64 +
	- pi^2 * cos(3*pi/4) / sin(3*pi/4)^2 / 16 +
	+ (digamma(3/4) - digamma(1/4)) * # == pi
		(digamma(7/8) - digamma(3/8)) / 64 +
	+ (digamma(5/8) - digamma(1/8)) * log(2) / 32;


### I( x^3 * log(1-x) / (x^4 + 1) )
integrate(\(x) x^3 * log(1-x) / (x^4+1), 0, 1)
- (pracma::psi(1, 1) - pracma::psi(1, 1/2)) / 64 +
	- 5/(3*64) * pi^2 + log(2)^2 / 32 +
	- (digamma(5/8) - digamma(1/8)) *
		(digamma(7/8) - digamma(3/8)) / 64 +
	+ (digamma(3/4) - digamma(1/4))^2 / 128;


### I( log(1-x) / x / (x^4+1) )
integrate(\(x) log(abs(1-x)) / x / (x^4+1), 0, 1)$value +
integrate(\(x) log(abs(1-x)) / x / (x^4+1), 1, Inf)$value;
- pi^2 * 31 / (3*64);

# on [0, 1]
integrate(\(x) log(abs(1-x)) / x / (x^4+1), 0, 1)
(pracma::psi(1, 1) - pracma::psi(1, 1/2)) / 64 +
	- pi^2 * 9/64 - log(2)^2 / 32 +
	+ (digamma(5/8) - digamma(1/8)) *
		(digamma(7/8) - digamma(3/8)) / 64 +
	- (digamma(3/4) - digamma(1/4))^2 / 128;


##################
### Log(1 + x) ###

### on [0, Inf] & [0. 1]

### I( log(1+x) / (x^4 + 1) )
integrate(\(x) log(1+x) / (x^4 + 1), 0, Inf)
- (pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 64 +
	- (digamma(3/4) - digamma(1/4)) / sin(pi/4) * pi / 16 +
	+ (digamma(3/4) - digamma(1/4)) / sin(3*pi/4) * pi / 32 +
	- (digamma(7/8) - digamma(3/8)) * pi / 32 +
	+ log(2) / sin(pi/4) * pi / 16;

# on [0, 1]
# TODO


### I( x * log(1+x) / (x^4 + 1) )
integrate(\(x) x * log(1+x) / (x^4+1), 0, Inf)
(pracma::psi(1, 3/4) - pracma::psi(1, 1/4)) / 64 +
	+ (digamma(3/4) - digamma(1/4)) * log(2) / 16 + # pi*log(2)
	- (digamma(7/8) - digamma(3/8))^2 / 64 +
	+ (digamma(5/8) - digamma(1/8))^2 / 64;

# on [0, 1]
integrate(\(x) x * log(1+x) / (x^4+1), 0, 1)
(pracma::psi(1, 3/4) - pracma::psi(1, 1/4)) / 64 +
	+ (digamma(3/4) - digamma(1/4)) * log(2) / 32 + # pi*log(2)
	- (digamma(7/8) - digamma(3/8))^2 / 128 +
	+ (digamma(5/8) - digamma(1/8))^2 / 128;


### I( x^2 * log(1+x) / (x^4 + 1) )
integrate(\(x) x^2 * log(1+x) / (x^4 + 1), 0, Inf)
- (pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 64 +
	- pi^2 * cos(3*pi/4) / sin(3*pi/4)^2 / 16 +
	- (digamma(3/4) - digamma(1/4)) / sin(pi/4) * pi / 16 +
	+ (digamma(3/4) - digamma(1/4)) / sin(3*pi/4) * pi / 32 +
	- (digamma(7/8) - digamma(3/8)) * pi / 32 +
	+ log(2) / sin(pi/4) * pi / 16;


### I( x^3 * log(x+1) / (x^4+1) )
integrate(\(x) x^3 * log(x+1) / (x^4+1), 0, 1)
(pracma::psi(1, 1) - pracma::psi(1, 1/2)) / 64 +
	+ log(2)^2 / 32 +
	+ (digamma(5/8) - digamma(1/8)) *
		(digamma(7/8) - digamma(3/8)) / 64 +
	- (digamma(3/4) - digamma(1/4))^2 / 128;


### I( log(1+x) / x / (x^4+1) )
integrate(\(x) log(1+x) / x / (x^4+1), 0, Inf)
pi^2 / (6*32) + pi^2 / 12;
pi^2 * 17/(3*64);

### on [0, 1]
integrate(\(x) log(1+x) / x / (x^4+1), 0, 1)
- (pracma::psi(1, 1) - pracma::psi(1, 1/2)) / 64 +
	+ pi^2 / 12 - log(2)^2 / 32 +
	- (digamma(5/8) - digamma(1/8)) *
		(digamma(7/8) - digamma(3/8)) / 64 +
	+ (digamma(3/4) - digamma(1/4))^2 / 128;


### Helper:
integrate(\(x) log(1+x) / x, 0, 1)
pi^2 / 12

# Limit:
e = 1E-4
integrate(\(x) log(x) / x / (x^4+1), 1, Inf)
(pracma::psi(1, e/8) - pracma::psi(1, (e+4)/8)) / 64 +
	- pi^2 * cos(e*pi/4) / sin(e*pi/4)^2 / 16;
(16/e^2 - pi^2 * cos(e*pi/4) / sin(e*pi/4)^2) / 32;
pi^2 / (6*32)


###################

### I( log(x^2 + 1) / (x^4 + 1) )

### I( log(x^2 + 1) / (x^4 + 1) )
integrate(\(x) log(x^2+1) / (x^4+1), 0, 1)
integrate(\(x) -2*atan(x) / (x^4+1), 0, 1)$value +
	+ (pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 32 +
	+ (digamma(5/8) - digamma(1/8)) * pi/8 +
	- (digamma(3/4) - digamma(1/4)) *
		(digamma(7/8) - digamma(3/8)) / 32 +
	+ (digamma(5/8) - digamma(1/8)) * log(2) / 16;
# TODO

### I( x^2 * log(x^2 + 1) / (x^4 + 1) )
integrate(\(x) x^2 * log(x^2+1) / (x^4+1), 0, 1)
integrate(\(x) 2*x^2 * atan(x) / (x^4+1), 0, 1)$value +
	+ (pracma::psi(1, 7/8) - pracma::psi(1, 3/8)) / 32 +
	- (digamma(7/8) - digamma(3/8)) * pi / 8 +
	+ (digamma(3/4) - digamma(1/4)) *
		(digamma(5/8) - digamma(1/8)) / 32 +
	+ (digamma(7/8) - digamma(3/8)) * log(2) / 16;
# TODO

### I( (x^2 + 1) * log(x^2 + 1) / (x^4 + 1) )
integrate(\(x) (x^2 + 1) * log(x^2 + 1) / (x^4 + 1), 0, 1)
sqrt(2)*pi/8 * (pi/2 + log(2)) +
	+ 2*(pracma::psi(1, 7/8) - pracma::psi(1, 3/8)) / 64 +
	- pi/2 * (atan(exp(1i*pi/4))*exp(1i*pi/4) +
		+ atan(exp(-1i*pi/4))*exp(-1i*pi/4));


# Note:
(digamma(3/4) - digamma(1/4)) # == pi;
#
sum(atan(exp(c(1i,-1i)*pi/4))*exp(c(1i,-1i)*pi/4)) # ==
(digamma(7/8) - digamma(3/8)) / 4
#
sum(atan(exp(c(1i,-1i)*pi/4))*exp(c(-1i,1i)*pi/4)) # ==
(digamma(5/8) - digamma(1/8)) / 4


#############

##########
### Other:

### I( log(1+x^n) / x / (x^(4*n)+1) )
n = sqrt(3);
integrate(\(x) log(1+x^n) / x / (x^(4*n)+1), 0, 1)
((pracma::psi(1, 1/2) - pracma::psi(1, 1)) / 64 +
	+ pi^2 / 12 - log(2)^2 / 32 +
	- (digamma(5/8) - digamma(1/8)) *
		(digamma(7/8) - digamma(3/8)) / 64 +
	+ (digamma(3/4) - digamma(1/4))^2 / 128 ) / n;


### I( log(x) / ((x-1)^4 + 1) )
integrate(\(x) log(x) / ((x-1)^4 + 1), 1, Inf)
integrate(\(x) log(1+x) / (x^4 + 1), 0, Inf)
# Note: see previous integral;
# Variation =>
integrate(\(x) - x^2 * log(x) / (x^4 + (x-1)^4), 0, 1)

###
integrate(\(x) log(x) / ((x-1)^4 + 1), 0, 1)
integrate(\(x) log(1-x) / (x^4+1), 0, 1)
(pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 64 +
	+ (digamma(5/8) - digamma(1/8)) * log(2) / 32 +
	+ (digamma(7/8) - digamma(3/8)) *
		(digamma(3/4) - digamma(1/4)) / 64;


###############
### Derivation:

### Base:
# see file: Integrals.Log.Fractions.R;
b = sqrt(pi)
integrate(\(x) log(x^2 + b^2) / (x^4 + 1), 0, Inf)
sqrt(2)*pi/8 * (2*atan(1/b^2) + log(b^4 + 1)) +
	- pi/2 * (atan(1/b*exp(1i*pi/4))*exp(1i*pi/4) +
		+ atan(1/b*exp(-1i*pi/4))*exp(-1i*pi/4));

###
b = sqrt(5)
integrate(\(x) x^2 * log(x^2 + b^2) / (x^4 + 1), 0, Inf)
sqrt(2)*pi/8 * (2*atan(b^2) + log(b^4 + 1)) +
	- pi^2*cos(3*pi/4)/sin(3*pi/4)^2 / 8 +
	- pi/2 * (atan(b*exp(1i*pi/4))*exp(1i*pi/4) +
		+ atan(b*exp(-1i*pi/4))*exp(-1i*pi/4));


### on [0, Inf]
# from log(x+b) / (x^4 + 1)
integrate(\(x) 1/(x+b) / (1 + x^4), 0, Inf)
integrate(\(x) (1/(x+b) - (x^3-b*x^2+b^2*x-b^3)/ (1 + x^4)) / (b^4+1), 0, Inf)
(-log(b) + pi*(b/sin(3*pi/4) - b^2 + b^3/sin(pi/4))/4) / (b^4+1)

### on [0, 1]
# from log(x^2+b^2) / (x^4 + 1)
# Note: does NOT include factor * 2*b;
b = sqrt(5)
integrate(\(x) 1/(x^2 + b^2) / (x^4 + 1), 0, 1)
integrate(\(x) (1/(x^2 + b^2) - (x^2-b^2)/ (x^4 + 1)) / (b^4+1), 0, 1)
(pi/(2*b) - atan(b)/b - (digamma(7/8) - digamma(3/8))/8 +
	+ b^2 * (digamma(5/8) - digamma(1/8)) / 8) / (b^4+1);

### Helper
integrate(\(x) x * log(1-x^4) / (x^4+1), 0, 1)
pi*log(2)/8 - Catalan / 2

#
integrate(\(x) x * log(1-x^2) / (x^4+1), 0, 1)
pi*log(2)/16 - Catalan/2
# Diff( x*log(1-x^2), x*log(1+x) )
# => I( x*log(1-x) / ... );

#
integrate(\(x) x^3 * log(1-x^2) / (x^4+1), 0, 1)
-5/(3*64) * pi^2 + log(2)^2 / 16;


# from x * log(x - b) / (x^4 + 1)
b = 1 / sqrt(3)
integrate(\(x) x * 1/(x - b) / (x^4+1), 1, Inf)
integrate(\(x) x * (1/(x - b) - (x+b)*(x^2+b^2)/ (x^4+1)) / (b^4+1), 1, Inf)
integrate(\(x) (b/(x - b) - (b*x^3+b^2*x^2+b^3*x-1)/ (x^4+1)) / (b^4+1), 1, Inf)
(- b*(log(1-b) - log(2)/4) +
	- b^2*pi/sin(3*pi/4)/4 + b^2*(digamma(7/8) - digamma(3/8)) / 8 +
	- b^3*pi/4 + b^3*(digamma(3/4) - digamma(1/4)) / 8 +
	+ pi/sin(pi/4)/4 - (digamma(5/8) - digamma(1/8)) / 8) / (b^4+1);


#####################
#####################

### Variants: Pow = 4

### I( log(x^4 - x^2 + 1) / (x^4 + 1) )
integrate(\(x) log(x^4 - x^2 + 1) / (x^4 + 1), 0, Inf)
b = cos(pi/3) + c(1i,-1i)*sin(pi/3);
sum(sqrt(2)*pi/8 * (2*atan(1/b^2) + log(b^4 + 1)) +
	- pi/2 * (atan(1/b*exp(1i*pi/4))*exp(1i*pi/4) +
		+ atan(1/b*exp(-1i*pi/4))*exp(-1i*pi/4)) );


### Gen: I( log(x^4 + 2*cos(a)*x^2 + 1) / (x^4 + 1) )
a = 1/pi;
integrate(\(x) log(x^4 + 2*cos(a)*x^2 + 1) / (x^4 + 1), 0, Inf)
b = cos(a/2) + c(1i,-1i)*sin(a/2);
sum(sqrt(2)*pi/8 * (2*atan(1/b^2) + log(b^4 + 1)) +
	- pi/2 * (atan(1/b*exp(1i*pi/4))*exp(1i*pi/4) +
		+ atan(1/b*exp(-1i*pi/4))*exp(-1i*pi/4)) );


### I( log(Poly(x^2)) / (x^4 + 1) )
integrate(\(x) log(x^10 + 9*x^8 + 28*x^6 + 35*x^4 + 15*x^2 + 1) / (x^4 + 1), 0, Inf)
cs = 2*cos(2*seq(5)*pi/11);
b  = abs(cs);
sum(sqrt(2)*pi/8 * (2*atan(1/b^2) + log(b^4 + 1)) +
	- pi/2 * (atan(1/b*exp(1i*pi/4))*exp(1i*pi/4) +
		+ atan(1/b*exp(-1i*pi/4))*exp(-1i*pi/4)) );
# Test:
x = - cs^2;
x^5 + 9*x^4 + 28*x^3 + 35*x^2 + 15*x + 1 # == 0

