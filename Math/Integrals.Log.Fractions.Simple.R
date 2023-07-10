


### Helper Functions

# I( x^p / (x^n + 1) ) on [0, 1]
int.FrU01 = function(n, p=0) {
	(digamma(((p+1)/n + 1)/2) - digamma((p+1)/n/2)) / (2*n);
}
# I( x^p * log(x) / (x^n + 1) ) on [0, 1]
int.LogFrU01 = function(n, p=0) {
	(pracma::psi(1, (p+1)/n) - pracma::psi(1, (p+1)/(2*n))/2) / n^2;
}

zeta = function(n) {
	(1 + 1/(2^(n-1) - 1)) * pracma::eta(n);
}

### Fractions: Polynomials
int.FrUInf = function(n, p=0, pow=1, coeff=1) {
	k = 1/pow;
	tmp = sapply(p, function(p) {
		gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) / n;
	});
	tmp = sum(coeff * tmp);
	return(tmp);
}


###############

### Simple Log:
n = sqrt(5)
integrate(\(x) log(x^n + 1), 0, 1)
log(2) - n * int.FrU01(n, p = n)


###
n = 3
integrate(\(x) log(1 - x^n), 0, 1)
log(n) - n - integrate(\(x) 1/(1 - x) - n / (1 - x^n), 0, 1)$value
- n - digamma(1/n) - Euler

#
n = sqrt(7)
integrate(\(x) log(1 - x^n), 0, 1)
- n - digamma(1/n) - Euler


##################

### Basic Log:

### I( log(x)^s )
s = log(5)
integrate(\(x) abs(log(x))^s, 0, 1)
gamma(s+1)


### I( x^p * log(x)^s )
s = log(5); p = sqrt(7);
integrate(\(x) x^p * abs(log(x))^s, 0, 1)
gamma(s+1) / (p+1)^(s+1)


### Shifted

### I( log(x+1)^s )

integrate(\(x) log(x+1), 0, 1)
2*log(2) - 1
#
integrate(\(x) log(x+1)^2, 0, 1)
2*log(2)^2 - 2*(2*log(2) - 1)
#
integrate(\(x) log(x+1)^3, 0, 1)
2*log(2)^3 - 3*(2*log(2)^2 - 2*(2*log(2) - 1))
#
integrate(\(x) log(x+1)^4, 0, 1)
2*log(2)^4 - 4*(2*log(2)^3 - 3*(2*log(2)^2 - 2*(2*log(2) - 1)))


########################
########################

########################
### Simple Fractions ###
########################

### I( log(x) / (x + 1)^p ) on [0, Inf]

# qncubed3: Complex Analysis: Integral of log(x)/(x+1)^2
# https://www.youtube.com/watch?v=tPveHNdBWR8

# Note:
# - easy pole with higher multiplicity;


### Full generalisation: I( log(x) / (x+1)^(p+1) )
p = sqrt(3)/7
integrate(function(x) log(x) / (x+1)^(p+1), 0, Inf, rel.tol=1E-8)
- (digamma(p) + Euler)/p

### Series: x^Integer * ...;
# TODO: Genralisation of Formula;
p = 1 + sqrt(3)/7
integrate(function(x) x * log(x) / (x+1)^(p+1), 0, Inf, rel.tol=1E-8)
- (digamma(p-1) + Euler - 1) / (p*(p-1))
#
- (digamma(p-1) + Euler)/(p-1) + (digamma(p) + Euler)/p
digamma(p)/p - digamma(p-1)/(p-1) - Euler/(p*(p-1))

#
p = sqrt(5)
integrate(function(x) x^2 * log(x) / (x+1)^(p+2), 0, Inf, rel.tol=1E-8)
#
- (digamma(p-1) + Euler)/(p-1) - (digamma(p+1) + Euler)/(p+1) + 2*(digamma(p) + Euler)/(p);
- (digamma(p) + Euler) * (1/(p-1) + 1/(p+1) - 2/p) +
	+ 1/(p-1)^2 - 1/(p*(p+1))


### Variants: Examples
p = sqrt(7) - sqrt(2)/2; # p > 1 (but should be >= 1.7 for numerical reasons)
integrate(function(x) x^2 * log(x) / (x+1)^(p+2), 0, Inf, rel.tol=1E-8)
integrate(function(x) 1/9 * log(x) / (x^(1/3) + 1)^(p+2), 0, Inf, rel.tol=1E-8)
integrate(function(x) 1/4  * x^(1/2) * log(x) / (x^(1/2) + 1)^(p+2), 0, Inf, rel.tol=1E-8)
integrate(function(x) 4/25 * x^(1/5) * log(x) / (x^(2/5) + 1)^(p+2), 0, Inf, rel.tol=1E-8)
- (digamma(p) + Euler) * (1/(p-1) + 1/(p+1) - 2/p) +
	+ 1/(p-1)^2 - 1/(p*(p+1))


### Series: x^(1/2)
integrate(function(x) x^(1/2) * log(x) / (x+1)^2, 0, Inf, rel.tol=1E-8)
pi
#
integrate(function(x) x^(1/2) * log(x) / (x+1)^3, 0, Inf, rel.tol=1E-8)
0
#
integrate(function(x) x^(1/2) * log(x) / (x+1)^4, 0, Inf, rel.tol=1E-8)
- pi/24
#
integrate(function(x) x^(1/2) * log(x) / (x+1)^5, 0, Inf, rel.tol=1E-8)
- pi/24
#
integrate(function(x) x^(1/2) * log(x) / (x+1)^6, 0, Inf, rel.tol=1E-8)
- 71/(30*64) * pi


### Initial Workout:
integrate(function(x) log(x)/(x+1)^2, 0, Inf)
# == 0

###
integrate(function(x) log(x)/(x+1)^(3/2), 0, Inf, rel.tol=1E-8)
4*log(2)

###
integrate(function(x) log(x)/(x+1)^3, 0, Inf)
-1/2

###
integrate(function(x) log(x)/(x+1)^(5/2), 0, Inf, rel.tol=1E-8)
4*log(2)/3 - 4/3

###
integrate(function(x) log(x)/(x+1)^4, 0, Inf)
-1/2

###
integrate(function(x) log(x)/(x+1)^(7/2), 0, Inf)
4*log(2)/5 - 16/15

### n = 5
integrate(function(x) log(x)/(x+1)^5, 0, Inf)
- 11 / gamma(5)

###
integrate(function(x) log(x)/(x+1)^(9/2), 0, Inf)
4*log(2)/7 - 92/(3*5*7)

###
integrate(function(x) log(x)/(x+1)^(11/2), 0, Inf, rel.tol=1E-8)
4*log(2)/9 - 704/(3*5*7*9)

### n = 6
integrate(function(x) log(x)/(x+1)^6, 0, Inf)
- 50 / gamma(6)

# x = exp(1i*pi);
# - 4i*pi*I + 4*pi^2/(n-1) = 2i*pi*(- 100/x^5 + 48*log(x)/x^5)/gamma(6)


### n = 7
integrate(function(x) log(x)/(x+1)^7, 0, Inf)
- 274 / gamma(7)

# x = exp(1i*pi);
# - 4i*pi*I + 4*pi^2/(n-1) = 2i*pi*(2*274 - 240*log(x)/x^6)/gamma(7)
# - 2i*I + 2*pi/(n-1) = 1i*(2*274 - 240*log(x)/x^6)/gamma(7)


### Pow = 1/n Series
integrate(function(x) log(x)/(x+1)^(4/3), 0, Inf)
- (digamma(1/3) + Euler)*3

###
integrate(function(x) log(x)/(x+1)^(6/5), 0, Inf, rel.tol=1E-8)
- (digamma(1/5) + Euler)*5


########################
########################

### I( log(x) / (x + 1)^p ) on [0, 1]

# TODO: truncated Polylog function;

###
integrate(function(x) log(x)/(x+1)^2, 0, 1)
- log(2)

###
integrate(function(x) log(x)/(x+1)^3, 0, 1)
- log(2)/2 + (1/2 - 1)/2

###
integrate(function(x) log(x)/(x+1)^4, 0, 1)
- log(2)/3 + (1/2 + 1/4/2 - 1 - 1/2)/3

###
integrate(function(x) log(x)/(x+1)^5, 0, 1)
- log(2)/4 + (1/2 + 1/4/2 + 1/8/3 - 1 - 1/2 - 1/3)/4

###
integrate(function(x) log(x)/(x+1)^6, 0, 1)
- log(2)/5 + (1/2 + 1/4/2 + 1/8/3 + 1/16/4 - 1 - 1/2 - 1/3 - 1/4)/5


### I( log(x) * sqrt(x + 1) )
integrate(\(x) log(x) * sqrt(x + 1), 0, 1)
- 2/3 * (2*log(2) - 2*asinh(1) + 2*sqrt(2) - 2 + 2/3 * (2^(3/2) - 1))
4/3 * (asinh(1) - log(2) - (1 + 2/3)*sqrt(2) + 1 + 1/3)

# Derivation:
integrate(\(x) (sqrt(x + 1) - 1)/x, 0, 1)
2*log(2) - 2*asinh(1) + 2*sqrt(2) - 2
#
integrate(\(x) ((x+1)*sqrt(x + 1) - 1)/x, 0, 1)
2*log(2) - 2*asinh(1) + 2*sqrt(2) - 2 + 2/3 * (2^(3/2) - 1)


### I( log(x) / sqrt(x + 1) )
integrate(\(x) log(x) / sqrt(x + 1), 0, 1)
4 * (asinh(1) - log(2) - sqrt(2) + 1)


### Power 1/3

### I( log(x) * (x + 1)^(1/3) )
integrate(\(x) log(x) * (x + 1)^(1/3), 0, 1)
-9/16*(2^(4/3) - 1 + 4*2^(1/3) - 4) + 9/8*log(2^(2/3) + 2^(1/3) + 1) - 9/8*log(3) +
	+ 9/8 / sin(pi/3) * atan((2^(1/3) + cos(pi/3)) / sin(pi/3)) +
	- 9/8 / sin(pi/3) * atan((1 + cos(pi/3)) / sin(pi/3));

### I( log(x) * (x + 1)^(2/3) )
integrate(\(x) log(x) * (x + 1)^(2/3), 0, 1)
integrate(\(x) - 3/5 * ((x + 1)^(5/3) - 1) / x, 0, 1)
-9*(2^(5/3)/25 + 2^(2/3)/10) + 9*(1/25 + 1/10) +
	+ 9/10 * log(2^(2/3) + 2^(1/3) + 1) - 9/10*log(3) +
	- 9/10 / sin(pi/3) * atan((2^(1/3) + cos(pi/3)) / sin(pi/3)) +
	+ 9/10 / sin(pi/3) * atan((1 + cos(pi/3)) / sin(pi/3));


### I( log(x) / (x + 1)^(1/3) )
integrate(\(x) log(x) / (x + 1)^(1/3), 0, 1)
5/2 * integrate(\(x) log(x) * (x + 1)^(2/3), 0, 1)$value +
	+ 3*2*2^(2/3)*(1/2 - 1/5) - 3*(1/2 - 1/5);
# see above for I( log(x) * (x+1)^(2/3) );


### I( log(x) / (x + 1)^(2/3) )
integrate(\(x) log(x) / (x + 1)^(2/3), 0, 1)
4 * integrate(\(x) log(x) * (x + 1)^(1/3), 0, 1)$value +
	+ 3*2*2^(1/3)*(1 - 1/4) - 3*(1 - 1/4);
# see above for I( log(x) * (x+1)^(1/3) );

# Derivations:
integrate(\(x) - 3/4 * ((x+1)^(4/3) - 1) / x, 0, 1)
integrate(\(x) - 9/4 * (x^6 - x^2) / (x^3 - 1), 1, 2^(1/3))
-9/16*(2^(4/3) - 1) + integrate(\(x) - 9/4 * (x^3 - x^2) / (x^3 - 1), 1, 2^(1/3))$value
-9/16*(2^(4/3) - 1 + 4*2^(1/3) - 4) +
	+ integrate(\(x) 9/4 * (x + 1) / (x^2 + x + 1), 1, 2^(1/3))$value
-9/16*(2^(4/3) - 1 + 4*2^(1/3) - 4) +
	+ integrate(\(x) 9/4 * (1/2*(2*x+1) + 1/2) / (x^2 + x + 1), 1, 2^(1/3))$value;
	# + 9/8 / sin(pi/3) * atan((x + cos(pi/3)) / sin(pi/3))
	# + integrate(\(x) 9/8 / (x^2 + x + 1), 1, 2^(1/3))$value

#
integrate(\(x) - 3/5 * ((x + 1)^(5/3) - 1) / x, 0, 1)
integrate(\(x) - 9/5 * (x^7 - x^2)/ (x^3 - 1), 1, 2^(1/3))
-9/25 * (2^(5/3) - 1) - 9/10*(2^(2/3) - 1) +
	+ integrate(\(x) 9/5 * x / (x^2 + x + 1), 1, 2^(1/3))$value
-9*(2^(5/3)/25 + 2^(2/3)/10) + 9*(1/25 + 1/10) +
	+ 9/10 * (log(2^(2/3) + 2^(1/3) + 1) - log(3)) +
	- 9/10 / sin(pi/3) * atan((2^(1/3) + cos(pi/3)) / sin(pi/3)) +
	+ 9/10 / sin(pi/3) * atan((1 + cos(pi/3)) / sin(pi/3));
	# + integrate(\(x) 9/5 * x / (x^2 + x + 1), 1, 2^(1/3))$value


# Derivation:
integrate(\(x) ((x+1)^(1/3) - 1)/x, 0, 1)
integrate(\(x) 3*(x^3 - x^2) / (x^3 - 1), 1, 2^(1/3))
#
integrate(\(x) ((x+1)^(4/3) - 1)/x, 0, 1)
integrate(\(x) 3*(x^6 - x^2) / (x^3 - 1), 1, 2^(1/3))

#
integrate(\(x) ((x+1)^(2/3) - 1)/x, 0, 1)
integrate(\(x) 3*(x^4 - x^2) / (x^3 - 1), 1, 2^(1/3))
#
integrate(\(x) ((x+1)^(5/3) - 1)/x, 0, 1)
integrate(\(x) 3*(x^7 - x^2) / (x^3 - 1), 1, 2^(1/3))


##################
##################

### I( log(1 - x) / x ) on [0, 1/2]
# see special values of polylog function:
# https://en.wikipedia.org/wiki/Polylogarithm

integrate(\(x) log(1-x) / x, 0, 1/2)
- pi^2/12 + log(2)^2 / 2

###
integrate(\(x) log(x) / (x^2 - x), 1, 2)
integrate(\(x) log(x + 1) / (x^2 + x), 0, 1)
pi^2/12 - log(2)^2 / 2


###
integrate(\(x) log(1 + x) / x, 0, 1)
integrate(\(x) - log(1 - x) / (2*x), 0, 1)
integrate(\(x) - log(1 - x^2) / x, 0, 1)
pi^2 / 12

#
integrate(\(x) log(1 + x) / x^(1/2), 0, 1)
2*log(2) - 4*(1 - pi/4)

#
integrate(\(x) log(1 + x) / x^(1/3), 0, 1)
integrate(\(x) 3 * x * log(1 + x^3), 0, 1)
3/2*log(2) - 9/4 + 9/2 * integrate(\(x) x / (1 + x^3), 0, 1)$value
3/2*log(2) - 9/4 + 9/2 * int.FrU01(3, p=1)

### Generalization

### p < 1
p = 1/5; # but can be negative;
integrate(\(x) log(1 + x) / x^p, 0, 1)
log(2)/(1 - p) - 1/(p*(1-p)) * int.FrU01(1/p, p = 2/p - 2)

###
p = 1/sqrt(5)
integrate(\(x) log(1 + x) / x^p, 0, 1)
log(2)/(1 - p) - 1/(p*(1-p)) * int.FrU01(1/p, p = 2/p - 2)

# Derivation:
integrate(\(x) 1/p * x^(1/p - 2) * log(1 + x^(1/p)), 0, 1)
1/(1 - p)*log(2) - 1/(1-p)/p * integrate(\(x) x^(2/p - 2) / (1 + x^(1/p)), 0, 1)$value
1/(1 - p)*log(2) - 1/(1-p)/p * int.FrU01(1/p, p = 2/p - 2)


### I( log(x^p + 1) / x^p )
integrate(\(x) log(x^2 + 1) / x^2, 0, 1)
- log(2) + integrate(\(x) 2 / (x^2 + 1), 0, 1)$value
pi/2 - log(2)

###
integrate(\(x) log(x^3 + 1) / x^3, 0, 1)
- log(2)/2 + integrate(\(x) 3/2 / (x^3 + 1), 0, 1)$value
- log(2)/2 + 3/2 * int.FrU01(3, p = 0)

### I( log(x^p + 1) / x^p )
p = sqrt(5)
integrate(\(x) log(x^p + 1) / x^p, 0, 1)
- log(2)/(p-1) + p/(p-1) * int.FrU01(p, p = 0)


### Other:
integrate(\(x) log( (1 + x^2)/(1 - x^2) ) / x, 0, 1)
integrate(\(x) log( (1 + x)/(1 - x) ) / (2*x), 0, 1)
integrate(\(x) log( (1 + x)/(x - 1) ) / (2*x), 1, Inf)
pi^2 / 8


###
integrate(\(x) log(x) / (x - 1), 0, 1/2)
pi^2 / 12 + log(2)^2 / 2
integrate(\(x) log(x) / (x - 1), 1/2, 1)
pi^2 / 12 - log(2)^2 / 2

# Derivation:
integrate(\(x) log(1 + x) / x, 0, 1)
integrate(\(x) log(x) / (x^2 - x), 1/2, 1)
pi^2 / 12
# =>
integrate(\(x) log(x) / (x - 1), 1/2, 1)
pi^2 / 12 - log(2)^2 / 2


######################
######################


### I( log(x)^s / (x^n - 1) ) on [0, 1]
n = sqrt(5); s = 3; # s = integer;
integrate(\(x) log(x)^s / (x^n - 1), 0, 1)
pracma::psi(s, 1/n) / n^(s+1)

### I( log(x)^s / (x^n + 1) ) on [0, 1]
n = sqrt(5); s = 3; # s = integer;
integrate(\(x) log(x)^s / (x^n + 1), 0, 1)
pracma::psi(s, 1/n) / n^(s+1) - 2*pracma::psi(s, 1/(2*n)) / (2*n)^(s+1)


### I( log(x) / (x^n - 1) ) on [0, 1]
n = sqrt(5)
integrate(\(x) log(x)/(x^n - 1), 0, 1)
pracma::psi(1, 1/n) / n^2


### Examples:
integrate(\(x) log(x)/(x - 1), 0, 1)
pi^2/6

### Gen: I( log(x)^s / (x - 1) )
s = 5
integrate(\(x) log(x)^s / (x - 1), 0, 1)
gamma(s+1) * zeta(s+1)

#
s = sqrt(5)
integrate(\(x) - abs(log(x))^s / (x - 1), 0, 1)
gamma(s+1) * zeta(s+1)


### Gen: I( log(x)^s / (x + 1) )
s = sqrt(5)
integrate(\(x) - abs(log(x))^s / (x + 1), 0, 1)
- gamma(s+1) * zeta(s+1) * (1 - 1/2^s)

#
s = sqrt(3)
integrate(\(x) - abs(log(x))^s / (x^2 - 1), 0, 1)
gamma(s+1) * zeta(s+1) * (1 - 1/2^(s+1))


### Derived Variants:
integrate(\(x) log(x)/(x^3 - 1), 0, 1)
pracma::psi(1, 1/3)/9
# =>
integrate(\(x) log(x) * (x+2) / (x^2 + x + 1), 0, 1)
pi^2/6 - pracma::psi(1, 1/3)/3

###
integrate(\(x) log(x) / (x^2 + x + 1), 0, 1)
- (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 9

### =>
integrate(\(x) log(x) * x / (x^2 + x + 1), 0, 1)
pi^2/6 - (pracma::psi(1, 1/3) + 2*pracma::psi(1, 2/3))/9
pi^2/(2*27) - pracma::psi(1, 2/3)/9


### I( x^p * log(x) / (x^n + 1) ) on [0, 1]
n = sqrt(5)
integrate(\(x) log(x)/(x^n + 1), 0, 1)
pracma::psi(1, 1/n) / n^2 - pracma::psi(1, 1/(2*n)) / (2*n^2)

### Full:
n = sqrt(5); p = sqrt(3);
integrate(\(x) x^p * log(x) / (x^n + 1), 0, 1)
(pracma::psi(1, (p+1)/n) - pracma::psi(1, (p+1)/(2*n))/2) / n^2

