########################
###
### Leonard Mada
### [the one and only]
###
### Integrals: Logarithms
### Simple Fractions
###
### draft v.0.2f



### History

# - [refactor] moved Radicals to new file:
#   Integrals.Log.Fractions.Radicals.R;


####################

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
# I( x^p / (x^n + 1)^pow ) on [0, Inf]
int.FrUInf = function(n, p=0, pow=1, coeff=1) {
	if(length(n) != 1) stop("Only numerator can be polynomial!");
	k = 1/pow;
	tmp = sapply(p, function(p) {
		gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) / n;
	});
	tmp = sum(coeff * tmp);
	return(tmp);
}


###############

### Simple Log:

### I( x^p * log(x^n + 1) )
n = sqrt(5)
integrate(\(x) log(x^n + 1), 0, 1)
log(2) - n * int.FrU01(n, p = n)

###
n = sqrt(11);
p = sqrt(3);
integrate(\(x) x^p * log(x^n + 1), 0, 1)
- (digamma((p+1)/(2*n) + 1) - digamma((p+1)/n + 1)) / (p+1);

### Special Case: p = -1
n = sqrt(11);
integrate(\(x) log(x^n + 1) / x, 0, 1)
pi^2 / (12*n)


### I( x^p * log(1 - x^n) )
n = 3
integrate(\(x) log(1 - x^n), 0, 1)
log(n) - n - integrate(\(x) 1/(1 - x) - n / (1 - x^n), 0, 1)$value
- n - digamma(1/n) - Euler

#
n = sqrt(7)
integrate(\(x) log(1 - x^n), 0, 1)
- n - digamma(1/n) - Euler

###
n = sqrt(11);
p = sqrt(3);
integrate(\(x) x^p * log(1 - x^n), 0, 1)
- (digamma((n+p+1)/n) + Euler) / (p+1);


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


###
s = 7
integrate(\(x) log(x+1)^s, 0, 1)
integrate(\(x) x^s * exp(x), 0, log(2))

### All:
s = sqrt(7)
integrate(\(x) log(x+1)^s, 0, 1)
integrate(\(x) x^s * exp(x), 0, log(2))

### All:
s = sqrt(7); up = sqrt(5);
integrate(\(x) log(x+1)^s, 0, up)
integrate(\(x) x^s * exp(x), 0, log(up + 1))

# TODO:
# - probably based on exponentials;


####################

### I( log(1 - x^n)^s / x )

### Integer:
s = 2
integrate(\(x) log(1 - x)^s / x, 0, 1, rel.tol=1E-8)
(-1)^s * gamma(s+1) * pracma::zeta(s+1)

###
s = 3
integrate(\(x) log(1 - x)^s / x, 0, 1, rel.tol=1E-8)
(-1)^s * gamma(s+1) * pracma::zeta(s+1)

### Gen 1:
s = sqrt(5)
integrate(\(x) abs(log(1 - x))^s / x, 0, 1, rel.tol=1E-8)
gamma(s+1) * pracma::zeta(s+1)

### Gen 2:
s = sqrt(5); n = sqrt(7)
integrate(\(x) abs(log(1 - x^n))^s / x, 0, 1, rel.tol=1E-8)
gamma(s+1) * pracma::zeta(s+1) / n


### I( log(1 - x)^s / x^2 )

### x-Pow = 2
s = sqrt(7)
integrate(\(x) abs(log(1 - x))^s / x^2, 0, 1, rel.tol=1E-8)
gamma(s+1) * pracma::zeta(s)


### x-Pow = 3
s = 3;
integrate(\(x) log(1 - x)^3 / x^3, 0, 1, rel.tol=1E-8)
- gamma(4) * (pracma::zeta(3) + pracma::zeta(2)) / 2

###
integrate(\(x) log(1 - x)^4 / x^3, 0, 1, rel.tol=1E-8)
gamma(5) * (pracma::zeta(4) + pracma::zeta(3)) / 2

### Gen: I( log(1 - x)^s / x^3 )
s = sqrt(11)
integrate(\(x) abs(log(1 - x))^s / x^3, 0, 1, rel.tol=1E-8)
gamma(s+1) * (pracma::zeta(s) + pracma::zeta(s-1)) / 2


### x-Pow = 4
s = sqrt(13)
integrate(\(x) abs(log(1 - x))^s / x^4, 0, 1, rel.tol=1E-8)
gamma(s+1) * (pracma::zeta(s) + 3/2*pracma::zeta(s-1) + pracma::zeta(s-2)/2) / 3

# Derivation:
integrate(\(x) s/3 * abs(log(1 - x))^(s-1) / x^3 / (1-x) +
	- s/3 * abs(log(1 - x))^(s-1)/(1-x), 0, 1, rel.tol=1E-8)
integrate(\(x) s/3 * abs(log(1 - x))^(s-1) * (1/x^3 + 1/x^2 + 1/x + 1/(1-x)) +
	- s/3 * abs(log(1 - x))^(s-1)/(1-x), 0, 1, rel.tol=1E-8)
integrate(\(x) s/3 * abs(log(1 - x))^(s-1) * (1/x^3 + 1/x^2 + 1/x), 0, 1, rel.tol=1E-8)


########################
########################

########################
### Simple Fractions ###
########################

### I( x^p * log(x) / (x^n + 1)^k )

# - started initially from a much simpler example;
# - full generalisation is below;


### I( log(x) / (x + 1)^k ) on [0, Inf]

# qncubed3: Complex Analysis: Integral of log(x)/(x+1)^2
# https://www.youtube.com/watch?v=tPveHNdBWR8

# Note:
# - easy pole with higher multiplicity;


### Simple Gen: I( log(x) / (x + 1)^k )
k = sqrt(3)/7
integrate(function(x) log(x) / (x+1)^(k+1), 0, Inf, rel.tol=1E-8)
- (digamma(k) + Euler)/k


### Full Gen: I( x^p * log(x) / (x^n + 1)^k )
p = sqrt(3); n = sqrt(7); k = sqrt(5);
integrate(function(x) x^p * log(x) / (x^n + 1)^k, 0, Inf, rel.tol=1E-8)
gamma((p+1)/n) * gamma(k - (p+1)/n) *
	(digamma((p+1)/n) - digamma(k - (p+1)/n)) / gamma(k) / n^2;


### [old]

### Series: x^Integer * ...;
# [DONE] Generalisation of Formula;

### Special Case: n = 1
# - using 1/k => k;
p = 1/sqrt(3); k = sqrt(5);
integrate(function(x) x^p * log(x) / (x + 1)^k, 0, Inf, rel.tol=1E-8)
gamma(p+1) * gamma(k - (p+1)) * (digamma((p+1)) - digamma(k - (p+1))) / gamma(k);


# [old]
p = 1 + sqrt(3)/7
integrate(function(x) x * log(x) / (x+1)^(p+1), 0, Inf, rel.tol=1E-8)
- (digamma(p-1) + Euler - 1) / (p*(p-1))
#
- (digamma(p-1) + Euler)/(p-1) + (digamma(p) + Euler)/p
digamma(p)/p - digamma(p-1)/(p-1) - Euler/(p*(p-1))

#
p = sqrt(5)
integrate(function(x) x^2 * log(x) / (x+1)^(p+2), 0, Inf, rel.tol=1E-8)
# ???
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

### Derivation: log(x) / (x^2+1)^(5/2)
integrate(function(x) log(x) / (x^2+1)^(3/2), 0, Inf, rel.tol=1E-8)
integrate(function(x) 1/4 * x^(-1/2) * log(x) / (x+1)^(3/2), 0, Inf, rel.tol=1E-8)
- log(2)
#
integrate(function(x) (log(x) - 1)*x^2 / (x^2+1)^(5/2), 0, Inf, rel.tol=1E-8)
- log(2)/3
# =>
integrate(function(x) log(x) / (x^2+1)^(5/2), 0, Inf, rel.tol=1E-8)
-2/3 * log(2) - int.FrUInf(2,0,3/2) + int.FrUInf(2,0,5/2)

### I( log(x) / (x^2+1)^(7/2) )
integrate(function(x) log(x) / (x^2+1)^(7/2), 0, Inf, rel.tol=1E-8)
-2/3*4/5 * log(2) - (4*int.FrUInf(2,0,3/2) + int.FrUInf(2,0,5/2)) / 5 +
	+ int.FrUInf(2,0,7/2)


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

###
n = 1/sqrt(7)
integrate(function(x) log(x)/(x+1)^(n+1), 0, Inf, rel.tol=1E-8)
- (digamma(n) + Euler) / n;


########################
########################

#################
### on [0, 1] ###

### Gen: Simple
p = - 1/3; q = 1/5; b = 4;
integrate(function(x) x^p * (b - x)^q * log(x), lower=0, upper=b)
b^(p+q+1) * gamma(p+1)*gamma(q+1) / gamma(p+q+2) *
	(digamma(p+1) - digamma(p+q+2) + log(b));

### Base:
p = 1/3; q = 1/5; b = 4;
integrate(function(x) x^p * (b - x)^q, lower=0, upper=b)
b^(p+q+1) * gamma(p+1)*gamma(q+1) / gamma(p+q+2)


### True [0,1]-Series

### I( x * log(x+1) / (x+1)^s )
s = 1/sqrt(3)
integrate(\(x) x * log(x+1) / (x+1)^s, 0, 1)
2^(2-s)* log(2)/(2-s) - 2^(1-s)*log(2)/(1-s) +
	- (2^(2-s)-1)/(2-s)^2 + (2^(1-s)-1)/(1-s)^2;


### I( x^2 * log(x+1) / (x+1)^s )
s = 1/sqrt(3)
integrate(\(x) x^2 * log(x+1) / (x+1)^s, 0, 1)
(4/(3-s) - 4/(2-s) + 1/(1-s)) * 2^(1-s)*log(2) +
	- (2^(3-s)-1)/(3-s)^2 + (2^(3-s)-2)/(2-s)^2 - (2^(1-s)-1)/(1-s)^2;


### Radicals: I( log(x) * (x+1)^s )
# - have been moved to file:
#   Integrals.Log.Fractions.Radicals.R
# - moved also: Pow = Integer;


##################
##################

# Note:
# - in-extenso coverage:
#   see section further below;

### I( x^p * log(x) / (x - 1) )
p = sqrt(3)
integrate(\(x) x^p * log(x) / (x - 1), 0, 1)
pracma::psi(1, p+1)


### I( x^p * log(x)^2 / (x - 1) )
p = sqrt(3)
integrate(\(x) x^p * log(x)^2 / (x - 1), 0, 1)
pracma::psi(2, p+1)


### I( x^p * log(x) / (x + 1) )
p = sqrt(3)
integrate(\(x) x^p * log(x) / (x + 1), 0, 1)
pracma::psi(1, p+1) - pracma::psi(1, (p+1)/2) / 2


### I( x^p * log(x)^2 / (x + 1) )
p = sqrt(3)
integrate(\(x) x^p * log(x)^2 / (x + 1), 0, 1)
pracma::psi(2, p+1) - pracma::psi(2, (p+1)/2) / 4

#
p = sqrt(3)
integrate(\(x) 4*x^(2*p+1) * log(x) / (x^2 - 1), 0, 1)
integrate(\(x) 2*x^(2*p+1) * log(x) * (1/(x - 1) - 1/(x+1)), 0, 1)
pracma::psi(1, p+1)


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


### Other:
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

### p < 1 (p < 2)
p = 1/5; # but can be negative;
integrate(\(x) log(1 + x) / x^p, 0, 1)
log(2)/(1 - p) - 1/(p*(1-p)) * int.FrU01(1/p, p = 2/p - 2)
# for p = 1: pi^2 / 12;

###
p = 1/sqrt(5)
integrate(\(x) log(1 + x) / x^p, 0, 1)
log(2)/(1 - p) - 1/(p*(1-p)) * int.FrU01(1/p, p = 2/p - 2)

# Derivation:
integrate(\(x) 1/p * x^(1/p - 2) * log(1 + x^(1/p)), 0, 1)
1/(1 - p)*log(2) - 1/(1-p)/p * integrate(\(x) x^(2/p - 2) / (1 + x^(1/p)), 0, 1)$value
1/(1 - p)*log(2) - 1/(1-p)/p * int.FrU01(1/p, p = 2/p - 2)


### I( log(x^p + 1) / x^p )

###
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


### I( log(1 - x^n) / x^p )

###
p = sqrt(5)
integrate(\(x) log(1 - x^p) / x^p, 0, 1)
(digamma(1/p) + Euler) / (p-1)


### I( x^p * log(1 - x^n) )
n = sqrt(5); p = sqrt(3);
integrate(\(x) x^p * log(1 - x^n), 0, 1)
- (digamma((n+p+1)/n) + Euler) / (p+1)

### Special Cases:
integrate(\(x) log(1 - x) / x, 0, 1)
- pi^2 / 6
# Note: Dp(digamma(1 + (p+1)/n)) => 1/n;
n = sqrt(5)
integrate(\(x) log(1 - x^n) / x, 0, 1)
- pi^2 / (6*n);


### Other:
integrate(\(x) log( (1 + x^2)/(1 - x^2) ) / x, 0, 1)
integrate(\(x) log( (1 + x)/(1 - x) ) / (2*x), 0, 1)
integrate(\(x) log( (1 + x)/(x - 1) ) / (2*x), 1, Inf)
pi^2 / 8

###
integrate(\(x) log( (1 + x^2)/(1 - x^2) ) / x^2, 0, 1)
pi/2 + log(2)
# based on:
integrate(\(x) log(1 + x^2)/x^2, 0, 1)
pi/2 - log(2)
integrate(\(x) log(1 - x^2)/x^2, 0, 1)
- 2*log(2)


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

##################
### Elementary ###

### I( x^p * log(x)^s / (x^n - 1) )
### I( x^p * log(x)^s / (x^n + 1) )
# - on [0, 1];
# - Note: pracma::psi handles only integer orders;

###
s = 3; # s = integer;
n = sqrt(5); p = sqrt(7);
integrate(\(x) x^p * log(x)^s / (x^n - 1), 0, 1)
pracma::psi(s, (p+1)/n) / n^(s+1)

###
s = 3; # s = integer;
n = sqrt(5); p = sqrt(7);
integrate(\(x) x^p * log(x)^s / (x^n + 1), 0, 1)
pracma::psi(s, (p+1)/n) / n^(s+1) - 2*pracma::psi(s, (p+1)/(2*n)) / (2*n)^(s+1)


### Simpler Versions:

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
integrate(\(x) log(x) / (x^n - 1), 0, 1)
pracma::psi(1, 1/n) / n^2


### Examples:
integrate(\(x) log(x) / (x - 1), 0, 1)
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

# p = 0:
n = sqrt(5)
integrate(\(x) log(x) / (x^n + 1), 0, 1)
pracma::psi(1, 1/n) / n^2 - pracma::psi(1, 1/(2*n)) / (2*n^2)

### Full:
n = sqrt(5); p = sqrt(3);
integrate(\(x) x^p * log(x) / (x^n + 1), 0, 1)
(pracma::psi(1, (p+1)/n) - pracma::psi(1, (p+1)/(2*n))/2) / n^2

