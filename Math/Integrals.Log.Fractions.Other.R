########################
##
## Leonard Mada
## [the one and only]
##
## Integrals: Logarithms
## Log-Fractions: Other
##
## draft v.0.3h


##################
### Logarithms ###
##################

# - definite integrals;
# - various types of Logarithms combined with fractions;


### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;
# Note:
# Catalan = - I(log(x)/(x^2 + 1), lower=0, upper=1)


### History:

# [refactor] Integrals of type log(...) * log(...)
# moved to new file: Integrals.Log.Prod.R;
# [refactor] Integrals of type atan(...) * log(...)
# moved to new file: Integrals.Trig.Tan.Atan.R;
# [refactor] Integrals of type log(...) / (x^4 + 1)
# moved to new file: Integrals.Log.Fractions.P4.R;


###################
###################

### I( log(1-x)^2 / x )
integrate(\(x) log(1-x)^2 / x, 0, 1, rel.tol=1E-8)
2 * pracma::zeta(3)

### I( log(1-x)^3 / x )
integrate(\(x) log(1-x)^3 / x, 0, 1, rel.tol=1E-8)
gamma(4) * pracma::zeta(4)


### Gen: I( log(1-x)^p / x )
p = sqrt(5)
integrate(\(x) abs(log(1-x))^p / x, 0, 1, rel.tol=1E-8)
gamma(p+1) * pracma::zeta(p+1)


### Other Intervals

### on [0, 1/2]
integrate(\(x) log(1 - x) / x, 0, 1/2)
- pi^2/12 + log(2)^2/2

### on [0, 1/3]
integrate(\(x) log(1 - x) / x, 0, 1/3)
- pracma::polylog(1/3, 2)

# TODO: closed formula;


####################

### I( |log( (1-x) / (1+x) )|^n ) on [0, 1]
# Maths 505: This integral is actually one of your favorite constants
# https://www.youtube.com/watch?v=83mUOaF7G9A


### Gen:
n = sqrt(3)
integrate(\(x) abs(log((1-x)/(1+x)))^n, 0, 1)
integrate(\(x) log((1+x)/(1-x))^n, 0, 1)
2*gamma(n+1)*(1 - 2^(1-n)) * pracma::zeta(n)


### I( x * log((1+x)/(1-x))^2 / (1-x^2) )
integrate(\(x) x * log((1+x)/(1-x))^2 / (1-x^2) - 1/2 * log((1-x)/2)^2/(1-x), 0, 1)
- 1/3 * gamma(4)*(1 - 1/4) * pracma::zeta(3) + log(2)^3 / 6

# library(Rmpfr)
x = integrate(\(x) {
	x = mpfr(x, 240);
	y = x * log((1+x)/(1-x))^(n-1) / (1-x^2) - 1/2*log((1-x)/2)^2/(1-x);
	as.numeric(y); }, 0, 1, rel.tol=1E-8)
x$value + 1/3 * gamma(4)*(1 - 1/4) * pracma::zeta(3) - log(2)^3 / 6


### I( log((x - 1) / (x + 1)) / x ) on [1, Inf]
# 1. Maths 505: A RIDICULOUSLY AWESOME LOG INTEGRAL!!!
#    https://www.youtube.com/watch?v=tTu6hedSlm0
# 2. Maths 505: A satisfying integral
#    https://www.youtube.com/watch?v=QXptzDZveTQ

# numerical issue:
integrate(\(x) log((exp(x) - 1) / (exp(x) + 1)), 0, Inf)
integrate(\(x) log((exp(x) - 1) / (exp(x) + 1)), 0, 100)
# robust:
integrate(\(x) log((1 - exp(-x)) / (1 + exp(-x))), 0, Inf, rel.tol=1E-9)
# alternative:
integrate(\(x) log((x - 1) / (x + 1)) / x, 1, Inf)
integrate(\(x) log((1 - x) / (x + 1)) / x, 0, 1)
- pi^2/4


### Pow = 2
integrate(\(x) log((exp(x)+1) / (exp(x)-1))^2, 0, Inf) # numerical issues
integrate(\(x) log((1+exp(-x)) / (1 - exp(-x)))^2, 0, Inf, rel.tol=1E-9)
integrate(\(x) log((1-x)/(1+x))^2 / x, 0, 1)
integrate(\(x) 2 * log(x)^2 / (1 - x^2), 0, 1)
7/2 * pracma::zeta(3)

### Pow = 3
# see Fractions.Unity: => psi(n, 1/2);
integrate(\(x) log((1-x)/(1+x))^3 / x, 0, 1)
integrate(\(x) 2 * log(x)^3 / (1 - x^2), 0, 1)
- 45/4 * pracma::zeta(4);
pracma::psi(3, 1/2) / 2^3;


### Gen: I( log((1-x)/(1+x))^n / x )
n = sqrt(5)
integrate(\(x) abs(log((1-x)/(1+x)))^n / x, 0, 1)
integrate(\(x) 2 * log(1/x)^n / (1 - x^2), 0, 1)
gamma(n+1) * (2-1/2^n) * pracma::zeta(n+1)

# Note: pracma::psi requires n = Integer;
# - may be correct only for integers anyway!
pracma::psi(n, 1/2) / 2^n;


### Simple:
integrate(\(x) log((1 - x) / (x + 1)), 0, 1)
- 2*log(2)

###
integrate(\(x) log(x) * log((1 - x) / (x + 1)), 0, 1)
2*log(2) - pi^2/12


###################
###################

# Note:
# - the polylog2 function is in file:
#   Integrals.Polylog.Helper.R;


### Gen: I( x^p * log(1-x^2) / (1+x) )
p = sqrt(5);
integrate(\(x) x^p * log(1-x^2) / (1+x), 0, 1, rel.tol=1E-13)
((pracma::psi(1, (p+2)/2) - pracma::psi(1, (p+1)/2)) +
	+ (digamma((p+1)/2) - digamma((p+2)/2)) * Euler * 2 +
	+ (digamma((p+1)/2)^2 - digamma((p+2)/2)^2)) / 4;


### I( log(x^2+x+1) / (x+1) )
integrate(\(x) log(x^2+x+1) / (x+1), 0, 1)
- polylog2(-2) - pi^2/9;


### I( log(x^2-x+1) / (x+1) )
integrate(\(x) log(x^2-x+1) / (x+1), 0, 1)
polylog2(-2) + pi^2/18 + log(2)*log(3);


### Varia:

###
integrate(\(x) log((1+x)/x) / (1 + 3*x), 0, Inf, rel.tol=1E-12)
pi^2/18 - polylog2(-2)/3;


#######################

### Types: log(x +/- 1) / (x^2+1)

# I( x * log(1+x) / (x^2 + 1) ) on [0, 1]
# I( x * log(1-x) / (x^2 + 1) ) on [0, 1]
# 1. Cipher: Integrate xln(1-x)/(1 + x^2)dx from 0 to 1
#    https://www.youtube.com/watch?v=CMRvFM2N7sw
# 2. Maths 505: ONE BALLER INTEGRAL: int(0, π/4) ln(1-tan(θ))/tan(θ)
#    https://www.youtube.com/watch?v=mqJYMoABUZs

###
integrate(\(x) log(1+x) / (x^2+1), 0, 1)
pi * log(2)/8

###
integrate(\(x) log(1-x) / (x^2+1), 0, 1)
pi*log(2)/8 - Catalan

### I( 1/x * log(1+x) / (x^2+1) )
integrate(\(x) 1/x * log(1+x) / (x^2+1), 0, 1)
integrate(\(x) log(1 + x)/x -  log(1 + x) * x / (1 + x^2), 0, 1)
integrate(\(x) pi^2/12 -  log(1 + x) * x / (1 + x^2), 0, 1)
pi^2 * 7/96 - log(2)^2 / 8;

### I( 1/x * log(1-x) / (x^2+1) )
integrate(\(x) 1/x * log(1-x) / (x^2+1), 0, 1)
- pi^2 * 11/96 - log(2)^2 / 8;

#
integrate(\(x) 1/x * log(1-x^2) / (x^2+1), 0, 1, rel.tol=1E-12)
- pi^2/24 - log(2)^2/4

###
integrate(\(x) x * log(1+x) / (x^2+1), 0, 1)
pi^2 / 96 + log(2)^2 / 8

###
integrate(\(x) x * log(1-x) / (x^2+1), 0, 1)
- 5/96 * pi^2 + log(2)^2 / 8


###########
### Pow = 3

### I( x^p * log(1-x^3) / (1-x^3) )
p = sqrt(3); n = 3;
integrate(\(x) x^p * (1-x) * log(1-x^3) / (1-x^3), 0, 1)
# TODO: ugly limits;
q = 1E-5;
(1/2*(pracma::psi(1, (p+1)/n) - pracma::psi(1, (p+2)/n)) +
	+ Euler * (digamma((p+1)/n) - digamma((p+2)/n)) +
	+ 1/2 * (digamma((p+1)/n)^2 - digamma((p+2)/n)^2)) / n +
	# TODO
	- (pracma::psi(1, (p+1)/n) * gamma((p+2)/n+q) * gamma((p+1)/n) +
		- pracma::psi(1, (p+2)/n) * gamma((p+1)/n+q) * gamma((p+2)/n)
	) / gamma((p+1)/n) / gamma((p+2)/n) / n;


### Gen: I( x^p * log(1-x^3) / (1 - x^3) )
p = sqrt(5);
integrate(\(x) x^p * log(1-x^3) / (1-x^3) - 1/3 * log(3-3*x) / (1-x), 0, 1)
(pracma::psi(1, 1) - pracma::psi(1, (p+1)/3) +
	+ (digamma((p+1)/3) - digamma(1)) * Euler * 2 +
	+ (digamma((p+1)/3)^2 - digamma(1)^2) - log(3)^2 ) / 6;

### Gen: I( x^p * log(1-x^3) / (x^2 + x + 1) )
p = sqrt(3);
integrate(\(x) x^p * log(1-x^3) * (1-x) / (1-x^3), 0, 1)
((pracma::psi(1, (p+2)/3) - pracma::psi(1, (p+1)/3)) +
	+ (digamma((p+1)/3) - digamma((p+2)/3)) * Euler * 2 +
	+ (digamma((p+1)/3)^2 - digamma((p+2)/3)^2)) / 6;

### I( log(1-x^3) / (x^2 + x + 1) )
integrate(\(x) log(1-x^3) * (1-x) / (1-x^3), 0, 1)
((pracma::psi(1, 2/3) - pracma::psi(1, 1/3)) +
	+ (digamma(1/3) - digamma(2/3)) * Euler * 2 +
	+ (digamma(1/3)^2 - digamma(2/3)^2)) / 6;

### I( x * log(1-x^3) / (x^2 + x + 1) )
integrate(\(x) x * log(1-x^3) * (1-x) / (1-x^3), 0, 1)
((pracma::psi(1, 1) - pracma::psi(1, 2/3)) +
	+ (digamma(2/3) - digamma(1)) * Euler * 2 +
	+ (digamma(2/3)^2 - digamma(1)^2)) / 6;

### I( x^2 * log(1-x^3) / (x^2 + x + 1) )
integrate(\(x) x^2 * log(1-x^3) * (1-x) / (1-x^3), 0, 1)
((pracma::psi(1, 4/3) - pracma::psi(1, 3/3)) +
	+ (digamma(3/3) - digamma(4/3)) * Euler * 2 +
	+ (digamma(3/3)^2 - digamma(4/3)^2)) / 6;


### Components:

### I( log(1-x) / (1 + x^3) )
integrate(\(x) log(1-x) / (1 + x^3), 0, 1)
- pi^2 / 18 + log(2)^2 / 6 +
	- pracma::psi(1, 1/3) / 18 + pracma::psi(1, 2/3) / 9;


### I( log(1-x) / (x^2 + x + 1) )
integrate(\(x) log(1-x) / (x^2+x+1), 0, 1)
integrate(\(x) 1/3 * log(1-x^3) / (x^2+x+1), 0, 1)$value +
	- (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 18;
((pracma::psi(1, 2/3) - pracma::psi(1, 1/3)) * 2 +
	+ (digamma(1/3) - digamma(2/3)) * Euler * 2 +
	+ (digamma(1/3)^2 - digamma(2/3)^2)) / 18;

### I( log(1-x) / (x^2 - x + 1) )
integrate(\(x) log(1-x) / (x^2-x+1), 0, 1)
- (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 6;


### I( log(1+x) / (x^2 + x + 1) )
integrate(\(x) log(1+x) / (x^2+x+1), 0, 1)
(pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 36;

### I( log(1+x) / (x^2 - x + 1) )
integrate(\(x) log(1+x) / (x^2-x+1), 0, 1)
- pracma::psi(1, 1/3) / 18 +
	+ (digamma(5/6) - digamma(1/3) + digamma(2/3) - digamma(1/6)) *
		(digamma(1) - digamma(1/3)) / 18;


### I( x * log(1-x) / ... )

### I( x * log(1-x) / (x^2 - x + 1) )
integrate(\(x) x * log(1-x) / (x^2-x+1), 0, 1)
pi^2/18 - (pracma::psi(1, 1/3) - pracma::psi(1, 3/3)) / 6;
pi^2/12 - pracma::psi(1, 1/3) / 6;

### I( x * log(1-x) / (x^2 + x + 1) )
integrate(\(x) x * log(1-x) / (x^2+x+1), 0, 1)
integrate(\(x) - x * log(1-x) / (x^2-x+1), 0, 1)$value +
+ ((pracma::psi(1, 1) - 2*pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) +
	+ (digamma(2/3) - digamma(1)) * Euler * 2 +
	+ (digamma(2/3)^2 - digamma(1)^2)) / 18;
(pracma::psi(1, 1/3) - pracma::psi(1, 2/3) +
	+ (digamma(2/3) - digamma(1)) * Euler * 2 +
	+ (digamma(2/3)^2 - digamma(1)^2) ) / 18 - pi^2 * 2/27;


### I( x * log(1+x) / (x^2 + x + 1) )
integrate(\(x) x * log(1+x) / (x^2+x+1), 0, 1)
integrate(\(x) - 1/2 * log(x^2+x+1) / (x+1), 0, 1)$value +
	log(2)*log(3)/2 - (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 72;
integrate(\(x) - 1/2 * log(1-x^3) / (x+1), 0, 1)$value +
	+ log(2)*log(3)/2 - (pi^2/12 - log(2)^2/2)/2 +
	- (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 72;
# TODO

# Solution using Li2:
# diff(Li2(c(-1/3, 1/3))) = 0.67524635646487196;
(pracma::psi(1, 1/3) + 2*pracma::psi(1, 2/3)) / 36 +
	- (pracma::polylog(1/3, 2) - pracma::polylog(-1/3, 2)) / 2;


### I( x * log(1+x) / (x^2 - x + 1) )
integrate(\(x) x * log(1+x) / (x^2-x+1), 0, 1)
integrate(\(x) - x * log(1+x) / (x^2+x+1), 0, 1)$value +
	- (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 36 +
	+ (digamma(6/6) - digamma(2/6) + digamma(5/6) - digamma(3/6)) *
		(digamma(1) - digamma(2/3)) / 9 +
	- (digamma(5/6) - digamma(2/6) + digamma(4/6) - digamma(1/6)) *
		(digamma(1) - digamma(1/3)) / 36 + pi^2 / 12;
# TODO


### I( x^2 * log(1-x) / ... )

### I( x^2 * log(1-x) / (x^2-x+1) )
integrate(\(x) x^2 * log(1-x) / (x^2-x+1), 0, 1)
pi^2/12 - pracma::psi(1, 2/3) / 6 - 1;


### Pow = 6

### I( log(1-x^6) * (1-x) / (1-x^6) )
integrate(\(x) log(1-x^6) * (1-x) / (1-x^6), 0, 1, rel.tol=1E-13)
((pracma::psi(1, 2/6) - pracma::psi(1, 1/6)) +
	+ (digamma(1/6) - digamma(2/6)) * Euler * 2 +
	+ (digamma(1/6)^2 - digamma(2/6)^2)) / 12;


### Gen: I( x^p * log(1-x^6) / (1-x^6) )
p = 1/sqrt(5);
integrate(\(x) x^p * log(1-x^6) / (1-x^6) - 1/6 * log(6-6*x)/(1-x), 0, 1)
(pracma::psi(1, 1) - pracma::psi(1, (p+1)/6) +
	+ (digamma((p+1)/6) - digamma(1)) * Euler * 2 +
	+ (digamma((p+1)/6)^2 - digamma(1)^2) - log(6)^2 ) / 12;

### Gen: I( x^p * log(1-x^6) / (1-x^3) )
p = 1/sqrt(5);
integrate(\(x) x^p * log(1-x^6) / (1-x^3) - 1/3 * log(6-6*x)/(1-x), 0, 1)
(2*pracma::psi(1, 1) - pracma::psi(1, (p+1)/6) - pracma::psi(1, (p+4)/6) +
	+ (digamma((p+1)/6) - digamma(1)) * Euler * 2 +
	+ (digamma((p+4)/6) - digamma(1)) * Euler * 2 +
	+ (digamma((p+1)/6)^2 - digamma(1)^2) - log(6)^2 +
	+ (digamma((p+4)/6)^2 - digamma(1)^2) - log(6)^2 ) / 12;


### Gen: I( x^p * log(1-x^6) / (1+x^3) )
p = sqrt(2);
integrate(\(x) x^p * log(1-x^6) / (1+x^3), 0, 1, rel.tol=1E-13)
((pracma::psi(1, (p+4)/6) - pracma::psi(1, (p+1)/6)) +
	+ (digamma((p+1)/6) - digamma((p+4)/6)) * Euler * 2 +
	+ (digamma((p+1)/6)^2 - digamma((p+4)/6)^2)) / 12;

### I( log(1-x^6) / (1+x^3) )
integrate(\(x) log(1-x^6) / (1+x^3), 0, 1, rel.tol=1E-13)
((pracma::psi(1, 4/6) - pracma::psi(1, 1/6)) +
	+ (digamma(1/6) - digamma(4/6)) * Euler * 2 +
	+ (digamma(1/6)^2 - digamma(4/6)^2)) / 12;

### I( x * log(1-x^6) / (1+x^3) )
integrate(\(x) x * log(1-x^6) / (1+x^3), 0, 1, rel.tol=1E-13)
((pracma::psi(1, 5/6) - pracma::psi(1, 2/6)) +
	+ (digamma(2/6) - digamma(5/6)) * Euler * 2 +
	+ (digamma(2/6)^2 - digamma(5/6)^2)) / 12;


### Gen: I( x^p * log(1-x^6) / (1+x) )
p = sqrt(3);
integrate(\(x) x^p * log(1-x^6) / (1+x), 0, 1, rel.tol=1E-13)
( (pracma::psi(1, (p+4)/6) - pracma::psi(1, (p+1)/6)) +
- (pracma::psi(1, (p+5)/6) - pracma::psi(1, (p+2)/6)) +
+ (pracma::psi(1, (p+6)/6) - pracma::psi(1, (p+3)/6)) +
	+ (digamma((p+1)/6) - digamma((p+4)/6)) * Euler * 2 +
	- (digamma((p+2)/6) - digamma((p+5)/6)) * Euler * 2 +
	+ (digamma((p+3)/6) - digamma((p+6)/6)) * Euler * 2 +
	+ (digamma((p+1)/6)^2 - digamma((p+4)/6)^2) +
	- (digamma((p+2)/6)^2 - digamma((p+5)/6)^2) +
	+ (digamma((p+3)/6)^2 - digamma((p+6)/6)^2) ) / 12;


### I( x^p * log(1+x^3) / (1-x^3) )
p = 1/sqrt(5);
integrate(\(x) x^p * log(1+x^3) / (1-x^3) - 1/3 * log(2)/(1-x), 0, 1)
(2*pracma::psi(1, (p+1)/3) - pracma::psi(1, (p+1)/6) - pracma::psi(1, (p+4)/6) +
	+ (digamma((p+1)/6) - digamma(1)) * Euler * 2 +
	+ (digamma((p+4)/6) - digamma(1)) * Euler * 2 +
	- (digamma((p+1)/3) - digamma(1)) * Euler * 4 +
	+ (digamma((p+1)/6)^2 - digamma(1)^2) +
	+ (digamma((p+4)/6)^2 - digamma(1)^2) + 2*log(3)^2 - 2*log(6)^2 +
	- (digamma((p+1)/3)^2 - digamma(1)^2) * 2 ) / 12;

### Special Case:

### I( x^2 * log(1+x^3) / (1-x^3) )
integrate(\(x) x^2 * log(1+x^3) / (1-x^3) - 1/3 * log(2)/(1-x), 0, 1)
sum(digamma(c(2,4,6,6)/6) * c(-1,-1,1,1)) * c(digamma(1/2) + Euler) / 18 +
	- pi^2/36 + log(2)^2/6;
log(3) * c(digamma(1/2) + Euler) / 6 - pi^2/36 + log(2)^2/6;


###
integrate(\(x) log(1 + x^3) / (x+1), 0, 1)
integrate(\(x) log(1 - x^6) * (1-x^3)*(x^2-x+1) / (1-x^6), 0, 1)$value +
integrate(\(x) - log(1 - x^3) * (1-x^3)*(x^2-x+1) / (1-x^6), 0, 1)$value


### I( log(1+x) / (1 - x^3) )

### I( log(1+x) / (1 - x^3) )
integrate(\(x) log(1+x) / (1-x^3) - 1/3 * log(2)/(1-x), 0, 1)
# TODO

### I( x * log(1+x) / (1 - x^3) )
integrate(\(x) x * log(1+x) / (1-x^3) - 1/3 * log(2)/(1-x), 0, 1)
integrate(\(x) log(1+x)/(1-x^3) - 1/3 * log(2)/(1-x), 0, 1)$value +
	- (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 36;
# TODO

### I( x^2 * log(1+x) / (1 - x^3) )
integrate(\(x) x^2 * log(1+x) / (1-x^3) - 1/3 * log(2)/(1-x), 0, 1)
integrate(\(x) -2 * log(1+x)/(1-x^3) + 2/3 * log(2)/(1-x), 0, 1)$value +
	+ (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 36 - pi^2/12 + log(2)^2/2;
# TODO

# Other Relations:
integrate(\(x) x^2 * log(1+x) / (1-x^3) - 1/3 * log(2)/(1-x), 0, 1)
integrate(\(x) -2*x * log(1+x) / (1-x^3) + 2/3 * log(2)/(1-x), 0, 1)$value +
	- (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 36 +
	- pi^2/12 + log(2)^2/2;
# using Li2:
- (pi^2/9 + log(3) * log(2) + (pi^2/6 - log(2)^2) / 2 + polylog2(-2)) / 3;
# TODO

# [Variant]
integrate(\(x) x^2 * log(1+x) / (1 - x^3) - log(2)/3/(1-x), 0, 1)
integrate(\(x) 1/3 * log(1-x^3) / (1+x), 0, 1)$value - log(2)*log(3)/3;
# [Variant]
integrate(\(x) x^2 * (log(1+x) - log(2)) / (1 - x^3), 0, 1)
integrate(\(x) 1/3 * log(1-x^3) / (1+x), 0, 1, rel.tol=1E-13)


### Composite: 
# I( (1-x) * log(1+x) / (1 - x^3) )
integrate(\(x) (1-x) * log(1+x) / (1-x^3), 0, 1)
(pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 36;


# polylog2: see file Integrals.Polylog.Helper.R; (but NOT complex)


### I( log(1-x) / (x^3 + 1) )
# - copied also above;
integrate(\(x) log(1-x) / (x^3 + 1), 0, 1)
- pi^2 / 18 + log(2)^2 / 6 +
	- pracma::psi(1, 1/3) / 18 + pracma::psi(1, 2/3) / 9;

### I( x * log(1-x) / (x^3 + 1) )
integrate(\(x) x * log(1-x) / (x^3 + 1), 0, 1)
(pracma::psi(1, 2/3) - 2*pracma::psi(1, 1/3)) / 18 +
	+ pi^2 / 18 - log(2)^2 / 6;

### I( x^2 * log(1-x) / (x^3 + 1) )
integrate(\(x) x^2 * log(1-x) / (x^3 + 1), 0, 1)
- (pracma::psi(1, 1/3) + pracma::psi(1, 2/3)) / 18 +
	+ pi^2 / 36 + log(2)^2 / 6;


### I( x * log(1-x) / (x^2 - x + 1) )
integrate(\(x) x * log(1-x) / (x^2-x+1), 0, 1)
pi^2/12 - pracma::psi(1, 1/3) / 6;


### LOG(1+x) / (1 + x^3)

### I( log(1+x) / (x^3 + 1) )
integrate(\(x) log(1+x) / (x^3 + 1), 0, 1)
integrate(\(x) (log(1+x) - log(2)) / (1 - x^3), 0, 1)$value +
(pracma::psi(1, 2/3) - pracma::psi(1, 1/6)) / 36 +
	+ pracma::psi(1, 1/3) / 9 +
	- (digamma(5/6) - digamma(1/3)) *
		(digamma(2/3) - digamma(1/3)) / 18 +
	+ (digamma(2/3) - digamma(1/6)) *
		(digamma(1) - digamma(1/3)) / 18;
# TODO

### I( x * log(1+x) / (x^3 + 1) )
integrate(\(x) x * log(1+x) / (x^3 + 1), 0, 1)
integrate(\(x) (log(1+x) - log(2)) / (x^3 - 1), 0, 1)$value +
(pracma::psi(1, 5/6) - pracma::psi(1, 2/3)) / 36 +
	- pracma::psi(1, 2/3) / 9 +
	- (digamma(2/3) - digamma(1/3)) * log(2) / 9 +
	+ (digamma(5/6) - digamma(1/3)) *
		(digamma(1) - digamma(1/3)) / 18;
# TODO

### I( x^2 * log(1+x) / (x^3 + 1) )
integrate(\(x) x^2 * log(1+x) / (x^3 + 1), 0, 1)
integrate(\(x) (x^2*log(1+x) - log(2)) / (1 - x^3), 0, 1)$value +
	+ (pracma::psi(1, 1) - pracma::psi(1, 1/2)) / 36 +
	+ pracma::psi(1, 1) / 9 +
	+ (digamma(1) - digamma(1/3)) * log(2) / 9 +
	+ (digamma(2/3) - digamma(1/6)) *
		(digamma(2/3) - digamma(1/3)) / 18;
# => I( log(1+x) / (1-x^3) )
integrate(\(x) -2 * log(1+x)/(1-x^3) + 2/3 * log(2)/(1-x), 0, 1)$value +
	+ (5*pracma::psi(1, 1) - pracma::psi(1, 1/2)) / 36 +
	+ (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 36 +
	+ (digamma(2/3) - digamma(1/6)) *
		(digamma(2/3) - digamma(1/3)) / 18 +
	- (digamma(1) - digamma(1/3)) * log(2) * 2/9 +
	- pi^2/12 + log(2)^2/2 + log(2)*log(3) / 3;
# TODO

### Transform:
integrate(\(x) (log(1+x) - log(2)) / (x^3 - 1), 0, 1)
integrate(\(x) log(1+x)/(x^3 - 1) - log(2)/3/(x-1), 0, 1)$value +
	+ (digamma(1) - digamma(1/3) - log(3)) * log(2) / 3;
# TODO


### I( log(1+x) / x / (x^3 + 1) )
# on [0, Inf]
integrate(\(x) log(1+x) / x / (x^3 + 1), 0, Inf)
pi^2 * 5 / (2*27);


### I( log(1+x) / x / (x^2 + x + 1) )
integrate(\(x) log(x+1) / x / (x^2+x+1), 0, 1)
integrate(\(x) x * log(1+x) / (x^2-x+1), 0, 1)$value +
- (digamma(6/6) - digamma(2/6) + digamma(5/6) - digamma(3/6)) *
	(digamma(1) - digamma(2/3)) / 9 +
+ (digamma(5/6) - digamma(2/6) + digamma(4/6) - digamma(1/6)) *
	(digamma(1) - digamma(1/3)) / 36;
# TODO

### Alternative:
integrate(\(x) log(1+x) / x / (x^2+x+1), 0, 1)
integrate(\(x) - x * log(1+x) / (x^2+x+1), 0, 1)$value +
	+ pi^2 / 12 - (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 36;
# Variant:
integrate(\(x) -1/2 * log(x^2-x+1) / (x+1), 0, 1)$value +
	+ (pracma::psi(1, 1/3) / 2 + pracma::psi(1, 2/3)) / 18;
# TODO


### I( (x-x^2) * log(1+x) / (1 - x^3) )
integrate(\(x) x * log(1+x) / (x^2+x+1), 0, 1)
integrate(\(x) (x-x^2) * log(1+x) / (1-x^3), 0, 1)
integrate(\(x) log(x+1)/x - log(x+1)/(x^2+x+1) - log(x+1)/x/(x^2+x+1), 0, 1)
integrate(\(x) - log(x+1)/x/(x^2+x+1), 0, 1)$value +
	+ pi^2 / 12 - (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 36;
# TODO


### I( log(1-x) / x / (x^3 + 1) )
# on [0, Inf]
integrate(\(x) log(abs(1-x)) / x / (x^3 + 1), 0, 1)$value +
integrate(\(x) log(abs(1-x)) / x / (x^3 + 1), 1, Inf)$value;
- pi^2 * 17 / (4*27);


### I( log(x^2+x+1) / (1-x) )
integrate(\(x) (log(x^2+x+1) - log(3)) / (1-x), 0, 1)
integrate(\(x) 2 * x * log(1-x) / (x^2+x+1), 0, 1, rel.tol=1E-13)$value +
	+ ((pracma::psi(1, 2/3) - pracma::psi(1, 1/3)) * 2 +
	+ (digamma(1/3) - digamma(2/3)) * Euler * 2 +
	+ (digamma(1/3)^2 - digamma(2/3)^2)) / 18;
(digamma(1/3) + digamma(2/3) - 2*digamma(1)) * Euler / 9 +
	+ (digamma(1/3)^2 + digamma(2/3)^2 - 2*digamma(1)^2) / 18 +
	- pi^2 * 4/27;

### I( log(x^2+x+1) / (1+x) )
integrate(\(x) log(x^2+x+1) / (1+x), 0, 1)
integrate(\(x) -2 * x * log(1+x) / (x^2+x+1), 0, 1)$value +
	+ (pracma::psi(1, 2/3) - pracma::psi(1, 1/3)) / 36 + log(2)*log(3);
# TODO

# Solution using Li2:
- (pracma::psi(1, 1/3) + pracma::psi(1, 2/3)) / 12 + # - pi^2 / 9
	+ (pracma::polylog(1/3, 2) - pracma::polylog(-1/3, 2)) + log(2)*log(3);


### I( log(x^2-x+1) / (1+x) )
integrate(\(x) log(x^2-x+1) / (1+x), 0, 1)
integrate(\(x) - log(x^2+x+1) / (1+x), 0, 1)$value +
	- pi^2/18 + log(2)*log(3);
# Solution using Li2:
- (pracma::polylog(1/3, 2) - pracma::polylog(-1/3, 2)) + pi^2/18;


### I( log(1-x^3) / (1+x) )
integrate(\(x) log(1-x^3) / (1+x), 0, 1)
# using Li2:
integrate(\(x) 1/2 * log(1-x^6) / (1+x), 0, 1, rel.tol=1E-13)$value + # OK: see above
	- (polylog2(-2) - pracma::polylog(-1/2, 2)) / 2 +
	- pi^2 / 24 + (log(3/2)^2 - log(3)^2) / 4;

### I( log(1+x^3) / (1+x) )
integrate(\(x) log(1+x^3) / (1+x), 0, 1)
# using Li2:
integrate(\(x) 1/2 * log(1-x^6) / (1+x), 0, 1, rel.tol=1E-13)$value + # OK: see above
	+ (polylog2(-2) - pracma::polylog(-1/2, 2)) / 2 +
	+ pi^2 / 24 - (log(3/2)^2 - log(3)^2) / 4;

### Variants:
integrate(\(x) (log(1+x^3) - log(1-x^3)) / (1+x), 0, 1)
polylog2(-2) - pracma::polylog(-1/2, 2) + pi^2/12 - (log(3/2)^2 - log(3)^2)/2;

# Note:
polylog2(-2) # ==
- pracma::polylog(-1/2, 2) - (pi^2/6 + log(2)^2/2);


### Helper

###
integrate(\(x) log(x) / (x^3 - 1), 1, Inf)
pracma::psi(1, 2/3) / 9

###
integrate(\(x) log(1+x^3) / x / (x^3 + 1), 0, Inf)
pi^2/18

###
integrate(\(x) log(abs(1-x^3)) / x / (x^3 + 1), 0, 1)$value +
integrate(\(x) log(abs(1-x^3)) / x / (x^3 + 1), 1, Inf)$value;
pi^2/36

### Li2(1/2)
integrate(\(x) log(1+x) / (1-x) - log(2)/(1-x), 0, 1)
- pi^2/12 + log(2)^2/2
#
integrate(\(x) log(1-x) / (x+1), 0, 1)
- pi^2/12 + log(2)^2/2

#
integrate(\(x) log(1+x) / (1+x), 0, 1)
log(2)^2 / 2;

#
integrate(\(x) log(1+x) * (x^2 + 2) / (1-x^3) - log(2)/(1-x), 0, 1)
(pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 36 +
	- pi^2/12 + log(2)^2/2;


### Derivation:

# from log(x + b) / (x^3 + 1)
b = sqrt(5)
integrate(\(x) 1 / (x+b) / (x^3+1), 0, 1)
integrate(\(x) ((x^2-b*x+b^2)/(x^3+1) - 1/(x+b)) / (b^3-1), 0, 1)
(log(b) - log(b+1) + log(2) / 3 +
	- (digamma(5/6) - digamma(1/3)) * b / 6 +
	+ (digamma(2/3) - digamma(1/6)) * b^2 / 6) / (b^3 - 1);

# from x * log(x + b) / (x^3 + 1)
b = sqrt(5)
integrate(\(x) x / (x+b) / (x^3+1), 0, 1)
integrate(\(x) x * ((x^2-b*x+b^2)/(x^3+1) - 1/(x+b)) / (b^3-1), 0, 1)
integrate(\(x) (b/(x+b) - (b*x^2-b^2*x+1)/(x^3+1)) / (b^3-1), 0, 1)
(b*log(b+1) - b*log(b) - b*log(2) / 3 +
	+ (digamma(5/6) - digamma(1/3)) * b^2 / 6 +
	- (digamma(2/3) - digamma(1/6)) / 6) / (b^3 - 1);

# from x^2 * log(x + b) / (x^3 + 1)
b = sqrt(5)
integrate(\(x) x^2 / (x+b) / (x^3+1), 0, 1)
integrate(\(x) x^2 * ((x^2-b*x+b^2)/(x^3+1) - 1/(x+b)) / (b^3-1), 0, 1)
integrate(\(x) ((b^2*x^2-x+b)/(x^3+1) - b^2/(x+b)) / (b^3-1), 0, 1)
(b^2*log(b) - b^2*log(b+1) + b^2*log(2) / 3 +
	- (digamma(5/6) - digamma(1/3)) / 6 +
	+ (digamma(2/3) - digamma(1/6)) * b / 6) / (b^3 - 1);


# from log(x^3+b^3) / (x^3+1)
# Note: does NOT include factor * 3*b^2;
b = sqrt(5)
integrate(\(x) 1/(x^3+b^3) / (x^3+1), 0, 1)
integrate(\(x) - (1/(x^3+b^3) - 1/(x^3+1)) / (b^3-1), 0, 1)
((digamma(2/3) - digamma(1/6)) * b^2 / 2 +
	- log(b+1) + log(b^2-b+1)/2 - sin(2*pi/3)*pi/3 +
	- 2*sin(2*pi/3)*atan((1/b + cos(2*pi/3))/sin(2*pi/3))) / (3*b^2) / (b^3-1)

# Helper:
integrate(\(x) 1/(x^3+b^3), 0, 1)
integrate(\(x) 1/(x^3+1) / b^2, 0, 1/b)
(log(b+1) + cos(2*pi/3)*log(b^2-b+1) +
	- log(b) - 2*cos(2*pi/3)*log(b) +
	+ 2*sin(2*pi/3)*atan((1/b + cos(2*pi/3))/sin(2*pi/3)) +
	+ sin(2*pi/3)*pi/3) / (3*b^2)

### I( log(x+1) / (x^3-1) )
integrate(\(x) (log(x+1) - log(2)) / (x^3-1), 0, 1)
integrate(\(x) log(x+1) / (x^3-1) - log(2)/3/(x-1), 0, 1)$value +
	+ (digamma(1) - digamma(1/3) - log(3)) * log(2) / 3;
integrate(\(x) -1/3 * x * log(x+1) / (x^2+x+1), 0, 1)$value +
	- (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 54 +
	+ (digamma(1) - digamma(1/3) - log(3)) * log(2) / 3 +
	+ (pi^2/12 - log(2)^2/2) / 3;
# TODO: ?


#
integrate(\(x) log(x+1) / (x-1) - log(2)/(x-1), 0, 1)
pi^2/12 - log(2)^2/2 # Li(1/2)

#
integrate(\(x) Im(log(x+1) / (x + 1/2+1i*sqrt(3)/2)), 0, 1)
- (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) * sqrt(3) / 72;


###########
### Pow = 4

# Note:
# - moved to file: Integrals.Log.Fractions.P4.R;

###
integrate(\(x) x^3 * log(1+x^4) / (1-x^4) - 1/4 * log(2)/(1-x), 0, 1)
log(4) * (digamma(1/2) + Euler) / 8 - (pi^2/6 - log(2)^2)/8;

### Gen: I( x^(p-1) * log(1+x^p) / (1-x^p) )
p = sqrt(5)
integrate(\(x) x^(p-1) * log(1+x^p) / (1-x^p) - 1/p * log(2)/(1-x), 0, 1)
log(p) * (digamma(1/2) + Euler) / (2*p) - (pi^2/6 - log(2)^2)/(2*p);
- log(2)*log(p) / p - (pi^2/6 - log(2)^2)/(2*p);


##################
##################

### I( log((b - x)/(b + x)) / (x * sqrt(1 - x^2)) ) on [0, 1]
# Maths 505: DESTROYING a MONSTER integral using Feynman's technique
# https://www.youtube.com/watch?v=PWtzeKvOVEo

b = sqrt(3)
integrate(\(x) log((b + x)/(b - x)) / (x * sqrt(1 - x^2)), 0, 1)
pi^2 / 2 - pi*acos(1/b)

### Derived:
b = sqrt(3)
integrate(\(x) 1 / ((b^2 - x^2) * sqrt(1 - x^2)), 0, 1)
pi / sqrt(b^2 - 1) / (2*b)


##################

### I( log(1 + (x/(1-x))^phi) / x )
# Maths 505: A golden ratio integral
# https://www.youtube.com/watch?v=xSfPruD0S5o
# (x/(1-x))^phi = y & split int [0,1] + [1, Inf];

phi = (1 + sqrt(5))/2;
integrate(\(x) log(1 + (x/(1-x))^phi) / x, 0, 1)
integrate(\(x) log(1 + (x/(1-x))^(phi-1)) / x, 0, 1)
pi^2 * (phi^2 + 1)/(12*phi);


######################
######################

### I( log(sqrt(x^2 + 1) + x) / (x^2 + 1) )
# Maths 505: An AWESOME log integral with a surprising result
# https://www.youtube.com/watch?v=a2On9rAvqks

integrate(\(x) log(sqrt(x^2 + 1) + x) / (x^2 + 1), 0, Inf)
integrate(\(x) log(sqrt(x^2 + 1) + 1) / (x^2 + 1), 0, Inf)
integrate(\(x) - log(sqrt(x^2 + 1) - x) / (x^2 + 1), 0, Inf)
integrate(\(x) - log(sqrt(x^2 + 1) - 1) / (x^2 + 1), 0, Inf)
2*Catalan


### Derived / Sub-Integrals:
integrate(\(x) log((1 + cos(x)) / cos(x)), 0, pi/2)
2*Catalan

###
# - see also file: Integrals.Log.Trig.R
integrate(\(x) log(1 + cos(x)), 0, pi/2)
integrate(\(x) 4*log(cos(x)), 0, pi/4)$value + pi*log(2)/2
- pi*log(2)/2 + 2*Catalan


#######################
#######################

### Varia:

###
n = sqrt(5)
integrate(function(x) log(2 - x^n) / (2 - x^n)^(1+1/n), lower=0, upper=1)
integrate(function(x) - n*(n+1) * x^n * log(x) / (2 - x^n)^(2+1/n), lower=0, upper=1)
# TODO


#######################
#######################

### 1 / LOG()

###
integrate(\(x) 1/log(x) + 1/(1-x), 0, 1)
Euler


### I( (x^n - 1) / log(x) )
# Michael Penn: My favorite way "around" Feynman's trick -- and a challenge for you!
# https://www.youtube.com/watch?v=EseFhG1QARM

n = sqrt(5)
integrate(\(x) (x^n - 1) / log(x), 0, 1)
log(n+1)

###
n = -1/2
integrate(\(x) x^n / log(x) + 1/(1-x), 0, 1, rel.tol=1E-8)
log(n+1) + Euler

###
n = sqrt(5); p = sqrt(3)
integrate(\(x) x^p * (x^n - 1) / log(x), 0, 1)
log(n+p+1) - log(p+1)


### I( 1 / ((x+1) * log(x)) )
integrate(\(x) 1 / ((x+1) * log(x)) + 1/(1-x^2), 0, 1, rel.tol=1E-8)
- log(gamma(1/2)) - digamma(1/2)/2

### I( (x-1) / ((x+1) * log(x)) )
integrate(\(x) (x-1) / ((x+1) * log(x)), 0, 1, rel.tol=1E-8)
log(pi/2)


### I( x^p * (x^n - 1)^2 / log(x)^2 )
# Maths 505: How Richard Feynman would solve this awesome golden integral
# https://www.youtube.com/watch?v=g2NPdw4ig5M

n = sqrt(5)
integrate(\(x) (x^n - 1)^2 / log(x)^2, 0, 1)
(2*n+1)*log(2*n+1) - (2*n+2)*log(n+1)


### I( 1 / log(x)^2 )
integrate(\(x) 1 / log(x)^2 - 1/(1-x)^2 + 1/(1-x), 0, 1, rel.tol=1E-6)
Euler - 1/2


### I( x^p / log(x)^2 )
p = sqrt(3)
integrate(\(x) x^p / log(x)^2 - 1/log(x)^2 - p/log(x), 0, 1)
(p+1)*log(p+1) - p;

### I( x^p * (x^n - 1)^2 / log(x)^2 )
n = sqrt(5); p = sqrt(3)
integrate(\(x) x^p * (x^n - 1)^2 / log(x)^2, 0, 1)
(2*n+p+1)*log(2*n+p+1) - 2*(n+p+1)*log(n+p+1) + (p+1)*log(p+1)


### Pow = 3

### I( (x^n - 1)^3 / log(x)^3 )
n = sqrt(5)
integrate(\(x) (x^n - 1)^3 / log(x)^3, 0, 1)
1/2*(3*n+1)^2*log(3*n+1) - 3/2*(2*n+1)^2*log(2*n+1) + 3/2*(n+1)^2*log(n+1)

### I( x^p * (x^n - 1)^3 / log(x)^3 )
n = sqrt(5); p = sqrt(3)
integrate(\(x) x^p * (x^n - 1)^3 / log(x)^3, 0, 1)
1/2*(3*n+p+1)^2*log(3*n+p+1) - 3/2*(2*n+p+1)^2*log(2*n+p+1) +
	+ 3/2*(n+p+1)^2*log(n+p+1) - 1/2*(p+1)^2*log(p+1)

# Note:
# shortcut: integrate (f(3*n) - 3*f(2*n) + 3*f(n)) dp;


### I( x^p * (x^n - 1)^4 / log(x)^4 )
n = sqrt(5); p = sqrt(3)
integrate(\(x) x^p * (x^n - 1)^4 / log(x)^4, 0, 1)
1/6*(4*n+p+1)^3*log(4*n+p+1) - 2/3*(3*n+p+1)^3*log(3*n+p+1) +
	+ (2*n+p+1)^3*log(2*n+p+1) - 2/3*(n+p+1)^3*log(n+p+1) +
	+ 1/6*(p+1)^3*log(p+1);


### I( x^p * (x^n - 1)^5 / log(x)^5 )
n = sqrt(5); p = sqrt(3)
integrate(\(x) x^p * (x^n - 1)^5 / log(x)^5, 0, 1)
((5*n+p+1)^4*log(5*n+p+1) - 5*(4*n+p+1)^4*log(4*n+p+1) +
	+ 10*(3*n+p+1)^4*log(3*n+p+1) - 10*(2*n+p+1)^4*log(2*n+p+1) +
	+ 5*(n+p+1)^4*log(n+p+1) - (p+1)^4*log(p+1)) / 24;


##################
##################

### I( 1 / (1 - x*log(x)) ) on [0, 1]
# Flammable Maths: I bet BPRP can not solve this Integral
# https://www.youtube.com/watch?v=p1uJkif5zE0

###
id = seq(0, 30)
integrate(\(x) 1 / (1 - x*log(x)), 0, 1)
sum( (-1)^id * gamma(id + 1) / (id + 1)^(id + 1) )

###
integrate(\(x) 1 / (1 + x*log(x)), 0, 1)
sum( gamma(id + 1) / (id + 1)^(id + 1) )

# TODO: any closed formulas?


##################
##################

### I( 1 / ((x+1) * (|log(x)|^n + b)) )
# Maths 505: A fun little integral exploration
# https://www.youtube.com/watch?v=Hl_ko1RSO3I

b = sqrt(3)
n = 6
integrate(\(x) 1 / ((x+1) * (abs(log(x))^n + b)), 0, Inf)
pi/n * b^(1/n-1) / sin(pi/n)

#
b = sqrt(5)
n = sqrt(7)
# accuracy issues with small n;
integrate(\(x) 1 / ((x+1) * (abs(log(x))^n + b)), 0, Inf, rel.tol = 1E-8, subdivisions = 1025)
pi/n * b^(1/n-1) / sin(pi/n)

#
b = 5
n = sqrt(7)
integrate(\(x) 1 / ((x+1) * (abs(log(x))^n + b)), 0, Inf, rel.tol = 1E-8, subdivisions = 1025)
pi/n * b^(1/n-1) / sin(pi/n)

