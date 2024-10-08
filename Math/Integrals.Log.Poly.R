

### Integrals: Log(Poly)


################
################

###

# Maths 505: A RIDICULOUSLY AWESOME INTEGRAL:
#   int (0,infty) ln^2(x) (1+x^2)/(1+x^4)
# https://www.youtube.com/watch?v=9DV0R5W9feA
# - this is actually a trivial integral;
#   (I misinterpreted the formula in the thumbnail)
# - below are a few more interesting types;


###
integrate(\(x) log(x^2*(x^2+1) / (x^4 + 1))^2, 0, Inf)
# TODO

### Varia:

### I( log(x^2*(x^2+1) / (x^4 + 1)) )
integrate(\(x) log(x^2*(x^2+1) / (x^4 + 1)), 0, Inf)
gamma(1/4)*gamma(-1/4)/4 - gamma(1/2)*gamma(-1/2)/2;
pi*(1/sin(pi/2) - 1/sin(pi/4))

### Limit: k -> 0
k = 1E-4
gamma(1/2) * gamma(k - 1/2) / gamma(k) / 2 *
	(digamma(1/2) - digamma(k - 1/2)) +
gamma(1/2) * gamma(k - 1/2) / gamma(k) / 2 *
	(digamma(k) - digamma(k - 1/2)) +
- gamma(1/4) * gamma(k - 1/4) / gamma(k) / 4 *
	(digamma(k) - digamma(k - 1/4));


### I( log(x^3 * (x^3 + 1) / (x^6 + 1)) )
integrate(\(x) log(x^3*(x^3+1) / (x^6 + 1)), 0, Inf)
gamma(1/6)*gamma(-1/6)/6 - gamma(1/3)*gamma(-1/3)/3;
pi*(1/sin(pi/3) - 1/sin(pi/6))


### I( log(x^n * (x^n + 1) / (x^(2*n) + 1)) )
n = sqrt(3)
integrate(\(x) log(x^n * (x^n+1) / (x^(2*n) + 1)), 0, Inf)
pi*(1/sin(pi/n) - 1/sin(pi/(2*n)))


### I( log((x^2+1)^2 / (x^4 + 1)) )
integrate(\(x) log((x^2+1)^2 / (x^4 + 1)), 0, Inf)
gamma(1/4) * gamma(- 1/4)/4 - gamma(1/2) * gamma(- 1/2)
pi*(2/sin(pi/2) - 1/sin(pi/4))


### I( log((x^n+1)^2 / (x^(2*n) + 1)) )
n = sqrt(3)
integrate(\(x) log((x^n + 1)^2 / (x^(2*n) + 1)), 0, Inf)
pi*(2/sin(pi/n) - 1/sin(pi/(2*n)))


### Derivations:

### Limit: k -> 0
# Note: n = 3 for both log(x) & log(x^3 + 1);
k = 1E-4
# library(Rmpfr); k = mpfr("1E-8", 240);
# Note: fractions need proper computation as well!
gamma(1/3) * gamma(k - 1/3) / gamma(k) / 3 *
	(digamma(1/3) - digamma(k - 1/3)) +
gamma(1/3) * gamma(k - 1/3) / gamma(k) / 3 *
	(digamma(k) - digamma(k - 1/3)) +
- gamma(1/6) * gamma(k - 1/6) / gamma(k) / 6 *
	(digamma(k) - digamma(k - 1/6));


### Helper:
p = sqrt(3); n = sqrt(7); k = sqrt(5);
integrate(\(x) x^p / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n
# D(p)
integrate(\(x) x^p * log(x) / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) * (digamma((p+1)/n) - digamma(k - (p+1)/n)) / gamma(k) / n^2
# D2(p)
integrate(\(x) x^p * log(x)^2 / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n^3 *
	(pracma::psi(1, (p+1)/n) + pracma::psi(1, k - (p+1)/n) +
	+ (digamma((p+1)/n) - digamma(k - (p+1)/n))^2)
# D3(p)
integrate(\(x) x^p * log(x)^3 / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n^4 *
	(pracma::psi(2, (p+1)/n) - pracma::psi(2, k - (p+1)/n) +
	3*(pracma::psi(1, (p+1)/n) + pracma::psi(1, k - (p+1)/n)) *
	(digamma((p+1)/n) - digamma(k - (p+1)/n)) +
	+ (digamma((p+1)/n) - digamma(k - (p+1)/n))^3)

# D(k)
integrate(\(x) x^p * log(x^n + 1) / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n *
	(digamma(k) - digamma(k - (p+1)/n))
# Case: p = 0
gamma(1/n) * gamma(k - 1/n) / gamma(k) / n *
	(digamma(k) - digamma(k - 1/n))


################
################

### I( log(x^2 + 2*sin(a)*x + 1) / x ) on [0,1]
# Maths 505: A beautifully unusual integral calculation
# https://www.youtube.com/watch?v=Na5kpJY0s9g

a = asin(sqrt(3)/4)
integrate(\(x) log(x^2 - 2*sin(a)*x + 1) / x, 0, 1)
- a^2/2 - pi/2*a + pi^2 / 24

###
a = asin(1/pi)
integrate(\(x) log(x^2 - 2*sin(a)*x + 1) / x, 0, 1)
- a^2/2 - pi/2*a + pi^2 / 24

###
a = asin(1/pi)
integrate(\(x) log(x^2 + 2*sin(a)*x + 1) / x, 0, 1)
- a^2/2 + pi/2*a + pi^2 / 24

###
a = asin(pi + 0i)
integrate(\(x) log(x^2 + 2*Re(sin(a))*x + 1) / x, 0, 1)
- a^2/2 + pi/2*a + pi^2 / 24


### Special: x^n + 1
# see file: Integrals.Log.Fractions.R
n = sqrt(5); p = sqrt(3);
integrate(function(x) log(x^n + 1) / x^(p+1), 0, 1)
(digamma(((n-p)/n + 1)/2) - digamma((n-p)/(2*n))) / (2*p) - log(2)/p;

# p -> 0
p = 1E-4
(digamma(1 - p/(2*n)) - digamma(1/2 - p/(2*n))) / (2*p) - log(2)/p;
pi^2/(12*n)

###
integrate(\(x) log(x^3 + 1) / x, 0, 1)
a = 2*pi/3
- a^2 / 2 + pi^2 / 4

# Deriv:
- (pi/2 - a)^2/2 + pi/2*(pi/2 - a) + pi^2 / 24 + pi^2/12

###
n = sqrt(11)
integrate(\(x) log(x^n + 1) / x, 0, 1)
pi^2 / (12*n)

