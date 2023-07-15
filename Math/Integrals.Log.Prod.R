########################
###
### Leonard Mada
### [the one and only]
###
### Integrals: Logarithms
### Log-Products
###
### draft v.0.1c


### Constants

Catalan = 0.915965594177219015054603514;


################
################

### I( log(x) * log(1 - x^n) / x )
# Maths 505: A beautiful log-trig integral featuring an important constant
# https://www.youtube.com/watch?v=2bwvuWsQSEY
# [the intermediate integral]

###
integrate(\(x) log(x) * log(1 - x^2) / x, 0, 1)
pracma::zeta(3) / 4

### Gen:
n = sqrt(5)
integrate(\(x) log(x) * log(1 - x^n) / x, 0, 1)
pracma::zeta(3) / n^2

###
integrate(\(x) log(x) * log(1 + x) / x, 0, 1)
- pracma::zeta(3) * 3/4

### Gen:
n = sqrt(5)
integrate(\(x) log(x) * log(1 + x^n) / x, 0, 1)
- pracma::zeta(3) * 3/4 / n^2

### Derived:
integrate(\(x) log(x) * log(x^2 + x + 1) / x, 0, 1)
- pracma::zeta(3) * 8/9
#
integrate(\(x) log(x) * log(x^2 - x + 1) / x, 0, 1)
pracma::zeta(3) * 2/3


### Series: Pow = 2
# P() / x^2

###
integrate(\(x) log(x) * log(1 - x^2) / x^2, 0, 1)
pi^2 / 4 - 2*log(2)

###
integrate(\(x) log(x) * log(1 + x^2) / x^2, 0, 1)
pi/2 - log(2) - 2*Catalan


### Gen: I( log(x) * log(1 - x^n) / x^p )
n = sqrt(7); p = sqrt(3);
integrate(\(x) log(x) * log(1 - x^n) / x^p, 0, 1)
pracma::psi(1, (n-p+1)/n) / ((p-1)*n) + (digamma((n-p+1)/n) + Euler) / (p-1)^2


###
n = sqrt(3)
integrate(\(x) log(x) * log(1 - x^(2*n)) / x^(n+1), 0, 1)
(pi^2 / 4 - 2*log(2)) / n^2

# Partial Variants:
# =>
integrate(\(x) log(x) * log(1 + x) / x^2 - log(x)/x + log(x)/2, 0, 1)
pi^2/12 - 2*log(2) + 0.5
integrate(\(x) log(x) * log(1 - x) / x^2 + log(x)/x + log(x)/2, 0, 1)
pi^2/6 - 1.5


### Base-Variants:

### Other: Simple
integrate(function(x) log(x) * log(1 - x), 0, 1)
2 - pi^2/6
#
integrate(function(x) log(x) * log(1 - x), 0, 1/2)
integrate(function(x) log(x) * log(1 - x), 1/2, 1)
integrate(function(x) log(1/2 - x) * log(1/2 + x), 0, 1/2)
1 - pi^2/12


### Base: Pow = 1
integrate(\(x) log(x) * log(1 - x) / x, 0, 1)
pracma::zeta(3)
# =>
integrate(\(x) log(x) * log(1 + x) / x, 0, 1)
-3/4 * pracma::zeta(3)


### Higher Power

# =>
integrate(\(x) log(x)^3 * log(1 - x) / x, 0, 1)
6 * pracma::zeta(5)
#
integrate(\(x) log(x)^3 * log(1 + x) / x, 0, 1)
-45/8 * pracma::zeta(5)

### Generalized:
p = sqrt(5)
integrate(\(x) abs(log(x))^p * log(1 - x) / x, 0, 1)
- gamma(p+1) * pracma::zeta(p+2)


### I( |log(x)|^p * log(1 - x^n) / x )
p = sqrt(3); n = sqrt(7);
integrate(\(x) abs(log(x))^p * log(1 - x^n) / x, 0, 1, rel.tol=1E-8)
- gamma(p+1) * pracma::zeta(p+2) / n^(p+1)


### I( |log(x)|^p * log(1 + x^n) / x )
# workout of log(1 - x^2) =>
p = sqrt(7)
integrate(\(x) abs(log(x))^p * log(1 + x) / x, 0, 1, rel.tol=1E-8)
gamma(p+1) * pracma::zeta(p+2) * (1 - 1/2^(p+1))
#
n = sqrt(5)
integrate(\(x) abs(log(x))^p * log(1 + x^n) / x, 0, 1, rel.tol=1E-8)
gamma(p+1) * pracma::zeta(p+2) * (1 - 1/2^(p+1)) / n^(p+1)

#
integrate(\(x) abs(log(x))^p * log(1 + x + x^2) / x, 0, 1, rel.tol=1E-8)
gamma(p+1) * pracma::zeta(p+2) * (1 - 1/3^(p+1))
#
integrate(\(x) abs(log(x))^p * log(1 - x + x^2) / x, 0, 1, rel.tol=1E-8)
- gamma(p+1) * pracma::zeta(p+2) * (1 - 1/2^(p+1)) * (1 - 1/3^(p+1))


### Variants:
# I( log(x)^s * log((1-x)/(1+x) / x )

### Gen:
s = sqrt(3)
integrate(\(x) abs(log(x))^s * log((1 - x) / (1 + x)) / x, 0, 1, rel.tol=1E-8)
- gamma(s+1) * pracma::zeta(s+2) * (2 - 1/2^(s+1))


###
integrate(\(x) log(x) * log((1 - x) / (1 + x)) / x, 0, 1)
7/4 * pracma::zeta(3)

###
integrate(\(x) log(x)^2 * log((1 - x) / (1 + x)) / x, 0, 1, rel.tol = 1E-8)
- 15/4 * pracma::zeta(4)

###
integrate(\(x) log(x)^3 * log((1 - x) / (1 + x)) / x, 0, 1, rel.tol = 1E-8)
93/8 * pracma::zeta(5)

###
integrate(\(x) log(x)^4 * log((1 - x) / (1 + x)) / x, 0, 1, rel.tol = 1E-8)
- 189/4 * pracma::zeta(6)


### Simple:
# see also file;
# Integrals.Log.Fractions.Other.R;
integrate(\(x) log((1 - x) / (x + 1)), 0, 1)
- 2*log(2)

###
integrate(\(x) log(x) * log((1 - x) / (x + 1)), 0, 1)
2*log(2) - pi^2/12


####################
####################

### I( log(x) * log(x+1) )
# Maths 505: A surprisingly wonderful log integral
# https://www.youtube.com/watch?v=v6iNE-heksU

integrate(\(x) log(x) * log(x+1), 0, 1)
2 - 2*log(2) - pi^2/12

###
a = sqrt(5)
integrate(\(x) log(x) * log(x+a), 0, a)
integrate(\(x) a*log(a*x) * log(a*x+a), 0, 1)
a*(2 + 2*log(2)*log(a) - 2*log(2) - 2*log(a) + log(a)^2 - pi^2/12)


###
a = 1/sqrt(3)
integrate(\(x) log(x) * log(x+1), 0, a)
# TODO:
id = seq(10000)
log(a) * sum((-a)^(id+1) / (id*(id+1))) +
	sum((-1)^id * a^(id+1) / (id*(id+1)^2))

###
integrate(\(x) log(x) * log(x+1), 0, 1/2)
# TODO: find closed form?
Li2 = - integrate(\(x) x / (2*exp(x) + 1), 0, Inf)$value;
(1/2 - 3/2*log(3/2))*log(2) - 3/2*log(3/2) + Li2 + 1
# alternative for Li2:
id = seq(10000)
Li2 = sum((-1/2)^id / id^2);


###
integrate(\(x) log(x) * log(1-x) / x^2, 1/2, 1)
pi^2/12 + 2*log(2)^2 - 2*log(2)

#
integrate(\(x) log(x)^2 / x^2, 1/2, 1)
2*log(2)^2 - 4*log(2) + 2

