########################
###
### Leonard Mada
### [the one and only]
###
### Integrals: Logarithms
### Log-Products
###
### draft v.0.2g


### Constants

Euler   = 0.57721566490153286060651209008240243079;
Catalan = 0.915965594177219015054603514;


### Note:
# abs(log(x)) = log(1/x) for x in [0, 1];


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


### Varia:
integrate(\(x) log(1+x) * log(1-x) / x, 0, 1)
- pracma::zeta(3) * 5/8


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
pracma::psi(1, 1 - (p-1)/n) / ((p-1)*n) + (digamma(1 - (p-1)/n) + Euler) / (p-1)^2


### Gen: I( log(x) * log(1 + x^n) / x^p )
n = sqrt(7); p = sqrt(3);
integrate(\(x) log(x) * log(1 + x^n) / x^p, 0, 1)
(pracma::psi(1, 1 - (p-1)/(2*n)) - 2*pracma::psi(1, 1 - (p-1)/n)) / (2*(p-1)*n) +
	+ (digamma(1 - (p-1)/(2*n)) - digamma(1 - (p-1)/n)) / (p-1)^2;


### Derivation:
# - Feynman technique;
# - for I( log(1 - x^n) / x^p ) see file:
#   Integrals.Log.Fractions.Simple.R;


### Special Cases: p = 1:
# - see section above;
#   I( log(x) * log(1 - x^n) / x )
# - Generalization of Special Case:
#   I( log(x)^s * log(1 - x^n) / x )
#   see section further below;


### Pow = 2:

### I( log(x)^2 * log(1 - x^n) / x^p )
# - using Feynman technique;
p = sqrt(3); n = sqrt(5);
integrate(\(x) log(x)^2 * log(1 - x^n) / x^p, 0, 1)
pracma::psi(2, 1 - (p-1)/n) / ((p-1)*n^2) +
	+ 2*pracma::psi(1, 1 - (p-1)/n) / (n*(p-1)^2) +
	+ 2*(digamma(1 - (p-1)/n) + Euler) / (p-1)^3


### I( log(x)^2 * log(1 + x^n) / x^p )
p = sqrt(3); n = sqrt(5);
integrate(\(x) log(x)^2 * log(1 + x^n) / x^p, 0, 1)
(pracma::psi(2, 1 - (p-1)/(2*n)) - 4*pracma::psi(2, 1 - (p-1)/n)) / (4*(p-1)*n^2) +
	+ 2*(pracma::psi(1, 1 - (p-1)/(2*n)) +
		- 2*pracma::psi(1, 1 - (p-1)/n)) / (2*n*(p-1)^2) +
	+ 2*(digamma(1 - (p-1)/(2*n)) - digamma(1 - (p-1)/n)) / (p-1)^3;


### [old] Variants:
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


### Log: Higher Power

# see section further below for details;
# (example starting from a double Integral)
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

# - a very special case;

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


### on [0, 1/2]
integrate(\(x) log(x) * log(x+1), 0, 1/2)
# TODO: find closed form?
Li2 = - integrate(\(x) x / (2*exp(x) + 1), 0, Inf)$value;
(1/2 - 3/2*log(3/2))*log(2) - 3/2*log(3/2) + Li2 + 1
# alternative for Li2:
id = seq(10000)
Li2 = sum((-1/2)^id / id^2);

###
integrate(\(x) log(x) * log(1-x), 0, 1/2)
integrate(\(x) log(x) * log(1-x), 1/2,1)
1 - pi^2/12


###
integrate(\(x) log(x) * log(1-x) / x^2, 1/2, 1)
pi^2/12 + 2*log(2)^2 - 2*log(2)

#
integrate(\(x) log(x)^2 / x^2, 1/2, 1)
2*log(2)^2 - 4*log(2) + 2


######################
######################

### I( log(x)^s * log(1 - x^n) / x )
# Maths 505: AN EPIC DOUBLE INTEGRAL!!!
# https://www.youtube.com/watch?v=xBbCiOFSqSA&t=192s

###
s = sqrt(5)
integrate(\(x) abs(log(x))^s * log(1-x) / x, 0, 1)
- gamma(s+1) * pracma::zeta(s+2)

# Note: example starts with the Double Integral
# II( log(x)^s / (1 - x*y) ) on [0,1]x[0,1]


### Gen 1: I( log(x)^s * log(1 - x^n) / x )
s = sqrt(5); n = sqrt(3)
integrate(\(x) abs(log(x))^s * log(1 - x^n) / x, 0, 1)
- gamma(s+1) * pracma::zeta(s+2) / n^(s+1)


### Gen 2: I( log(x)^s * log(1 + x^n) / x )
s = sqrt(5); n = sqrt(3)
integrate(\(x) abs(log(x))^s * log(1 + x^n) / x, 0, 1)
gamma(s+1) * pracma::zeta(s+2) * (1 - 1/2^(s+1)) / n^(s+1)


### Simple: I( log(1 - x)^s / x )
s = sqrt(5)
integrate(\(x) abs(log(1-x))^s / x, 0, 1)
gamma(s+1) * pracma::zeta(s+1)


### Reverse:
integrate(\(x) log(1 - x)^2 * log(x) / x, 0, 1, rel.tol=1E-8)
- pracma::zeta(4) / 2

### I( log(1 - x^n)^2 * log(x) / x )
n = sqrt(3)
integrate(\(x) log(1 - x^n)^2 * log(x) / x, 0, 1, rel.tol=1E-8)
- pracma::zeta(4) / (2*n^2)


### I( log(1 - x) * log(1 + x) * log(x) / x )
integrate(\(x) log(1 - x) * log(1 + x) * log(x) / x, 0, 1, rel.tol=1E-8)
# TODO: ???


### I( log(1 - x)^3 * log(x) / x )
integrate(\(x) log(1 - x)^3 * log(x) / x, 0, 1, rel.tol=1E-8)
integrate(\(x) 3/2 * log(1 - x)^2 * log(x)^2 / x, 0, 1, rel.tol=1E-8)
(pi^2*pracma::psi(2,1) - pracma::psi(4,1))/2


### I( log(1 - x) * log(1 + x) )
integrate(\(x) log(1 - x) * log(1 + x), 0, 1)
log(2)^2 - 2*log(2) - pi^2/6 + 2

###  I( log(1 - x) * log(1 + x) / x )
integrate(\(x) log(1 - x) * log(1 + x) / x, 0, 1)
- pracma::zeta(3) * 5/8


####################

### I( x^p * (1 - x)^q * log(1-x) * log(x) )

###
p = 1/3; q = 1/5;
integrate(function(x) x^p * (1 - x)^q, 0, 1)
gamma(p+1)*gamma(q+1) / gamma(p+q+2)

###
p = 1/3; q = 1/5;
integrate(function(x) x^p * (1 - x)^q * log(1-x), 0, 1)
gamma(p+1) * gamma(q+1) * (digamma(q+1) - digamma(p+q+2)) / gamma(p+q+2)

### Lim: q = 0
p = -1/5;
integrate(function(x) x^p * log(1-x), 0, 1)
- (digamma(p+2) + Euler) / (p+1)


###
p = - 1/3; q = 1/5;
integrate(function(x) x^p * (1 - x)^q * log(1-x) * log(x), 0, 1)
gamma(p+1) * gamma(q+1) *
	( (digamma(p+q+2) - digamma(q+1)) * (digamma(p+q+2) - digamma(p+1)) +
		- pracma::psi(1, p+q+2) ) / gamma(p+q+2);

### Special Case: q = 0;
p = - 1/5;
integrate(function(x) x^p * log(1-x) * log(x), 0, 1)
( (digamma(p+2) + Euler) / (p+1) - pracma::psi(1, p+2) ) / (p+1);
# Lim: p -> -1
integrate(function(x) log(1-x) * log(x) / x, 0, 1)
pracma::zeta(3)
# x => y^n =>
n = sqrt(5)
integrate(function(x) log(1 - x^n) * log(x) / x, 0, 1)
pracma::zeta(3) / n^2

# - for a classic approach, see:
#   Flammable Math: Simple Trig Subs won't Help you here...
#   https://www.youtube.com/watch?v=Y6yYSS3YbD8


### I( x^p * (1 - x)^q * log(1-x)^2 * log(x) )
p = - 1/3; q = 1/5;
integrate(function(x) x^p * (1 - x)^q * log(1-x)^2 * log(x), 0, 1)
gamma(p+1) * gamma(q+1) * (
	- (digamma(p+q+2) - digamma(q+1))^2 * (digamma(p+q+2) - digamma(p+1)) +
	- pracma::psi(2, p+q+2) +
	+ pracma::psi(1, p+q+2) * (3*digamma(p+q+2) - 2*digamma(q+1) - digamma(p+1)) +
	- pracma::psi(1, q+1) * (digamma(p+q+2) - digamma(p+1))
	) / gamma(p+q+2);

# Special Case: q = 0
integrate(function(x) x^p * log(1-x)^2 * log(x), 0, 1)
( - (digamma(p+2) + Euler)^2 - pi^2/6 +
	- pracma::psi(2, p+2) * (p+1) +
	+ pracma::psi(1, p+2) * (2*digamma(p+2) + 2*Euler) * (p+1) +
	+ pracma::psi(1, p+2)
	) / (p+1)^2;

### I( x^p * log(1-x)^2 * log(x)^2 )
p = 1/7
integrate(function(x) x^p * log(1-x)^2 * log(x)^2, 0, 1)
(   - pracma::psi(3, p+2) +
	+ pracma::psi(2, p+2) * (digamma(p+2) + Euler) * 2 +
	+ pracma::psi(1, p+2) * pracma::psi(1, p+2) * 2
	) / (p+1) +
	- 2*( - (digamma(p+2) + Euler)^2 +
	- pracma::psi(2, p+2) * (p+1) +
	+ pracma::psi(1, p+2) * (2*digamma(p+2) + 2*Euler) * (p+1) +
	+ pracma::psi(1, p+2) - pi^2/6
	) / (p+1)^3;


### Lim: x -> -1
integrate(function(x) log(1-x)^2 * log(x)^2 / x, 0, 1)
pracma::psi(1,1)*pracma::psi(2,1)*2 - pracma::psi(4,1)/3
(pi^2*pracma::psi(2,1) - pracma::psi(4,1))/3


####################

### I( (1 + x)^p * (1 - x)^q * log(1+x) * log(1-x) )

###
p = 1/5; q = - 1/7;
# p = q = 1/7
integrate(function(x) (1 + x)^p * (1 - x)^q, lower=-1, upper=1)
2^(p+q+1) * gamma(p+1)*gamma(q+1) / gamma(p+q+2)

###
p = 1/5; q = - 1/7;
# p = q = 1/7
integrate(function(x) (1 + x)^p * (1 - x)^q * log(1+x), -1, 1)
2^(p+q+1) * gamma(p+1)*gamma(q+1) *
	(log(2) + digamma(p+1) - digamma(p+q+2)) / gamma(p+q+2)

###
p = 1/5; q = - 1/7;
# p = q = 1/7
integrate(function(x) (1 + x)^p * (1 - x)^q * log(1+x) * log(1-x), -1, 1)
- 2^(p+q+1) * gamma(p+1)*gamma(q+1) * pracma::psi(1, p+q+2) / gamma(p+q+2) +
2^(p+q+1) * gamma(p+1)*gamma(q+1) *
	(log(2) + digamma(q+1) - digamma(p+q+2)) *
	(log(2) + digamma(p+1) - digamma(p+q+2)) / gamma(p+q+2)

# Special Case: p = q =>
2^(2*p+1) * gamma(p+1)^2 *
	((log(2) + digamma(p+1) - digamma(2*p+2))^2 - pracma::psi(1, 2*p+2)) / gamma(2*p+2)
# p = q = 0 =>
2 * ((log(2) - 1)^2 - pracma::psi(1, 2))
2 * (log(2)^2 - 2*log(2) - pi^2/6 + 2)


####################
####################

### I( x^p * log(x) * log(x^n + 1) / (x^n + 1)^k )
p = sqrt(3); n = sqrt(5); k = sqrt(2) + sqrt(3);
integrate(function(x) x^p * log(x) * log(x^n + 1) / (x^n + 1)^k, lower=0, upper=Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) *
	(pracma::psi(1, k - (p+1)/n) + (digamma(k - (p+1)/n) - digamma((p+1)/n)) *
	(digamma(k - (p+1)/n) - digamma(k)) ) / gamma(k) / n^2;


### Special Cases:

### Lim: p -> -1
n = sqrt(5); k = sqrt(2) + sqrt(3);
integrate(function(x) log(x) * log(x^n + 1) / (x * (x^n + 1)^k), lower=0, upper=Inf)
- (pracma::psi(2, k)/2 + pracma::psi(1, k)*(digamma(k) + Euler)) / n^2


### Lim: p = -1; n = k = 1;
# Michael Penn: An integral with many logs
# https://www.youtube.com/watch?v=5blBvR7T24o
integrate(\(x) log(x) * log(x + 1) / (x*(x+1)), 0, Inf)
pracma::zeta(3)


### p = -2
n = sqrt(5); k = sqrt(2) + sqrt(3);
integrate(function(x) log(x) * log(x^n + 1) / (x^2 * (x^n + 1)^k), lower=0, upper=Inf)
gamma(-1/n) * gamma(k + 1/n) *
	(pracma::psi(1, k + 1/n) + (digamma(k + 1/n) - digamma(-1/n)) *
	(digamma(k + 1/n) - digamma(k)) ) / gamma(k) / n^2;
# Lim: k -> 0
integrate(function(x) log(x) * log(x^n + 1) / x^2, 0, Inf)
- gamma(-1/n) * gamma(1/n) * (pi /tan(pi/n) + n) / n^2
pi/sin(pi/n) * (pi / tan(pi/n) + n) / n

### I( log(x) * log(x^n + 1) / x^p )
# p > 1; n >= floor(p);
n = sqrt(5); p = sqrt(7);
integrate(function(x) log(x) * log(x^n + 1) / x^p, 0, Inf)
pi/sin(pi*(p-1)/n) * (pi / tan(pi*(p-1)/n) + n/(p-1)) / ((p-1)*n)


###
# Maths 505: A savage logarithmic integral!
# https://www.youtube.com/watch?v=gTBARM502JY

integrate(function(x) log(x) * log(x^2 + 1) / x^2, lower=0, upper=1)
pi/2 - log(2) - 2*Catalan;

#
integrate(\(x) log(1-x)^2 / x, 0, 1)
2*pracma::zeta(3)


####################
####################

### Other

dzeta = function(x, dx=1E-8) {
	(pracma::zeta(x + dx) - pracma::zeta(x)) / dx;
}

### I( log(1 - x)^s * log(- log(1 - x)) / x^3 )
s = sqrt(11); # s > 3
integrate(\(x) abs(log(1 - x))^s * log(-log(1 - x)) / x^3, 0, 1, rel.tol=1E-8)
(digamma(s+1) * (pracma::zeta(s) + pracma::zeta(s-1)) +
	+ dzeta(s) + dzeta(s-1)) * gamma(s+1) / 2


### I( x^p * log(x)^s * log(- log(x)) )
s = sqrt(5); p = sqrt(7);
integrate(\(x) x^p * abs(log(x))^s * log(-log(x)), 0, 1)
gamma(s+1) * (digamma(s+1) - log(p+1)) / (p+1)^(s+1)


### Special Case: p = 0
s = sqrt(7)
integrate(\(x) abs(log(x))^s * log(-log(x)), 0, 1)
gamma(s+1) * digamma(s+1)

